import mod


from collections import Counter
from overlay_graphs.canonicalisation import CanonicalGraph, CanonicalRule, GraphCanonicaliser
from overlay_graphs.filtered_rule import FilteredRule
from overlay_graphs.label_parser import abstract_vertex_term_details, is_term
from overlay_graphs.mechanism import Mechanism, Step
from overlay_graphs.rule_builder import RuleBuilder
from typing import Iterable, List, Optional, Set, Tuple


def _report_difference(direction: str, difference: Counter[CanonicalGraph, int], index: int, verbosity: int):
    if len(difference) <= 0:
        return

    report = {str(graph): count for graph, count in difference.items()}
    if verbosity >= 2:
        print(f"Found a {direction} diff in step {index}")
        print(f"{report}")


class ExtendableCanonicalRule(CanonicalRule):
    def __init__(self, rule: mod.Rule, canonicaliser: GraphCanonicaliser):
        super().__init__(rule, canonicaliser)

        self._extra_context: List[CanonicalGraph] = []

    @property
    def left(self) -> Tuple[CanonicalGraph]:
        return tuple(sorted(list(super().left) + self._extra_context, key=lambda g: g.canonical_smiles))

    @property
    def right(self) -> Tuple[CanonicalGraph]:
        return tuple(sorted(list(super().right) + self._extra_context, key=lambda g: g.canonical_smiles))

    def _to_rule_builder(self) -> RuleBuilder:
        rule_builder = RuleBuilder.from_rule(self.rule)

        for graph in self._extra_context:
            rule_builder.add_context_graph(graph.graph)

        return rule_builder

    def add_context(self, graphs: Iterable[CanonicalGraph]):
        self._extra_context.extend(graphs)

    def to_gml(self) -> str:
        return self._to_rule_builder().to_gml()

    def to_mod_rule(self) -> mod.Rule:
        return self._to_rule_builder().to_mod_rule()


class MechanismSanitiser:
    def __init__(self, small_molecules: Iterable[mod.Graph] = tuple(), preserve_peptide_chain_positions: bool = False,
                 ignore_dative_bonds: bool = False):
        self._label_settings: mod.LabelSettings = mod.LabelSettings(mod.LabelType.String, mod.LabelRelation.Isomorphism)
        self._canonicaliser: GraphCanonicaliser = GraphCanonicaliser()

        self._small_molecules: Set[CanonicalGraph] = set()
        self.add_small_molecules(small_molecules)

        self._preserve_peptide_chain_positions: bool = preserve_peptide_chain_positions
        self._ignore_dative_bonds: bool = ignore_dative_bonds

    def _get_preserve_peptide_chain_positions(self) -> bool:
        return self._preserve_peptide_chain_positions

    def _set_preserve_peptide_chain_positions(self, value: bool):
        self._preserve_peptide_chain_positions = value

    def _get_ignore_dative_bonds(self) -> bool:
        return self._ignore_dative_bonds

    def _set_ignore_dative_bonds(self, value: bool):
        self._ignore_dative_bonds = value

    preserve_peptide_chain_positions = property(_get_preserve_peptide_chain_positions,
                                                _set_preserve_peptide_chain_positions)

    ignore_dative_bonds = property(_get_ignore_dative_bonds, _set_ignore_dative_bonds)

    @staticmethod
    def _is_hydrogen_on_abstract_atom(vertex: mod.Rule.Vertex) -> bool:
        if vertex.left.stringLabel != "H":
            return False

        return any(not edge.left.isNull() and (is_term(edge.source.left.stringLabel) or
                                               is_term(edge.target.left.stringLabel)) for edge in vertex.incidentEdges)

    @staticmethod
    def _is_dative_bond(edge: mod.Rule.Edge) -> bool:
        if not edge.left.isNull() and edge.left.stringLabel == ">":
            return True

        if not edge.right.isNull() and edge.right.stringLabel == ">":
            return True

        return False

    @staticmethod
    def _graph_tuple_difference(first_tuple: Tuple[CanonicalGraph], second_tuple: Tuple[CanonicalGraph]) ->\
            Counter[CanonicalGraph, int]:
        return Counter(first_tuple) - Counter(second_tuple)

    @staticmethod
    def _rule_difference(first_rule: ExtendableCanonicalRule, second_rule: ExtendableCanonicalRule) ->\
            Counter[CanonicalGraph, int]:
        return MechanismSanitiser._graph_tuple_difference(first_rule.right, second_rule.left)

    @staticmethod
    def _rule_codifference(first_rule: ExtendableCanonicalRule, second_rule: ExtendableCanonicalRule) ->\
            Counter[CanonicalGraph, int]:
        return MechanismSanitiser._graph_tuple_difference(second_rule.left, first_rule.right)

    def _sanitise_step(self, step: Step) -> (mod.Rule, mod.Rule):
        sanitised_rule = FilteredRule(step.rule)
        abstract_rule = FilteredRule(step.rule)

        for vertex in step.rule.vertices:
            label = abstract_vertex_term_details(vertex.left.stringLabel, self.preserve_peptide_chain_positions)

            if self._is_hydrogen_on_abstract_atom(vertex):
                continue

            sanitised_rule.include_vertex(vertex, left_label=label)
            abstract_rule.include_vertex(vertex, left_label=label)

        for edge in step.rule.edges:
            if edge.source not in sanitised_rule.used_vertices or edge.target not in sanitised_rule.used_vertices:
                continue

            sanitised_rule.include_edge(edge)

            if self.ignore_dative_bonds and self._is_dative_bond(edge):
                continue

            abstract_rule.include_edge(edge)

        return abstract_rule.to_mod_rule(), sanitised_rule.to_mod_rule()

    def _canonicalise_steps(self, steps: List[Step], verbosity: int) ->\
            Iterable[Tuple[ExtendableCanonicalRule, ExtendableCanonicalRule]]:
        for step in steps:
            if step.rule is None:
                return

            if verbosity >= 2:
                print(f"sanitizing rule {step.rule}")

            abstract_rule, sanitised_rule = self._sanitise_step(step)
            yield ExtendableCanonicalRule(abstract_rule, self._canonicaliser),\
                  ExtendableCanonicalRule(sanitised_rule, self._canonicaliser)

    def _sanitise_steps(self, steps: List[Step], verbosity: int) -> Optional[List[ExtendableCanonicalRule]]:
        canonical_rules = list(self._canonicalise_steps(steps, verbosity))

        for index in range(1, len(canonical_rules)):
            first_abstract, first_rule = canonical_rules[index - 1]
            second_abstract, second_rule = canonical_rules[index]

            difference = self._rule_difference(first_abstract, second_abstract)

            if len(difference) > 0:
                codifference = self._rule_codifference(first_abstract, second_abstract)

                if len(codifference) > 0 and not set(difference).issubset(self._small_molecules) and\
                        not set(codifference).issubset(self._small_molecules):
                    if verbosity >= 2:
                        print(f"Failed to sanitise step {index}")
                    return None

                _report_difference("forward", difference, index, verbosity)

            second_abstract.add_context(difference.elements())
            second_rule.add_context(difference.elements())

        for index in range(1, len(canonical_rules)):
            first_abstract, first_rule = canonical_rules[len(canonical_rules) - index]
            second_abstract, second_rule = canonical_rules[len(canonical_rules) - index - 1]

            difference = self._rule_codifference(second_abstract, first_abstract)

            _report_difference("backward", difference, index, verbosity)

            second_abstract.add_context(difference.elements())
            second_rule.add_context(difference.elements())

        if verbosity >= 2:
            print(str(sum(graph.graph.numVertices for graph in canonical_rules[0][1].left) if
                      len(canonical_rules) > 0 else 0))

        return [sanitised_rule for abstract_rule, sanitised_rule in canonical_rules]

    def add_small_molecule(self, molecule: mod.Graph):
        self._small_molecules.add(self._canonicaliser.canonicalise_graph(molecule))

    def add_small_molecules(self, molecules: Iterable[mod.Graph]):
        for molecule in molecules:
            self.add_small_molecule(molecule)

    def sanitise_mechanism(self, mechanism: Mechanism, verbosity: int = 0) -> Optional[Mechanism]:
        canonical_steps = self._sanitise_steps(mechanism.steps, verbosity)

        if canonical_steps is None:
            return None

        return Mechanism(mechanism.entry, mechanism.number,
                         [Step(mechanism.entry, mechanism.number, index + 1, canonical_step.to_mod_rule())
                          for index, canonical_step in enumerate(canonical_steps)])

    def sanitise_mechanisms(self, mechanisms: Iterable[Mechanism], verbosity: int = 0) -> Iterable[Mechanism]:
        for mechanism in mechanisms:
            sanitised_mechanism = self.sanitise_mechanism(mechanism, verbosity)

            if sanitised_mechanism is not None:
                yield sanitised_mechanism
