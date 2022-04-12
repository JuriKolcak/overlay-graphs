import mod
import networkx as nx
import re


from collections import Counter
from networkx.algorithms.isomorphism import GraphMatcher
from overlay_graphs.networkx_converter import get_component_graphs, get_rule_component_graphs_with_nx
from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple


_coordinating_nonmetals = re.compile(r"^[NOS][^a-z]*$")
_label_settings = mod.LabelSettings(mod.LabelType.String, mod.LabelRelation.Isomorphism)


def _identify_metal_ions(graphs: Dict[mod.Graph, nx.Graph]) -> Dict[mod.Graph, Counter[str]]:
    metal_ions = {}

    for graph, nx_graph in graphs.items():
        metal_ions[graph] = Counter(data["label"] for node, data in nx_graph.nodes(data=True) if
                                    re.match(_coordinating_nonmetals, data["label"]) is None and
                                    any(nx_graph.edges[edge]["label"] == ":" for edge in nx_graph.edges(node)))

    return metal_ions


def _decompose_coordination_bonds(graphs: Dict[mod.Graph, nx.Graph]) -> Dict[mod.Graph, Dict[mod.Graph, nx.Graph]]:
    decompositions = {}

    for graph, nx_graph in graphs.items():
        decompositions[graph] = get_component_graphs(
            nx.subgraph_view(nx_graph, filter_edge=lambda source, target: nx_graph
                             .edges[source, target]["label"] != ":"),
            lambda n: n.id)

    return decompositions


def _compute_matches(first_graphs: Dict[mod.Graph, nx.Graph], second_graphs: Dict[mod.Graph, nx.Graph],
                     reaction_center: Set[int]) -> Iterable['Isomorphism']:
    sorted_first = sorted(first_graphs)

    matched = False
    for index, first_graph in enumerate(sorted_first):
        matched = False

        for second_graph, second_nx_graph in second_graphs.items():
            if first_graph.isomorphism(second_graph, labelSettings=_label_settings) <= 0:
                continue

            matched = True
            matcher = GraphMatcher(first_graphs[first_graph], second_nx_graph,
                                   lambda node1, node2: node1["label"] == node2["label"],
                                   lambda edge1, edge2: edge1["label"] == edge2["label"])

            isomorphism = Isomorphism.from_match(reaction_center, first_graph, first_graphs[first_graph], second_graph,
                                                 second_nx_graph, {first_node.id: second_node.id for
                                                                   first_node, second_node in
                                                                   next(matcher.isomorphisms_iter()).items()})

            remaining_first_graphs = {graph: first_graphs[graph] for graph in sorted_first[index + 1:]}

            remaining_second_graphs = dict(second_graphs)
            del remaining_second_graphs[second_graph]

            yield from [isomorphism + remaining_match for remaining_match in
                        _compute_matches(remaining_first_graphs, remaining_second_graphs, reaction_center)]

        if matched:
            break

    if not matched:
        yield Isomorphism(reaction_center)


def _match_metal_ions(first_metal_ions: Dict[mod.Graph, Counter[str]],
                      second_metal_ions: Dict[mod.Graph, Counter[str]]) -> Iterable[Dict[mod.Graph, List[mod.Graph]]]:
    if all(len(ions) == 0 for graph, ions in first_metal_ions.items()) or\
            all(len(ions) == 0 for graph, ions in second_metal_ions.items()):
        yield {}
        return

    sorted_first = sorted(first_metal_ions)

    for index, first_graph in enumerate(sorted_first):
        first_ions = first_metal_ions[first_graph]
        if len(first_ions) == 0:
            continue

        for second_graph, second_ions in second_metal_ions.items():
            if first_ions - second_ions == first_ions:
                continue

            remaining_first_ions = {graph: first_metal_ions[graph] for graph in sorted_first[index + 1:]}
            remaining_first_ions[first_graph] = first_ions - second_ions

            remaining_second_ions = dict(second_metal_ions)
            remaining_second_ions[second_graph] = second_ions - first_ions

            for ion_match in _match_metal_ions(remaining_first_ions, remaining_second_ions):
                match = dict(ion_match)
                if first_graph not in match:
                    match[first_graph] = []

                match[first_graph].append(second_graph)
                yield match


def _complete_match(match: 'Isomorphism', first_graphs: Dict[mod.Graph, nx.Graph],
                    second_graphs: Dict[mod.Graph, nx.Graph], reaction_center: Set[int]) -> Iterable['Isomorphism']:
    unmatched_first = {graph: nx_graph for graph, nx_graph in first_graphs.items() if graph not in match.first}
    unmatched_second = {graph: nx_graph for graph, nx_graph in second_graphs.items() if graph not in match.second}

    if len(unmatched_first) == 0 and len(unmatched_second) == 0:
        yield match
        return

    first_metal_ions = _identify_metal_ions(unmatched_first)
    second_metal_ions = _identify_metal_ions(unmatched_second)

    first_coordination_decompositions = _decompose_coordination_bonds(unmatched_first)
    second_coordination_decompositions = _decompose_coordination_bonds(unmatched_second)

    for ion_match in _match_metal_ions(first_metal_ions, second_metal_ions):
        ion_aware_matches = [match]
        for first_graph, ion_matched_second_graphs in ion_match.items():
            second_subgraphs = {}
            for second_graph in ion_matched_second_graphs:
                second_subgraphs.update(second_coordination_decompositions[second_graph])

            new_matches = []

            for submatch in _compute_matches(first_coordination_decompositions[first_graph], second_subgraphs,
                                             reaction_center):
                new_matches.extend(ion_aware_match + submatch for ion_aware_match in ion_aware_matches if
                                   ion_aware_match.compatible(submatch))
            if len(new_matches) > 0:
                ion_aware_matches = new_matches

        for ion_aware_match in ion_aware_matches:
            still_unmatched_first = {}
            for graph, decomposition in first_coordination_decompositions.items():
                still_unmatched_first.update({component: nx_component for component, nx_component in
                                              decomposition.items() if component not in ion_aware_match.first})

            still_unmatched_second = {}
            for graph, decomposition in second_coordination_decompositions.items():
                still_unmatched_second.update({component: nx_component for component, nx_component in
                                              decomposition.items() if component not in ion_aware_match.second})

            yield from [ion_aware_match + remaining_match for remaining_match in
                        _compute_matches(still_unmatched_first, still_unmatched_second, reaction_center)]


def _compute_sample_isomorphisms(first: mod.Rule.RightGraph, second: mod.Rule.LeftGraph, reaction_center: Set[int]) ->\
        Iterable['Isomorphism']:
    first_graphs = get_rule_component_graphs_with_nx(first)
    second_graphs = get_rule_component_graphs_with_nx(second)

    for globals_match in _compute_matches(first_graphs, second_graphs, reaction_center):
        for match in _complete_match(globals_match, first_graphs, second_graphs, reaction_center):
            if not match.is_complete(vertex.id for vertex in first.vertices):
                continue

            yield match


def _permute_isomorphisms(seed: Iterable['Isomorphism'], reaction_center: Iterable[int]) -> Iterable['Isomorphism']:
    isomorphisms = set(seed)
    for vertex in reaction_center:
        new_isomorphisms = set(isomorphisms)

        for isomorphism in isomorphisms:
            new_isomorphisms.update(isomorphism.permutations(vertex))

        isomorphisms = new_isomorphisms

    return isomorphisms


def _get_automorphisms(graph: mod.Graph, nx_graph: nx.Graph, vertex: int) -> List[Dict[int, int]]:
    group = graph.aut(_label_settings)

    automorphisms = []
    for generator in group.gens:
        automorphism = {}
        for node in nx_graph.nodes:
            target_id = next(target.id for target in nx_graph.nodes if
                             generator[graph.getVertexFromExternalId(node.id)] ==
                             graph.getVertexFromExternalId(target.id))

            if node.id != target_id:
                automorphism[node.id] = target_id

        automorphisms.append(automorphism)

    relevant_ids = {vertex}
    new = 1
    while new > 1:
        new_ids = set(relevant_ids)
        for automorphism in automorphisms:
            if len(relevant_ids.intersection(automorphism)) > 0:
                new_ids.update(automorphism)

        new = len(new_ids) - len(relevant_ids)
        relevant_ids = new_ids

    return [automorphism for automorphism in automorphisms if len(relevant_ids.intersection(automorphism))]


def _apply_automorphisms_to_one_side(seed: 'Isomorphism', automorphisms: List[Dict[int, int]],
                                     apply_call: Callable[['Isomorphism', Dict[int, int]], 'Isomorphism']) ->\
        Iterable['Isomorphism']:
    isomorphisms = {seed}
    new = 1
    while new > 0:
        new_isomorphisms = set(isomorphisms)
        for automorphism in automorphisms:
            new_isomorphisms.update(apply_call(isomorphism, automorphism) for isomorphism in isomorphisms)

        new = len(new_isomorphisms) - len(isomorphisms)
        isomorphisms = new_isomorphisms

    yield from isomorphisms


def _apply_automorphisms(seed: 'Isomorphism', automorphisms: List[Dict[int, int]], first: bool = True) ->\
        Iterable['Isomorphism']:

    if first:
        return _apply_automorphisms_to_one_side(seed, automorphisms, lambda isomorphism, automorphism: isomorphism.
                                                apply(first_automorphism=automorphism))
    else:
        return _apply_automorphisms_to_one_side(seed, automorphisms, lambda isomorphism, automorphism: isomorphism.
                                                apply(second_automorphism=automorphism))


class Isomorphism:
    def __init__(self, reaction_center: Iterable[int]):
        self._reaction_center: Set[int] = set(reaction_center)

        self._first_graphs: Dict[mod.Graph, nx.Graph] = {}
        self._second_graphs: Dict[mod.Graph, nx.Graph] = {}

        self._atom_map: Dict[int, int] = {}
        self._key: Optional[Tuple[int]] = None

    def __add__(self, other: 'Isomorphism') -> 'Isomorphism':
        result = Isomorphism(self._reaction_center.union(other._reaction_center))
        result._first_graphs.update(self._first_graphs)
        result._first_graphs.update(other._first_graphs)

        result._second_graphs.update(self._second_graphs)
        result._second_graphs.update(other._second_graphs)

        result._atom_map.update(self._atom_map)
        result._atom_map.update(other._atom_map)

        return result

    def __eq__(self, other: 'Isomorphism') -> bool:
        return self.key == other.key

    def __getitem__(self, index: int) -> int:
        return self._atom_map[index]

    def __hash__(self) -> int:
        return hash(self.key)

    def __ne__(self, other: 'Isomorphism') -> bool:
        return not self == other

    @property
    def first(self) -> List[mod.Graph]:
        return list(self._first_graphs)

    @property
    def second(self) -> List[mod.Graph]:
        return list(self._second_graphs)

    @property
    def key(self) -> Tuple[int]:
        if self._key is None:
            self._key = tuple(sorted(Counter({source: target for source, target in self._atom_map.items() if
                                              source in self._reaction_center}).elements()))

        return self._key

    @staticmethod
    def from_match(reaction_center: Iterable[int], first_graph: mod.Graph, first_nx_graph: nx.Graph,
                   second_graph: mod.Graph, second_nx_graph: nx.Graph, atom_map: Dict[int, int]) -> 'Isomorphism':
        isomorphism = Isomorphism(reaction_center)

        isomorphism._first_graphs[first_graph] = first_nx_graph
        isomorphism._second_graphs[second_graph] = second_nx_graph

        isomorphism._atom_map.update(atom_map)

        return isomorphism

    def _first_with_vertex(self, vertex: int) -> mod.Graph:
        return next(graph for graph, nx_graph in self._first_graphs.items() if
                    any(node.id == vertex for node in nx_graph.nodes))

    def _second_with_vertex(self, vertex: int) -> mod.Graph:
        return next(graph for graph, nx_graph in self._second_graphs.items() if
                    any(node.id == vertex for node in nx_graph.nodes))

    def compatible(self, other: 'Isomorphism') -> bool:
        return len(set(self.first).intersection(other.first)) == 0 and\
               len(set(self.second).intersection(other.second)) == 0

    def is_complete(self, vertex_ids: Iterable[int]) -> bool:
        return tuple(sorted(self._atom_map)) == tuple(sorted(vertex_ids))

    def apply(self, first_automorphism: Optional[Dict[int, int]] = None,
              second_automorphism: Optional[Dict[int, int]] = None) -> 'Isomorphism':
        if first_automorphism is None:
            first_automorphism = {}
        if second_automorphism is None:
            second_automorphism = {}

        result = Isomorphism(self._reaction_center)

        result._first_graphs = dict(self._first_graphs)
        result._second_graphs = dict(self._second_graphs)

        for source, target in self._atom_map.items():
            if source in first_automorphism:
                source = first_automorphism[source]
            if target in second_automorphism:
                target = second_automorphism[target]

            result._atom_map[source] = target

        return result

    def permutations(self, vertex: int) -> Iterable['Isomorphism']:
        first_graph = self._first_with_vertex(vertex)
        second_graph = self._second_with_vertex(vertex)

        first_automorphisms = _get_automorphisms(first_graph, self._first_graphs[first_graph], vertex)
        second_automorphisms = _get_automorphisms(second_graph, self._second_graphs[second_graph], vertex)

        for isomorphism in _apply_automorphisms(self, first_automorphisms):
            yield from _apply_automorphisms(isomorphism, second_automorphisms, False)


class IsomorphismCacheEntry:
    def __init__(self, value: Iterable[Isomorphism], parent: Optional['IsomorphismCacheEntry']):
        self._value: List[Isomorphism] = list(value)
        self._parent: Optional[IsomorphismCacheEntry] = parent

        self._cache: Dict[int, IsomorphismCacheEntry] = {}

    def get_isomorphisms(self, reaction_center: Tuple[int], verbosity: int = 0) -> Iterable[Isomorphism]:
        yield from self._value

        if len(reaction_center) == 0:
            return

        minimal_vertex = reaction_center[0]

        if minimal_vertex not in self._cache:
            if verbosity > 4:
                print(f"\t\t\t#\tCreating new cache entry for reaction center {reaction_center}.")

            known_isomorphisms = list(self.get_parent_isomorphisms())
            new_isomorphisms = []
            for permutation in _permute_isomorphisms(known_isomorphisms, [minimal_vertex]):
                if permutation not in known_isomorphisms:
                    new_isomorphisms.append(permutation)

            if verbosity > 5:
                print(f"\t\t\t#\tExpanding on {len(known_isomorphisms)} known isomorphisms.")
                print(f"\t\t\t#\tFound {len(new_isomorphisms)} new isomorphisms.")

            self._cache[minimal_vertex] = IsomorphismCacheEntry(new_isomorphisms, self)

        yield from self._cache[minimal_vertex].get_isomorphisms(reaction_center[1:], verbosity)

    def get_parent_isomorphisms(self) -> Iterable[Isomorphism]:
        yield from self._value

        if self._parent is not None:
            yield from self._parent.get_parent_isomorphisms()


class IsomorphismCache:
    def __init__(self):
        self._cache: Dict[mod.Rule, Dict[mod.Rule, IsomorphismCacheEntry]] = {}

    def get_isomorphisms(self, first: mod.Rule, second: mod.Rule, reaction_center: Tuple[int], verbosity: int = 0) ->\
            Iterable[Isomorphism]:
        if first not in self._cache:
            self._cache[first] = {}

        if second not in self._cache[first]:
            self._cache[first][second] = IsomorphismCacheEntry(_compute_sample_isomorphisms(first.right, second.left,
                                                                                            set(reaction_center)),
                                                               None)

        return self._cache[first][second].get_isomorphisms(reaction_center, verbosity)
