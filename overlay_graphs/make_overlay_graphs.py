import io
import json
import mod
import networkx as nx


from collections import Counter
from overlay_graphs.canonicalisation import GraphCanonicaliser
from overlay_graphs.draw import print_overlay_graph
from overlay_graphs.mechanism import Mechanism, Step
from overlay_graphs.networkx_converter import get_component_graphs, graph_to_nx_graph
from overlay_graphs.overlay_graph import OverlayGraph
from overlay_graphs.reaction_sequence_tracking import IsomorphismCache
from overlay_graphs.rule_builder import EdgeTuple
from typing import Dict, Iterable, List, Optional, Tuple, Set, Union


_edge_valence_symbols = ["?", "-", "=", "#"]


def _parse_atom_maps(mechanisms: List[Mechanism], atom_map_file: str = "manual_atom_maps.json") ->\
        Dict[Mechanism, List[Dict[int, int]]]:
    with open(atom_map_file, "r") as file:
        atom_map_json = json.load(file)

    atom_maps = {}
    for entry in atom_map_json:
        mechanism = next((m for m in mechanisms if m.entry == entry["entry"] and m.number == entry["mechanism"]), None)
        if mechanism is None:
            continue

        atom_maps[mechanism] = [{}] * len(mechanism)

        for map in entry["atom_maps"]:
            atom_map = {}
            for atom in map["map"]:
                atom_map[atom["atom"]] = atom["original"]

            atom_maps[mechanism][map["step"] - 1] = atom_map

    return atom_maps


def _extend_overlay_graph(canonicaliser: GraphCanonicaliser, isomorphism_cache: 'IsomorphismCache',
                          marking: 'OverlayMarking', atom_map: Dict[int, int], last_rule: mod.Rule,
                          mechanism: List[Step], atom_maps: List[Dict[int, int]], verbosity: int = 0) ->\
        Iterable[OverlayGraph]:
    if len(mechanism) == 0:
        yield OverlayGraph(marking.host_graph, marking.to_dictionary())
        return

    if verbosity > 1:
        print(f"\t#\tExtending overlay graph after rule {last_rule}. {len(mechanism)} steps to go.")

    canonical_overlay_graphs = set()
    action_atoms: Set[int] = {atom_map[vertex] for vertex in marking.action}
    reaction_center = tuple(sorted(action_atoms))

    for index, isomorphism in enumerate(isomorphism_cache.get_isomorphisms(last_rule, mechanism[0].rule,
                                                                           reaction_center, verbosity)):
        if verbosity > 2:
            print(f"\t\t#\tFound isomorphism number {index} for the OG extension after {last_rule}.")

        new_atom_map = {original_id: isomorphism[last_id] for original_id, last_id in atom_map.items()}

        violating_atom = next((atom for atom, original in atom_maps[0].items() if new_atom_map[original] != atom), None)
        if violating_atom is not None:
            if verbosity > 3:
                print(f"\t\t#\tIntermediary overlay graph after {last_rule} violates prescribed atom map on atom "
                      f"'{violating_atom}'. Discarding...")
            continue

        new_marking = marking.update_from_rule(mechanism[0].rule, {new_id: original_id for original_id, new_id in
                                                                   new_atom_map.items()})

        intermediary_og = OverlayGraph(new_marking.host_graph, new_marking.to_dictionary())
        canonical_og = canonicaliser.canonicalise_nx_graph(intermediary_og.to_labelled_graph("L_+_-", "L_+_-"))
        if canonical_og in canonical_overlay_graphs:
            if verbosity > 3:
                print(f"\t\t#\tIntermediary overlay graph after {last_rule} isomorphic to previous OG. Discarding...")
            continue

        canonical_overlay_graphs.add(canonical_og)

        yield from _extend_overlay_graph(canonicaliser, isomorphism_cache, new_marking, new_atom_map, mechanism[0].rule,
                                         mechanism[1:], atom_maps[1:], verbosity)


def compute_overlay_graphs(canonicaliser: GraphCanonicaliser, isomorphism_cache: 'IsomorphismCache',
                            mechanism: Mechanism, atom_maps: List[Dict[int, int]], verbosity: int = 0) ->\
        Iterable[OverlayGraph]:
    host_graph = graph_to_nx_graph(mechanism[0].rule.left, use_indices=True)

    atom_map = {node: node for node in host_graph.nodes}

    marking = OverlayMarking(host_graph).update_from_rule(mechanism[0].rule, atom_map)

    if verbosity > 0:
        print(f"\t#\tCreated initial marking on {len(atom_map)} vertices.")

    yield from _extend_overlay_graph(canonicaliser, isomorphism_cache, marking, atom_map, mechanism[0].rule,
                                     list(mechanism)[1:], atom_maps[1:], verbosity)


def overlay_graphs_for_mechanisms(mechanisms: List[Mechanism], output_name: str = "overlay_graphs",
                                  known_atom_maps: Optional[Dict[Mechanism, List[Dict[int, int]]]] = None,
                                  verbosity: int = 0):
    if known_atom_maps is None:
        known_atom_maps = {}

    canonicaliser = GraphCanonicaliser()

    with open(f"{output_name}.json", "w") as file:
        file.write("[\n]\n")

    for mechanism in mechanisms:
        if verbosity >= 1:
            print(f"#\n#\tComputing overlay graphs for mechanism '{mechanism}'.\n#")

        if mechanism.entry == -1 or len(mechanism) == 0 or any(step.rule is None for step in mechanism):
            continue

        isomorphism_cache = IsomorphismCache()

        canonical_overlay_graphs = {}

        if mechanism not in known_atom_maps:
            known_atom_maps[mechanism] = [{}] * len(mechanism)

        for index, overlay_graph in enumerate(compute_overlay_graphs(canonicaliser, isomorphism_cache, mechanism,
                                                                     list(known_atom_maps[mechanism]), verbosity)):
            canonical_overlay_graph = canonicaliser.canonicalise_nx_graph(overlay_graph.to_labelled_graph("L_+_-", "L_+_-"))
            if canonical_overlay_graph in canonical_overlay_graphs:
                continue

            canonical_overlay_graphs[canonical_overlay_graph] = overlay_graph

        mod.postChapter(f"{mechanism.entry}_{mechanism.number} [{len(canonical_overlay_graphs)}]")

        for index, overlay_graph in enumerate(canonical_overlay_graphs.values()):
            if index > 10:
                break

            mod.postSection(f"OG {index}")
            print_overlay_graph(overlay_graph)

        with open(f"{output_name}.json", "r+") as file:
            file.seek(0, io.SEEK_END)
            pos = file.tell()
            file.seek(pos - 2)
            if pos > 4:
                file.write(",\n")
            file.write(json.dumps({"mechanism": mechanism.serialise(),
                                   "overlay_graphs": [graph.serialise() for graph in canonical_overlay_graphs.values()]}))
            file.write("]\n")

        if verbosity >= 1:
            print(f"#\tFound {len(canonical_overlay_graphs)} unique overlay graphs.\n")

    with open(f"{output_name}.json", "r") as file:
        data = json.load(file)

    with open(f"{output_name}.json", "w") as file:
        json.dump(data, file, indent=2)


class OverlayMarking:
    def __init__(self, host_graph: nx.Graph):
        self._host_graph: nx.Graph = host_graph.copy()
        self._elements: Set[Union[int, EdgeTuple]] = set(EdgeTuple((source, target)) for source, target
                                                         in host_graph.edges).union(host_graph.nodes)

        self._electrons_donated: Counter[Union[int, EdgeTuple]] = Counter({element: 0 for element in self._elements})
        self._electrons_received: Counter[Union[int, EdgeTuple]] = Counter({element: 0 for element in self._elements})

    @property
    def host_graph(self) -> nx.Graph:
        clean_host = nx.Graph()

        for graph, nx_graph in get_component_graphs(self._host_graph).items():
            if all(node not in self.action for node in nx_graph.nodes) and\
                all(self._electrons_donated[EdgeTuple((source, target))] == 0 and
                    self._electrons_received[EdgeTuple((source, target))] == 0 for source, target in
                    nx_graph.edges):
                continue

            for node, data in nx_graph.nodes(data=True):
                clean_host.add_node(node, label=data["label"])

            for source, target, data in nx_graph.edges(data=True):
                label = data["label"] if data["label"] != ":" else ">"
                clean_host.add_edge(source, target, label=label)

        return clean_host

    @property
    def action(self) -> Iterable[int]:
        for element in self._elements:
            if self._electrons_donated[element] == 0 and self._electrons_received[element] == 0:
                continue

            if isinstance(element, int):
                yield element
            else:
                yield element[0]
                yield element[1]

    def _ensure_element_exists(self, element: Union[int, EdgeTuple]):
        if element in self._elements:
            return

        self._elements.add(element)
        if not isinstance(element, int):
            self._host_graph.add_edge(element[0], element[1], label="?")

        self._electrons_donated[element] = 0
        self._electrons_received[element] = 0

    def copy(self) -> 'OverlayMarking':
        copy = OverlayMarking(self._host_graph)

        copy._elements = set(self._elements)
        copy._electrons_donated = dict(self._electrons_donated)
        copy._electrons_received = dict(self._electrons_received)

        return copy

    def add_electron_donated(self, element: Union[int, Tuple[int, int]]):
        if not isinstance(element, int):
            element = EdgeTuple(element)

        self._ensure_element_exists(element)

        self._electrons_donated[element] += 1

    def add_electron_received(self, element: Union[int, Tuple[int, int]]):
        if not isinstance(element, int):
            element = EdgeTuple(element)

        self._ensure_element_exists(element)

        self._electrons_received[element] += 1

    def update_from_rule(self, rule: mod.Rule, atom_map: Dict[int, int]) -> 'OverlayMarking':
        result = self.copy()

        for vertex in rule.vertices:
            if vertex.left.stringLabel == vertex.right.stringLabel:
                continue

            if int(vertex.left.charge) > int(vertex.right.charge):
                result.add_electron_received(atom_map[vertex.id])
            elif int(vertex.left.charge) < int(vertex.right.charge):
                result.add_electron_donated(atom_map[vertex.id])

        for edge in rule.edges:
            if edge.left.isNull():
                if edge.right.stringLabel == ":":
                    continue

                result.add_electron_received((atom_map[edge.source.id], atom_map[edge.target.id]))
            elif edge.right.isNull():
                if edge.left.stringLabel == ":":
                    continue

                result.add_electron_donated((atom_map[edge.source.id], atom_map[edge.target.id]))
            elif edge.left.stringLabel == ":" or edge.right.stringLabel == ":":
                continue
            elif _edge_valence_symbols.index(edge.left.stringLabel) <\
                    _edge_valence_symbols.index(edge.right.stringLabel):
                result.add_electron_received((atom_map[edge.source.id], atom_map[edge.target.id]))
            elif _edge_valence_symbols.index(edge.left.stringLabel) >\
                    _edge_valence_symbols.index(edge.right.stringLabel):
                result.add_electron_donated((atom_map[edge.source.id], atom_map[edge.target.id]))

        return result

    def to_dictionary(self) -> Dict[Union[int, EdgeTuple], Tuple[int, int]]:
        return {element: (self._electrons_received[element], self._electrons_donated[element]) for element in
                self._elements if self._electrons_received[element] > 0 or self._electrons_donated[element] > 0}
