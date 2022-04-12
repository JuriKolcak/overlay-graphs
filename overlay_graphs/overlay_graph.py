import mod
import networkx as nx
import re


from overlay_graphs.networkx_converter import get_component_graphs
from overlay_graphs.rule_builder import EdgeTuple
from typing import Any, Dict, Iterable, Optional, Tuple, Union


_atom_label_pattern: re.Pattern = re.compile(r"^[\s]*([a-zA-Z]+)([0-9]*)([+\-]*)[\s]*$")


_edge_valence_symbols = ["?", "-", "=", "#", ":", ">"]


def _parse_atom_label(label: str) -> (str, int):
    match = _atom_label_pattern.match(label)
    if match is None:
        return label, 0

    atom_type = match.group(1)

    if match.group(3) != "":
        if match.group(2) != "":
            charge = int(match.group(2))
        else:
            charge = 1

        if match.group(3) == "-":
            charge = -charge
    else:
        charge = 0

    return atom_type, charge


def _compose_atom_label(atom_type: str, charge: int) -> str:
    if charge != 0:
        if charge < 0:
            charge_label = "-"
        else:
            charge_label = "+"

        if abs(charge) > 1:
            charge_label = f"{abs(charge)}{charge_label}"
    else:
        charge_label = ""

    return f"{atom_type}{charge_label}"


class OverlayGraph:
    def __init__(self, graph: nx.Graph, marking: Dict[Union[int, Tuple[int, int]], Tuple[int, int]]):
        self._host_graph: nx.Graph = graph
        self._marking: Dict[Union[int, EdgeTuple], Tuple[int, int]] =\
            {self._canonicalise_element(key): value for key, value in marking.items()}

        self._product_graph: Optional[nx.Graph] = None

    def __copy__(self) -> 'OverlayGraph':
        return OverlayGraph(self.host_graph, dict(self._marking))

    @property
    def host_graph(self) -> nx.Graph:
        return self._host_graph

    @property
    def product_graph(self) -> nx.Graph:
        if self._product_graph is None:
            self._product_graph = nx.Graph()

            for node, data in self.host_graph.nodes(data=True):
                atom_type, charge = _parse_atom_label(data["label"])

                electron_balance = self.electrons_received(node) - self.electrons_donated(node)

                self._product_graph.add_node(node, label=_compose_atom_label(atom_type, charge - electron_balance))

            for source, target, data in self.host_graph.edges(data=True):
                edge_tuple = EdgeTuple((source, target))

                valence = _edge_valence_symbols.index(data["label"])
                electron_balance = self.electrons_received(edge_tuple) - self.electrons_donated(edge_tuple)

                self._product_graph.add_edge(source, target, label=_edge_valence_symbols[valence + electron_balance])

        return self._product_graph

    @property
    def action(self) -> Iterable[Union[int, EdgeTuple]]:
        return self._marking

    @staticmethod
    def _canonicalise_element(element: Union[int, Tuple[int, int]]) -> Union[int, EdgeTuple]:
        if isinstance(element, int):
            return element

        return EdgeTuple(element)

    @staticmethod
    def deserialise(og_json: Dict[str, Any]) -> 'OverlayGraph':
        host_graph = nx.Graph()
        marking = {}

        for node in og_json["nodes"]:
            donated = node["electrons_donated"]
            received = node["electrons_received"]

            if donated > 0 or received > 0:
                marking[node["id"]] = (received, donated)

            host_graph.add_node(node["id"], label=node["label"], electrons_donated=donated, electrons_received=received)

        for edge in og_json["edges"]:
            donated = edge["electrons_donated"]
            received = edge["electrons_received"]

            if donated > 0 or received > 0:
                marking[(edge["src"], edge["tar"])] = (received, donated)

            host_graph.add_edge(edge["src"], edge["tar"], label=edge["label"], electrons_donated=donated,
                                electrons_received=received)

        return OverlayGraph(host_graph, marking)

    def _og_marking(self, element: Union[int, EdgeTuple], add: bool) -> Tuple[int, int]:
        if add:
            if element not in self._marking:
                self._marking[element] = (0, 0)

            return self._marking[element]

        return self.electrons_received(element), self.electrons_donated(element)

    def _apply_label_pattern(self, label_pattern: str, original_label: str, element: Union[int, EdgeTuple],
                             include_blue: bool = True) -> str:
        electrons_donated, electrons_received = self._og_marking(element, False)
        if not include_blue and electrons_donated == electrons_received:
            electrons_donated, electrons_received = 0, 0

        label = re.sub(r"L", original_label, label_pattern)
        label = re.sub(r"-", str(electrons_donated), label)
        label = re.sub(r"\+", str(electrons_received), label)
        return re.sub(r"b", str(electrons_received - electrons_donated), label)

    def marking_key(self) -> frozenset[Tuple[Union[int, EdgeTuple], Tuple[int, int]]]:
        return frozenset([(k, v) for k, v in self._marking.items()])

    def electrons_received(self, element: Union[int, Tuple[int, int]]) -> int:
        element = self._canonicalise_element(element)
        return self._marking[element][0] if element in self._marking else 0

    def electrons_donated(self, element: Union[int, Tuple[int, int]]) -> int:
        element = self._canonicalise_element(element)
        return self._marking[element][1] if element in self._marking else 0

    def add_received(self, element: Union[int, Tuple[int, int]]):
        element = self._canonicalise_element(element)
        received, donated = self._og_marking(element, True)

        self._marking[element] = (received + 1, donated)

    def add_donated(self, element: Union[int, Tuple[int, int]]):
        element = self._canonicalise_element(element)
        received, donated = self._og_marking(element, True)

        self._marking[element] = (received, donated + 1)

    def remove_received(self, element: Union[int, Tuple[int, int]]):
        element = self._canonicalise_element(element)
        received, donated = self._og_marking(element, False)

        assert (received > 0)

        self._marking[element] = (received - 1, donated)

    def remove_donated(self, element: Union[int, Tuple[int, int]]):
        element = self._canonicalise_element(element)
        received, donated = self._og_marking(element, False)

        assert (donated > 0)

        self._marking[element] = (received, donated - 1)

    def serialise(self) -> Dict[str, Any]:
        return {
            "nodes": [{
                "id": node,
                "label": data["label"],
                "electrons_donated": self.electrons_donated(node),
                "electrons_received": self.electrons_received(node)
            } for node, data in self.host_graph.nodes(data=True)],
            "edges": [{
                "src": source,
                "tar": target,
                "label": data["label"],
                "electrons_donated": self.electrons_donated((source, target)),
                "electrons_received": self.electrons_received((source, target))
            } for source, target, data in self.host_graph.edges(data=True)]
        }

    def catalytic_vertices(self):
        label_settings = mod.LabelSettings(mod.LabelType.String, mod.LabelRelation.Isomorphism)

        left = nx.subgraph_view(self.host_graph,
                                filter_edge=lambda src, tar: self.host_graph.edges[src, tar]["label"] != "?")
        right = nx.subgraph_view(self.product_graph,
                                 filter_edge=lambda src, tar: self.product_graph.edges[src, tar]["label"] != "?")
        left_graphs = get_component_graphs(left)
        right_graphs = get_component_graphs(right)

        left_enzyme_components = []
        right_enzyme_components = []
        for left_graph, left_nx in left_graphs.items():
            isomorphic_right_graphs = {graph: [node for node in left_nx.nodes if node in right_graphs[graph].nodes]
                                       for graph in right_graphs
                                       if graph.isomorphism(left_graph, labelSettings=label_settings) != 0 and
                                       graph not in right_enzyme_components}
            print(left_nx.nodes)
            print(len(isomorphic_right_graphs))
            print(right_enzyme_components)

            for graph in list(isomorphic_right_graphs):
                if len(isomorphic_right_graphs[graph]) == 0:
                    del isomorphic_right_graphs[graph]

            if len(isomorphic_right_graphs) == 0:
                continue

            for right_graph in sorted(isomorphic_right_graphs, key=lambda g: len(isomorphic_right_graphs[g]),
                                      reverse=True):
                if len(isomorphic_right_graphs[right_graph]) > left_graph.numVertices / 2 or \
                        any(left_nx.nodes[node]["label"] != "H" for node in isomorphic_right_graphs[right_graph]):

                    left_enzyme_components.append(left_graph)
                    right_enzyme_components.append(right_graph)
                    break
        left_enzyme_vertices = set()
        for component in left_enzyme_components:
            left_enzyme_vertices.update(node for node in left_graphs[component].nodes)
        print("LEFT: ", left_enzyme_vertices)

        right_enzyme_vertices = set()
        for component in right_enzyme_components:
            right_enzyme_vertices.update(node for node in right_graphs[component].nodes)
        print("RIGHT: ", left_enzyme_vertices)

        return left_enzyme_vertices.intersection(right_enzyme_vertices)

    def reindex(self) -> 'OverlayGraph':
        new_graph = nx.Graph()

        marking = {}
        node2node = {}

        for index, (node, data) in enumerate(self.host_graph.nodes(data=True)):
            new_graph.add_node(index, label=data["label"])
            node2node[node] = index
            og_marking = self._og_marking(node, False)
            marking[index] = og_marking

        for source, target, data in self.host_graph.edges(data=True):
            new_edge = EdgeTuple((node2node[source], node2node[target]))
            new_graph.add_edge(new_edge[0], new_edge[1], label=data["label"])
            og_marking = self._og_marking(EdgeTuple((source, target)), False)
            marking[new_edge] = og_marking

        return OverlayGraph(new_graph, marking)

    def to_labelled_graph(self, vertex_label_pattern: str, edge_label_pattern: str) -> nx.Graph:
        labelled_graph = self.host_graph.copy()

        for node, data in self.host_graph.nodes(data=True):
            labelled_graph.nodes[node]["label"] = self._apply_label_pattern(vertex_label_pattern, data["label"], node)

        for source, target, data in self.host_graph.edges(data=True):
            labelled_graph.edges[source, target]["label"] = self._apply_label_pattern(edge_label_pattern, data["label"],
                                                                                      EdgeTuple((source, target)))

        return labelled_graph

    def to_gml(self, vertex_label_pattern: str, edge_label_pattern: str, include_blue: bool = True) -> str:
        vertices = set()
        edges = set()

        for edge in self.host_graph.edges:
            edge_tuple = EdgeTuple(edge)
            if edge_tuple in self._marking and \
                    (include_blue or self.electrons_donated(edge) != self.electrons_received(edge)):
                vertices.update(edge_tuple)
                edges.add(edge_tuple)

        for node in self.host_graph.nodes:
            if node in self._marking and \
                    (include_blue or self.electrons_donated(node) != self.electrons_received(node)):
                vertices.add(node)

        out = ["graph", "["]
        for index in vertices:
            label = self._apply_label_pattern(vertex_label_pattern, self.host_graph.nodes[index]["label"], index,
                                              include_blue)
            out.append(f"\tnode [ id {index} label \"{label}\" ]")

        for indices in edges:
            label = self._apply_label_pattern(edge_label_pattern, self.host_graph.edges[indices]["label"], indices,
                                              include_blue)

            out.append(f"\tedge [ source {indices[0]} target {indices[1]} label \"{label}\" ]")
        out.append("]")

        return "\n".join(out)

    def to_mod_graph(self, vertex_label_pattern: str, edge_label_pattern: str, include_blue: bool = True) ->\
            Optional[mod.Graph]:
        try:
            return mod.graphGMLString(self.to_gml(vertex_label_pattern, edge_label_pattern, include_blue))
        except mod.InputError:
            return None
