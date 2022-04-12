import json
import mod
import networkx as nx
import re


from overlay_graphs.networkx_converter import get_component_graphs, nx_graph_to_gml
from overlay_graphs.rule_builder import EdgeTuple, RuleBuilder
from overlay_graphs.overlay_graph import OverlayGraph
from overlay_graphs.mcsadb import MCSAEntry
from typing import Dict, Iterable, Tuple


_non_and_dative = re.compile(r"^[?>]$")


def _filter_graph_edges(graph: nx.Graph, label_patter: re.Pattern) -> nx.Graph:
    return nx.subgraph_view(graph, filter_edge=lambda source, target: re.
                            match(label_patter, graph.edges[source, target]["label"]) is None)


def _remove_ghost_edges(graph: nx.Graph) -> nx.Graph:
    new_graph = graph.copy()

    for source, target, data in graph.edges(data=True):
        if data["label"] == "?":
            new_graph.remove_edge(source, target)

    return new_graph


def _make_substrate_rule(overlay_graph: OverlayGraph, rule_id: str, label_settings: mod.LabelSettings, disable_catalytic: bool) ->\
        Tuple[str, Dict[mod.Graph, nx.Graph], Dict[mod.Graph, nx.Graph]]:
    left = _filter_graph_edges(overlay_graph.host_graph, _non_and_dative)
    right = _filter_graph_edges(overlay_graph.product_graph, _non_and_dative)

    left_graphs = get_component_graphs(left)
    right_graphs = get_component_graphs(right)

    left_enzyme_components = []
    right_enzyme_components = []
    for left_graph, left_nx in left_graphs.items():
        isomorphic_right_graphs = {graph: [node for node in left_nx.nodes if node in right_graphs[graph].nodes]
                                   for graph in right_graphs
                                   if graph.isomorphism(left_graph, labelSettings=label_settings) != 0 and
                                   graph not in right_enzyme_components}

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

    right_enzyme_vertices = set()
    for component in right_enzyme_components:
        right_enzyme_vertices.update(node for node in right_graphs[component].nodes)

    rule_builder = RuleBuilder(rule_id)

    for index, data in left.nodes(data=True):
        if index not in overlay_graph.action:
            continue

        if disable_catalytic or index not in left_enzyme_vertices:
            rule_builder.add_left_vertex(index, data["label"])
        if disable_catalytic or index not in right_enzyme_vertices:
            rule_builder.add_right_vertex(index, right.nodes[index]["label"])

    for source, target, data in overlay_graph.host_graph.edges(data=True):
        edge_tuple = EdgeTuple((source, target))

        if edge_tuple not in overlay_graph.action:
            continue

        for index in edge_tuple:
            if index not in rule_builder.vertices:
                if disable_catalytic or index not in left_enzyme_vertices:
                    rule_builder.add_left_vertex(index, left.nodes[index]["label"])
                if disable_catalytic or index not in right_enzyme_vertices:
                    rule_builder.add_right_vertex(index, right.nodes[index]["label"])

        if edge_tuple in left.edges and (
                disable_catalytic or all(index not in left_enzyme_vertices for index in edge_tuple)):
            rule_builder.add_left_edge(source, target, data["label"])

        if edge_tuple in right.edges and (
                disable_catalytic or all(index not in right_enzyme_vertices for index in edge_tuple)):
            rule_builder.add_right_edge(source, target, right.edges[edge_tuple]["label"])

    for source, target, data in overlay_graph.host_graph.edges(data=True):
        edge_tuple = EdgeTuple((source, target))

        if rule_builder.has_vertex(source) and rule_builder.has_vertex(target) and \
                not rule_builder.has_edge(source, target):
            if edge_tuple in left.edges:
                rule_builder.add_left_edge(source, target, data["label"])
            if edge_tuple in right.edges:
                rule_builder.add_right_edge(source, target, right.edges[edge_tuple]["label"])

            continue

    return (rule_builder.to_gml(), {comp: left_graphs[comp] for comp in left_enzyme_components},
            {comp: right_graphs[comp] for comp in right_enzyme_components})


def _make_substrate_rules(entry: MCSAEntry, rule_id: str, label_settings: mod.LabelSettings, disable_catalytic: bool) ->\
        Iterable[Tuple[str, Dict[mod.Graph, nx.Graph], Dict[mod.Graph, nx.Graph]]]:
    for index, overlay_graph in enumerate(entry.overlay_graphs):
        yield _make_substrate_rule(overlay_graph, f"{rule_id}_{index}", label_settings, disable_catalytic)


def substrate_rules_for_entries(entries: Iterable[MCSAEntry]):
    label_settings = mod.LabelSettings(mod.LabelType.String, mod.LabelRelation.Isomorphism)

    rules = {}
    disable_catalytic = False

    for entry in entries:
        rule_id = f"OG_{entry.entry}_{entry.mechanism}"

        for index, substrate_rule_tuple in enumerate(_make_substrate_rules(entry, rule_id, label_settings, disable_catalytic)):
            rules[f"{rule_id}_{index}"] = substrate_rule_tuple

    with open("substrate_rules.json", "w") as file:
        json.dump({
            "rules": [{"name": rule_id, "gml": gml,
                       "left_catalysts": [{"name": catalyst.name, "gml": nx_graph_to_gml(nx)}
                                          for catalyst, nx in left_catalysts.items()],
                       "right_catalysts": [{"name": catalyst.name, "gml": nx_graph_to_gml(nx)}
                                           for catalyst, nx in right_catalysts.items()]}
                      for rule_id, (gml, left_catalysts, right_catalysts) in rules.items()]
        }, file, indent=2)

    printer = mod.GraphPrinter()
    mod.config.stereo.silenceDeductionWarnings = True
    for rule_id, (gml, left, right) in rules.items():
        mod.postSection(rule_id)
        mod.ruleGMLString(gml).print(printer)


