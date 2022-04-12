import mod


from typing import Dict, Optional, List, Set, Tuple, Union


class RuleEdgeTuple(Tuple[mod.Rule.Vertex, mod.Rule.Vertex]):
    def __new__(cls, edge: mod.Rule.Edge):
        return super().__new__(cls, sorted([edge.source, edge.target], key=lambda v: v.id))


class RuleGraphContainer:
    def __init__(self):
        self._vertices: List[Union[mod.Rule.LeftGraph.Vertex, mod.Rule.RightGraph.Vertex]] = []
        self._edges: List[Union[mod.Rule.LeftGraph.Edge, mod.Rule.RightGraph.Edge]] = []

    @property
    def vertices(self) -> List[Union[mod.Rule.LeftGraph.Vertex, mod.Rule.RightGraph.Vertex]]:
        return list(self._vertices)

    @property
    def edges(self) -> List[Union[mod.Rule.LeftGraph.Edge, mod.Rule.RightGraph.Edge]]:
        return list(self._edges)

    def add_vertex(self, vertex: Union[mod.Rule.LeftGraph.Vertex, mod.Rule.RightGraph.Vertex]):
        self._vertices.append(vertex)

    def add_edge(self, edge: Union[mod.Rule.LeftGraph.Edge, mod.Rule.RightGraph.Edge]):
        self._edges.append(edge)


class FilteredRule:
    def __init__(self, rule: mod.Rule):
        self._rule: mod.Rule = rule

        self._used_vertices: Set[mod.Rule.Vertex] = set()
        self._used_edges: Dict[RuleEdgeTuple, mod.Rule.Edge] = {}

        self._protected_vertices: Set[mod.Rule.Vertex] = set()

        self._relabels: Dict[Union[mod.Rule.LeftGraph.Vertex, mod.Rule.RightGraph.Vertex], str] = {}

    @property
    def rule(self) -> mod.Rule:
        return self._rule

    @property
    def used_vertices(self) -> Set[mod.Rule.Vertex]:
        return set(self._used_vertices)

    @property
    def used_edges(self) -> Dict[RuleEdgeTuple, mod.Rule.Edge]:
        return dict(self._used_edges)

    def uses_vertex(self, vertex: mod.Rule.Vertex):
        return vertex in self.used_vertices

    def uses_edge(self, edge: mod.Rule.Edge):
        return RuleEdgeTuple(edge) in self.used_edges

    def include_all(self):
        for vertex in self.rule.vertices:
            self.include_vertex(vertex)

        for edge in self.rule.edges:
            self.include_edge(edge)

    def include_vertex(self, vertex: mod.Rule.Vertex, protected: bool = False, left_label: str = None,
                       right_label: str = None):
        if vertex.rule != self.rule:
            return

        self._used_vertices.add(vertex)

        if protected:
            self._protected_vertices.add(vertex)

        if left_label is not None:
            assert(not vertex.left.isNull())
            assert(vertex.left not in self._relabels)
            self._relabels[vertex.left] = left_label

        if right_label is not None:
            assert(not vertex.right.isNull())
            assert(vertex.right not in self._relabels)
            self._relabels[vertex.right] = right_label

    def include_edge(self, edge: mod.Rule.Edge, protected: bool = False):
        if edge.rule != self.rule:
            return

        self._used_edges[RuleEdgeTuple(edge)] = edge

        self.include_vertex(edge.source, protected)
        self.include_vertex(edge.target, protected)

    def include_reaction_center(self, only_include_changed_bonds: bool = False, protected: bool = False):
        for vertex in self.rule.vertices:
            if vertex.left.stringLabel != vertex.right.stringLabel:
                self.include_vertex(vertex, protected)

        for edge in self.rule.edges:
            if edge.left.isNull() or edge.right.isNull() or edge.left.stringLabel != edge.right.stringLabel:
                self.include_edge(edge, protected)

        if only_include_changed_bonds:
            return

        for edge in self.rule.edges:
            if edge.source in self.used_vertices and edge.target in self.used_vertices:
                self.include_edge(edge, protected)

    def remove_vertex(self, vertex: mod.Rule.Vertex):
        if not self.uses_vertex(vertex) or vertex in self._protected_vertices:
            return

        self._used_vertices.remove(vertex)
        for edge in vertex.incidentEdges:
            self.remove_edge(edge)

    def remove_edge(self, edge: mod.Rule.Edge):
        if not self.uses_edge(edge) or\
                (edge.source in self._protected_vertices and edge.target in self._protected_vertices):
            return

        del self._used_edges[RuleEdgeTuple(edge)]

    def to_gml(self, name: Optional[str] = None, unlabelled_vertices: Optional[Set[mod.Rule.LeftGraph.Vertex]] = None):
        if name is None:
            name = self.rule.name

        left = RuleGraphContainer()
        context = RuleGraphContainer()
        right = RuleGraphContainer()

        for vertex in self.used_vertices:
            if vertex.left.isNull():
                right.add_vertex(vertex.right)
            elif vertex.right.isNull():
                left.add_vertex(vertex.left)
            elif vertex.left.stringLabel != vertex.right.stringLabel:
                left.add_vertex(vertex.left)
                right.add_vertex(vertex.right)
            else:
                context.add_vertex(vertex.left)

        for edge in self.used_edges.values():
            if edge.left.isNull():
                right.add_edge(edge.right)
            elif edge.right.isNull():
                left.add_edge(edge.left)
            elif edge.left.stringLabel != edge.right.stringLabel:
                left.add_edge(edge.left)
                right.add_edge(edge.right)
            else:
                context.add_edge(edge.left)

        rule_sections = (('left', left), ('context', context), ('right', right))

        out = [f'rule [ ruleID "{name}"']

        for name, container in rule_sections:
            out.append(f'\t{name} [')

            for vertex in container.vertices:
                label: str = vertex.stringLabel
                if name == "context" and unlabelled_vertices is not None and vertex in unlabelled_vertices:
                    label = "*"
                if vertex in self._relabels:
                    label = self._relabels[vertex]

                out.append(f'\t\tnode [ id {vertex.id} label "{label}" ]')

            for edge in container.edges:
                out.append(f'\t\tedge [ source {edge.source.id} target {edge.target.id} label "{edge.stringLabel}" ]')

            out.append('\t]')

        out.append(']')

        return '\n'.join(out)

    def to_mod_rule(self) -> mod.Rule:
        return mod.ruleGMLString(self.to_gml())
