import json


from overlay_graphs.mechanism import ECNumber, Mechanism
from overlay_graphs.overlay_graph import OverlayGraph
from typing import Any, Dict, Iterable, List, Optional


class MCSAEntry:
    def __init__(self, mechanism: Mechanism, overlay_graphs: Iterable[OverlayGraph]):
        self._mechanism: Mechanism = mechanism

        self._overlay_graphs = list(overlay_graphs)

    def __eq__(self, other: 'MCSAEntry') -> bool:
        return self._mechanism == other._mechanism

    def __hash__(self) -> int:
        return hash(self._mechanism)

    def __ne__(self, other: 'MCSAEntry') -> bool:
        return not self == other

    @property
    def entry(self) -> int:
        return self._mechanism.entry

    @property
    def mechanism(self) -> int:
        return self._mechanism.number

    def _get_ec(self) -> Optional[ECNumber]:
        return self._mechanism.ec

    def _set_ec(self, value: ECNumber):
        self._mechanism.ec = value

    ec = property(fget=_get_ec, fset=_set_ec)

    @property
    def overlay_graphs(self) -> List[OverlayGraph]:
        return list(self._overlay_graphs)

    @staticmethod
    def deserialize(entry_json: Dict[str, Any]):
        return MCSAEntry(Mechanism.deserialise(entry_json["mechanism"]),
                         (OverlayGraph.deserialise(og_json) for og_json in entry_json["overlay_graphs"]))


class MCSADB:
    def __init__(self, path: str, limit: int = 0):
        with open(path) as file:
            data = json.load(file)

        if limit <= 0:
            limit = len(data)

        self.entries: List[MCSAEntry] = [MCSAEntry.deserialize(entry_data) for entry_data in data[:limit]]
