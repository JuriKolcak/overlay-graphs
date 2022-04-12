import json
import subprocess as sp


from overlay_graphs.mechanism import Mechanism, Step
from typing import Dict, Iterable, List, Tuple


def convert_svg(_svg: str, _pdf: str):
    sp.run(['rsvg-convert', '-f', 'pdf', '-o', _pdf, _svg])


def load_mechanisms(mechanism_file: str = "mechanisms.json") -> Iterable[Mechanism]:
    with open(mechanism_file, "r") as file:
        mechanisms_json = json.load(file)

    mechanisms: Dict[Tuple[int, int], List[Step]] = {}

    for step_json in mechanisms_json:
        step_json["gml"] = step_json["gml"].replace(">", ":")

        entry = step_json["entry"]
        mechanism = step_json["proposal"]
        step = Step.deserialise(step_json)

        if (entry, mechanism) not in mechanisms:
            mechanisms[(entry, mechanism)] = []

        mechanisms[(entry, mechanism)].append(step)

    for (entry, number), steps in mechanisms.items():
        yield Mechanism(entry, number, steps)
