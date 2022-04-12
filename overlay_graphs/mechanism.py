import mod


from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple


class ECNumber:
    def __init__(self, ec_id: str):
        self._levels: Tuple[str] = tuple(ec_id.split("."))

    def __eq__(self, other: 'ECNumber') -> bool:
        return self._levels == other._levels

    def __hash__(self) -> int:
        return hash(self._levels)

    def __ge__(self, other: 'ECNumber') -> bool:
        return self > other or self == other

    def __gt__(self, other: 'ECNumber') -> bool:
        return other < self

    def __le__(self, other: 'ECNumber') -> bool:
        return self < other or self == other

    def __lt__(self, other: 'ECNumber') -> bool:
        for i in range(0, 4):
            if self._levels[i] == other._levels[i]:
                continue

            if self._levels[i] == "-":
                return True

            if other._levels[i] == "-":
                return False

            if self._levels[i].isnumeric():
                if not other._levels[i].isnumeric():
                    return True
                if int(self._levels[i]) != int(other._levels[i]):
                    return int(self._levels[i]) < int(other._levels[i])

            if other._levels[i].isnumeric():
                return False

            return self._levels[i] < other._levels[i]

        return False

    def __iter__(self) -> Iterator[str]:
        return iter(self._levels)

    def __getitem__(self, key: int) -> str:
        return self._levels[key]

    def __contains__(self, ec: 'ECNumber') -> bool:
        for i in range(0, 4):
            if self._levels[i] != ec._levels[i] and self._levels[i] != "-":
                return False

        return True

    def __str__(self) -> str:
        return ".".join(self._levels)

    def abstract(self, level: int) -> 'ECNumber':
        return ECNumber(".".join(list(self._levels[:level]) + ["-"] * (4 - level)))


class Step:
    def __init__(self, entry: int, mechanism: int, number: int, rule: Optional[mod.Rule] = None):
        self._entry: int = entry
        self._mechanism: int = mechanism
        self._number: int = number
        self._rule: Optional[mod.Rule] = rule

        self._components: List[str] = []

    def __str__(self):
        rule_name = self.rule.name if self.rule is not None else "None"
        return f'entry: {self.entry}, mechanism: {self.mechanism}, step: {self.number}, rule: {rule_name}'

    @property
    def entry(self) -> int:
        return self._entry

    @property
    def mechanism(self) -> int:
        return self._mechanism

    @property
    def number(self) -> int:
        return self._number

    @property
    def rule(self) -> Optional[mod.Rule]:
        return self._rule

    @property
    def components(self):
        return self._components

    @staticmethod
    def deserialise(json_object: Dict[str, Any], load_rule: bool = True, verbosity: int = 0) -> 'Step':
        step = Step(json_object["entry"], json_object["proposal"], json_object["step"])

        if load_rule and "gml" in json_object:
            try:
                step._rule = mod.ruleGMLString(json_object["gml"], add=False)
            except mod.InputError as error:
                if verbosity > 0:
                    print(f"#\tError loading gml rule for Step {step.entry}_{step.mechanism}_{step.number}: '{error}'")

        if "components" in json_object:
            for component in json_object["components"]:
                step.add_component(component)

        return step

    def add_component(self, component: str):
        self._components.append(component)

    def serialise(self) -> Dict[str, Any]:
        json_object = {
            'entry': self.entry,
            'proposal': self.mechanism,
            'step': self.number
        }

        if self.rule is not None:
            json_object["gml"] = self.rule.getGMLString()

        if len(self._components) > 0:
            json_object["components"] = self._components

        return json_object


class Mechanism:
    def __init__(self, entry: int, number: int, steps: Iterable[Step]):
        self._entry: int = entry
        self._number: int = number
        self._steps: List[Step] = list(sorted(steps, key=lambda s: s.number))

        self._ec: Optional[ECNumber] = None

    def __eq__(self, other: 'Mechanism') -> bool:
        return self.entry == other.entry and self.number == other.number

    def __getitem__(self, index: int) -> Step:
        return self._steps[index]

    def __hash__(self) -> int:
        return hash(pow(self.entry, self.number))

    def __iter__(self) -> Iterator[Step]:
        return iter(self._steps)

    def __len__(self) -> int:
        return len(self._steps)

    def __ne__(self, other: 'Mechanism') -> bool:
        return not self == other

    def __next__(self) -> Step:
        return next(iter(self))

    def __str__(self):
        return f'Proposal(entry: {self.entry}, mechanism: {self.number}, #steps: {len(self.steps)})'

    @property
    def entry(self) -> int:
        return self._entry

    @property
    def number(self) -> int:
        return self._number

    @property
    def steps(self) -> List[Step]:
        return list(self._steps)

    def _get_ec(self) -> Optional[ECNumber]:
        return self._ec

    def _set_ec(self, value: ECNumber):
        self._ec = value

    ec = property(fget=_get_ec, fset=_set_ec)

    @staticmethod
    def deserialise(mechanism_json: Dict[str, Any]) -> 'Mechanism':
        mechanism = Mechanism(mechanism_json['entry'], mechanism_json['proposal'],
                              [Step.deserialise(step_json) for step_json in mechanism_json['steps']])

        if "ec" in mechanism_json:
            mechanism.ec = ECNumber(mechanism_json["ec"])

        return mechanism

    def serialise(self) -> Dict[str, Any]:
        mechanism_json = {
            'entry': self.entry,
            'proposal': self.number,
            'steps': [step.serialise() for step in self.steps]
        }

        if self.ec is not None:
            mechanism_json["ec"] = str(self.ec)

        return mechanism_json
