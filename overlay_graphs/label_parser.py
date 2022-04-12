import re


_amino_pattern: re.Pattern = re.compile(r"^[\s]*Amino\([\s]*([a-zA-Z0-9+\-]+)[\s]*,[\s]*([a-zA-Z]+)[\s]*,"
                                        r"[\s]*([0-9]+)[^)]+\)[\s]*$")
_alias_pattern: re.Pattern = re.compile(r"^[\s]*Alias\([\s]*([a-zA-Z0-9+\-]+)[\s]*,[^)]+\)[\s]*$")


def is_amino_term(label: str) -> bool:
    return _amino_pattern.match(label) is not None


def is_alias_term(label: str) -> bool:
    return _alias_pattern.match(label) is not None


def is_term(label: str) -> bool:
    return is_amino_term(label) or is_alias_term(label)


def abstract_vertex_term_details(label: str, preserve_chain_position: bool = False) -> str:
    amino_acid_substitute = r"Amino(\1, \2, \3, *)" if preserve_chain_position else r"Amino(\1, \2, *, *)"

    label = _amino_pattern.sub(amino_acid_substitute, label)
    label = _alias_pattern.sub(r"Alias(\1, *)", label)

    return label


def remove_vertex_term(label: str) -> str:
    label = _amino_pattern.sub(r"\1", label)
    label = _alias_pattern.sub(r"\1", label)

    return label

