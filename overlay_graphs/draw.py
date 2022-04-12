import mod
import re


from enum import Enum
from overlay_graphs.util import convert_svg
from overlay_graphs.overlay_graph import OverlayGraph
from rdkit import Chem


_atom_charge_label_pattern: re.Pattern = re.compile(r"[0-9]*[+\-]")


class ChangeType(Enum):
    ADD = 1
    REMOVE = 2
    CONTEXT = 3
    MODIFIED = 4


def _get_type(element, overlay_graph: OverlayGraph):
    donated: int = overlay_graph.electrons_donated(element)
    received: int = overlay_graph.electrons_received(element)

    if donated == 0 and received == 0:
        return ChangeType.CONTEXT
    elif donated != 0 and received == 0:
        return  ChangeType.REMOVE
    elif donated == 0 and received != 0:
        return ChangeType.ADD
    else:
        return ChangeType.MODIFIED


def _get_charge(v_data):
    lbl: str = v_data["label"]
    neg_charge = lbl.count("-")
    pos_charge = lbl.count("+")
    if neg_charge > 0:
        return -neg_charge
    elif pos_charge > 0:
        return pos_charge
    else:
        return 0


def draw_overlay_graph(graph: OverlayGraph, path,
                   add_atom_indices: bool = False,
                   add_electron_count: bool = False,
                   add_context: bool = True,
                   abstract_hydrogens: bool = True,
                   include_vertices=None):
    amino_pattern: re.Pattern = re.compile(r"[\s]*Amino\([\s]*([a-zA-Z+\-]+)[\s]*,[\s]*([a-zA-Z]+)[^)]+\)[\s]*")
    alias_pattern: re.Pattern = re.compile(r"[\s]*Alias\([\s]*([a-zA-Z+\-]+)[\s]*,[^)]+\)[\s]*")
    nxgraph = graph.host_graph
    periodic_table = Chem.GetPeriodicTable()

    def label_converter(label: str):
        amino_match = amino_pattern.match(label)
        if amino_match is not None:
            return True, amino_pattern.sub(r"\2", label)

        alias_match = alias_pattern.match(label)
        if alias_match is not None:
            return True, alias_pattern.sub(r"\1", label)

        ignore_labels = ["R", "U", "1*", "2*", "Mg2"]
        if label.replace("-", "").replace("+", "") in ignore_labels:
            return True, label

        return False, label

    label2bond = {
        "-": Chem.rdchem.BondType.SINGLE,
        "?": Chem.rdchem.BondType.SINGLE,
        "=": Chem.rdchem.BondType.DOUBLE,
        "#": Chem.rdchem.BondType.TRIPLE,
        ">": Chem.rdchem.BondType.HYDROGEN,
        ":": Chem.rdchem.BondType.SINGLE
    }
    changetype2color = {
        ChangeType.ADD: (0, 1, 0),
        ChangeType.REMOVE: (1, 0, 0),
        ChangeType.MODIFIED: (0, 0, 1)
    }

    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    gid2molid = {}
    mol = Chem.EditableMol(Chem.Mol())
    highligt_atoms = []
    highlight_bonds = []
    highligt_atom_colors = {}
    highlight_bond_colors = {}
    for v, v_data in nxgraph.nodes(data=True):
        if include_vertices is not None and v not in include_vertices:
            continue

        v_type = _get_type(v, graph)
        if abstract_hydrogens and (not add_context or v_data["label"] == "H") and v_type == ChangeType.CONTEXT:
            is_only_context = True
            for e in nxgraph.edges(v):
                if _get_type(e, graph) != ChangeType.CONTEXT:
                    is_only_context = False
                    break

            if is_only_context:
                continue
        is_alias, lbl = label_converter(v_data["label"])
        if is_alias:
            atom = Chem.Atom("U")
            atom.SetProp("atomLabel", lbl)
        else:
            atom = Chem.Atom(re.sub(_atom_charge_label_pattern, "", lbl))
        atom.SetNoImplicit(True)
        atom_notes = []
        if add_atom_indices:
            atom_notes.append(str(v))
        if add_electron_count and v_type != ChangeType.CONTEXT:
            atom_notes.append(f"({graph.electrons_donated(v)}, {graph.electrons_received(v)})")
        if len(atom_notes) > 0:
            atom.SetProp("atomNote", ",".join(atom_notes))

        chg = _get_charge(v_data)
        if chg != 0:
            atom.SetFormalCharge(chg)
        idx = mol.AddAtom(atom)
        gid2molid[v] = idx
        if v_type != ChangeType.CONTEXT:
            highligt_atoms.append(idx)
            highligt_atom_colors[idx] = changetype2color[v_type]

    bond_notes = {}
    for src, tar in nxgraph.edges:
        e_type = _get_type((src, tar), graph)
        if not add_context and e_type == ChangeType.CONTEXT:
            continue

        e_data = nxgraph.edges[(src, tar)]
        if src not in gid2molid or tar not in gid2molid:
            continue
        bond_type = label2bond[e_data["label"]]
        # for some reason, EditableMol is 1 indexed, while GetMol() is 0  indexed...
        idx = mol.AddBond(gid2molid[src], gid2molid[tar], bond_type) - 1
        if e_type != ChangeType.CONTEXT:
            if add_electron_count:
                bond_notes[idx] = f"({graph.electrons_donated((src, tar))}, {graph.electrons_received((src, tar))})"
            highlight_bonds.append(idx)
            highlight_bond_colors[idx] = changetype2color[e_type]

    static_mol: Chem.Mol = mol.GetMol()
    for idx, note in bond_notes.items():
        static_mol.GetBondWithIdx(idx).SetProp("bondNote", note)

    static_mol.UpdatePropertyCache(strict=False)
    d = rdMolDraw2D.MolDraw2DSVG(500, 500)
    opts: Draw.MolDrawOptions = d.drawOptions()
    opts.useBWAtomPalette()
    rdMolDraw2D.PrepareAndDrawMolecule(d, static_mol,
                                       highlightAtoms=highligt_atoms,
                                       highlightBonds=highlight_bonds,
                                       highlightAtomColors=highligt_atom_colors,
                                       highlightBondColors=highlight_bond_colors)
    d.FinishDrawing()
    with open(path, "w") as f:
        f.write(d.GetDrawingText())
    return gid2molid


def print_overlay_graph(overlay_graph: OverlayGraph, add_atom_indices: bool = False, add_electron_count: bool = False,
                    add_context: bool = True, abstract_hydrogens: bool = True):
    prefix = mod.makeUniqueFilePrefix()

    draw_overlay_graph(overlay_graph, f"{prefix}og.svg", add_atom_indices, add_electron_count, add_context,
                   abstract_hydrogens)

    convert_svg(f'{prefix}og.svg', f'{prefix}og.pdf')

    with open(f"{prefix}og.tex", "w") as file:
        file.write("\centering\n")
        file.write(r"\includegraphics[width=0.6\textwidth]{" + f"{prefix}og.pdf" + "}")

    mod.post(f"\\summaryInput {prefix}og.tex")
