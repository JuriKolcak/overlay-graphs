import mod


from overlay_graphs.make_overlay_graphs import overlay_graphs_for_mechanisms
from overlay_graphs.mechanism import Mechanism
from overlay_graphs.sanitiser import MechanismSanitiser
from overlay_graphs.util import load_mechanisms
from typing import Iterable, List


def sanitize_mechanisms(mechanisms: List[Mechanism], verbosity: int = 0) -> Iterable[Mechanism]:
    sanitiser = MechanismSanitiser(
        [mod.smiles("O", name="Water", add=False),
         mod.smiles("[OH3+]", name="Hydronium", add=False),
         mod.smiles("[OH-]", name="Hydroxide", add=False),
         mod.smiles("O=O", name="Oxygen Molecule", add=False),
         mod.smiles("N#N", name="Nitrogen Molecule", add=False),
         mod.smiles("O=C=O", name="Carbon Dioxide", add=False)],
        preserve_peptide_chain_positions=True,
        ignore_dative_bonds=True
    )

    sanitised_mechanisms = list(sanitiser.sanitise_mechanisms(mechanisms, verbosity))

    unsanitised_mechanisms = [mechanism for mechanism in mechanisms if all(mechanism.entry != sanitised_mechanism.entry or
                                                                           mechanism.number != sanitised_mechanism.number
                                                                           for sanitised_mechanism in sanitised_mechanisms)]
    if verbosity >= 1:
        print(f"Number of incompatible proposals {len(unsanitised_mechanisms)} out of {len(mechanisms)}")
        print(f"One step proposals {len([mechanism for mechanism in sanitised_mechanisms if len(mechanism.steps) == 1])}")
        print(f"Processed a total of {sum(len(mechanism.steps) for mechanism in sanitised_mechanisms)} reaction steps")

        if verbosity >= 3:
            print("\nIncompatible proposals: "+" ".join(map(str, unsanitised_mechanisms)))

    return sanitised_mechanisms


def _compute_overlay_graphs():
    mechanisms = list(load_mechanisms())[:8]

    sanitised_mechanisms = list(sanitize_mechanisms(mechanisms))

    overlay_graphs_for_mechanisms(sanitised_mechanisms)


if __name__ == "__main__":
    _compute_overlay_graphs()
