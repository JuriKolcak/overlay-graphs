from overlay_graphs.substrate_rules import substrate_rules_for_entries
from overlay_graphs.mcsadb import MCSADB


def _compute_substrate_rules():
    db: MCSADB = MCSADB("overlay_graphs.json")

    substrate_rules_for_entries(db.entries)


if __name__ == "__main__":
    _compute_substrate_rules()
