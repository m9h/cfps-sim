"""Normalize a cached competition inventory without declaring it authoritative."""
import argparse
import hashlib
import json
import re
from pathlib import Path

ALIASES = {"cell_lysate": "lysate"}
STOICHIOMETRY = {
    "magnesium_glutamate": {"magnesium": 1, "glutamate": 2},
    "potassium_glutamate": {"potassium": 1, "glutamate": 1},
    "ammonium_glutamate": {"ammonium": 1, "glutamate": 1},
    "potassium_phosphate_dibasic": {"potassium": 2, "phosphate": 1},
    "potassium_phosphate_monobasic": {"potassium": 1, "phosphate": 1},
}
SUPPORTED = {"atp": "ATP", "ctp": "CTP", "gtp": "GTP", "utp": "UTP",
             "tyrosine": "aa_Y", "cysteine": "aa_C"}


def normalize(source, draft):
    raw = json.loads(source.read_text())
    stocks, conflicts = {}, []
    for group in raw["reagentGroups"]:
        for reagent in group["reagents"]:
            name = ALIASES.get(reagent["id"], reagent["id"])
            label = reagent.get("concentration", "")
            matched = re.fullmatch(r"([0-9.]+)M", label)
            concentration = float(matched[1]) * 1000 if matched else None
            # Mixed-stock labels do not establish each constituent's concentration.
            components = (STOICHIOMETRY.get(name, {name: 1}) if concentration is not None
                          and "mix" not in name and name != "aa_mix_17" else {})
            if name == "dna_template":
                ids = ["template_sfgfp", "template_petase", "template_reteplase"]
            else:
                ids = [name]
            for stock_id in ids:
                cached_price = reagent.get("costPerMl")
                draft_price = draft["usd_per_ml"].get(stock_id)
                if stock_id in draft["usd_per_ml"] and draft_price != cached_price:
                    conflicts.append({"stock": stock_id, "cached_usd_per_ml": cached_price,
                                      "draft_usd_per_ml": draft_price})
                stocks[stock_id] = {
                    "name": reagent["name"], "source_id": reagent["id"],
                    "stock_label": label, "stock_mM": concentration,
                    "components_per_formula_unit": components,
                    "components_status": "formula_accounting" if components else "unresolved",
                    "model_species": SUPPORTED.get(name),
                    "model_coverage": "initial_pool_only" if name in SUPPORTED else "unsupported",
                    "cached_in_stock": reagent.get("in_stock"),
                    "availability_status": "unconfirmed_cached_snapshot",
                    "cached_usd_per_ml": cached_price, "draft_usd_per_ml": draft_price,
                }
    return {"schema_version": 1, "status": "draft_requires_reconciliation",
            "source": {"url": "https://coopetition.htgaa.org/cfps", "artifact": source.name,
                       "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                       "status": "local_cache_capture_date_unknown"},
            "constraints": {"reaction_nl": 20000, "transfer_increment_nl": 25,
                            "lysate_nl": 5000, "base_buffer_nl": 2000, "template_nl": 2000,
                            "status": "existing_project_assumptions_require_current_confirmation"},
            "base_buffer_components": None,
            "unknowns": ["current_inventory", "operative_price_ledger", "base_buffer_composition",
                         "lysate_background_pools", "mixed_stock_composition", "target_assay_conversion"],
            "price_conflicts": conflicts, "stocks": stocks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("--draft", type=Path, default=Path("examples/prices.draft.json"))
    parser.add_argument("--output", type=Path, default=Path("data/registry/inventory.draft.json"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(normalize(args.source, json.loads(args.draft.read_text())), indent=2)+"\n")
