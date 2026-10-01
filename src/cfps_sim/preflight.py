"""Source-aware recipe chemistry and readiness gates; no implied calibration."""
from .economics import number, cost_recipe
from .recipe import validate


def audit_recipe(recipe, registry, prices=None, base_final_mM=None):
    checked = validate(recipe)
    totals, unsupported, unresolved, absent = {}, [], [], []
    for name, volume in checked["volumes_nl"].items():
        if not volume:
            continue
        stock = registry["stocks"].get(name)
        if stock is None:
            absent.append(name)
            continue
        concentration = stock.get("stock_mM")
        components = stock.get("components_per_formula_unit", {})
        if concentration is None or not components:
            unresolved.append(name)
        else:
            for chemical, coefficient in components.items():
                amount = number(concentration, name) * number(coefficient, chemical) * volume / 20000
                totals[chemical] = totals.get(chemical, 0) + amount
        if stock.get("model_coverage") == "unsupported" and name not in {
                "lysate", "base_buffer", "nuclease_free_water", "template_sfgfp",
                "template_petase", "template_reteplase"}:
            unsupported.append(name)
    base = base_final_mM if base_final_mM is not None else registry.get("base_buffer_components")
    if base is not None:
        for chemical, value in base.items():
            totals[chemical] = totals.get(chemical, 0) + number(value, chemical)
    conflicts = [row for row in registry.get("price_conflicts", [])
                 if checked["volumes_nl"].get(row["stock"], 0)]
    blockers = []
    if registry.get("status") != "confirmed" or registry.get("constraints", {}).get("status") != "confirmed":
        blockers.append("inventory_and_execution_constraints_unconfirmed")
    if base is None:
        blockers.append("base_buffer_composition_unknown")
    if conflicts:
        blockers.append("conflicting_price_sources")
    if absent:
        blockers.append("stocks_absent_from_registry")
    if set(unresolved) - {"lysate", "base_buffer", "nuclease_free_water", "template_sfgfp", "template_petase", "template_reteplase"}:
        blockers.append("stock_or_constituent_concentrations_unresolved")
    ledger = cost_recipe(recipe, prices) if prices is not None else None
    if ledger is None or not ledger["complete"]:
        blockers.append("complete_price_ledger_missing")
    if prices is None or prices.get("provenance", {}).get("status") != "confirmed":
        blockers.append("operative_price_ledger_unconfirmed")
    target = next(name.removeprefix("template_") for name in checked["volumes_nl"] if name.startswith("template_") and checked["volumes_nl"][name])
    return {"target": target, "volumes_nl": checked["volumes_nl"], "known_final_mM": totals,
            "base_composition_supplied": base is not None,
            "unresolved_stocks": sorted(unresolved), "absent_stocks": sorted(absent),
            "unsupported_model_effects": sorted(unsupported), "price_conflicts": conflicts,
            "cost": ledger, "submission_ready": not blockers, "blockers": blockers,
            "note": "Analytical totals, not free-ion concentrations. Unknown lysate pools are excluded."
            " Formula stoichiometry does not model binding, pH, or activity."}


def model_initial_conditions(audit, lysate_background_mM, other_initial):
    """Map supported pools only, requiring explicit backgrounds and other states.

    Supply all 20 amino-acid constituents in the audited composition/background.
    Fuel is an explicitly assumed lumped quantity, never inferred from glucose.
    Unsupported effects remain visible even when initial pools can be mapped.
    """
    from .construct import AMINO_ACIDS
    if audit["target"] != "sfgfp":
        raise ValueError("Dynamic initial-pool mapping currently supports sfGFP only")
    aliases = {"atp": "ATP", "gtp": "GTP", "ctp": "CTP", "utp": "UTP"}
    aa_names = dict(zip(AMINO_ACIDS, [
        "alanine", "cysteine", "aspartate", "glutamate", "phenylalanine", "glycine",
        "histidine", "isoleucine", "lysine", "leucine", "methionine", "asparagine",
        "proline", "glutamine", "arginine", "serine", "threonine", "valine",
        "tryptophan", "tyrosine"]))
    aliases.update({name: f"aa_{aa}" for aa, name in aa_names.items()})
    if not audit["base_composition_supplied"]:
        raise ValueError("Explicit base composition is required")
    unresolved = set(audit["unresolved_stocks"]) - {"lysate", "base_buffer", "nuclease_free_water", "template_sfgfp", "template_petase", "template_reteplase"}
    if unresolved or audit["absent_stocks"]:
        raise ValueError("Resolve all dispensed chemical stocks before mapping initial pools")
    missing = set(aliases) - set(lysate_background_mM)
    if missing:
        raise ValueError(f"Explicit lysate backgrounds required: {sorted(missing)}")
    required = {"ADP", "GDP", "Fuel", "mRNA", "Protein", "GFP", "Machinery"}
    if set(other_initial) != required:
        raise ValueError(f"Explicit other_initial must contain exactly {sorted(required)}")
    result = {state: audit["known_final_mM"].get(chemical, 0) + number(lysate_background_mM[chemical], chemical)
              for chemical, state in aliases.items()}
    result.update({name: number(value, name) for name, value in other_initial.items()})
    return {"initial": result, "model_status": "uncalibrated_partial_chemistry",
            "unsupported_model_effects": audit["unsupported_model_effects"],
            "recipe_blockers": audit["blockers"], "note": "Unknown concentrations must not be encoded as zero."
            " Explicit zero backgrounds are assumptions unless measured."}
