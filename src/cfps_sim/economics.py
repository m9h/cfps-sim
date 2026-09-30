"""Auditable stock-volume costs and target-specific production economics."""
import math

from .recipe import FIXED, TEMPLATES, validate


def number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f"{name}: requires a finite nonnegative number")
    return float(value)


def cost_recipe(recipe, prices):
    """Charge every dispensed stock, including fixed inputs and explicit water.

    Missing prices stay missing. The subtotal must not be used as a full cost.
    Mixed stock prices refer to the entire dispensed stock; do not also charge
    its constituents. This ledger reports reagent cost, excluding other costs.
    """
    if prices.get("currency") != "USD" or prices.get("price_basis") != "per_mL_stock":
        raise ValueError("Prices require USD and per_mL_stock basis")
    if prices.get("cost_scope") != "dispensed_reagents_only":
        raise ValueError("Declare cost_scope: dispensed_reagents_only")
    if not prices.get("provenance", {}).get("source"):
        raise ValueError("Price provenance.source is required")
    rates = prices["usd_per_ml"]
    for name, value in rates.items():
        if value is not None:
            number(value, f"price of {name}")
    checked = validate(recipe)
    lines, missing = [], []
    for name, volume in checked["volumes_nl"].items():
        if volume == 0:
            continue
        rate = rates.get(name)
        cost = None if rate is None else number(rate, f"price of {name}") * volume / 1_000_000
        if cost is None:
            missing.append(name)
        lines.append({"reagent": name, "volume_nl": volume, "usd_per_ml_stock": rate,
                      "cost_usd": cost, "fixed_input": name in FIXED or name in TEMPLATES})
    subtotal = sum(row["cost_usd"] for row in lines if row["cost_usd"] is not None)
    complete = not missing
    return {"currency": "USD", "cost_scope": prices["cost_scope"], "complete": complete,
            "reaction_volume_l": sum(checked["volumes_nl"].values()) / 1_000_000_000,
            "known_subtotal_usd": subtotal, "total_cost_usd": subtotal if complete else None,
            "missing_prices": sorted(missing), "line_items": lines,
            "price_provenance": prices["provenance"],
            "note": "Dispensed reagent costs only. No labor, plates, assays, purification, compute, or scale-up costs are included. Missing prices are never treated as zero."}


def score_recipe(recipe, prices, observation):
    report = cost_recipe(recipe, prices)
    if observation.get("data_kind") not in {"synthetic", "experimental"}:
        raise ValueError("Yield data_kind must be synthetic or experimental")
    if not observation.get("provenance", {}).get("source") or not observation.get("readout"):
        raise ValueError("Yield source provenance and readout are required")
    template = next(name for name in TEMPLATES if recipe["volumes_nl"].get(name, 0))
    target = template.removeprefix("template_")
    if observation.get("target", "").lower() != target:
        raise ValueError("Yield target must match the recipe template")
    value = number(observation["value"], "yield")
    units = observation["units"]
    if units == "g/L":
        amount, amount_units, score_units = value * report["reaction_volume_l"], "g", "USD/g"
    elif units == "U/mL":
        amount, amount_units, score_units = value * report["reaction_volume_l"] * 1000, "U", "USD/U"
    else:
        raise ValueError("Use g/L or U/mL with an explicit validated assay conversion; no implicit apparent-g/L conversion")
    if report["complete"] and report["total_cost_usd"] <= 0:
        raise ValueError("Scoring requires positive total reaction cost")
    scorable = report["complete"] and amount > 0
    report.update({"observation": observation, "product_amount": amount, "product_units": amount_units,
                   "score_units": score_units, "scorable": scorable,
                   "specific_cost": report["total_cost_usd"] / amount if scorable else None,
                   "unscorable_reason": None if scorable else "missing_prices" if not report["complete"] else "zero_product",
                   "objective_direction": "minimize_specific_cost"})
    return report


def pareto_indices(costs, yields):
    """Indices not dominated by lower/equal cost and higher/equal yield."""
    if len(costs) != len(yields):
        raise ValueError("Costs and yields must have equal length")
    costs = [number(value, "cost") for value in costs]
    yields = [number(value, "yield") for value in yields]
    return [i for i, (cost, output) in enumerate(zip(costs, yields))
            if not any(c <= cost and y >= output and (c < cost or y > output)
                       for c, y in zip(costs, yields))]
