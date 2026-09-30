"""Volume accounting; chemical composition must be supplied explicitly."""
import math

FIXED = {"lysate": 5000, "base_buffer": 2000}
TEMPLATES = {"template_sfgfp", "template_petase", "template_reteplase"}


def validate(recipe):
    volumes = recipe["volumes_nl"]
    for name, volume in volumes.items():
        if isinstance(volume, bool) or not isinstance(volume, (int, float)):
            raise ValueError(f"{name}: volume must be numeric")
        if not math.isfinite(volume) or volume < 0 or volume % 25:
            raise ValueError(f"{name}: volume must be a nonnegative multiple of 25 nL")
    for name, volume in FIXED.items():
        if volumes.get(name) != volume:
            raise ValueError(f"{name}: requires {volume} nL")
    templates = [name for name in TEMPLATES if volumes.get(name, 0) > 0]
    if len(templates) != 1 or volumes[templates[0]] != 2000:
        raise ValueError("Exactly one DNA template at 2000 nL is required")
    water = 20000 - sum(v for k, v in volumes.items() if k != "nuclease_free_water")
    if water < 0:
        raise ValueError("Recipe exceeds 20000 nL")
    if "nuclease_free_water" in volumes and volumes["nuclease_free_water"] != water:
        raise ValueError(f"Water must be {water} nL")
    stocks = recipe.get("stocks_mM", {})
    for name, concentration in stocks.items():
        if not math.isfinite(concentration) or concentration < 0:
            raise ValueError(f"{name}: invalid stock concentration")
    additions = {name: stocks[name] * v / 20000 for name, v in volumes.items() if name in stocks}
    base = recipe.get("base_final_mM", {})
    for name, concentration in base.items():
        if not math.isfinite(concentration) or concentration < 0:
            raise ValueError(f"{name}: invalid base concentration")
    known = {name: additions.get(name, 0) + base.get(name, 0) for name in additions.keys() | base.keys()}
    return {"volumes_nl": dict(volumes, nuclease_free_water=water),
            "supplement_mM": additions, "known_final_mM": known,
            "unmapped_reagents": sorted(set(volumes) - set(stocks) - set(FIXED) - TEMPLATES - {"nuclease_free_water"}),
            "note": "Known additions only; endogenous lysate pools are excluded. Base composition is user supplied."}
