from copy import deepcopy
import json
from pathlib import Path

import pytest

from cfps_sim.economics import cost_recipe, pareto_indices, score_recipe


def inputs():
    recipe = json.loads(Path("examples/recipe.json").read_text())
    prices = json.loads(Path("examples/prices.draft.json").read_text())
    observation = json.loads(Path("examples/yield.synthetic.json").read_text())
    return recipe, prices, observation


def test_price_completeness_and_units():
    recipe, prices, observation = inputs()
    incomplete = score_recipe(recipe, prices, observation)
    assert incomplete["missing_prices"] == ["base_buffer"]
    assert incomplete["specific_cost"] is None
    assert incomplete["known_subtotal_usd"] == pytest.approx(0.02448395)
    prices["usd_per_ml"]["base_buffer"] = 0.1  # Deliberate test assumption.
    score = score_recipe(recipe, prices, observation)
    assert score["reaction_volume_l"] == pytest.approx(20e-6)
    assert score["product_amount"] == pytest.approx(20e-6)
    assert score["total_cost_usd"] == pytest.approx(0.02468395)
    assert score["specific_cost"] == pytest.approx(1234.1975)
    doubled = deepcopy(observation)
    doubled["value"] *= 2
    assert score_recipe(recipe, prices, doubled)["specific_cost"] == pytest.approx(score["specific_cost"]/2)
    observation["value"] = 0
    assert score_recipe(recipe, prices, observation)["unscorable_reason"] == "zero_product"


def test_activity_target_and_invalid_inputs():
    recipe, prices, observation = inputs()
    prices["usd_per_ml"]["base_buffer"] = 0.1
    recipe["volumes_nl"]["template_petase"] = recipe["volumes_nl"].pop("template_sfgfp")
    with pytest.raises(ValueError, match="target"):
        score_recipe(recipe, prices, observation)
    observation.update(target="petase", value=10.0, units="U/mL", readout="declared activity assay")
    score = score_recipe(recipe, prices, observation)
    assert score["product_amount"] == pytest.approx(0.2)
    assert score["score_units"] == "USD/U"
    prices["usd_per_ml"]["atp"] = float("nan")
    with pytest.raises(ValueError, match="price"):
        cost_recipe(recipe, prices)
    prices["usd_per_ml"]["atp"] = 1.539
    observation["units"] = "apparent g/L"
    with pytest.raises(ValueError, match="conversion"):
        score_recipe(recipe, prices, observation)


def test_cost_yield_tradeoff_and_pareto():
    # The maximum-yield point can have worse economics; remove only domination.
    costs, yields = [0.05, 0.08, 0.04, 0.06], [1.0, 1.2, 0.9, 0.8]
    assert min(range(4), key=lambda i: costs[i]/yields[i]) == 2
    assert max(range(4), key=lambda i: yields[i]) == 1
    assert pareto_indices(costs, yields) == [0, 1, 2]
