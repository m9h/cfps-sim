import json
from pathlib import Path
import numpy as np
import pytest
from cfps_sim.recipe import validate
from cfps_sim.model import simulate


def test_recipe_volume_and_concentration():
    recipe = json.loads(Path("examples/recipe.json").read_text())
    result = validate(recipe)
    assert sum(result["volumes_nl"].values()) == 20000
    assert result["known_final_mM"]["magnesium_glutamate"] == 6.25
    recipe["volumes_nl"]["atp"] = 301
    with pytest.raises(ValueError):
        validate(recipe)


@pytest.mark.gpu
def test_no_substrate_no_expression_and_resource_balance():
    rr, trajectory = simulate()
    assert np.isfinite(trajectory).all()
    assert np.min(trajectory[:, 1:]) >= -1e-9
    assert trajectory[-1, 9] > 0
    # GTP + GDP + nucleotide incorporated into RNA remains balanced,
    # with degraded RNA leaving the modeled pool.
    assert rr["GTP"] + rr["GDP"] <= 1.5 + 1e-8
    assert rr["AA"] + 240*(rr["Protein"] + rr["GFP"]) == pytest.approx(40, abs=1e-7)
    assert rr["ATP"] + rr["ADP"] <= 2 + 1e-8
    rr.resetAll()
    for species in ["ATP", "ADP", "GTP", "CTP", "UTP"]:
        rr[species] = 0
    rr.simulate(0, 240, 10)
    assert rr["GFP"] == pytest.approx(0, abs=1e-12)
