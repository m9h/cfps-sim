import csv
import json

import pytest

from cfps_sim.construct import read_construct
from cfps_sim.experiments import calibrate, make_synthetic_study, run_suite
from cfps_sim.model import load_model


@pytest.mark.gpu
def test_challenges_and_override_export(tmp_path):
    report = run_suite("constructs/competition/sfgfp.json", tmp_path/"suite")
    assert report["all_checks_pass"]
    assert report["initial_units"]["Machinery"] == "dimensionless"
    assert report["parameter_units"]["regen"] == "mM/min"
    model = load_model(read_construct("constructs/competition/sfgfp.json"), {"k_tl": 0.3}, {"GTP": 0.4})
    assert model["k_tl"] == 0.3
    assert model["GTP"] == 0.4
    with pytest.raises(ValueError, match="absent"):
        load_model(initial={"magnesium": 10})
    with pytest.raises(ValueError, match="positive"):
        load_model(parameters={"K_nt": 0})


@pytest.mark.gpu
def test_synthetic_recovery_holdout_and_no_truth_leakage(tmp_path):
    study_path = make_synthetic_study("constructs/competition/sfgfp.json", tmp_path/"study")
    study = json.loads(study_path.read_text())
    truth = study.pop("synthetic_truth")
    study_path.write_text(json.dumps(study))  # Fitting must not read generating values.
    result = calibrate(study_path, study_path.parent/"observations.csv", tmp_path/"fit")
    assert result["optimizer_success"], result
    assert result["optimizer"]["jacobian_step_relative_difference"] < 1e-3
    assert result["jacobian_rank"] == 2
    for name, value in truth.items():
        assert result["fitted_parameters"][name] == pytest.approx(value, rel=0.1)
    heldout = next(row for row in result["metrics"] if row["split"] == "test" and row["observable"] == "GFP")
    assert heldout["rmse_mM"] < heldout["baseline_rmse_mM"]
    # Withheld measured values cannot influence fitting.
    with (study_path.parent/"observations.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        if row["split"] == "test":
            row["value_mM"] = "1000"
    changed = tmp_path/"changed.csv"
    with changed.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    refit = calibrate(study_path, changed, tmp_path/"refit")
    assert refit["fitted_parameters"] == result["fitted_parameters"]
    rows[0]["split"] = "test"
    with changed.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="entire conditions"):
        calibrate(study_path, changed, tmp_path/"invalid")
