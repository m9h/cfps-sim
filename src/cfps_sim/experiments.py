"""Reproducible in-silico challenges and conditional parameter calibration.

The demo uses synthetic observations generated from this same reduced model.
It tests numerical plumbing, not biological accuracy or model adequacy.
"""
from copy import deepcopy
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from .construct import AMINO_ACIDS, demands, read_construct
from .model import PARAMETERS, load_model, simulate


def write_rows(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_suite(construct_path, output, seed=2026):
    construct = read_construct(construct_path)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rr = load_model(construct)
    baseline_parameters = {name: float(rr[name]) for name in sorted(PARAMETERS)}
    baseline_initial = {name: float(rr[name]) for name in ["ATP", "ADP", "GTP", "GDP", "CTP", "UTP", "Fuel", "mRNA", "Protein", "GFP", "Machinery"] + [f"aa_{aa}" for aa in AMINO_ACIDS]}
    cases = [
        ("reference", {}, {}),
        ("no_template", {"dna_template": 0}, {}),
        ("no_transcription", {"k_tx": 0}, {}),
        ("no_translation", {"k_tl": 0}, {}),
        ("no_maturation", {"k_mat": 0}, {}),
        ("no_regeneration", {"regen": 0}, {}),
        ("no_fuel", {}, {"Fuel": 0}),
        ("no_gtp", {}, {"GTP": 0}),
        ("no_leucine", {}, {"aa_L": 0}),
        ("no_machinery", {}, {"Machinery": 0}),
    ]
    rows, manifest_cases = [], []
    max_balance_error = 0.0
    min_concentration = 0.0
    stats = demands(construct)
    for name, parameters, initial in cases:
        model, trajectory = simulate(construct=construct, parameters=parameters, initial=initial)
        protein = float(model["Protein"] + model["GFP"])
        for aa, count in stats["aa_per_protein"].items():
            pool = initial.get(f"aa_{aa}", 2.0)
            max_balance_error = max(max_balance_error, abs(model[f"aa_{aa}"] + count*protein-pool))
        min_concentration = min(min_concentration, float(trajectory[:, 1:].min()))
        rows.append({"case": name, "gfp_mM": float(model["GFP"]), "total_protein_mM": protein,
                     "mrna_mM": float(model["mRNA"]), "atp_mM": float(model["ATP"]), "gtp_mM": float(model["GTP"])})
        manifest_cases.append({"case": name, "parameters": parameters, "initial_values": initial})
        np.savetxt(output/f"{name}.csv", trajectory, delimiter=",", header=",".join(trajectory.colnames), comments="")
    reference = rows[0]["gfp_mM"]
    sensitivity = []
    for parameter in ["k_tx", "k_tl", "regen", "k_decay", "k_mrna", "k_mat"]:
        endpoints = []
        for factor in [0.9, 1.1]:
            model, _ = simulate(construct=construct, parameters={parameter: baseline_parameters[parameter]*factor})
            endpoints.append(float(model["GFP"]))
        elasticity = np.log(endpoints[1]/endpoints[0])/np.log(1.1/0.9)
        sensitivity.append({"parameter": parameter, "log_endpoint_elasticity": float(elasticity)})
    grid = []
    for atp in [0.25, 0.75, 1.5, 3.0]:
        for gtp in [0.25, 0.75, 1.5, 3.0]:
            model, _ = simulate(construct=construct, initial={"ATP": atp, "GTP": gtp})
            grid.append({"atp_initial_mM": atp, "gtp_initial_mM": gtp, "gfp_mM": float(model["GFP"])})
    capacity = []
    for tx_scale in [0.25, 0.5, 1.0, 2.0, 4.0]:
        for tl_scale in [0.25, 0.5, 1.0, 2.0, 4.0]:
            model, _ = simulate(construct=construct, parameters={"k_tx": baseline_parameters["k_tx"]*tx_scale, "k_tl": baseline_parameters["k_tl"]*tl_scale})
            capacity.append({"tx_scale": tx_scale, "tl_scale": tl_scale, "gfp_mM": float(model["GFP"]), "mrna_mM": float(model["mRNA"])})
    rng = np.random.default_rng(seed)
    ensemble = []
    for index in range(64):
        parameters = {name: baseline_parameters[name]*float(np.exp(rng.normal(0, 0.2))) for name in ["k_tx", "k_tl", "regen", "k_decay", "k_mrna", "k_mat"]}
        initial = {name: baseline_initial[name]*float(np.exp(rng.normal(0, 0.02))) for name in ["ATP", "GTP", "CTP", "UTP", "Fuel"]}
        model, _ = simulate(construct=construct, parameters=parameters, initial=initial)
        ensemble.append(dict(draw=index, gfp_mM=float(model["GFP"]), **parameters, **{f"initial_{name}_mM": value for name, value in initial.items()}))
    # An uncertainty scenario, not a newly designed competition template.
    shortened = deepcopy(construct)
    shortened["id"] += "-assumed-cds-end"
    shortened["transcript"]["end"] = shortened["cds"]["end"]
    shortened["transcript"]["basis"] = "Sensitivity scenario: assumed transcript ends at CDS stop; no experimental verification"
    short_model, _ = simulate(construct=shortened)
    refined = load_model(construct)
    refined.integrator.relative_tolerance = 1e-10
    refined.integrator.absolute_tolerance = 1e-14
    refined.simulate(0, 240, 481)
    convergence = abs(float(refined["GFP"])-reference)/reference
    zero_cases = {"no_template", "no_transcription", "no_translation", "no_maturation", "no_gtp", "no_leucine", "no_machinery"}
    checks = {"finite_endpoints": bool(all(np.isfinite(row["gfp_mM"]) for row in rows)),
              "nonnegative_trajectories": min_concentration >= -1e-9,
              "individual_aa_balance": max_balance_error < 1e-7,
              "required_resource_controls": all(abs(row["gfp_mM"]) < 1e-12 for row in rows if row["case"] in zero_cases),
              "solver_convergence": convergence < 1e-5}
    report = {"status": "uncalibrated_computational_challenges", "construct": stats,
              "baseline_parameters": baseline_parameters, "baseline_initial_values": baseline_initial,
              "initial_units": {name: "dimensionless" if name == "Machinery" else "mM" for name in baseline_initial}, "duration_minutes": 240,
              "parameter_units": {name: "mM" if name == "dna_template" or name.startswith("K_") else "mM/min" if name == "regen" else "1/min" for name in baseline_parameters},
              "cases": manifest_cases, "checks": checks, "all_checks_pass": all(checks.values()),
              "maximum_aa_balance_error_mM": max_balance_error, "solver_relative_endpoint_difference": convergence,
              "ensemble": {"seed": seed, "n": 64, "gfp_mM_quantiles_05_50_95": np.quantile([row["gfp_mM"] for row in ensemble], [0.05, 0.5, 0.95]).tolist(),
                           "assumption": "Independent illustrative lognormal perturbations: log-SD 0.2 for six kinetics, 0.02 for five initial pools. Scenario spread, not calibrated prediction intervals or measured dispensing noise."},
              "transcript_scenarios": [{"boundary": "supplied_sequence_end", "ntp_counts": stats["ntp_per_transcript"], "gfp_mM": reference},
                                       {"boundary": "assumed_cds_end", "ntp_counts": demands(shortened)["ntp_per_transcript"], "gfp_mM": float(short_model["GFP"])}],
              "note": "Initial pools and parameter changes are model-space interventions, not feasible competition recipes. Magnesium, pH, oxygen, and enzyme-specific folding are absent; no predictions for them are claimed."}
    write_rows(output/"challenges.csv", rows)
    write_rows(output/"sensitivity.csv", sensitivity)
    write_rows(output/"atp_gtp_grid.csv", grid)
    write_rows(output/"tx_tl_grid.csv", capacity)
    write_rows(output/"ensemble.csv", ensemble)
    (output/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    return report


def make_synthetic_study(construct_path, output, seed=2026):
    construct = read_construct(construct_path)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    truth = {"k_tl": 0.26, "k_decay": 0.003}
    conditions = [{"id": "reference", "initial_mM": {}, "split": "train"},
                  {"id": "low_gtp", "initial_mM": {"GTP": 0.4}, "split": "train"},
                  {"id": "low_fuel", "initial_mM": {"Fuel": 1}, "split": "train"},
                  {"id": "withheld_recipe", "initial_mM": {"GTP": 0.8, "Fuel": 4}, "split": "test"}]
    rows = []
    for condition in conditions:
        model = load_model(construct, truth, condition["initial_mM"])
        values = model.simulate(times=list(range(0, 241, 30)), selections=["time", "GFP", "mRNA"])
        for column, observable in [(1, "GFP"), (2, "mRNA")]:
            sigma = float(max(values[:, column].max()*0.05, 1e-12))
            for index in range(1, len(values)):
                for replicate in range(1, 3):
                    rows.append({"condition": condition["id"], "split": condition["split"], "replicate": replicate,
                                 "time_min": float(values[index, 0]), "observable": observable,
                                 "value_mM": float(values[index, column]+rng.normal(0, sigma)), "sigma_mM": sigma})
    write_rows(output/"observations.csv", rows)
    study = {"schema_version": 1, "data_kind": "synthetic", "seed": seed,
             "construct": construct, "conditions": conditions,
             "fit_parameters": {"k_tl": {"lower": 0.05, "upper": 0.6}, "k_decay": {"lower": 0.0005, "upper": 0.008}},
             "synthetic_truth": truth,
             "provenance": {"source": "cfps-sim reduced model; not Horvath's experimental data", "noise": "Independent Gaussian; known SD = 5% of channel maximum; background-corrected values may be negative"},
             "note": "Same-model recovery exercise. All other parameters, initial pools, transcript bounds, and observation scale are fixed. This is not biological calibration."}
    (output/"study.json").write_text(json.dumps(study, indent=2)+"\n")
    return output/"study.json"


def calibrate(study_path, observations_path, output):
    study = json.loads(Path(study_path).read_text())
    if study.get("schema_version") != 1 or study.get("data_kind") not in {"synthetic", "experimental"} or not study.get("provenance", {}).get("source"):
        raise ValueError("Study requires schema v1, data kind, and source provenance")
    with Path(observations_path).open() as handle:
        observations = list(csv.DictReader(handle))
    if not observations:
        raise ValueError("No observations supplied")
    conditions = {c["id"]: c for c in study["conditions"]}
    if len(conditions) != len(study["conditions"]):
        raise ValueError("Condition identifiers must be unique")
    if any("Machinery" in c["initial_mM"] for c in conditions.values()):
        raise ValueError("Machinery is dimensionless and cannot be supplied in initial_mM")
    models = {key: load_model(study["construct"], initial=c["initial_mM"]) for key, c in conditions.items()}
    names = list(study["fit_parameters"])
    if not names or set(names)-PARAMETERS:
        raise ValueError("Fit only explicitly supported parameters")
    bounds = np.array([[study["fit_parameters"][name][key] for name in names] for key in ["lower", "upper"]], dtype=float)
    if not np.isfinite(bounds).all() or np.any(bounds[0] <= 0) or np.any(bounds[0] >= bounds[1]):
        raise ValueError("Log-space bounds must be finite, positive, and ordered")
    for row in observations:
        if row["condition"] not in conditions or row["split"] != conditions[row["condition"]]["split"] or row["split"] not in {"train", "test"}:
            raise ValueError("Assign entire conditions to train/test; do not split replicates")
        for key in ["time_min", "value_mM", "sigma_mM"]:
            row[key] = float(row[key])
            if not np.isfinite(row[key]):
                raise ValueError("Measurements, times, and errors must be finite")
        if row["time_min"] <= 0 or row["sigma_mM"] <= 0 or row["observable"] not in {"GFP", "mRNA", "ATP", "GTP", "CTP", "UTP"}:
            raise ValueError("Positive times/errors and supported concentration observables are required")
    train = np.array([r["split"] == "train" for r in observations])
    if train.sum() <= len(names):
        raise ValueError("Insufficient training observations")
    defaults = np.array([float(next(iter(models.values()))[name]) for name in names])

    def predict(parameters):
        lookup = {}
        for key, model in models.items():
            requested = [r for r in observations if r["condition"] == key]
            if not requested:
                continue
            times = sorted({r["time_min"] for r in requested})
            observables = sorted({r["observable"] for r in requested})
            model.resetAll()
            for name, value in zip(names, parameters):
                model[name] = float(value)
            values = model.simulate(times=[0]+times, selections=["time"]+observables)
            for i, time in enumerate(times, 1):
                for j, observable in enumerate(observables, 1):
                    lookup[(key, time, observable)] = float(values[i, j])
        return np.array([lookup[(r["condition"], r["time_min"], r["observable"])] for r in observations])

    measured = np.array([r["value_mM"] for r in observations])
    sigma = np.array([r["sigma_mM"] for r in observations])
    baseline = predict(defaults)
    residual = lambda log_parameters: ((predict(np.exp(log_parameters))-measured)/sigma)[train]
    starts = [np.clip(defaults, bounds[0]*1.001, bounds[1]*0.999), np.sqrt(bounds[0]*bounds[1])]
    fits = [least_squares(residual, np.log(start), bounds=np.log(bounds), diff_step=1e-4, max_nfev=200) for start in starts]
    fit = min(fits, key=lambda result: result.cost)
    predictions = predict(np.exp(fit.x))
    singular_values = np.linalg.svd(fit.jac, compute_uv=False)
    rank = int(np.sum(singular_values > singular_values[0]*1e-8))
    intervals = None
    if rank == len(names):
        covariance = np.linalg.inv(fit.jac.T@fit.jac)
        sd = np.sqrt(np.maximum(np.diag(covariance), 0))
        intervals = {name: np.exp([value-1.96*error, value+1.96*error]).tolist() for name, value, error in zip(names, fit.x, sd)}
    metrics = []
    for split in ["train", "test"]:
        for observable in sorted({r["observable"] for r in observations}):
            mask = np.array([r["split"] == split and r["observable"] == observable for r in observations])
            if mask.any():
                metrics.append({"split": split, "observable": observable, "n": int(mask.sum()),
                                "rmse_mM": float(np.sqrt(np.mean((predictions[mask]-measured[mask])**2))),
                                "baseline_rmse_mM": float(np.sqrt(np.mean((baseline[mask]-measured[mask])**2))),
                                "standardized_rmse": float(np.sqrt(np.mean(((predictions[mask]-measured[mask])/sigma[mask])**2)))})
    for row, predicted, base in zip(observations, predictions, baseline):
        row.update(predicted_mM=float(predicted), baseline_mM=float(base))
    report = {"data_kind": study["data_kind"], "optimizer_success": bool(fit.success),
              "fitted_parameters": dict(zip(names, np.exp(fit.x).tolist())), "conditional_approximate_95pct_intervals": intervals,
              "jacobian_rank": rank, "jacobian_singular_values": singular_values.tolist(),
              "jacobian_condition_number": float(singular_values[0]/singular_values[-1]) if singular_values[-1] else None,
              "metrics": metrics, "provenance": study["provenance"],
              "note": "Intervals are local Gaussian approximations conditional on fixed chemistry, known independent errors, and model correctness. They exclude structural discrepancy, batch effects, uncertain pools, and transcript boundaries. No interval-coverage claim is made."}
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    write_rows(output/"predictions.csv", observations)
    (output/"calibration.json").write_text(json.dumps(report, indent=2)+"\n")
    return report
