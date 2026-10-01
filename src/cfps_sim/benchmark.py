"""GPU retrospective baseline on measured CFPS data with disjoint recipes."""
import csv
import hashlib
import json
import math
from pathlib import Path

from .backend import cuda_device, metadata

FEATURES = ["nad", "folinic_acid", "coa", "nucleo_mix", "spermidin", "pga",
            "aa", "trna", "mg_gluta", "camp", "K_gluta"]


def read_rows(path):
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key in FEATURES + ["yield", "yield_std"]:
            row[key] = float(row[key])
            if not math.isfinite(row[key]) or row[key] < 0:
                raise ValueError(f"Invalid {key} in measured table")
        row["plate"] = int(row["plate"])
        if row["plate"] not in range(1, 11):
            raise ValueError("Benchmark requires plate identifiers 1 through 10")
        row["recipe_key"] = tuple(row[key] for key in FEATURES)
    return rows


def partition(rows):
    """Later blocks own shared recipes, preventing replication across splits."""
    groups = {"train": [], "validation": [], "calibration": [], "test": []}
    seen, purged = set(), {}
    for name, plates in [("test", {10}), ("calibration", {9}), ("validation", {8}),
                         ("train", set(range(1, 8)))]:
        candidates = [row for row in rows if row["plate"] in plates]
        groups[name] = [row for row in candidates if row["recipe_key"] not in seen]
        purged[name] = len(candidates) - len(groups[name])
        seen.update(row["recipe_key"] for row in groups[name])
        if not groups[name]:
            raise ValueError(f"No independent {name} recipes remain after purging")
    return groups, purged


def run_benchmark(observations, provenance, output):
    import torch
    device = cuda_device()
    rows = read_rows(observations)
    groups, purged = partition(rows)
    source = json.loads(Path(provenance).read_text())
    checksum = hashlib.sha256(Path(observations).read_bytes()).hexdigest()
    if checksum != source["normalized_sha256"]:
        raise ValueError("Normalized measurement checksum mismatch")
    def tensor(block):
        return (torch.tensor([[row[key] for key in FEATURES] for row in block], device=device, dtype=torch.float64),
                torch.tensor([row["yield"] for row in block], device=device, dtype=torch.float64))
    data = {name: tensor(block) for name, block in groups.items()}
    def fit_predict(training, destinations, alpha, quadratic):
        x, y = training
        mean, scale = x.mean(0), x.std(0, correction=0).clamp_min(1e-12)
        def expand(raw):
            z = (raw-mean)/scale
            if quadratic:
                pair = torch.triu_indices(z.shape[1], z.shape[1], device=device)
                z = torch.cat([z, z[:, pair[0]] * z[:, pair[1]]], dim=1)
            return z
        features = expand(x)
        fm, fs = features.mean(0), features.std(0, correction=0).clamp_min(1e-12)
        def design(raw):
            normalized = (expand(raw)-fm)/fs
            return torch.cat([torch.ones((len(raw), 1), device=device, dtype=raw.dtype), normalized], dim=1)
        matrix = design(x)
        penalty = torch.eye(matrix.shape[1], device=device, dtype=x.dtype)*alpha
        penalty[0, 0] = 0
        weights = torch.linalg.solve(matrix.T@matrix+penalty, matrix.T@y)
        return {name: (design(raw)@weights).clamp_min(0) for name, raw in destinations.items()}
    candidates = []
    for quadratic in [False, True]:
        for alpha in [0.1, 1., 10., 100.]:
            prediction = fit_predict(data["train"], {"validation": data["validation"][0]}, alpha, quadratic)["validation"]
            mae = float((prediction-data["validation"][1]).abs().mean())
            candidates.append({"model": "quadratic_ridge" if quadratic else "linear_ridge", "alpha": alpha, "validation_mae": mae})
    selected = min(candidates, key=lambda record: record["validation_mae"])
    training = tuple(torch.cat([data["train"][i], data["validation"][i]]) for i in [0, 1])
    predictions = fit_predict(training, {name: data[name][0] for name in ["calibration", "test"]},
                              selected["alpha"], selected["model"] == "quadratic_ridge")
    residuals = (predictions["calibration"]-data["calibration"][1]).abs().sort().values
    coverage_target = .9
    order = math.ceil((len(residuals)+1)*coverage_target)
    if order > len(residuals):
        raise ValueError("Too few calibration observations for finite conformal interval")
    radius = float(residuals[order-1])
    observed, predicted = data["test"][1], predictions["test"]
    lower, upper = (predicted-radius).clamp_min(0), predicted+radius
    baseline = training[1].mean().expand_as(observed)
    def ranks(values):
        unique, inverse, counts = torch.unique(values, sorted=True, return_inverse=True, return_counts=True)
        return (counts.cumsum(0)-counts + (counts-1)/2)[inverse].to(torch.float64)
    rank_observed, rank_predicted = ranks(observed), ranks(predicted)
    correlation = torch.corrcoef(torch.stack([rank_observed, rank_predicted]))[0, 1]
    metrics = {"mae": float((predicted-observed).abs().mean()),
               "rmse": float(((predicted-observed)**2).mean().sqrt()),
               "mean_baseline_mae": float((baseline-observed).abs().mean()),
               "spearman": float(correlation) if torch.isfinite(correlation) else None,
               "interval_coverage": float(((observed >= lower) & (observed <= upper)).double().mean()),
               "mean_interval_width": float((upper-lower).mean()), "interval_radius": radius}
    report = {"data_kind": "experimental", "execution": metadata(), "source": source["manifest"],
              "normalized_sha256": checksum, "feature_names": FEATURES,
              "split": {name: {"rows": len(block), "unique_recipes": len({r['recipe_key'] for r in block}),
                               "plates": sorted({r['plate'] for r in block}), "purged_rows": purged[name]}
                        for name, block in groups.items()},
              "recipe_overlap_across_splits": 0, "selection": selected, "validation_candidates": candidates,
              "metrics": metrics, "interval_target": coverage_target,
              "limitations": source["manifest"]["limitations"] + [
                  "chronological plate holdout is distribution-shifted, not an independent lysate-batch test",
                  "split-conformal coverage assumes exchangeability; later-plate coverage is empirical only",
                  "intervals concern processed formulation means, not individual wells",
                  "negative regression outputs clipped to zero; no assay-to-mass conversion",
                  "no cost objective, mechanistic comparison, or Ginkgo performance claim"]}
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    (out/"benchmark.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    with (out/"heldout_predictions.csv").open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["plate", "source_row", "observed_relative_yield", "predicted_relative_yield", "lower", "upper"])
        for row, actual, prediction, lo, hi in zip(groups["test"], observed.tolist(), predicted.tolist(), lower.tolist(), upper.tolist()):
            writer.writerow([row["plate"], row["source_row"], actual, prediction, lo, hi])
    return report
