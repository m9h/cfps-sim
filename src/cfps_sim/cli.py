import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="CFPS starter simulator")
    commands = parser.add_subparsers(dest="command", required=True)
    validation = commands.add_parser("validate")
    validation.add_argument("recipe", type=Path)
    simulation = commands.add_parser("simulate")
    simulation.add_argument("--output", type=Path, default=Path("results/demo"))
    simulation.add_argument("--minutes", type=float, default=240)
    simulation.add_argument("--construct", type=Path, help="Versioned sfGFP construct JSON")
    construct = commands.add_parser("construct", help="Inspect or exchange sequence artifacts")
    actions = construct.add_subparsers(dest="construct_action", required=True)
    inspect = actions.add_parser("inspect")
    inspect.add_argument("input", type=Path)
    inspect.add_argument("--aa-pools", type=Path, help="JSON mapping one-letter amino acids to available mM")
    export = actions.add_parser("export-genbank")
    export.add_argument("input", type=Path)
    export.add_argument("--output", type=Path, required=True)
    imported = actions.add_parser("import-genbank")
    imported.add_argument("input", type=Path)
    imported.add_argument("--output", type=Path, required=True)
    imported.add_argument("--id")
    imported.add_argument("--target")
    imported.add_argument("--transcript-start", type=int)
    imported.add_argument("--transcript-end", type=int)
    commands.add_parser("bo-demo")
    experiments = commands.add_parser("experiments", help="Run model challenges or a synthetic calibration study")
    experiments.add_argument("--construct", type=Path, default=Path("constructs/competition/sfgfp.json"))
    experiments.add_argument("--output", type=Path, default=Path("results/experiments"))
    experiments.add_argument("--seed", type=int, default=2026)
    calibration = commands.add_parser("calibrate", help="Fit supported concentration observations with declared study provenance")
    calibration.add_argument("study", type=Path)
    calibration.add_argument("observations", type=Path)
    calibration.add_argument("--output", type=Path, default=Path("results/calibration"))
    args = parser.parse_args()
    if args.command == "validate":
        from .recipe import validate
        print(json.dumps(validate(json.loads(args.recipe.read_text())), indent=2))
    elif args.command == "simulate":
        import numpy as np
        from .model import simulate
        from .construct import read_construct, demands
        construct = read_construct(args.construct) if args.construct else None
        rr, values = simulate(args.minutes, construct=construct)
        args.output.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.output / "trajectory.csv", values, delimiter=",", header=",".join(values.colnames), comments="")
        (args.output / "model.xml").write_text(rr.getSBML())
        metadata = {"model_status": "uncalibrated", "minutes": args.minutes,
                    "construct": demands(construct) if construct else None,
                    "initial_aa_pools": "2 mM each; illustrative" if construct else "40 mM pooled; illustrative"}
        (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
        print(f"Saved uncalibrated demonstration to {args.output}")
    elif args.command == "construct":
        from .construct import read_construct, demands, export_genbank, import_genbank
        if args.construct_action == "inspect":
            pools = json.loads(args.aa_pools.read_text()) if args.aa_pools else None
            print(json.dumps(demands(read_construct(args.input), pools), indent=2))
        elif args.construct_action == "export-genbank":
            export_genbank(read_construct(args.input), args.output)
        else:
            record = import_genbank(args.input, args.id, args.target, args.transcript_start, args.transcript_end)
            args.output.write_text(json.dumps(record, indent=2)+"\n")
    elif args.command == "experiments":
        from .experiments import run_suite, make_synthetic_study, calibrate
        checks = run_suite(args.construct, args.output, args.seed)
        study = make_synthetic_study(args.construct, args.output/"synthetic", args.seed)
        fitted = calibrate(study, study.parent/"observations.csv", args.output/"calibration")
        print(json.dumps({"output": str(args.output), "checks": checks["checks"], "calibration_data": fitted["data_kind"],
                          "fitted_parameters": fitted["fitted_parameters"], "metrics": fitted["metrics"]}, indent=2))
        if not checks["all_checks_pass"] or not fitted["optimizer_success"]:
            raise SystemExit("Computational challenge or fitting check failed")
    elif args.command == "calibrate":
        from .experiments import calibrate
        report = calibrate(args.study, args.observations, args.output)
        print(json.dumps(report, indent=2))
        if not report["optimizer_success"]:
            raise SystemExit("Calibration optimizer did not converge")
    else:
        from .bo import demo
        print(json.dumps(demo(), indent=2))
