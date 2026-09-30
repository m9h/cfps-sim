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
    else:
        from .bo import demo
        print(json.dumps(demo(), indent=2))
