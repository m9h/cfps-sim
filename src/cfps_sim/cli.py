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
    commands.add_parser("bo-demo")
    args = parser.parse_args()
    if args.command == "validate":
        from .recipe import validate
        print(json.dumps(validate(json.loads(args.recipe.read_text())), indent=2))
    elif args.command == "simulate":
        import numpy as np
        from .model import simulate
        rr, values = simulate(args.minutes)
        args.output.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.output / "trajectory.csv", values, delimiter=",", header=",".join(values.colnames), comments="")
        (args.output / "model.xml").write_text(rr.getSBML())
        print(f"Saved uncalibrated demonstration to {args.output}")
    else:
        from .bo import demo
        print(json.dumps(demo(), indent=2))
