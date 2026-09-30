"""Run the same package on Modal. Launch from the repository root."""
from pathlib import Path
import modal

ROOT = Path(__file__).parent.resolve()
app = modal.App("cfps-sim")
image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("uv")
    .add_local_file(ROOT / "pyproject.toml", "/opt/cfps-sim/pyproject.toml", copy=True)
    .add_local_file(ROOT / "uv.lock", "/opt/cfps-sim/uv.lock", copy=True)
    .add_local_dir(ROOT / "src", "/opt/cfps-sim/src", copy=True)
    .run_commands("cd /opt/cfps-sim && uv sync --frozen --extra bo --no-dev")
    .env({"PATH": "/opt/cfps-sim/.venv/bin:/usr/local/bin:/usr/bin:/bin"})
)


@app.function(image=image, cpu=2, memory=4096, timeout=600)
def run_demo():
    import sys
    # Modal's runtime interpreter loads packages from the locked project environment.
    sys.path.insert(0, "/opt/cfps-sim/.venv/lib/python3.11/site-packages")
    sys.path.insert(0, "/opt/cfps-sim/src")
    from cfps_sim.model import simulate
    from cfps_sim.bo import demo
    rr, values = simulate()
    return {"bo": demo(), "model": "uncalibrated demonstration", "final_gfp_mM": float(values[-1, 9]),
            "columns": list(values.colnames), "trajectory": values.tolist(), "sbml": rr.getSBML()}


@app.local_entrypoint()
def main():
    import json
    result = run_demo.remote()
    out = ROOT / "results/modal"
    out.mkdir(parents=True, exist_ok=True)
    (out / "demo.json").write_text(json.dumps(result, indent=2))
    (out / "model.xml").write_text(result["sbml"])
    print(json.dumps({"bo": result["bo"], "final_gfp_mM": result["final_gfp_mM"], "output": str(out)}, indent=2))
