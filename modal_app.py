"""Run the same package on Modal. Launch from the repository root."""
from pathlib import Path
import modal

ROOT = Path(__file__).parent.resolve()
app = modal.App("cfps-sim")
image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("curl", "ca-certificates", "libfontconfig1", "fonts-dejavu-core",
                 "libgraphite2-3", "libharfbuzz0b", "libicu72")
    .pip_install("uv")
    .run_commands("curl -L --fail https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%400.15.0/tectonic-0.15.0-x86_64-unknown-linux-gnu.tar.gz -o /tmp/tectonic.tar.gz", "tar -xzf /tmp/tectonic.tar.gz -C /usr/local/bin tectonic")
    .run_commands("tectonic --version")
    .add_local_file(ROOT / "pyproject.toml", "/opt/cfps-sim/pyproject.toml", copy=True)
    .add_local_file(ROOT / "uv.lock", "/opt/cfps-sim/uv.lock", copy=True)
    .run_commands("cd /opt/cfps-sim && uv sync --frozen --extra bo --extra docs --no-install-project")
    .add_local_dir(ROOT / "src", "/opt/cfps-sim/src", copy=True)
    .add_local_dir(ROOT / "tests", "/opt/cfps-sim/tests", copy=True)
    .add_local_dir(ROOT / "constructs", "/opt/cfps-sim/constructs", copy=True)
    .add_local_dir(ROOT / "examples", "/opt/cfps-sim/examples", copy=True)
    .add_local_dir(ROOT / "docs", "/opt/cfps-sim/docs", copy=True)
    .add_local_dir(ROOT / "scripts", "/opt/cfps-sim/scripts", copy=True)
    .run_commands("cd /opt/cfps-sim && uv sync --frozen --extra bo --extra docs")
    .env({"PATH": "/opt/cfps-sim/.venv/bin:/usr/local/bin:/usr/bin:/bin"})
)


@app.function(image=image, gpu="A100", memory=8192, timeout=3600)
def run_demo(construct=None, verify: bool = False, build_docs: bool = False):
    import sys
    # Modal's runtime interpreter loads packages from the locked project environment.
    sys.path.insert(0, "/opt/cfps-sim/.venv/lib/python3.11/site-packages")
    sys.path.insert(0, "/opt/cfps-sim/src")
    from cfps_sim.model import simulate
    from cfps_sim.bo import demo, cost_demo
    from cfps_sim.construct import demands
    if verify or build_docs:
        from cfps_sim.backend import cuda_device
        cuda_device()
        import subprocess
        checks = {}
        if verify:
            subprocess.run(["/opt/cfps-sim/.venv/bin/python", "-m", "pytest", "-q"], cwd="/opt/cfps-sim", check=True)
        if verify and not build_docs:
            from cfps_sim.experiments import run_suite, make_synthetic_study, calibrate
            output = Path("/opt/cfps-sim/results/verification")
            checks["challenges"] = run_suite("/opt/cfps-sim/constructs/competition/sfgfp.json", output)
            study = make_synthetic_study("/opt/cfps-sim/constructs/competition/sfgfp.json", output/"synthetic")
            checks["calibration"] = calibrate(study, study.parent/"observations.csv", output/"calibration")
        if build_docs:
            subprocess.run(["/opt/cfps-sim/.venv/bin/python", "scripts/build_docs.py"], cwd="/opt/cfps-sim", check=True)
        if build_docs:
            import json
            output = Path("/opt/cfps-sim/results/article-experiments")
            checks["challenges"] = json.loads((output/"report.json").read_text())
            checks["calibration"] = json.loads((output/"calibration/calibration.json").read_text())
        if not checks["challenges"]["all_checks_pass"] or not checks["calibration"]["optimizer_success"]:
            raise RuntimeError("GPU verification failed")
        return {"verification": checks, "article_pdf": (Path("/opt/cfps-sim/docs/article.pdf").read_bytes() if build_docs else None)}
    rr, values = simulate(construct=construct)
    return {"bo": demo(), "bo_cost": cost_demo(), "model": "uncalibrated demonstration", "final_gfp_mM": float(values[-1, 9]),
            "execution": rr.execution_metadata(), "construct": demands(construct) if construct else None,
            "columns": list(values.colnames), "trajectory": values.tolist(), "sbml": rr.getSBML()}


@app.local_entrypoint()
def main(construct: str = "", verify: bool = False, build_docs: bool = False,
         background: bool = False, retrieve: str = ""):
    import json
    from cfps_sim.construct import read_construct
    artifact = read_construct(construct) if construct else None
    out = ROOT / "results/modal"
    out.mkdir(parents=True, exist_ok=True)
    if retrieve:
        try:
            result = modal.FunctionCall.from_id(retrieve).get(timeout=0)
        except TimeoutError:
            print("GPU job is still running")
            return
    elif background:
        call = run_demo.spawn(artifact, verify=verify, build_docs=build_docs)
        (out / "function_call.json").write_text(json.dumps({"id": call.object_id})+"\n")
        print(f"Started GPU job {call.object_id}; retrieve with --retrieve {call.object_id}")
        return
    else:
        result = run_demo.remote(artifact, verify=verify, build_docs=build_docs)
    if "verification" in result:
        pdf = result.pop("article_pdf")
        if pdf is not None:
            (ROOT / "docs/article.pdf").write_bytes(pdf)
        (out / "verification.json").write_text(json.dumps(result, indent=2)+"\n")
        print(json.dumps({"all_challenge_checks_pass": result["verification"]["challenges"]["all_checks_pass"],
                          "optimizer_success": result["verification"]["calibration"]["optimizer_success"],
                          "article_rebuilt": pdf is not None, "output": str(out)}, indent=2))
        return
    (out / "demo.json").write_text(json.dumps(result, indent=2))
    (out / "model.xml").write_text(result["sbml"])
    print(json.dumps({"bo": result["bo"], "bo_cost": result["bo_cost"], "final_gfp_mM": result["final_gfp_mM"], "output": str(out)}, indent=2))
