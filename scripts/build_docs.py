"""Weave executable documentation, then compile with Tectonic."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"

if __name__ == "__main__":
    if not shutil.which("tectonic"):
        raise SystemExit("Install Tectonic first: brew install tectonic")
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    # Make the installed project interpreter visible to Pweave's kernel launcher.
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env["PATH"]
    subprocess.run(["pweave", "-f", "tex", "article.Pnw"], cwd=DOCS, env=env, check=True)
    if "Traceback" in (DOCS / "article.tex").read_text():
        raise SystemExit("A documentation code chunk failed; inspect docs/article.tex")
    for name in ["intro-resources.pdf", "intro-workflow.pdf", "trajectory.pdf", "computational-study.pdf"]:
        if not (DOCS / "figures" / name).is_file():
            raise SystemExit(f"Pweave did not generate {name}")
    subprocess.run(["tectonic", "article.tex"], cwd=DOCS, env=env, check=True)
    print(f"Built {DOCS / 'article.pdf'}")
