# cfps-sim

A first computational stack for lysate-based cell-free protein synthesis (CFPS):
recipe accounting, a reduced Antimony/RoadRunner model with SBML export, and a
BoTorch experiment-selection demonstration. Runs locally or on Modal.

**Status:** runnable research scaffold. The kinetic constants are illustrative,
the model is uncalibrated, and the BoTorch example uses synthetic observations.
It does not predict competition winners. Recipe validation and simulation are
currently separate commands; recipe-to-model mapping is a planned milestone.

## Install with uv

On macOS, install uv with Homebrew if needed:

```sh
brew install uv
git clone https://github.com/m9h/cfps-sim.git
cd cfps-sim
uv sync --frozen --extra bo --extra cloud
```

uv uses the Python 3.11 version specified in `.python-version`; it can download
that interpreter if necessary. `uv.lock` pins dependencies. Native wheels are
used where available; no global Python packages are needed. COPASI is an optional
GUI for inspecting the exported SBML, not a runtime requirement.
On Linux, RoadRunner's wheel requires the shared Python 3.11 library. The CI
workflow provisions Python with `actions/setup-python` before `uv sync`; Modal's
system Python already supplies it. If a standalone interpreter reports missing
`libpython3.11.so.1.0`, use a Python installation with its shared library available
to the dynamic loader, then recreate the environment with `uv sync --python /path/to/python3.11`.

```sh
uv run cfps-sim validate examples/recipe.json
uv run cfps-sim simulate --output results/demo
uv run --extra bo cfps-sim bo-demo
uv run pytest
```

## Executable article

`docs/article.Pnw` documents the model, assumptions, recipe constraints,
optimization, and Modal workflow. Its figure and numerical recipe summary are
generated from the package with Pweave; Tectonic builds the PDF.
The introduction also covers CFPS history, a curated modern reading guide,
and Joongoo Lee's POSTECH i-POPFLEX work presented in the Build-a-Cell series.

```sh
brew install tectonic
uv sync --frozen --extra bo --extra cloud --extra docs
uv run --extra docs python scripts/build_docs.py
```

The output is `docs/article.pdf`. The initial build may download TeX support
files. Generated TeX and figures are ignored; the article source and PDF are tracked.
`AGENTS.md` requires updating the article and rebuilding its PDF when scientific
behavior, assumptions, workflows, or validation change.

The simulator writes `trajectory.csv` and `model.xml`. Time is in minutes and
concentrations in mM. `Protein` is immature protein and `GFP` is mature protein;
neither is calibrated fluorescence. Without `--construct`, sequence counts remain
illustrative placeholders. With a supplied sfGFP construct, the model uses its
CDS-derived amino-acid counts and assumed-transcript NTP counts.

## DNA constructs and Popcorn exchange

The three competition plasmids are saved as versioned JSON and GenBank artifacts
in `constructs/competition/`, with source provenance and sequence checksums. Their
CDS translations match the proteins in the draft rules, including tags.
[Construct documentation](constructs/README.md) specifies the schema and the
assumed transcript boundaries. These are supplied templates; custom DNA is not
an established competition submission variable.

```sh
uv run cfps-sim construct inspect constructs/competition/sfgfp.json
uv run cfps-sim simulate --construct constructs/competition/sfgfp.json --output results/sfgfp
uv run cfps-sim construct export-genbank constructs/competition/sfgfp.json --output results/sfgfp.gb
uv run cfps-sim construct import-genbank results/sfgfp.gb --output results/sfgfp.json
```

Keep `m9h/popcorn` in its own environment and exchange GenBank artifacts. Import
requires a complete expression construct with one CDS and an mRNA annotation,
or explicit zero-based `--transcript-start` and `--transcript-end` coordinates.
The existing Popcorn example's CDS-only output and random-sequence fallback
are not verified competition constructs. No sequence-score-to-rate calibration
or DNA optimization is implemented here.

`construct inspect --aa-pools pools.json` accepts a map of canonical one-letter
amino-acid codes to available mM and calculates a static stoichiometric ceiling.
It excludes synthesis, recycling, and unknown lysate pools. Inspection works for
all three targets; simulation currently accepts sfGFP only because enzyme folding
and activity observation models are still required. Simulation also writes
`metadata.json` with sequence demands and explicit model assumptions.

## Apple Silicon and BoTorch

BoTorch runs on Apple Silicon CPUs with PyTorch `torch.double` tensors. This
project explicitly uses CPU float64. BoTorch recommends double precision for
Gaussian-process stability, while PyTorch's MPS/Metal backend cannot allocate
float64 tensors. GPU acceleration is therefore not the default on a Mac.

- [BoTorch double-precision guidance](https://github.com/meta-pytorch/botorch/discussions/1444)
- [PyTorch MPS limitations](https://docs.pytorch.org/docs/main/notes/mps.html)

The demo fits a GP to 12 synthetic observations and chooses a candidate using
log expected improvement from a discrete normalized grid. This grid demonstrates
selection only; it is not a dispensing-feasible competition candidate set.

## Run on Modal

```sh
uv run --extra cloud modal setup
uv run --extra cloud modal run modal_app.py
uv run --extra cloud modal run modal_app.py --construct constructs/competition/sfgfp.json
```

The first command authenticates your account. The second builds a Linux Python
3.11 image using the same lockfile, runs the simulator and BoTorch remotely on
2 CPU cores with 4 GiB RAM, and saves results under `results/modal/`. This is an
on-demand run with a 10-minute timeout. It does not require a persistent deployment.
Modal compute is billed to your account.

For larger workloads, add parameter-ensemble or campaign functions and distribute
independent runs. Consider CUDA GPU execution only after benchmarking GP fitting
and acquisition workloads; small datasets often fit well on CPU.

## Recipe accounting

The example uses saved competition constraints: 20,000 nL total, 25 nL increments,
5,000 nL lysate, 2,000 nL base buffer, and one 2,000 nL DNA template. Reconfirm
these rules before submission at [the competition](https://coopetition.htgaa.org/cfps).

`stocks_mM` defines molar stocks explicitly. `base_final_mM` contains already
diluted base contributions, keyed by the same chemical identifiers. The result
excludes endogenous lysate pools. Mass concentrations, mixed amino-acid stocks,
polymers, activity units, counterions, and volume fractions need explicit mappings;
they must not be silently interpreted as molar stocks. Unknown reagents are listed
in the output rather than inferred.

## Model boundaries

The reduced model includes individual NTP consumption in transcription, mRNA
decay, amino-acid and ATP/GTP use in translation, fuel-dependent ATP regeneration,
machinery decay, and GFP maturation. The default demo pools amino acids; the
construct-aware model tracks all 20 individually and uses the smallest required
amino-acid saturation as an illustrative translation limitation. Initial pools
are 2 mM each, not measured lysate concentrations. Charging ATP use is
lumped into ADP equivalents; AMP, pyrophosphate, phosphate, recycling, magnesium,
pH, oxygen transfer, and detailed lysate metabolism are not represented yet.

ATP regeneration does not regenerate GTP, CTP, or UTP. These choices make resource
limits explicit but prevent detailed biochemical interpretation of trajectories.
Parameters and default placeholder sequence composition live in `src/cfps_sim/model.py`.
Construct-aware simulations replace sequence demands, while retaining uncalibrated
kinetics. Their SBML exports retain the initial conditions and actual sequence counts.

## Implementation milestones

1. Audit well-level recipes, batch/replicate identities, target sequences, assay
   timing, stock semantics, base buffer, and platform constraints.
2. Reconcile reference recipe calculations and map supported chemistry to model
   initial conditions; document unknown lysate contributions separately.
3. Refine the supplied constructs' assumed transcript boundaries and source
   kinetic parameters with uncertainty ranges. Select an initial energy pathway
   based on available data.
4. Fit identifiable effective parameters and compare against an empirical baseline
   using held-out recipes and batches, with uncertainty coverage checks.
5. Build quantized feasible candidate sets, mechanistic residual models, and
   batch acquisition with controls and replicates; compare campaign strategies.
6. Add separate PETase/reteplase folding and assay observation models.

## Papers

The competition also builds on the OpenAI–Ginkgo autonomous-lab study. The
article introduction connects that work to the Fall 2026 HTBAA course and the
current multi-target competition.

- Smith et al., *Using a GPT-5-driven autonomous lab to optimize the cost and
  titer of cell-free protein synthesis* (2026).
  [bioRxiv v1](https://www.biorxiv.org/content/10.64898/2026.02.05.703998v1),
  [official PDF](https://cdn.openai.com/pdf/5a12a3bc-96b7-4e07-9386-db6ee5bb2ed9/using-a-gpt-5-driven-autonomous-lab-to-optimize-the-cost-and-titer-of-cell-free-protein-synthesis.pdf).

- Horvath et al., *Toward a genome scale sequence specific dynamic model of
  cell-free protein synthesis in Escherichia coli* (2020).
  [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7136494/),
  [DOI](https://doi.org/10.1016/j.mec.2019.e00113).
- Chae et al., *Automated, modular assembly of reconstituted cell-free systems
  from in vitro-produced components* (2026), Joongoo Lee's group.
  [DOI](https://doi.org/10.1016/j.tibtech.2026.05.002).
- Nishimori et al., *Uncertainty-Driven Experiment Design in Cell-Free Protein
  Synthesis with Bayesian Optimization* (2026).
  [DOI](https://doi.org/10.1021/acssynbio.6c00179).

This repository contains citation links; publisher PDFs and conversation exports
are not bundled. Local downloads can be stored in the ignored `papers/downloads`
directory. Paper review should extract equations, units, parameter ranges,
data, and code availability before incorporating mechanisms.
