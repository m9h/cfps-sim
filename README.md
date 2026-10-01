# cfps-sim

A first computational stack for lysate-based cell-free protein synthesis (CFPS):
recipe accounting, a reduced CUDA PyTorch model with SBML export, and a
BoTorch experiment-selection demonstration. Numerical execution requires a CUDA GPU; Modal provides a remote GPU runner.

**Status:** runnable research scaffold. The kinetic constants are illustrative,
the model is uncalibrated, and the BoTorch example uses synthetic observations.
It does not predict competition winners. Recipe validation and simulation are
separate workflows. Initial-pool mapping now requires explicit base chemistry and lysate backgrounds; a calibrated recipe-response model remains future work.

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
Simulation, calibration, and BoTorch acquisition require CUDA-enabled PyTorch
and a visible NVIDIA GPU. Recipe accounting, sequence inspection, and SBML
serialization can run on the development host. The numerical backend raises
an error when CUDA is unavailable; it does not silently change execution devices.

```sh
uv run cfps-sim validate examples/recipe.json
uv run cfps-sim simulate --output results/demo
uv run --extra bo cfps-sim bo-cost-demo
uv run pytest
```

## Pre-experiment evidence and readiness

The source-aware registry in `data/registry/inventory.draft.json` normalizes 66
cached stock entries and preserves 12 price conflicts with the draft ledger.
Neither source is designated the operative inventory. Mass/activity stocks,
mixed-stock constituents, base chemistry, availability and assay conversions
remain unresolved where source evidence is insufficient. Salt accounting reports
analytical totals, not free ions or pH. No substrate is silently converted into
the model's lumped Fuel.

```sh
uv run cfps-sim preflight examples/recipe.json --prices examples/prices.draft.json
uv run python scripts/prepare_review.py
uv run python scripts/fetch_benchmark.py
uv run cfps-sim benchmark                         # requires CUDA
uv run --extra cloud modal run modal_app.py --benchmark
```

`preflight` intentionally exits 2 while readiness gates are blocked.
`model_initial_conditions` requires explicit backgrounds for all supported NTP
and amino-acid pools, explicit other initial states, and resolved dispensed
chemical stocks. It returns missing-mechanism warnings with the mapping; this
is not a calibrated recipe response. `reports/pre_experiment_review.json`
contains prioritized diagnostic hypotheses and shortlist acceptance criteria,
with candidate rankings withheld pending in-domain evidence and reconciled costs.

The measured-data benchmark downloads ten pinned Borkowski author tables
(1,017 formulation means) into ignored `data/raw/borkowski2020`, verifies each
Git blob hash, and records SHA-256 provenance. Raw tables are not redistributed.
The benchmark uses source feature scales and relative yield, without an assumed
conversion to mass or competition reagents. Plates 1–7 train the model; plate 8
selects linear/quadratic ridge and regularization; plate 9 calibrates intervals;
plate 10 tests predictions. Shared recipes are purged from earlier splits, and
all preprocessing is fitted on training data. Final fitting uses training and
validation plates only.

The A100 result in `reports/borkowski2020.gpu.json` has MAE 0.581 versus 1.011
for a training-mean baseline, but Spearman correlation −0.315 and only 61.8%
coverage for nominal 90% intervals. This warns against using this baseline for
recipe ranking. Later active-learning plates are distribution-shifted; interval
coverage is an empirical check, not a guarantee. The benchmark neither reproduces
the authors' original policy nor validates Ginkgo predictions or cost savings.
Fetch the data before a remote benchmark or executable article build.

## Production economics

The competition objective is **minimize specific production cost**, not maximize
yield alone: `USD/g = reaction cost / (titer in g/L × reaction volume in L)`.
The ledger charges all dispensed stocks, including fixed lysate, base buffer,
DNA, and explicit water. Missing prices are reported and block a complete score.
It currently covers dispensed reagent costs, not labor, assays, purification,
plates, or compute. Confirm the official scoring scope and full price ledger.

```sh
uv run cfps-sim cost examples/recipe.json examples/prices.draft.json
uv run cfps-sim score examples/recipe.json examples/prices.illustrative.json \
  examples/yield.synthetic.json
uv run --extra bo cfps-sim bo-cost-demo
```

The first command intentionally exits with status 2: the partial cached draft
has no verified base-buffer price. The complete illustrative ledger assumes
that price; its yield file is synthetic. Neither example is an official score
or a demonstrated productive recipe. `score` accepts explicit, source-matched
`g/L` or `U/mL` observations and keeps `USD/g` separate from `USD/U`. Apparent
enzyme g/L requires an audited assay conversion; raw fluorescence or activity
must not be silently treated as mass. Zero product is unscorable.

`bo-cost-demo` models synthetic log titer and applies known toy reaction costs
to both observed and candidate posterior samples. It selects expected improvement
in `log(g/USD)` using noisy Monte Carlo acquisition. This accounts for production
cost within the objective; the original `bo-demo` remains a yield-surface smoke
test. Neither BO example uses experimental competition observations or an
implemented recipe-to-model mapping. Highest yield and lowest USD/g can select
different points. Report both raw yield and cost, retain the cost/yield Pareto
front, and rank competition candidates by their target-specific economic score.

## Executable article

`docs/article.Pnw` documents the model, assumptions, recipe constraints,
optimization, and GPU workflow. Its figure and numerical recipe summary are
generated from the package with Pweave; Tectonic builds the PDF.
The introduction also covers CFPS history, a curated modern reading guide,
and Joongoo Lee's POSTECH i-POPFLEX work presented in the Build-a-Cell series.
The concise abstract summarizes the implemented scope and remaining validation gaps.
The article also distinguishes measured CFPS formulation improvements from
Illinois iBioFAB's enzyme-activity and metabolic-production results, and proposes
assay calibration, module-swap experiments, and prospective comparisons of search
policies under equal experimental budgets. These experimental steps remain planned.

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

## GPU computation and numerical accuracy

The implemented backend uses PyTorch CUDA for the ODE solver, parameter fitting,
Gaussian processes, and acquisition evaluation. State tensors and GP calculations
use `float64`; `torch.compile` fuses per-step flux calculations on the GPU.
Adaptive Dormand–Prince integration uses `rtol=1e-10` and
`atol=1e-14`; the challenge suite repeats the endpoint calculation at
`rtol=1e-12`, `atol=1e-16`. These are numerical error controls, not guarantees
of predictive accuracy. Conservation, nonnegativity, missing-resource controls,
and tolerance refinement must pass on the target GPU.

Calibration uses bounded damped Gauss–Newton in log-parameter space. Residuals,
finite-difference Jacobians, and linear solves stay on CUDA. Two starts check
local sensitivity; Jacobians are compared at two perturbation sizes, and
conditional covariance uses an SVD. Convergence uses scaled gradient, objective,
or step criteria; Jacobian disagreement above 0.1% rejects fit success. Artifacts record GPU identity, PyTorch/CUDA
versions, solver tolerances, and fit diagnostics. Host code handles inputs,
reporting, and serialization.

The BO demo fits twelve synthetic observations and selects from a normalized
grid. This demonstrates acquisition only; it is not a dispensing-feasible
competition candidate set. Better biological accuracy requires measured time
courses, an audited assay mapping, identifiable parameters, and validation on
whole held-out recipes and lysate batches. GPU execution alone does not supply
that evidence.

[PyTorch numerical accuracy](https://docs.pytorch.org/docs/stable/notes/numerical_accuracy.html)
and [torchdiffeq](https://github.com/rtqichen/torchdiffeq) describe the numerical
backend. A JAX GPU/TPU implementation would be a separate future backend; none
is implemented or validated here.

## Run on Modal GPU

```sh
uv run --extra cloud modal setup
uv run --extra cloud modal run modal_app.py
uv run --extra cloud modal run modal_app.py --construct constructs/competition/sfgfp.json
# Run GPU tests, challenge/calibration checks, and rebuild the executable article:
uv run --extra cloud modal run --detach modal_app.py --verify --build-docs --background
# Retrieve the printed function-call ID once complete:
uv run --extra cloud modal run modal_app.py --retrieve fc-...
```

The runner requests an A100 GPU and builds the locked Linux environment.
Results go to `results/modal/`; `--build-docs` also replaces `docs/article.pdf`
with the GPU-generated article. The on-demand function has a one-hour timeout.
[Modal GPU compute](https://modal.com/docs/guide/gpu) is billed to the authenticated
account. A CUDA-enabled RunPod host can run the same locked package and commands;
no provider-specific RunPod integration is included.

All 18 tests passed on a Modal A100 GPU on September 30, 2026, including
conservation, solver refinement, parameter recovery, holdout isolation, and
acquisition checks. These validate software behavior; measured biological
accuracy remains unestablished. Hosted GitHub CI runs host checks and skips GPU
tests; use the GPU verification command for the complete suite.

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

## Computational experiments and calibration

```sh
uv run cfps-sim experiments --output results/experiments --seed 2026
uv run cfps-sim calibrate results/experiments/synthetic/study.json \
  results/experiments/synthetic/observations.csv --output results/refit
```

This runs resource/process knockouts, conservation and solver checks, local
sensitivity, ATP/GTP and transcription/translation grids, transcript-boundary
scenarios, and a 64-member uncertainty exercise. The synthetic calibration fits
two effective parameters to three conditions and predicts a withheld condition;
it does not fit real competition observations. Outputs retain inputs, provenance,
units, seeds, predictions, and diagnostic statistics.
[Experiment instructions](experiments/README.md) document the data schema and
limitations. The article includes an executable worked example and a published
calibration precedent from Horvath et al.

[Research priorities](examples/research_priorities.json) distinguish draft-catalog
ingredients from proposed requests and mechanisms still absent from the model.
The article also covers T7 resource allocation, chemical/kinetic databases, and
lessons from Ben Blaiszik's Argonne/Chicago data ecosystem and the October
21–22, 2026 [Open Scientific Intelligence Hackathon](https://llmhackathon.github.io/).

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

- Borkowski et al., *Large scale active-learning-guided exploration for in vitro
  protein production optimization* (2020). Reports up to 34-fold improvement in
  relative sfGFP yield over the tested homemade lysate's initial buffer.
  [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7170859/),
  [DOI](https://doi.org/10.1038/s41467-020-15798-5).
- HamediRad et al., *Towards a fully automated algorithm driven platform for
  biosystems design* (2019). Illinois BioAutomata/iBioFAB; cellular lycopene
  optimization compared with random screening.
  [DOI](https://doi.org/10.1038/s41467-019-13189-z).
- Singh et al., *A generalized platform for artificial intelligence-powered
  autonomous enzyme engineering* (2025). Illinois iBioFAB; enzyme activity
  improvements rather than CFPS protein-mass yield.
  [DOI](https://doi.org/10.1038/s41467-025-61209-y).

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
