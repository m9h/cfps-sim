# Computational challenges and calibration

Run from the repository root on a CUDA GPU host:

```sh
uv run cfps-sim experiments --output results/experiments --seed 2026
uv run cfps-sim calibrate results/experiments/synthetic/study.json \
  results/experiments/synthetic/observations.csv --output results/refit
```

The suite uses the supplied sfGFP construct. It produces:

| Experiment | Output | Question |
|---|---|---|
| Reference plus nine resource/process interventions | `challenges.csv`, individual trajectories | Do missing required resources prevent expression? |
| Individual amino-acid balances and solver refinement | `report.json` | Is the computation finite, nonnegative, conservative, and converged? |
| Six kinetic parameters perturbed by ±10% | `sensitivity.csv` | Which local changes affect mature GFP at 240 minutes? |
| 4 × 4 initial ATP/GTP grid | `atp_gtp_grid.csv` | How do shared transcription/translation resources interact? |
| 5 × 5 transcription/translation capacity grid | `tx_tl_grid.csv` | Does greater effective transcription always improve expression? |
| Supplied sequence end versus assumed CDS end | `report.json` | How much do uncertain transcript bounds affect resource predictions? |
| 64 perturbed kinetic/pool scenarios | `ensemble.csv` | How broad is scenario spread under declared uncertainty assumptions? |
| Two-parameter synthetic recovery with one withheld condition | `synthetic/`, `calibration/` | Can fitting recover generating parameters without using withheld measurements? |

Pool interventions are not dispensing-feasible recipes. The ensemble assumes
independent lognormal perturbations with log-SD 0.2 for six rates and 0.02 for
five starting pools; these are illustrative assumptions, not measured error or
calibrated prediction intervals. Local elasticities and response grids are
properties of this reduced model, not evidence about which real reagent wins.
`no_fuel` can still express protein because initial NTPs remain available.

## Runnable calibration example

The demonstration generates GFP and mRNA time courses for three training
conditions and one withheld condition, each with eight times and two replicates.
The generating values are `k_tl = 0.26 /min` and `k_decay = 0.003 /min`.
Independent Gaussian measurement noise has known SD equal to 5% of the maximum
of each condition/channel. All other chemistry and measurement scaling are fixed.
This is intentionally a same-model recovery exercise, not experimental calibration.

The fitter uses GPU bounded damped Gauss–Newton least squares in log-parameter space, with residuals
scaled by the supplied measurement SD. Two initializations check basic local
optimizer sensitivity. Residuals, Jacobians, and linear solves execute on CUDA.
It compares Jacobians at two perturbation sizes and uses SVD for conditional
covariance. It reports per-observable train/test errors, Jacobian
singular values and rank, and approximate conditional 95% parameter intervals
when full rank. Those intervals exclude uncertain pools, model discrepancy,
correlated errors, batch effects, and transcript uncertainty. They do not replace
an ensemble/posterior or coverage assessment. Test measurements never enter the
fitting objective; the tests change them to verify that fitted values stay fixed.

For experimental inputs, preserve the generated schema:

- `study.json`: schema version 1; `data_kind: experimental`; full construct;
  unique condition IDs, `initial_mM`, and whole-condition `split` (`train`/`test`);
  positive ordered bounds for supported `fit_parameters`; `provenance.source`.
- CSV columns: `condition,split,replicate,time_min,observable,value_mM,sigma_mM`.
  Supported observables: GFP, mRNA, ATP, GTP, CTP, UTP concentrations. Times and
  measurement SD must be positive. Convert time and concentration units before
  import and retain original data and conversion provenance.

Raw fluorescence, activity units, and apparent enzyme g/L cannot be inserted as
these concentrations without a validated observation mapping. No recipe-to-model
conversion, unknown lysate-pool inference, assay fitting, or batch model is added
by this example. Initial conditions must be audited separately. Replicates of
the same condition stay together, and lysate batches need an additional explicit
held-out design. A successful fit cannot establish model adequacy.

## Published calibration precedent

[Horvath et al., *Toward a genome scale sequence specific dynamic model of
cell-free protein synthesis in Escherichia coli*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7136494/)
fits dynamic protein and metabolite data in a glucose-powered CAT system using
bounded parameter ensembles. Its [author code](https://github.com/varnerlab/Kinetic-CFPS-Model-Publication-Code)
includes raw measurements in `src/data/` and fitted/random ensembles.

Use it to learn how to align model observables with measurements, constrain
parameters from literature, normalize differently scaled channels, retain
multiple compatible parameter sets, and challenge pathways through knockouts.
The next real-data reproduction should run the authors' original model against
its original observations first, preserving its Julia environment and units.
The CAT/glucose study is not the competition's sfGFP/T7 recipe dataset, and our
synthetic example does not reproduce the paper's fitted model or results.
