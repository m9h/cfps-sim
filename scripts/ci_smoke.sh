#!/usr/bin/env bash
# Host-side CLI smoke test: everything here must work without a CUDA GPU.
# Run locally with: bash scripts/ci_smoke.sh
set -euo pipefail

out="$(mktemp -d)"
trap 'rm -rf "$out"' EXIT
run() { echo "+ cfps-sim $*" >&2; uv run --frozen cfps-sim "$@"; }
field() { python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])' "$1" "$2"; }

run --help >/dev/null
run validate examples/recipe.json >"$out/validate.json"
run cost examples/recipe.json examples/prices.draft.json >"$out/cost.json"
run score examples/recipe.json examples/prices.illustrative.json \
  examples/yield.synthetic.json >"$out/score.json"

if [ "$(field "$out/cost.json" complete)" != "True" ]; then
  echo "draft price ledger no longer covers the example recipe" >&2
  exit 1
fi

# preflight must keep exiting 2 while readiness gates are blocked.
status=0
run preflight examples/recipe.json --prices examples/prices.draft.json \
  >"$out/preflight.json" || status=$?
if [ "$status" -ne 2 ]; then
  echo "preflight exited $status, expected 2 (blocked readiness gates)" >&2
  exit 1
fi
if [ "$(field "$out/preflight.json" submission_ready)" != "False" ]; then
  echo "preflight exited 2 but reports submission_ready" >&2
  exit 1
fi

# Every competition construct inspects, and survives a GenBank round trip.
for target in sfgfp petase reteplase; do
  src="constructs/competition/$target.json"
  run construct inspect "$src" >"$out/$target.inspect.json"
  run construct export-genbank "$src" --output "$out/$target.gb" >/dev/null
  run construct import-genbank "$out/$target.gb" --output "$out/$target.json" >/dev/null
  if [ "$(field "$src" sequence_sha256)" != "$(field "$out/$target.json" sequence_sha256)" ]; then
    echo "$target: sequence checksum changed across GenBank round trip" >&2
    exit 1
  fi
done

# Numerical commands must refuse to run without CUDA rather than fall back to CPU.
if ! uv run --frozen python -c 'import torch,sys; sys.exit(torch.cuda.is_available())'; then
  echo "CUDA present: skipping the no-GPU refusal check" >&2
else
  status=0
  run simulate --output "$out/demo" >/dev/null 2>"$out/simulate.err" || status=$?
  if [ "$status" -eq 0 ] || ! grep -q "requires a CUDA GPU" "$out/simulate.err"; then
    echo "simulate did not refuse to run without CUDA" >&2
    cat "$out/simulate.err" >&2
    exit 1
  fi
fi

echo "CLI smoke test passed"
