# Project documentation

Keep `docs/article.Pnw` synchronized with project changes that affect scientific
behavior, inputs, assumptions, installation, workflows, or validation. Update its
prose and executable examples, run `uv run --extra docs python scripts/build_docs.py`,
and include the rebuilt `docs/article.pdf` with those changes. The article is the
project's primary scientific documentation; updating the README alone is not
sufficient.

Distinguish source facts from assumptions and uncalibrated simulations. Preserve
construct provenance and version changed sequences. Keep experimental
observations separate from synthetic or mechanistic predictions.
