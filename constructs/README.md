# Construct artifacts

`competition/` contains the three supplied plasmids transcribed from the
[competition draft rules](https://docs.google.com/document/d/1hFs2_vGpHP_hZ7UB3lRBsfVlN4xVPOoLy1Hu2INrBoo/edit)
on 2026-09-29. Each CDS translates exactly to its accompanying source protein,
including tags. This verifies consistency with the document, not experimental
identity or activity. JSON and GenBank contain the same sequence and provenance.

| Target | Source DNA ID | Plasmid bp | Protein aa |
|---|---|---:|---:|
| sfGFP | m9134733 | 2511 | 240 |
| PETase/TfCut2 | m9307404 | 2696 | 309 |
| reteplase/vtPA | m9307400 | 2966 | 399 |

The transcript starts immediately after the upstream T7 promoter motif and ends
at the supplied sequence boundary. **Both are modeling assumptions**, especially
the termination boundary; circular templates may produce other transcript
lengths. NTP counts depend on these bounds. Amino-acid counts follow the matched
CDS and do not depend on transcript termination. Change assumptions in a new
version and retain their basis; do not imply experimental verification.

## Schema v1

JSON requires `schema_version: 1`, `id`, `target`, `template_format` (`circular`
or `linear`), unambiguous `sequence`, its uppercase SHA-256 `sequence_sha256`,
`provenance.source`, and `cds`/`transcript` objects. Coordinates are **zero-based,
half-open**; only contiguous forward-strand (`strand: 1`) features are supported.
The complete CDS includes its start and terminal stop codons. The transcript
contains the CDS and specifies `boundary_status` (`assumed` or `annotated`) and
`basis`. Annotation alone does not establish experimental verification.
An optional `protein_reference` must match translation using bacterial table 11.

GenBank export includes CDS, mRNA, topology, translation, and a base64-encoded
JSON metadata comment for lossless provenance exchange. Import requires exactly
one CDS and either one mRNA annotation or explicit transcript coordinates. It
rejects ambiguous bases, reverse/compound CDS locations, and CDS-only files
without transcript bounds. Schema v1 does not support circular origin-crossing
features or multiple expression units.

To regenerate from a fresh text export of the same rules:

```sh
uv run python scripts/import_competition_constructs.py /path/to/rules.txt
```

The importer checks the three expected lengths and exact protein matches. It
records the source file checksum. Review source changes before replacing a
versioned artifact; the source remains a draft.
