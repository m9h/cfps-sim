"""Extract published constructs from a locally downloaded competition draft.

Usage: uv run python scripts/import_competition_constructs.py draft.txt
DNA/AA sequences are source facts; transcript boundaries remain explicit assumptions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
from Bio.Seq import Seq
from cfps_sim.construct import digest, export_genbank, validate_construct

SOURCE = "https://docs.google.com/document/d/1hFs2_vGpHP_hZ7UB3lRBsfVlN4xVPOoLy1Hu2INrBoo/edit"
TARGETS = [("sfgfp", "m9134733", 2511, 240), ("petase", "m9307404", 2696, 309), ("reteplase", "m9307400", 2966, 399)]


def extract(path, output):
    raw = Path(path).read_bytes()
    text = raw.decode("utf-8-sig")
    dna = re.findall(r"^\s*([acgt]{1000,})\s*$", text, re.M)
    proteins = re.findall(r"^\s*([ACDEFGHIKLMNPQRSTVWY]{200,})\s*$", text, re.M)
    if len(dna) != 3 or len(proteins) != 3:
        raise ValueError("Expected exactly three DNA and three protein sequences in source draft")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for (target, source_id, bp, aa), sequence, protein in zip(TARGETS, dna, proteins):
        sequence = sequence.upper()
        if len(sequence) != bp or len(protein) != aa:
            raise ValueError(f"Source length mismatch for {target}")
        matches = []
        for start in range(len(sequence)-3*(len(protein)+1)+1):
            if sequence[start:start+3] != "ATG":
                continue
            end = start+3*(len(protein)+1)
            try:
                if str(Seq(sequence[start:end]).translate(table=11, cds=True)) == protein:
                    matches.append((start, end))
            except Exception:
                pass
        if len(matches) != 1:
            raise ValueError(f"Expected unique CDS matching reference protein for {target}")
        start, end = matches[0]
        promoter = sequence.rfind("TAATACGACTCACTATA", 0, start)
        if promoter < 0:
            raise ValueError("Cannot find expected T7 promoter")
        record = dict(schema_version=1, id=f"competition-{source_id}-v1", target=target,
                      sequence=sequence, sequence_sha256=digest(sequence), template_format="circular",
                      protein_reference=protein,
                      cds={"start": start, "end": end, "strand": 1},
                      transcript={"start": promoter+17, "end": len(sequence), "strand": 1,
                                  "boundary_status": "assumed",
                                  "basis": "T7 +1 inferred after TAATACGACTCACTATA; transcript assumed through the supplied terminal terminator region to plasmid sequence end. Actual termination/readthrough is unmeasured."},
                      provenance={"source": SOURCE, "source_status": "competition_draft", "retrieved": "2026-09-29",
                                  "source_file_sha256": hashlib.sha256(raw).hexdigest(), "ginkgo_dna_id": source_id,
                                  "verification": "DNA and protein lengths match draft; CDS translation exactly matches draft protein; no experimental sequence verification"})
        validate_construct(record)
        (output/f"{target}.json").write_text(json.dumps(record, indent=2)+"\n")
        export_genbank(record, output/f"{target}.gb")
        print(target, bp, aa, "CDS", start, end, "transcript bounds assumed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, default=Path("constructs/competition"))
    args = parser.parse_args()
    extract(args.source, args.output)
