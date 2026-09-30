"""Versioned DNA artifacts with explicit boundaries and sequence provenance.

Coordinates are zero-based, half-open, on the forward strand. Transcript
termination must be annotated explicitly; it is never inferred as verified.
"""
from collections import Counter
import base64
import hashlib
import json
import math
import textwrap
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def digest(sequence):
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


def validate_construct(record):
    if record.get("schema_version") != 1:
        raise ValueError("Unsupported construct schema_version; expected 1")
    for name in ["id", "target", "provenance", "template_format"]:
        if not record.get(name):
            raise ValueError(f"Missing construct {name}")
    if record["template_format"] not in {"circular", "linear"}:
        raise ValueError("Template format must be circular or linear")
    sequence = record.get("sequence", "").upper()
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("DNA must contain only unambiguous A/C/G/T bases")
    if record.get("sequence_sha256") != digest(sequence):
        raise ValueError("DNA checksum does not match; create a new version for modified sequences")
    for name in ["cds", "transcript"]:
        feature = record[name]
        start, end = feature["start"], feature["end"]
        if type(start) is not int or type(end) is not int or not 0 <= start < end <= len(sequence):
            raise ValueError(f"Invalid {name} coordinates")
        if feature.get("strand") != 1:
            raise ValueError("Only contiguous forward-strand features are supported in schema v1")
    cds, transcript = record["cds"], record["transcript"]
    if not transcript["start"] <= cds["start"] < cds["end"] <= transcript["end"]:
        raise ValueError("Transcript must contain the complete CDS")
    if transcript.get("boundary_status") not in {"assumed", "annotated"} or not transcript.get("basis"):
        raise ValueError("Transcript needs an explicit boundary status and basis")
    if not isinstance(record["provenance"], dict) or not record["provenance"].get("source"):
        raise ValueError("Construct provenance needs a source")
    coding = sequence[cds["start"]:cds["end"]]
    if len(coding) % 3 or not coding.startswith("ATG"):
        raise ValueError("CDS must start with ATG and contain complete codons")
    try:
        protein = str(Seq(coding).translate(table=11, cds=True))
    except Exception as exc:
        raise ValueError(f"Invalid bacterial CDS: {exc}") from exc
    if not protein or set(protein) - set(AMINO_ACIDS):
        raise ValueError("Protein must contain only canonical amino acids")
    if record.get("protein_reference") and protein != record["protein_reference"]:
        raise ValueError("Translated CDS does not match the reference protein")
    return record


def demands(record, available_aa_mM=None):
    validate_construct(record)
    sequence = record["sequence"].upper()
    cds, transcript = record["cds"], record["transcript"]
    protein = str(Seq(sequence[cds["start"]:cds["end"]]).translate(table=11, cds=True))
    rna = sequence[transcript["start"]:transcript["end"]].replace("T", "U")
    counts = Counter(protein)
    result = {"construct_id": record["id"], "sequence_sha256": record["sequence_sha256"],
              "template_format": record["template_format"], "plasmid_bp": len(sequence),
              "protein_aa": len(protein), "protein_sequence": protein,
              "transcript_nt": len(rna), "transcript_boundary_status": transcript["boundary_status"],
              "transcript_basis": transcript["basis"],
              "ntp_per_transcript": {base: rna.count(base) for base in "ACGU"},
              "aa_per_protein": {aa: counts[aa] for aa in AMINO_ACIDS},
              "translation_atp_equivalents_per_protein": 2*len(protein),
              "translation_gtp_per_protein": 2*len(protein),
              "note": "Sequence demands only; charging/elongation costs are reduced-model assumptions. NTP demands depend on supplied transcript boundaries."}
    if available_aa_mM is not None:
        if set(available_aa_mM) - set(AMINO_ACIDS):
            raise ValueError("Amino-acid pools must use canonical one-letter codes")
        required = {aa for aa, n in counts.items() if n}
        if not required <= available_aa_mM.keys():
            raise ValueError("All required amino-acid pools must be supplied for a ceiling calculation")
        for value in available_aa_mM.values():
            if not math.isfinite(value) or value < 0:
                raise ValueError("Amino-acid concentrations must be finite and nonnegative")
        ceilings = {aa: available_aa_mM[aa]/counts[aa] for aa in required}
        limiting = min(ceilings, key=ceilings.get)
        result["aa_only_ceiling_protein_mM"] = ceilings[limiting]
        result["limiting_aa"] = limiting
        result["ceiling_note"] = "Static stoichiometric ceiling; excludes recycling, synthesis, endogenous pools, and competing consumption."
    return result


def read_construct(path):
    return validate_construct(json.loads(Path(path).read_text()))


def export_genbank(record, path):
    validate_construct(record)
    seq = SeqRecord(Seq(record["sequence"].upper()), id=record["id"], name=record["target"],
                    description=f"CFPS construct {record['id']}",
                    annotations={"molecule_type": "DNA", "topology": record["template_format"],
                                 "comment": "cfps-sim-metadata-base64:\n" + "\n".join(textwrap.wrap(base64.b64encode(json.dumps({k: v for k, v in record.items() if k not in {"sequence", "cds", "transcript"}}).encode()).decode(), 60))})
    for name, kind in [("cds", "CDS"), ("transcript", "mRNA")]:
        feature = record[name]
        qualifiers = {"label": [name]}
        if name == "cds":
            qualifiers["translation"] = [demands(record)["protein_sequence"]]
            qualifiers["transl_table"] = ["11"]
        else:
            qualifiers.update(boundary_status=[feature["boundary_status"]], basis=[feature["basis"]])
        seq.features.append(SeqFeature(FeatureLocation(feature["start"], feature["end"], strand=1), type=kind, qualifiers=qualifiers))
    SeqIO.write(seq, path, "genbank")


def import_genbank(path, construct_id=None, target=None, transcript_start=None, transcript_end=None):
    seq = SeqIO.read(path, "genbank")
    cds_features = [f for f in seq.features if f.type == "CDS"]
    transcripts = [f for f in seq.features if f.type == "mRNA"]
    if len(cds_features) != 1:
        raise ValueError("Import requires exactly one CDS; export/select a single expression construct first")
    cds = cds_features[0]
    if len(cds.location.parts) != 1 or cds.location.strand != 1:
        raise ValueError("Only contiguous forward-strand CDS features are supported")
    metadata = {}
    comment = seq.annotations.get("comment", "")
    if comment.startswith("cfps-sim-metadata-base64:"):
        metadata = json.loads(base64.b64decode("".join(comment.split(":", 1)[1].split())))
    if transcript_start is not None or transcript_end is not None:
        if transcript_start is None or transcript_end is None:
            raise ValueError("Both transcript-start and transcript-end are required")
        transcript = dict(start=transcript_start, end=transcript_end, strand=1,
                          boundary_status="assumed", basis="User-supplied import boundaries; termination not experimentally verified")
    elif len(transcripts) == 1 and len(transcripts[0].location.parts) == 1:
        feature = transcripts[0]
        transcript = dict(start=int(feature.location.start), end=int(feature.location.end), strand=feature.location.strand,
                          boundary_status=feature.qualifiers.get("boundary_status", ["annotated"])[0],
                          basis=feature.qualifiers.get("basis", ["GenBank mRNA annotation; experimental verification not established"])[0])
    else:
        raise ValueError("Supply transcript boundaries or exactly one mRNA annotation; CDS is not the whole transcript")
    sequence = str(seq.seq).upper()
    if metadata.get("sequence_sha256") and metadata["sequence_sha256"] != digest(sequence):
        raise ValueError("GenBank sequence differs from its metadata checksum; supply a newly versioned artifact with updated provenance")
    record = dict(metadata, schema_version=1, id=construct_id or metadata.get("id") or seq.id,
                  target=target or metadata.get("target") or seq.name,
                  sequence=sequence, sequence_sha256=digest(sequence),
                  template_format=seq.annotations.get("topology", "linear"),
                  cds={"start": int(cds.location.start), "end": int(cds.location.end), "strand": 1},
                  transcript=transcript)
    record.setdefault("provenance", {"source": str(Path(path)), "file_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest()})
    if "translation" in cds.qualifiers:
        record["protein_reference"] = cds.qualifiers["translation"][0]
    return validate_construct(record)
