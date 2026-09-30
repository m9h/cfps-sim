from copy import deepcopy

import numpy as np
import pytest
import roadrunner
from Bio import SeqIO
from Bio.Seq import Seq

from cfps_sim.construct import AMINO_ACIDS, demands, export_genbank, import_genbank, read_construct, validate_construct
from cfps_sim.model import load_model, simulate


@pytest.mark.parametrize("target,bp,length", [("sfgfp", 2511, 240), ("petase", 2696, 309), ("reteplase", 2966, 399)])
def test_source_sequences_and_exchange(target, bp, length, tmp_path):
    record = read_construct(f"constructs/competition/{target}.json")
    counts = demands(record)
    assert counts["plasmid_bp"] == bp
    assert counts["protein_aa"] == length
    assert counts["protein_sequence"] == record["protein_reference"]
    assert sum(counts["aa_per_protein"].values()) == length
    assert sum(counts["ntp_per_transcript"].values()) == counts["transcript_nt"]
    assert counts["transcript_boundary_status"] == "assumed"
    output = tmp_path / "construct.gb"
    export_genbank(record, output)
    imported = import_genbank(output)
    assert imported == record


def test_reject_modified_sequence_and_implicit_transcript(tmp_path):
    record = read_construct("constructs/competition/sfgfp.json")
    modified = deepcopy(record)
    modified["sequence"] = "A" + modified["sequence"][1:]
    with pytest.raises(ValueError, match="checksum"):
        validate_construct(modified)
    seq = SeqIO.read("constructs/competition/sfgfp.gb", "genbank")
    seq.seq = Seq("A" + str(seq.seq)[1:])
    altered = tmp_path / "modified.gb"
    SeqIO.write(seq, altered, "genbank")
    with pytest.raises(ValueError, match="checksum"):
        import_genbank(altered)
    seq = SeqIO.read("constructs/competition/sfgfp.gb", "genbank")
    seq.annotations.pop("comment")
    seq.features = [f for f in seq.features if f.type == "CDS"]
    output = tmp_path / "cds_only.gb"
    SeqIO.write(seq, output, "genbank")
    with pytest.raises(ValueError, match="transcript boundaries"):
        import_genbank(output)


def test_sequence_balance_limitation_and_sbml_replay():
    record = read_construct("constructs/competition/sfgfp.json")
    derived = demands(record, {aa: 2 for aa in AMINO_ACIDS})
    expected = min(2/n for n in derived["aa_per_protein"].values() if n)
    assert derived["aa_only_ceiling_protein_mM"] == expected
    rr, trajectory = simulate(construct=record)
    assert np.isfinite(trajectory).all()
    assert np.min(trajectory[:, 1:]) >= -1e-9
    for aa, count in derived["aa_per_protein"].items():
        assert rr[f"aa_{aa}"] + count*(rr["Protein"] + rr["GFP"]) == pytest.approx(2, abs=1e-7)
    replay = roadrunner.RoadRunner(rr.getSBML())
    assert replay["aa_A"] == 2  # Export retains initial conditions, not final state.
    for base, count in derived["ntp_per_transcript"].items():
        assert replay[f"nt_{base.lower()}"] == count
    replay.simulate(0, 240, 241)
    assert replay["GFP"] == pytest.approx(rr["GFP"], rel=1e-5)
    depleted = load_model(record)
    depleted["aa_L"] = 0
    depleted.simulate(0, 240, 20)
    assert depleted["GFP"] == pytest.approx(0, abs=1e-12)
    with pytest.raises(ValueError, match="sfGFP only"):
        load_model(read_construct("constructs/competition/petase.json"))
