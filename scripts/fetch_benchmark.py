"""Download pinned public measurements locally; verify Git blob hashes."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
from urllib.request import urlopen


def fetch(manifest_path, output):
    manifest = json.loads(manifest_path.read_text())
    output.mkdir(parents=True, exist_ok=True)
    rows, hashes, columns = [], [], None
    for entry in manifest["files"]:
        url = f'https://raw.githubusercontent.com/brsynth/active_learning_cell_free/{manifest["revision"]}/{entry["path"]}'
        with urlopen(url, timeout=60) as response:
            data = response.read()
        blob = hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest()
        if blob != entry["git_blob_sha1"]:
            raise ValueError(f"Source hash mismatch: {entry['path']}")
        (output / f"plate_{entry['plate']}.csv").write_bytes(data)
        reader = csv.DictReader(io.StringIO(data.decode()))
        if columns is not None and columns != reader.fieldnames:
            raise ValueError("Inconsistent source table schema")
        columns = reader.fieldnames
        for index, row in enumerate(reader):
            rows.append(dict(row, plate=str(entry["plate"]), source_row=str(index+2)))
        hashes.append({"url": url, "sha256": hashlib.sha256(data).hexdigest(), "git_blob_sha1": blob})
    with (output / "observations.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns+["plate", "source_row"])
        writer.writeheader()
        writer.writerows(rows)
    (output / "provenance.json").write_text(json.dumps({"manifest": manifest, "downloads": hashes,
        "rows": len(rows), "normalized_sha256": hashlib.sha256((output / "observations.csv").read_bytes()).hexdigest()}, indent=2)+"\n")
    print(f"Downloaded and verified {len(rows)} measured formulation rows from {len(hashes)} plates")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("benchmarks/borkowski2020.json"))
    parser.add_argument("--output", type=Path, default=Path("data/raw/borkowski2020"))
    args = parser.parse_args()
    fetch(args.manifest, args.output)
