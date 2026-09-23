#!/usr/bin/env python3
"""Export an explicit working-tree selection with content hashes, without Git writes."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def build(output: Path) -> None:
    spec = json.loads((ROOT / "config/source_manifest.json").read_text())
    selected = {name: name for name in spec["files"]}
    selected.update({name: name for name in spec.get("optional_files", [])
                     if (ROOT / name).is_file()})
    for name in spec["trees"]:
        for path in sorted((ROOT / name).rglob("*")):
            relative = path.relative_to(ROOT)
            if any(part in spec["excluded_parts"] for part in relative.parts):
                continue
            if path.suffix in spec["excluded_suffixes"] or not path.is_file():
                continue
            selected[str(relative)] = str(relative)
    entries = {}
    for source, target in selected.items():
        path = ROOT / source
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"release source is missing or a symlink: {source}")
        entries[target] = path.read_bytes()
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit = None
    provenance = {
        "base_commit": commit,
        "snapshot": "selected working-tree contents, including uncommitted files",
        "files": {
            name: {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
            for name, data in sorted(entries.items())
        },
    }
    entries["SOURCE_MANIFEST.json"] = (json.dumps(provenance, indent=2) + "\n").encode()
    output.parent.mkdir(parents=True, exist_ok=True)
    # Refuse accidental replacement of an existing release artifact.
    with zipfile.ZipFile(output, "x", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(entries.items()):
            info = zipfile.ZipInfo("rmi-solver/" + name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(output.suffix + ".sha256").write_text(f"{digest}  {output.name}\n")
    print(f"{output}: {len(entries)} files; SHA-256 {digest}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    build(parser.parse_args().output)
