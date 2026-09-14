#!/usr/bin/env python3
"""Compare a suspect local BGZF tail to GEO without downloading it again."""

import argparse
import hashlib
import json
from pathlib import Path
import struct
from urllib.request import Request, urlopen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default="/hpc/group/xielab/xj58/xverse-m-data")
    parser.add_argument("--filename", default="GSM8443612_rep2_atac_fragments.tsv.gz")
    args = parser.parse_args()
    root = Path(args.data_root) / "perturb_multiome_gse274113"
    plan = json.loads((root / "download_plan.json").read_text())
    item = next(f for f in plan["files"] if Path(f["path"]).name == args.filename)
    path = root / item["path"]
    if not path.exists():
        path = path.with_name(path.name + ".part")
    size = path.stat().st_size
    with path.open("rb") as handle:
        handle.seek(max(0, size - 131072))
        local = handle.read()
    request = Request(item["url"], headers={"Range": f"bytes={size-len(local)}-{size-1}"})
    with urlopen(request, timeout=120) as response:
        expected_range = f"bytes {size-len(local)}-{size-1}/{item['bytes']}"
        if response.status != 206 or response.headers.get("Content-Range") != expected_range:
            raise ValueError("Server did not honor the requested byte range")
        remote = response.read(len(local) + 1)
    report = dict(file=args.filename, url=item["url"], local_bytes=size,
                  expected_bytes=item["bytes"], tail_matches_source=(local == remote),
                  tail_sha256=hashlib.sha256(local).hexdigest())
    # Locate a standard BGZF header. This is a tail diagnostic, not a full validator.
    offset = local.rfind(b"\x1f\x8b\x08\x04")
    if offset >= 0 and offset + 18 <= len(local) and local[offset+12:offset+16] == b"BC\x02\x00":
        report["last_bgzf_block_available_bytes"] = len(local) - offset
        report["last_bgzf_block_expected_bytes"] = struct.unpack("<H", local[offset+16:offset+18])[0] + 1
    destination = root / "metadata" / (args.filename + ".tail_diagnostic.json")
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
