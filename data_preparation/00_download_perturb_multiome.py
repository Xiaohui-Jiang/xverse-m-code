#!/usr/bin/env python3
"""Download the original Science Perturb-multiome GEO supplementary data."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
from urllib.request import Request, urlopen


SERIES = "GSE274113"
BASE = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/{SERIES}"
DEFAULT_ROOT = "/hpc/group/xielab/xj58/xverse-m-data"


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def fetch(url):
    with urlopen(url, timeout=120) as response:
        return response.read()


def discover(root):
    """Snapshot official metadata and resolve individual files instead of TAR."""
    metadata = root / "metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    listing = fetch(BASE + "/suppl/")
    (metadata / "geo_suppl_listing.html").write_bytes(listing)
    filelist = fetch(BASE + "/suppl/filelist.txt")
    (metadata / "filelist.txt").write_bytes(filelist)
    soft = fetch(BASE + f"/soft/{SERIES}_family.soft.gz")
    gzip.decompress(soft)
    (metadata / f"{SERIES}_family.soft.gz").write_bytes(soft)
    files = []
    for name in sorted(set(re.findall(r'href="(GSE274113_[^"/]+)"', listing.decode()))):
        if name.endswith(".tar"):
            continue
        url = BASE + "/suppl/" + name
        with urlopen(Request(url, method="HEAD"), timeout=120) as response:
            size = int(response.headers["Content-Length"])
        files.append(dict(url=url, path="raw/series/" + name, bytes=size))
    for row in filelist.decode().splitlines():
        columns = row.split("\t")
        if columns[0] != "File":
            continue
        name, size = columns[1], int(columns[3])
        gsm = name.split("_")[0]
        if not re.fullmatch(r"GSM\d+", gsm) or Path(name).name != name:
            raise ValueError(f"Unexpected GEO sample filename: {name}")
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}"
        files.append(dict(url=url, path="raw/samples/" + name, bytes=size))
    if not files or not any("fragments.tsv.gz" in f["path"] for f in files):
        raise ValueError("Incomplete GEO inventory: no fragment files")
    plan = dict(series=SERIES, discovered_at=now(), paper_doi="10.1126/science.ads7951",
                scope="All series supplementary files except duplicate TAR; all TAR members downloaded individually",
                files=files, total_bytes=sum(f["bytes"] for f in files))
    atomic_json(root / "download_plan.json", plan)
    return plan


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def validate(path):
    """Check gzip CRC, HDF5 structure, or XLSX archive integrity."""
    if path.name.endswith(".gz"):
        with gzip.open(path, "rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass
    elif path.name.endswith(".h5"):
        import h5py
        with h5py.File(path, "r") as handle:
            matrix = handle["matrix"]
            shape = matrix["shape"][:]
            if len(matrix["barcodes"]) != shape[1]:
                raise ValueError(f"Invalid barcode dimension in {path}")
            if len(matrix["indptr"]) != shape[1] + 1:
                raise ValueError(f"Invalid CSC pointer dimension in {path}")
            if len(matrix["data"]) != len(matrix["indices"]):
                raise ValueError(f"Invalid CSC data dimension in {path}")
    elif path.name.endswith(".xlsx"):
        import zipfile
        with zipfile.ZipFile(path) as archive:
            if archive.testzip() is not None:
                raise ValueError(f"Invalid XLSX archive: {path}")


def download_one(root, item, attempts):
    target = root / item["path"]
    target.parent.mkdir(parents=True, exist_ok=True)
    receipt = target.with_name(target.name + ".receipt.json")
    if target.exists():
        if target.stat().st_size != item["bytes"]:
            raise ValueError(f"Existing file size mismatch; inspect and move aside: {target}")
        checksum = digest(target)
        if receipt.exists():
            previous = json.loads(receipt.read_text())
            if previous.get("sha256") != checksum or previous.get("url") != item["url"]:
                raise ValueError(f"Existing receipt mismatch: {target}")
            print(f"[Verified existing] {target.name}", flush=True)
            return previous
        validate(target)
    else:
        partial = target.with_name(target.name + ".part")
        print(f"[Downloading] {target.name} ({item['bytes']} bytes)", flush=True)
        if not partial.exists() or partial.stat().st_size != item["bytes"]:
            subprocess.run([
                "curl", "--fail", "--location", "--show-error", "--silent",
                "--retry", str(attempts), "--retry-all-errors", "--retry-delay", "5",
                "--connect-timeout", "60", "--speed-time", "180", "--speed-limit", "1024",
                "--continue-at", "-", "--output", str(partial), item["url"],
            ], check=True)
        if partial.stat().st_size != item["bytes"]:
            raise ValueError(f"Downloaded size mismatch: {partial}")
        # Use the original suffix for format validation before atomic publication.
        if target.name.endswith(".gz"):
            with gzip.open(partial, "rb") as handle:
                while handle.read(8 * 1024 * 1024):
                    pass
        checksum = digest(partial)
        partial.replace(target)
        validate(target) if not target.name.endswith(".gz") else None
    record = dict(item, sha256=checksum, verified_at=now(), status="complete")
    atomic_json(receipt, record)
    print(f"[Complete] {target.name}", flush=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=DEFAULT_ROOT)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--attempts", type=int, default=5)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--refresh-plan", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.attempts < 1:
        parser.error("workers and attempts must be positive")
    root = Path(args.data_root).expanduser().resolve() / "perturb_multiome_gse274113"
    root.mkdir(parents=True, exist_ok=True)
    for directory in ("processed", "qc", "logs"):
        (root / directory).mkdir(exist_ok=True)
    with (root / ".download.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        plan_path = root / "download_plan.json"
        plan = (json.loads(plan_path.read_text()) if plan_path.exists() and not args.refresh_plan
                else discover(root))
        print(f"[Plan] {len(plan['files'])} files; {plan['total_bytes'] / 1e9:.3f} GB; {root}", flush=True)
        if args.plan_only:
            return
        manifest = dict(series=SERIES, started_at=now(), status="running", files=[], errors=[])
        manifest_path = root / "download_manifest.json"
        atomic_json(manifest_path, manifest)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(download_one, root, item, args.attempts): item for item in plan["files"]}
            for future in as_completed(futures):
                try:
                    manifest["files"].append(future.result())
                except Exception as error:
                    manifest["errors"].append(dict(path=futures[future]["path"], error=str(error)))
                    print(f"[ERROR] {futures[future]['path']}: {error}", flush=True)
                atomic_json(manifest_path, manifest)
        manifest["finished_at"] = now()
        manifest["status"] = "failed" if manifest["errors"] else "complete"
        atomic_json(manifest_path, manifest)
        if manifest["errors"]:
            raise RuntimeError("Some files failed; inspect manifest and rerun to resume")
        print(f"[Complete] All {len(manifest['files'])} files verified", flush=True)


if __name__ == "__main__":
    main()
