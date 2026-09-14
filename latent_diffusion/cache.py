"""Stream deterministic frozen embeddings from xVERSE compiled sparse data."""
import json
import subprocess
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset

from .xverse import encode, imports, load_frozen, sha256


def read_genes(path):
    genes = Path(path).read_text().splitlines()
    if not genes or len(set(genes)) != len(genes) or any(not g for g in genes):
        raise ValueError("Gene IDs must be nonempty and unique")
    return genes


def extract(args):
    device = torch.device(args.device)
    checkpoint_digest = sha256(args.checkpoint)
    model = load_frozen(args.xverse_repo, args.checkpoint, device)
    module, data = imports(args.xverse_repo)
    root = Path(args.compiled_root).resolve()
    genes = read_genes(args.gene_ids)
    if genes != read_genes(root / "global_gene_ids.txt") or len(genes) != model.total_gene:
        raise ValueError("Compiled gene order must exactly match the checkpoint's supplied gene list")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    manifest = {
        "format": "xverse_latents_v1", "checkpoint": str(Path(args.checkpoint).resolve()),
        "checkpoint_sha256": checkpoint_digest, "xverse_repo": str(Path(args.xverse_repo).resolve()),
        "source_code_sha256": {name: sha256(Path(args.xverse_repo) / name)
                               for name in ("main/utils_model.py", "main/data.py")},
        "compiled_root": str(root), "compiled_manifest_sha256": sha256(root / "manifest.json"),
        "gene_ids_sha256": sha256(args.gene_ids), "dim": model.hidden_dim,
        "normalize_z_bio": model.normalize_z_bio, "seed": args.seed,
        "selection": "uniform_without_replacement_within_existing_split", "splits": {},
    }
    try:
        manifest["xverse_commit"] = subprocess.check_output(
            ["git", "-C", args.xverse_repo, "rev-parse", "HEAD"], text=True).strip()
    except subprocess.CalledProcessError:
        manifest["xverse_commit"] = None
    (out / "gene_ids.txt").write_text("\n".join(genes) + "\n")
    # Completion manifest is written last. Incomplete directories are never reused.
    for split, cap in (("train", args.train_cells), ("val", args.val_cells)):
        dataset = data.CompiledShardDataset(str(root), split, max_cached_shards=4)
        if not dataset.has_observed_panel_masks:
            raise ValueError("True measured-panel masks are mandatory; no nonzero-only fallback")
        n = len(dataset) if cap == 0 else min(cap, len(dataset))
        if n < 2:
            raise ValueError(f"Need at least two cells in {split}")
        rng = np.random.default_rng(args.seed + (split == "val"))
        indices = np.arange(n) if n == len(dataset) else np.sort(rng.choice(len(dataset), n, replace=False))
        np.save(out / f"{split}_indices.npy", indices)
        latents = np.lib.format.open_memmap(out / f"{split}.npy", mode="w+", dtype="float32", shape=(n, model.hidden_dim))
        metadata = np.lib.format.open_memmap(out / f"{split}_labels.npy", mode="w+", dtype="int64", shape=(n, 2))
        loader = DataLoader(Subset(dataset, indices), batch_size=args.batch_size,
                            num_workers=args.workers, pin_memory=device.type == "cuda",
                            collate_fn=data.CompiledSparseBatchCollator(dataset, model.total_gene))
        offset = 0
        total = np.zeros(model.hidden_dim, dtype=np.float64)
        squares = total.copy()
        with torch.inference_mode():
            for batch in loader:
                sample, celltype, counts, mask = module._unpack_batch(batch, model.total_gene, device)
                z = encode(model, counts, mask).float().cpu().numpy()
                if not np.isfinite(z).all():
                    raise FloatingPointError("Non-finite frozen embeddings")
                latents[offset:offset + len(z)] = z
                metadata[offset:offset + len(z)] = np.column_stack((sample.cpu().numpy(), celltype.cpu().numpy()))
                total += z.sum(axis=0, dtype=np.float64)
                squares += np.square(z.astype(np.float64)).sum(axis=0)
                offset += len(z)
                if offset % (args.batch_size * 20) == 0:
                    print(f"extract {split}: {offset}/{n}", flush=True)
        latents.flush(); metadata.flush()
        if offset != n:
            raise RuntimeError("Incomplete extraction")
        manifest["splits"][split] = {"cells": n, "sha256": sha256(out / f"{split}.npy")}
        if split == "train":
            mean = total / n
            std = np.sqrt(np.maximum(squares / n - mean**2, 0)).clip(1e-4)
            np.savez(out / "normalization.npz", mean=mean.astype("float32"), std=std.astype("float32"))
        print(f"extract {split} complete: {n} cells", flush=True)
    if sha256(args.checkpoint) != manifest["checkpoint_sha256"]:
        raise RuntimeError("Source checkpoint changed during extraction; cache is incomplete")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
