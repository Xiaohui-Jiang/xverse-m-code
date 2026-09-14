"""Generate latent states, decode sparse counts, and report latent diagnostics."""
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.stats import wasserstein_distance
import torch
import torch.nn.functional as F

from .model import Denoiser, Diffusion
from .train import seed_all
from .xverse import decode, load_frozen, sha256


def compare(a, b, seed=0):
    """Descriptive latent metrics, not proof of biological generation fidelity."""
    rng = np.random.default_rng(seed)
    directions = rng.standard_normal((a.shape[1], 64))
    directions /= np.linalg.norm(directions, axis=0)
    pa, pb = a @ directions, b @ directions
    return {"mean_rmse": float(np.sqrt(np.mean((a.mean(0) - b.mean(0))**2))),
            "covariance_rmse": float(np.sqrt(np.mean((np.cov(a, rowvar=False) - np.cov(b, rowvar=False))**2))),
            "sliced_wasserstein_64": float(np.mean([wasserstein_distance(pa[:, j], pb[:, j]) for j in range(64)]))}


def sample(args):
    seed_all(args.seed)
    device = torch.device(args.device)
    checkpoint = torch.load(args.diffusion_checkpoint, map_location="cpu", weights_only=False)
    if checkpoint.get("format") != "xverse_diffusion_v1":
        raise ValueError("Unsupported diffusion checkpoint")
    provenance = checkpoint["cache_manifest"]
    if sha256(args.checkpoint) != provenance["checkpoint_sha256"]:
        raise ValueError("xVERSE checkpoint differs from the embedding source")
    for name, digest in provenance["source_code_sha256"].items():
        if sha256(Path(args.xverse_repo) / name) != digest:
            raise ValueError(f"xVERSE source changed: {name}")
    cache = Path(args.cache)
    if json.loads((cache / "manifest.json").read_text()) != provenance:
        raise ValueError("Cache provenance differs from diffusion checkpoint")
    model = Denoiser(**checkpoint["config"]).to(device)
    model.load_state_dict(checkpoint["ema"], strict=True); model.eval()
    diffusion = Diffusion(checkpoint["diffusion_steps"], device)
    frozen = load_frozen(args.xverse_repo, args.checkpoint, device)
    mean = checkpoint["normalization"]["mean"].to(device)
    std = checkpoint["normalization"]["std"].to(device)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    generator = torch.Generator(device=device).manual_seed(args.seed)
    latents, raw_norms, shards = [], [], []
    for start in range(0, args.cells, args.batch_size):
        n = min(args.batch_size, args.cells - start)
        standardized = diffusion.sample(model, n, len(mean), args.sampling_steps, generator)
        z = standardized * std + mean
        raw_norms.append(z.norm(dim=-1).cpu().numpy())
        # The frozen decoder was trained on this sphere. Preserve preprojection
        # diagnostics; this correction is not a learned manifold constraint.
        if frozen.normalize_z_bio:
            z = F.normalize(z, dim=-1)
        expected, counts = decode(frozen, z)
        filename = f"counts_{start:09d}.npz"
        sparse.save_npz(out / filename, sparse.csr_matrix(counts.cpu().numpy()))
        if args.save_means:
            np.save(out / f"mean_{start:09d}.npy", expected.cpu().numpy())
        shards.append({"path": filename, "start": start, "cells": n})
        latents.append(z.cpu().numpy())
    z = np.concatenate(latents)
    np.save(out / "latents.npy", z)
    norms = np.concatenate(raw_norms)
    np.save(out / "preprojection_norms.npy", norms)
    (out / "gene_ids.txt").write_text((cache / "gene_ids.txt").read_text())
    train_values = np.load(cache / "train.npy", mmap_mode="r")
    val_values = np.load(cache / "val.npy", mmap_mode="r")
    rng = np.random.default_rng(args.seed)
    n = min(len(z), len(val_values), len(train_values), 2048)
    real = np.asarray(val_values[rng.choice(len(val_values), n, replace=False)])
    bootstrap = np.asarray(train_values[rng.choice(len(train_values), n, replace=True)])
    gaussian = rng.normal(size=(n, z.shape[1])) * std.cpu().numpy() + mean.cpu().numpy()
    if frozen.normalize_z_bio:
        gaussian /= np.maximum(np.linalg.norm(gaussian, axis=1, keepdims=True), 1e-12)
    report = {"cells_per_comparison": n, "diffusion_vs_val": compare(z[:n], real),
              "train_bootstrap_vs_val": compare(bootstrap, real),
              "diagonal_gaussian_vs_val": compare(gaussian, real),
              "preprojection_norm_mean": float(norms.mean()),
              "preprojection_norm_sd": float(norms.std()),
              "sphere_projection": frozen.normalize_z_bio,
              "scope": "Descriptive latent diagnostics on validation, not held-out biological evidence."}
    (out / "diagnostics.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "manifest.json").write_text(json.dumps({
        "format": "xverse_diffusion_samples_v1", "shards": shards,
        "diffusion_checkpoint_sha256": sha256(args.diffusion_checkpoint),
        "xverse_checkpoint_sha256": provenance["checkpoint_sha256"],
        "source_cache": str(cache.resolve()), "seed": args.seed,
        "likelihood": frozen.recon_distribution, "decoder": "mu_bio",
        "sphere_projection": frozen.normalize_z_bio,
        "sampling_steps": args.sampling_steps}, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
