"""Train only the latent denoiser; xVERSE is not loaded by this module."""
import copy
import json
import random
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

from .model import Denoiser, Diffusion
from .xverse import sha256


class Latents(Dataset):
    def __init__(self, root, split):
        self.values = np.load(Path(root) / f"{split}.npy", mmap_mode="r")
        with np.load(Path(root) / "normalization.npz") as stats:
            self.mean, self.std = stats["mean"], stats["std"]

    def __len__(self):
        return len(self.values)

    def __getitem__(self, index):
        return torch.from_numpy((self.values[index] - self.mean) / self.std)


def seed_all(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def validate(model, diffusion, loader, device, seed):
    model.eval()
    generator = torch.Generator(device=device).manual_seed(seed)
    total = 0.0
    with torch.inference_mode():
        for x in loader:
            total += diffusion.loss(model, x.to(device), generator).item() * len(x)
    return total / len(loader.dataset)


def train(args):
    seed_all(args.seed)
    device = torch.device(args.device)
    root = Path(args.cache)
    manifest = json.loads((root / "manifest.json").read_text())
    if manifest["format"] != "xverse_latents_v1":
        raise ValueError("Unsupported cache format")
    for split in ("train", "val"):
        if sha256(root / f"{split}.npy") != manifest["splits"][split]["sha256"]:
            raise ValueError(f"Corrupt {split} latent cache")
    config = {"dim": manifest["dim"], "width": args.width, "depth": args.depth}
    model = Denoiser(**config).to(device)
    ema = copy.deepcopy(model).requires_grad_(False).eval()
    diffusion = Diffusion(args.diffusion_steps, device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    train_loader = DataLoader(Latents(root, "train"), batch_size=args.batch_size, shuffle=True,
                              num_workers=args.workers, pin_memory=device.type == "cuda")
    val_loader = DataLoader(Latents(root, "val"), batch_size=args.batch_size, shuffle=False,
                            num_workers=args.workers)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    print(f"Denoiser parameters: {sum(p.numel() for p in model.parameters()):,}", flush=True)
    initial = validate(ema, diffusion, val_loader, device, args.seed + 100)
    best, stale = float("inf"), 0
    history = [{"epoch": 0, "val_velocity_mse": initial}]
    (out / "history.json").write_text(json.dumps(history, indent=2))
    for epoch in range(1, args.epochs + 1):
        model.train()
        total = 0.0
        for x in train_loader:
            optimizer.zero_grad(set_to_none=True)
            loss = diffusion.loss(model, x.to(device))
            if not torch.isfinite(loss):
                raise FloatingPointError("Non-finite training loss")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0, error_if_nonfinite=True)
            optimizer.step()
            with torch.no_grad():
                for target, source in zip(ema.parameters(), model.parameters()):
                    target.lerp_(source, 1 - args.ema_decay)
            total += loss.item() * len(x)
        score = validate(ema, diffusion, val_loader, device, args.seed + 100)
        if not np.isfinite(score):
            raise FloatingPointError("Non-finite validation loss")
        row = {"epoch": epoch, "train_velocity_mse": total / len(train_loader.dataset),
               "val_velocity_mse": score}
        history.append(row); print(json.dumps(row), flush=True)
        improved = score < best
        best = min(best, score)
        stale = 0 if improved else stale + 1
        payload = {"format": "xverse_diffusion_v1", "config": config,
                   "diffusion_steps": args.diffusion_steps, "model": model.state_dict(),
                   "ema": ema.state_dict(), "optimizer": optimizer.state_dict(),
                   "epoch": epoch, "val_velocity_mse": score, "best_val": best,
                   "args": vars(args), "cache_manifest": manifest,
                   "normalization": {"mean": torch.from_numpy(train_loader.dataset.mean.copy()),
                                     "std": torch.from_numpy(train_loader.dataset.std.copy())}}
        torch.save(payload, out / "last.pt")
        if improved:
            torch.save(payload, out / "best.pt")
        (out / "history.json").write_text(json.dumps(history, indent=2) + "\n")
        if args.patience > 0 and stale >= args.patience:
            break
    print(f"Initial validation={initial:.6f}; best={best:.6f}", flush=True)
