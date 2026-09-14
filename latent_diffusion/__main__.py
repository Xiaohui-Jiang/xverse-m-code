"""Entry points for frozen xVERSE latent diffusion."""
import argparse

import torch


def positive(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("Must be positive")
    return value


def nonnegative(value):
    value = int(value)
    if value < 0:
        raise argparse.ArgumentTypeError("Must be nonnegative")
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("extract", "train", "sample"):
        p = sub.add_parser(name)
        p.add_argument("--output", required=True, help="New output directory; existing paths are rejected")
        p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
        p.add_argument("--seed", type=int, default=42)
        p.add_argument("--batch-size", type=positive, default=256)
        if name != "train":
            p.add_argument("--xverse-repo", required=True)
            p.add_argument("--checkpoint", required=True, help="Trusted xVERSE checkpoint with training_args")
        if name != "sample":
            p.add_argument("--workers", type=nonnegative, default=0)
        if name != "extract":
            p.add_argument("--cache", required=True)
        if name == "extract":
            p.add_argument("--compiled-root", required=True)
            p.add_argument("--gene-ids", required=True, help="Exact gene list used to train the checkpoint")
            p.add_argument("--train-cells", type=nonnegative, default=100000, help="0 extracts all train cells")
            p.add_argument("--val-cells", type=nonnegative, default=10000, help="0 extracts all validation cells")
        elif name == "train":
            p.add_argument("--width", type=positive, default=512)
            p.add_argument("--depth", type=positive, default=6)
            p.add_argument("--diffusion-steps", type=positive, default=1000)
            p.add_argument("--epochs", type=positive, default=100)
            p.add_argument("--lr", type=float, default=2e-4)
            p.add_argument("--ema-decay", type=float, default=0.999)
            p.add_argument("--patience", type=nonnegative, default=20)
        else:
            p.add_argument("--diffusion-checkpoint", required=True)
            p.add_argument("--cells", type=positive, default=1024)
            p.add_argument("--sampling-steps", type=positive, default=100)
            p.add_argument("--save-means", action="store_true")
    args = parser.parse_args()
    if args.command == "extract":
        from .cache import extract
        extract(args)
    elif args.command == "train":
        if not 0 <= args.ema_decay < 1 or args.lr <= 0:
            parser.error("Require 0 <= ema-decay < 1 and positive learning rate")
        from .train import train
        train(args)
    else:
        if args.cells < 2:
            parser.error("At least two generated cells are needed for diagnostics")
        from .sample import sample
        sample(args)


if __name__ == "__main__":
    main()
