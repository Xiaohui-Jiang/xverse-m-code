"""Regression tests for diffusion algebra and frozen xVERSE integration."""
import argparse
import json
import os
from pathlib import Path
import tempfile
import unittest

import numpy as np
import torch

from latent_diffusion.cache import extract
from latent_diffusion.model import Denoiser, Diffusion
from latent_diffusion.sample import sample
from latent_diffusion.train import train
from latent_diffusion.xverse import decode, encode, imports, load_frozen, sha256

REPO = os.environ.get("XVERSE_REPO", "/Users/xiaohui/LocalFiles/Codes/xverse-code")
torch.set_num_threads(1)


class DiffusionTests(unittest.TestCase):
    def test_velocity_roundtrip_and_sampling(self):
        d = Diffusion(32)
        x, noise = torch.randn(10, 8), torch.randn(10, 8)
        t = torch.arange(10) * 3
        noisy, velocity = d.corrupt(x, t, noise)
        recovered, recovered_noise = d.recover(noisy, velocity, t)
        torch.testing.assert_close(x, recovered)
        torch.testing.assert_close(noise, recovered_noise)
        m = Denoiser(8, 16, 2)
        loss = d.loss(m, x)
        loss.backward()
        self.assertGreater(m.output[-1].weight.grad.abs().sum().item(), 0)
        a = d.sample(m, 10, 8, 8, torch.Generator().manual_seed(2))
        b = d.sample(m, 10, 8, 8, torch.Generator().manual_seed(2))
        torch.testing.assert_close(a, b)
        self.assertTrue(torch.isfinite(a).all())
        with self.assertRaises(ValueError):
            d.sample(m, 10, 8, 33)


@unittest.skipUnless(Path(REPO, "main/utils_model.py").exists(), "Set XVERSE_REPO to the reference checkout")
class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.module, self.data = imports(REPO)

    def tearDown(self):
        self.tmp.cleanup()

    def checkpoint(self, distribution="nb", mode="cell_gene"):
        config = dict(normalize_z_bio=True, recon_distribution=distribution,
                      nb_dispersion_mode=mode, expr_encoder_depth=2,
                      mask_encoder_depth=2, decoder_depth=1)
        torch.manual_seed(4)
        m = self.module.XVerseModel(total_gene=6, hidden_dim=8, num_samples=2,
                                   total_cell_type=3, expr_hidden_dim=16,
                                   mask_hidden_dim=16, dec_hidden_dim=16, **config).eval()
        path = self.root / f"{distribution}_{mode}.pth"
        torch.save({"model_state_dict": m.state_dict(), "training_args": config}, path)
        return path, m

    def fixture(self):
        root = self.root / "compiled"
        root.mkdir()
        (root / "global_gene_ids.txt").write_text("\n".join(f"g{i}" for i in range(6)) + "\n")
        meta = root / "panel.npz"
        np.savez(meta, gene_ids=np.array(["g0", "g1", "g2", "g3"]))
        manifest = {"format": "xverse_train_v1", "global_num_genes": 6,
                    "source_pairs": [{"pair_id": 0, "meta_path": str(meta)}], "splits": {}}
        for split, n in (("train", 32), ("val", 16)):
            directory = root / split
            directory.mkdir()
            for key, value in {
                "cell_ptr": np.arange(n + 1), "gene_idx": np.zeros(n, dtype="int32"),
                "gene_val": np.arange(n) % 8 + 1, "celltype_id": np.zeros(n, dtype="int64"),
                "sample_id": np.arange(n) % 2, "source_pair_id": np.zeros(n, dtype="int64"),
            }.items():
                np.save(directory / f"{key}.npy", value)
            manifest["splits"][split] = {"num_cells": n, "shards": [
                {"path": str(directory), "global_cell_start": 0, "global_cell_end": n}]}
        (root / "manifest.json").write_text(json.dumps(manifest))
        return root

    def test_strict_load_and_likelihoods(self):
        for distribution, mode in [("poisson", "gene"), ("zinb", "gene")] + [
                ("nb", m) for m in ("global", "gene", "cell", "factorized", "cell_gene", "lowrank")]:
            with self.subTest(distribution=distribution, mode=mode):
                path, original = self.checkpoint(distribution, mode)
                frozen = load_frozen(REPO, path)
                self.assertFalse(frozen.training)
                self.assertTrue(all(not p.requires_grad for p in frozen.parameters()))
                x = torch.ones(4, 6); mask = torch.ones_like(x, dtype=torch.bool)
                z = encode(frozen, x, mask)
                torch.testing.assert_close(z, encode(original, x, mask))
                mean, counts = decode(frozen, z)
                self.assertEqual(counts.shape, x.shape)
                self.assertTrue(torch.isfinite(mean).all())
                self.assertTrue((counts >= 0).all())
                torch.testing.assert_close(counts, counts.round())
        ckpt = torch.load(path, weights_only=False)
        del ckpt["model_state_dict"]["decoder.out.bias"]
        torch.save(ckpt, path)
        with self.assertRaises(RuntimeError):
            load_frozen(REPO, path)

    def test_pipeline_preserves_frozen_weights_and_masks(self):
        path, _ = self.checkpoint()
        before = sha256(path)
        root = self.fixture()
        ds = self.data.CompiledShardDataset(root)
        batch = self.data.CompiledSparseBatchCollator(ds, 6)([ds[0]])
        _, _, counts, mask = self.module._unpack_batch(batch, 6, "cpu")
        self.assertEqual(counts[0, 1], 0)
        self.assertTrue(mask[0, 1])  # Measured zero.
        self.assertFalse(mask[0, 5])  # Unmeasured gene.
        cache, run, samples = (self.root / name for name in ("cache", "run", "samples"))
        args = argparse.Namespace(xverse_repo=REPO, checkpoint=str(path), device="cpu",
                                  compiled_root=str(root), gene_ids=str(root / "global_gene_ids.txt"),
                                  output=str(cache), train_cells=32, val_cells=16, seed=2,
                                  batch_size=8, workers=0)
        extract(args)
        self.assertEqual(np.load(cache / "train.npy").shape, (32, 8))
        with np.load(cache / "normalization.npz") as stats:
            np.testing.assert_allclose(stats["mean"], np.load(cache / "train.npy").mean(0), atol=1e-6)
        train(argparse.Namespace(cache=str(cache), device="cpu", seed=2, width=16, depth=2,
                                 diffusion_steps=32, lr=1e-3, batch_size=8, workers=0,
                                 output=str(run), epochs=3, ema_decay=0.9, patience=0))
        sample(argparse.Namespace(cache=str(cache), device="cpu", seed=2, xverse_repo=REPO,
                                  checkpoint=str(path), diffusion_checkpoint=str(run / "best.pt"),
                                  output=str(samples), cells=8, batch_size=4, sampling_steps=8,
                                  save_means=False))
        z = np.load(samples / "latents.npy")
        np.testing.assert_allclose(np.linalg.norm(z, axis=1), 1, atol=1e-6)
        self.assertEqual(before, sha256(path))
        self.assertTrue((samples / "manifest.json").exists())
        with self.assertRaises(FileExistsError):
            extract(args)
        (root / "global_gene_ids.txt").write_text("g1\ng0\ng2\ng3\ng4\ng5\n")
        # A different supplied order must fail, despite an equal gene count.
        reference = self.root / "genes.txt"
        reference.write_text("g0\ng1\ng2\ng3\ng4\ng5\n")
        args.gene_ids = str(reference)
        with self.assertRaises(ValueError):
            extract(args)


if __name__ == "__main__":
    unittest.main()
