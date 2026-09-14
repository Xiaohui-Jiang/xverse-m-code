# Frozen xVERSE latent diffusion

This first RNA prototype learns the distribution of the existing xVERSE
`z_bio` embeddings. The pretrained encoder, decoder, dispersion heads, and
sample embeddings are frozen. No new autoencoder is trained and no Gaussian
KL penalty is added to the embeddings. This is currently unconditional
population generation, not template-conditioned biological-copy generation.

## Data contract

The original `main/data.py` supports legacy sparse NPZ blocks and compiled
`xverse_train_v1` memory-mapped shards. This prototype reuses the latter:
`CompiledShardDataset`, `CompiledSparseBatchCollator`, and the original sparse
batch unpacker. Counts become dense only for the current encoder batch.
`source_pair_id.npy` and source metadata reconstruct the full measured panel,
including measured zeros. Missing panel metadata causes a hard failure.

Pass the exact gene list used for checkpoint training. Its complete order is
compared with `global_gene_ids.txt`, not merely its length. The checkpoint must
contain architecture metadata in `training_args`; all tensors load strictly.
The old generic CLI's permissive checkpoint loader is not used.

Extraction uses one deterministic full-panel encoder pass under eval mode and
inference mode, without running the decoder or auxiliary heads. The existing
train/val splits are preserved. By default a fixed uniform cell sample is drawn
without replacement from each split, then sorted for shard locality. This is
cell-weighted population sampling, not xVERSE's sample-balanced pretraining
sampler. It includes all assay panels present in the compiled dataset; there
is no RNA-platform-only filter. A learned density may retain panel/batch effects.

Latents, sample/cell-type IDs, and original split indices are cached separately.
Labels are retained for auditing only and are not denoiser inputs. Mean and
coordinate scales are fitted on training embeddings only. Standardization
changes units; it does not Gaussianize the distribution. Output directories
are immutable, and a completion manifest is written only after extraction.
Source checkpoint, source code, compiled manifest and latent hashes are saved.
The original sparse source files are read-only, not copied or rehashed in full.

## Diffusion

The denoiser is a time-modulated residual MLP (default width 512, six blocks),
using sinusoidal time embeddings, LayerNorm, SiLU, and scale/shift conditioning.
Training uses cosine DDPM noise, velocity prediction, AdamW, gradient clipping,
and an exponential moving average of denoiser weights. Validation uses fixed
noise and timesteps at each epoch for reproducible checkpoint selection.
Sampling uses DDIM with the EMA weights. Gaussian starting noise is an algorithmic
base distribution, not a Gaussian assumption on the final latent distribution.

The xVERSE baseline has unit-norm `z_bio`. Diffusion operates in standardized
Euclidean coordinates; after reversing standardization, generated states are
projected to the unit sphere before decoding when the source model normalizes
`z_bio`. Preprojection norms are saved to expose manifold mismatch. This is a
practical constrained-decoder baseline, not intrinsic spherical diffusion.

Generation uses the original biological `mu_bio` decoder with its learned
library size. NB uses the model's actual predicted dispersion, including
cell-by-gene dispersion; Poisson and ZINB sampling are also supported.
Counts are saved as compressed CSR shards in the original model gene order.
`--save-means` additionally saves dense expected-count shards (including the
zero-inflation factor for ZINB). No target donor embeddings are introduced.
Whole-dictionary output is model prediction, not validation of unmeasured genes.

## Usage

Run from the repository root with the dependencies in `requirements.txt`.
The existing DCC `SpaRest` environment can be used without installation.
Only load trusted PyTorch checkpoints and original source metadata.

```bash
python -m latent_diffusion extract \
  --xverse-repo /path/to/xverse-code \
  --checkpoint /path/to/best_model.pth \
  --compiled-root /path/to/compiled_train_v1_all \
  --gene-ids /path/to/training_gene_ids.txt \
  --output outputs/pilot/cache --train-cells 100000 --val-cells 10000

python -m latent_diffusion train \
  --cache outputs/pilot/cache --output outputs/pilot/diffusion

python -m latent_diffusion sample \
  --cache outputs/pilot/cache \
  --diffusion-checkpoint outputs/pilot/diffusion/best.pt \
  --xverse-repo /path/to/xverse-code --checkpoint /path/to/best_model.pth \
  --output outputs/pilot/samples --cells 1024
```

Use `--device cpu` for local tests. Extraction caps of zero process an entire
existing split. This can require substantial I/O and storage: 71 million
384-dimensional float32 latents alone are about 109 GB. Establish pilot
behavior before increasing the cap. Training currently starts fresh; `last.pt`
records optimizer state but automatic resume is not implemented.

DCC bounded integration pilot (4,096 train, 1,024 validation, 20 epochs):

```bash
sbatch bashfiles/01_latent_diffusion_pilot.sh
```

The launcher snapshots the original checkpoint and uses a fresh job-specific
output directory. It requests one generic GPU for at most one hour. It does not
submit a full-atlas training job. `TRAIN_CELLS`, `VAL_CELLS`, `EPOCHS`, `WIDTH`,
`DEPTH`, and the documented path variables may be overridden before submission.

## Validation and interpretation

```bash
XVERSE_REPO=/path/to/xverse-code python -m unittest discover -s tests -v
```

Tests cover diffusion velocity algebra, deterministic sampling, strict frozen
checkpoint loading, all eight likelihood/dispersion configurations, measured
zeros versus missing genes, training-only normalization, gene-order rejection,
and extraction-to-training-to-decoding without changing source weights.

`history.json` reports train/validation velocity MSE; this is not count
likelihood or evidence of biological quality. `diagnostics.json` compares latent
means, covariances, and sliced Wasserstein distance with validation embeddings,
including diagonal-Gaussian and empirical train-bootstrap baselines. For a
normalized checkpoint the Gaussian baseline receives the same sphere projection.
These are descriptive validation diagnostics, not independent test results.
Full count-distribution fidelity, cell-type coverage, independent-donor transfer,
and downstream augmentation require additional experiments.
