# Previous xVERSE Version: Technical Summary

This document records only the technical designs that are present in the codebase at `/Users/xiaohui/LocalFiles/Codes/xverse-code`.

## 1. Model Architecture

The core model is `XVerseModel` in `main/utils_model.py`:

- The input is a dense count matrix in a fixed global gene universe, together with an `observed_mask` for each cell.
- `FiLMMaskEncoder` separately encodes `log1p(counts)` and the observed-gene mask.
- The mask encoder produces FiLM gamma and beta values that modulate the expression representation.
- The encoder outputs `z_bio_raw`; optional L2 normalization produces the public `z_bio`.
- `DenseExpressionDecoder` maps `z_bio` to gene logits. `library_head` predicts library size, and a softmax produces the nonnegative `mu_bio`.
- When sample embeddings are enabled, the decoder can additionally produce a sample-conditioned `mu`.

## 2. Observed and Random Masks

`FiLMMaskEncoder._apply_random_mask()` randomly hides already observed genes during training:

- For cells with few observed genes, it hides a subset while retaining at least one observed gene.
- For cells with many observed genes, it hides a stochastic fraction ranging from 0 to 70 percent.
- For sufficiently large panels, some views intersect the observed genes with a configured panel mask to simulate spatial panels.

The input `observed_mask` and the encoder's random mask have different roles: the former determines which genes contribute to the likelihood, while the latter determines which genes are visible in a particular encoder view.

## 3. Count Likelihoods

The code implements three reconstruction losses:

- Poisson negative log-likelihood.
- Negative Binomial negative log-likelihood.
- Zero-Inflated Negative Binomial negative log-likelihood.

NB dispersion supports six parameterizations: `global`, `gene`, `cell`, `factorized`, `cell_gene`, and `lowrank`. Numerical likelihood calculations run in fp32, with finite-value handling and bounds applied to counts, means, and dispersion values.

## 4. Auxiliary Objectives and Training

`pretrain_one_epoch()` performs two independent forward passes for each batch, creating two random-mask views. The training objective is:

$$
\mathcal{L} =
\lambda_{sample}\mathcal{L}_{recon}(\mu)
+ \lambda_{bio}\mathcal{L}_{recon}(\mu_{bio})
+ \lambda_{celltype}\mathcal{L}_{celltype}
+ \lambda_{contrast}\mathcal{L}_{contrast}.
$$

- Reconstruction losses for `mu` and `mu_bio` are computed only on observed genes.
- Cell types use cross-entropy with label smoothing; cells without labels are ignored.
- The projections from the two views use a bidirectional InfoNCE contrastive loss.
- Optimization uses Adam, AMP, gradient clipping, and `ReduceLROnPlateau`.
- Checkpoints are primarily selected using validation `loss_nb_bio`.

The code can replace cell-type classification with cosine-style matching to precomputed text embeddings through an auxiliary prediction head.

## 5. Data Loading

`main/data.py` provides two data paths.

### Legacy NPZ Blocks

`FastXVerseBatchDataset` lazily loads sparse CSR blocks, maps local gene lists to the global gene order, and uses an LRU cache to control memory. `SparseBatchCollator` converts sparse rows into a batch payload.

### Compiled `xverse_train_v1`

`CompiledShardDataset` uses memory-mapped shards and a manifest. It stores global cell ranges, sample IDs, cell-type IDs, and sparse gene/value arrays. When `source_pair_id` and source metadata are available, the loader recovers the complete measured panel and distinguishes measured zeros from unmeasured genes. Otherwise, it falls back to nonzero gene indices and prints a warning.

## 6. Sampling Strategy

`BalancedSampleSampler` and `CompiledBalancedSampler` sample by sample ID and support:

- A fixed number of sampled cells per sample.
- Reproducible ordering within each sample.
- Training-data-fraction ablations with a fixed number of training steps.
- Shard locality, active-shard limits, and reorder windows.

This design changes cell diversity in data-volume experiments while keeping optimizer steps and the validation split as consistent as possible.

## 7. Inference and Fine-tuning

`main/cli_xverse.py` provides embedding and generation tasks:

- Embedding writes `z_bio` to `adata.obsm["xVerse"]`.
- Generation writes `mu_bio` and NB samples to a new `.h5ad` file.

`XVerseFineTuneModel` reuses the base model's encoder, decoder, library head, and dispersion, and can add `sample_emb_ft`. It outputs `mu_bio` and, when the sample ID is valid, an additional sample-conditioned `mu`. The fine-tuning loop uses reconstruction and contrastive losses, a validation split, a scheduler, and early stopping.

The class defines `freeze_base_model()`, but the CLI fine-tuning loop does not call it automatically. By default, the reused base encoder and decoder parameters can therefore be updated together with the sample embedding when sample-specific fine-tuning is enabled.

## 8. Technical Principles Encoded in the Code

1. Use an explicit observed mask for different gene panels instead of treating panel absence as a true zero.
2. Use two randomly masked views to learn panel-robust representations.
3. Use `z_bio` to generate a biological mean and sample embeddings to express sample-specific variation.
4. Use count likelihoods rather than ordinary MSE to constrain expression generation.
5. Use sparse loading, memory mapping, and sample-balanced sampling to support large-scale training.
6. Let one model serve embedding, imputation, virtual-cell generation, and downstream augmentation.
