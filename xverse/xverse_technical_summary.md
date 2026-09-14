# xVERSE: brief background

Based on the xVERSE repository at commit
`fab26c3d7dd14486a24735edeb6ae66cdac36b22` (reviewed 2026-09-14).
This note summarizes earlier work only; it does not define the new project.

## Purpose

xVERSE combines transferable transcriptomic representations with
template-conditioned virtual-cell generation. A real cell supplies the input
profile; pretrained parameters help predict a count distribution from which
new expression profiles are sampled. Applications include gene imputation and
augmentation for analyses with limited observed cells.

## Model

- A mask-aware autoencoder separately processes expression and the measured-gene
  mask, combines them through FiLM, and produces the biological embedding
  `z_bio` (384 dimensions in the baseline).
- Observation masks distinguish measured zeros from unmeasured genes. Paired
  random gene masking and contrastive learning encourage panel robustness.
- A shared decoder predicts gene proportions and library size, producing
  sample-independent `mu_bio` or sample-conditioned `mu`.
- Poisson, negative-binomial, and zero-inflated negative-binomial likelihoods
  are supported, with several dispersion parameterizations. Generation samples
  counts from the predicted distribution; the encoder has no stochastic latent
  prior.
- Pretraining combines biological and sample-conditioned reconstruction,
  panel contrast, and cell-type text alignment. Fine-tuning adapts to target
  samples and supports the base model's likelihood and dispersion heads.

## Evaluation

Current Figure 2 protocols cover representation quality and measured-panel
population generation. The embedding benchmark uses five single-cell/spatial
dataset families and official scIB scores. Generation compares zero-shot,
fine-tuned, and from-scratch xVERSE with scVI, using donor-balanced A/B splits,
fixed target-cell budgets, and marginal and multivariate distribution metrics.
The current outside-panel pipeline still awaits finalized ground-truth pairing.

Earlier manuscript work also explored spatial imputation, rare-cell and DEG
recovery, and cross-modality prediction. Historical results and future plans
should be distinguished from the current numerical-validation protocols.
Virtual cells are model-derived augmentations, not independent biological
replicates.

## Reference

Source: `/Users/xiaohui/LocalFiles/Codes/xverse-code`.
Read `main/model-structure.md` for architecture and the experiment-specific
READMEs under `figure_codes/figure2/` for current evaluation details.
Some parent documents retain older protocols. The legacy CLI also has narrower
checkpoint-loading and generation-export support than the model itself.
