# Upstream attribution and adaptation boundary

Source: Martin-Rufino et al., *Transcription factor networks disproportionately
enrich for heritability of blood cell phenotypes*, Science (2025).
[Paper](https://doi.org/10.1126/science.ads7951) ·
[Official repository](https://github.com/sankaranlab/perturb_multiome)

Pinned revision: `7455af0957569b7fc74809f71fa05cef282956d0`.
The original BSD 2-Clause notice is retained in [LICENSE.txt](LICENSE.txt).
Copyright (c) 2024, jmartinrufino. Python adaptations were added by the xVERSE-M
project in 2026. The guide mapping is extracted verbatim from notebook 2's
`recoding_vector` and `recoding_vector_target` (63 sequences).

| Official source | Python adaptation |
| --- | --- |
| [1: RNA object](https://github.com/sankaranlab/perturb_multiome/blob/7455af0957569b7fc74809f71fa05cef282956d0/1_GITHUB_processing_RNA_object.ipynb) | `02_prepare_perturb_multiome.py`: feature-type selection, per-sample barcode intersection, sample prefixes, maximum-count guide assignment, first-feature tie handling, zero-count unassignment |
| [2: multiome object](https://github.com/sankaranlab/perturb_multiome/blob/7455af0957569b7fc74809f71fa05cef282956d0/2_GITHUB_processing_perturb_multiome_object.ipynb) | Guide sequence/name/target recoding in `guide_map.tsv`; AAVS1 and NT guides share target NT; `--cell-selection assigned` keeps recognized targets |
| [3a: TF-sensitive element computation](https://github.com/sankaranlab/perturb_multiome/blob/7455af0957569b7fc74809f71fa05cef282956d0/3_TF_sensitive_element_computation/3a_GITHUB_TF_sensitive_element_computation.ipynb) | Exclusion of rep11 and rep15, described upstream as wetting failures; published annotations are joined by cell ID |

## Deliberate differences and steps not ported

This is a runnable port of the data-loading and guide-assignment steps, not a
full reproduction of the paper's analysis. The GEO release has 14 samples;
upstream notebooks initially load 16 and exclude rep11/rep15 in notebook 3a.
GEO filenames are resolved to their observed modalities rather than copied
from the author's local paths.

- Outputs are three aligned sparse AnnData files per sample, with raw counts
  in `X`, stable feature IDs in `var_names`, and symbols/sequences in `var.name`.
  The dialout's Gene Expression assay is not a second population of cells.
- Default selection uses the deposited annotated cell IDs. `intersection`
  reproduces notebook 1's barcode selection; `assigned` adds notebook 2's
  recognized-guide selection **without its ATAC QC**. All selection modes omit
  rep11/rep15. These modes are not equivalent to rebuilding final Seurat QC.
- All-empty CSV records are removed. Individual missing fields are retained;
  literal `NA` cell-type annotations are preserved. Recomputed guide labels
  and `published_*` annotations remain separate and disagreements are counted.
- ATAC uses each sample's native deposited peaks. Notebook 2 instead reduces
  peaks across samples (20 < width < 10,000), recounts fragments with Signac,
  and applies `min.cells=10`, `min.features=200`. **That recount/QC is not
  implemented here. Do not concatenate these ATAC matrices as a common
  feature space or call them EpiAgent cCRE inputs.** The rep2 source fragment
  file is currently truncated, so a complete fragment-based recreation also
  needs a corrected input. Native rep2 H5 counts remain usable.
- Notebook 3a's TF-IDF/LSI, RNA PCA, clustering, nearest-control perturbation
  signatures and Mixscale scores are not implemented. Its ATAC scoring relies
  on a modified Mixscale installation with an empty GitHub repository string;
  this adapter does not substitute an unverified Python approximation.
- No doublet removal, uncertainty calibration, perturbation-effect estimation,
  normalization, pretraining-overlap audit, or train/test split is implied.

The extra guide fraction/tie diagnostics, input provenance, metadata comparison
and atomic output writes are additions in this repository, not paper methods.
