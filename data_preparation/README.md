# Science Perturb-multiome data acquisition

## Source and scope

Martin-Rufino et al., *Transcription factor networks disproportionately enrich
for heritability of blood cell phenotypes*, Science (2025),
[DOI: 10.1126/science.ads7951](https://doi.org/10.1126/science.ads7951).
The paired single-cell RNA/ATAC study is
[GEO GSE274113](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274113).
The experiment profiles CRISPR TF perturbations in primary human hematopoietic
cells at days 7, 9, 11 and 14. RNA and ATAC are paired within a cell; control
and perturbed cells are not longitudinal measurements of the same cell.

This downloader obtains the original GEO supplementary release, not the
MultiFlow-derived H5MU. It retains the original matrices, annotations and
ATAC fragments needed to build an EpiAgent-compatible cCRE matrix later.
It does not download SRA FASTQs or the separate GSE274110 bulk experiment.
No filtering, normalization, cell matching or train/test split is performed.

The layout follows xVERSE's `figure_codes/figure2/00_data_preparation/`:
numbered Python acquisition scripts, separate `bashfiles/` Slurm launchers,
an external `--data-root`, and per-dataset raw/processed/QC locations.

## Paths

```text
/hpc/group/xielab/xj58/xverse-m-data/
  perturb_multiome_gse274113/
    metadata/                  # Official HTML listing, filelist and family SOFT
    raw/series/                # H5 matrices, annotated CSV and fragment inventory XLSX
    raw/samples/               # ATAC fragments, indexes and peak annotations
    processed/                 # Reserved for future format conversion
    qc/                        # Reserved for future QC
    logs/
    download_plan.json         # Frozen URLs and expected byte counts
    download_manifest.json     # Progress, successful files and failures
```

As of 2026-09-14, GEO exposes 28 H5 matrices, one annotated metadata CSV,
one fragment inventory spreadsheet and 42 sample files (14 fragment files with
indexes and peak annotations). The two H5 filename families are preserved;
their modalities must be determined from HDF5 features and sample metadata,
not inferred from the filename. Sample labels such as `rep1` must not be
assumed to identify independent donors.

Inspection of the downloaded rep1 files confirms:

- `GSE274113_filtered_feature_bc_matrix_1.h5`: `Gene Expression` and
  `CRISPR Guide Capture`, 36,664 features by 14,525 barcodes.
- `GSE274113_rep1_filtered_feature_bc_matrix.h5`: `Gene Expression` and
  `Peaks`, 194,891 features by 13,168 barcodes.

Do not concatenate these files as independent cells or assume identical
barcode sets. Join by sample-prefixed barcode and inspect `feature_type`.

The original annotated CSV contains **137,604 nonempty, unique cell IDs**
and **137,604 all-empty records**. The raw file is preserved unchanged;
future preprocessing must explicitly discard all-empty records. Its unnamed
first column contains IDs such as `rep1_AAACAGCCAACAGCCT-1`. Important fields
are `replicate`, `Timepoint`, `perturbation` (guide sequence),
`perturbation_name` (guide name), `target`, `new_CellType` and
`annotation_simplified`. The target field has 19 TFs plus `NT`; do not infer
the experimental control type from that aggregate label without checking
guide-level annotations. There are 14 sample labels, not an established
count of biological donors.

| Timepoint | Annotated cells |
| --- | ---: |
| day 7 | 43,723 |
| day 9 | 50,116 |
| day 11 | 20,791 |
| day 14 | 22,974 |

These counts describe the deposited annotation, not a new QC selection.

GEO's approximately 45.4 GB `GSE274113_RAW.tar` contains the 42 sample files.
The script downloads those members directly from their official GSM URLs,
avoiding duplicate TAR storage and extraction. The spreadsheet is retained
as provenance and contains GSM accessions and placeholder filenames. The
actual fragments are available directly from GEO. Exact totals are in the plan.
The 2026-09-14 plan contains 72 data files totaling 49,422,011,482 bytes
(49.42 GB decimal), excluding the small provenance snapshots and receipts.

## Run

Requires Python 3.9+, `h5py`, and a curl version supporting
`--retry-all-errors` (7.71+). The DCC SpaRest environment supplies Python/h5py.
The download uses CPU only; it does not modify model checkpoints.

```bash
cd /hpc/group/xielab/xj58/xverse-m-code
PYTHON=/hpc/group/xielab/xj58/miniconda/envs/SpaRest/bin/python
$PYTHON data_preparation/00_download_perturb_multiome.py --plan-only
mkdir -p /hpc/group/xielab/xj58/sbatch_output
sbatch bashfiles/00_download_perturb_multiome.sh
```

When biostat CPU nodes are occupied, use the authorized fallback:
`sbatch -p common -A xielab bashfiles/00_download_perturb_multiome.sh`.

Offline integrity checks:
`python -m unittest discover -s data_preparation -p 'test_*.py'`.

Direct execution on another machine:

```bash
python data_preparation/00_download_perturb_multiome.py \
  --data-root /path/to/xverse-m-data --workers 3
```

The launcher accepts `XVERSE_M_DATA_ROOT` and `XVERSE_M_PYTHON` overrides.
Allow at least 60 GB free storage for this release. The script snapshots GEO
metadata on the first run and reuses the plan on subsequent runs. Use
`--refresh-plan` only when deliberately adopting a newer upstream release.

## Integrity and resuming

- HTTP transfers retry and resume from `.part` files.
- Sizes are checked against exact GEO filelist bytes (sample files) or HTTP
  Content-Length (series files).
- Gzip files are fully decompressed for CRC validation without retaining
  uncompressed copies. H5 files receive basic sparse-matrix structural checks;
  the spreadsheet receives ZIP integrity validation.
- Each successful file has a SHA-256 receipt. These are locally calculated
  checksums, not publisher-supplied checksums. Reruns verify existing receipts.
- One writer per dataset is enforced with a filesystem lock. The manifest is
  replaced atomically after each completed file; `status=complete` appears
  only after every planned file succeeds.
- A mismatched existing file causes an error rather than silent replacement.
  Inspect and move it aside before retrying. A corrupt full-size `.part` must
  likewise be moved aside before retrying. Interrupted partial files otherwise
  resume automatically.

Monitor Slurm output in `/hpc/group/xielab/xj58/sbatch_output/xvm_science_data_JOBID.out`
and `.err`; check `download_manifest.json` for verified-file progress.
Pretraining overlap for xVERSE and EpiAgent has not been audited by this download.

## Upstream integrity issue discovered on 2026-09-14

DCC acquisition job `55461528` ran on `common` / `xielab` for 12 minutes
14 seconds. All 72 planned files reached their advertised byte counts;
71 passed validation (including all 28 H5 matrices, metadata, and 13 of 14
fragment files). The job intentionally exited with code 1 and the manifest
has `status=failed` because the remaining fragment file is corrupt upstream.
The other verified files are available for use. No QC or model fitting was run.

`GSM8443612_rep2_atac_fragments.tsv.gz` has the exact GEO-advertised length
(3,052,699,648 bytes), but full gzip validation fails with an unexpected EOF.
An independent HTTP Range request reproduced the identical last 128 KiB.
The last BGZF block declares 14,165 bytes but only 1,663 remain in the file.
This indicates truncation in the currently served source, not simply a short
local transfer. The downloader deliberately retains this as `.part` and
reports failure rather than treating it as usable. Existing rep2 H5 matrices
are separate assets and are not invalidated by this fragment-file defect.

Reproduce the tail diagnostic (JSON is saved under `metadata/`):

```bash
python data_preparation/01_check_fragment_tail.py
```

Do not repair this by appending a gzip footer, silently discarding the final
block, or using the truncated data as a complete sample. A corrected source
or reprocessing of the corresponding raw sequencing reads would be required
for a fully validated rep2 fragment input. Other files continue downloading
even when one file fails.
