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
    raw/series/                # H5 matrices, annotated CSV and fragment XLSX
    raw/samples/               # ATAC fragments, indexes and peak annotations
    processed/                 # Reserved for future format conversion
    qc/                        # Reserved for future QC
    logs/
    download_plan.json         # Frozen URLs and expected byte counts
    download_manifest.json     # Progress, successful files and failures
```

As of 2026-09-14, GEO exposes 28 H5 matrices, one annotated metadata CSV,
one fragment-link spreadsheet and 42 sample files (14 fragment files with
indexes and peak annotations). The two H5 filename families are preserved;
their modalities must be determined from HDF5 features and sample metadata,
not inferred from the filename. Sample labels such as `rep1` must not be
assumed to identify independent donors.

GEO's approximately 45.4 GB `GSE274113_RAW.tar` contains the 42 sample files.
The script downloads those members directly from their official GSM URLs,
avoiding duplicate TAR storage and extraction. The spreadsheet is retained
as provenance; its external links are unnecessary because fragments are now
available directly from GEO. Exact totals are recorded in the download plan.

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
