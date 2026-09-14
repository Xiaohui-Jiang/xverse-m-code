#!/usr/bin/env bash
#SBATCH --job-name=xvm_science_data
#SBATCH --partition=biostat
#SBATCH --account=biostat
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/hpc/group/xielab/xj58/sbatch_output/xvm_science_data_%j.out
#SBATCH --error=/hpc/group/xielab/xj58/sbatch_output/xvm_science_data_%j.err

set -euo pipefail
REPO_ROOT=/hpc/group/xielab/xj58/xverse-m-code
DATA_ROOT=${XVERSE_M_DATA_ROOT:-/hpc/group/xielab/xj58/xverse-m-data}
PYTHON_BIN=${XVERSE_M_PYTHON:-/hpc/group/xielab/xj58/miniconda/envs/SpaRest/bin/python}
cd "$REPO_ROOT"
git rev-parse HEAD
"$PYTHON_BIN" -u data_preparation/00_download_perturb_multiome.py \
    --data-root "$DATA_ROOT" --workers 3 "$@"
