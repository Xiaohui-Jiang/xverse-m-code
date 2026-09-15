#!/bin/bash
#SBATCH --job-name=xvm_prepare_multiome
#SBATCH -p biostat
#SBATCH -A biostat
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /hpc/group/xielab/xj58/sbatch_output/xvm_prepare_multiome_%j.out
#SBATCH -e /hpc/group/xielab/xj58/sbatch_output/xvm_prepare_multiome_%j.err
set -euo pipefail
cd /hpc/group/xielab/xj58/xverse-m-code
PYTHON=${XVERSE_M_PYTHON:-/hpc/group/xielab/xj58/miniconda/envs/SpaRest/bin/python}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
"$PYTHON" -m unittest discover -s data_preparation -p 'test_*.py'
"$PYTHON" data_preparation/02_prepare_perturb_multiome.py \
  --data-root "${XVERSE_M_DATA_ROOT:-/hpc/group/xielab/xj58/xverse-m-data}" "$@"
