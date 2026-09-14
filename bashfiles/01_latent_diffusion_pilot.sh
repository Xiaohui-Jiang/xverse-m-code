#!/bin/bash
#SBATCH --job-name=xverse_m_ldm_pilot
#SBATCH --partition=biostat-gpu
#SBATCH --account=biostat
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/hpc/group/xielab/xj58/sbatch_output/%x_%j.out
#SBATCH --error=/hpc/group/xielab/xj58/sbatch_output/%x_%j.err
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/hpc/group/xielab/xj58/xverse-m-code}"
XVERSE_REPO="${XVERSE_REPO:-/hpc/group/xielab/xj58/xverse-code}"
PYTHON="${PYTHON:-/hpc/group/xielab/xj58/miniconda/envs/SpaRest/bin/python}"
COMPILED_ROOT="${COMPILED_ROOT:-/hpc/group/xielab/xj58/xVerseAtlas/compiled_train_v1_all}"
GENE_IDS="${GENE_IDS:-/hpc/group/xielab/xj58/xVerseAtlas/npz_tissue_dataset_donor/ensg_keys_high_quality.txt}"
BASE_CKPT="${BASE_CKPT:-/hpc/group/xielab/xj58/pretrain_model_celltype/xverse_count_ablation0820/nb_cell_gene/best_model.pth}"
OUTPUT_ROOT="${OUTPUT_ROOT:-/hpc/group/xielab/xj58/xverse-m-results/latent_diffusion/pilot_${SLURM_JOB_ID}}"
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
cd "$REPO_ROOT"
mkdir -p "$(dirname "$OUTPUT_ROOT")"
mkdir "$OUTPUT_ROOT"
git rev-parse HEAD > "$OUTPUT_ROOT/code_commit.txt"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
# Snapshot trusted weights so an upstream best-checkpoint update cannot change
# the encoder between extraction and sampling. The upstream file is untouched.
cp "$BASE_CKPT" "$OUTPUT_ROOT/xverse_frozen.pth"
"$PYTHON" -m latent_diffusion extract \
  --xverse-repo "$XVERSE_REPO" --checkpoint "$OUTPUT_ROOT/xverse_frozen.pth" \
  --compiled-root "$COMPILED_ROOT" --gene-ids "$GENE_IDS" \
  --output "$OUTPUT_ROOT/cache" --device cuda --batch-size 256 --workers 0 \
  --train-cells "${TRAIN_CELLS:-4096}" --val-cells "${VAL_CELLS:-1024}"
"$PYTHON" -m latent_diffusion train \
  --cache "$OUTPUT_ROOT/cache" --output "$OUTPUT_ROOT/diffusion" \
  --device cuda --batch-size 256 --width "${WIDTH:-256}" --depth "${DEPTH:-4}" \
  --epochs "${EPOCHS:-20}" --ema-decay 0.95 --patience 0
"$PYTHON" -m latent_diffusion sample \
  --cache "$OUTPUT_ROOT/cache" --diffusion-checkpoint "$OUTPUT_ROOT/diffusion/best.pt" \
  --xverse-repo "$XVERSE_REPO" --checkpoint "$OUTPUT_ROOT/xverse_frozen.pth" \
  --output "$OUTPUT_ROOT/samples" --device cuda --batch-size 128 \
  --cells 256 --sampling-steps 100
