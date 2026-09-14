# Repository Instructions

## Language

All writing in this repository — code comments, docstrings, commit messages,
and any Markdown/documentation files — must be in English, regardless of the
language used to discuss the work.

## Local and remote repositories

- Use `/Users/xiaohui/LocalFiles/Codes/xverse-m-code` as the local repository.
- The corresponding HPC repository is `/hpc/group/xielab/xj58/xverse-m-code`.
- Connect to the HPC system with the SSH host alias `dcc`.

## Project scope

The current prototype fits latent diffusion to a frozen pretrained xVERSE RNA
encoder/decoder. Further project scope will be specified by the user.
The short summary under `xverse/` is background on the earlier xVERSE project;
it does not prescribe this repository's architecture, experiments, or roadmap.
The reference repository is `/Users/xiaohui/LocalFiles/Codes/xverse-code`
(locally) and `/hpc/group/xielab/xj58/xverse-code` (on DCC).

## Git and synchronization workflow

After modifying tracked repository files:

1. Validate the changed files locally with checks appropriate to the change.
2. Review `git status` and the relevant diff. Preserve unrelated or pre-existing
   user changes.
3. When making a commit, include only files belonging to the current task;
   never include unrelated working-tree changes.

Commits, pushes, and pulls should be performed when requested or when they are
useful for the current workflow. They are not mandatory after every edit.
Remote pull is optional rather than a completion requirement.

Ordinary untracked files such as `slurm-*.out` do not block normal Git
operations. When synchronizing the HPC checkout, avoid destructive operations
and do not overwrite tracked remote changes or resolve real conflicts without
the user's direction.

Use Git for tracked-code synchronization. Do not substitute `rsync` or `scp`
unless the user explicitly requests it.

## DCC Slurm partitions and accounts

The combinations below were checked against DCC with `sinfo`, `scontrol`,
`sacctmgr`, historical jobs for user `xj58`, and `sbatch --test-only` on
2026-09-07. Cluster availability and policy can change; re-run the diagnostic
commands at the end of this section if a submission is rejected.

### Recommended project partitions

| Partition | Account | GPU | Maximum time | Intended use |
| --- | --- | --- | --- | --- |
| `biostat` | `biostat` | No | 30 days | Default for preprocessing, downloads, CellTypist, PCA/Harmony, scIB, aggregation, and other CPU jobs. |
| `biostat-gpu` | `biostat` | Yes; heterogeneous RTX 5000 Ada and A6000 nodes | 7 days | Ordinary GPU training when a specific H200 is not required. Request `--gres=gpu:1`. |
| `h200-hp` | `xielab_h200` | NVIDIA H200, eight GPUs per node | 45 days | Default for large-model inference/training that needs an H200 or must be hardware-comparable. Request the H200 type explicitly. |
| `h200-hp` | `xielab_h200_r` | NVIDIA H200, eight GPUs per node | 45 days | Valid alternate/recharge H200 account; use only when requested or when its allocation is intended. |
| `scavenger-h200` | `scavenger-h200` | NVIDIA H200, eight GPUs per node | 1 day for user `xj58` | Opportunistic H200 fallback. The partition allows 7 days, but the current user association imposes `MaxWall=1-00:00:00`; it can queue much longer and should be treated as preemptible/requeueable. |
| `scavenger-gpu` | `xielab` | Yes; mixed 2080, RTX 5000 Ada, RTX 6000 Ada, A5000, and A6000 | 7 days | Opportunistic generic-GPU work that does not require identical hardware. Request `--gres=gpu:1`; do not assume the GPU model. |
| `scavenger` | `xielab` | No | 30 days | Opportunistic CPU fallback; default time is one day, so always set `-t`. |

Current association-level limits include 8 GPUs for `xielab_h200`, 2 GPUs
for `xielab_h200_r`, and 2 GPUs for `scavenger-h200`. The `biostat`, `xielab`,
and `chsi` associations each currently have group caps of 400 CPUs and 1.5 TB
memory. Large arrays can therefore show `AssocGrpMemLimit` even when nodes are
otherwise available; request realistic memory rather than the node maximum.

### General DCC fallback partitions

These combinations also passed `sbatch --test-only`, but prefer the project
partitions above for routine xVERSE work.

| Partition | Account | GPU | Maximum time | Notes |
| --- | --- | --- | --- | --- |
| `common` | `xielab` | Treat as CPU-only | 30 days | General shared CPU pool; it is DCC's default partition. Historical `xj58` jobs also used `common` with account `chsi`, but use that account only for work belonging to that allocation. |
| `gpu-common` | `xielab` | Mixed 2080, RTX 5000 Ada, and A5000 | 2 days | General shared GPU fallback; request `--gres=gpu:1`. |
| `interactive` | `xielab` | No | 1 day | Short interactive CPU sessions; not the default for batch pipelines. |

### Other verified account associations

These resources are available to `xj58`, but they are tied to other projects.
Do not charge xVERSE work to them unless the work belongs to that allocation or
the user explicitly requests it.

| Partition | Account | GPU | Maximum time | Notes |
| --- | --- | --- | --- | --- |
| `chsi` | `chsi` | No | 30 days | CHSI project CPU allocation. |
| `chsi-gpu` | `chsi` | NVIDIA 2080, four GPUs per observed node | 7 days | CHSI project GPU allocation; request `--gres=gpu:1`. |
| `common` | `chsi` | Treat as CPU-only | 30 days | General CPU fallback charged to the CHSI allocation; this combination has also appeared in historical `xj58` jobs. |

The `coursesf26cbb914` account association exists, but there is currently no
same-named partition and a test submission to that name is invalid. Do not use
it without current course-specific submission instructions. The `chsi-high`
partition does not currently accept the `chsi` account.

Do not use `h200-shared`: its allow-list does not include the xVERSE accounts.
The `h200ea` association currently has `gres/gpu=0` and is not a usable GPU
allocation. Visibility in `sinfo` alone does not imply permission to use a
lab-owned partition.

### SBATCH templates

Default CPU job:

```bash
#SBATCH -p biostat
#SBATCH -A biostat
#SBATCH -c 20
#SBATCH --mem=80G
#SBATCH -t 20:00:00
```

Biostat generic-GPU job:

```bash
#SBATCH -p biostat-gpu
#SBATCH -A biostat
#SBATCH --gres=gpu:1
#SBATCH -c 20
#SBATCH --mem=200G
#SBATCH -t 20:00:00
```

Hardware-controlled H200 job:

```bash
#SBATCH -p h200-hp
#SBATCH -A xielab_h200
#SBATCH --gres=gpu:h200:1
#SBATCH -c 20
#SBATCH --mem=200G
#SBATCH -t 20:00:00
```

Opportunistic generic-GPU job:

```bash
#SBATCH -p scavenger-gpu
#SBATCH -A xielab
#SBATCH --gres=gpu:1
#SBATCH -c 20
#SBATCH --mem=200G
#SBATCH -t 20:00:00
```

Opportunistic H200 job:

```bash
#SBATCH -p scavenger-h200
#SBATCH -A scavenger-h200
#SBATCH --gres=gpu:h200:1
#SBATCH -c 20
#SBATCH --mem=200G
#SBATCH -t 20:00:00
```

Use a CPU partition unless the program actually executes GPU kernels. For a
heterogeneous GPU partition, record the allocated hardware inside the job:

```bash
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
```

Useful read-only checks:

```bash
sacctmgr -nP show assoc where user=xj58 \
  format=Cluster,Account,Partition,QOS,GrpTRES
sinfo -h -o '%P|%a|%l|%D|%c|%m|%G' | sort -u
scontrol show partition <partition> -o
sbatch --test-only -p <partition> -A <account> [--gres=gpu[:type]:1] \
  -c 1 --mem=1G -t 00:01:00 --wrap=true
```

`sbatch --test-only` validates scheduling without submitting a real job.
