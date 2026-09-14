# Project data locations

Large data are stored outside Git, following the xVERSE repository convention.

- DCC data root: `/hpc/group/xielab/xj58/xverse-m-data`.
- Science Perturb-multiome: `<data-root>/perturb_multiome_gse274113/`.
- Download scripts and source documentation: [data preparation](../data_preparation/README.md).

The tracked `data/` directory is a location index, not a second copy of the data.
Use `--data-root` to select another filesystem; no data symlink is required.
