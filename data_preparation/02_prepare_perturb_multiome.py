#!/usr/bin/env python3
"""Adapt official Perturb-multiome loading/guide assignment to sparse AnnData.

Source: sankaranlab/perturb_multiome, commit 7455af0957569b7fc74809f71fa05cef282956d0,
notebooks 1, 2 (guide recoding), and 3a (replicate exclusions).
Copyright (c) 2024, jmartinrufino; BSD-2-Clause. See sources/perturb_multiome/.
This is a data-loading adaptation, not the complete Seurat/Signac/Mixscale pipeline.
"""
import argparse
import json
import os
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

SOURCE_COMMIT = '7455af0957569b7fc74809f71fa05cef282956d0'
SAMPLES = tuple(list(range(1, 11)) + [12, 13, 14, 16])


def read_10x(path):
    """Read feature x barcode CSC as barcode x feature CSR without densifying."""
    with h5py.File(path, 'r') as handle:
        group = handle['matrix']
        matrix = sparse.csc_matrix((group['data'][:], group['indices'][:],
                                    group['indptr'][:]), shape=tuple(group['shape'][:])).T.tocsr()
        features = pd.DataFrame({key: group['features'][key].asstr()[:]
                                 for key in ('id', 'name', 'feature_type')})
        barcodes = pd.Index(group['barcodes'].asstr()[:], name='barcode')
    if not barcodes.is_unique:
        raise ValueError(f'Duplicate barcodes in {path}')
    if matrix.shape != (len(barcodes), len(features)) or np.any(matrix.data < 0):
        raise ValueError(f'Invalid count matrix in {path}')
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    return matrix, barcodes, features


def read_metadata(path):
    # keep_default_na=False preserves meaningful literal annotations such as "NA".
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    empty = frame.eq('').all(axis=1)
    removed = int(empty.sum())
    frame = frame.loc[~empty].rename(columns={frame.columns[0]: 'cell_id'})
    if frame.cell_id.eq('').any() or frame.cell_id.duplicated().any():
        raise ValueError('Nonempty metadata records must have unique, nonempty cell IDs')
    return frame.set_index('cell_id'), removed


def assign_guides(counts, sequences, mapping):
    """Match R which.max: first maximum in feature order; zero totals unassigned."""
    total = np.asarray(counts.sum(axis=1)).ravel()
    maximum = np.asarray(counts.max(axis=1).toarray()).ravel()
    index = np.asarray(counts.argmax(axis=1)).ravel()
    seq = np.asarray(sequences, dtype=object)[index].copy()
    seq[total == 0] = ''
    result = pd.DataFrame({'perturbation': seq, 'nCount_cropseq': total,
                           'guide_max_count': maximum})
    result['perturbation_name'] = result.perturbation.map(mapping.perturbation_name).fillna('NP')
    result['target'] = result.perturbation.map(mapping.target).fillna('')
    result['guide_max_fraction'] = np.divide(maximum, total, out=np.zeros(len(total), dtype=float), where=total > 0)
    result['guide_max_ties'] = [int(np.count_nonzero(counts.data[counts.indptr[i]:counts.indptr[i+1]] == maximum[i]))
                               if total[i] else 0 for i in range(len(total))]
    return result


def prepare_sample(root, output, sample, metadata, mapping, selection, limit_cells):
    rep = f'rep{sample}'
    multi_path = root / 'raw/series' / f'GSE274113_{rep}_filtered_feature_bc_matrix.h5'
    guide_path = root / 'raw/series' / f'GSE274113_filtered_feature_bc_matrix_{sample}.h5'
    multi, barcodes, features = read_10x(multi_path)
    dialout, guide_barcodes, guide_features = read_10x(guide_path)
    shared = barcodes[barcodes.isin(guide_barcodes)]
    ids = pd.Index([f'{rep}_{barcode}' for barcode in shared], name='cell_id')
    guide_mask = guide_features.feature_type.eq('CRISPR Guide Capture').to_numpy()
    if not guide_mask.any():
        raise ValueError(f'No guide features in {guide_path}')
    guides = dialout[guide_barcodes.get_indexer(shared)][:, guide_mask].tocsr()
    guide_vars = guide_features.loc[guide_mask].copy()
    obs = assign_guides(guides, guide_vars['name'], mapping)
    obs.index = ids
    obs['barcode'] = shared.to_numpy()
    obs['replicate'] = rep
    obs['in_published_metadata'] = ids.isin(metadata.index)
    source = metadata.reindex(ids)
    for column in source:
        obs[f'published_{column}'] = source[column].fillna('')
    mismatches = {}
    for column in ('perturbation', 'perturbation_name', 'target'):
        compared = obs.in_published_metadata & source[column].ne('')
        mismatches[column] = int((obs.loc[compared, column] != source.loc[compared, column]).sum())
    keep = np.ones(len(obs), dtype=bool)
    if selection == 'published':
        keep = obs.in_published_metadata.to_numpy()
    elif selection == 'assigned':
        keep = obs.target.ne('').to_numpy()
    selected = np.flatnonzero(keep)
    before_limit = len(selected)
    if limit_cells:
        selected = selected[:limit_cells]
    if not len(selected):
        raise ValueError(f'No selected cells for {rep}')
    dest = output / rep
    dest.mkdir(parents=True, exist_ok=False)
    obs = obs.iloc[selected].copy()
    provenance = {'upstream_commit': SOURCE_COMMIT, 'cell_selection': selection,
                  'atac_feature_space': 'native_sample_peaks', 'normalized': False,
                  'input_multiome': str(multi_path), 'input_dialout': str(guide_path),
                  'limit_cells': limit_cells}
    matrix_rows = barcodes.get_indexer(shared[selected])
    for modality, feature_type in [('rna', 'Gene Expression'), ('atac', 'Peaks'), ('guides', 'CRISPR Guide Capture')]:
        if modality == 'guides':
            matrix, variables = guides[selected], guide_vars.copy()
        else:
            mask = features.feature_type.eq(feature_type).to_numpy()
            if not mask.any():
                raise ValueError(f'Missing {feature_type} in {multi_path}')
            matrix, variables = multi[matrix_rows][:, mask], features.loc[mask].copy()
        variables.index = pd.Index(variables['id'], name='feature_id')
        if not variables.index.is_unique:
            raise ValueError(f'Duplicate {modality} feature IDs')
        obj = ad.AnnData(X=matrix.tocsr(), obs=obs.copy(), var=variables)
        obj.uns['provenance'] = provenance
        temporary = dest / f'{modality}.tmp.h5ad'
        obj.write_h5ad(temporary, compression='gzip')
        os.replace(temporary, dest / f'{modality}.h5ad')
    report = {'sample': rep, 'multiome_barcodes': len(barcodes), 'dialout_barcodes': len(guide_barcodes),
              'intersected_cells': len(shared), 'selected_before_limit': before_limit,
              'written_cells': len(obs), 'published_label_mismatches_before_selection': mismatches,
              'published_cells_missing_from_intersection': int((~metadata.index[metadata.index.str.startswith(rep + '_')].isin(ids)).sum()),
              'target_counts': obs.target.value_counts().to_dict(), 'provenance': provenance}
    (dest / 'qc.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=Path('/hpc/group/xielab/xj58/xverse-m-data'))
    parser.add_argument('--output', type=Path, help='New output directory; existing paths are refused')
    parser.add_argument('--samples', type=int, nargs='+', default=list(SAMPLES), choices=SAMPLES)
    parser.add_argument('--cell-selection', choices=['published', 'assigned', 'intersection'], default='published')
    parser.add_argument('--limit-cells', type=int, default=0, help='First N selected cells per sample for smoke tests; 0 means all')
    args = parser.parse_args()
    if args.limit_cells < 0 or len(set(args.samples)) != len(args.samples):
        parser.error('Use a nonnegative cell limit and unique sample numbers')
    root = args.data_root / 'perturb_multiome_gse274113'
    output = args.output or root / 'processed/python_native'
    output.mkdir(parents=True, exist_ok=False)
    metadata, removed = read_metadata(root / 'raw/series/GSE274113_annotated_metadata.csv.gz')
    mapping = pd.read_csv(Path(__file__).parent / 'sources/perturb_multiome/guide_map.tsv', sep='\t').set_index('sequence')
    reports = []
    for sample in args.samples:
        report = prepare_sample(root, output, sample, metadata, mapping, args.cell_selection, args.limit_cells)
        reports.append(report)
        print(json.dumps(report), flush=True)
    (output / 'manifest.json').write_text(json.dumps({'status': 'complete', 'empty_metadata_records_removed': removed,
        'published_metadata_cells': len(metadata), 'samples': reports}, indent=2) + '\n')


if __name__ == '__main__':
    main()
