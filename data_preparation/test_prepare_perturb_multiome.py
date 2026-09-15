"""Regression checks for alignment and the official guide-assignment convention."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

spec = importlib.util.spec_from_file_location('prepare', Path(__file__).with_name('02_prepare_perturb_multiome.py'))
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class PreparationTests(unittest.TestCase):
    def test_ties_zero_and_control_mapping(self):
        mapping = pd.DataFrame({'perturbation_name': ['AAVS1_1', 'GATA1_674'],
                                'target': ['NT', 'GATA1']}, index=['AAA', 'CCC'])
        result = prepare.assign_guides(sparse.csr_matrix([[3, 3], [0, 0], [1, 4]]), ['AAA', 'CCC'], mapping)
        self.assertEqual(result.perturbation.tolist(), ['AAA', '', 'CCC'])
        self.assertEqual(result.target.tolist(), ['NT', '', 'GATA1'])
        self.assertEqual(result.guide_max_ties.tolist(), [2, 0, 1])
        np.testing.assert_allclose(result.guide_max_fraction, [0.5, 0, 0.8])

    def test_metadata_empty_records_and_literal_na(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'metadata.csv'
            path.write_text(',target,annotation_simplified\n,,\nrep1_A,NT,NA\nrep1_B,,\n')
            frame, removed = prepare.read_metadata(path)
            self.assertEqual(removed, 1)
            self.assertEqual(frame.index.tolist(), ['rep1_A', 'rep1_B'])
            self.assertEqual(frame.loc['rep1_A', 'annotation_simplified'], 'NA')
            path.write_text(',target\nrep1_A,NT\nrep1_A,NT\n')
            with self.assertRaises(ValueError):
                prepare.read_metadata(path)

    def test_sample_alignment_and_sparse_roundtrip(self):
        import h5py
        import anndata as ad
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / 'raw/series'
            raw.mkdir(parents=True)
            def write(name, counts, barcodes, ids, types):
                matrix = sparse.csc_matrix(counts)
                with h5py.File(raw / name, 'w') as handle:
                    group = handle.create_group('matrix')
                    for key, value in [('data', matrix.data), ('indices', matrix.indices), ('indptr', matrix.indptr), ('shape', matrix.shape)]:
                        group[key] = value
                    group['barcodes'] = np.array(barcodes, dtype='S')
                    features = group.create_group('features')
                    for key, value in [('id', ids), ('name', ids), ('feature_type', types)]:
                        features[key] = np.array(value, dtype='S')
            write('GSE274113_rep1_filtered_feature_bc_matrix.h5', [[1, 2, 3], [4, 5, 6]],
                  ['A', 'B', 'C'], ['gene', 'chr1:0-100'], ['Gene Expression', 'Peaks'])
            write('GSE274113_filtered_feature_bc_matrix_1.h5', [[9, 8, 7]],
                  ['C', 'A', 'D'], ['AAA'], ['CRISPR Guide Capture'])
            metadata = pd.DataFrame({'perturbation': ['AAA', 'AAA'], 'perturbation_name': ['NT_1', 'NT_1'],
                                     'target': ['NT', 'NT']}, index=['rep1_C', 'rep1_A'])
            mapping = pd.DataFrame({'perturbation_name': ['NT_1'], 'target': ['NT']}, index=['AAA'])
            out = root / 'out'
            report = prepare.prepare_sample(root, out, 1, metadata, mapping, 'published', 0)
            rna, atac, guides = [ad.read_h5ad(out / 'rep1' / f'{name}.h5ad') for name in ['rna', 'atac', 'guides']]
            self.assertEqual(rna.obs_names.tolist(), ['rep1_A', 'rep1_C'])
            self.assertTrue(rna.obs_names.equals(atac.obs_names) and rna.obs_names.equals(guides.obs_names))
            np.testing.assert_array_equal(rna.X.toarray().ravel(), [1, 3])
            np.testing.assert_array_equal(atac.X.toarray().ravel(), [4, 6])
            np.testing.assert_array_equal(guides.X.toarray().ravel(), [8, 9])
            self.assertEqual(report['published_label_mismatches_before_selection']['target'], 0)


if __name__ == '__main__':
    unittest.main()
