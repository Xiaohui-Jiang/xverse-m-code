"""Small offline integrity tests for resumable data acquisition."""

import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "downloader", Path(__file__).with_name("00_download_perturb_multiome.py"))
downloader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(downloader)


class DownloadIntegrityTest(unittest.TestCase):
    def test_completed_partial_is_validated_and_existing_corruption_rejected(self):
        payload = gzip.compress(b"chr1\t1\t5\tcell\t1\n")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = dict(path="raw/fragments.tsv.gz", bytes=len(payload), url="https://example.org/data")
            target = root / item["path"]
            target.parent.mkdir()
            target.with_name(target.name + ".part").write_bytes(payload)
            with patch.object(downloader.subprocess, "run") as transfer:
                first = downloader.download_one(root, item, 1)
                second = downloader.download_one(root, item, 1)
                transfer.assert_not_called()
            self.assertEqual(first["sha256"], second["sha256"])
            target.write_bytes(b"x" * len(payload))
            with self.assertRaisesRegex(ValueError, "receipt mismatch"):
                downloader.download_one(root, item, 1)

    def test_corrupt_full_partial_is_not_published(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            item = dict(path="bad.gz", bytes=7, url="https://example.org/data")
            (root / "bad.gz.part").write_bytes(b"corrupt")
            with self.assertRaises(gzip.BadGzipFile):
                downloader.download_one(root, item, 1)
            self.assertFalse((root / "bad.gz").exists())


if __name__ == "__main__":
    unittest.main()
