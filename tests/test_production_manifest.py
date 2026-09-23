from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from src.production_config import load_production_config
from src.production_manifest import build_cases, bundle_records, read_manifest, write_manifest


CONFIG = Path(__file__).parents[1] / "config/production_sweep_v1.json"


class ProductionManifestTests(unittest.TestCase):
    def setUp(self):
        self.config = load_production_config(CONFIG)

    def test_expected_case_counts_ids_and_shapes(self):
        group_a = build_cases(self.config, "A")
        self.assertEqual(len(group_a), 343)
        self.assertEqual(group_a[0].case_id, "a_r00_m00_g00")
        self.assertEqual(group_a[-1].case_id, "a_r06_m06_g06")
        self.assertTrue(all(len(case.axes["Ek"]) == 501 for case in group_a))
        for group in ("B_rho", "B_mu", "B_G"):
            cases = build_cases(self.config, group)
            self.assertEqual(len(cases), 49)
            self.assertTrue(all(tuple(len(case.axes[name]) for name in case.swept) == (501, 501) for case in cases))

    def test_manifest_generation_is_byte_deterministic(self):
        cases = build_cases(self.config, "B_rho")
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "one.jsonl", Path(directory) / "two.jsonl"
            write_manifest(first, cases)
            write_manifest(second, cases)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(len(read_manifest(first)), 49)

    def test_bundle_selection_preserves_case_order(self):
        records = tuple(case.to_json() for case in build_cases(self.config, "A"))
        selected = bundle_records(records, 2, 7)
        self.assertEqual([row["case_index"] for row in selected], list(range(14, 21)))


if __name__ == "__main__":
    unittest.main()
