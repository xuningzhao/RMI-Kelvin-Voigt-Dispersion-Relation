from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from src.production_config import load_production_config


CONFIG = Path(__file__).parents[1] / "config/production_sweep_v1.json"


class ProductionConfigTests(unittest.TestCase):
    def test_default_configuration_is_versioned_and_candidate_storage_is_off(self):
        config = load_production_config(CONFIG)
        self.assertEqual(config.data["schema_version"], 1)
        self.assertFalse(config.save_all_candidates)
        self.assertEqual(len(config.levels), 7)
        self.assertEqual(config.data["scans"]["Ek"]["samples"], 501)
        self.assertEqual(config.data["scans"]["contrast"]["samples"], 501)
        self.assertEqual(config.data["execution"]["case_workers"], 49)
        self.assertEqual(config.data["execution"]["point_workers"], 1)

    def test_configuration_hash_is_independent_of_json_formatting(self):
        first = load_production_config(CONFIG)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.json"
            path.write_text(json.dumps(first.data, indent=None), encoding="utf-8")
            second = load_production_config(path)
        self.assertEqual(first.sha256, second.sha256)


if __name__ == "__main__":
    unittest.main()
