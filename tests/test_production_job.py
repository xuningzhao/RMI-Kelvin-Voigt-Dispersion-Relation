from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.run_production_bundle import _run_lane
from src.production_config import load_production_config
from src.production_manifest import build_cases


BASE_CONFIG = Path(__file__).parents[1] / "config/production_sweep_v1.json"


class ProductionJobTests(unittest.TestCase):
    def test_lane_runs_cases_sequentially_and_isolates_case_failure(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            data = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
            data["scans"]["Ek"]["samples"] = 3
            data["scans"]["contrast"]["samples"] = 3
            data["execution"]["case_workers"] = 1
            data["execution"]["group_a_cases_per_worker"] = 2
            config_path = directory / "config.json"
            config_path.write_text(json.dumps(data), encoding="utf-8")
            config = load_production_config(config_path)
            cases = tuple(case.to_json() for case in build_cases(config, "A")[:2])
            complete = {"case_id": cases[1]["case_id"], "complete": True}
            with patch("scripts.run_production_bundle.run_case",
                side_effect=(RuntimeError("forced"), complete)) as mocked:
                statuses = _run_lane((cases, config_path, directory, None))
            self.assertEqual(mocked.call_count, 2)
            self.assertFalse(statuses[0]["complete"])
            self.assertTrue(statuses[1]["complete"])
            status_path = directory / cases[0]["output_path"] / "status.json"
            self.assertEqual(json.loads(status_path.read_text())["status"], "case_worker_failed")


if __name__ == "__main__":
    unittest.main()
