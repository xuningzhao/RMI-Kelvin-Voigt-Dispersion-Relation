from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from src.production_config import load_production_config
from src.production_manifest import build_cases
from src.production_runner import run_case


BASE_CONFIG = Path(__file__).parents[1] / "config/production_sweep_v1.json"


def small_config(directory: Path):
    data = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    data["scans"]["Ek"]["samples"] = 3
    data["scans"]["contrast"]["samples"] = 3
    data["execution"]["group_b_tile_shape"] = [2, 2]
    path = directory / "config.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return load_production_config(path)


class ProductionRunnerTests(unittest.TestCase):
    def test_manifest_configuration_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            config = small_config(directory)
            case = build_cases(config, "A")[0].to_json()
            case["config_sha256"] = "wrong"
            with self.assertRaisesRegex(ValueError, "hashes differ"):
                run_case(case, config, directory)

    def test_group_a_is_single_payload_and_restarts(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            config = small_config(directory)
            case = build_cases(config, "A")[171].to_json()
            first = run_case(case, config, directory)
            second = run_case(case, config, directory)
            self.assertEqual(first["completed_units"], 1)
            self.assertEqual(second["skipped_units"], 1)
            output = directory / case["output_path"]
            self.assertTrue((output / "spectrum.npz").is_file())
            self.assertFalse((output / "tiles").exists())

    def test_group_b_tiles_restart_and_corruption_recovery(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            config = small_config(directory)
            case = build_cases(config, "B_G")[24].to_json()
            first = run_case(case, config, directory)
            output = directory / case["output_path"]
            payloads = sorted((output / "tiles").glob("*.npz"))
            self.assertEqual(first["expected_units"], 4)
            self.assertEqual(len(payloads), 4)
            payloads[0].write_bytes(b"corrupt")
            second = run_case(case, config, directory)
            self.assertEqual(second["completed_units"], 1)
            self.assertEqual(second["skipped_units"], 3)

    def test_partial_group_b_run_is_restartable(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            config = small_config(directory)
            case = build_cases(config, "B_G")[24].to_json()
            partial = run_case(case, config, directory, maximum_new_units=2)
            self.assertFalse(partial["complete"])
            self.assertEqual(partial["completed_units"], 2)
            resumed = run_case(case, config, directory)
            self.assertTrue(resumed["complete"])
            self.assertEqual(resumed["skipped_units"], 2)
            self.assertEqual(resumed["completed_units"], 2)

    def test_point_failures_are_recorded_without_terminating_case(self):
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            config = small_config(directory)
            case = build_cases(config, "A")[171].to_json()
            with patch("src.production_runner.analyze_parameter_point", side_effect=RuntimeError("forced")):
                status = run_case(case, config, directory)
            self.assertEqual(status["status"], "complete_with_failures")
            self.assertEqual(status["failure_count"], 3)
            failures = json.loads((directory / case["output_path"] / "failures.json").read_text())
            self.assertEqual(len(failures["failures"]), 3)
            retried = run_case(case, config, directory)
            self.assertEqual(retried["completed_units"], 1)
            self.assertEqual(retried["failure_count"], 0)
            self.assertFalse((directory / case["output_path"] / "failures.json").exists())


if __name__ == "__main__":
    unittest.main()
