from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from src.production_storage import tile_bounds, validate_payload, write_payload


def record(index: int, *, rejected: bool = False):
    candidate = {"real": -1.0 - index, "imag": 0.5}
    verification = {
        "candidate": candidate, "domain_valid": True,
        "absolute_residual": 1e-12, "relative_residual": 1e-13,
        "verification_threshold": 1e-8, "mathematically_genuine": not rejected,
        "high_precision_refinement": None,
    }
    physical = {
        "candidate": candidate, "physically_admissible": not rejected,
        "q_plus": {"real": 1.0, "imag": 0.0},
        "q_minus": {"real": 1.0, "imag": 0.0},
    }
    result = {
        "counts": {"total_polynomial_candidates": 1, "physically_admissible_roots": int(not rejected)},
        "polynomial_candidates": [candidate], "verification_results": [verification],
        "physical_results": [physical],
    }
    return {"index": index, "coordinates": {"Ek": 0.0}, "parameters": {},
        "classification": "only_spurious_candidates" if rejected else "physically_admissible_roots",
        "runtime_seconds": 0.1, "error": None, "result": result}


class ProductionStorageTests(unittest.TestCase):
    def test_default_storage_omits_rejected_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            data, sidecar = Path(directory) / "data.npz", Path(directory) / "data.json"
            write_payload(data, sidecar, [record(0), record(1, rejected=True)],
                case_id="case", config_sha256="abc", save_all_candidates=False)
            self.assertTrue(validate_payload(data, sidecar, config_sha256="abc", case_id="case"))
            with np.load(data, allow_pickle=False) as loaded:
                self.assertEqual(len(loaded["root_real"]), 1)
                self.assertFalse(bool(loaded["save_all_candidates"]))

    def test_candidate_storage_can_be_enabled(self):
        with tempfile.TemporaryDirectory() as directory:
            data, sidecar = Path(directory) / "data.npz", Path(directory) / "data.json"
            write_payload(data, sidecar, [record(0), record(1, rejected=True)],
                case_id="case", config_sha256="abc", save_all_candidates=True)
            with np.load(data, allow_pickle=False) as loaded:
                self.assertEqual(len(loaded["root_real"]), 2)

    def test_corrupted_payload_is_not_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            data, sidecar = Path(directory) / "data.npz", Path(directory) / "data.json"
            write_payload(data, sidecar, [record(0)], case_id="case",
                config_sha256="abc", save_all_candidates=False)
            data.write_bytes(b"corrupt")
            self.assertFalse(validate_payload(data, sidecar, config_sha256="abc", case_id="case"))

    def test_failed_or_wrong_schema_payload_is_not_restart_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            data, sidecar = Path(directory) / "data.npz", Path(directory) / "data.json"
            write_payload(data, sidecar, [record(0)], case_id="case",
                config_sha256="abc", save_all_candidates=False)
            metadata = json.loads(sidecar.read_text(encoding="utf-8"))
            metadata["failure_count"] = 1
            sidecar.write_text(json.dumps(metadata), encoding="utf-8")
            self.assertFalse(validate_payload(data, sidecar, config_sha256="abc", case_id="case"))
            metadata["failure_count"] = 0
            metadata["schema_version"] = 999
            sidecar.write_text(json.dumps(metadata), encoding="utf-8")
            self.assertFalse(validate_payload(data, sidecar, config_sha256="abc", case_id="case"))

    def test_tiles_cover_grid_once(self):
        covered = np.zeros((7, 8), dtype=int)
        for r0, r1, c0, c1 in tile_bounds((7, 8), (3, 5)):
            covered[r0:r1, c0:c1] += 1
        np.testing.assert_array_equal(covered, np.ones((7, 8), dtype=int))

    def test_production_tile_layout_has_121_restart_units(self):
        bounds = list(tile_bounds((501, 501), (50, 50)))
        self.assertEqual(len(bounds), 121)
        self.assertEqual(bounds[0], (0, 50, 0, 50))
        self.assertEqual(bounds[-1], (500, 501, 500, 501))


if __name__ == "__main__":
    unittest.main()
