from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from src.production_index import build_global_index


class ProductionIndexTests(unittest.TestCase):
    def test_index_collects_case_status_without_shared_task_writes(self):
        with tempfile.TemporaryDirectory() as directory_name:
            root = Path(directory_name)
            for group, case in (("group_a", "a_r00_m00_g00"), ("group_b_rho", "brho_m00_g00")):
                path = root / group / case
                path.mkdir(parents=True)
                (path / "status.json").write_text(json.dumps({
                    "case_id": case, "group": group, "status": "complete",
                    "complete": True, "runtime_seconds": 1.0,
                    "output_path": str(path), "failure_count": 0,
                }), encoding="utf-8")
            index = build_global_index(root)
            self.assertEqual(index["case_count"], 2)
            self.assertEqual(index["complete_case_count"], 2)
            self.assertTrue((root / "global_index.json").is_file())


if __name__ == "__main__":
    unittest.main()
