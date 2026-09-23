"""Tests for Sprint 7 sweep-summary extraction."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from src.parameter_sweep import (
    CLASSIFICATIONS,
    cartesian_grid_definition,
    explicit_sweep_definition,
    run_parameter_sweep,
)
from src.sweep_summary import (
    CLASSIFICATION_CODES,
    SweepSummaryError,
    read_sweep_summary,
    write_summary_outputs,
)


GENERIC = {"Arho": 0.2, "Amu": -0.3, "AG": 0.4, "Lambda": 1.7}
ONLY_SPURIOUS = {"Arho": 0.5, "Amu": 0.5, "AG": 0.0, "Lambda": 81.0}
OUTSIDE_EK = {"Arho": 0.2, "Amu": -0.3, "AG": 0.4, "Ek": 1.0}


def _make_sweep_jsonl(path: Path, definition) -> None:
    run_parameter_sweep(definition, checkpoint_path=path)


class SweepSummaryTests(unittest.TestCase):
    def test_cartesian_grid_summary_and_npz_export(self):
        definition = cartesian_grid_definition(
            fixed={"Arho": 0.2, "Amu": -0.3},
            sweep={"AG": [0.4, 0.0], "Lambda": [1.7, 81.0]},
        )
        with tempfile.TemporaryDirectory() as tmp:
            jsonl = Path(tmp) / "sweep.jsonl"
            output = Path(tmp) / "summary"
            _make_sweep_jsonl(jsonl, definition)

            summary = read_sweep_summary(jsonl)
            paths = write_summary_outputs(summary, output)

            self.assertEqual(len(summary.point_rows), 4)
            self.assertEqual(summary.coordinate_system, "Lambda")
            self.assertEqual(summary.parameter_ordering, ("Arho", "Amu", "AG", "Lambda"))
            self.assertEqual(
                summary.aggregate_classification_counts["physically_admissible_roots"],
                3,
            )
            self.assertIn("grid_data_npz", paths)

            grid = np.load(paths["grid_data_npz"], allow_pickle=False)
            self.assertEqual(grid["classification_code"].shape, (2, 2))
            self.assertTrue(np.array_equal(grid["coord_AG"], np.array([0.4, 0.0])))
            self.assertTrue(np.array_equal(grid["coord_Lambda"], np.array([1.7, 81.0])))

            with Path(paths["point_summary_csv"]).open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 4)

    def test_explicit_list_summary_has_flat_tables_only(self):
        definition = explicit_sweep_definition([GENERIC, ONLY_SPURIOUS])
        with tempfile.TemporaryDirectory() as tmp:
            jsonl = Path(tmp) / "explicit.jsonl"
            output = Path(tmp) / "summary"
            _make_sweep_jsonl(jsonl, definition)

            summary = read_sweep_summary(jsonl)
            paths = write_summary_outputs(summary, output)

            self.assertEqual(summary.sweep_definition["kind"], "explicit")
            self.assertNotIn("grid_data_npz", paths)
            self.assertTrue(Path(paths["point_summary_csv"]).exists())
            self.assertTrue(Path(paths["admissible_roots_csv"]).exists())
            self.assertTrue(Path(paths["metadata_json"]).exists())

    def test_resumed_sweep_records_are_handled_consistently(self):
        definition = explicit_sweep_definition([GENERIC, ONLY_SPURIOUS])
        with tempfile.TemporaryDirectory() as tmp:
            jsonl = Path(tmp) / "resume.jsonl"
            run_parameter_sweep(definition, checkpoint_path=jsonl)
            with patch("src.parameter_sweep.analyze_parameter_point") as analyzer:
                run_parameter_sweep(definition, checkpoint_path=jsonl, resume=True)

            analyzer.assert_not_called()
            summary = read_sweep_summary(jsonl)
            self.assertEqual(len(summary.metadata_records), 2)
            self.assertEqual(len(summary.point_rows), 2)
            self.assertEqual(
                summary.validation_report["root_counts_match_point_rows"],
                True,
            )

    def test_complex_root_reconstruction_and_root_count_consistency(self):
        definition = explicit_sweep_definition([GENERIC])
        with tempfile.TemporaryDirectory() as tmp:
            jsonl = Path(tmp) / "roots.jsonl"
            _make_sweep_jsonl(jsonl, definition)
            summary = read_sweep_summary(jsonl)

            point = summary.point_rows[0]
            self.assertEqual(
                len(summary.admissible_root_rows),
                point["physically_admissible_root_count"],
            )
            first_root = summary.admissible_root_rows[0]
            self.assertIsInstance(first_root["s_real"], float)
            self.assertIsInstance(first_root["s_imag"], float)
            self.assertIsInstance(first_root["q_plus_real"], float)
            self.assertLess(first_root["absolute_residual"], 1e-8)

    def test_classification_code_consistency_and_Ek_preservation(self):
        self.assertEqual(
            CLASSIFICATION_CODES,
            {classification: index for index, classification in enumerate(CLASSIFICATIONS)},
        )
        definition = explicit_sweep_definition([OUTSIDE_EK])
        with tempfile.TemporaryDirectory() as tmp:
            jsonl = Path(tmp) / "ek.jsonl"
            _make_sweep_jsonl(jsonl, definition)
            summary = read_sweep_summary(jsonl)

            row = summary.point_rows[0]
            self.assertEqual(row["Ek"], 1.0)
            self.assertEqual(row["Lambda"], 0.0)
            self.assertEqual(row["final_classification"], "outside_theorem_scope")

    def test_malformed_and_incomplete_inputs_raise(self):
        with tempfile.TemporaryDirectory() as tmp:
            malformed = Path(tmp) / "bad.jsonl"
            malformed.write_text("{not-json}\n", encoding="utf-8")
            with self.assertRaises(SweepSummaryError):
                read_sweep_summary(malformed)

            no_metadata = Path(tmp) / "no_metadata.jsonl"
            no_metadata.write_text(
                json.dumps({"record_type": "point", "index": 0}) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(SweepSummaryError):
                read_sweep_summary(no_metadata)

    def test_duplicate_point_detection(self):
        definition = explicit_sweep_definition([GENERIC])
        with tempfile.TemporaryDirectory() as tmp:
            jsonl = Path(tmp) / "duplicate.jsonl"
            _make_sweep_jsonl(jsonl, definition)
            records = jsonl.read_text(encoding="utf-8").splitlines()
            jsonl.write_text("\n".join(records + [records[-1]]) + "\n", encoding="utf-8")

            with self.assertRaises(SweepSummaryError):
                read_sweep_summary(jsonl)

    def test_incomplete_point_set_detection(self):
        definition = explicit_sweep_definition([GENERIC, ONLY_SPURIOUS])
        with tempfile.TemporaryDirectory() as tmp:
            jsonl = Path(tmp) / "incomplete.jsonl"
            _make_sweep_jsonl(jsonl, definition)
            records = jsonl.read_text(encoding="utf-8").splitlines()
            jsonl.write_text("\n".join(records[:-1]) + "\n", encoding="utf-8")

            with self.assertRaises(SweepSummaryError):
                read_sweep_summary(jsonl)


if __name__ == "__main__":
    unittest.main()
