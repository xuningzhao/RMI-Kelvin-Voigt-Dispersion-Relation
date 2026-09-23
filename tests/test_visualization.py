"""Tests for Sprint 8 visualization layer."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from src.sweep_summary import CLASSIFICATION_CODES, POINT_SUMMARY_FIELDS, ROOT_SUMMARY_FIELDS
from src.visualization import (
    CLASSIFICATION_LABELS,
    VisualizationError,
    _require_two_dimensional_grid,
    load_summary_data,
    plot_admissible_roots_complex_plane,
    plot_classification_map,
    plot_grid_scalar_map,
)


def _write_csv(path: Path, fieldnames, rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _synthetic_summary_dir(base: Path, *, grid: bool = True) -> Path:
    directory = base / ("grid_summary" if grid else "explicit_summary")
    directory.mkdir()
    point_rows = [
        {
            "point_id": 0,
            "Arho": 0.2,
            "Amu": -0.3,
            "AG": 0.0,
            "Lambda": 2.25,
            "Ek": -0.2,
            "theorem_scope_applies": True,
            "theorem_scope_failed_conditions": "",
            "final_classification": "physically_admissible_roots",
            "classification_code": CLASSIFICATION_CODES["physically_admissible_roots"],
            "effective_polynomial_degree": 14,
            "reduced_degree": False,
            "polynomial_candidate_count": 14,
            "domain_valid_candidate_count": 14,
            "mathematically_genuine_root_count": 1,
            "physically_admissible_root_count": 1,
            "marginal_root_count": 0,
            "physically_nonadmissible_root_count": 0,
            "numerically_unresolved": False,
            "runtime_seconds": 0.01,
            "warning_count": 0,
        },
        {
            "point_id": 1,
            "Arho": 0.2,
            "Amu": -0.3,
            "AG": 0.0,
            "Lambda": 0.0,
            "Ek": 1.0,
            "theorem_scope_applies": False,
            "theorem_scope_failed_conditions": "Lambda_positive",
            "final_classification": "outside_theorem_scope",
            "classification_code": CLASSIFICATION_CODES["outside_theorem_scope"],
            "effective_polynomial_degree": "",
            "reduced_degree": False,
            "polynomial_candidate_count": 0,
            "domain_valid_candidate_count": 0,
            "mathematically_genuine_root_count": 0,
            "physically_admissible_root_count": 0,
            "marginal_root_count": 0,
            "physically_nonadmissible_root_count": 0,
            "numerically_unresolved": False,
            "runtime_seconds": 0.01,
            "warning_count": 0,
        },
        {
            "point_id": 2,
            "Arho": 0.2,
            "Amu": -0.3,
            "AG": 0.4,
            "Lambda": 2.25,
            "Ek": -0.2,
            "theorem_scope_applies": True,
            "theorem_scope_failed_conditions": "",
            "final_classification": "only_spurious_candidates",
            "classification_code": CLASSIFICATION_CODES["only_spurious_candidates"],
            "effective_polynomial_degree": 14,
            "reduced_degree": False,
            "polynomial_candidate_count": 14,
            "domain_valid_candidate_count": 14,
            "mathematically_genuine_root_count": 0,
            "physically_admissible_root_count": 0,
            "marginal_root_count": 0,
            "physically_nonadmissible_root_count": 0,
            "numerically_unresolved": False,
            "runtime_seconds": 0.01,
            "warning_count": 0,
        },
        {
            "point_id": 3,
            "Arho": 0.2,
            "Amu": -0.3,
            "AG": 0.4,
            "Lambda": 0.0,
            "Ek": 1.0,
            "theorem_scope_applies": False,
            "theorem_scope_failed_conditions": "Lambda_positive",
            "final_classification": "numerically_unresolved",
            "classification_code": CLASSIFICATION_CODES["numerically_unresolved"],
            "effective_polynomial_degree": "",
            "reduced_degree": False,
            "polynomial_candidate_count": 0,
            "domain_valid_candidate_count": 0,
            "mathematically_genuine_root_count": 0,
            "physically_admissible_root_count": 0,
            "marginal_root_count": 0,
            "physically_nonadmissible_root_count": 0,
            "numerically_unresolved": True,
            "runtime_seconds": 0.01,
            "warning_count": 0,
        },
    ]
    root_rows = [
        {
            "point_id": 0,
            "root_index": 0,
            "s_real": -1.0,
            "s_imag": 2.0,
            "temporal_status": "decaying",
            "q_plus_real": 1.2,
            "q_plus_imag": -0.4,
            "q_minus_real": 0.7,
            "q_minus_imag": 0.5,
            "absolute_residual": 1e-13,
            "relative_residual": 1e-14,
            "admissibility_status": "physically_admissible",
        }
    ]
    _write_csv(directory / "point_summary.csv", POINT_SUMMARY_FIELDS, point_rows)
    _write_csv(directory / "admissible_roots.csv", ROOT_SUMMARY_FIELDS, root_rows)

    metadata = {
        "source_jsonl_path": "/tmp/synthetic.jsonl",
        "workflow_version": "p14-root-workflow-sprint6",
        "sweep_definition": {
            "kind": "cartesian_grid" if grid else "explicit",
            "fixed_values": {"Arho": 0.2, "Amu": -0.3},
            "swept_values": {"AG": [0.0, 0.4], "Ek": [-0.2, 1.0]} if grid else {},
            "coordinate_system": "Ek",
            "point_count": 4,
            "points": [],
        },
        "analysis_options": {},
        "creation_time": 0.0,
        "parameter_ordering": ["Arho", "Amu", "AG", "Ek"],
        "coordinate_system": "Ek",
        "classification_codes": CLASSIFICATION_CODES,
    }
    (directory / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    if grid:
        np.savez(
            directory / "grid_data.npz",
            grid_parameter_names=np.asarray(["AG", "Ek"]),
            coord_AG=np.asarray([0.0, 0.4]),
            coord_Ek=np.asarray([-0.2, 1.0]),
            classification_code=np.asarray([[7, 0], [4, 8]]),
            mathematically_genuine_root_count=np.asarray([[1, 0], [0, 0]]),
            physically_admissible_root_count=np.asarray([[1, 0], [0, 0]]),
            marginal_root_count=np.asarray([[0, 0], [0, 0]]),
            physically_nonadmissible_root_count=np.asarray([[0, 0], [0, 0]]),
            effective_polynomial_degree=np.asarray([[14, -1], [14, -1]]),
            polynomial_candidate_count=np.asarray([[14, 0], [14, 0]]),
            domain_valid_candidate_count=np.asarray([[14, 0], [14, 0]]),
            runtime_seconds=np.asarray([[0.01, 0.01], [0.01, 0.01]]),
        )
    return directory


class VisualizationTests(unittest.TestCase):
    def test_classification_map_generation_and_colorbar_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_dir = _synthetic_summary_dir(Path(tmp))
            output = Path(tmp) / "classification.png"
            with patch("matplotlib.axes.Axes.set_yticklabels") as set_labels:
                path = plot_classification_map(summary_dir, output, dpi=80)

            self.assertTrue(Path(path).exists())
            label_args = [call.args[0] for call in set_labels.call_args_list if call.args]
            self.assertTrue(any(list(CLASSIFICATION_LABELS.values()) == list(args) for args in label_args))

    def test_axis_orientation_uses_stored_grid_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_dir = _synthetic_summary_dir(Path(tmp))
            data = load_summary_data(summary_dir)
            y_name, x_name, grid = _require_two_dimensional_grid(data)

            self.assertEqual((y_name, x_name), ("AG", "Ek"))
            self.assertEqual(grid.shape, (2, 2))

    def test_root_count_and_degree_maps(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_dir = _synthetic_summary_dir(Path(tmp))
            for field in (
                "mathematically_genuine_root_count",
                "physically_admissible_root_count",
                "marginal_root_count",
                "effective_polynomial_degree",
            ):
                output = Path(tmp) / f"{field}.png"
                self.assertTrue(
                    Path(plot_grid_scalar_map(summary_dir, field, output, dpi=80)).exists()
                )

    def test_complex_plane_root_plot(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_dir = _synthetic_summary_dir(Path(tmp))
            output = Path(tmp) / "roots.png"
            path = plot_admissible_roots_complex_plane(
                summary_dir,
                output,
                color_by="Ek",
                dpi=80,
            )
            self.assertTrue(Path(path).exists())

    def test_irregular_sweep_rejected_for_grid_plots(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_dir = _synthetic_summary_dir(Path(tmp), grid=False)
            with self.assertRaises(VisualizationError):
                plot_classification_map(summary_dir, Path(tmp) / "bad.png")

    def test_missing_or_malformed_summary_files_raise(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(VisualizationError):
                load_summary_data(Path(tmp))

            summary_dir = _synthetic_summary_dir(Path(tmp))
            (summary_dir / "grid_data.npz").write_text("not an npz", encoding="utf-8")
            with self.assertRaises(VisualizationError):
                load_summary_data(summary_dir)

    def test_no_solver_or_root_verification_functions_are_called(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_dir = _synthetic_summary_dir(Path(tmp))
            output = Path(tmp) / "classification.png"
            with patch("src.workflow.analyze_parameter_point") as analyze:
                with patch("src.root_verification.verify_candidate_roots") as verify:
                    plot_classification_map(summary_dir, output, dpi=80)

            analyze.assert_not_called()
            verify.assert_not_called()


if __name__ == "__main__":
    unittest.main()

