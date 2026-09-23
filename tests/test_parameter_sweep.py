"""Tests for batch parameter-sweep infrastructure."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from src.parameter_sweep import (
    cartesian_grid_definition,
    completed_indices_from_checkpoint,
    explicit_sweep_definition,
    lambda_from_ek,
    load_checkpoint_records,
    parameter_point_to_json,
    run_parameter_sweep,
)
from src.workflow import analyze_parameter_point


GENERIC = {"Arho": 0.2, "Amu": -0.3, "AG": 0.4, "Lambda": 1.7}
ONLY_SPURIOUS = {"Arho": 0.5, "Amu": 0.5, "AG": 0.0, "Lambda": 81.0}
OUTSIDE = {"Arho": 0.2, "Amu": -0.3, "AG": 0.4, "Lambda": 0.0}


class ParameterSweepTests(unittest.TestCase):
    def test_small_serial_sweep_and_mixed_classifications(self):
        definition = explicit_sweep_definition([GENERIC, ONLY_SPURIOUS, OUTSIDE])
        result = run_parameter_sweep(definition)

        self.assertEqual(len(result.points), 3)
        self.assertEqual(
            [point.classification for point in result.points],
            [
                "physically_admissible_roots",
                "only_spurious_candidates",
                "outside_theorem_scope",
            ],
        )
        self.assertEqual(result.aggregate_classification_counts["physically_admissible_roots"], 1)
        self.assertEqual(result.aggregate_classification_counts["only_spurious_candidates"], 1)
        self.assertEqual(result.aggregate_classification_counts["outside_theorem_scope"], 1)

    def test_deterministic_ordering_for_cartesian_grid(self):
        definition = cartesian_grid_definition(
            fixed={"Arho": 0.2, "Amu": -0.3},
            sweep={"AG": [0.4, 0.0], "Lambda": [1.7, 81.0]},
        )
        result = run_parameter_sweep(definition)

        self.assertEqual([point.index for point in result.points], [0, 1, 2, 3])
        self.assertEqual(
            [point.coordinates["AG"] for point in result.points],
            [0.4, 0.4, 0.0, 0.0],
        )
        self.assertEqual(
            [point.coordinates["Lambda"] for point in result.points],
            [1.7, 81.0, 1.7, 81.0],
        )

    def test_sweep_result_agrees_with_direct_analyze_call(self):
        definition = explicit_sweep_definition([GENERIC])
        sweep = run_parameter_sweep(definition)
        direct = analyze_parameter_point(**GENERIC)

        self.assertEqual(sweep.points[0].classification, direct.classification)
        self.assertEqual(
            sweep.points[0].result.counts.mathematically_genuine_roots,
            direct.counts.mathematically_genuine_roots,
        )

    def test_one_failed_point_does_not_terminate_batch(self):
        definition = explicit_sweep_definition([GENERIC, ONLY_SPURIOUS])

        def flaky(Arho, Amu, AG, Lambda, **kwargs):
            if Lambda == ONLY_SPURIOUS["Lambda"]:
                raise RuntimeError("intentional failure")
            return analyze_parameter_point(Arho, Amu, AG, Lambda, **kwargs)

        with patch("src.parameter_sweep.analyze_parameter_point", side_effect=flaky):
            result = run_parameter_sweep(definition)

        self.assertEqual(len(result.points), 2)
        self.assertEqual(result.points[0].classification, "physically_admissible_roots")
        self.assertEqual(result.points[1].classification, "numerically_unresolved")
        self.assertEqual(result.unresolved_points, (1,))
        self.assertEqual(result.points[1].error, "RuntimeError")

    def test_checkpoint_creation_and_resume_skip_completed_points(self):
        definition = explicit_sweep_definition([GENERIC, ONLY_SPURIOUS])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sweep.jsonl"
            first = run_parameter_sweep(definition, checkpoint_path=path)

            self.assertTrue(path.exists())
            self.assertEqual(completed_indices_from_checkpoint(path), {0, 1})

            with patch("src.parameter_sweep.analyze_parameter_point") as analyzer:
                resumed = run_parameter_sweep(
                    definition,
                    checkpoint_path=path,
                    resume=True,
                )

            analyzer.assert_not_called()
            self.assertEqual(len(resumed.points), 2)
            self.assertEqual([point.index for point in resumed.points], [0, 1])
            self.assertTrue(all(point.skipped_existing for point in resumed.points))
            self.assertTrue(all(point.checkpoint_record is not None for point in resumed.points))
            self.assertEqual(
                [point.classification for point in resumed.points],
                ["physically_admissible_roots", "only_spurious_candidates"],
            )
            records = load_checkpoint_records(path)
            self.assertTrue(any(record.get("record_type") == "metadata" for record in records))
            self.assertEqual(first.checkpoint_path, str(path))

    def test_serialization_preserves_complex_roots(self):
        result = analyze_parameter_point(**GENERIC)
        payload = parameter_point_to_json(result)

        self.assertTrue(payload["polynomial_candidates"])
        first = payload["polynomial_candidates"][0]
        self.assertIn("real", first)
        self.assertIn("imag", first)

        encoded = json.dumps(payload)
        decoded = json.loads(encoded)
        self.assertEqual(decoded["classification"], result.classification)

    def test_parallel_execution_optional(self):
        definition = explicit_sweep_definition([GENERIC, ONLY_SPURIOUS])
        result = run_parameter_sweep(definition, max_workers=2)

        self.assertEqual([point.index for point in result.points], [0, 1])
        self.assertEqual(len(result.points), 2)

    def test_Ek_to_Lambda_conversion_consistency(self):
        ek = -0.2
        expected = ((1.0 - ek) / (1.0 + ek)) ** 2
        self.assertAlmostEqual(lambda_from_ek(ek), expected)

        definition = explicit_sweep_definition(
            [{"Arho": 0.2, "Amu": -0.3, "AG": 0.4, "Ek": ek}]
        )
        point = definition.points[0]
        self.assertIn("Ek", point.coordinates)
        self.assertIn("Lambda", point.coordinates)
        self.assertAlmostEqual(point.parameters["Lambda"], expected)

    def test_Ek_endpoint_one_maps_to_outside_theorem_scope(self):
        definition = explicit_sweep_definition(
            [{"Arho": 0.2, "Amu": -0.3, "AG": 0.4, "Ek": 1.0}]
        )
        result = run_parameter_sweep(definition)

        self.assertEqual(result.points[0].parameters["Lambda"], 0.0)
        self.assertEqual(result.points[0].classification, "outside_theorem_scope")


if __name__ == "__main__":
    unittest.main()
