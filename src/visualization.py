"""Visualization helpers for Sprint 7 sweep-summary outputs.

The functions in this module consume only compact Sprint 7 summary artifacts:

- ``point_summary.csv``
- ``admissible_roots.csv``
- ``grid_data.npz``
- ``metadata.json``

They do not call the root solver, polynomial evaluator, dispersion evaluator,
or admissibility filters.  This layer is presentation-only.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping

_MPL_CACHE = Path(tempfile.gettempdir()) / "rmi_dispersion_matplotlib"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np


CLASSIFICATION_LABELS: dict[int, str] = {
    0: "outside theorem scope",
    1: "identically zero polynomial",
    2: "constant nonzero polynomial",
    3: "only domain-invalid candidates",
    4: "only spurious candidates",
    5: "genuine roots, nonadmissible",
    6: "marginal spatial decay",
    7: "physically admissible roots",
    8: "numerically unresolved",
}
"""Human-readable labels for categorical classification codes."""


CLASSIFICATION_COLORS: dict[int, str] = {
    0: "#d9d9d9",
    1: "#6a51a3",
    2: "#9e9ac8",
    3: "#fdae6b",
    4: "#f16913",
    5: "#9ecae1",
    6: "#fee391",
    7: "#31a354",
    8: "#de2d26",
}
"""Central categorical color policy.

Colors are categorical identifiers only.  They do not imply an ordinal ranking
of classifications.
"""


COUNT_MAP_FIELDS: dict[str, str] = {
    "mathematically_genuine_root_count": "Mathematically genuine root count",
    "physically_admissible_root_count": "Physically admissible root count",
    "marginal_root_count": "Marginal spatial-decay root count",
    "effective_polynomial_degree": "Effective polynomial degree",
}


class VisualizationError(ValueError):
    """Raised when summary artifacts are missing or inconsistent."""


@dataclass(frozen=True)
class SummaryPlotData:
    """Loaded Sprint 7 summary artifacts."""

    summary_dir: str
    point_rows: tuple[dict[str, Any], ...]
    admissible_root_rows: tuple[dict[str, Any], ...]
    metadata: dict[str, Any]
    grid_data: dict[str, np.ndarray] | None


def load_summary_data(summary_dir: str | Path) -> SummaryPlotData:
    """Load and validate a Sprint 7 summary directory."""

    directory = Path(summary_dir)
    point_path = directory / "point_summary.csv"
    root_path = directory / "admissible_roots.csv"
    metadata_path = directory / "metadata.json"
    for path in (point_path, root_path, metadata_path):
        if not path.exists():
            raise VisualizationError(f"Missing summary file: {path}")

    point_rows = tuple(_read_csv_rows(point_path))
    root_rows = tuple(_read_csv_rows(root_path))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    grid_path = directory / "grid_data.npz"
    grid_data = _load_npz(grid_path) if grid_path.exists() else None

    _validate_root_counts(point_rows, root_rows)
    if grid_data is not None:
        _validate_grid_data(point_rows, metadata, grid_data)

    return SummaryPlotData(
        summary_dir=str(directory),
        point_rows=point_rows,
        admissible_root_rows=root_rows,
        metadata=metadata,
        grid_data=grid_data,
    )


def _read_csv_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _load_npz(path: Path) -> dict[str, np.ndarray]:
    try:
        with np.load(path, allow_pickle=False) as data:
            return {key: data[key] for key in data.files}
    except Exception as exc:  # noqa: BLE001 - give a clean user-facing error
        raise VisualizationError(f"Could not read grid_data.npz: {exc}") from exc


def _as_int(value: Any, default: int = 0) -> int:
    if value in {None, ""}:
        return default
    return int(float(value))


def _as_float(value: Any) -> float:
    if value in {None, ""}:
        return float("nan")
    return float(value)


def _validate_root_counts(
    point_rows: tuple[dict[str, Any], ...],
    root_rows: tuple[dict[str, Any], ...],
) -> None:
    counts: dict[int, int] = {}
    for row in root_rows:
        point_id = _as_int(row["point_id"])
        counts[point_id] = counts.get(point_id, 0) + 1
    for row in point_rows:
        point_id = _as_int(row["point_id"])
        expected = _as_int(row["physically_admissible_root_count"])
        actual = counts.get(point_id, 0)
        if expected != actual:
            raise VisualizationError(
                f"Root-count mismatch for point {point_id}: "
                f"summary says {expected}, root table has {actual}"
            )


def _validate_grid_data(
    point_rows: tuple[dict[str, Any], ...],
    metadata: Mapping[str, Any],
    grid_data: Mapping[str, np.ndarray],
) -> None:
    required = {"grid_parameter_names", "classification_code"}
    missing = sorted(required - set(grid_data))
    if missing:
        raise VisualizationError("grid_data.npz missing: " + ", ".join(missing))
    names = _grid_parameter_names(grid_data)
    if len(names) < 1:
        raise VisualizationError("grid_data.npz contains no grid axes")
    for name in names:
        if f"coord_{name}" not in grid_data:
            raise VisualizationError(f"grid_data.npz missing coord_{name}")
    expected_points = len(point_rows)
    actual_points = int(np.asarray(grid_data["classification_code"]).size)
    if expected_points != actual_points:
        raise VisualizationError(
            f"Grid has {actual_points} point(s), point_summary has {expected_points}"
        )
    codes = np.asarray(grid_data["classification_code"])
    if np.isnan(codes.astype(float)).any() or np.any(codes < 0):
        raise VisualizationError("classification_code grid contains missing values")
    metadata_codes = metadata.get("classification_codes") or {}
    valid_codes = {int(value) for value in metadata_codes.values()} or set(
        CLASSIFICATION_LABELS
    )
    unknown = sorted(set(int(code) for code in np.ravel(codes)) - valid_codes)
    if unknown:
        raise VisualizationError(
            "classification_code grid contains unknown code(s): "
            + ", ".join(str(code) for code in unknown)
        )


def _grid_parameter_names(grid_data: Mapping[str, np.ndarray]) -> tuple[str, ...]:
    return tuple(str(name) for name in np.asarray(grid_data["grid_parameter_names"]))


def _require_two_dimensional_grid(data: SummaryPlotData) -> tuple[str, str, np.ndarray]:
    if data.grid_data is None:
        raise VisualizationError("Grid plot requested for an irregular/non-grid summary")
    names = _grid_parameter_names(data.grid_data)
    if len(names) != 2:
        raise VisualizationError(
            f"Grid plot requires exactly two swept parameters; found {len(names)}"
        )
    codes = np.asarray(data.grid_data["classification_code"])
    if codes.ndim != 2:
        raise VisualizationError(f"Grid arrays must be 2D; found shape {codes.shape}")
    return names[0], names[1], codes


def _coordinate_edges(centers: np.ndarray) -> np.ndarray:
    values = np.asarray(centers, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise VisualizationError("Grid coordinate arrays must be nonempty 1D arrays")
    if values.size == 1:
        return np.asarray([values[0] - 0.5, values[0] + 0.5])
    mids = 0.5 * (values[:-1] + values[1:])
    first = values[0] - (mids[0] - values[0])
    last = values[-1] + (values[-1] - mids[-1])
    return np.concatenate([[first], mids, [last]])


def _axis_label(name: str, overrides: Mapping[str, str] | None = None) -> str:
    if overrides and name in overrides:
        return overrides[name]
    labels = {
        "Arho": r"$A_\rho$",
        "Amu": r"$A_\mu$",
        "AG": r"$A_G$",
        "Lambda": r"$\Lambda$",
        "Ek": r"$E_k$",
    }
    return labels.get(name, name)


def _new_figure(figsize: tuple[float, float]) -> tuple[Any, Any]:
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


def plot_classification_map(
    summary_dir: str | Path,
    output_path: str | Path,
    *,
    title: str | None = None,
    figsize: tuple[float, float] = (6.0, 4.8),
    dpi: int = 150,
    axis_labels: Mapping[str, str] | None = None,
) -> str:
    """Plot a categorical final-classification phase diagram."""

    data = load_summary_data(summary_dir)
    y_name, x_name, codes = _require_two_dimensional_grid(data)
    grid = data.grid_data
    assert grid is not None

    y = np.asarray(grid[f"coord_{y_name}"], dtype=float)
    x = np.asarray(grid[f"coord_{x_name}"], dtype=float)
    z = np.asarray(codes, dtype=int)

    fig, ax = _new_figure(figsize)
    ordered_codes = sorted(CLASSIFICATION_LABELS)
    cmap = ListedColormap([CLASSIFICATION_COLORS[code] for code in ordered_codes])
    norm = BoundaryNorm(
        np.arange(min(ordered_codes) - 0.5, max(ordered_codes) + 1.5),
        cmap.N,
    )
    mesh = ax.pcolormesh(
        _coordinate_edges(x),
        _coordinate_edges(y),
        z,
        cmap=cmap,
        norm=norm,
        shading="flat",
    )
    colorbar = fig.colorbar(mesh, ax=ax, ticks=ordered_codes)
    colorbar.ax.set_yticklabels([CLASSIFICATION_LABELS[code] for code in ordered_codes])
    colorbar.set_label("Final classification")
    ax.set_xlabel(_axis_label(x_name, axis_labels))
    ax.set_ylabel(_axis_label(y_name, axis_labels))
    ax.set_title(title or "Parameter classification")
    fig.tight_layout()
    return _save_and_close(fig, output_path, dpi)


def plot_grid_scalar_map(
    summary_dir: str | Path,
    field: str,
    output_path: str | Path,
    *,
    title: str | None = None,
    colorbar_label: str | None = None,
    figsize: tuple[float, float] = (6.0, 4.8),
    dpi: int = 150,
    axis_labels: Mapping[str, str] | None = None,
) -> str:
    """Plot one scalar grid field from ``grid_data.npz``."""

    data = load_summary_data(summary_dir)
    y_name, x_name, _ = _require_two_dimensional_grid(data)
    grid = data.grid_data
    assert grid is not None
    if field not in grid:
        raise VisualizationError(f"grid_data.npz does not contain field: {field}")

    y = np.asarray(grid[f"coord_{y_name}"], dtype=float)
    x = np.asarray(grid[f"coord_{x_name}"], dtype=float)
    z = np.asarray(grid[field], dtype=float)
    z = np.where(z < 0, np.nan, z)

    fig, ax = _new_figure(figsize)
    mesh = ax.pcolormesh(
        _coordinate_edges(x),
        _coordinate_edges(y),
        z,
        shading="flat",
        cmap="viridis",
    )
    colorbar = fig.colorbar(mesh, ax=ax)
    label = colorbar_label or COUNT_MAP_FIELDS.get(field, field)
    colorbar.set_label(label)
    ax.set_xlabel(_axis_label(x_name, axis_labels))
    ax.set_ylabel(_axis_label(y_name, axis_labels))
    ax.set_title(title or label)
    fig.tight_layout()
    return _save_and_close(fig, output_path, dpi)


def plot_standard_grid_maps(
    summary_dir: str | Path,
    output_dir: str | Path,
    *,
    image_format: str = "png",
    dpi: int = 150,
    figsize: tuple[float, float] = (6.0, 4.8),
) -> dict[str, str]:
    """Generate standard count/degree grid maps."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    for field, label in COUNT_MAP_FIELDS.items():
        name = field.replace("_root_count", "").replace("_", "-")
        path = output / f"{name}.{image_format}"
        paths[field] = plot_grid_scalar_map(
            summary_dir,
            field,
            path,
            title=label,
            colorbar_label=label,
            dpi=dpi,
            figsize=figsize,
        )
    return paths


def plot_admissible_roots_complex_plane(
    summary_dir: str | Path,
    output_path: str | Path,
    *,
    color_by: str | None = None,
    filters: Mapping[str, tuple[float, float]] | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (6.0, 4.8),
    dpi: int = 150,
) -> str:
    """Plot physically admissible roots in the complex ``s`` plane."""

    data = load_summary_data(summary_dir)
    point_by_id = {_as_int(row["point_id"]): row for row in data.point_rows}
    rows = [
        row
        for row in data.admissible_root_rows
        if _point_passes_filters(point_by_id[_as_int(row["point_id"])], filters)
    ]

    fig, ax = _new_figure(figsize)
    ax.axhline(0.0, color="0.6", linewidth=0.8)
    ax.axvline(0.0, color="0.6", linewidth=0.8)
    if rows:
        s_real = np.asarray([_as_float(row["s_real"]) for row in rows])
        s_imag = np.asarray([_as_float(row["s_imag"]) for row in rows])
        if color_by is None:
            ax.scatter(s_real, s_imag, s=28, color="#3182bd", edgecolor="none")
        else:
            color_values = np.asarray(
                [_as_float(point_by_id[_as_int(row["point_id"])].get(color_by)) for row in rows]
            )
            scatter = ax.scatter(
                s_real,
                s_imag,
                c=color_values,
                s=28,
                cmap="viridis",
                edgecolor="none",
            )
            colorbar = fig.colorbar(scatter, ax=ax)
            colorbar.set_label(_axis_label(color_by))
    ax.set_xlabel(r"$\operatorname{Re}(s)$")
    ax.set_ylabel(r"$\operatorname{Im}(s)$")
    ax.set_title(title or "Physically admissible roots")
    fig.tight_layout()
    return _save_and_close(fig, output_path, dpi)


def _point_passes_filters(
    point_row: Mapping[str, Any],
    filters: Mapping[str, tuple[float, float]] | None,
) -> bool:
    if not filters:
        return True
    for name, (lower, upper) in filters.items():
        value = _as_float(point_row.get(name))
        if not (lower <= value <= upper):
            return False
    return True


def plot_root_component_maps(
    summary_dir: str | Path,
    output_dir: str | Path,
    *,
    image_format: str = "png",
    dpi: int = 150,
    figsize: tuple[float, float] = (6.0, 4.8),
) -> dict[str, str]:
    """Conditionally plot root-component maps when no branch choice is needed.

    Sprint 8 does not perform branch tracking.  This helper only supports the
    trivial unambiguous case where every parameter point has exactly one
    physically admissible root.  Otherwise it writes no files and returns an
    explanatory ``{"skipped": ...}`` entry.
    """

    data = load_summary_data(summary_dir)
    if data.grid_data is None:
        raise VisualizationError("Root-component maps require a regular grid summary")
    if any(_as_int(row["physically_admissible_root_count"]) != 1 for row in data.point_rows):
        return {
            "skipped": (
                "Root-component maps require exactly one admissible root at every "
                "parameter point; branch tracking belongs to a later sprint."
            )
        }

    y_name, x_name, _ = _require_two_dimensional_grid(data)
    grid = data.grid_data
    assert grid is not None
    shape = np.asarray(grid["classification_code"]).shape
    re_grid = np.full(shape, np.nan)
    im_grid = np.full(shape, np.nan)
    for row in data.admissible_root_rows:
        index = np.unravel_index(_as_int(row["point_id"]), shape, order="C")
        re_grid[index] = _as_float(row["s_real"])
        im_grid[index] = _as_float(row["s_imag"])

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    for field_name, z, label in (
        ("root_real_component", re_grid, r"$\operatorname{Re}(s)$"),
        ("root_imag_component", im_grid, r"$\operatorname{Im}(s)$"),
    ):
        fig, ax = _new_figure(figsize)
        mesh = ax.pcolormesh(
            _coordinate_edges(np.asarray(grid[f"coord_{x_name}"], dtype=float)),
            _coordinate_edges(np.asarray(grid[f"coord_{y_name}"], dtype=float)),
            z,
            shading="flat",
            cmap="viridis",
        )
        fig.colorbar(mesh, ax=ax).set_label(label)
        ax.set_xlabel(_axis_label(x_name))
        ax.set_ylabel(_axis_label(y_name))
        ax.set_title(label)
        fig.tight_layout()
        paths[field_name] = _save_and_close(fig, output / f"{field_name}.{image_format}", dpi)
    return paths


def _save_and_close(fig: Any, output_path: str | Path, dpi: int) -> str:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    return str(path)
