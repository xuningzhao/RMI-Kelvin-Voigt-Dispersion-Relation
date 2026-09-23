#!/usr/bin/env python
"""Generate Sprint 8 figures from Sprint 7 sweep-summary outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.visualization import (  # noqa: E402
    COUNT_MAP_FIELDS,
    VisualizationError,
    plot_admissible_roots_complex_plane,
    plot_classification_map,
    plot_grid_scalar_map,
    plot_root_component_maps,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_dir", type=Path, help="Sprint 7 summary directory")
    parser.add_argument("--output-dir", type=Path, help="Figure output directory")
    parser.add_argument("--format", default="png", choices=("png", "pdf", "svg"))
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--figsize", nargs=2, type=float, default=(6.0, 4.8), metavar=("W", "H"))
    parser.add_argument("--title", help="Optional title for single selected plots")
    parser.add_argument("--classification-map", action="store_true")
    parser.add_argument("--root-count-maps", action="store_true")
    parser.add_argument("--effective-degree-map", action="store_true")
    parser.add_argument("--complex-plane", action="store_true")
    parser.add_argument("--root-component-maps", action="store_true")
    parser.add_argument(
        "--color-by",
        default=None,
        help="Point-summary column used to color complex-plane roots, e.g. Ek or Lambda",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Generate classification, root-count, degree, and complex-plane figures",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or args.summary_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figsize = (float(args.figsize[0]), float(args.figsize[1]))

    make_classification = args.all or args.classification_map
    make_counts = args.all or args.root_count_maps
    make_degree = args.all or args.effective_degree_map
    make_complex = args.all or args.complex_plane
    make_components = args.root_component_maps
    if not any((make_classification, make_counts, make_degree, make_complex, make_components)):
        make_classification = make_counts = make_degree = make_complex = True

    created: list[str] = []
    try:
        if make_classification:
            created.append(
                plot_classification_map(
                    args.summary_dir,
                    output_dir / f"classification_map.{args.format}",
                    title=args.title,
                    dpi=args.dpi,
                    figsize=figsize,
                )
            )
        if make_counts:
            for field in (
                "mathematically_genuine_root_count",
                "physically_admissible_root_count",
                "marginal_root_count",
            ):
                created.append(
                    plot_grid_scalar_map(
                        args.summary_dir,
                        field,
                        output_dir / f"{field}.{args.format}",
                        title=COUNT_MAP_FIELDS[field],
                        colorbar_label=COUNT_MAP_FIELDS[field],
                        dpi=args.dpi,
                        figsize=figsize,
                    )
                )
        if make_degree:
            field = "effective_polynomial_degree"
            created.append(
                plot_grid_scalar_map(
                    args.summary_dir,
                    field,
                    output_dir / f"{field}.{args.format}",
                    title=COUNT_MAP_FIELDS[field],
                    colorbar_label=COUNT_MAP_FIELDS[field],
                    dpi=args.dpi,
                    figsize=figsize,
                )
            )
        if make_complex:
            created.append(
                plot_admissible_roots_complex_plane(
                    args.summary_dir,
                    output_dir / f"admissible_roots_complex_plane.{args.format}",
                    color_by=args.color_by,
                    title=args.title,
                    dpi=args.dpi,
                    figsize=figsize,
                )
            )
        if make_components:
            component_paths = plot_root_component_maps(
                args.summary_dir,
                output_dir,
                image_format=args.format,
                dpi=args.dpi,
                figsize=figsize,
            )
            skipped = component_paths.pop("skipped", None)
            if skipped:
                print(f"Root-component maps skipped: {skipped}")
            created.extend(component_paths.values())
    except VisualizationError as exc:
        raise SystemExit(f"Plotting failed: {exc}") from exc

    print("Created figures")
    for path in created:
        print(f"  {path}")


if __name__ == "__main__":
    main()
