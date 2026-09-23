"""Shared physical layout and Matplotlib style for PRE manuscript figures.

Consistency is defined at the subplot level.  Every main axes created here has
the same physical dimensions; colorbars occupy additional, separate axes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


SUBPLOT_WIDTH_IN = 5.00
SUBPLOT_HEIGHT_IN = 5.00
LEFT_MARGIN_IN = 1.25
RIGHT_MARGIN_IN = 0.22
BOTTOM_MARGIN_IN = 0.78
TOP_MARGIN_IN = 0.65
COLUMN_GAP_IN = 1.80
ROW_GAP_IN = 1.00
COLORBAR_PAD_IN = 0.27
COLORBAR_HEIGHT_IN = 0.18

AXIS_LABEL_FONTSIZE = 21
TICK_LABEL_FONTSIZE = 18.75
LEGEND_FONTSIZE = 18.75
TITLE_FONTSIZE = 21
PANEL_LABEL_FONTSIZE = 21
CRITICAL_LABEL_FONTSIZE = 14.25

BRANCH_LINEWIDTH = 1.5
CRITICAL_LINEWIDTH = 1.0
AXIS_LINEWIDTH = 0.8
MAJOR_GRID_LINEWIDTH = 0.55
MINOR_GRID_LINEWIDTH = 0.35
MAJOR_TICK_LENGTH = 4.0
MINOR_TICK_LENGTH = 2.3
MARKER_SIZE = 3.5
SHADE_ALPHA = 0.10
DEFAULT_DPI = 300


RC_PARAMS = {
    "font.family": "sans-serif",
    "font.size": TICK_LABEL_FONTSIZE,
    "axes.labelsize": AXIS_LABEL_FONTSIZE,
    "axes.titlesize": TITLE_FONTSIZE,
    "axes.linewidth": AXIS_LINEWIDTH,
    "xtick.labelsize": TICK_LABEL_FONTSIZE,
    "ytick.labelsize": TICK_LABEL_FONTSIZE,
    "xtick.major.size": MAJOR_TICK_LENGTH,
    "ytick.major.size": MAJOR_TICK_LENGTH,
    "xtick.minor.size": MINOR_TICK_LENGTH,
    "ytick.minor.size": MINOR_TICK_LENGTH,
    "xtick.major.width": AXIS_LINEWIDTH,
    "ytick.major.width": AXIS_LINEWIDTH,
    "xtick.minor.width": AXIS_LINEWIDTH,
    "ytick.minor.width": AXIS_LINEWIDTH,
    "legend.fontsize": LEGEND_FONTSIZE,
    "lines.linewidth": BRANCH_LINEWIDTH,
    "lines.markersize": MARKER_SIZE,
    "savefig.dpi": DEFAULT_DPI,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


@dataclass(frozen=True)
class PublicationLayout:
    figure: Figure
    axes: tuple[Axes, ...]
    colorbar_axes: tuple[Axes | None, ...]
    figsize: tuple[float, float]


def apply_publication_style() -> None:
    """Apply the project-wide PRE plotting style."""
    mpl.rcParams.update(RC_PARAMS)


def publication_figure(
    nrows: int,
    ncols: int,
    *,
    colorbars: Sequence[bool] | None = None,
    row_gap_in: float = ROW_GAP_IN,
) -> PublicationLayout:
    """Create fixed-size square subplot axes and optional external colorbars."""
    if nrows < 1 or ncols < 1:
        raise ValueError("nrows and ncols must be positive")
    panel_count = nrows * ncols
    flags = tuple(colorbars or (False,) * panel_count)
    if len(flags) != panel_count:
        raise ValueError("colorbars must have one flag per subplot")

    row_extras = []
    for row in range(nrows):
        has_colorbar = any(flags[row * ncols + column] for column in range(ncols))
        row_extras.append(COLORBAR_PAD_IN + COLORBAR_HEIGHT_IN if has_colorbar else 0.0)
    width = (
        LEFT_MARGIN_IN
        + ncols * SUBPLOT_WIDTH_IN
        + (ncols - 1) * COLUMN_GAP_IN
        + RIGHT_MARGIN_IN
    )
    height = (
        BOTTOM_MARGIN_IN
        + nrows * SUBPLOT_HEIGHT_IN
        + sum(row_extras)
        + (nrows - 1) * row_gap_in
        + TOP_MARGIN_IN
    )
    figure = plt.figure(figsize=(width, height))
    axes: list[Axes] = []
    colorbar_axes: list[Axes | None] = []
    column_lefts = []
    cursor = LEFT_MARGIN_IN
    for column in range(ncols):
        column_lefts.append(cursor)
        cursor += SUBPLOT_WIDTH_IN + COLUMN_GAP_IN

    row_bottoms = [0.0] * nrows
    cursor = BOTTOM_MARGIN_IN
    for row in range(nrows - 1, -1, -1):
        row_bottoms[row] = cursor
        cursor += SUBPLOT_HEIGHT_IN + row_extras[row] + row_gap_in
    for row in range(nrows):
        bottom = row_bottoms[row]
        for column in range(ncols):
            index = row * ncols + column
            left = column_lefts[column]
            axes.append(
                figure.add_axes(
                    [
                        left / width,
                        bottom / height,
                        SUBPLOT_WIDTH_IN / width,
                        SUBPLOT_HEIGHT_IN / height,
                    ]
                )
            )
            if flags[index]:
                colorbar_axes.append(
                    figure.add_axes(
                        [
                            left / width,
                            (bottom + SUBPLOT_HEIGHT_IN + COLORBAR_PAD_IN) / height,
                            SUBPLOT_WIDTH_IN / width,
                            COLORBAR_HEIGHT_IN / height,
                        ]
                    )
                )
            else:
                colorbar_axes.append(None)
    return PublicationLayout(figure, tuple(axes), tuple(colorbar_axes), (width, height))


def style_parameter_axes(ax: Axes, *, equal: bool = False, grid: bool = True) -> None:
    """Apply common limits-independent styling to a parameter/spectrum axes."""
    if grid:
        ax.minorticks_on()
        ax.grid(True, which="major", color="0.86", linewidth=MAJOR_GRID_LINEWIDTH)
        ax.grid(True, which="minor", color="0.92", linewidth=MINOR_GRID_LINEWIDTH)
    if equal:
        ax.set_aspect("equal", adjustable="box")


def panel_heading(
    ax: Axes,
    label: str,
    subtitle: str | None = None,
    *,
    above_colorbar: bool = False,
) -> None:
    """Place the panel label at the outside upper-left, followed by a subtitle."""
    y_position = (
        1.0
        + (COLORBAR_PAD_IN + COLORBAR_HEIGHT_IN + 0.40) / SUBPLOT_HEIGHT_IN
        if above_colorbar
        else 1.055
    )
    ax.text(
        -0.24,
        y_position,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=PANEL_LABEL_FONTSIZE,
        clip_on=False,
    )
    if subtitle:
        ax.text(
            0.5,
            y_position,
            subtitle,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=TITLE_FONTSIZE,
            clip_on=False,
        )


def contrast_ticks(ax: Axes, *, x: bool = False, y: bool = False) -> None:
    """Use manuscript-standard 0.5 major spacing for A_G and E_k axes."""
    from matplotlib.ticker import MultipleLocator

    if x:
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
    if y:
        ax.yaxis.set_major_locator(MultipleLocator(0.5))


def style_horizontal_colorbar(colorbar: object) -> None:
    """Place horizontal colorbar ticks above the bar at shared font scale."""
    colorbar.ax.xaxis.set_ticks_position("top")
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.ax.tick_params(axis="x", labelsize=TICK_LABEL_FONTSIZE, pad=2)


apply_publication_style()
