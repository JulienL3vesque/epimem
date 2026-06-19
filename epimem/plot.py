"""Optional plotting for epimem: the MEM surveillance chart.

Needs matplotlib, an optional dependency - install with `pip install epimem[plot]`. matplotlib is
imported lazily inside the function, so the core package stays numpy/scipy-only and a plain
`import epimem` never pulls it in.

    from epimem import mem_model
    from epimem.plot import mem_chart
    model = mem_model(season_matrix)
    mem_chart(season_matrix, model, season_labels=[...], ylabel="ILI share of 811 calls (%)")
"""
import numpy as np

from .model import MemThresholds

# The four threshold reference lines, in increasing severity: label, MemThresholds attribute, colour.
_THRESHOLD_LINES = (
    ("epidemic onset", "epidemic_onset", "#378ADD"),
    ("medium", "medium", "#639922"),
    ("high", "high", "#BA7517"),
    ("very high", "very_high", "#E24B4A"),
)

PAST_SEASON_COLOUR = "#B4B2A9"
CURRENT_SEASON_COLOUR = "#534AB7"


def mem_chart(seasons, model: MemThresholds, *, current_index: int | None = None,
              season_labels=None, week_numbers=None, ylabel: str = "value",
              title: str | None = None, ax=None):
    """Draw the MEM surveillance chart and return the matplotlib Axes.

    `seasons` is 2-D (rows = surveillance weeks, columns = seasons); each column is drawn as a
    line, past seasons in grey and the `current_index` column (default: the last) bold. `model`
    supplies the epidemic-onset and medium / high / very-high levels, drawn as labelled horizontal
    reference lines. Pass `ax` to draw onto an existing Axes; otherwise a new figure is created.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - only reachable without the [plot] extra
        raise ImportError("mem_chart needs matplotlib; install epimem[plot]") from exc

    matrix = np.asarray(seasons, dtype=float)
    if matrix.ndim == 1:
        matrix = matrix[:, None]
    week_count, season_count = matrix.shape
    weeks = np.arange(week_count) if week_numbers is None else np.asarray(week_numbers)
    if current_index is None:
        current_index = season_count - 1

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 5))

    for season in range(season_count):
        is_current = season == current_index
        ax.plot(
            weeks, matrix[:, season],
            color=CURRENT_SEASON_COLOUR if is_current else PAST_SEASON_COLOUR,
            linewidth=2.6 if is_current else 1.3,
            zorder=3 if is_current else 2,
            label=None if season_labels is None else season_labels[season],
        )

    for name, attribute, colour in _THRESHOLD_LINES:
        level = getattr(model, attribute)
        ax.axhline(level, color=colour, linestyle="--", linewidth=1.2, zorder=1)
        ax.text(weeks[-1], level, f" {name} {level:.2f}", color=colour, va="center", fontsize=8)

    ax.set_xlabel("surveillance week")
    ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    if season_labels is not None:
        ax.legend(loc="upper left", fontsize=8)
    ax.margins(x=0)
    return ax
