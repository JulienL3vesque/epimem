"""Threshold diagnostics across seasons - mem's memstability and memevolution.

mem_stability refits the model on the last 2, last 3, ..., all seasons, so comparing rows shows
whether the thresholds settle as history grows. mem_evolution refits season by season - on the
prior seasons ("sequential", a real-time backtest) or on every other season ("cross",
leave-one-out) - so you can see how the thresholds would have looked at each point in time. Both
only re-run mem_model on subsets; deterministic, reproducing R mem.
"""

from dataclasses import dataclass

import numpy as np

from .goodness import _nearest_other_seasons
from .model import mem_model

# Smallest window worth refitting: a single season has no spread to report.
MINIMUM_WINDOW = 2

# The 14 reported columns, in R's order. Each row of Stability.data / Evolution.data holds these.
STABILITY_COLUMNS = (
    "durationll", "duration", "durationul", "startll", "start", "startul",
    "percentagell", "percentage", "percentageul",
    "epidemic", "postepidemic", "medium", "high", "veryhigh",
)
COLUMN_COUNT = len(STABILITY_COLUMNS)

# The two memevolution windowing rules.
EVOLUTION_METHODS = ("sequential", "cross")


@dataclass(frozen=True)
class Stability:
    """One row per expanding-window fit."""
    counts: np.ndarray        # Seasons kept in each fit (the row labels).
    data: np.ndarray          # Shape (fit_count, COLUMN_COUNT); see STABILITY_COLUMNS.
    seasons_used: np.ndarray  # Shape (fit_count, season_count) bool: which seasons each fit used.
    columns: tuple = STABILITY_COLUMNS


@dataclass(frozen=True)
class Evolution:
    """One row per season-by-season refit; the last row is the upcoming ("next") season."""
    counts: np.ndarray        # Seasons kept in each fit.
    data: np.ndarray          # Shape (row_count, COLUMN_COUNT); see STABILITY_COLUMNS.
    seasons_used: np.ndarray  # Shape (row_count, season_count) bool: which seasons each fit used.
    method: str               # "sequential" or "cross".
    columns: tuple = STABILITY_COLUMNS


def mem_stability(seasons, **model_kwargs) -> Stability:
    """Refit on the last 2, 3, ..., all seasons and tabulate the outputs (mem's memstability).

    `**model_kwargs` are forwarded to mem_model (e.g. type_threshold, param). The recency cap
    is forced off, since each window is already the chosen set, matching R's i.seasons = NA.
    """
    season_matrix = np.asarray(seasons, dtype=float)
    if season_matrix.ndim != 2:
        raise ValueError("`seasons` must be 2-D: rows = weeks, columns = seasons")
    season_count = season_matrix.shape[1]
    if season_count < MINIMUM_WINDOW:
        return Stability(
            np.array([]),
            np.empty((0, COLUMN_COUNT)),
            np.empty((0, season_count), dtype=bool),
        )

    rows, counts, membership = [], [], []
    for window in range(MINIMUM_WINDOW, season_count + 1):
        # Take the `window` most recent seasons (the rightmost columns) for this fit.
        used_seasons = np.arange(season_count - window, season_count)
        model = mem_model(season_matrix[:, used_seasons], max_seasons=None, **model_kwargs)
        rows.append([
            *model.duration_confidence_interval,
            *model.start_confidence_interval,
            *model.percent_confidence_interval,
            model.epidemic_onset, model.post_epidemic,
            model.medium, model.high, model.very_high,
        ])
        counts.append(model.n_seasons)
        season_was_used = np.zeros(season_count, dtype=bool)
        season_was_used[used_seasons] = True
        membership.append(season_was_used)

    return Stability(
        counts=np.array(counts),
        data=np.array(rows, dtype=float),
        seasons_used=np.array(membership),
    )


def mem_evolution(seasons, *, evolution_seasons: int | None = 10, method: str = "sequential",
                  **model_kwargs) -> Evolution:
    """How the thresholds would have evolved season by season (mem's memevolution).

    For each season the thresholds are rebuilt from a window of other seasons, set by `method`:

      "sequential"  use only the seasons BEFORE it - a real-time backtest. Rows run from the third
                    season onward, since two prior seasons are the minimum mem_model needs.
      "cross"       use every other season (leave-one-out), nearest seasons first; one row per
                    season.

    Either way the final row, the upcoming ("next") season, uses the most recent
    `evolution_seasons` seasons. `evolution_seasons` caps how many seasons feed each fit (R's
    i.evolution.seasons; None means use all available). Columns match mem_stability
    (STABILITY_COLUMNS); `**model_kwargs` pass through to mem_model.
    """
    season_matrix = np.asarray(seasons, dtype=float)
    if season_matrix.ndim != 2:
        raise ValueError("`seasons` must be 2-D: rows = weeks, columns = seasons")
    if method not in EVOLUTION_METHODS:
        raise ValueError(f"method must be one of {EVOLUTION_METHODS}")
    season_count = season_matrix.shape[1]
    if season_count < MINIMUM_WINDOW:
        return Evolution(
            np.array([]),
            np.empty((0, COLUMN_COUNT)),
            np.empty((0, season_count), dtype=bool),
            method,
        )

    cap = season_count if evolution_seasons is None else evolution_seasons

    if method == "sequential":
        # Each season from the third on is modelled from the seasons before it; the last is "next".
        windows = [np.arange(max(0, upto - 1 - cap), upto - 1)
                   for upto in range(MINIMUM_WINDOW + 1, season_count + 2)]
    else:
        # Each season is modelled from every other season (leave-one-out), nearest first; then "next".
        windows = [_nearest_other_seasons(season_count, target, cap)
                   for target in range(season_count)]
        windows.append(np.arange(max(0, season_count - cap), season_count))

    rows, counts, membership = [], [], []
    for used_seasons in windows:
        model = mem_model(season_matrix[:, used_seasons], max_seasons=None, **model_kwargs)
        rows.append([
            *model.duration_confidence_interval,
            *model.start_confidence_interval,
            *model.percent_confidence_interval,
            model.epidemic_onset, model.post_epidemic,
            model.medium, model.high, model.very_high,
        ])
        counts.append(model.n_seasons)
        season_was_used = np.zeros(season_count, dtype=bool)
        season_was_used[used_seasons] = True
        membership.append(season_was_used)

    return Evolution(
        counts=np.array(counts),
        data=np.array(rows, dtype=float),
        seasons_used=np.array(membership),
        method=method,
    )
