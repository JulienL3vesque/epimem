"""Threshold stability across an expanding set of seasons (mem's memstability).

Refit the model on the last 2, last 3, ..., all seasons and tabulate the key outputs for each
fit. Comparing rows shows whether the thresholds settle as history grows, i.e. how many seasons
are needed before the thresholds are stable. Pure diagnostic: it only re-runs mem_model on
larger subsets. Deterministic, reproduces R mem exactly.
"""

from dataclasses import dataclass

import numpy as np

from .model import mem_model

# Smallest window worth refitting: a single season has no spread to report.
MINIMUM_WINDOW = 2

# The 14 reported columns, in R's order. Each row of Stability.data holds these values.
STABILITY_COLUMNS = (
    "durationll", "duration", "durationul", "startll", "start", "startul",
    "percentagell", "percentage", "percentageul",
    "epidemic", "postepidemic", "medium", "high", "veryhigh",
)
COLUMN_COUNT = len(STABILITY_COLUMNS)


@dataclass(frozen=True)
class Stability:
    """One row per expanding-window fit."""
    counts: np.ndarray        # Seasons kept in each fit (the row labels).
    data: np.ndarray          # Shape (fit_count, COLUMN_COUNT); see STABILITY_COLUMNS.
    seasons_used: np.ndarray  # Shape (fit_count, season_count) bool: which seasons each fit used.
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
