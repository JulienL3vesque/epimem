"""Choose the slope parameter from a season matrix - mem's roc.analysis and optimum.by.inspection.

roc_analysis scores each candidate slope by leave-one-season-out cross-validation (mem_goodness).
optimum_by_inspection scores each candidate slope against epidemics an analyst marked by hand
instead of by cross-validation. Both sweep a grid of slopes and, for each of eight criteria,
report the slope that scores best. Deterministic, so they match the R reference exactly.
"""

from dataclasses import dataclass

import numpy as np
from scipy.stats import rankdata

from .goodness import METRIC_NAMES, ConfusionMatrix, Goodness, mem_goodness
from .timing import map_curve, optimum_criterion

# Default slope grid, R's seq(1, 5, 0.1).
SLOPE_GRID_START = 1.0
SLOPE_GRID_STOP = 5.0
SLOPE_GRID_STEP = 0.1

# A floor of zero filters nothing, so every slope is kept.
NO_FLOOR = 0.0

# Sweep-table columns: the slope, then the fifteen Goodness metrics.
ROC_COLUMNS = ("value", *METRIC_NAMES)

# The eight criteria a winning slope is reported for.
CRITERIA = (
    "pos_likelihood",
    "neg_likelihood",
    "additive",
    "multiplicative",
    "mixed",
    "percent",
    "matthews",
    "youden",
)


@dataclass(frozen=True)
class SweepResult:
    """Result of a slope sweep: the winning slope per criterion, the ranking scores, and the table."""

    optimum: dict          # criterion -> best slope (NaN if the criterion is undecidable)
    rankings: dict         # criterion -> ranking score per surviving row
    table: np.ndarray      # one row per surviving slope; columns are ROC_COLUMNS
    columns: tuple = ROC_COLUMNS

    def best(self, criterion: str = "youden") -> float:
        """The chosen slope under one criterion (Youden's index by default)."""
        return self.optimum[criterion]


def _default_slope_grid() -> np.ndarray:
    """The default slope grid: 1.0, 1.1, ..., 5.0 (R's seq(1, 5, 0.1))."""
    value_count = int(round((SLOPE_GRID_STOP - SLOPE_GRID_START) / SLOPE_GRID_STEP)) + 1
    return SLOPE_GRID_START + SLOPE_GRID_STEP * np.arange(value_count)


def _validated_season_matrix(seasons, min_seasons: int) -> np.ndarray:
    """Drop empty seasons and check that enough remain to tune on."""
    season_matrix = np.asarray(seasons, dtype=float)
    if season_matrix.ndim != 2:
        raise ValueError("`seasons` must be 2-D: rows = weeks, columns = seasons")
    season_matrix = season_matrix[:, np.nansum(season_matrix, axis=0) > 0]
    season_count = season_matrix.shape[1]
    if season_count <= 2:
        raise ValueError("need at least three seasons of non-zero data")
    if season_count < min_seasons:
        raise ValueError(f"need at least min_seasons={min_seasons} valid seasons")
    return season_matrix


def _floor(value) -> float:
    """A missing sensitivity/specificity floor means no filtering, so return zero (as in R)."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return NO_FLOOR
    return float(value)


def _rank_best_first(values: np.ndarray) -> np.ndarray:
    """Rank so a larger value gets a smaller rank (best = rank 1); missing values rank last.
    Port of R's rank(-values, na.last = TRUE)."""
    negated = -np.asarray(values, dtype=float)
    return rankdata(np.where(np.isnan(negated), np.inf, negated), method="average")


def _rank_smallest_first(values: np.ndarray) -> np.ndarray:
    """Rank so a smaller value gets a smaller rank; missing values rank last (port of R's rank())."""
    array = np.asarray(values, dtype=float)
    return rankdata(np.where(np.isnan(array), np.inf, array), method="average")


def _winning_slope(ranking: np.ndarray, slope_grid: np.ndarray) -> float:
    """The slope with the best (lowest) rank, ties broken by the smallest slope.
    The ranking is over the possibly filtered rows but the grid it indexes is the full one; this
    matches R, which has the same mismatch. With the default floors nothing is filtered, so the
    two line up.
    """
    return float(slope_grid[int(np.argmin(ranking))])


def _criteria_optima(table: np.ndarray, slope_grid: np.ndarray) -> tuple[dict, dict]:
    """The winning slope and ranking score under each of the eight criteria."""
    column = {name: table[:, index] for index, name in enumerate(ROC_COLUMNS)}
    sensitivity = column["sensitivity"]
    specificity = column["specificity"]
    have_sensitivity_and_specificity = np.any(np.isfinite(sensitivity)) and np.any(np.isfinite(specificity))

    rankings: dict = {}
    optimum: dict = {}

    def record(name: str, ranking) -> None:
        rankings[name] = ranking
        optimum[name] = np.nan if ranking is None else _winning_slope(ranking, slope_grid)

    def record_single_metric(name: str, metric: str) -> None:
        values = column[metric]
        record(name, _rank_best_first(values) if np.any(np.isfinite(values)) else None)

    record(
        "additive",
        _rank_best_first(sensitivity) + _rank_best_first(specificity) if have_sensitivity_and_specificity else None,
    )
    record(
        "multiplicative",
        _rank_best_first(sensitivity * specificity) if have_sensitivity_and_specificity else None,
    )
    if have_sensitivity_and_specificity:
        # The "mixed" criterion rewards a slope whose sensitivity and specificity are both high
        # and close together. It sums three terms to minimize: the gap between them, the shortfall
        # from a perfect 1 + 1, and the distance to the top-left corner of an ROC plot.
        gap = np.abs(sensitivity - specificity)
        shortfall = 2 - sensitivity - specificity
        distance_to_corner = (1 - sensitivity) ** 2 + (1 - specificity) ** 2
        record(
            "mixed",
            _rank_smallest_first(gap) + _rank_smallest_first(shortfall) + _rank_smallest_first(distance_to_corner),
        )
    else:
        record("mixed", None)

    record_single_metric("pos_likelihood", "positive_likelihood_ratio")
    record_single_metric("neg_likelihood", "negative_likelihood_ratio")
    record_single_metric("percent", "percent_agreement")
    record_single_metric("matthews", "matthews_corrcoef")
    record_single_metric("youden", "youden_index")
    return optimum, rankings


def _sweep_row(slope: float, report: Goodness) -> list:
    """One sweep-table row: the slope, then the fifteen Goodness metrics in ROC_COLUMNS order."""
    return [slope, *(getattr(report, metric) for metric in METRIC_NAMES)]


def roc_analysis(
    seasons,
    *,
    param_values=None,
    min_seasons: int = 6,
    min_sensitivity=None,
    min_specificity=None,
    **goodness_kwargs,
) -> SweepResult:
    """Sweep the slope and report the best value under eight criteria (mem's roc.analysis).
    `**goodness_kwargs` pass through to mem_goodness (e.g. detection_values, type_threshold)."""
    season_matrix = _validated_season_matrix(seasons, min_seasons)
    slope_grid = _default_slope_grid() if param_values is None else np.asarray(param_values, dtype=float)

    table = np.array([
        _sweep_row(slope, mem_goodness(season_matrix, param=float(slope), min_seasons=min_seasons, **goodness_kwargs))
        for slope in slope_grid
    ])

    # Drop slopes below the sensitivity / specificity floors. A NaN fails the comparison, so
    # those rows drop too, which matches R's filter.
    column = {name: table[:, index] for index, name in enumerate(ROC_COLUMNS)}
    survives_floor = (column["specificity"] >= _floor(min_specificity)) & (column["sensitivity"] >= _floor(min_sensitivity))
    table = table[survives_floor]

    optimum, rankings = _criteria_optima(table, slope_grid)
    return SweepResult(optimum=optimum, rankings=rankings, table=table)


def _epidemic_window(start_week: int, end_week: int, week_count: int) -> np.ndarray:
    """Boolean mask of the weeks inside an epidemic [start, end]; out-of-range clamps to the edges."""
    first_week = start_week if 1 <= start_week <= week_count else 1
    last_week = end_week if 1 <= end_week <= week_count else week_count
    week_number = np.arange(week_count) + 1
    return (week_number >= first_week) & (week_number <= last_week)


def _inspection_confusion(season_values, marked_timing, automatic_timing):
    """Compare one season's automatic epidemic window against the analyst's marked one.
    Returns the season's week count, its non-missing week count, and the confusion matrix. The
    marked epidemic is the truth, the automatic timing is the prediction.
    """
    weekly_values = np.asarray(season_values, dtype=float).ravel()
    week_count = weekly_values.size
    week_has_data = np.isfinite(weekly_values)
    marked_epidemic = _epidemic_window(marked_timing[0], marked_timing[1], week_count)
    automatic_epidemic = _epidemic_window(automatic_timing[0], automatic_timing[1], week_count)

    confusion = ConfusionMatrix(
        true_positives=int(np.sum(week_has_data & marked_epidemic & automatic_epidemic)),
        false_positives=int(np.sum(week_has_data & ~marked_epidemic & automatic_epidemic)),
        true_negatives=int(np.sum(week_has_data & ~marked_epidemic & ~automatic_epidemic)),
        false_negatives=int(np.sum(week_has_data & marked_epidemic & ~automatic_epidemic)),
    )
    return week_count, int(week_has_data.sum()), confusion


def optimum_by_inspection(seasons, inspection_timings, *, param_values=None) -> SweepResult:
    """Tune the slope against analyst-marked epidemics (mem's optimum.by.inspection).
    `inspection_timings` is one (start_week, end_week) pair per season, 1-indexed, marked by eye.
    R reads these from mouse clicks; here they are a required argument. For each slope, the
    automatic epidemic timing is compared to the marked one, the confusion matrices are pooled
    across seasons, and the slopes are ranked the same eight ways as roc_analysis.
    """
    season_matrix = np.asarray(seasons, dtype=float)
    if season_matrix.ndim != 2:
        raise ValueError("`seasons` must be 2-D: rows = weeks, columns = seasons")
    season_count = season_matrix.shape[1]
    if len(inspection_timings) != season_count:
        raise ValueError("need one (start, end) inspection timing per season")
    slope_grid = _default_slope_grid() if param_values is None else np.asarray(param_values, dtype=float)

    # A season's MAP curve does not depend on the slope, so compute it once per season.
    season_map_curves = [map_curve(season_matrix[:, season]) for season in range(season_count)]

    rows = []
    for slope in slope_grid:
        total_weeks = 0
        total_non_missing_weeks = 0
        pooled = ConfusionMatrix()
        for season in range(season_count):
            _, automatic_start, automatic_end = optimum_criterion(season_map_curves[season], float(slope))
            weeks, non_missing_weeks, confusion = _inspection_confusion(
                season_matrix[:, season], inspection_timings[season], (automatic_start, automatic_end)
            )
            total_weeks += weeks
            total_non_missing_weeks += non_missing_weeks
            pooled += confusion
        report = Goodness.from_pooled(total_weeks, total_non_missing_weeks, pooled)
        rows.append(_sweep_row(slope, report))

    table = np.array(rows)
    optimum, rankings = _criteria_optima(table, slope_grid)
    return SweepResult(optimum=optimum, rankings=rankings, table=table)
