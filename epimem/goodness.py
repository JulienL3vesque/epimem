"""Cross-validated goodness scores for a threshold setting (mem's memgoodness).

Leave-one-season-out evaluation of how well a threshold setting separates epidemic from
non-epidemic weeks:

    1. Leave one season out.
    2. Build thresholds from the remaining seasons.
    3. Classify every week of the left-out season as epidemic or not, using those thresholds.
    4. Compare that against the season's own gold standard, its MAP optimum.
    5. Sum every season's confusion matrix and read the diagnostic scores off the total.

Deterministic, so the numbers match the R reference exactly.
"""

from dataclasses import asdict, dataclass
from typing import NamedTuple

import numpy as np

from .model import mem_model
from .timing import map_curve, optimum_criterion

# Default detection grid, R's seq(1, 5, 0.1).
DETECTION_START = 1.0
DETECTION_STOP = 5.0
DETECTION_STEP = 0.1

# The fifteen reported values, in order; these are the scalar fields of Goodness.
METRIC_NAMES = (
    "weeks",
    "non_missing_weeks",
    "true_positives",
    "false_positives",
    "true_negatives",
    "false_negatives",
    "sensitivity",
    "specificity",
    "positive_predictive_value",
    "negative_predictive_value",
    "positive_likelihood_ratio",
    "negative_likelihood_ratio",
    "percent_agreement",
    "matthews_corrcoef",
    "youden_index",
)


class DiagnosticScores(NamedTuple):
    """Nine diagnostic-test scores derived from a confusion matrix."""

    sensitivity: float
    specificity: float
    positive_predictive_value: float
    negative_predictive_value: float
    positive_likelihood_ratio: float
    negative_likelihood_ratio: float
    percent_agreement: float
    matthews_corrcoef: float
    youden_index: float


@dataclass(frozen=True)
class ConfusionMatrix:
    """A 2x2 count of epidemic calls. The positive class is an epidemic week."""

    true_positives: int = 0      # epidemic weeks correctly flagged
    false_positives: int = 0     # non-epidemic weeks wrongly flagged
    true_negatives: int = 0      # non-epidemic weeks correctly left unflagged
    false_negatives: int = 0     # epidemic weeks missed

    def __add__(self, other: "ConfusionMatrix") -> "ConfusionMatrix":
        return ConfusionMatrix(
            self.true_positives + other.true_positives,
            self.false_positives + other.false_positives,
            self.true_negatives + other.true_negatives,
            self.false_negatives + other.false_negatives,
        )

    def scores(self) -> DiagnosticScores:
        """The nine diagnostic scores derived from these counts."""
        true_positives = np.float64(self.true_positives)
        false_positives = np.float64(self.false_positives)
        true_negatives = np.float64(self.true_negatives)
        false_negatives = np.float64(self.false_negatives)

        epidemic_weeks = true_positives + false_negatives    # weeks that were epidemic
        quiet_weeks = true_negatives + false_positives       # weeks that were non-epidemic
        alerts_fired = true_positives + false_positives       # weeks the thresholds flagged
        all_clears = true_negatives + false_negatives         # weeks the thresholds did not flag

        # A zero denominator gives NaN, or an infinite likelihood ratio; the float cast keeps
        # the arithmetic well defined instead of raising.
        with np.errstate(divide="ignore", invalid="ignore"):
            sensitivity = true_positives / epidemic_weeks
            specificity = true_negatives / quiet_weeks
            return DiagnosticScores(
                sensitivity=float(sensitivity),
                specificity=float(specificity),
                positive_predictive_value=float(true_positives / alerts_fired),
                negative_predictive_value=float(true_negatives / all_clears),
                positive_likelihood_ratio=float(sensitivity / (1 - specificity)),
                negative_likelihood_ratio=float((1 - sensitivity) / specificity),
                percent_agreement=float((true_positives + true_negatives) / (epidemic_weeks + quiet_weeks)),
                matthews_corrcoef=float(
                    (true_positives * true_negatives - false_positives * false_negatives)
                    / np.sqrt(alerts_fired * epidemic_weeks * quiet_weeks * all_clears)
                ),
                youden_index=float(sensitivity + specificity - 1),
            )


class SeasonResult(NamedTuple):
    """One held-out season's week counts and the confusion matrix it produced."""

    weeks: int
    non_missing_weeks: int
    confusion: ConfusionMatrix


@dataclass(frozen=True)
class Goodness:
    """Pooled, cross-validated scores for one threshold setting."""

    weeks: int
    non_missing_weeks: int
    true_positives: int
    false_positives: int
    true_negatives: int
    false_negatives: int
    sensitivity: float
    specificity: float
    positive_predictive_value: float
    negative_predictive_value: float
    positive_likelihood_ratio: float
    negative_likelihood_ratio: float
    percent_agreement: float
    matthews_corrcoef: float
    youden_index: float
    per_season: tuple = ()       # one ConfusionMatrix per held-out season

    @classmethod
    def from_pooled(cls, weeks: int, non_missing_weeks: int, confusion: ConfusionMatrix,
                    per_season: tuple = ()) -> "Goodness":
        """Build a report from pooled week counts and a pooled confusion matrix."""
        return cls(
            weeks=weeks,
            non_missing_weeks=non_missing_weeks,
            **asdict(confusion),
            **confusion.scores()._asdict(),
            per_season=per_season,
        )


def _default_detection_values() -> np.ndarray:
    """The default detection grid: 1.0, 1.1, ..., 5.0 (R's seq(1, 5, 0.1))."""
    value_count = int(round((DETECTION_STOP - DETECTION_START) / DETECTION_STEP)) + 1
    return DETECTION_START + DETECTION_STEP * np.arange(value_count)


def _nearest_other_seasons(season_count: int, target: int, pool_size: int) -> np.ndarray:
    """The pool_size seasons nearest the target, target excluded, ties broken toward the earlier
    season."""
    distance = np.arange(season_count) - target
    nearest_first = np.lexsort((distance, np.abs(distance)))   # by distance, then earlier season
    return np.sort(nearest_first[1: pool_size + 1])


def _predicted_epidemic_weeks(weekly_values, onset_threshold, post_epidemic_threshold) -> np.ndarray:
    """Which weeks the thresholds call epidemic, splitting the season at its peak.

    Before the peak, a week is epidemic if it is strictly above the onset threshold. From the
    peak on, if it is at or above the post-epidemic threshold.
    """
    peak_week = int(np.nanargmax(weekly_values)) + 1   # 1-indexed first peak
    is_epidemic = np.zeros(weekly_values.size, dtype=bool)
    is_epidemic[:peak_week] = weekly_values[:peak_week] > onset_threshold
    is_epidemic[peak_week:] = weekly_values[peak_week:] >= post_epidemic_threshold
    return is_epidemic


def _confusion_for_season(
    season_values, onset_threshold, post_epidemic_threshold, detection_values
) -> SeasonResult:
    """Score one held-out season against its own MAP optimum, pooled over the detection grid."""
    weekly_values = np.asarray(season_values, dtype=float).ravel()
    week_has_data = np.isfinite(weekly_values)
    week_number = np.arange(weekly_values.size) + 1   # 1-indexed week numbers

    # A missing onset threshold means nothing counts as epidemic before the peak.
    if not np.isfinite(onset_threshold):
        onset_threshold = np.inf
    if post_epidemic_threshold is None or not np.isfinite(post_epidemic_threshold):
        post_epidemic_threshold = onset_threshold

    predicted_epidemic = _predicted_epidemic_weeks(weekly_values, onset_threshold, post_epidemic_threshold)
    map_curve_table = map_curve(weekly_values)

    confusion = ConfusionMatrix()
    for detection_value in detection_values:
        _, start_week, end_week = optimum_criterion(map_curve_table, float(detection_value))
        gold_epidemic = (week_number >= start_week) & (week_number <= end_week)
        confusion += ConfusionMatrix(
            true_positives=int(np.sum(week_has_data & gold_epidemic & predicted_epidemic)),
            false_positives=int(np.sum(week_has_data & ~gold_epidemic & predicted_epidemic)),
            true_negatives=int(np.sum(week_has_data & ~gold_epidemic & ~predicted_epidemic)),
            false_negatives=int(np.sum(week_has_data & gold_epidemic & ~predicted_epidemic)),
        )

    return SeasonResult(weekly_values.size, int(week_has_data.sum()), confusion)


def mem_goodness(
    seasons,
    *,
    max_seasons: int | None = 10,
    type_threshold: int = 5,
    level_threshold: float = 0.95,
    tails_threshold: int = 1,
    type_intensity: int = 6,
    level_intensity=(0.40, 0.90, 0.975),
    tails_intensity: int = 1,
    method: int = 2,
    param: float = 2.8,
    n_values: int | None = None,
    use_t: bool = False,
    detection_values=None,
    min_seasons: int = 6,
) -> Goodness:
    """Leave-one-season-out goodness of a threshold setting (mem's memgoodness).

    Fit the thresholds on every season but one, score the season left out, then read the
    metrics off the pooled confusion matrix rather than by averaging per-season metrics.
    """
    season_matrix = np.asarray(seasons, dtype=float)
    if season_matrix.ndim != 2:
        raise ValueError("`seasons` must be 2-D: rows = weeks, columns = seasons")

    season_matrix = season_matrix[:, np.nansum(season_matrix, axis=0) > 0]   # drop empty seasons
    season_count = season_matrix.shape[1]
    if season_count <= 2:
        raise ValueError("memgoodness needs at least three seasons of non-zero data")
    if season_count < min_seasons:
        raise ValueError(f"memgoodness needs at least min_seasons={min_seasons} valid seasons")

    pool_size = season_count if max_seasons is None else max_seasons
    grid = _default_detection_values() if detection_values is None else detection_values
    def score_held_out(held_out: int) -> SeasonResult:
        model_seasons = _nearest_other_seasons(season_count, held_out, pool_size)
        thresholds = mem_model(
            season_matrix[:, model_seasons], max_seasons=pool_size,
            n_values=n_values, method=method, param=param,
            type_threshold=type_threshold, level_threshold=level_threshold, tails_threshold=tails_threshold,
            type_intensity=type_intensity, level_intensity=level_intensity, tails_intensity=tails_intensity,
            use_t=use_t,
        )
        return _confusion_for_season(
            season_matrix[:, held_out], thresholds.epidemic_onset, thresholds.post_epidemic, grid
        )

    results = [score_held_out(held_out) for held_out in range(season_count)]
    pooled = sum((result.confusion for result in results), ConfusionMatrix())

    return Goodness.from_pooled(
        weeks=sum(result.weeks for result in results),
        non_missing_weeks=sum(result.non_missing_weeks for result in results),
        confusion=pooled,
        per_season=tuple(result.confusion for result in results),
    )
