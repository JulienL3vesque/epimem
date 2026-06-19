"""Per-season epidemic timing (mem's fill.missing, calcular.map, calcular.optimo.criterio,
memtiming).

Given one season of weekly values, find the epidemic period and return the top
pre-epidemic, epidemic, and post-epidemic values that the threshold model pools.
"""

from dataclasses import dataclass
from typing import cast

import numpy as np
from scipy.optimize import brentq

# mem's smoother (R's sm.regression) is a Gaussian local-linear fit. The optimum uses a fixed
# bandwidth of one week. Gap-filling instead selects the bandwidth that gives the smoother
# SMOOTHER_TARGET_DF effective degrees of freedom.
DEFAULT_SMOOTHING_BANDWIDTH = 1.0
SMOOTHER_TARGET_DF = 6.0

# The slope of the percent-captured curve must fall below this before the epidemic is declared
# over. mem calls it the criterion; 2.8 is its default.
DEFAULT_OPTIMUM_CRITERION = 2.8

# Bracketing interval for the bandwidth root-find, as multiples of the data span.
BANDWIDTH_SEARCH_FLOOR = 0.3
BANDWIDTH_SEARCH_CEILING_FACTOR = 3.0
BANDWIDTH_SEARCH_TOLERANCE = 1e-8

# Column positions in a MAP table row: [duration, percent_captured, peak_sum, start, end].
DURATION_COLUMN = 0
PERCENT_CAPTURED_COLUMN = 1
PEAK_SUM_COLUMN = 2
START_WEEK_COLUMN = 3
END_WEEK_COLUMN = 4


def max_n_values(values, n_max: int = 1) -> np.ndarray:
    """The n_max highest finite values, descending, padded with NaN (mem's maxnvalores)."""
    finite_values = np.asarray(values, dtype=float).ravel()
    finite_values = finite_values[np.isfinite(finite_values)]
    descending = np.sort(finite_values)[::-1]
    result = np.full(n_max, np.nan)
    take = min(n_max, descending.size)
    result[:take] = descending[:take]
    return result


def _smoother_matrix(eval_points, sample_points, bandwidth: float) -> np.ndarray:
    """Weight matrix for the local-linear smoother (port of sm::sm.weight).

    smoothed = S @ observations, and trace(S) is the smoother's degrees of freedom.
    """
    offset = (
        np.asarray(eval_points, dtype=float)[:, None]
        - np.asarray(sample_points, dtype=float)[None, :]
    )
    kernel = np.exp(-0.5 * (offset / bandwidth) ** 2)

    # Gaussian-weighted moments of the offsets, the local-linear normal equations.
    weight_sum = kernel.sum(1)[:, None]
    weighted_offset = (kernel * offset).sum(1)[:, None]
    weighted_offset_squared = (kernel * offset ** 2).sum(1)[:, None]

    return (
        kernel
        * (weighted_offset_squared - offset * weighted_offset)
        / (weighted_offset_squared * weight_sum - weighted_offset ** 2)
    )


def _bandwidth_for_df(sample_points, target_df: float = SMOOTHER_TARGET_DF) -> float:
    """Bandwidth whose smoother has target_df degrees of freedom (R's sm::h.select 'df' rule, df=6)."""
    positions = np.asarray(sample_points, dtype=float)

    # Degrees of freedom = trace(S); it decreases as the bandwidth widens, so solve for the target.
    def trace_gap(bandwidth: float) -> float:
        return float(np.trace(_smoother_matrix(positions, positions, bandwidth))) - target_df

    data_span = positions.max() - positions.min()
    return cast(float, brentq(
        trace_gap,
        BANDWIDTH_SEARCH_FLOOR,
        data_span * BANDWIDTH_SEARCH_CEILING_FACTOR,
        xtol=BANDWIDTH_SEARCH_TOLERANCE,
    ))


def smooth_local_linear(values, bandwidth: float) -> np.ndarray:
    """Local-linear smoother (mem's suavizado, R's sm.regression).

    For each index, fit a line through the nearby points weighted by a Gaussian kernel and take
    its value, which reduces short-term noise. Negative results are set to 0. The optimum uses
    bandwidth 1.
    """
    weekly_values = np.asarray(values, dtype=float)
    positions = np.arange(1, weekly_values.size + 1, dtype=float)
    observed = np.isfinite(weekly_values)
    if observed.sum() <= 1:
        return weekly_values.copy()
    fitted = _smoother_matrix(positions, positions[observed], bandwidth) @ weekly_values[observed]
    fitted[fitted < 0] = 0.0
    return fitted


def fill_missing(values) -> np.ndarray:
    """Fill interior gaps with the kernel smoother (mem's fill.missing).

    Interior NaNs are filled with sm.regression at an auto-selected bandwidth (df=6).
    Leading and trailing NaNs and complete seasons are left unchanged.
    """
    weekly_values = np.asarray(values, dtype=float).copy()
    observed_indices = np.where(np.isfinite(weekly_values))[0]
    if observed_indices.size < 2:
        return weekly_values

    # An interior gap is a NaN between the first and last observed week. Gaps before the first
    # or after the last observation are left unchanged.
    first_observed, last_observed = observed_indices[0], observed_indices[-1]
    is_interior_gap = np.zeros(weekly_values.size, dtype=bool)
    interior_span = slice(first_observed, last_observed + 1)
    is_interior_gap[interior_span] = ~np.isfinite(weekly_values[interior_span])
    if not is_interior_gap.any():
        return weekly_values

    observed_positions = (observed_indices + 1).astype(float)         # 1-indexed, as in R
    bandwidth = _bandwidth_for_df(observed_positions)
    all_positions = np.arange(1, weekly_values.size + 1, dtype=float)
    fitted = _smoother_matrix(all_positions, observed_positions, bandwidth) @ weekly_values[observed_indices]
    weekly_values[is_interior_gap] = fitted[is_interior_gap]
    return weekly_values


def map_curve(values) -> np.ndarray:
    """MAP curve (mem's calcular.map): for each window length, the heaviest consecutive run.

    One row per duration 0..n_weeks: [duration, pct_of_total, max_sum, start, end], 1-indexed.
    """
    weekly_values = np.asarray(values, dtype=float)
    week_count = weekly_values.size
    season_total = np.nansum(weekly_values)

    rows = []
    for duration in range(1, week_count + 1):
        best_sum = -np.inf
        best_start = 1
        # Slide a window of this length across the season and keep the heaviest run, breaking
        # ties toward the earliest window.
        for start_week in range(1, week_count - duration + 2):            # 1-indexed start weeks
            window = weekly_values[start_week - 1: start_week - 1 + duration]
            window_sum = np.nansum(window)
            if window_sum > best_sum:
                best_sum = window_sum
                best_start = start_week
        percent_captured = 100 * best_sum / season_total if season_total != 0 else 0.0
        rows.append([duration, percent_captured, best_sum, best_start, best_start + duration - 1])

    zero_row = [0, 0, 0, 0, 0]                                            # R prepends a zero row
    return np.vstack([zero_row, np.array(rows, dtype=float)])


def optimum_criterion(map_table: np.ndarray, criterion: float = DEFAULT_OPTIMUM_CRITERION) -> tuple[int, int, int]:
    """Optimal epidemic duration (mem's calcular.optimo.criterio, the default method).

    Take the first duration where the smoothed slope of the percent-captured curve drops below
    `criterion`, minus one. Returns (duration, start, end), 1-indexed.
    """
    durations = map_table[:, DURATION_COLUMN]
    percent_captured = map_table[:, PERCENT_CAPTURED_COLUMN]

    # Each added epidemic week captures less of the season than the previous one. Once that
    # marginal gain falls below the criterion the epidemic is treated as over.
    smoothed = smooth_local_linear(percent_captured, bandwidth=DEFAULT_SMOOTHING_BANDWIDTH)
    marginal_gain = np.diff(smoothed)
    marginal_gain_durations = durations[1:]

    gain_below_criterion = marginal_gain < criterion
    if np.any(gain_below_criterion):
        first_below = int(np.min(marginal_gain_durations[gain_below_criterion]))
        duration = max(first_below - 1, 1)
    else:
        duration = int(np.max(marginal_gain_durations))

    chosen_row = map_table[map_table[:, DURATION_COLUMN] == duration][0]
    return duration, int(chosen_row[START_WEEK_COLUMN]), int(chosen_row[END_WEEK_COLUMN])


@dataclass(frozen=True)
class EpidemicTiming:
    """One season's epidemic period and the extreme values the thresholds use."""
    start: int                  # 1-indexed first epidemic week
    end: int                    # 1-indexed last epidemic week
    duration: int
    percent: float              # percent of the season's total inside the epidemic window
    pre_epi: np.ndarray         # top pre-epidemic values
    epi: np.ndarray             # top epidemic values
    post_epi: np.ndarray        # top post-epidemic values


def mem_timing(season_values, n_values: int = 5, method: int = 2, param: float = DEFAULT_OPTIMUM_CRITERION) -> EpidemicTiming:
    """Find one season's epidemic period (mem's memtiming)."""
    if method != 2:
        raise NotImplementedError("only the default criterion optimum (method=2) is ported")

    weekly_values = fill_missing(season_values)
    curve = map_curve(weekly_values)
    duration, start, end = optimum_criterion(curve, param)

    duration_row = curve[:, DURATION_COLUMN] == duration
    return EpidemicTiming(
        start=start,
        end=end,
        duration=duration,
        percent=float(curve[duration_row, PERCENT_CAPTURED_COLUMN][0]),
        pre_epi=max_n_values(weekly_values[: start - 1], n_values),
        epi=max_n_values(weekly_values[start - 1: end], n_values),
        post_epi=max_n_values(weekly_values[end:], n_values),
    )
