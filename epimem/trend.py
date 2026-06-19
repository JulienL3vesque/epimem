"""Real-time rising / falling signal, mem's memtrend.

From the week-over-week change inside each past epidemic, pool the rises and the falls
separately and put a one-sided confidence interval on each. That yields two cut-points for
the latest observed change:

    ascending   a rise above this means the epidemic is clearly starting or climbing.
    descending  a fall below this (a negative number) means it is clearly declining.

Deterministic on the default arithmetic confidence interval (type 1), so it matches the R
reference exactly.
"""

from dataclasses import dataclass

import numpy as np

from .confidence import confidence_interval

from .timing import mem_timing

# A week-over-week change of exactly zero is neither a clear rise nor a clear fall, so it is
# counted in both pools; these bounds keep that "0 belongs to both sides" rule explicit.
RISE_FLOOR = 0.0
FALL_CEILING = 0.0

# Both thresholds come from one-sided confidence intervals: we only care about how far the
# change could plausibly reach in a single direction.
ONE_SIDED = 1


@dataclass(frozen=True)
class TrendThresholds:
    """Week-over-week change cut-points (NaN if a pool of that sign is empty)."""

    ascending: float       # A rise above this means the epidemic is clearly starting or rising.
    descending: float      # A change below this (negative) means it is clearly declining.


def mem_trend(
    seasons,
    *,
    confidence_interval_type: int = 1,
    level: float = 0.95,
    use_t: bool = False,
    method: int = 2,
    param: float = 2.8,
) -> TrendThresholds:
    """Ascending / descending week-over-week thresholds (mem's memtrend).

    For each season, keep the first differences that fall inside its epidemic window,
    [start .. min(week_count, end + 1)], then split them into rises (>= 0) and falls (<= 0);
    a change of exactly zero counts in both. The ascending threshold is the lower bound of the
    rises' one-sided confidence interval, the descending threshold the upper bound of the
    falls'.
    """
    season_matrix = np.asarray(seasons, dtype=float)
    if season_matrix.ndim != 2:
        raise ValueError("`seasons` must be 2-D: rows = weeks, columns = seasons")
    week_count, season_count = season_matrix.shape

    # First differences down each column. The first week has no prior week to compare against,
    # so its difference is left missing.
    weekly_change = np.full_like(season_matrix, np.nan)
    weekly_change[1:] = season_matrix[1:] - season_matrix[:-1]

    # Keep only the changes that occur inside each season's own epidemic window; everything
    # outside the epidemic is irrelevant to a rising / falling signal.
    epidemic_change = np.full_like(season_matrix, np.nan)
    for season in range(season_count):
        timing = mem_timing(season_matrix[:, season], method=method, param=param)
        first_row = timing.start - 1                       # 1-indexed start to 0-indexed first row.
        last_row = min(week_count, timing.end + 1)         # Include the first post-epidemic drop.
        epidemic_change[first_row:last_row, season] = weekly_change[first_row:last_row, season]

    rises = np.where(epidemic_change >= RISE_FLOOR, epidemic_change, np.nan)
    falls = np.where(epidemic_change <= FALL_CEILING, epidemic_change, np.nan)
    ascending = confidence_interval(
        rises.ravel(), level=level, confidence_interval_type=confidence_interval_type, tails=ONE_SIDED, use_t=use_t
    ).lower
    descending = confidence_interval(
        falls.ravel(), level=level, confidence_interval_type=confidence_interval_type, tails=ONE_SIDED, use_t=use_t
    ).upper
    return TrendThresholds(ascending=float(ascending), descending=float(descending))
