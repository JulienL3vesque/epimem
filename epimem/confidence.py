"""Confidence intervals (a centre estimate plus a margin) used to set thresholds.

Each function returns ConfidenceInterval(lower, centre, upper); MEM reads .upper for its
thresholds. The type numbers match R's mem:

    1  arithmetic mean confidence interval         (margin = sd / sqrt(n))
    2  geometric mean confidence interval
    3  median confidence interval                  (two-sided order-statistic / sign-test)
    5  arithmetic prediction confidence interval   (margin = sd)
    6  geometric prediction confidence interval

Type 4 (bootstrap) is not ported: it depends on R's random-number stream, so it cannot be
reproduced exactly, and no default uses it.
"""

from typing import NamedTuple

import numpy as np
from scipy.stats import binom, norm, t as t_distribution

# The median is the 50th percentile, the only quantile this family is asked to bracket.
MEDIAN_PROBABILITY = 0.5

# A geometric interval works on the logs, and a zero value has no log. Shift every value up by
# one before logging and undo the shift after exponentiating, as mem does.
ZERO_SAFE_LOG_SHIFT = 1


class ConfidenceInterval(NamedTuple):
    lower: float
    centre: float
    upper: float


def _finite_values(values) -> np.ndarray:
    """Non-NaN, non-Inf values, flattened to a 1-D float array. mem drops NA/Inf too."""
    array = np.asarray(values, dtype=float).ravel()
    return array[np.isfinite(array)]


def _margin_factor(level: float, tails: int, value_count: int, use_t: bool) -> float:
    """Margin width in standard errors: a normal quantile, or a t quantile when use_t."""
    upper_tail_probability = 1 - (1 - level) / tails
    if use_t:
        prediction_inflation = np.sqrt(1 + 1 / value_count)
        return float(t_distribution.ppf(upper_tail_probability, value_count - 1) * prediction_inflation)
    return float(norm.ppf(upper_tail_probability))


def _arithmetic(values, level: float, tails: int, use_t: bool, prediction: bool) -> ConfidenceInterval:
    """mean +/- factor * spread. The spread is the full sd for a prediction interval (type 5) or
    the standard error sd / sqrt(n) for an interval on the mean (type 1)."""
    data = _finite_values(values)
    if data.size == 0:
        return ConfidenceInterval(np.nan, np.nan, np.nan)
    # A single value has no spread to estimate, so the interval collapses onto it, as in R.
    if data.size < 2:
        return ConfidenceInterval(np.nan, float(data.mean()), np.nan)
    factor = _margin_factor(level, tails, data.size, use_t)
    mean = float(data.mean())
    standard_deviation = float(data.std(ddof=1))
    spread = standard_deviation if prediction else standard_deviation / np.sqrt(data.size)
    return ConfidenceInterval(mean - factor * spread, mean, mean + factor * spread)


def _geometric(values, level: float, tails: int, use_t: bool, prediction: bool) -> ConfidenceInterval:
    """The arithmetic interval computed on the logs, then mapped back through exp().

    Suited to skewed counts that bunch near zero. Uses log(x), or log(x + 1) - 1 when any value
    is zero, as mem does.
    """
    data = _finite_values(values)
    if data.size == 0:
        return ConfidenceInterval(np.nan, np.nan, np.nan)
    shift = 0 if np.all(data != 0) else ZERO_SAFE_LOG_SHIFT
    lower, centre, upper = _arithmetic(np.log(data + shift), level, tails, use_t, prediction)
    return ConfidenceInterval(np.exp(lower) - shift, np.exp(centre) - shift, np.exp(upper) - shift)


def _median(values, level: float, tails: int, use_t: bool, prediction: bool) -> ConfidenceInterval:
    """Two-sided confidence interval for the median, from order statistics (the sign test).

    The lower and upper limits are the values whose ranks come from the binomial(n, 0.5)
    distribution, the standard nonparametric median interval. This replaces R's interpolated
    EnvStats version with a few lines, at the cost of not exactly matching R on the
    typical-season summary columns, which nothing downstream depends on. use_t and prediction
    do not apply.
    """
    if tails != 2:
        raise NotImplementedError("median confidence interval (type 3) is two-sided only")
    data = _finite_values(values)
    if data.size == 0:
        return ConfidenceInterval(np.nan, np.nan, np.nan)

    ordered_values = np.sort(data)
    value_count = ordered_values.size
    centre = float(np.median(ordered_values))

    # The limits are a symmetric pair of order statistics; their ranks come from the binomial.
    tail_probability = (1 - level) / 2
    lower_rank = max(1, int(binom.ppf(tail_probability, value_count, MEDIAN_PROBABILITY)))
    upper_rank = min(value_count, value_count + 1 - lower_rank)
    return ConfidenceInterval(
        float(ordered_values[lower_rank - 1]), centre, float(ordered_values[upper_rank - 1])
    )


# Each type maps to the function that computes it and whether it is a prediction interval.
_CONFIDENCE_INTERVAL_TYPES = {
    1: (_arithmetic, False),
    2: (_geometric, False),
    3: (_median, False),
    5: (_arithmetic, True),
    6: (_geometric, True),
}


def confidence_interval(values, level: float = 0.95, confidence_interval_type: int = 1,
                        tails: int = 2, use_t: bool = False) -> ConfidenceInterval:
    """Select a confidence-interval type by its number (see the module docstring)."""
    if not 0 < level < 1:
        raise ValueError(f"level must be in (0, 1), got {level}")
    if confidence_interval_type not in _CONFIDENCE_INTERVAL_TYPES:
        raise ValueError(
            f"confidence_interval_type must be one of {sorted(_CONFIDENCE_INTERVAL_TYPES)} "
            "(4=bootstrap not ported)"
        )
    compute, prediction = _CONFIDENCE_INTERVAL_TYPES[confidence_interval_type]
    return compute(values, level, tails, use_t, prediction)
