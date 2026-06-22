"""Cross-season MEM thresholds - mem's memmodel (the threshold outputs).

From several past seasons, learn the epidemic-onset and medium/high/very-high thresholds.
MemThresholds.level_of() then grades a weekly value against them.
"""

from dataclasses import dataclass

import numpy as np

from .confidence import confidence_interval
from .timing import mem_timing

# mem pools roughly thirty extreme values across all seasons. When the caller does not give a
# per-season count, divide that total over the seasons available.
TARGET_POOLED_VALUE_COUNT = 30

# The three intensity levels are ordered medium, high, very high. Name those positions so the
# thresholds are never read from `level_intensity` by a bare index.
MEDIUM_LEVEL = 0
HIGH_LEVEL = 1
VERY_HIGH_LEVEL = 2
INTENSITY_LEVEL_COUNT = 3

# R defaults for the three intensity thresholds: 40%, 90% and 97.5% of the epidemic peaks.
DEFAULT_INTENSITY_LEVELS = (0.40, 0.90, 0.975)

# The typical-duration / start-week / %-covered confidence intervals are always two-sided.
TWO_SIDED = 2


@dataclass(frozen=True)
class MemThresholds:
    """The thresholds for one indicator, plus the settings they were built from."""
    epidemic_onset: float
    post_epidemic: float
    medium: float
    high: float
    very_high: float
    n_seasons: int
    n_values_per_season: int
    epidemic_durations: list[int]
    intensity_levels: tuple = DEFAULT_INTENSITY_LEVELS
    # duration / start / percent intervals: (lower, centre, upper) triples, filled in by mem_model.
    duration_confidence_interval: tuple[float, ...] = ()
    start_confidence_interval: tuple[float, ...] = ()
    percent_confidence_interval: tuple[float, ...] = ()

    def level_of(self, value: float) -> str:
        """Activity level for one weekly value."""
        if not np.isfinite(value):
            return "no data"
        if value >= self.very_high:
            return "very high"
        if value >= self.high:
            return "high"
        if value >= self.medium:
            return "medium"
        if value >= self.epidemic_onset:
            return "low (epidemic started)"
        return "baseline"


def mem_model(
    seasons,
    *,
    n_values: int | None = None,            # i.n.max ; None -> round(30 / n_seasons)
    method: int = 2,                        # i.method
    param: float = 2.8,                     # i.param
    type_threshold: int = 5,                # i.type.threshold
    level_threshold: float = 0.95,          # i.level.threshold
    tails_threshold: int = 1,               # i.tails.threshold
    type_intensity: int = 6,                # i.type.intensity
    level_intensity = DEFAULT_INTENSITY_LEVELS,   # i.level.intensity
    tails_intensity: int = 1,               # i.tails.intensity
    type_other: int = 3,                    # i.type.other (duration/start/% CIs; 3 = median)
    level_other: float = 0.95,              # i.level.other
    use_t: bool = False,                    # i.use.t
    max_seasons: int | None = 10,           # i.seasons
) -> MemThresholds:
    """Learn the thresholds from several seasons (mem memmodel).

    `seasons` is 2-D: rows = surveillance weeks, columns = seasons. Defaults match R, so the
    result lines up with memmodel's epidemic / intensity thresholds.
    """
    season_matrix = np.asarray(seasons, dtype=float)
    if season_matrix.ndim != 2:
        raise ValueError("`seasons` must be 2-D: rows = weeks, columns = seasons")
    if len(level_intensity) != INTENSITY_LEVEL_COUNT:
        raise ValueError("`level_intensity` must have three levels (medium, high, very high)")

    # Drop seasons with no signal, then keep only the most recent max_seasons.
    season_matrix = season_matrix[:, np.nansum(season_matrix, axis=0) > 0]
    if max_seasons and season_matrix.shape[1] > max_seasons:
        season_matrix = season_matrix[:, -max_seasons:]

    n_seasons = season_matrix.shape[1]
    if n_seasons < 2:
        raise ValueError("mem_model needs at least two seasons of non-zero data")

    if n_values is None:
        n_per_season = max(1, int(round(TARGET_POOLED_VALUE_COUNT / n_seasons)))
    else:
        n_per_season = int(n_values)

    season_timings = [
        mem_timing(season_matrix[:, season], n_values=n_per_season, method=method, param=param)
        for season in range(n_seasons)
    ]

    # Pool every season's top pre-epidemic, epidemic and post-epidemic values into one sample each.
    pre_epidemic_values = np.concatenate([timing.pre_epi for timing in season_timings])
    post_epidemic_values = np.concatenate([timing.post_epi for timing in season_timings])
    epidemic_values = np.concatenate([timing.epi for timing in season_timings])

    # The onset and post-epidemic thresholds are the upper end of a one-sided confidence interval
    # on those pools.
    onset_threshold = confidence_interval(
        pre_epidemic_values, level_threshold, type_threshold, tails_threshold, use_t
    ).upper
    post_epidemic_threshold = confidence_interval(
        post_epidemic_values, level_threshold, type_threshold, tails_threshold, use_t
    ).upper

    # Each intensity threshold is the upper end of a one-sided confidence interval on the pooled
    # epidemic peaks.
    intensity_thresholds = [
        confidence_interval(epidemic_values, level, type_intensity, tails_intensity, use_t).upper
        for level in level_intensity
    ]

    # Summarise the typical epidemic length, start week and %-of-season-covered, each as a
    # two-sided median confidence interval across the seasons. mem rounds the start interval to
    # whole weeks (floor / round / ceil); length and percent stay raw.
    durations = np.array([timing.duration for timing in season_timings], dtype=float)
    start_weeks = np.array([timing.start for timing in season_timings], dtype=float)
    percents_covered = np.array([timing.percent for timing in season_timings], dtype=float)
    duration_interval = confidence_interval(durations, level=level_other, confidence_interval_type=type_other, tails=TWO_SIDED)
    start_interval = confidence_interval(start_weeks, level=level_other, confidence_interval_type=type_other, tails=TWO_SIDED)
    percent_interval = confidence_interval(percents_covered, level=level_other, confidence_interval_type=type_other, tails=TWO_SIDED)

    return MemThresholds(
        epidemic_onset=onset_threshold,
        post_epidemic=post_epidemic_threshold,
        medium=intensity_thresholds[MEDIUM_LEVEL],
        high=intensity_thresholds[HIGH_LEVEL],
        very_high=intensity_thresholds[VERY_HIGH_LEVEL],
        n_seasons=n_seasons,
        n_values_per_season=n_per_season,
        epidemic_durations=[timing.duration for timing in season_timings],
        intensity_levels=tuple(level_intensity),
        duration_confidence_interval=(duration_interval.lower, duration_interval.centre, duration_interval.upper),
        start_confidence_interval=(
            float(np.floor(start_interval.lower)),
            float(round(start_interval.centre)),
            float(np.ceil(start_interval.upper)),
        ),
        percent_confidence_interval=(percent_interval.lower, percent_interval.centre, percent_interval.upper),
    )


def _intensity_label(level: float) -> str:
    """'40', '90', '97.5' - R's as.character(round(level*100, 1)), no trailing '.0'."""
    return f"{round(level * 100, 1):g}"


def mem_intensity(model: MemThresholds) -> dict[str, float]:
    """The four labelled cut-points (mem memintensity).

    Packages the onset plus the three intensity uppers already in `model`, labelled by level
    (e.g. 'Medium (40%)'). To grade a value, use model.level_of.
    """
    intensity_names = ("Medium", "High", "Very high")
    intensity_values = (model.medium, model.high, model.very_high)
    cut_points = {"Epidemic": model.epidemic_onset}
    for name, level, value in zip(intensity_names, model.intensity_levels, intensity_values):
        cut_points[f"{name} ({_intensity_label(level)}%)"] = value
    return cut_points
