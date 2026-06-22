"""epimem - Python port of the Moving Epidemic Method, from the R `mem` package.

    from epimem import mem_model
    th = mem_model(season_matrix)   # rows = weeks, cols = seasons
    th.level_of(value)              # 'baseline' ... 'very high'
"""
from .confidence import ConfidenceInterval, confidence_interval
from .timing import (
    EpidemicTiming,
    mem_timing,
    map_curve,
    optimum_criterion,
    smooth_local_linear,
    fill_missing,
    max_n_values,
)
from .model import MemThresholds, mem_model, mem_intensity
from .goodness import ConfusionMatrix, DiagnosticScores, Goodness, mem_goodness, METRIC_NAMES
from .tuning import SweepResult, roc_analysis, optimum_by_inspection, CRITERIA, ROC_COLUMNS
from .trend import TrendThresholds, mem_trend
from .stability import Evolution, Stability, mem_evolution, mem_stability, STABILITY_COLUMNS

__all__ = [
    "mem_model",
    "MemThresholds",
    "mem_intensity",
    "mem_trend",
    "TrendThresholds",
    "mem_stability",
    "Stability",
    "STABILITY_COLUMNS",
    "mem_evolution",
    "Evolution",
    "mem_timing",
    "EpidemicTiming",
    "confidence_interval",
    "ConfidenceInterval",
    "map_curve",
    "optimum_criterion",
    "smooth_local_linear",
    "fill_missing",
    "max_n_values",
    "mem_goodness",
    "Goodness",
    "ConfusionMatrix",
    "DiagnosticScores",
    "METRIC_NAMES",
    "roc_analysis",
    "optimum_by_inspection",
    "SweepResult",
    "CRITERIA",
    "ROC_COLUMNS",
]

__version__ = "0.1.0"
