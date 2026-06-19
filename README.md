# epimem

A faithful, readable Python port of the **core Moving Epidemic Method (MEM)** from the R
package [`lozalojo/mem`](https://github.com/lozalojo/mem) (José E. Lozano).

MEM learns **epidemic-onset** and **intensity** thresholds (medium / high / very high) for a
seasonal indicator — e.g. weekly influenza rates, or an 811 respiratory call share — from a set
of past complete seasons, then lets you grade the current season against them.

It also ports the **auto-tuner** and the surveillance companions, so you can pick the best
detection setting for your own seasons and run the usual reporting outputs.

| Function | What it does |
|---|---|
| `mem_model` | the epidemic-onset + intensity thresholds (the core) |
| `mem_intensity` | the four headline cut-points, labelled (`Epidemic`, `Medium (40%)`, …) |
| `mem_trend` | a real-time rising / falling signal from week-over-week changes |
| `mem_goodness` | how good a setting is, by leave-one-season-out cross-validation |
| `roc_analysis` | the **auto-tuner**: sweep the slope and pick the best value for your seasons |
| `optimum_by_inspection` | tune the slope against epidemics you marked by eye |
| `mem_stability` | how much the thresholds move as seasons accumulate |

Every one of these is verified to reproduce R `mem` on its own `flucyl` data (see Fidelity below).

> Not affiliated with or endorsed by the original author. This is an independent
> re-implementation for use where R is not available.

## Install

Needs **Python 3.11+**. From a clone of this repo:

```bash
pip install -e .            # core: numpy + scipy only
pip install -e ".[plot]"    # + matplotlib, for the surveillance chart
```

The core depends only on `numpy` and `scipy`; charting (`epimem.plot`) is an optional extra.

## Quickstart

```python
import numpy as np
from epimem import mem_model

# season_matrix: one ROW per surveillance week, one COLUMN per season.
# Three toy seasons (52 weeks, a winter bump) — use your real data instead.
weeks = np.arange(52)
season_matrix = np.column_stack([5 + amp * np.exp(-((weeks - 20) ** 2) / 30) for amp in (18, 25, 15)])

model = mem_model(season_matrix)             # -> MemThresholds (a frozen dataclass)
model.epidemic_onset                         # onset: the season has started
model.medium, model.high, model.very_high    # intensity levels
model.level_of(season_matrix[20, -1])        # grade a week -> 'baseline' / 'low (epidemic started)' / … / 'very high'
```

Tune the detection slope to your own seasons, then read the live trend signal:

```python
from epimem import roc_analysis, mem_trend

# roc_analysis (and mem_goodness) need >= 6 seasons by default; mem_model needs only 2.
# Pass min_seasons to tune on fewer, as here with 3.
tuning = roc_analysis(season_matrix, min_seasons=3)   # sweep the slope, cross-validated
best_slope = tuning.best("youden")                    # recommended value under Youden's index
model = mem_model(season_matrix, param=best_slope)    # refit with the tuned slope (overrides param=2.8)

trend = mem_trend(season_matrix)             # week-over-week rising / falling cut-points
trend.ascending, trend.descending
```

## Plotting (optional)

With the `plot` extra installed (`pip install -e ".[plot]"`), draw the standard MEM surveillance
chart — past seasons in grey, the current one bold, with the threshold lines:

```python
from epimem import mem_model
from epimem.plot import mem_chart

model = mem_model(season_matrix)
ax = mem_chart(season_matrix, model,
               season_labels=["2022-23", "2023-24", "2024-25"],
               ylabel="ILI share of 811 calls (%)")
ax.figure.savefig("mem.svg")
```

matplotlib is imported lazily, so the core package never requires it.

## Parameters (same names/defaults as R `memmodel`)

| Python (`mem_model`) | R (`memmodel`) | Default | Meaning |
|---|---|---|---|
| `method` | `i.method` | 2 | optimum method (only 2 = criterion ported) |
| `param` | `i.param` | 2.8 | slope criterion (% added per week) |
| `n_values` | `i.n.max` | `round(30/n_seasons)` | highest values per season used |
| `type_threshold` | `i.type.threshold` | 5 | onset interval: arithmetic prediction |
| `level_threshold` | `i.level.threshold` | 0.95 | onset confidence level |
| `tails_threshold` | `i.tails.threshold` | 1 | one-sided |
| `type_intensity` | `i.type.intensity` | 6 | intensity interval: geometric prediction |
| `level_intensity` | `i.level.intensity` | (0.40, 0.90, 0.975) | medium / high / very-high levels |
| `tails_intensity` | `i.tails.intensity` | 1 | one-sided |
| `use_t` | `i.use.t` | False | t vs normal quantiles |

Confidence-interval **types** match R: 1 = arithmetic mean, 2 = geometric mean, 3 = nonparametric
median (two-sided), 5 = arithmetic prediction, 6 = geometric prediction.

## Fidelity & limitations (read these)

- **Verified against R, not just claimed.** Every function above reproduces R `mem` on the
  package's own `flucyl` data to machine precision — `mem_model` thresholds (as-is *and* with
  interior gaps poked in, so the missing-week smoother is exercised), `mem_goodness`, the
  `roc_analysis` / `optimum_by_inspection` sweeps, `mem_intensity`, `mem_trend`, and the
  `mem_stability` thresholds. The checks live in `tests/test_equivalence.py`; the R-generated
  reference numbers in `tests/reference/`, each with the `generate_*.R` script that produced it.
- **Faithfully ported, not approximated:** the mean / prediction / geometric confidence intervals,
  the MAP curve, the slope-criterion optimum, per-season timing, cross-season pooling, the
  leave-one-season-out goodness metrics, the kernel smoother + its bandwidth selection
  (`sm.regression` plus `h.select`'s df=6 rule), and the rank-based tuner.
- **Deliberately simpler than R:** the median confidence interval (type 3) uses the standard
  order-statistic / sign-test interval rather than R's fiddly interpolated one. It feeds only the
  typical-duration / start-week / %-covered summary, so those few numbers don't bit-match R — every
  actual threshold still does.
- **Not ported** (intentional scope): confidence-interval type 4 (bootstrap) — it can't be made
  bit-exact against R's RNG and no default path uses it (`confidence_interval` raises a clear error
  if you ask for it); the one-sided median interval; optimum methods 1, 3, 4; and the modelled
  typical curve. The defaults don't use these.

## Reference

Vega T, Lozano JE, Ortiz de Lejarazu R, et al. *Influenza surveillance in Europe: establishing
epidemic thresholds by the moving epidemic method.* Influenza Other Respir Viruses. 2013;7(4):546-58.
Original R package: https://github.com/lozalojo/mem
