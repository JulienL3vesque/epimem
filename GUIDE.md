# epimem guide

What `epimem` is, how to use it, its limits, and how it was verified against the R original.
Glossary at the end.

---

## 1. In one paragraph

`epimem` turns past flu seasons into a small set of warning lines ("the epidemic has started",
"activity is medium, high, or very high") and grades each new week against them. It is a Python port
of the Moving Epidemic Method (MEM), used by health agencies across Europe and beyond. It does not
predict the future. It places this week against a normal season.

---

## 2. The problem it solves

Respiratory illness rises and falls every winter, and a surveillance team has to answer the same
question each week: is this normal for the time of year, or an epidemic, and if so how bad?

This week's number alone cannot answer it. A reading of 8% might be ordinary in January and alarming
in August. A yardstick built from past seasons is needed, one that is standard and defensible, so
that "activity is HIGH this week" means something specific. That is what `epimem` builds.

---

## 3. The big idea

MEM compares the weekly number (say, weekly influenza rates, or the share of visits that are
flu-like) against warning lines learned from past seasons. The number can be anything tracked weekly,
a rate per 100,000, a percent share, or a raw count; the thresholds come back in the same unit. Below
the first line is normal; above it the season has started; the higher lines grade intensity.

Two rules define the method:

1. The lines come from past, finished seasons, never the current one (that would be circular). Last
   year's seasons draw the lines; this year's data is read against them.
2. The lines are fixed for the season. They are set each autumn and read against all winter.

---

## 4. How it works, step by step

`mem_model` runs all five steps; none need to be done by hand. Knowing them makes the output readable
rather than opaque. "The model" is the warning lines plus the steps that drew them.

1. Line up the seasons. A table with one column per past season and one row per week, aligned so week
   1 is, say, early August every year.
2. Find each season's epidemic. A season runs low and flat in autumn, rises to a winter peak, and
   falls back. The epidemic is the steep middle; the flat ends are the shoulders. MEM tries every
   width of middle, measures how much of the season's activity each captures, and stops where
   widening it further adds almost nothing. (See `map_curve` and `optimum_criterion`.)
3. Pool the extremes across seasons. The highest pre-epidemic weeks from every season set the onset
   line: higher than the calm weeks normally get means the season has begun. The highest in-epidemic
   weeks set the medium, high, and very-high lines.
4. Draw each line with a margin. Each pool is a handful of numbers, not one. MEM takes a typical
   value plus a margin (a confidence range) and uses the top of the range. The margin reflects having
   few seasons.
5. Grade each week. The band this week's value lands in is the reported level.

Everything else in the package supports these steps, tunes one of them, or checks the result.

---

## 5. The functions

The names in `this typeface` are the importable names. All come from `epimem`; the chart is in
`epimem.plot`.

### 5.1 Everyday

- `mem_model`. Given a table of past seasons, returns the warning lines. Defaults match the R
  original.
- `MemThresholds`. What `mem_model` returns. Holds the lines (`.epidemic_onset`, `.medium`, `.high`,
  `.very_high`, and `.post_epidemic`), summary ranges, and a grader.
- `MemThresholds.level_of(value)`. Takes this week's number, returns a label: `'baseline'`,
  `'low (epidemic started)'`, `'medium'`, `'high'`, `'very high'`, or `'no data'`.
- `mem_intensity(model)`. The four lines as a labelled table (`Epidemic`, `Medium (40%)`,
  `High (90%)`, `Very high (97.5%)`) for a report or legend. The percentages mark how high each line
  sits, not a count.
- `mem_chart(...)` (in `epimem.plot`). The standard surveillance chart: past seasons grey, current one
  bold, threshold lines across. Needs matplotlib (`pip install -e ".[plot]"`).

### 5.2 Tuning

MEM has one dial, `param` (default 2.8): how steep the week-over-week rise must be before it counts
as the epidemic. The functions below let the software pick a value from the seasons instead. The
cross-validated tuners (`roc_analysis`, `mem_goodness`) need six or more seasons and refuse fewer, so
with three or four seasons the default stands; `optimum_by_inspection` works on fewer.

- `roc_analysis`. Tries a range of settings, scores each by cross-validation, and reports the winner
  under eight definitions of "best". Returns a `SweepResult`; `.best("youden")` gives the recommended
  value. (ROC and Youden are standard ways to score hits against false alarms.)
- `optimum_by_inspection`. For when an analyst has marked each past epidemic's start and end by eye.
  Picks the setting that best reproduces those marks. Needs one (start, end) pair per season.
- `mem_goodness`. Scores a setting by leave-one-season-out cross-validation: hide a season, build the
  lines from the rest, check the hidden one, repeat. Returns sensitivity (epidemics caught),
  specificity (quiet weeks left alone), and related scores. Can be run alone to check the current
  settings.
- `SweepResult`, `CRITERIA`, `ROC_COLUMNS`. The tuner's result container, the eight scoring criteria,
  and the result-table column names.

### 5.3 In-season

- `mem_trend(seasons)`. Learns how big a week-over-week jump means a real climb and how big a drop
  means a real decline. Returns `.ascending` and `.descending` cut-offs, in the same unit as the
  indicator: a week is rising fast when `this_week - last_week` exceeds `.ascending`.

### 5.4 Diagnostics

- `mem_stability(seasons)`. Rebuilds the lines on the last 2 seasons, then 3, and so on, showing
  whether they have settled or are still moving. Useful for judging how many seasons are needed.
- `mem_evolution(seasons)`. Rebuilds the lines as they would have looked each past year: from only the
  seasons before each one (`method="sequential"`, a real-time back-test), or leaving each season out
  (`method="cross"`). Shows how one big season pulls the lines around.
- `mem_timing(one_season)`. Step 2 for a single season: its epidemic window and peak values. Returns
  start, end, duration, share of season, and the before, during, and after peaks.

### 5.5 Internals (rarely called directly)

- `confidence_interval(values, ...)`. The low/typical/high summary from step 4. Several forms: plain;
  one for skewed counts; one around the median; and wider ones that bracket a single future week. MEM
  reads the high end (`.upper`). The bootstrap form is left out (see section 9).
- `map_curve(one_season)`. The table behind step 2: the busiest run of N consecutive weeks held X% of
  the season.
- `optimum_criterion(table, ...)`. Picks the epidemic length from that table.
- `smooth_local_linear(values, ...)`. Smooths a weekly series so noise does not fool the
  epidemic-finder.
- `fill_missing(values)`. Fills missing weeks inside a season; leaves gaps at the start or end alone.
- `max_n_values(values, n_max)`. The top n numbers from a list.
- `ConfidenceInterval`, `EpidemicTiming`, `Goodness`, `Stability`, `TrendThresholds`. Result
  containers.
- `METRIC_NAMES`, `STABILITY_COLUMNS`. Fixed name lists so tables line up.

---

## 6. How to use it

Two rhythms.

Once a year, draw the lines. In early autumn, after the last season has finished, refit on every
complete past season:

```python
from epimem import mem_model, mem_trend

# season_matrix: one column per past season, one row per surveillance week
model = mem_model(season_matrix)     # the four warning lines, fixed for the year
trend = mem_trend(season_matrix)     # the rise/fall cut-offs
```

Every week, read the band:

```python
this_week = current_week_value                     # the indicator this week (a rate, a %, ...)
level  = model.level_of(this_week)                 # 'baseline' ... 'very high'
rising = (this_week - last_week) > trend.ascending # climbing faster than usual?
```

A bulletin then has two parts:

1. Onset. The first week the level leaves `'baseline'`, the season has started. A syndromic signal
   (call lines, visits, prescriptions) is often available the same day, while lab confirmation lags
   one to two weeks, so this call can come earlier than the lab system's.
2. Intensity. `'high'`, `'very high'`, and so on, defensible because it is the published MEM method.

Two rules: the signal should be a sharp single-disease one, not an "all respiratory" mix (section 8);
and the lines are re-drawn each year, folding the season just graded into the next year's set.

### A worked example: the bundled `flucyl` data

`epimem` ships with R mem's `flucyl` dataset (eight past influenza seasons, weekly rates per
100,000). Learning the thresholds from all eight gives onset 53, medium 250, high 442, very high
568. Grading a few example weeks against them:

| this week's rate | level |
|---|---|
| 30 | baseline |
| 120 | low (epidemic started) |
| 300 | medium |
| 480 | high |
| 650 | very high |

Each week is one reading, looked up against the lines; the first reading above onset (53) is the
onset alert. The whole thing is reproducible: `example_seasons()` returns the data and
`mem_model(...).level_of(value)` does the grading. With only a few seasons the top lines can still
move as history grows; `mem_stability` and `mem_evolution` (5.4) show how settled they are.

---

## 7. It is not a forecast

MEM measures the present against the past: where things stand now versus a normal season.

The nearest thing to a forecast it offers is a typical-season profile (climatology): a normal flu
season starts around late November, lasts about 9 to 10 weeks, and peaks in the high or very-high
band. That is a planning baseline, the same every year until more seasons are added. It does not say
whether next season will be mild or severe.

Real prediction (a nowcast, estimating this week's not-yet-reported lab value, or a forecast of next
week) is a separate, supervised model. MEM is the "how bad is it now" layer; forecasting sits on top.
On a few seasons of history, any forecast should be kept simple and benchmarked hard against a naive
guess.

---

## 8. Lessons from real data

A few practical lessons:

- A sharp, single-disease indicator works best, not a broad mixed bucket. MEM needs a signal that
  spikes well above its off-season floor. A flu-specific signal does that and gives clean thresholds.
  A broad "all respiratory" signal barely rises (its denominator is seasonal too), and MEM misbehaves
  on it: the onset line can land above the medium line. Best run per disease.
- One peak per season is the norm, and what flu gives. Two peaks usually means viruses mixed together
  (COVID and flu peak at different times), or a rarer two-strain year. MEM frames one epidemic window
  per season, a known limit.
- Few seasons means rough lines. On three or four seasons they are usable but keep moving as history
  grows; `mem_stability` shows how much.

---

## 9. How it was verified

The R `mem` package ships with a demonstration dataset (`flucyl`, eight seasons). Every user-facing
function was run in both the R original and this port on that data and the numbers compared:

- The thresholds, the cross-validated scores, both tuners, the trend cut-offs, the intensity table,
  and the stability and evolution tables match R to about one part in a million, the core thresholds
  closer. A version with gaps poked into it (to exercise the missing-week fill) matches to about one
  part in a million.

The R comparison is the evidence. The port was written by reading the open R source and checking the
output against real R numbers until they matched, and those checks rerun on every change.

The checks are automated tests (`tests/test_equivalence.py`) that run without R installed, against
R-generated reference numbers in the repo. The `generate_*.R` scripts reproduce every reference file.

The one exception is the bootstrap confidence form (internally type 4): it relies on random
resampling, and R's and Python's randomness cannot be made identical, so it could not be verified
exactly. No default uses it, and asking for it raises an error.

---

## 10. References

MEM was developed by Tomás Vega, José E. Lozano, and colleagues in Spain, starting around 2001 (the
Health Sentinel Network of Castilla y León) and validated for Europe in 2013 (sensitivity about 72%,
specificity about 95.5% for detecting epidemic periods).

It is widely used: the ECDC (in TESSy), the WHO Regional Office for Europe (EuroFlu, since 2012), the
US CDC (for national flu-severity intensity thresholds, though not for onset, which the US sets from
ILINet baselines), the UK Health Security Agency (for influenza thresholds, though it has been moving
several indicators to a mean-and-standard-deviation method for 2024-25 onward), national institutes
such as Montenegro's and the Netherlands' RIVM (for RSV), and various local programmes.

1. Vega T, Lozano JE, Meerhoff T, et al. *Influenza surveillance in Europe: establishing epidemic
   thresholds by the moving epidemic method.* Influenza Other Respir Viruses. 2013;7(4):546-558.
   https://doi.org/10.1111/j.1750-2659.2012.00422.x
2. Vega T, Lozano JE, Meerhoff T, et al. *Influenza surveillance in Europe: comparing intensity levels
   calculated using the moving epidemic method.* Influenza Other Respir Viruses. 2015;9(5):234-246.
   https://doi.org/10.1111/irv.12330
3. Rakocevic B, et al. *Determining the epidemic threshold for influenza by using the Moving Epidemic
   Method (MEM), Montenegro, 2010/11 to 2017/18.* Euro Surveill. 2019;24(12):1800042.
   https://doi.org/10.2807/1560-7917.ES.2019.24.12.1800042
4. The `mem` R package by José E. Lozano, the reference implementation this port follows. On CRAN and
   at https://github.com/lozalojo/mem.

The R `mem` package is licensed GPL (2 or later); `epimem` is GPL-3.0, an independent port not
endorsed by or connected to the original authors. Published thresholds should cite the Vega/Lozano
work as the method, noting that the figures came from an independent Python implementation verified
against the R original.

---

## 11. Glossary

- Indicator. The weekly number being tracked (e.g. weekly influenza rates, or the share of visits
  that are flu-like).
- ILI (influenza-like illness). The flu-like syndrome.
- Season. One respiratory year, here August to July, so the winter peak sits in the middle.
- Shoulder. The calm weeks on either side of the winter peak.
- Threshold / line. A learned value each week is compared against.
- Onset threshold. The line that, once crossed, marks the start of the season.
- Intensity thresholds. The medium, high, and very-high lines.
- Confidence range. A low/typical/high summary of a few numbers; MEM uses the high end.
- Sensitivity. Of the truly-epidemic weeks, the fraction flagged.
- Specificity. Of the truly-quiet weeks, the fraction left unflagged.
- Cross-validation (leave-one-season-out). Hide a season, build from the rest, check the hidden one,
  repeat.
- Climatology. The typical season averaged from history; a planning baseline, not a forecast.
- Nowcast / forecast. Estimating the not-yet-reported present, or predicting the future; separate
  tools.
- Peaky. A signal with a sharp winter spike (good for MEM), not a broad gentle bump.
