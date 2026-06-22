# epimem — a plain-English companion guide

This guide explains, with as little jargon as possible, what `epimem` is, the idea behind it,
what every piece of the code does, how you use it week to week, what it can and
cannot tell you, how we know it is correct, and where the method comes from.

It is written for a non-programmer — a public-health analyst, not a software engineer.
Wherever a technical word is unavoidable, it is explained the first time it appears, and there
is a glossary at the end.

---

## 1. The one-paragraph version

`epimem` turns several past flu seasons into a small set of **warning lines** — "the epidemic
has started," "activity is medium / high / very high" — and then lets you place each new week's
number against those lines to say, in plain words, how bad things are right now. It is a faithful
**Python** copy (Python is a popular, free programming language) of a well-established method (the
**Moving Epidemic Method**, or **MEM**) that health agencies around the world use for this purpose. It does not predict the future; it tells you
where you are today compared with a normal season.

---

## 2. The problem it solves

Every winter, respiratory illness rises and falls. A surveillance team has to answer a simple
question every week:

> "Is what we're seeing right now **normal for the time of year**, or is this an epidemic — and
> if so, how bad?"

You cannot answer that from this week's number alone. 8% of calls being about
breathing problems might be ordinary in January and alarming in August. You need a
yardstick built from history — "what does a normal season look like, and what counts as
high?" — and it needs to be standardized and defensible, so that when you tell a
decision-maker "activity is HIGH this week," that word means something specific and consistent.

That yardstick is what `epimem` builds.

---

## 3. The big idea (MEM in everyday terms)

A useful comparison is a thermometer with a fever scale printed on it.

- The thermometer is your weekly number (say, the share of 811 calls that are flu-like).
- The fever scale — "37° normal, 38° fever, 40° dangerous" — is the set of warning lines.

MEM's job is to print the right fever scale for your indicator, learned from your own
past seasons. Once the scale is printed, reading the thermometer each week is straightforward: you
see which band the needle is in.

Two points to keep in mind, because they are central to the method:

1. **The scale is drawn from past, finished seasons.** You never use the season you are
   currently in to draw its own scale, because that would be circular. Last year's seasons
   draw the lines; this year's data is the needle moving across them.
2. **The scale stays fixed for the whole season.** You do not redefine "fever" every time you
   take a temperature. You set the lines once (each autumn) and read against them all winter.

---

## 4. How it works, step by step

Here is the whole pipeline in plain language. You never have to run these steps by hand —
`mem_model` does all of it — but understanding them makes the output trustworthy rather than
opaque.

Two words used throughout: **"the model"** means the bundle of warning lines plus the steps
that drew them; nothing is being predicted, only summarized from past data. And when the method
**"learns"** something, that means it *calculates it from your history* — nothing more.

**Step 1 — Line up your seasons.** You hand the model a table: one **column per past season**,
one **row per week** of the season (aligned so that, say, week 1 is early August every year).

**Step 2 — Find the epidemic in each season.** Picture one season's weekly numbers as a **hill**:
low and flat through autumn, rising to a winter peak, then falling back. The epidemic is the steep
middle of the hill; the flat ends are the calm weeks before and after (the "shoulders"). To find
where the steep middle begins and ends, the method tries **every possible width of middle** and
measures how much of the season's total activity each width captures. A narrow middle misses a lot;
a very wide one captures almost everything but starts including the calm shoulders too. It settles
on the width right where **widening it further barely adds anything** — the natural edge of the
hill. (Internally, the table of "how much each width captures" is the *concentration table*, and the
rule for "where it stops being worth widening" is the *optimum*; see `map_curve` and
`optimum_criterion` below.)

**Step 3 — Pool the extremes across seasons.** Now the method gathers, from every season:
- the **highest off-season (pre-epidemic) values** → these set the **"epidemic has started"
  line** (the *onset threshold*). The logic: if this week is higher than the calm weeks normally
  get, the season has begun.
- the **highest in-epidemic values** → these set the **medium / high / very-high** intensity
  lines. The logic: a normal epidemic peak reaches roughly these heights, so we can grade severity.

**Step 4 — Turn each pool into a single line with a sensible margin.** A pool of "the highest
pre-season weeks from each year" is a handful of numbers, not one number. MEM summarizes each
pool as a **typical value plus a margin** (a *confidence range*) and uses the top of that range
as the line. That margin reflects having only a few seasons of history.

**Step 5 — Grade each week.** With the four lines drawn, grading is a lookup: take this week's
value, see which band it lands in, report the plain-English level.

Everything else in the package either supports these five steps, tunes one dial in
them, or checks how trustworthy the result is.

---

## 5. The tools — what each function does

Grouped by what you would actually be trying to do. A "function" here means a **named command you
give the software**: you type its name, hand it your data, and it hands back an answer. The names in
`this typeface` (like `mem_model`) are exactly what you would type. Every name here is available from
`epimem`; the optional chart tool lives in `epimem.plot`.

### 5.1 The ones you'll use almost every time

- **`mem_model`** — the headline tool. Give it your table of past seasons; it runs the whole
  five-step pipeline and hands back the warning lines. This is almost always the first thing you
  run. Defaults match the original R tool, so results line up with it.
- **`MemThresholds`** — the object `mem_model` gives back. It holds the lines
  (`.epidemic_onset`, `.medium`, `.high`, `.very_high`, plus an end-of-epidemic line
  `.post_epidemic`), some summary ranges about typical epidemics, and a built-in grader.
- **`MemThresholds.level_of(value)`** — the everyday move. Hand it this week's number; it
  returns a plain label: `'baseline'` (quiet), `'low (epidemic started)'`, `'medium'`, `'high'`,
  `'very high'`, or `'no data'`. This is what produces your weekly headline.
- **`mem_intensity(model)`** — tidies the four headline lines into a labelled table
  (`Epidemic`, `Medium (40%)`, `High (90%)`, `Very high (97.5%)`) for a report or chart legend.
  It does no new maths — it selects and labels. (The percentages are labels marking how high
  each line sits — a higher % means a higher, rarer level — not a count of anything.)
- **`mem_chart`** (in `epimem.plot`) — draws the standard surveillance chart: past seasons in
  grey, the current one bold, the threshold lines across. Optional — needs matplotlib
  (`pip install epimem[plot]`).

### 5.2 Choosing the best settings for *your* data (tuning)

MEM has one main dial — roughly, "how steep must the rise be before we call it the epidemic."
The standard value is **2.8**, and that is a reasonable default. If you want the software to pick the
best value for your own seasons instead, use these. (Two names to define up front: **"ROC"** is
a standard way of scoring a test's correct hits against its false alarms, and **"Youden"** is
one particular all-round score that balances catching real epidemics against avoiding false alarms.
You do not need the details of either.)

**Worth knowing before you reach for this:** tuning only pays off with roughly **six or more
seasons** — `roc_analysis` and `mem_goodness` need that many and refuse below it (they judge a
setting by hiding one season at a time and checking that the rest predict it, so a handful of
seasons is too little to learn from). With only three or four seasons, **stick with the default
2.8** — which is exactly what the worked example in §6 does. Come back to tuning once your history
has grown.

- **`roc_analysis`** — the automatic tuner. It tries a whole range of settings, scores each one
  by honest testing (see `mem_goodness`), and reports the winner under **eight different
  definitions of "best"** (for example, best balance of catching epidemics vs avoiding false alarms).
  Returns a `SweepResult`; call `.best("youden")` to get the single recommended value.
- **`optimum_by_inspection`** — an alternative tuner for when an analyst has **marked each past
  epidemic's start and end by eye**. It picks the setting that best reproduces your hand-marked
  calls. You must supply one (start week, end week) pair per season.
- **`mem_goodness`** — the report card behind the tuner. It grades a setting honestly using
  **"leave-one-season-out" testing**: hide one season, build the lines from the rest, check the
  hidden season week-by-week, repeat for every season, and tally how often the alerts were right.
  Returns a `Goodness` object with **sensitivity** (when there truly was an epidemic, how often
  did we catch it?), **specificity** (when there wasn't, how often did we correctly stay quiet?),
  plus a few related scores (for example, when an alert does fire, how often it is a true alert). You
  can run it on its own to ask "how trustworthy are my current settings?"
- **`SweepResult` / `CRITERIA` / `ROC_COLUMNS`** — the container the tuners return, the list of
  the eight "best" yardsticks, and the column names of the results table.

### 5.3 In-season extras

- **`mem_trend(seasons)`** — learns, from past epidemics, **how big a week-over-week jump signals
  a real climb** and how big a drop signals a real decline. Returns `TrendThresholds` with
  `.ascending` (the rise cut-off) and `.descending` (the fall cut-off). Useful for an "↑ rising /
  ↓ falling / → stable" arrow on a dashboard, on top of the level label.

### 5.4 Diagnostics — "should I trust this yet?"

- **`mem_stability(seasons)`** — rebuilds the lines using the **last 2 seasons, then the last 3,
  and so on**, so you can see whether the numbers move around (too little history) or settle
  down (enough history). Returns a `Stability` table. Use it when standing up a new indicator to
  judge how many seasons you need before the thresholds are reliable.
- **`mem_evolution(seasons)`** — the season-by-season cousin of `mem_stability`: it rebuilds the
  lines **as they would have looked each past year** — using only the seasons *before* each one
  (`method="sequential"`, a real-time back-test), or leaving each season out in turn
  (`method="cross"`). Returns an `Evolution` table. Use it to see how one big season (a severe flu
  year, say) pulls the lines around, and which lines are stable versus driven by a single year.
- **`mem_timing(one_season)`** — runs Step 2 for a **single** season: finds that season's
  epidemic window and pulls out its peak values. Returns an `EpidemicTiming` (start, end,
  duration, share-of-season, and the before/during/after peaks). Useful to ask "when did the
  epidemic run that year?"

### 5.5 Internal helpers (you rarely call these directly)

These are the building blocks the steps above use. You can call them if you are curious, but day
to day you will not need to.

- **`confidence_interval(values, ...)`** — the "low / typical / high" summarizer from Step 4. It
  comes in a few forms: a straightforward one; one tuned for numbers that pile up near zero with
  an occasional big value (common with case counts); one built around the middle value (the
  *median*) instead of the average; and wider ones that try to bracket a single future week rather
  than the typical level. MEM almost always reads the **high end** (`.upper`) to set a line. (One
  form — a "bootstrap" that relies on random resampling — is deliberately left out; see §9.)
- **`map_curve(one_season)`** — builds the **concentration table** from Step 2 ("the busiest run of
  N *consecutive* weeks held X% of the season").
- **`optimum_criterion(table, ...)`** — the **"where does the climb flatten?"** rule that picks
  the epidemic length from that table.
- **`smooth_local_linear(values, ...)`** — **smooths a weekly series** so short-term noise does
  not fool the epidemic-finder.
- **`fill_missing(values)`** — **patches over missing weeks** in the middle of a season by
  estimating what they probably were (gaps at the very start or end are left alone).
- **`max_n_values(values, n_max)`** — "give me the top few numbers from this list."
- **`ConfidenceInterval` / `EpidemicTiming` / `Goodness` / `Stability` / `TrendThresholds`** —
  small labelled containers that keep related numbers together so nothing gets mixed up.
- **`METRIC_NAMES` / `STABILITY_COLUMNS`** — fixed lists of names so result tables line up.

---

## 6. How to actually use it (the operational playbook)

There are two rhythms, and keeping them straight is the main thing.

**Once a year — draw the lines (this is the only time you "run MEM").**
In early autumn, after last season has fully finished, refit on all your complete past seasons:

```python
from epimem import mem_model, mem_trend

# season_matrix: one column per past season, one row per surveillance week
model = mem_model(season_matrix)     # -> the four warning lines, frozen for the year
trend = mem_trend(season_matrix)     # -> the rise/fall cut-offs
```

**Every week — read the needle (just a lookup, instant).**

```python
this_week = flu_calls_this_week / all_calls_this_week * 100   # your indicator
level  = model.level_of(this_week)                            # 'baseline' ... 'very high'
rising = (this_week - last_week) > trend.ascending            # is it climbing fast?
```

In plain words: work out this week's percentage; ask the model which band it falls in; and check
whether it jumped up, since last week, by more than the "clearly climbing" cut-off.

From those two lines you get the two things a bulletin needs:

1. **An onset alert** — the first week the level leaves `'baseline'`, you can announce *"the flu
   season has started."* Because an 811 call signal is available the **same day** while
   laboratory results lag one to two weeks, you can often make that call **earlier than the lab
   system** does.
2. **An intensity label** — `'high'`, `'very high'`, etc. for the weekly report, defensible
   because it is the published MEM method rather than an ad-hoc cut-off.

**Two rules that matter:**

- **Feed it a sharp, single-disease signal, not a mixed "all respiratory" bucket.** This is the
  most important practical lesson — see §8.
- **Re-draw the lines once a year.** When this season finishes, fold it into the pile that draws
  next year's lines. The season you just graded becomes part of the yardstick for the year after.

### A worked example — your 2025-26 flu season

Lines drawn from the three prior seasons (your real ILI — *influenza-like illness*, i.e. flu-like —
data): onset **2.4%**, medium **3.5%**, high **5.6%**, very high **6.9%** of calls. Watching the
2025-26 season cross them:

| When | ILI share | Level it reported |
|---|---|---|
| Aug – mid-Nov | ~0.5–1.9% | baseline (quiet) |
| **24 Nov** | 3.0% | **low — epidemic started** ← the onset alert |
| 1 Dec | 4.3% | medium |
| 8 Dec → 29 Dec | 7–13% | **very high** (peak 13.4% on 15 Dec) |
| Jan | dropping | medium → low → back to baseline |
| Feb onward | ~2% | baseline — season over |

Each new week is one new reading; whichever band it falls in is your headline, and the first
reading above the onset line is your early warning. No judgement call, no arbitrary cut-off.

**Folding 2025-26 in for next year.** Now that the 2025-26 season is finished, it joins the pile
that draws *next* year's lines (the "re-draw once a year" rule above). With all four seasons, the
2026-27 thresholds come out higher — onset **2.5%**, medium **4.1%**, high **7.7%**, very high
**10.1%** — because that severe season is now part of the yardstick. One honest caveat: with only
four seasons, the top two lines (high, very high) are essentially set by 2025-26 *alone* — drop it
and they fall by about a third — so treat them as "best estimate for now." The onset line barely
moves and is the one to trust. (`mem_evolution`, §5.4, is the tool that shows this season by season.)

---

## 7. What it is **not**: this is not a forecast

This trips up even professionals, so it is worth stating plainly: **MEM does not predict the future.**
It is a yardstick built from the past, held up against the present. It answers *"where are we
now, compared with a normal season?"* — never *"what will happen next?"*

The closest thing MEM offers to "what to expect" is a **typical-season profile** (a *climatology*):
from your seasons it can say a normal flu season **starts around late November, lasts roughly 9–10
weeks, and a normal peak reaches the high/very-high band**. That is useful for planning, but read
it for what it is — an **average of the past**, the same every year until you add new seasons. It
does not say next season will be mild or severe, or peak on a particular date.

If you genuinely want prediction — "estimate this week's not-yet-reported lab value" (a *nowcast*),
or "what will next week look like" (a *forecast*) — that is a **different kind of tool** (a
supervised model), and a separate piece of work. MEM is the "how bad is it now, has it started"
layer; forecasting is a complementary layer on top. With only a few seasons of history, any
forecast carries wide uncertainty and should be kept simple and benchmarked hard against a naive
guess.

---

## 8. Lessons from running it on real data

These came out of testing `epimem` on the actual 811 series, and they shape how you should use it.

- **Use a *peaky*, single-disease indicator — not the broad "all respiratory" bucket.** MEM was
  built for signals that **spike sharply** above their off-season floor. The flu-like (ILI) call
  share rises about **11×** in winter and gives clean, sensible thresholds. The *broad*
  respiratory share only **doubles** (because the denominator — all calls — is seasonal too,
  which flattens the bump), and on it MEM misbehaves: the "epidemic started" line can come out
  *above* the "medium" line, and the epidemic-finder grabs noise. The maths is right; the broad
  signal just is not epidemic-shaped. **Run MEM per disease (flu, COVID) on its sharp share.**
- **One clean peak per season is normal and ideal.** Flu arrives as a single dominant wave each
  winter, so one peak is what you should see — and the shape MEM is designed for.
  You would only see "a big peak plus smaller ones" if you **mix viruses together** (each peaks on
  its own schedule — COVID can peak in summer/fall, flu in midwinter), or in a rarer two-strain flu
  year. If you ever do have two genuine waves, note that MEM frames only one epidemic window per
  season — a known limit of the method, not a bug.
- **A few seasons is thin.** With only three or four seasons, the lines are usable but rough, and
  they will keep shifting as you add history. `mem_stability` is there precisely to show you how
  much they are still moving.

---

## 9. How we know it's correct

`epimem` is an independent re-implementation, so the obvious question is "does it actually match the
real method?" It does.

The original R `mem` package ships with a standard demonstration dataset (`flucyl`, eight seasons
of flu data). We ran **every user-facing tool** in both the R original and this Python port on that
data and compared the numbers (the under-the-hood helpers are covered indirectly, through those
end-to-end checks):

- The core thresholds, the cross-validated report card, both auto-tuners, the trend cut-offs, the
  intensity table, and the full stability and evolution tables all match R to **very tight tolerances** — the
  cross-validated, auto-tuner and companion outputs agree to about one part in a million, and the
  core thresholds even closer. A deliberately gap-poked version of the data (to exercise the
  missing-week patcher) agrees to within about a hundredth of a unit — still far inside rounding
  noise. The result is the same answers as R, not an approximation.

These checks are committed as automated tests (`tests/test_equivalence.py`) that run without needing
R (the statistics software the original was written in) installed, by comparing against
R-generated reference numbers stored in the repo. The committed R generator scripts (`generate_*.R`)
reproduce every reference file, so anyone can regenerate and re-verify them.

**One deliberate gap.** The random-resampling "bootstrap" confidence form (one specific style of
the low/typical/high summary, internally numbered type 4) is **not** ported. It relies on random
number generation, and R's and Python's randomness cannot be made
identical, so it could never be verified to match exactly — and no default setting uses it. Asking
for it returns a clear error rather than a silently-different number. Everything on the standard
path is exact.

---

## 10. Where the method comes from (references)

`epimem` is a Python port of the **Moving Epidemic Method (MEM)**, developed by **Tomás Vega and
José E. Lozano** and colleagues in Spain. Development began around 2001 (the Health Sentinel
Network of Castilla y León), with the first published account in 2003; the method was then
formally established and validated for Europe in 2013 — reporting, in that study, a sensitivity of
about 72% and a specificity of about 95–96% (95.5%) for detecting epidemic periods.

It is not a fringe technique. It is used by, among others:

- the **European Centre for Disease Prevention and Control (ECDC)**, in its European surveillance
  system (TESSy);
- the **World Health Organization Regional Office for Europe**, in its EuroFlu platform (since
  2012);
- the **US Centers for Disease Control and Prevention (CDC)**, which uses MEM to set its national
  flu-**severity** intensity thresholds (low / moderate / high / very high) for influenza-like
  illness, hospitalizations and deaths — though for epidemic *onset* the US relies on ILINet
  baseline percentages rather than MEM, so MEM is not its onset method;
- the **UK Health Security Agency (UKHSA)**, which has used MEM for influenza activity thresholds
  (note: as of the 2024-25 / 2025-26 seasons UKHSA has been moving several flu indicators, and
  COVID-19, to a different "mean ± standard deviation" method);
- national institutes such as **Montenegro's** (for flu thresholds) and the **Netherlands' RIVM**
  (for RSV);
- and various local / programme-level users (e.g. Los Angeles County; England's syndromic
  surveillance).

### Key references

1. **Vega T, Lozano JE, Meerhoff T, et al.** *Influenza surveillance in Europe: establishing
   epidemic thresholds by the moving epidemic method.* Influenza Other Respir Viruses.
   2013;7(4):546–558. — the foundational paper that established and validated MEM.
   https://doi.org/10.1111/j.1750-2659.2012.00422.x
2. **Vega T, Lozano JE, Meerhoff T, et al.** *Influenza surveillance in Europe: comparing intensity
   levels calculated using the moving epidemic method.* Influenza Other Respir Viruses.
   2015;9(5):234–246. — extends MEM to compare intensity levels across countries and seasons.
   https://doi.org/10.1111/irv.12330
3. **Rakocevic B, et al.** *Determining the epidemic threshold for influenza by using the Moving
   Epidemic Method (MEM), Montenegro, 2010/11 to 2017/18.* Euro Surveill. 2019;24(12):1800042. —
   a worked national application. https://doi.org/10.2807/1560-7917.ES.2019.24.12.1800042
4. **The `mem` R package**, by José E. Lozano — the reference implementation this port follows.
   On CRAN and at https://github.com/lozalojo/mem (a companion web-app, `memapp`, provides a
   point-and-click GUI).

### License and attribution

The original R `mem` package is licensed **GPL (≥ 2)**. `epimem` is released under **GPL-3.0**
(compatible with the upstream license), which means you are free to use, modify, and share it,
provided you keep it open under the same terms. `epimem` is an **independent, unaffiliated** port —
it is **not** endorsed by or connected to the original authors. If you publish thresholds produced
with it, cite the Vega/Lozano work above as the method, and note that the figures were produced by
an independent Python implementation verified against the R original.

---

## 11. Glossary

- **Indicator** — the weekly number you track (e.g. the share of 811 calls that are flu-like).
- **ILI (influenza-like illness)** — the flu-like syndrome; here, the share of calls flagged as
  flu-like.
- **Season** — one respiratory year, here defined August → July, so the winter peak sits in the
  middle.
- **Shoulder** — the calm, flat weeks on either side of the winter peak (the ends of the "hill").
- **Function** — a named command you give the software: you type its name, hand it your data, and
  it hands back an answer (e.g. `mem_model`).
- **R / Python / port** — R is the statistics software the original `mem` was written in; Python is
  the (different) programming language this version is written in; a *port* is the same method
  rebuilt in another language.
- **Threshold / line** — a learned value you compare each week against.
- **Epidemic (onset) threshold** — the line that, once crossed, says "the season has started."
- **Intensity thresholds** — the medium / high / very-high lines that grade how bad it is.
- **Confidence range** — a "low / typical / high" summary of a handful of numbers that is honest
  about small samples; MEM uses the high end as a line.
- **Sensitivity** — of the weeks that truly were epidemic, the fraction the thresholds flagged.
- **Specificity** — of the weeks that truly were quiet, the fraction the thresholds left unflagged.
- **Cross-validation (leave-one-season-out)** — an honest test: hide one season, build from the
  rest, check the hidden one; repeat for each season.
- **Climatology** — the "typical season" averaged from history; a planning baseline, **not** a
  forecast.
- **Nowcast / forecast** — estimating the present-but-not-yet-reported value / predicting a future
  value; both are separate tools, not part of MEM.
- **Peaky** — a signal with a sharp, tall winter spike (good for MEM) versus a broad, gentle bump
  (poor for MEM).
