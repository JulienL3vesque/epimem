"""Tests for epimem. Run: python tests/test_epimem.py  (or: pytest -q)

Component checks verify each ported piece against the R formula. The flucyl test is the
end-to-end equivalence check vs R; it's skipped until someone runs R once (see its note).
"""
import numpy as np
from scipy.stats import norm

import epimem


def test_prediction_interval():
    # Type 5: mean +/- z*sd. Check lower / centre / upper, not just upper.
    v = np.array([10.0, 12.0, 14.0])
    ci = epimem.confidence_interval(v, 0.95, confidence_interval_type=5, tails=1)
    margin = norm.ppf(0.95) * v.std(ddof=1)
    assert np.isclose(ci.upper, 12 + margin)
    assert np.isclose(ci.centre, 12.0)
    assert np.isclose(ci.lower, 12 - margin)


def test_mean_interval():
    # Type 1: mean + z*sd/sqrt(n).
    v = np.array([10.0, 12.0, 14.0])
    got = epimem.confidence_interval(v, 0.95, confidence_interval_type=1, tails=1).upper
    assert np.isclose(got, 12 + norm.ppf(0.95) * v.std(ddof=1) / np.sqrt(3))


def test_geometric_interval():
    # Type 6: exp of the type-5 interval computed on the logs.
    v = np.array([10.0, 12.0, 14.0])
    expected = np.exp(np.log(v).mean() + norm.ppf(0.90) * np.log(v).std(ddof=1))
    assert np.isclose(epimem.confidence_interval(v, 0.90, confidence_interval_type=6, tails=1).upper, expected)


def test_unsupported_interval_type_raises():
    try:
        epimem.confidence_interval([1.0, 2.0], confidence_interval_type=4)   # bootstrap isn't ported
    except ValueError:
        return
    raise AssertionError("expected ValueError for confidence_interval_type=4")


def test_example_seasons():
    seasons = epimem.example_seasons()
    assert seasons.shape == (33, 8)              # the bundled flucyl data
    assert np.isfinite(seasons).all()
    model = epimem.mem_model(seasons)            # must feed straight into the model
    assert model.epidemic_onset < model.high


def test_max_n_values():
    assert list(epimem.max_n_values([1, 5, 3, np.nan, 2], 3)[:3]) == [5.0, 3.0, 2.0]


def test_map_curve():
    # tiny season, total = 18
    table = epimem.map_curve(np.array([0.0, 0, 5, 10, 3, 0]))
    assert np.isclose(table[1, 1], 100 * 10 / 18) and int(table[1, 3]) == 4    # best 1 week
    assert np.isclose(table[2, 1], 100 * 15 / 18)                              # best 2 weeks
    assert int(table[2, 3]) == 3 and int(table[2, 4]) == 4


def test_mem_model():
    rng = np.random.default_rng(3)

    def season(amp, peak=22):
        w = np.arange(52)
        return np.clip(5 + rng.normal(0, 0.5, 52) + amp * np.exp(-((w - peak) ** 2) / (2 * 5 ** 2)), 0, None)

    th = epimem.mem_model(np.column_stack([season(20), season(28), season(18)]))
    assert th.medium < th.high < th.very_high
    assert th.epidemic_onset < th.high
    assert th.n_values_per_season == 10          # round(30 / 3)
    assert th.level_of(th.very_high + 1) == "very high"
    assert th.level_of(0) == "baseline"


def test_mem_evolution():
    rng = np.random.default_rng(7)

    def season(amp):
        w = np.arange(40)
        return np.clip(5 + rng.normal(0, 0.5, 40) + amp * np.exp(-((w - 18) ** 2) / (2 * 5 ** 2)), 0, None)

    seasons = np.column_stack([season(a) for a in (16, 22, 14, 25, 19, 28)])   # 40 weeks x 6 seasons
    n = seasons.shape[1]
    column = {name: i for i, name in enumerate(epimem.STABILITY_COLUMNS)}

    # sequential: a row per season from the 3rd on, plus "next". cross: a row per season, plus "next".
    seq = epimem.mem_evolution(seasons, method="sequential")
    crs = epimem.mem_evolution(seasons, method="cross")
    assert seq.data.shape == (n - 1, len(epimem.STABILITY_COLUMNS))
    assert crs.data.shape == (n + 1, len(epimem.STABILITY_COLUMNS))

    # The final "next" row reproduces the all-seasons model on every threshold.
    full = epimem.mem_model(seasons)
    for name, attr in [("epidemic", "epidemic_onset"), ("postepidemic", "post_epidemic"),
                       ("medium", "medium"), ("high", "high"), ("veryhigh", "very_high")]:
        assert np.isclose(seq.data[-1, column[name]], getattr(full, attr))
        assert np.isclose(crs.data[-1, column[name]], getattr(full, attr))

    # evolution_seasons caps how many seasons feed each fit: uncapped leave-one-out uses every
    # other season, capped uses only the nearest few.
    assert crs.seasons_used[0].sum() == n - 1
    capped = epimem.mem_evolution(seasons, method="cross", evolution_seasons=3)
    assert capped.seasons_used[0].sum() == 3

    # an unknown method is rejected.
    try:
        epimem.mem_evolution(seasons, method="bogus")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for an unknown method")


def test_edge_cases():
    # One value: the interval centre is the value, the bounds are NaN (matches R, no warning).
    interval = epimem.confidence_interval([7.0], confidence_interval_type=5)
    assert interval.centre == 7 and np.isnan(interval.upper) and np.isnan(interval.lower)
    # fewer than two seasons -> ValueError
    try:
        epimem.mem_model(np.array([[1.0], [2.0], [3.0]]))   # one season
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for <2 seasons")
    # all-zero season: percentages are all 0, nothing blows up
    table = epimem.map_curve(np.zeros(4))
    assert np.all(table[:, 1] == 0) and np.all(np.isfinite(table))


def test_median_interval():
    # Type 3 is our own order-statistic (sign-test) median interval, not R's interpolated version.
    seven = epimem.confidence_interval([2, 2, 5, 9, 9, 12, 15], confidence_interval_type=3, tails=2)
    assert (seven.lower, seven.centre, seven.upper) == (2.0, 9.0, 15.0)
    five = epimem.confidence_interval([1, 2, 3, 4, 5], confidence_interval_type=3, tails=2)
    assert (five.lower, five.centre, five.upper) == (1.0, 3.0, 5.0)
    # Larger n: the order-statistic ranks are no longer the plain min/max, so these exercise the
    # binomial rank math that the small-n cases (where the ranks collapse to 1 and n) cannot.
    twelve = epimem.confidence_interval(range(1, 13), confidence_interval_type=3, tails=2)
    assert (twelve.lower, twelve.centre, twelve.upper) == (3.0, 6.5, 10.0)
    fifteen = epimem.confidence_interval(range(1, 16), confidence_interval_type=3, tails=2)
    assert (fifteen.lower, fifteen.centre, fifteen.upper) == (4.0, 8.0, 12.0)
    # The one-sided median interval is intentionally not supported.
    try:
        epimem.confidence_interval([1, 2, 3, 4, 5], confidence_interval_type=3, tails=1)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("expected NotImplementedError for a one-sided median interval")


def test_confidence_branches():
    # Branches the equivalence tests never reach: use_t, the geometric zero-shift, and two-sided tails.
    # use_t=True prediction uses the t quantile * sqrt(1 + 1/n), not the normal quantile.
    use_t = epimem.confidence_interval([2.0, 4, 6], confidence_interval_type=5, tails=1, use_t=True).upper
    assert np.isclose(use_t, 10.7434178)
    # A geometric interval containing a zero exercises the log(x + 1) shift and its undo.
    geo = epimem.confidence_interval([0.0, 2, 4, 6], confidence_interval_type=6, tails=1)
    assert np.isclose(geo.centre, 2.2010859) and np.isclose(geo.upper, 11.9628730)
    # Two-sided uses the wider 0.975 quantile; one-sided the 0.95 one.
    one = epimem.confidence_interval([2.0, 4, 6], confidence_interval_type=1, tails=1).upper
    two = epimem.confidence_interval([2.0, 4, 6], confidence_interval_type=1, tails=2).upper
    assert np.isclose(one, 4 + norm.ppf(0.95) * 2 / np.sqrt(3))
    assert np.isclose(two, 4 + norm.ppf(0.975) * 2 / np.sqrt(3))


def test_fill_missing_boundary():
    # Only interior gaps are reconstructed; leading and trailing NaNs are left as NaN.
    weeks = np.arange(20)
    series = 5 + 10 * np.exp(-((weeks - 10) ** 2) / 8)   # a smooth bump, enough points to smooth
    series[0:2] = np.nan        # leading gap
    series[18:20] = np.nan      # trailing gap
    series[10] = np.nan         # interior gap, at the peak
    filled = epimem.fill_missing(series)
    assert np.isnan(filled[0]) and np.isnan(filled[1])        # leading left untouched
    assert np.isnan(filled[18]) and np.isnan(filled[19])      # trailing left untouched
    assert np.isfinite(filled[10])                            # interior gap filled


def test_confusion_matrix():
    # Addition combines the counts field by field, and sum() starts from an empty matrix.
    a = epimem.ConfusionMatrix(true_positives=3, false_positives=1, true_negatives=5, false_negatives=2)
    b = epimem.ConfusionMatrix(1, 1, 1, 1)
    total = a + b
    assert (total.true_positives, total.false_positives, total.true_negatives, total.false_negatives) == (4, 2, 6, 3)
    assert sum([a, b], epimem.ConfusionMatrix()) == total

    # The diagnostic scores match the textbook definitions on a hand-checked 2x2.
    scores = epimem.ConfusionMatrix(true_positives=8, false_positives=2, true_negatives=80, false_negatives=10).scores()
    assert np.isclose(scores.sensitivity, 8 / 18)
    assert np.isclose(scores.specificity, 80 / 82)
    assert np.isclose(scores.percent_agreement, 0.88)
    assert np.isclose(scores.youden_index, 8 / 18 + 80 / 82 - 1)

    # A zero denominator is reported as NaN rather than raising.
    assert np.isnan(epimem.ConfusionMatrix().scores().sensitivity)


def test_mem_chart():
    # The optional epimem.plot module renders the MEM chart without error. Skipped (-> [SKIP]) when
    # matplotlib is not installed, since it is only an extra (epimem[plot]).
    try:
        import matplotlib
    except ImportError:
        raise NotImplementedError("matplotlib not installed (epimem[plot])")
    matplotlib.use("Agg")                              # headless: render to a buffer, no window
    import matplotlib.pyplot as plt

    from epimem.plot import mem_chart

    weeks = np.arange(30)
    seasons = np.column_stack([
        5 + amplitude * np.exp(-((weeks - peak) ** 2) / 20.0)
        for amplitude, peak in [(10, 15), (12, 16), (9, 14)]
    ])
    ax = mem_chart(seasons, epimem.mem_model(seasons), season_labels=["s1", "s2", "s3"], ylabel="%")
    # three season curves + four threshold reference lines.
    assert len(ax.get_lines()) == 3 + 4
    plt.close("all")


# The full mem-equivalence proof (epimem vs R on flucyl, clean and gapped) lives in
# tests/test_equivalence.py - it checks committed R-generated reference numbers.


if __name__ == "__main__":
    import sys

    failures = []
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"[OK]   {name}")
            except NotImplementedError:
                print(f"[SKIP] {name}")
            except Exception as exc:
                failures.append(name)
                print(f"[FAIL] {name}: {exc}")
    print(f"\n{len(failures)} failed" if failures else "\nall passed")
    sys.exit(1 if failures else 0)
