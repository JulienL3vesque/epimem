"""Tests for epimem. Run: python tests/test_epimem.py  (or: pytest -q)

Component checks verify each ported piece against the R formula. The flucyl test is the
end-to-end equivalence check vs R; it's skipped until someone runs R once (see its note).
"""
import numpy as np
from scipy.stats import norm

import epimem


def test_prediction_interval():
    # Type 5: mean + z*sd.
    v = np.array([10.0, 12.0, 14.0])
    got = epimem.confidence_interval(v, 0.95, confidence_interval_type=5, tails=1).upper
    assert np.isclose(got, 12 + norm.ppf(0.95) * v.std(ddof=1))


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
    # The one-sided median interval is intentionally not supported.
    try:
        epimem.confidence_interval([1, 2, 3, 4, 5], confidence_interval_type=3, tails=1)
    except NotImplementedError:
        pass
    else:
        raise AssertionError("expected NotImplementedError for a one-sided median interval")


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
