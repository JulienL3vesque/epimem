"""Equivalence tests: epimem vs R `mem`, checked against committed R-generated numbers.

The reference CSVs in reference/ were produced by R `mem` itself (see the generate_*.R scripts
there). These tests need no R - they check epimem against those committed numbers. To refresh
them, re-run the generators in tests/reference/. Files are read with the standard-library `csv`
module + numpy, so the test suite needs nothing beyond the package's own numpy / scipy.
"""
import csv
from pathlib import Path

import numpy as np

import epimem

REFERENCE = Path(__file__).parent / "reference"
THRESHOLD_NAMES = ["epidemic_onset", "post_epidemic", "medium", "high", "very_high"]


def _num(cell):
    """A CSV cell as float; R's '' / 'NA' become NaN, while 'Inf' / '-Inf' / 'NaN' parse directly."""
    cell = cell.strip()
    return np.nan if cell in ("", "NA") else float(cell)


def _matrix(name):
    """A numeric reference table (one header row, then numbers) as a float array."""
    with open(REFERENCE / name, newline="") as f:
        body = list(csv.reader(f))[1:]
    return np.array([[_num(cell) for cell in row] for row in body], dtype=float)


def _rows(name):
    """A labelled reference table as a list of {column: value} dicts."""
    with open(REFERENCE / name, newline="") as f:
        return list(csv.DictReader(f))


def _flucyl():
    return _matrix("flucyl.csv")


# memmodel: the epidemic thresholds every other layer is built on.

def _expected_thresholds():
    table = {}
    for row in _rows("reference_thresholds.csv"):
        table.setdefault(row["case"], {})[row["name"]] = float(row["value"])
    return table


def _thresholds(csv_name):
    th = epimem.mem_model(_matrix(csv_name))
    return {name: getattr(th, name) for name in THRESHOLD_NAMES}


def test_thresholds_clean():
    expected, got = _expected_thresholds()["clean"], _thresholds("flucyl.csv")
    for name in THRESHOLD_NAMES:
        assert abs(got[name] - expected[name]) <= 1e-4, f"clean/{name}: {got[name]} vs {expected[name]}"


def test_thresholds_gappy():
    # The gapped series is reconstructed by fill_missing before thresholds are built, and R and
    # numpy round that reconstruction slightly differently, so allow a looser absolute tolerance
    # (1e-2) than the clean case (1e-4). Still far tighter than any threshold difference that
    # would change a public-health call.
    expected, got = _expected_thresholds()["gappy"], _thresholds("flucyl_gappy.csv")
    for name in THRESHOLD_NAMES:
        assert abs(got[name] - expected[name]) <= 1e-2, f"gappy/{name}: {got[name]} vs {expected[name]}"


# memgoodness: leave-one-season-out diagnostic scores.

def _goodness_expected(case):
    return {row["name"]: float(row["value"])
            for row in _rows("goodness_reference.csv") if row["grid"] == case}


def test_goodness_equivalence():
    flucyl = _flucyl()
    for case, detection in {"default": None, "coarse": 1.0 + 0.5 * np.arange(9)}.items():
        g = epimem.mem_goodness(flucyl, detection_values=detection)
        expected = _goodness_expected(case)
        for name in epimem.METRIC_NAMES:
            got = getattr(g, name)
            assert abs(got - expected[name]) <= 1e-6, f"goodness {case}/{name}: {got} vs {expected[name]}"


# roc.analysis: pick the slope by cross-validation.

_ROC_RENAME = {"pos.likehood": "pos_likelihood", "neg.likehood": "neg_likelihood", "aditive": "additive"}


def test_roc_equivalence():
    roc = epimem.roc_analysis(_flucyl(), param_values=np.array([2.0, 2.5, 3.0]))
    # equal_nan handles the Inf/NaN cells (e.g. +LR = sens/0 when specificity hits 1)
    assert np.isclose(roc.table, _matrix("roc_reference.csv"), rtol=1e-6, atol=1e-9, equal_nan=True).all()
    for row in _rows("roc_optimum.csv"):
        key = _ROC_RENAME.get(row["criterion"], row["criterion"])
        assert roc.optimum[key] == float(row["value"]), \
            f"roc optimum {key}: {roc.optimum[key]} vs {row['value']}"


# optimum.by.inspection: pick the slope against analyst-marked epidemics.

def test_inspection_equivalence():
    timings = [(int(float(row["start"])), int(float(row["end"])))
               for row in _rows("inspection_timings.csv")]
    insp = epimem.optimum_by_inspection(_flucyl(), timings, param_values=2.0 + 0.1 * np.arange(11))

    # equal_nan handles the Inf/NaN cells (e.g. +LR = sens/0 when specificity hits 1)
    assert np.isclose(insp.table, _matrix("inspection_reference.csv"),
                      rtol=1e-6, atol=1e-9, equal_nan=True).all()
    for row in _rows("inspection_optimum.csv"):
        assert insp.optimum[row["criterion"]] == float(row["value"]), \
            f"inspection optimum {row['criterion']}: {insp.optimum[row['criterion']]} vs {row['value']}"


# memintensity, memtrend, memstability: the surveillance companions.

def test_intensity_equivalence():
    cut_points = epimem.mem_intensity(epimem.mem_model(_flucyl()))
    for row in _rows("intensity_reference.csv"):
        assert abs(cut_points[row["label"]] - float(row["value"])) <= 1e-4, row["label"]


def test_trend_equivalence():
    ref = {row["name"]: float(row["value"]) for row in _rows("trend_reference.csv")}
    tr = epimem.mem_trend(_flucyl())
    assert abs(tr.ascending - ref["ascending"]) <= 1e-6
    assert abs(tr.descending - ref["descending"]) <= 1e-6


def test_stability_equivalence():
    st = epimem.mem_stability(_flucyl())
    ref = _rows("stability_reference.csv")
    assert list(st.counts) == [int(float(row["count"])) for row in ref]

    column = {name: index for index, name in enumerate(epimem.STABILITY_COLUMNS)}

    # The threshold columns use mem's mean / geometric CIs and still match R exactly.
    threshold_columns = ["epidemic", "postepidemic", "medium", "high", "veryhigh"]
    got_values = st.data[:, [column[name] for name in threshold_columns]]
    ref_values = np.array([[float(row[name]) for name in threshold_columns] for row in ref])
    assert np.max(np.abs(got_values - ref_values)) <= 1e-6

    # The duration / start / %-covered columns use our own simpler median confidence interval
    # rather than R's interpolated one, so we only check they are internally sensible:
    # lower <= centre <= upper.
    for lower, centre, upper in [("durationll", "duration", "durationul"),
                                 ("startll", "start", "startul"),
                                 ("percentagell", "percentage", "percentageul")]:
        assert (st.data[:, column[lower]] <= st.data[:, column[centre]] + 1e-9).all()
        assert (st.data[:, column[centre]] <= st.data[:, column[upper]] + 1e-9).all()


# memevolution: how the thresholds would have evolved season by season (both validation methods).

def test_evolution_equivalence():
    flucyl = _flucyl()
    ref = _rows("evolution_reference.csv")
    column = {name: index for index, name in enumerate(epimem.STABILITY_COLUMNS)}
    threshold_columns = ["epidemic", "postepidemic", "medium", "high", "veryhigh"]

    for method in ("sequential", "cross"):
        ev = epimem.mem_evolution(flucyl, method=method)
        rows = sorted((row for row in ref if row["method"] == method), key=lambda r: int(r["row"]))
        assert ev.data.shape[0] == len(rows), f"{method}: {ev.data.shape[0]} rows vs {len(rows)}"

        # The "number" column (seasons used in each fit) matches R exactly.
        assert list(ev.counts) == [int(float(row["number"])) for row in rows], f"{method} counts"

        # The threshold columns use mem's mean / geometric CIs and still match R exactly.
        got_values = ev.data[:, [column[name] for name in threshold_columns]]
        ref_values = np.array([[float(row[name]) for name in threshold_columns] for row in rows])
        assert np.max(np.abs(got_values - ref_values)) <= 1e-6, f"{method} thresholds"

        # The duration / start / %-covered columns use our own simpler median CI rather than R's
        # interpolated one, so we only check they are internally sensible: lower <= centre <= upper.
        for lower, centre, upper in [("durationll", "duration", "durationul"),
                                     ("startll", "start", "startul"),
                                     ("percentagell", "percentage", "percentageul")]:
            assert (ev.data[:, column[lower]] <= ev.data[:, column[centre]] + 1e-9).all()
            assert (ev.data[:, column[centre]] <= ev.data[:, column[upper]] + 1e-9).all()


if __name__ == "__main__":
    import sys

    failures = []
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"[OK]   {name}")
            except Exception as exc:                       # noqa: BLE001 - test runner
                failures.append(name)
                print(f"[FAIL] {name}: {exc}")
    print(f"\n{len(failures)} failed" if failures else "\nall equivalence checks passed")
    sys.exit(1 if failures else 0)
