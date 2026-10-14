import numpy as np
import pandas as pd
from src.matching.match_truth import match_batch
from src.truthset.intervals import cluster_ro, reciprocal_overlap


def test_ro_basic():
    assert reciprocal_overlap(0, 100, 50, 150) == 0.5
    assert reciprocal_overlap(0, 100, 0, 100) == 1.0
    assert reciprocal_overlap(0, 100, 200, 300) == 0.0
    assert abs(reciprocal_overlap(0, 1000, 0, 100) - 0.1) < 1e-9


def test_cluster_components():
    chrom = np.array(["a", "a", "a", "b"])
    start = np.array([0, 10, 5000, 0])
    end = np.array([1000, 1010, 6000, 1000])
    ids = cluster_ro(chrom, start, end, 0.8)
    assert ids[0] == ids[1] and ids[0] != ids[2] and ids[0] != ids[3]


def test_cluster_requires_same_key():
    ids = cluster_ro(np.array(["a", "a"]), np.array([0, 0]), np.array([100, 100]), 0.8, same_key=np.array(["x", "y"]))
    assert ids[0] != ids[1]


def test_match_type_aware():
    truth = pd.DataFrame({"sv_id": ["s1", "s2"], "chrom": ["c", "c"], "start": [1000, 50000], "end": [2000, 60000],
                          "svtype": ["DEL", "DUP"]})
    calls = pd.DataFrame({"chrom": ["c", "c"], "start": [1010, 50000], "end": [2005, 60000], "svtype": ["DEL", "DEL"]})
    m = match_batch(truth, calls, "std")
    assert m.loc[0, "matched"] and not m.loc[1, "matched"]  # DUP not matched by a DEL call
    assert m.loc[0, "bp_err_start"] == 10 and m.loc[0, "bp_err_end"] == 5
