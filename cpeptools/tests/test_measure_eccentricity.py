from cpeptools.metrics import calculate_eccentricity
import mdtraj as md
import numpy as np

def test_calculate_eccentricity():
    traj = md.load("cpeptools/tests/data/test.h5")

    out = []
    for result in calculate_eccentricity(traj):
        out.append(result)
    assert np.isnan(sum([i[0] for i in out]))
    # assert len(out) == 100
