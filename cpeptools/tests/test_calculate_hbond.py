from cpeptools.metrics import calculate_hbond_occurences
import mdtraj as md

def test_calculate_hbond():
    traj = md.load("cpeptools/tests/data/test.h5")

    out =  calculate_hbond_occurences(traj)
    assert "occurences" in out and "labels" in out and "num_frames" in out
