from cpeptools.metrics import identify_interference

import pytest
import mdtraj as md


def test_one_frame():
    traj = md.load("cpeptools/tests/data/test.h5")[0:100]
    identify_interference(traj)
