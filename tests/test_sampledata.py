# This file is part of the perconet package
# (c) 2022 Eindhoven University of Technology
# Released under EUPL v1.2
# See LICENSE file for details
# Contributors:
# * Chiara Raffaelli
# * Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import perconet as pn
import numpy as np
from pn_test_helpers import initialize_test


def single_test_set(set_id):
    nloops = np.loadtxt(f"tests/testdata/bonds_{set_id}_loops.dat", dtype=int)
    nfiles = len(nloops)
    for ifile in range(nfiles):
        filename = f"tests/testdata/bonds_{set_id}_{ifile:02}.dat"
        bondlist = np.loadtxt(filename, dtype=int)
        network = initialize_test(bondlist)
        loopfinder = pn.LoopFinder(network, verbose=False)
        _, n_loops = loopfinder.get_independent_loops()
        # _, n_loops = loopfinder._compare_independence_methods()
        assert n_loops == nloops[ifile]
        print(n_loops)


def test_samples():
    sets = [0, 20, 40, 50, 60, 80]
    for set in sets:
        single_test_set(set)


if __name__ == "__main__":
    test_samples()
