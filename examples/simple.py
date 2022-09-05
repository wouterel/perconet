# Usage example for perconet package
# See https://github.com/wouterel/perconet

import perconet as pn
import numpy as np


def test_simple():
    # the example starts counting node numbers from 1,
    # but perconet counts from 0, so we define 5 nodes instead of 4.
    # The unused node "0" does not affect the percolation properties.
    number_of_nodes = 5
    max_coordination = 6
    testnet = pn.PeriodicNetwork(number_of_nodes,
                                 max_coordination,
                                 verbose=False)

    # first add the three internal edges connecting 1-2, 2-3, and 3-4.
    testnet.add_edge(1, 2, np.array([0, 0, 0]))
    testnet.add_edge(2, 3, np.array([0, 0, 0]))
    testnet.add_edge(3, 4, np.array([0, 0, 0]))

    # we now have a small network that doesn't do anything with the boundaries yet
    loopfinder = pn.LoopFinder(testnet, verbose=False)
    loops, Nloops = loopfinder.get_independent_loops()
    print(f"Found {Nloops} loops ( = 0 because no boundary-crossing bonds are defined).")

    print("Adding the boundary-crossing bonds")
    # Note the sign of the boundary crossing for an edge between i and j
    # is determined by the direction in which you go if you follow the edge from i to j
    # add a bond between nodes 1 and 3 that crosses the x-boundary
    testnet.add_edge(1, 3, np.array([1, 0, 0]))
    # add a bond between nodes 1 and 4 that crosses the negative y-boundary
    testnet.add_edge(1, 4, np.array([0, -1, 0]))
    # add a bond between nodes 2 and 4 that crosses the negative y-boundary
    testnet.add_edge(2, 4, np.array([0, -1, 0]))
    # add a bond between nodes 3 and 4 that crosses the negative x-boundary
    testnet.add_edge(3, 4, np.array([-1, 0, 0]))

    # now the network percolates across x and y boundaries.
    loops, Nloops = loopfinder.get_independent_loops()
    print(f"Found {Nloops} independent loops.")
    for loop in loops:
        print(f"Loop: {loop}")


if __name__ == "__main__":
    test_simple()
