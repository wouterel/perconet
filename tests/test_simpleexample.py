# import sys
# sys.path.append("..")
import perconet as pn
import numpy as np


def test_simple():

    number_of_nodes = 4
    max_coordination = 6
    testnet = pn.PeriodicNetwork(number_of_nodes,
                                 max_coordination,
                                 verbose=False)

    # add a bond between nodes 0 and 1 that crosses no boundary
    testnet.add_edge(0, 1, np.array([0, 0, 0]))
    # add a bond between nodes 2 and 1 that crosses no boundary
    testnet.add_edge(2, 1, np.array([0, 0, 0]))
    # add a bond between nodes 2 and 3 that crosses no boundary
    testnet.add_edge(2, 3, np.array([0, 0, 0]))
    # add a bond between nodes 0 and 3 that crosses no boundary
    testnet.add_edge(0, 3, np.array([0, 0, 0]))

    # we now have a square that doesn't do anything

    if not testnet.crosses_boundaries:
        print("Network does not cross periodic boundaries in any direction, "
              "therefore it does not percolate")

    # add a bond between nodes 1 and 0 that crosses the x-boundary
    testnet.add_edge(1, 0, np.array([1, 0, 0]))
    # add a bond between nodes 1 and 0 that crosses the x-boundary
    testnet.add_edge(2, 3, np.array([1, 0, 0]))
    # add a bond between nodes 3 and 0 that crosses the negative y-boundary
    testnet.add_edge(3, 0, np.array([0, -1, 0]))

    if testnet.needs_reducing:
        print("Reducing network to simpler form...")
        my_reduced_network = testnet.get_reduced_network()
        # print(my_reduced_network)
        myloops = pn.LoopFinder(my_reduced_network, verbose=False)
    else:
        myloops = pn.LoopFinder(testnet, verbose=False)
    loops = myloops.get_independent_loops()

    assert len(loops[0]) == 2  # tests number of loops (2)
    loop1 = loops[0][0]
    loop2 = loops[0][1]
    assert loop1[2] == 0  # tests if loop1 is in the xy plane
    assert loop2[2] == 0  # tests if loop2 is in the xy plane

    def dotprod_integer(a, b):
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

    assert dotprod_integer(loop1, loop1) == 1  # tests if loop1 has length 1
    assert dotprod_integer(loop2, loop2) == 1  # tests if loop2 has length 1

    crossprod = loop1[1]*loop2[0]-loop1[0]*loop2[1]
    assert abs(crossprod) == 1  # tests if loops span the xy plane

    print("number of loops \n", len(loops[0]),
          " \n independent loops: \n", loops[0])


if __name__ == "__main__":
    test_simple()
