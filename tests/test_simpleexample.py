# import sys
# sys.path.append("..")
import perconet as pn
import numpy as np
import pytest


def test_empty():
    with pytest.raises(ValueError):
        pn.PeriodicNetwork(0, max_degree=3)


def test_maxdegree_single():
    # note: we save self-bonds only ones but technically the node has degree 2,
    # so technically perhaps we should have this fail for not satisfying max_degree
    testnet = pn.PeriodicNetwork(1, max_degree=1)
    assert testnet.add_edge(0, 0, [0, 0, 0])
    assert not testnet.add_edge(0, 0, [1, 0, 0])


def test_maxdegree_twonodes():
    testnet = pn.PeriodicNetwork(2, max_degree=2)
    assert testnet.add_edge(1, 1, [0, 0, 0])
    assert testnet.add_edge(0, 1, [1, 0, 0])
    assert not testnet.add_edge(0, 1, [0, 0, 0])


def test_loops_singlenode():
    testnet = pn.PeriodicNetwork(1, max_degree=1)
    testbc = np.random.randint(-5, 6, 3)
    # Is it bad form to use random numbers in a test?
    assert testnet.add_edge(0, 0, testbc)
    finder = pn.LoopFinder(testnet, verbose=False)
    loops, n_loops = finder.get_independent_loops()
    # the algorithm would also be correct if this returned loops = -testbc
    # at the time of writing of this test it always gives loops = testbc though
    assert n_loops == 1
    assert np.all(loops[0] == testbc) or np.all(loops[0] == -testbc)


def test_loops_reversewrap():
    testnet = pn.PeriodicNetwork(2, max_degree=2)
    assert testnet.add_edge(1, 0, [2, 0, 0])
    assert testnet.add_edge(0, 1, [-2, 0, 0])
    finder = pn.LoopFinder(testnet, verbose=False)
    loops, n_loops = finder.get_independent_loops()
    assert n_loops == 0


def test_loops_twonodes():
    testnet = pn.PeriodicNetwork(2, max_degree=3)
    assert testnet.add_edge(1, 0, [2, 0, 0])
    assert testnet.add_edge(0, 1, [0, 1, 0])
    finder = pn.LoopFinder(testnet)
    loops, n_loops = finder.get_independent_loops()
    # the algorithm would also be correct if this returned loops = -testbc
    # at the time of writing of this test it always gives loops = testbc though
    assert n_loops == 1
    checkloop = np.asarray([-2, -1, 0])
    assert np.all(loops[0] == checkloop) or np.all(loops[0] == -checkloop)
    assert testnet.add_edge(0, 1, [0, 0, 0])
    finder = pn.LoopFinder(testnet)
    loops, n_loops = finder.get_independent_loops()
    assert n_loops == 2
    l0 = loops[0]
    l1 = loops[1]
    # check if the two loops have no z-component
    assert l0[2] == 0
    assert l1[2] == 0
    # now check if the two loops span the xy plane
    assert l0[0]*l1[1] != l0[1]*l1[0]
    # During the writing of this test I found that adding an edge to a network after
    # initializing the loopfinder will not update the loopfinder's network!


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
    assert testnet.needs_reducing()
    assert not testnet.crosses_boundaries()

    # add a bond between nodes 1 and 0 that crosses the x-boundary
    testnet.add_edge(1, 0, np.array([1, 0, 0]))
    assert testnet.crosses_boundaries()
    # add a bond between nodes 1 and 0 that crosses the x-boundary
    testnet.add_edge(2, 3, np.array([1, 0, 0]))
    # add a bond between nodes 3 and 0 that crosses the negative y-boundary
    testnet.add_edge(3, 0, np.array([0, -1, 0]))

    assert testnet.needs_reducing
    my_reduced_network = testnet.get_reduced_network()
    assert my_reduced_network.crosses_boundaries()
    assert not my_reduced_network.needs_reducing()
    myloops = pn.LoopFinder(my_reduced_network, verbose=False)
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


def test_5d():
    number_of_nodes = 4
    max_coordination = 6
    testnet = pn.PeriodicNetwork(number_of_nodes,
                                 max_coordination,
                                 verbose=False,
                                 dim=5)
    # check if bc vector with wrong length is rejected
    # both for numpy array and python list
    assert not testnet.add_edge(0, 1, np.array([0, 0, 0, 1]))
    assert not testnet.add_edge(0, 1, [0, 1, 1])
    # check if bond indeed was not added
    assert testnet.get_number_of_edges() == 0

    testnet.add_edge(0, 1, np.array([0, 0, 0, 1, 0]))
    testnet.add_edge(1, 2, np.array([0, 0, 0, 0, -1]))
    testnet.add_edge(2, 3, np.array([0, 0, 1, 0, 0]))
    testnet.add_edge(3, 0, np.array([0, 0, 0, 1, 1]))
    testnet.add_edge(1, 3, np.array([1, 0, 0, 0, 0]))
    assert not testnet.needs_reducing()
    assert testnet.crosses_boundaries()

    myloops = pn.LoopFinder(testnet, verbose=False)
    loops, n_loops = myloops.get_independent_loops()
    compare = [np.array([1, 0, 0, 2, 1]),
               np.array([0, 0, 1, 2, 0])]
    # print(loops)
    # print(compare)
    for i_loop, loop in enumerate(loops):
        print(f"loop {i_loop}: {loop}")
        print(f"comparing to {compare[i_loop]}")
        assert np.all(loop == compare[i_loop]) or np.all(loop == -compare[i_loop])

    testnet.add_edge(0, 2, [0, 0, 0, 0, 0])
    assert testnet.needs_reducing()
    my_reduced_network = testnet.get_reduced_network()
    assert my_reduced_network.get_number_of_nodes() == 3
    assert my_reduced_network.crosses_boundaries()
    assert not my_reduced_network.needs_reducing()
    loops, n_loops = myloops.get_independent_loops()
    assert n_loops == 3
    # not adding an assertion on the loops themselves
    # because the ige result is not unique
    # using hermite normal form for this in the future will help
    # print(n_loops, loops)
    return


if __name__ == "__main__":
    test_5d()
