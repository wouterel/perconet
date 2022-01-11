
"""
Created on Fri Apr 10 18:04:33 2020

@authors: craffaelli,wouterel
"""

# import sys
# sys.path.append("..")
# removed this because within the poetry setup it is not necessary
import perconet as pn
import numpy as np
# sys.path.append("tests")
import networkdata_testing as ndata
from pn_test_helpers import initialize_test


def test_loops():
    for net in dir(ndata):
        if not net.startswith("testcase_"):
            continue
        print(net)
        edgelist, solution = eval(f"ndata.{net}")
        network = initialize_test(edgelist)
        loopfinder = pn.LoopFinder(network, verbose=False)
        loops, n_loops = loopfinder.get_independent_loops()
        print(loops)
        print(net)
        if len(solution) == 0:  # solution = []
            assert n_loops == 0
        elif len(solution) == 1 and len(solution[0]) == 0:  # solution = [[]]
            assert n_loops == 0
        else:
            assert n_loops > 0
            # next assert is too strict.
            # Future optimizations could change the precise list of loops
            # without changing their meaning (e.g. as linear combinations)
            # Fix later.
            assert np.all(np.asarray(solution) == np.asarray(loops))
    # assert False


def test_reduction_C():
    edgelist, solution = ndata.testcase_C
    assert len(edgelist) == 23
    network = initialize_test(edgelist)
    assert network.get_number_of_nodes() == 10
    assert network.get_number_of_edges() == len(edgelist)  # =23
    network = network.get_reduced_network()
    assert network.get_number_of_nodes() == 5
    assert network.get_number_of_edges() == 14
    loopfinder = pn.LoopFinder(network, verbose=False)
    loops, n_loops = loopfinder.get_independent_loops()
    assert n_loops == 3


def test_duplicate_removal_C():
    edgelist, solution = ndata.testcase_C
    assert len(edgelist) == 23
    edgelist = np.unique(edgelist, axis=0)
    # this used to test an exposed routine, now it just tests a numpy feature
    assert len(edgelist) == 22
    network = initialize_test(edgelist)
    assert network.get_number_of_nodes() == 10
    assert network.get_number_of_edges() == len(edgelist)  # =22
    network = network.get_reduced_network()
    assert network.get_number_of_nodes() == 5
    assert network.get_number_of_edges() == 14
    loopfinder = pn.LoopFinder(network, verbose=False)
    loops, n_loops = loopfinder.get_independent_loops()
    assert n_loops == 3


def oldstuff():
    dropped_list = ndata.edgelist  # use old "dropped_list" terminology below
    # dropped_list has description of boundary-crossing bonds in a reduced network
    # these fully characterize the network
    # deduce number of nodes from contents of dropped_list
    number_of_nodes = np.amax(dropped_list) + 1
    my_test_network = pn.PeriodicNetwork(number_of_nodes, number_of_nodes)

    print("#edges including duplicates:", len(dropped_list))
    dropped_list = pn.rows_uniq_elems(dropped_list)
    print("#edges after removing duplicates:", len(dropped_list))

    for edge_info in dropped_list:
        print(edge_info)
        my_test_network.add_edge(edge_info[0], edge_info[1], edge_info[2:])

    print("neighbors", my_test_network.neighbors)
    print("boundary crossing ", my_test_network.bond_is_across_boundary)
    for i in range(len(my_test_network.neighbors)):
        print("number of neighbours: ", my_test_network.get_number_of_neighbors(i))
        for n_index in range(my_test_network.get_number_of_neighbors(i)):
            # neigh = my_test_network.get_neighbor(i, n_index)
            edge = my_test_network.get_edge(i, n_index)
            print("edge", edge)

    if not my_test_network.crosses_boundaries:
        print("Network does not cross any periodic boundary, so it does not percolate")
        exit()

    # my_test_network.needs_reducing=0
    # the following lines could already be part of get_reduced_network
    if my_test_network.needs_reducing:
        print("Reducing network to simpler form...")
        my_reduced_network = my_test_network.get_reduced_network()
        # think of a way to return it already in the right format you would get with add_edges
        myloops = pn.LoopFinder(my_reduced_network)
    else:
        myloops = pn.LoopFinder(my_test_network)
    loops, n_loops = myloops.get_independent_loops()

    print("number of loops \n", n_loops, " \n independent loops: \n", loops)
