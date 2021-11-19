#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 18:04:33 2020

@authors: craffaelli,wouterel
"""

import sys
# sys.path.append("..")
# removed this because within the poetry setup it is not necessary
import perconet as pn
import numpy as np
sys.path.append("tests")
from networkdata_testing import edgelist

dropped_list = edgelist  # use old ropped_list" terminology below
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
        neigh = my_test_network.get_neighbor(i, n_index)
        edge = my_test_network.get_edge(i, n_index)
        print("edge", edge)

if not my_test_network.crosses_boundaries:
    print("Network does not cross periodic boundaries in any direction, so it does not percolate")
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
