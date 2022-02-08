# This file is part of the PercoNet package
# (c) 2022 Eindhoven University of Technology
# Released under the BSD 3-clause license
# See LICENSE file for details
# Authors: Chiara Raffaelli, Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import numpy as np
# import sympy
import perconet.looptools as looptools


class LoopFinder:
    """
    Class implementing a depth-first search to determine the percolation directions of the network.

    Args:
        network (:obj:`perconet.PeriodicNetwork`): A PeriodicNetwork object representing
            the graph to analyze.
        verbose (bool, optional): Generate verbose output to stdout (to be replaced by
            Logging in future release)
    """
    def __init__(self, network, verbose=True):
        self.network = network
        # allocate arrays and values for dfs stuff
        self.discovery_time_node = 0
        self.discovery_time_edge = 0
        self.visited_nodes = -1 * np.ones(network.number_of_nodes, dtype=int)
        self.visited_edges = -1 * np.ones(np.amax(network.edges_list) + 1, dtype=int)
        self.crossing_sum = [[0, 0, 0]]
        self.loops_temp = []
        self.loops_list = []
        self.verbose = verbose

    def __dfs(self, start):
        # in here, get network info whenever using self.network.get_blabla()
        # this function can be private (not callable from the outside)
        self.visited_nodes[start] = self.discovery_time_node

        # print ("crossing sum now: \n", self.crossing_sum , "\n visited nodes: ",
        #         self.visited_nodes, "\n visited edges: \n", self.visited_edges)
        if self.verbose:
            print("number of neighbours: ", self.network.get_number_of_neighbors(start),
                  "start: ", start)

        for n_index in range(self.network.get_number_of_neighbors(start)):  # fix this
            neigh = self.network.get_neighbor(start, n_index)
            edge = self.network.get_edge(start, n_index)
            if self.verbose:
                print("edge", edge)
            if edge == -1:
                break

            if self.visited_edges[edge] == -1:  # otherwhise neigh is a parent of start
                self.visited_edges[edge] = self.discovery_time_edge
                ("visited edges now: ", self.visited_edges)
                self.discovery_time_edge += 1
                timestep = self.visited_nodes[start]

                current_crossing = self.crossing_sum[timestep] + \
                    self.network.get_boundary_crossing(start, n_index)
                xC, yC, zC = current_crossing[:]
                # print(f"Current crossing is {self.network.get_boundary_crossing(start, n_index)}")
                # print(f"    Total crossing for {start}-{neigh} (index: {n_index}): {[xC,yC,zC]}.")

                if self.visited_nodes[neigh] == -1:
                    self.discovery_time_node += 1
                    self.crossing_sum.append([xC, yC, zC])
                    # print ("visited nodes: \n",visited_nodes, "\n crossing sum \n", crossing_sum)
                    # loops_temp, discovery_timeC, discovery_timeE = dfs (neigh, discovery_timeC,
                    #  discovery_timeE, neighbors, visited_nodes, current_crossing, crossing_sum,
                    # boundary_crossing, edges_list, visited_edges, loops_temp, Nclusters)
                    self.__dfs(neigh)  # find right one
                    # print ("crossing sum now: \n", crossing_sum, "loops temp is:", loops_temp)
                else:  # we found a loop!

                    loop_timestep = self.visited_nodes[neigh]

                    loop = current_crossing - self.crossing_sum[loop_timestep]

                    self.loops_temp.append(loop.tolist())
                    current_crossing -= self.network.get_boundary_crossing(start, n_index)
                    if self.verbose:
                        print("loop!  -> ", self.loops_temp, " loop_timestep: ", loop_timestep,
                              ", crossing sum: ", self.crossing_sum)
        return

    def get_loops(self):
        """
        Generate a raw list of loops. Most use cases will require get_independent_loops instead.

        Returns:
            Tuple[:obj:`List` of :obj:`List` of int, int]:
                (list, int) A tuple containing a raw list of loops,
                with elements of the form [Bx, By, Bz] and the length of that list.
        """

        original_network = self.network
        # print(f"needs red: {self.network.needs_reducing()}")
        self.network = self.network.get_reduced_network()
        # print(f"Number of nodes after reduction: {network.get_number_of_nodes()}.")

        # starting_nodes_with_loops = 0
        for node in range(self.network.number_of_nodes):
            # If this node has not been visited, this is a new cluster and we start a fresh search
            if self.visited_nodes[node] == -1:
                if self.verbose:
                    print("new origin of network: ", node)
                # reset variables used during search to 0
                self.crossing_sum = [[0, 0, 0]]
                self.loops_temp = []
                self.discovery_time_node = 0
                self.discovery_time_edge = 0
                # enter the recursive search
                self.__dfs(node)

                self.discovery_time_node += 1
                if self.loops_temp:
                    # starting_nodes_with_loops += 1
                    self.loops_list.extend(self.loops_temp)
        self.network = original_network  # this means no more extracting info after this point
        # print(f"Loops found from {starting_nodes_with_loops} starting points")
        # assert starting_nodes_with_loops < 2
        return self.loops_list

    def get_independent_loops(self):
        """
        Generate a list of linearly independent topologically nontrivial loops.

        Returns:
            Tuple[:obj:`List` of :obj:`List` of int, int]:
                (list, int) A tuple containing a list of the independent loops,
                with elements of the form [Bx, By, Bz] and the length of that list.
        """
        # loops info should include a cluster label for each loop to denote connected components
        # lump before clean would discard that info
        # clean before lump would allow having it
        # Note that we have not yet encountered data for which this difference mattered
        # and if all one wants to know is the weak directions of a material it
        # is irrelevant anyway
        myloops_list = np.asarray(self.get_loops())
        if len(myloops_list) == 0:
            # No loops found
            return myloops_list, 0
        if np.count_nonzero(myloops_list) == 0:
            # If only topologically trivial loops we treat this as no loops
            return np.asarray([]), 0
        if self.verbose:
            print("loops list: ", myloops_list)
            print(f"rank according to numpy.linalg {np.linalg.matrix_rank(myloops_list, tol=1e-8)}")
        independent_loops, _, Nloops = looptools.integer_gaussian_elimination(myloops_list)
        # the ige method has left the 0-rows in there so need to remove them now
        independent_loops = independent_loops[:Nloops]

        if self.verbose:
            print(f"Found {Nloops} loops")
            print(independent_loops)
        return independent_loops, Nloops

    def _compare_independence_methods(self):
        """
        Checks if ranks obtained with three methods agree.
        """
        import sympy  # put this here so sympy will be a development-only dependency
        myloops_list = np.asarray(self.get_loops())
        if len(myloops_list) == 0:
            # No loops found
            return myloops_list, 0
        if np.count_nonzero(myloops_list) == 0:
            # If only topologically trivial loops we treat this as no loops
            return np.asarray([]), 0
        nprank = np.linalg.matrix_rank(myloops_list, tol=1e-8)
        if self.verbose:
            print("loops list: ", myloops_list)
            print(f"rank according to numpy.linalg {nprank}")

        _, inds = sympy.Matrix(myloops_list).T.rref()
        Nloops = len(inds)
        assert Nloops == nprank
        if self.verbose:
            print(f"Found {Nloops} loops")
            print(list(inds))
            print(myloops_list)
        independent_loops = myloops_list[list(inds)]

        my_ind_loops, _, my_Nloops = looptools.integer_gaussian_elimination(myloops_list)
        # the ige method has left the 0-rows in there so need to remove them now
        my_ind_loops = my_ind_loops[:my_Nloops]
        assert my_Nloops == Nloops
        # print("independent loops according to sympy:")
        # print(independent_loops)
        # print("independent loops according to ige:")
        # print(my_ind_loops)
        # print("independent_loops index and loops: ", independent_loops, inds)
        # assert not break_here
        assert my_Nloops > 0
        if my_Nloops == 1:
            assert np.all(my_ind_loops == independent_loops)
        # next 4 lines could give false alarms but worked fine on existing test data
        elif my_Nloops == 2:
            print(my_ind_loops)
            for loop in my_ind_loops.tolist():
                assert loop in independent_loops.tolist()
        return independent_loops, Nloops
