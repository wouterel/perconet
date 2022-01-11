# This file is part of the PercoNet package
# (c) 2022 Eindhoven University of Technology
# Released under the BSD 3-clause license
# See LICENSE file for details
# Authors: Chiara Raffaelli, Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import numpy as np
import sympy


class LoopFinder:
    """
    Class implementing a depth-first search to determine the percolation directions of the network.

    Args:
        network (:obj:`perconet.PeriodicNetwork`): A PeriodicNetwork object representing the graph to analyze.
        verbose (bool, optional): Generate verbose output to stdout (to be replaced by Logging in future release)
    """
    def __init__(self,network,verbose=True):
        self.network=network
        #allocate arrays and values for dfs stuff
        self.discovery_time_node = 0
        self.discovery_time_edge = 0
        self.visited_nodes = -1 * np.ones (network.number_of_nodes, dtype = int)
        self.visited_edges = -1 * np.ones (np.amax(network.edges_list) + 1, dtype = int)
        #self.current_crossing = [0,0,0]
        self.crossing_sum = [[0,0,0]]
        self.loops_temp = []
        self.loops_list = []
        self.verbose = verbose
        
    def __dfs(self, start):
        #in here, get network info whenever using self.network.get_blabla()
        #this function can be private (not callable from the outside)
        self.visited_nodes [start] = self.discovery_time_node
        
        #print ("crossing sum now: \n", self.crossing_sum , "\n visited nodes: ", self.visited_nodes, "\n visited edges: \n", self.visited_edges)
        if self.verbose:
            print ("number of neighbours: ", self.network.get_number_of_neighbors(start), "start: ", start)
        
        for n_index in range(self.network.get_number_of_neighbors(start)): #fix this
            neigh = self.network.get_neighbor(start,n_index)
            edge = self.network.get_edge(start,n_index)
            if self.verbose:
                print ("edge", edge)
            if edge == -1:
                break
            
            if self.visited_edges[edge] == -1: # otherwhise neigh is a parent of start
                self.visited_edges[edge] = self.discovery_time_edge
                ("visited edges now: ", self.visited_edges)
                self.discovery_time_edge += 1
                timestep = self.visited_nodes [start]
                
                current_crossing = self.crossing_sum [timestep] + self.network.get_boundary_crossing(start, n_index)
                xC,yC,zC = current_crossing[:]
                #print ("current crossing sum for ",start,neigh,  ' (index: ', n_index, "), is ", xC,yC,zC, ", current crossing is ",self.network.get_boundary_crossing(start, n_index) )
                if self.visited_nodes[neigh] == -1:
                    self.discovery_time_node += 1
                    #self.visited_nodes[neigh] = self.discovery_time_node (no need to do this twice)
                    self.crossing_sum.append ([xC,yC,zC])
                    #print ("visited nodes: \n",visited_nodes, "\n crossing sum \n", crossing_sum)
                    #loops_temp, discovery_timeC, discovery_timeE = dfs (neigh, discovery_timeC, discovery_timeE, neighbors, visited_nodes, current_crossing, crossing_sum, boundary_crossing, edges_list, visited_edges, loops_temp, Nclusters)
                    self.__dfs(neigh ) #find right one
                    #print ("crossing sum now: \n", crossing_sum, "loops temp is:", loops_temp)
                else: #we found a loop!
                    
                    loop_timestep = self.visited_nodes [neigh]
                    
                    loop = current_crossing - self.crossing_sum [loop_timestep]

                    self.loops_temp.append (loop.tolist())
                    current_crossing -= self.network.get_boundary_crossing(start, n_index)
                    if self.verbose:
                        print("loop!  -> ", self.loops_temp, " loop_timestep: ", loop_timestep, ", crossing sum: ", self.crossing_sum )
        return
        
    def get_loops(self):
        #this is what used to be called launch_def
        #launch dfs and get unclean list
        #call some function that sweeps the list (function can be outside)
        #returns cleaned-up list of loops
        
        
        for node in range (self.network.number_of_nodes):
            if self.visited_nodes[node] == -1:           #should this be made into a method like is_node_visited (self, node) > return True/False or is this an overkill?
                if self.verbose:
                    print ("new origin of network: ", node)
                #reset stuff to 0
                self.crossing_sum = [[0,0,0]] #every time I start visiting a new network (that has never been visited before), I set the sum of crossings to 0
                #self.current_crossing = 0
                self.loops_temp = []
                self.discovery_time_node = 0
                self.discovery_time_edge = 0
                self.__dfs(node)
                
                self.discovery_time_node += 1
                if self.loops_temp: 
                
                    self.loops_list.extend(self.loops_temp)
                    
        return self.loops_list

    def get_independent_loops(self):
        """
        Generate a list of linearly independent topologically nontrivial loops.

        Returns:
            Tuple[:obj:`List` of :obj:`List` of int, int]:
                (list, int) A tuple containing a list of the independent loops,
                with element of the form [Bx, By, Bz] and the length of that list.
        """
        """
        ISSUE: The Matrix sweep used here is only valid if all loops have the same color in a full-graph coloring.
        FIX PLAN: Decompose the full graph into its connected subgraphs first, then determine loops for each color separately.
        NOTE: We should also start using a different term for the clusters as coloring has another meaning in graph theory.
        """
        #loops info should also provide a color for each loop to denote connected components
        # lump before clean would discard that info
        # clean before lump would allow having it
        myloops_list=np.asarray(self.get_loops())
        if self.verbose:
            print("loops list: ", myloops_list)
            print("rank according to numpy.linalg {}".format(np.linalg.matrix_rank(myloops_list,tol=1e-8)))
        _, inds = sympy.Matrix(myloops_list).T.rref()
        Nloops = len(inds)
        if self.verbose:
            print("Found {} loops".format(Nloops))
            print(list(inds))
            print(myloops_list)
        independent_loops = myloops_list[list(inds)] 
        #print ("independent_loops index and loops: ", independent_loops, inds)
        return independent_loops, Nloops



# # remove identical rows > should maybe be done somewhere in the code as well
# def rows_uniq_elems(a):
#     if len(a)==0:
#         return np.array([])
#     new_array = [tuple(row) for row in a]
#     uniques = np.unique(new_array, axis = 0)
#     #uniques = np.unique(a, axis = 0)
#     return uniques



 