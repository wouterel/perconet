# This file is part of the PercoNet package
# (c) 2022 Eindhoven University of Technology
# Released under the BSD 3-clause license
# See LICENSE file for details
# Authors: Chiara Raffaelli, Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import numpy as np


class PeriodicNetwork:
    """Store and analyze the topology of a periodic net.

    Periodic nets are graphs embedded in a periodic topology. This class stores
    the topology of such a graph for the case of a threedimensional periodic
    box (a 3-torus) and allows to determine its percolation properties.
    The key question to determine percolation is whether the graph contains any
    loops that are topologically nontrivial in that they allow to travel from a
    node to itself with a net nonzero number of boundary crossings.

    Args:
        n (int):
            The number of nodes of the graph.
        max_degree (int):
            The largest number of edges coming out of any node.
    """
    def __init__(self, n: int, max_degree=6, verbose=True):
        if n < 1:
            raise ValueError("Number of nodes must be a positive integer.")
        self.number_of_nodes = n
        self.max_degree = max_degree
        # allocate arrays for per-node info
        self.boundary_crossing = np.zeros((n, max_degree, 3), dtype=int)
        self.neighbors = -1 * np.ones((n, max_degree), dtype=int)
        self.edges_list = -1 * np.ones((n, max_degree), dtype=int)
        self.neighbors_counter = np.zeros(n, dtype=int)
        self.edges_counter = 0  # keeps track of the total number of edges while building edge list
        self.simple_edges_list = []  # is this duplicate info?
        self.simple_boundary_crossing = []
        self.bond_is_across_boundary = []
        self.needs_reducing = 0
        self.crosses_boundaries = 0  # should I set anything? is there a way to keep ot empty?
        self.verbose = verbose

    def get_number_of_nodes(self):
        return self.number_of_nodes

    def add_edge(self, node1: int, node2: int, boundary_vector, existing_boundary_crossing_flag=False):
        # by default(for now, boundary crossing information is taken from boundary vector)
        """
        Add an edge to the periodic network
        
        Args:
            node1 (int):
                The index of the first node of the pair that defines this edge. Must satisfy (add requirement)
            node2 (int):
                The index of the second node of the pair that defines this edge. Must satisfy (add requirement)
            boundary_vector ((3) int):
                Vector of three integers (bvx,bvy,bvz) denoting the number of times the edge wraps around the x, y, and z boundaries, respectively.
                The sign indicates the wrapping direction (e.g. (-1,0,0) indicates that the edge goes around the x-boundary in the negative x-direction
                when going from node1 to node2.
            existing_boundary_crossing_flag (bool, optional):
                Deprecated parameter may be removed in future version.

        Returns:
            (bool): True if succesful. False if an error occurred.
        """
        if(node1>=self.number_of_nodes):
            print("Error in add_edge(): node {} does not exist (number of nodes = {})".format(node1,self.number_of_nodes))
            return False
        if(node2>=self.number_of_nodes):
            print("Error in add_edge(): node {} does not exist (number of nodes = {})".format(node2,self.number_of_nodes))
            return False
        if(self.neighbors_counter[node1] == self.max_degree):
            print(f"Cannot add edge: node {node1} already has {self.max_degree} edges")
            return False
        if(self.neighbors_counter[node2] == self.max_degree):
            print(f"Cannot add edge: node {node2} already has {self.max_degree} edges")
            return False
        if self.verbose:
            print("boundary vector is: ", boundary_vector, self.neighbors_counter)
        if not existing_boundary_crossing_flag:     #Improve way of counting
            if any (boundary_vector):
                self.bond_is_across_boundary.append(True) #bond/edge at the same index as "bond_counter" (see edges_list) is flagged as crossing the boundary
                self.crosses_boundaries +=1
            else:
                self.bond_is_across_boundary.append(False) #bond/edge at the same index as "bond_counter" (see edges_list) is flagged as not crossing the boundary
                self.needs_reducing +=1
        self.simple_edges_list.append([node1,node2])
        #this list constructor makes the next line work regardless of whether boundary_vector is passed as a Python list or numpy array
        self.simple_boundary_crossing.append([x for x in boundary_vector])
        
        self.neighbors[node1,self.neighbors_counter[node1]] = node2
        self.edges_list [node1,self.neighbors_counter[node1]] = self.edges_counter  #whatis a smart way to count "edge_count"?
        self.boundary_crossing [node1,self.neighbors_counter[node1], :] = boundary_vector
        self.neighbors_counter [node1] += 1
        
        if node2 != node1:
        # do the same inverting nodes
            self.neighbors[node2,self.neighbors_counter[node2]] = node1
            self.edges_list [node2,self.neighbors_counter[node2]] = self.edges_counter  #whatis a smart way to count "edge_count"?
            self.boundary_crossing [node2,self.neighbors_counter[node2], :] = -np.asarray(boundary_vector)  # fix this
            self.neighbors_counter [node2] += 1
        
        
        # update edges_counter
        self.edges_counter += 1
        return True
    
    def get_number_of_neighbors(self,i):
        """
        Get the number of bonds of node i
        
        Args:
            i (int): node number

        Returns:
            int: The number of edges (bonds) involving node i
        """

        return np.count_nonzero(self.neighbors[i,:] != -1)
        
    def get_neighbors(self, i, padded = True):
        """
        Get array of neighbor indices of node i.
        
        Args:
            i (int): node number
            padded (bool, optional): If true (the default), the list will be padded
                with values -1 to the value of maximum_neighbors_per_node passed to
                the constructor. Otherwise the length will be the number of neighbors of i.

        Returns:
            :obj:`list` of int: list of neighbors (along bonds) of node i

        """
        neighbors_of_i = self.neighbors[i, :]
        if padded:
            return neighbors_of_i
        if not padded:
            stripped_neighbors_of_i = neighbors_of_i [neighbors_of_i != -1]
            return stripped_neighbors_of_i

    def get_neighbor(self,i,n_index):
        """
        Get n_index'th neighbor of node i

        Args:
            i (int): node number
            n_index (int): the position of the neighbor in the list returned by get_neighbors()

        Returns:
            int: The index of that neighbor (the value of get_neighbors(i)[n_index])

        """
        return self.neighbors[i, n_index]

    def get_edges(self,i):
        """ Returns edge corresponding to """
        return self.edges_list [i, :]

    def get_number_of_edges(self):
        """ Returns number of edges """
        return self.edges_counter

    def get_edge(self,node1,k):
        """
        Get the edge number of the k'th edge of node1.

        Args:
            node1 (int): node number
            k (int): the position of the edge in the list of edges of node1.

        Returns:
            int: The edge number of that edge (to be used as an index in arrays of edge properties)
        """
        return self.edges_list [node1,k]
    
    def get_boundary_crossing(self,i,n_index):
        """ 
        Get the boundary crossing vector of the n_index'th neighbor of node i.
        
        Args:
            i (int): node number
            n_index (int): index of neighbor in neighbor list of i
        
        Returns:
            int[3]: The vector of integers [bx, by, bz] denoting the number of times each boundary is crossed by this edge.
        """
        return self.boundary_crossing [i,n_index, :]
    
    def __coloring(self, start,  current_color,color):
        color[start] = current_color
        for index in range(len(self.neighbors[1])):
            neigh = self.neighbors[start, index]
            edge_is_outside = self.bond_is_across_boundary [self.edges_list[start, index]] #check if bond is "inside" 
            if not edge_is_outside and neigh != -1 and color[neigh] == -1: #empty: not connected 
                #color[neigh] = current_color #commented this out because the first line of the recursed self.__coloring does this anyway
                self.__coloring(neigh,current_color,color)
                

    def cluster_find(self):
        """
        Obtain the cluster decomposition of the network (using only internal bonds).

        Returns:
            Tuple[:obj:`List` of int, int]: A list with the cluster ID of each node and the number of clusters
        """
        #initiate with first cluster colour
        current_color = 0
        # initialise list of star colours (-1 == not coloured)
        color = -1* np.ones(self.number_of_nodes, dtype = int)
        for star in range(self.number_of_nodes):
            if color[star] == -1:
                #start_cluster = star #don't need this because could just use star below
                self.__coloring (star, current_color,color)
                current_color += 1
        Ncolors = np.amax(color) + 1
        #print("THESE SHOULD BE EQUAL, RIGHT? {} {}".format(Ncolors,current_color))
        return color, Ncolors

    def nodeid_to_clusterid (self,list_colors):
        if self.verbose:
            print (self.simple_edges_list, self.simple_boundary_crossing)
        reduced_network = []
        for i in range(len(self.simple_edges_list)):
            if self.bond_is_across_boundary[i]:
                item1 = self.simple_edges_list[i][0]
                item2 = self.simple_edges_list[i][1]
                
                color1 = list_colors[item1]
                color2 = list_colors[item2]
                #add_edge here
                
                reduced_network.append ([color1, color2, self.simple_boundary_crossing[i][0], self.simple_boundary_crossing[i][1],self.simple_boundary_crossing[i][2]])
        
        reduced_network = np.asarray(reduced_network) 
        reduced_network = np.unique(reduced_network, axis=0)
        return reduced_network
    

    #methods that reduce network
    def get_reduced_network (self):
        """
        Generate the reduced network with identical boundary crossing properties but no internal edges.

        Returns:
            :obj:`PeriodicNetwork`: The reduced network    
        """    
        # to do: make sure self.needs_reducing works and return self if it is False
        #colouring algorithm (could be a generic function/method outside) that returns a list_of_colors and Ncolors
        list_colors, Ncolors = self.cluster_find()
        if self.verbose:
            print ("color", list_colors, Ncolors)
        #use list_of_colours and dropped list to turn dropped list into a coloured based list (clst.molid_to_clusterid) - 
            #this is now the format of the dropped lists used so far in testing 
        reduced_network_list = self.nodeid_to_clusterid (list_colors)
        if len(reduced_network_list)==0:
            #return network without any bonds, but with the proper number of nodes
            return PeriodicNetwork(Ncolors,1)
        #next line finds the number of occurrences of the most-occuring number in the first two columns of reduced_network_list
        largest_functionality = np.max(np.unique(reduced_network_list[:,0:2],return_counts=True)[1])
        
        #reduced_network = PeriodicNetwork(len(reduced_network_list), len(reduced_network_list)) 
        #reduced_network = PeriodicNetwork(Ncolors,len(reduced_network_list))  #
        reduced_network = PeriodicNetwork(Ncolors,largest_functionality,verbose=self.verbose)  #
        if self.verbose:
            print ("reduced network list:", reduced_network_list) 
        for i in range(len(reduced_network_list)):
            edge_info = reduced_network_list[i,:]
            if self.verbose:
                print(edge_info)
    
            reduced_network.add_edge(edge_info[0], edge_info[1], edge_info[2:])
        # create new instance of PeriodicNetwork and use add_edges
        return reduced_network
