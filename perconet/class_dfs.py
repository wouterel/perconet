# -*- coding: utf-8 -*-
"""
authors: Chiara Raffaelli, Wouter G. Ellenbroek
Maintainer: w.g.ellenbroek@tue.nl
"""
import numpy as np
import sympy


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

    def add_edge(self, node1, node2, boundary_vector, existing_boundary_crossing_flag=False):
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
        
    def get_neighbors(self,i, padded = True):
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
        reduced_network = rows_uniq_elems(reduced_network)
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



# remove identical rows > should maybe be done somewhere in the code as well
def rows_uniq_elems(a):
    if len(a)==0:
        return np.array([])
    new_array = [tuple(row) for row in a]
    uniques = np.unique(new_array, axis = 0)
    return uniques



 