#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 18:04:33 2020

@author: craffaelli
"""
import numpy as np
from dropped_list_test import dropped_list  
import sympy


class PeriodicNetwork:
    #stores network structure
    #contains methods like get_something() to provide network info
    def __init__(self,n,maximum_neighbors_per_node):
        #allocate arrays for per-node info
        
        #(substitute function neighbour_clusters in clst)
        self.number_of_clusters = n
        self.boundary_crossing = np.zeros((n,maximum_neighbors_per_node,3), dtype = int)
        self.neighbors = -1 * np.ones((n,maximum_neighbors_per_node), dtype = int)
        self.edges_list = -1 * np.ones((n,maximum_neighbors_per_node), dtype = int)
        self.neighbors_counter = np.zeros(n, dtype = int)
        self.edges_counter = 0 #keeps track of the total number of edges while building the edges list
        self.simple_edges_list = []  #is this duplicate info?
        self.simple_boundary_crossing = []
        self.bond_is_across_boundary = []
        self.needs_reducing = 0
        self.crosses_boundaries = 0 #should I set anything? is there a way to keep ot empty?
        
    def add_edge(self,node1,node2, boundary_vector, existing_boundary_crossing_flag = False): #by default(for now, boundary crossing information is taken from boundary vector)
        """ Adds info for this edge """
        print("boundary vector is: ", boundary_vector, self.neighbors_counter)
        if not existing_boundary_crossing_flag:     #Improve way of counting
            if any (boundary_vector):
                self.bond_is_across_boundary.append(True) #bond/edge at the same index as "bond_counter" (see edges_list) is flagged as crossing the boundary
                #self.bond_is_across_boundary[node1,self.neighbors_counter[node1]] = self.bond_is_across_boundary[node2,self.neighbors_counter[node2]] = True
                self.crosses_boundaries +=1
            else:
                self.bond_is_across_boundary.append(False) #bond/edge at the same index as "bond_counter" (see edges_list) is flagged as not crossing the boundary
                #self.bond_is_across_boundary[node1,self.neighbors_counter[node1]] = self.bond_is_across_boundary[node2,self.neighbors_counter[node2]] = False
                self.needs_reducing +=1
        self.simple_edges_list.append([node1,node2])
        self.simple_boundary_crossing.append([x for x in boundary_vector])
        
        self.neighbors[node1,self.neighbors_counter[node1]] = node2
        self.edges_list [node1,self.neighbors_counter[node1]] = self.edges_counter  #whatis a smart way to count "edge_count"?
        self.boundary_crossing [node1,self.neighbors_counter[node1], :] = boundary_vector
        self.neighbors_counter [node1] += 1
        
        if node2 != node1:
        # do the same inverting nodes
            self.neighbors[node2,self.neighbors_counter[node2]] = node1
            self.edges_list [node2,self.neighbors_counter[node2]] = self.edges_counter  #whatis a smart way to count "edge_count"?
            self.boundary_crossing [node2,self.neighbors_counter[node2], :] = -boundary_vector
            self.neighbors_counter [node2] += 1
        
        
        # update edges_counter
        self.edges_counter += 1
    

    def get_number_of_neighbors(self,i):
        """ Returns number of neighbors of node i """
        return np.count_nonzero(self.neighbors[i,:] != -1)
        
    def get_neighbors(self,i, padded = True):
        """ Returns array of neighbor indices of node i (could be length maxN, padded with -1s, or length actual number)"""
        neighbors_of_i = self.neighbors[i, :]
        if padded:
            return neighbors_of_i
        if not padded:
            stripped_neighbors_of_i = neighbors_of_i [neighbors_of_i != -1]
            return stripped_neighbors_of_i
             
    def get_neighbor(self,i,n_index):
        """ Returns n_index'th neighbor of node i """
        return self.neighbors[i, n_index]
    
    def get_edges(self,i):
        """ Returns edge corresponding to """
        return self.edges_list [i, :]
    
    def get_edge(self,node1,node2):
        """ Returns edge between two nodes  ##check if edge exists! """
        return self.edges_list [node1,node2]
    
    def get_boundary_crossing(self,i,n_index):
        """ Returns the bc vector of the n_index'th neighbor of node i """
        return self.boundary_crossing [i,n_index, :]
    
    def coloring(self, start,  current_color,color):
        color[start] = current_color
        for index in range(len(self.neighbors[1])):
            neigh = self.neighbors[start, index]
            edge_is_outside = self.bond_is_across_boundary [self.edges_list[start, index]] #check if bond is "inside" 
            if not edge_is_outside and neigh != -1 and color[neigh] == -1: #empty: not connected 
                color[neigh] = current_color
                PeriodicNetwork.coloring(self, neigh,current_color,color)
                

    def cluster_find(self):  
        #initiate with first cluster colour
        current_color = 0
        # initialise list of star colours (-1 == not coloured)
        color = -1* np.ones(self.number_of_clusters, dtype = int)
        for star in range(self.number_of_clusters):
            if color[star] == -1:
                start_cluster = star
                PeriodicNetwork.coloring (self, start_cluster, current_color,color)
                current_color += 1
        Ncolors = np.amax(color) + 1
        return color, Ncolors

    def nodeid_to_clusterid (self,list_colors):
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
        
        #colouring algorithm (could be a generic function/method outside) that returns a list_of_colors and Ncolors
        list_colors, Ncolors = self.cluster_find()
        print ("color", list_colors, Ncolors)
        #use list_of_colours and dropped list to turn dropped list into a coloured based list (clst.molid_to_clusterid) - 
            #this is now the format of the dropped lists used so far in testing 
        reduced_network_list = self.nodeid_to_clusterid (list_colors)
        reduced_network = PeriodicNetwork(len(reduced_network_list), len(reduced_network_list)) 
        #reduced_network = PeriodicNetwork(Ncolors,len(reduced_network_list))  #
        print ("reduced network list:", reduced_network_list) 
        for i in range(len(reduced_network_list)):
            edge_info = reduced_network_list[i,:]
            print(edge_info)
    
            reduced_network.add_edge(edge_info[0], edge_info[1], edge_info[2:])
        # create new instance of PeriodicNetwork and use add_edges
        return reduced_network
        

class LoopFinder:
    def __init__(self,network):
        self.network=network
        #allocate arrays and values for dfs stuff
        self.discovery_time_node = 0
        self.discovery_time_edge = 0
        self.visited_nodes = -1 * np.ones (network.number_of_clusters, dtype = int)
        self.visited_edges = -1 * np.ones (np.amax(network.edges_list) + 1, dtype = int)
        #self.current_crossing = [0,0,0]
        self.crossing_sum = [[0,0,0]]
        self.loops_temp = []
        self.loops_list = []
        
    def dfs(self, start):
        #in here, get network info whenever using self.network.get_blabla()
        #this function can be private (not callable from the outside)
        self.visited_nodes [start] = self.discovery_time_node
        
        #print ("crossing sum now: \n", self.crossing_sum , "\n visited nodes: ", self.visited_nodes, "\n visited edges: \n", self.visited_edges)
        print ("number of neighbours: ", self.network.get_number_of_neighbors(start), "start: ", start)
        
        for n_index in range(self.network.get_number_of_neighbors(start)): #fix this
            neigh = self.network.get_neighbor(start,n_index)
            edge = self.network.get_edge(start,n_index)
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
                    LoopFinder.dfs(self,neigh ) #find right one
                    #print ("crossing sum now: \n", crossing_sum, "loops temp is:", loops_temp)
                else: #we found a loop!
                    
                    loop_timestep = self.visited_nodes [neigh]
                    
                    loop = current_crossing - self.crossing_sum [loop_timestep]

                    self.loops_temp.append (loop.tolist())
                    current_crossing -= self.network.get_boundary_crossing(start, n_index)
                    print("loop!  -> ", self.loops_temp, " loop_timestep: ", loop_timestep, ", crossing sum: ", self.crossing_sum )
        return
        
    
    
    def get_loops(self):
        #this is what used to be called launch_def
        #launch dfs and get unclean list
        #call some function that sweeps the list (function can be outside)
        #returns cleaned-up list of loops
        
        
        for node in range (self.network.number_of_clusters):
            if self.visited_nodes[node] == -1:           #should this be made into a method like is_node_visited (self, node) > return True/False or is this an overkill?
                print ("new origin of network: ", node)
                #reset stuff to 0
                self.crossing_sum = [[0,0,0]] #every time I start visiting a new network (that has never been visited before), I set the sum of crossings to 0
                #self.current_crossing = 0
                self.loops_temp = []
                self.discovery_time_node = 0
                self.discovery_time_edge = 0
                LoopFinder.dfs(self, node)
                
                self.discovery_time_node += 1
                if self.loops_temp: 
                
                    self.loops_list.extend(self.loops_temp)
                    
        return self.loops_list

# remove identical rows > should maybe be done somewhere in the code as well
def rows_uniq_elems(a):
    new_array = [tuple(row) for row in a]
    uniques = np.unique(new_array, axis = 0)
    return uniques

 #loops info should also provide a color for each loop to denote connected components
# lump before clean would discard that info
# clean before lump would allow having it

def linearly_independent (loops_list):
    loops_list = np.asarray(loops_list)
    print("loops list: ", loops_list)
    print("rank according to numpy.linalg: {}".format(np.linalg.matrix_rank(loops_list,tol=1e-8)))
    _, inds = sympy.Matrix(loops_list).T.rref()
    Nloops = len(inds)
    
    independent_loops = loops_list[list(inds)] 
    #print ("independent_loops index and loops: ", independent_loops, inds)
    return independent_loops, Nloops
  
number_of_clusters = np.amax(dropped_list) +  1
## NOTE: at this step, in my previous code, I have a dropped_list that only contains boundary crossing elements
my_test_network = PeriodicNetwork(number_of_clusters, number_of_clusters)  #fix this so you don't take boundary crossings i to account; should be fixed once we turn rev list into object


print(dropped_list)
dropped_list = rows_uniq_elems(dropped_list)
print(dropped_list)
for i in range(len(dropped_list)):
    edge_info = dropped_list[i,:]
    print(edge_info)
    
    my_test_network.add_edge(edge_info[0], edge_info[1], edge_info[2:])
    
print ("neighbors", my_test_network.neighbors) 
print ("boundary crossing ", my_test_network.bond_is_across_boundary)  
for i in range (len (my_test_network.neighbors)):
    print ("number of neighbours: ", my_test_network.get_number_of_neighbors(i))
    for n_index in range(my_test_network.get_number_of_neighbors(i)):
                neigh = my_test_network.get_neighbor(i,n_index)
                edge = my_test_network.get_edge(i,n_index)
                print ("edge", edge)
                

if not my_test_network.crosses_boundaries:
    print ("Network does not cross periodic boundaries in any direction, therefore it does not percolate")
    exit()

#the following lines could already be part of get_reduced_network
if my_test_network.needs_reducing:
    print ("Reducing network to simpler form...")
    my_reduced_network = my_test_network.get_reduced_network() #think of a way to return it already in the right format you would get with add_edges
    myloops=LoopFinder(my_reduced_network)
else:          
    myloops=LoopFinder(my_test_network)
loops=myloops.get_loops()


print ("number of loops \n" , len(linearly_independent(loops)[0]), " \n independent loops: \n", linearly_independent(loops)[0])


