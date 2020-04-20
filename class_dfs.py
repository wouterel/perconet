#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 18:04:33 2020

@author: craffaelli
"""
import numpy as np
from dropped_list_test import dropped_list  


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
        
        
    def add_edge(self,node1,node2,boundary_vector):
        """ Adds info for this edge """
        self.neighbors[node1,self.neighbors_counter[node1]] = node2
        self.edges_list [node1,self.neighbors_counter[node1]] = self.edges_counter  #whatis a smart way to count "edge_count"?
        self.boundary_crossing [node1,self.neighbors_counter[node1], :] = boundary_vector
        self.neighbors_counter [node1] += 1
        
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
        #returns array of neighbor indices of node i (could be length maxN, padded with -1s, or length actual number)
        neighbors_of_i = self.neighbors[i, :]
        if padded:
            return neighbors_of_i
        if not padded:
            stripped_neighbors_of_i = neighbors_of_i [neighbors_of_i != -1]
            return stripped_neighbors_of_i
             
    def get_neighbor(self,i,n_index):
        #returns n_index'th neighbor of node i
        return self.neighbors[i, n_index]
    
    def get_edges(self,i):
        #returns edge corresponding to 
        return self.edges_list [i, :]
    
    def get_edge(self,node1,node2):
        #returns edge between two nodes  ##check if edge exists!
        return self.edges_list [node1,node2]
    
    def get_boundary_crossing(self,i,n_index):
        #returns the bc vector of the n_index'th neighbor of node i
        return self.boundary_crossing [i,n_index, :]
                
number_of_clusters = np.amax(dropped_list) +  1
## NOTE: at this step, in my previous code, I have a dropped_list that only contains boundary crossing elements
my_test_network = PeriodicNetwork(number_of_clusters, number_of_clusters)  #fix this so you don't take boundary crossings i to account; should be fixed once we turn rev list into object

for i in range(len(dropped_list)):
    edge_info = dropped_list[i,:]
    my_test_network.add_edge(edge_info[0], edge_info[1], edge_info[2:])
    
    

class LoopFinder:
    def __init__(self,network):
        self.network=network
        #allocate arrays and values for dfs stuff
        self.start = 0
        self.discovery_time_node = 0
        self.discovery_time_edge = 0
        self.visited_nodes = -1 * np.ones (network.number_of_clusters, dtype = int)
        self.visited_edges = -1 * np.ones (np.amax(network.edges_list) + 1, dtype = int)
        #self.current_crossing = [0,0,0]
        self.crossing_sum = [[0,0,0]]
        self.loops_temp = []
        self.loops_list = []
        
    def dfs(self):
        #in here, get network info whenever using self.network.get_blabla()
        #this function can be private (not callable from the outside)
        self.visited_nodes [self.start] = self.discovery_time_node
        
        print ("crossing sum now: \n", self.crossing_sum , "\n visited nodes: ", self.visited_nodes, "\n visited edges: \n", self.visited_edges)
        
        for n_index in range(self.network.get_number_of_neighbors(self.start)):
            neigh = self.network.get_neighbor(self.start,n_index)
            edge = self.network.get_edge(self.start,n_index)
            
            if edge==-1:
                break
            
            if self.visited_edges[edge] == -1: # otherwhise neigh is a parent of start
                self.visited_edges[edge] = self.discovery_time_edge
                ("visited edges now: ", self.visited_edges)
                self.discovery_time_edge += 1
                timestep = self.visited_nodes [self.start]
                
                current_crossing = self.crossing_sum [timestep] + self.network.get_boundary_crossing(self.start, n_index)
                xC,yC,zC = current_crossing[:]
                #print ("current crossing sum for ",start,neigh,  ' (index: ', n_index, "), is ", xC,yC,zC, ", current crossing is ",boundary_crossing[start, n_index] )
                if self.visited_nodes[neigh] == -1:
                    self.discovery_time_node += 1
                    self.visited_nodes[neigh] = self.discovery_time_node
                    self.crossing_sum.append ([xC,yC,zC])
                    
                    #print ("visited nodes: \n",visited_nodes, "\n crossing sum \n", crossing_sum)
                    #loops_temp, discovery_timeC, discovery_timeE = dfs (neigh, discovery_timeC, discovery_timeE, neighbors, visited_nodes, current_crossing, crossing_sum, boundary_crossing, edges_list, visited_edges, loops_temp, Nclusters)
                    self.loops_temp, self.discovery_time_edge, self.discovery_time_node = myloops.dfs()
                    #print ("crossing sum now: \n", crossing_sum, "loops temp is:", loops_temp)
                else: #we found a loop!
                    
                    loop_timestep = self.visited_nodes [neigh]
                    
                    loop = current_crossing - self.crossing_sum [loop_timestep]

                    self.loops_temp.append (loop.tolist())
                    current_crossing -= self.boundary_crossing[self.start, n_index]
                    #print("loop!  -> ", loops_temp, " loop_timestep: ", loop_timestep, ", crossing sum: ", crossing_sum )
        
        
        
        return self.loops_temp, self.discovery_time_edge, self.discovery_time_node
    
    
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
                self.loops_temp, self.discovery_time_edge, self.discovery_time_node = myloops.dfs()
                
                self.discovery_time_node += 1
                if self.loops_temp: 
                
                    self.loops_list.extend(self.loops_temp)
                    
        return self.loops_list
    
 #loops info should also provide a color for each loop to denote connected components
# lump before clean would discard that info
# clean before lump would allow having it
  

        
myloops=LoopFinder(my_test_network)
loops=myloops.get_loops()

print (loops)