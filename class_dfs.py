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
        return
                
number_of_clusters = np.amax(dropped_list) +  1
## NOTE: at this step, in my previous code, I have a dropped_list that only contains boundary crossing elements
my_test_network = PeriodicNetwork(number_of_clusters, number_of_clusters)  #fix this so you don't take boundary crossings i to account; should be fixed once we turn rev list into object

for i in range(len(dropped_list)):
    edge_info = dropped_list[i,:]
    my_test_network.add_edge(edge_info[0], edge_info[1], edge_info[2:])
    
    

class LoopFinder:
    def __init__(self,network):
        self.network=network
        #allocate arrays for dfs stuff
        
    def dfs(self):
        #in here, get network info whenever using self.network.get_blabla()
        #this function can be private (not callable from the outside)
        return
    def get_loops():
        #this is what used to be called launch_def
        #launch dfs and get unclean list
        #call some function that sweeps the list (function can be outside)
        #returns cleaned-up list of loops
        return
 #loops info should also provide a color for each loop to denote connected components
# lump before clean would discard that info
# clean before lump would allow having it

        
        
myloops=LoopFinder(mynetwork)
loops=myloops.get_loops()

