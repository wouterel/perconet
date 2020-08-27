#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 10:23:55 2020

@author: wouterel
"""
import numpy as np
   

class testrecursive:
    def __init__(self,n):
        self.factor=n
        self.result=1
    def recurse(self):
        self.result=self.result*self.factor
        print(self.result)
        if (self.result<10000):
            self.recurse()
            
        
class PeriodicNetwork:
    #stores network structure
    #contains methods like get_something() to provide network info
    def __init__(self,n,maximum_neighbors_per_node):
        #allocate arrays for per-node info
        
        #(substitute function neighbour_clusters in clst)
        self.boundary_crossing = np.zeros((Nclusters,max_links,3), dtype = int)
        self.neighbors = -1 * np.ones((n,maximum_neighbors_per_node), dtype = int)
        self.edges_list = -1 * np.ones((n,maximum_neighbors_per_node), dtype = int)

        
    def add_edge(self,node1,node2,boundary_vector):
        #add info for this edge

    def get_number_of_neighbors(self,i):
        #returns number of neighbors of node i
        
    def get_neighbors(self,i):
        #returns array of neighbor indices of node i (could be length maxN, padded with -1s, or length actual number)
        
    def get_neighbor(self,i,n_index):
        #returns n_index'th neighbor of node i
        
    def get_edges(self,i):
        #returns array of edge indices 

    def get_boundary_crossing(self,i,n_index):
        #returns the bc vector of the n_index'th neighbor of node i
        
        
class LoopFinder:
    def __init__(self,network):
        self.network=network
        #allocate arrays for dfs stuff
        
    def dfs(self):
        #in here, get network info whenever using self.network.get_blabla()
        #this function can be private (not callable from the outside)
        
    def get_loops():
        #this is what used to be called launch_def
        #launch dfs and get unclean list
        #call some function that sweeps the list (function can be outside)
        #returns cleaned-up list of loops
        
 #loops info should also provide a color for each loop to denote connected components
# lump before clean would discard that info
# clean before lump would allow having it

#testnetwork: import testnetwork from dropped_list_test
from dropped_list_test import dropped_list
testnetwork =         

mynetwork=PeriodicNetwork(10,9)
#now add edges using mynetwork.add_edge(3,5,[1,1,0]) etc (this would happen in the lammps reader code)


myloops=LoopFinder(mynetwork)
loops=myloops.get_loops()
#myloops.dfs(bla) would give error because it's private?




        
        
#one class to store network topology + boundary crossing vectors
#one class to do the depth-first search 
        
        
test=testrecursive(2)
test.recurse()