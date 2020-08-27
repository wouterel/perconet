#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 10:33:40 2018

@author: craffaelli
"""

import numpy as np
import sympy 
#import pandas as pd


def dfs (start, discovery_timeC, discovery_timeE, neighbors, visited_nodes, current_crossing, crossing_sum, boundary_crossing, edges_list, visited_edges, loops_temp, Nclusters) :
    visited_nodes[ start ] = discovery_timeC
    print ("crossing sum now: \n", crossing_sum , "\n visited nodes: ", visited_nodes, "\n visited edges: \n", visited_edges)
    #print("neighbours and edges list are: \n ", neighbors,"\n", edges_list)
    #loops_temp = []
    #print ("len neighbour (0) is ",(len(neighbors[0])), neighbors, "N clusters is ", Nclusters )
    

    #Even better:
    #for n_index in network.get_number_of_neighbors(start):
    #   neigh=network.get_neighbor(start,n_index)
    #   edge=network.get_edge(start,n_index)
    #   if edge==-1:
    #       break
    #   bc=network.get_bc(start,n_index) #could put this further down to avoid unnecessary calls

    for n_index in range (len(neighbors[0])): #neighbors
        neigh = neighbors[start, n_index]
        edge = edges_list[start, n_index]

        print( "start and neigh are: ", start, neigh, "edge ",  edge, "index ",  n_index, "visited_edges[edge] ", visited_edges[edge])
        if edge != -1 : #if edge is -1 it's a loose end (equivalent to writing (ïf neigh != -1))
            
            if visited_edges[edge] == -1: # otherwhise neigh is a parent of start
                visited_edges[edge] = discovery_timeE
                ("visited edges now: ", visited_edges)
                discovery_timeE += 1
                timestep = visited_nodes [start]
                
                current_crossing = crossing_sum [timestep] + boundary_crossing[start, n_index]
                xC,yC,zC = current_crossing[:]
                #print ("current crossing sum for ",start,neigh,  ' (index: ', n_index, "), is ", xC,yC,zC, ", current crossing is ",boundary_crossing[start, n_index] )
                if visited_nodes[neigh] == -1:
                    discovery_timeC += 1
                    visited_nodes[neigh] = discovery_timeC
                    crossing_sum.append ([xC,yC,zC])
                    
                    #print ("visited nodes: \n",visited_nodes, "\n crossing sum \n", crossing_sum)
                    loops_temp, discovery_timeC, discovery_timeE = dfs (neigh, discovery_timeC, discovery_timeE, neighbors, visited_nodes, current_crossing, crossing_sum, boundary_crossing, edges_list, visited_edges, loops_temp, Nclusters)
                    #print ("crossing sum now: \n", crossing_sum, "loops temp is:", loops_temp)
                else: #we found a loop!
                    
                    loop_timestep = visited_nodes [neigh]
                    
                    loop = current_crossing - crossing_sum [loop_timestep]

                    loops_temp.append (loop.tolist())
                    current_crossing -= boundary_crossing[start, n_index]
                    #print("loop!  -> ", loops_temp, " loop_timestep: ", loop_timestep, ", crossing sum: ", crossing_sum )
            
        else:
            print ("NO more clusters bound to cluster ", start)
            break
            
            #print ("visited edges: \n", visited_edges)
            #print("loops now  ", loops_temp)
            #print (" visited nodes at the end: ", visited_nodes, "crossing sum now: \n", crossing_sum)
    return loops_temp , discovery_timeC, discovery_timeE
                
def linearly_independent (loops_list):
    loops_list = np.asarray(loops_list)
    print("loops list: ", loops_list)
    _, inds = sympy.Matrix(loops_list).T.rref()
    Nloops = len(inds)
    
    independent_loops = loops_list[list(inds)] 
    #print ("independent_loops index and loops: ", independent_loops, inds)
    return independent_loops, Nloops


def launch_dfs (neighbors, Nclusters, boundary_crossing, edges_list):
    
    loops_list = []
    visited_nodes = -1 * np.ones (Nclusters, dtype = int)
    visited_edges = -1 * np.ones (np.amax(edges_list) + 1, dtype = int)
    #discovery_timeC = 0
    #discovery_timeE = 0
    for node in range(Nclusters):
        #print ("node and visited node tag: ",node, visited_nodes[node])
        if visited_nodes[node] == -1:
            print ("new origin of network: ", node)
            crossing_sum = [[0,0,0]] #every time I start visiting a new network (that has never been visited before), I set the sum of crossings to 0
            current_crossing = 0 # I also set the current crossing to zero
            loops_temp = [] # and also the temporary list of loops found in the current network is set to 0
            discovery_timeC = 0
            discovery_timeE = 0
            loops, discovery_timeC, discovery_timeE = dfs (node, discovery_timeC, discovery_timeE, neighbors, visited_nodes, current_crossing, crossing_sum, boundary_crossing, edges_list, visited_edges, loops_temp, Nclusters)
            discovery_timeC += 1
            print('loops is ', loops)
            if loops: 
                
                loops_list.extend(loops)
                #print ('current loops found are: ', loops_list)
                
    
    #loops_list = np.asarray(loops_list)
    independent_loops, Nloops = linearly_independent (loops_list)
    #print (" visited nodes at the end: ", visited_nodes, "crossing sum now: \n", crossing_sum)
    return independent_loops, Nloops
