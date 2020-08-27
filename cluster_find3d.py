#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 16:06:12 2018

@author: craffaelli
"""

# import libraries:
import numpy as np
#import pandas as pd



# counts how many redundant bonds (two stars connected via more than 1 crosslink) there are
def count_redundant_bonds(crosslinks, armsN):
    """ Counts how many redundant bonds (two stars connected via more than 1 crosslink) there are. """
    unique, counts = np.unique(crosslinks, return_counts=True)
    N_multiples = dict(zip(unique, counts))
    useless_info =[0,1] #delete info for the recurrence of 0 and 1
    for e in useless_info:
        N_multiples.pop(e, None)
    # make 1 entry of the array: N_double links, N_triple links etc 
    multiple_crosslinks = np.zeros(armsN-1, dtype= int)
    for i in range (len(N_multiples)):
        multiple_crosslinks[i] = list(N_multiples.keys())[i]*list(N_multiples.values())[i] 
    return multiple_crosslinks

       
       # make array of plots 
       
# cluster finding algorithm

def store_neighbours(N_stars,armsN, crosslinks):
    counter = np.zeros(N_stars,dtype = int)
    #initialise array of nearest neighbours: -1 == empty
    neighbours = -1* np.ones((N_stars,armsN), dtype= int)
    #loop over non zero elements of the array, make list of bonded stuff
    for index,x in np.ndenumerate(crosslinks):  
        #index = (first bonded star, second bonded star), x (how many bonds)
        if x != 0:
            ind1,ind2 = index[0], index[1]
            neighbours[ind1,counter[ind1]] = ind2
            counter[ind1] += 1
            neighbours[ind2,counter[ind2]] = ind1
            counter[ind2] += 1
    return neighbours


def coloring(start, current_color, neighbours, color):
    color[start] = current_color
    for index in range(len(neighbours[1])):
        neigh = neighbours[start, index]
        if neigh != -1 and color[neigh] == -1: #empty: not connected 
            color[neigh] = current_color
            coloring(neigh, current_color, neighbours, color)        
        

def cluster_find(neighbours, N_stars):  
    #initiate with first cluster colour
    current_color = 0
    # initialise list of star colours (-1 == not coloured)
    color = -1* np.ones(N_stars, dtype = int)
    for star in range(N_stars):
        if color[star] == -1:
            start_cluster = star
            coloring (start_cluster, current_color, neighbours, color)
            current_color += 1
    Ncolors = np.amax(color) + 1
    return color, Ncolors

def cluster_distribution(list_colors, N_stars):
    #count size of each cluster
    unique, counts = np.unique(list_colors, return_counts=True)
    cluster_sizes = dict(zip(unique, counts))
    # make list of cluster sizes 
    size_list = []
    for size in cluster_sizes.values():
        size_list.append(size)
    # count occurrence of each size (is a dict {"cluster size": 'occurrence'}
    unique, counts = np.unique(size_list, return_counts=True)
    size_distr_dict = dict(zip(unique, counts))

    # store size distribution in 1D array
    size_distr = np.zeros(N_stars, dtype =int)
    idx = [int(key)-1 for key in size_distr_dict]
    val = [v for v in size_distr_dict.values()]
    size_distr[idx] = val
    return size_distr


def drop (crosslinks, crosslink_boundaries):
    """ removes all crosslinks across the boundaries and stores them in another list for later use"""
    #dropped = abs(crosslink_boundaries).sum (axis = 2)
    crosslinks_tmp = crosslinks
    #print ("dropped, crosslinks, tmp \n", dropped,"\n", crosslinks, "\n",crosslinks_tmp )
    dropped_list = []
    for i in range(len(crosslink_boundaries)):
        for j in range(len(crosslink_boundaries)):
            val = crosslink_boundaries[i,j,:]
            if np.any (val):
                dropped_list.append([i,j,val[0], val[1], val [2]])
                crosslinks_tmp [i,j] = 0
    dropped_list = np.asarray(dropped_list)

    #gets rid of information on multiply bonded stars:
    if dropped_list.any():
        bondsIDs = dropped_list[:, :2]
        boundary_crossing = np.sign(dropped_list[:,2:])
        dropped_list = np.concatenate((bondsIDs,boundary_crossing), axis = 1 )
    return crosslinks_tmp, dropped_list

    
def connected_clusters (dropped_list, list_colors, direction):
    """ Checks if clusters are connected to themselves through at least one linker. """
    percolation = False
    Nc = np.amax(list_colors) + 1
    #print (Nc, list_colors)
    Nc_matrix = np.zeros((Nc,Nc), dtype = int)
    for i in range(len(dropped_list)):
        item1 = dropped_list[i,0]
        item2 = dropped_list[i,1]
        
        color1 = list_colors[item1]
        color2 = list_colors[item2]
        
        Nc_matrix [color1, color2], Nc_matrix [color2, color1] = dropped_list[i, 2], -dropped_list[i, 2]
        if color1 == color2:
            percolation = True
            print ("system is percolating in direction ", direction, ", cluster ", color1 , ", 1st round")
    return percolation, Nc_matrix


def molid_to_clusterid (list_colors, dropped_list):
    for i in range(len(dropped_list)):
        item1 = dropped_list[i,0]
        item2 = dropped_list[i,1]
        
        color1 = list_colors[item1]
        color2 = list_colors[item2]
        dropped_list[i,:2] = color1, color2
    dropped_list = np.unique(dropped_list, axis=0)
    return dropped_list
        
def neighbour_clusters (dropped_list, Nclusters):
    edge_count = 0
    counter = np.zeros(Nclusters,dtype = int)
    #print ("counter, Nclusters: ", counter, Nclusters)
    max_links = len(dropped_list) * 2
    #print ("total number of linkers across boundaries: ", max_links)
    Cneighbors = -1 * np.ones((Nclusters,max_links), dtype = int)
    edges_list = -1 * np.ones((Nclusters,max_links), dtype = int)
    #print (Cneighbors)
    Ccrossings = np.zeros((Nclusters,max_links,3), dtype = int)
    
    for link in range(len(dropped_list)):
        ind1,ind2 = dropped_list[link, 0], dropped_list[link, 1]
        #print ("index1 , index2, counter, counter[ix1]: ", ind1,ind2, counter,counter[ind1] )
        Cneighbors[ind1,counter[ind1]] = ind2
        edges_list [ind1,counter[ind1]] = edge_count
        Ccrossings[ind1,counter[ind1], :] = dropped_list [link,2:]
        counter[ind1] += 1
        Cneighbors[ind2,counter[ind2]] = ind1
        edges_list [ind2,counter[ind2]] = edge_count
        Ccrossings[ind2,counter[ind2], :] = -1* dropped_list [link,2:]
        counter[ind2] += 1
        
        edge_count += 1
        # print ("edges list: \n", edges_list, "\n neighbours: \n", Cneighbors )
    return Cneighbors, Ccrossings, edges_list
    


def loose_ends (Nc_matrix, n,x):
    proceed = True
    # check that matrix is not empty
    #print("matrix sum = ", np.sum(np.absolute(Nc_matrix)))
    if np.sum(np.absolute(Nc_matrix)) != 0:
        valid=np.where(np.any(Nc_matrix,axis=0))[0]
        Nc_matrix = Nc_matrix[np.ix_(valid,valid)]
        Nc = len (Nc_matrix)
        deleted = 0
        for i in range(Nc):
            count = 0
            for j in range(Nc):
                if Nc_matrix[i,j] != 0:
                    count +=1
                    jstar= j
            #print ("count at stage ", n, "is ", count )
            if count == 1:
                Nc_matrix[i,jstar] = 0
                Nc_matrix [jstar,i] = 0
                deleted += 1
        n+=1
        #print("Nc_matrix at round", n,"in direction ", x, "\n",  pd.DataFrame(Nc_matrix))
        
        if deleted != 0:
            loose_ends(Nc_matrix, n, x)
    else:
        #print ("system does NOT percolate in direction ", x)
        proceed = False
    return Nc_matrix, proceed

# valid=np.where(np.any(Nc_matrix,axis=0))[0]
# NC2=Nc_matrix[np.ix_(valid,valid)]
#
#def burning (Nr_matrix, x):
#    Nr = len(Nr_matrix)
#    q_a = np.zeros((Nr,2), dtype = int)
#    t = 0
#    #q_a [0,0], q_a[0,1] = 0, 1 
#    
#    for b in range(Nr):
#        if q_a [b, 1] == 0:
#            q_a [b, 1] = 1  # assign value and 1 because it is burning
#            # run through neighbours 
#            for a in range (Nr):
#                burning_list = []
#                S_ab = Nr_matrix [a, b]
#                if S_ab != 0:
#                    burning_list.append(a)
#                    q_a [a, 0] = 
    
    
