#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 11:41:35 2018

@author: craffaelli
"""

 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 19:01:44 2018

@author: craffaelli
"""

from lammpstools import dump_reader
from lammpstools import block_data
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import copy
import sys
import os



import import_files as imp
import import_snapshot as imp0
import cluster_find3d as clst
import DFS_cycle_finder as dfs
import reversible_bonds as rev
import import_reversible as impRev

#general variables, paths to files; 
path = "/Users/craffaelli/projects/phd/percolation/code/test/C810Cr70S17/"
dumpname = path+"centerTraj_c8_s10_cRev70_s117eq.bin" #dump file with
print("dump", dumpname) 

revdumpname = path + "revEndTraj_c8_s10_cRev70_s117eq.bin"
snapshot="final_snapshot_c8_s10_cRev70_s117eq.txt"
snapshotname = path + snapshot
revBonds_count = "revBonds_c8_s10_cRev70_s117eq_fixed.txt"
revcountname = path + revBonds_count

name = 'test'
print("name",name)

#set the types that are used in LAMMPS
irreversibleBondType = 2
irreversibleBeadType = 4
armsN=4
revBondsCutoff = 0.5

# using Stefan's lammpstools
d = dump_reader.dump_reader( dumpname,
                            file_format = "BIN", dump_format = "LAMMPS" )

dr = dump_reader.dump_reader( revdumpname,
                            file_format = "BIN", dump_format = "LAMMPS" )
#di = dump_reader.dump_reader( irrdumpname,
#                            file_format = "BIN", dump_format = "LAMMPS" )
#
db = dump_reader.dump_reader( revcountname,
                              file_format = "PLAIN", dump_format = "LAMMPS",
                              is_local = True )

#d.no_block_data_copy = False

# You _have_ to tell it upfront which headers to expect in the dump file.
# This is because it needs those to map them to the internal arrays
# ids, mol, type and x. They should match your dump custom command
# but don't have to be in order. Also, you can omit "mol" if you're not
# interested in that.
d.set_column_headers( [ 'id', 'mol', 'type', 'x', 'y', 'z'] )
dr.set_column_headers( [ 'id', 'mol', 'type', 'x', 'y', 'z'] )

## Stefan: I think maybe the number of headers you pass has to match the                                                                                                      
## number in the file, not exactly sure.
id0='id'
isbound='c_revBonds'
db.set_column_headers( [id0, isbound] )

#set percolation as False by default - this will be changed to True if the network of irreversible bonds percolates
percolation = False


# extract total number of beads
with open(snapshotname) as f:
    imp.skip_lines(2, f)
    
    Natoms = imp.meta(1, int,f)
    
    imp.skip_lines(4, f)
    
    xlim = imp.meta(2, float ,f)
    ylim = imp.meta(2, float , f)
    zlim = imp.meta(2, float ,f)
    boxsize = [abs(xlim[1] - xlim[0]), abs(ylim[1] - ylim[0]), abs(zlim[1] - zlim[0])]
    
    # extract data about bonds
    imp.find_keyword("Atoms", f)
    #imp0.snapshot_maptest(irreversibleBeadType, Natoms, f)
    map_ends_to_star = imp0.snapshot_map(irreversibleBeadType, Natoms, f)
    print("map: \n", map_ends_to_star)
   
    imp.find_keyword("Bonds", f)
    irrev_bonds, Nbonds = imp0.snapshot_mapBonds(irreversibleBondType, f)
    #print("bonds number: ", Nbonds)
    

# make a dictionary that connects endbeads to centers via their molecule number; 
#    we use the first timestep of the reversible bonds dump file
br_0 = dr.next_block()
NbeadsPart = br_0.meta.N
map_RevEnds_to_star = imp.create_map(Natoms,NbeadsPart,br_0)

b_0 = d.next_block()
Ncenters = b_0.meta.N
Nstars= Ncenters
max_crosslinks = (Nstars * armsN) / float(2)
time_now = b_0.meta.t



# **** check if the system percolates if we only take irreversible bonds into account **** 

coordinates = imp.map_coords(Nstars, b_0)

#print (coordinates)
crosslinks, crosslink_boundaries = imp.store_crosslinks1(Nstars, Nbonds, map_ends_to_star, irrev_bonds, coordinates, boxsize)
#print ("crosslinks: \n", crosslinks, "\n crosslink boudaries: \n", crosslink_boundaries)
crosslink_tmp, dropped_list = clst.drop (crosslinks, crosslink_boundaries)
print (dropped_list)
if dropped_list.any(): # if false, no percolation in any direction for sure   
    neighbours = clst.store_neighbours(Nstars, armsN, crosslink_tmp)
    print ("neighbour list",neighbours)
    list_colors, Ncolors = clst.cluster_find(neighbours, Nstars)
    #print ("list colors: \n",list_colors, len(list_colors))
    #print ("dropped list: \n", dropped_list, len(dropped_list))
    dropped_list = clst.molid_to_clusterid (list_colors, dropped_list)
    #print ("dropped list: \n", dropped_list)
    Cneigh, boundary_crossing, edges_list = clst.neighbour_clusters (dropped_list, Ncolors)
    #print ("boundary crossing \n",boundary_crossing,"\n \n edges list \n", edges_list, "\n \n C_neigh and Ncolors \n ",Cneigh, Ncolors)
    #print ("Cneigh: ", Cneigh)
    #print ("n colours", Ncolors)
    independent_loops, Nloops = dfs.launch_dfs (Cneigh, Ncolors, boundary_crossing, edges_list)
    #print ("independent loops: ",independent_loops)
    #print ("network percolates in ", Nloops, " independent direction(s) -> ", independent_loops)
    percolation_directions = Nloops
    if Nloops == 3:
        percolation = True
    #print ("system percolates in ", Nloops, " independent direction(s)")

counter = 0
#print(time_now)


    

    
c = 0
timestep = 700010000

if percolation == False:    # if system does not percolate with irreversible crosslinks, then check percolation including reversible crosslinks
    
    map_allEnds_to_star = map_RevEnds_to_star + map_ends_to_star
    percolation_time = []
    while d and dr and db:
        # synchronise the dump files
        blocks = imp.synchro(d, dr, db)
        b = blocks[0]
        br = blocks[1]
        bb = blocks[2]
        if b is None or br is None or bb is None:
            break       
                # extract general information
        time_now = b.meta.t
        #print (time_now)
        c+= 1
        #if time_now == timestep:
        if (c%2500 == 0):
                # do stuff
                
        # 0. get reversible bonds 
            NbondsMax = bb.meta.N
            reversible_bonds, bound = impRev.import_rev (bb, NbondsMax, id0, isbound)
            revBeads_ixMax = max(bb.data_by_name(id0))
            revEnd_coordinates = impRev.map_coords (revBeads_ixMax,NbondsMax, br)
            
            #print (reversible_bonds)
            bond_ratio = bound/NbondsMax
            revList = rev.make_revbondlist (reversible_bonds, bound, revEnd_coordinates, boxsize, revBondsCutoff ** 2)
            #print ('reversible bonds list: ', revList)
            
            # join REVERSBLE and IRREVERSIBLE bonds lists
            allBondsList = np.concatenate((revList, irrev_bonds))
            NbondsTot = len(allBondsList)
            # 1. make map of 
            coordinates = imp.map_coords(Nstars, b)
            
            crosslinks, crosslink_boundaries = imp.store_crosslinks1(Nstars, NbondsTot, map_allEnds_to_star, allBondsList, coordinates, boxsize)
            #print ("crosslinks: \n", crosslinks, "\n crosslink boudaries: \n", crosslink_boundaries)
            #print ("crosslinks: \n", crosslinks)
            
            crosslink_tmp, dropped_list = clst.drop (crosslinks, crosslink_boundaries)
            #print (dropped_list)
            if dropped_list.any(): # no percolation in any direction for sure   
                neighbours = clst.store_neighbours(Nstars, armsN, crosslink_tmp)
                #print ("neighbour list",neighbours)
                list_colors, Ncolors = clst.cluster_find(neighbours, Nstars)
                #print ("list colors: \n",list_colors, len(list_colors))
                #print ("dropped list: \n", dropped_list, len(dropped_list))
                dropped_list = clst.molid_to_clusterid (list_colors, dropped_list)
                #print ("dropped list: \n", dropped_list)
                Cneigh, boundary_crossing, edges_list = clst.neighbour_clusters (dropped_list, Ncolors)
                #print ("boundary crossing \n",boundary_crossing,"\n \n edges list \n", edges_list, "\n \n C_neigh and Ncolors \n ",Cneigh, Ncolors)
                #print ("Cneigh: ", Cneigh)
                #print ("n colours", Ncolors)
                independent_loops, Nloops = dfs.launch_dfs (Cneigh, Ncolors, boundary_crossing, edges_list)
                #print ("independent loops: ",independent_loops)
                #print ("network percolates in ", Nloops, " independent direction(s) -> ", independent_loops)
                percolation_directions = Nloops
                if Nloops == 3:
                    percolation = True
                print ("system percolates in ", Nloops, " independent direction(s)")
                percolation_time.append ([time_now, Nloops])
            
            
            
            
            
            
            #break

            
    print ("percolation time: ", percolation_time)    
        
