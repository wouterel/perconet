#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 09:19:30 2020

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

#import class_dfs
#from class_dfs import PeriodicNetwork
#from class_dfs import LoopFinder

#general variables, paths to files; 
from localsettings import path

def is_file_old(list_files):
    for path_item in list_files:
        check = os.path.isfile(path_item) 
        if not check:
            print(path_item, " ERROR! file cannot be found!" ) 
            sys.exit()
    return(check)

def is_file(list_files):
    all_files_exist = True
    for path_item in list_files:
        check = os.path.isfile(path_item) 
        if not check:
            print(path_item, " ERROR! file cannot be found!" ) 
            #sys.exit()
            all_files_exist = False
    return(all_files_exist)            
        
def is_path(list_paths):
    for path_item in list_paths:
        check = os.path.isdir(path_item) 
        if not check:
            print(path_item, " ERROR! folder does NOT exist!" ) 
            sys.exit()

dumpname = os.path.join(path,"centerTraj_c8_s10_cRev70_s117eq.bin") #dump file with
print("dump", dumpname) 

revdumpname = os.path.join(path,"revEndTraj_c8_s10_cRev70_s117eq.bin")
snapshot="final_snapshot_c8_s10_cRev70_s117eq.txt"
snapshotname = os.path.join(path,snapshot)
revBonds_count = "revBonds_c8_s10_cRev70_s117eq_fixed.txt"
revcountname = os.path.join(path,revBonds_count)

name = 'test'

#set the types that are used in LAMMPS  /// how do we make this more general?
irreversibleBondType = 2
irreversibleBeadType = 4
armsN=4
revBondsCutoff = 0.5

# using Stefan's lammpstools
d = dump_reader.dump_reader( dumpname,
                            file_format = "BIN", dump_format = "LAMMPS" )

dr = dump_reader.dump_reader( revdumpname,
                            file_format = "BIN", dump_format = "LAMMPS" )

db = dump_reader.dump_reader( revcountname,
                              file_format = "PLAIN", dump_format = "LAMMPS",
                              is_local = True )

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

#lammpsreader class (import_files)


#make class for system and snapshot ///should only 
class SystemBox:
    """ Stores the data that is unchanged through the simulation"""
    def __init__(self):
        #define details
        #self.file = file
        self.xlim = []
        self.ylim = []
        self.zlim = []
        self.boxsize = []
        
        self.Natoms = 0
        self.Nbonds = 0
        
        self.NbeadsPart = 0
        self.Ncenters = 0
        self.max_crosslinks = 0
        
        self.map_ends_to_star = []
        self.map_bonds = []
        self.map_RevEnds_to_star = []
        
        
    def open_files(self, filename, d, dr):
        #is_file = os.path.isfile(filename) 
        if is_file([filename]): ## shuld find a way for this to work with all files
            with open(filename) as f:
                self.file=f
                #here you do all the things that have to be done in a certain order
                imp.skip_lines(2, f)
                self.Natoms = imp.meta(1, int,f)
        
                imp.skip_lines(4, f)
        
        
                self.xlim = imp.meta(2, float ,f)
                self.ylim = imp.meta(2, float , f)
                self.zlim = imp.meta(2, float ,f)
        
                self.boxsize = [abs(self.xlim[1] - self.xlim[0]), abs(self.ylim[1] - self.ylim[0]), abs(self.zlim[1] - self.zlim[0])]
                
                # extract data about bonds
                imp.find_keyword("Atoms", f)
                #imp0.snapshot_maptest(irreversibleBeadType, Natoms, f)
                self.make_map_ends_to_star()
                print("map: \n", self.map_ends_to_star)
               
                imp.find_keyword("Bonds", f)
                self.make_map_bonds()
                #irrev_bonds, Nbonds = imp0.snapshot_mapBonds(irreversibleBondType, f)
                print("bonds number: ", self.Nbonds)   
                
                # get info from other files
                br_0 = dr.next_block()
                self.NbeadsPart = br_0.meta.N
                self.map_RevEnds_to_star = imp.create_map(self.Natoms,self.NbeadsPart,br_0)
                
                b_0 = d.next_block()
                self.Nstars = b_0.meta.N
                self.max_crosslinks = (self.Nstars * armsN) / float(2)
                
        else:
            print(filename, " ERROR! file cannot be found!" ) 
        
        return is_file

    def get_boxsize(self):
        return self.boxsize
    
    def get():
        return
        
    def make_map_ends_to_star(self):
        self.map_ends_to_star = imp0.snapshot_map(irreversibleBeadType, self.Natoms, self.file) #actually import function here
        return
    
    def make_map_bonds(self):
        self.irrev_bonds, self.Nbonds = imp0.snapshot_mapBonds(irreversibleBondType, self.file)
    
    def make_map_RevEnds_to_star(self):
        return
      
class Snapshot:
    """ stores information that is typical of one timeframe"""
    def __init__ (self, t, contains_reversible = False):
        self.t = t
        
    def snapshot_to_PeriodicNetwork(self):
        """Turn this into same format as PeriodicNetwork class """
        network = 0
        return network        

#start getting information from the network 
# extract total number of beads

mysystem = SystemBox() #import data/blocks from files
general_info = mysystem.open_files(snapshotname, d ,dr)

    


br_0 = dr.next_block()
NbeadsPart = br_0.meta.N
# map_RevEnds_to_star = imp.create_map(Natoms,NbeadsPart,br_0)

      
#set percolation as False by default - this will be changed to True if the network of irreversible bonds percolates
percolation = False       