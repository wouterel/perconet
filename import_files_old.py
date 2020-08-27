#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 13:38:41 2018

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

# functions:
# synchronise dump files

#def synchro(d, dl):
#    bl = dl.next_block()
#    b = d.next_block()
#    #print ()
#    if b is None or bl is None:
#        return None, None
#    else:
#        while d and dl:
#            if b.meta.t == bl.meta.t:
#                print('synched! t= ', b.meta.t)
#                return b,bl
#            elif b.meta.t < bl.meta.t:
#                print('bond dump ahead, t=', b.meta.t, 'tb=', bl.meta.t)
#                b = d.next_block()
#            elif b.meta.t > bl.meta.t:
#                print('dump ahead, t=', b.meta.t, 'tb=', bl.meta.t)
#                bl = dl.next_block()
                
def synchro(d, dl, db):
    bl = dl.next_block()
    b = d.next_block()
    bb = db.next_block()
    #print ()
    if b is None or bl is None or bb is None:
        return None, None, None
    else:
        while d and dl and bb:
            if b.meta.t == bl.meta.t:
                if bb.meta.t == b.meta.t:
                    print('synched! t= ', b.meta.t)
                    return b,bl,bb
                elif b.meta.t < bb.meta.t:
                    print('bond dump ahead, t=', b.meta.t, 'tb=', bb.meta.t)
                    b = d.next_block()
                elif b.meta.t > bb.meta.t:
                    print('dump ahead, t=', b.meta.t, 'tb=', bb.meta.t)
                    bb = db.next_block()
            elif b.meta.t == bb.meta.t:
                if b.meta.t < bl.meta.t:
                    print('bond dump ahead, t=', b.meta.t, 'tb=', bl.meta.t)
                    b = d.next_block()
                elif b.meta.t > bl.meta.t:
                    print('dump ahead, t=', b.meta.t, 'tb=', bl.meta.t)
                    bl = dl.next_block()
            elif bl.meta.t == bb.meta.t:
                if b.meta.t < bl.meta.t:
                    print('bond dump ahead, t=', b.meta.t, 'tb=', bl.meta.t)
                    b = d.next_block()
                elif b.meta.t > bl.meta.t:
                    print('dump ahead, t=', b.meta.t, 'tb=', bl.meta.t)
                    bl = dl.next_block()
            elif bl.meta.t != bb.meta.t != b.meta.t:
                print ('bug! cannot synch files')
                return None, None, None
        
                
# skip lines
def skip_lines(n,f):
    for i in range(n):
        next(f)
        
def is_title(s,keyword):
    """ function to check if a line
         starts with some character.
    """
    # return true if a line starts with #
    return s.startswith(keyword)

# find keyword and skip one line
def find_keyword(keyword, f):
    """ Finds keyword and skip one line """
    check = False
    while check is False:
        s = f.readline().rstrip()
        check = is_title(s, keyword)
        #print(s)
    else:
        next(f) 
               

       
def meta(length, typ, f):
    l= f.readline()
    line = l.rstrip().split()
    if length == 1:
        val = typ(line[0])
    else:
        val = []
        for i in range(length):
            val.append(typ(line[i]))
    return val

# 2. create map to map the bond id to its star index
def create_map(NatomsTot, number_of_beads, b):
    """ Makes a map from bond id to star index. """
    map_id_to_star = np.zeros(NatomsTot, dtype = int)

    for index in range(number_of_beads):
        bead_id = int(b.ids[index])
        star_id = int(b.mol[index]) - 2  # molecule_id correction [star_id = int(b.mol[index]) -1 ]

        map_id_to_star[bead_id - 1] = star_id
        
    return map_id_to_star

#2.map the bond id to its center coordinates
    
def map_coords(Nstars, b):
    """ Makes a map from the start id to its center coordinates. """
    star_coords = np.zeros((Nstars ,3), dtype = float)
    for index in range(Nstars):
        
        x = b.x[index]
        star_id = int (b.mol[index]) -1 # molecule_id correction [star_id = int (b.mol[index]) ] 

        star_coords[star_id - 1, :] = x[:]
    return star_coords
    
# 3. make N_stars*N_stars array of zeros; add +1 everytime a star polymer is bonded to another
    
# function to store bounary crossing
    
def crossings(id1,id2, coord, boxsize):
    crossing = np.zeros(3, dtype =int)
    #id1,id2 = id1 -1, id2 -1
    #printt = False
    for i in range (3):
        dist = coord[id2, i] - coord[id1, i]
        if dist >= boxsize[i]/ float(2):
            crossing[i] += 1
            #printt =True
        elif dist <= - boxsize[i]/float(2):
            crossing[i] -= 1
            #printt = True
    
    #if id1 == 14 and id2 ==136:
        #print("id14 and id136: ", crossing)
    return crossing

def block_to_array(N_bonds,bl, id1,id2):
    """ Converts data about bonds in a block to an array 
        containing bond information (so that it is homogeneous
        with the snapshot import)."""
    crosslinks_list = []
    for bond in range(N_bonds):
        # for each bond, gather the ids of the two bonded beads
        first_bead, second_bead = int(bl.data_by_name(id1)[bond] -1), int(bl.data_by_name(id2)[bond] -1)
        crosslinks_list.append([first_bead, second_bead])
    crosslinks_list = np.asarray(crosslinks_list)
        
            
def store_crosslinks1(N_stars, N_bonds, map_id_to_star,crosslinks_list, coord, boxsize):
    """ Makes an N_stars*N_stars array with {c(i,j)}, where c(i,j) counts
        how many times molecule i is bonded to molecule j. """
       
    crosslinks = np.zeros((N_stars,N_stars), dtype = int)
    
    crosslink_boundaries = np.zeros((N_stars,N_stars, 3), dtype = int)
    
    for bond in range(N_bonds):
        # for each bond, gather the ids of the two bonded beads
        first_bead, second_bead = crosslinks_list [bond, :]
        #print (first_bead, second_bead)
        bonded_star1, bonded_star2 = map_id_to_star[first_bead], map_id_to_star[second_bead]
        bondedID = [bonded_star1, bonded_star2]
        bondedID.sort() #always smaller ID first
        #print (bondedID)
        #crosslinks[bondedID[0]-1, bondedID[1]-1] += 1
        crosslinks[bondedID[0], bondedID[1]] = 1 # DOUBLE BONDS ARE NOT REGISTERED ANYMORE
        crosslink_boundaries[bondedID[0], bondedID[1] , :] += crossings(bondedID[0], bondedID[1], coord,boxsize)[:]
    crosslink_boundaries = np.asarray(crosslink_boundaries)
    return crosslinks , crosslink_boundaries

def store_crosslinks(N_stars, N_bonds, map_id_to_star, bl, id1,id2, coord, boxsize):
    """ Makes an N_stars*N_stars array with {c(i,j)}, where c(i,j) counts
        how many times molecule i is bonded to molecule j. """
    crosslinks = np.zeros((N_stars,N_stars), dtype = int)
    
    crosslink_boundaries = np.zeros((N_stars,N_stars, 3), dtype = int)
    
    for bond in range(N_bonds):
        # for each bond, gather the ids of the two bonded beads
        first_bead, second_bead = int( bl.data_by_name(id1)[bond]), int( bl.data_by_name(id2)[bond])
        bonded_star1, bonded_star2 = map_id_to_star[first_bead -1], map_id_to_star[second_bead -1]
        bondedID = [bonded_star1, bonded_star2]
        bondedID.sort() #always smaller ID first
        #crosslinks[bondedID[0]-1, bondedID[1]-1] += 1
        crosslinks[bondedID[0]-1, bondedID[1]-1] = 1 # DOUBLE BONDS ARE NOT REGISTERED ANYMORE
        crosslink_boundaries[bondedID[0]-1, bondedID[1]-1 , :] += crossings(bondedID[0], bondedID[1], coord,boxsize)[:]
    crosslink_boundaries = np.asarray(crosslink_boundaries)
    return crosslinks , crosslink_boundaries
        
        # check if the bonds are further away than half the box
    