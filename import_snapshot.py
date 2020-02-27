#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 14:38:26 2019

@author: craffaelli
"""

from lammpstools import dump_reader
from lammpstools import block_data
import numpy as np
import pandas as pd

def snapshot_map(typeN, NatomsTot, f):
    """ Makes a map from the endbead id to the star id. """
    empty = False
    map_id_to_star = np.zeros(NatomsTot, dtype = int)
    while empty is False:
        l= f.readline()
        line = l.rstrip().split()
        if line:
            if line[2] == str(typeN):
                bead_id, star_id = int(line[0]), int(line[1])
                map_id_to_star[bead_id - 1] = star_id -2 # molecule_id correction  [map_id_to_star[bead_id - 1] = star_id -1] when bug is fixed
        else:
            empty = True
    return map_id_to_star





def snapshot_mapBonds(type_ix, f):
    empty = False
    crosslinks_list = []
    while empty is False:
        l= f.readline()
        line = l.rstrip().split()
        if line:
            if line[1] == str(type_ix):
                #print (line)
                id1, id2 = int(line[2]), int(line[3])
                crosslinks_list.append([id1-1, id2 -1])
        else:
            empty = True
    crosslinks_list = np.asarray(crosslinks_list)
    Nbonds = len(crosslinks_list)
    return crosslinks_list, Nbonds


#def snapshot_maptest(typeN, NatomsTot, f):
#    empty = False
#    map_id_to_star = np.zeros(NatomsTot, dtype = int)
#    while empty is False:
#        l= f.readline()
#        line = l.rstrip().split()
#        if line:
#            if line[2] == str(1):
#                print (line[1])
#        else:
#            empty = True
#    return map_id_to_star