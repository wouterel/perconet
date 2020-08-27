#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 09:58:09 2019

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




def import_rev (bb, maxBonds, id0, Nbonds):
    reversible_bonds = []
    for i in range (maxBonds):
        ix0, bond = int(bb.data_by_name(id0)[i] -1), int(bb.data_by_name(Nbonds)[i])
        if bond == 1:
            reversible_bonds.append ([ix0, -1])
        elif bond != 0:
            ('simulation problem! more than one bead is reversibly bound to this bead! index = ', ix0, ', n bonds = ', bond)
            
    reversible_bonds = np.asarray(reversible_bonds)
    return reversible_bonds, len(reversible_bonds)


#2.map the bond id to its coordinates
    
def map_coords(ixMax, N, b):
    """ Makes a map from the bead id to its center coordinates. """
    bead_coords = np.zeros((ixMax ,3), dtype = float)
    for index in range(N):
        x = b.x[index]
        bead_id = int (b.ids[index])

        bead_coords[bead_id - 1, :] = x[:]
    return bead_coords
