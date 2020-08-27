#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 09:46:16 2019

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


def check_distance(id1,id2,coords, boxsize, cutoffSq):
    x1, x2 = coords[id1,:],coords[id2,:] 
    dxSq = 0
    for i in range(3):
        dxi = abs (x1[i] - x2[i])
        if dxi >= boxsize[i]/float(2):
            dxi -= boxsize[i]
        dxSq += dxi ** 2
    #print (dxSq)
    if dxSq <= cutoffSq:
        #print (dxSq , cutoffSq)
        return True
    else:
        return False
    

def make_revbondlist (revList, boundN, revEnd_coordinates, boxsize, cutoffSq):


    for ix1 in range (boundN):
        if revList[ix1, 1] == -1:
            #print (revList[ix1,:])
            for ix2 in [x for x in range(boundN) if x != ix1]: 
                if revList[ix2,1] == -1:
                    id1, id2 = revList[ix1, 0], revList[ix2, 0]
                    is_bound = check_distance(id1,id2, revEnd_coordinates, boxsize, cutoffSq) 
                    if is_bound:
                       revList[ix1, 1] = id2
                       revList[ix2, 1] = id1
                       break
                
    return revList
                   
            
        
                    
                       
                      
                
    