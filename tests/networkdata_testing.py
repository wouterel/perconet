# This file is part of the perconet package
# (c) 2022 Eindhoven University of Technology
# Released under EUPL v1.2
# See LICENSE file for details
# Contributors:
# * Chiara Raffaelli
# * Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import numpy as np

edgelist_noloops = np.array([
        [0, 23,  0, -1,  0],
        [2, 56,  0, -1,  0],
        [5, 25,  1,  0,  0],
        [12, 26,  0,  0, -1],
        [15, 55,  0,  1,  0],
        [15, 75,  0,  1,  0],
        [21, 64,  0,  0, -1],
        [26, 42,  1,  0,  0],
        [28,  2,  0,  1,  0],
        [33, 58, -1,  0,  1],
        [38, 51, -1,  0,  0],
        [39, 19,  0,  0,  1],
        [45, 22,  0,  0, -1],
        [48, 65,  0,  0,  1],
        [50, 65,  0,  0,  1]])
testcase_noloops = (edgelist_noloops, [[]])

edgelist_example1 = np.array([
        [1, 8,  0, -1,  0],
        [2, 8,  0, -1,  0],
        [3, 4,  1,  0,  0],
        [5, 2,  0,  0, -1],
        [6, 2,  0,  0,  -1],
        [8, 6,  0,  0,  1],
        [8, 9,  0,  0, 1],
        [8, 10,  0,  0,  1],
        [9,  2,  0,  0,  -1]])


# specific cases (see powerpoint)
edgelist_A = np.array([
        [1, 2,  -1, 0,  0],
        [2, 3,  -1, 0,  0],
        [3,  1,  -1,  0,  0]])
solution_A = [[-3, 0, 0]]
testcase_A = (edgelist_A, solution_A)

edgelist_B = np.array([
        [1, 2,  1, 1,  0],
        [1, 4,  1, 0,  0],
        [1, 7,  1,  0,  0],
        [2, 3,  1,  1, 0],
        [2, 3,  -1,  -1,  0],
        [3, 4,  0,  1,  0],
        [4, 1,  0,  -1, 0],
        [4, 5,  0,  -1, -1],
        [5,  3,  0,  1,  0],
        [5, 6, 1,  1,  1],
        [6, 7, 0,  -1,  -1],
        [7, 2,  1,  0,  0],
        [7, 8,  0,  0, -1],
        [8, 5,  0,  -1,  0],
        [8, 9,  0,  0,  1]])
solution_B = [[2, 2, 0], [1, 3, 0], [0, 1, -1]]
testcase_B = (edgelist_B, solution_B)

edgelist_C = np.array([
        [1, 2,  1, 1,  0],
        [1, 4,  1, 0,  0],
        [1, 7,  1,  0,  0],
        [2, 3,  1,  1, 0],
        [2, 3,  -1,  -1,  0],
        [3, 4,  0,  1,  0],
        [4, 1,  0,  -1, 0],
        [4, 5,  0,  -1, -1],
        [5,  3,  0,  1,  0],
        [5, 6, 1,  1,  1],
        [6, 7, 0,  -1,  -1],
        [7, 2,  1,  0,  0],
        [7, 8,  0,  0, -1],
        [8, 5,  0,  -1,  0],
        [8, 9,  0,  0,  1],
        [1, 2,  0, 0,  0],
        [3, 4,  0, 0,  0],
        [1, 9,  0, 0,  0],
        [2, 4,  0, 0,  0],
        [1, 6,  0, 0,  0],
        [6, 3,  0, 0,  0],
        [1, 3,  0, 0,  0],
        [1, 3,  0, 0,  0]])
solution_C = [[2, 2, 0], [1, 3, 0], [0, 1, -1]]
testcase_C = (edgelist_C, solution_C)

edgelist_D = np.array([
        [1, 2,  1, 0, 0],
        [2, 3,  0, 0, 0],
        [3, 2,  0, 0, 1],
        [2, 1,  0, 1, 0]])
solution_D = [[1, 1, 1], [1, 1, 0]]
testcase_D = (edgelist_D, solution_D)


# keep next line for backwards compatibility
edgelist = edgelist_C
