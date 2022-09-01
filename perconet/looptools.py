# This file is part of the perconet package
# (c) 2022 Eindhoven University of Technology
# Released under EUPL v1.2
# See LICENSE file for details
# Contributors:
# * Chiara Raffaelli
# * Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import numpy as np


def swaprows(a: np.ndarray, i1: int, i2: int):
    """
    Swap (in-place) rows i1 and i2 in twodimensional :obj:`ndarray` a.

    Will throw an error if called on an array with fewer than 2 dimensions.
    """
    a[(i1, i2), :] = a[(i2, i1), :]


def multiply_row(mat, shadow=None, row=0, factor=1):
    """
    Mutiply row of matrix mat by factor. Note that the only nontrivial
    multiplication that is allowed is by a factor -1 when reducing a lattice basis.

    If set, apply the same transformation to the shadow matrix.
    This is useful to keep track of the transformation matrix representing
    the elimination result.
    """
    mat[row, :] = mat[row, :] * factor
    if isinstance(shadow, np.ndarray):
        shadow[row, :] = shadow[row, :] * factor


def reduce_row(mat, shadow=None, pivot_col=0, source_row=0, target_row=0):
    """
    Add integer multiple of source row of mat to target row of mat such that
    the pivot element mat[target_row, pivot_col] becomes as close to 0
    as possbible (rounds such that it ends up nonnegative).

    If set, apply the same transformation to the shadow matrix.
    This is useful to keep track of the transformation matrix representing
    the elimination result.
    """
    if target_row == source_row:
        return
    q = mat[target_row, pivot_col] // mat[source_row, pivot_col]  # integer division
    if q != 0:
        # update matrix a and keep track of the corresponding transformation matrix
        mat[target_row, :] -= q * mat[source_row, :]
        if isinstance(shadow, np.ndarray):
            shadow[target_row, :] -= q * shadow[source_row, :]


def hermite_normal_form(input: np.ndarray):
    """
    Construct the Hermite normal form of the input matrix using only
    * row swapping;
    * addition of integer multiples of other rows to rows;
    * multiplication of entire rows by -1;

    Returns:
        Tuple[:obj:`np.ndarray`, :obj:`np.ndarray`, int]:
            Tuple (r, u, rank) with r representing the Hermite normal form of input,
            u representing the orthogonal transformation matrix such that
            u*input=r and rank the number of nonzero rows of r which equals
            the rank of input.
    """
    a = input.copy()
    (m, n) = a.shape  # number of rows, columns

    # initialize the matrix that will contain the transformation representing the sweep
    transform = np.eye(m, dtype=int)
    col = 0  # pivot column
    row = 0  # current working row (rows with lower number are no longer involved)
    rank = m  # default: will be changed if zero-rows are detected later

    # if the matrix is full rank the loop will end when row reaches m
    # else the loop will end when col reaches n
    while col < n and row < m:
        # find on which rows (starting from "row") column "col" has nonzero elements
        nonzero_rows = np.nonzero(a[row:, col])[0]
        if len(nonzero_rows) == 0:
            # no nonzero elements here. stay on same row but move 1 column to the right
            col += 1
            if col == n:
                # reached the last column while not yet on last row. set rank variable
                rank = row
                break
            continue
        # find which row has (in this column) the smallest absolute-value entry
        smallest_abs = np.min(np.abs(a[row+nonzero_rows, col]))
        source_row = np.where(np.abs(a[row:, col]) == smallest_abs)[0][0]+row
        to_next_col = True  # default value. will be set to false later if nonzero elements remain
        for target_row in range(row, m):
            if source_row == target_row:
                continue
            reduce_row(a, shadow=transform, pivot_col=col, source_row=source_row,
                       target_row=target_row)
            if a[target_row, col] != 0:
                # at least one of the rows still has a nonzero element
                # this happens when not all elements in this column are multiples of the smallest
                # will need to do this column again with whichever row is now the smallest.
                to_next_col = False
        if to_next_col:
            # all relevant elements in this column are 0 now.
            # swap source_row to the top of active block and increase row variable so it will
            # no longer be touched
            if source_row != row:
                swaprows(a, source_row, row)
                swaprows(transform, source_row, row)
            # pivot element is now in "row"
            # HNF rule: pivot element must be positive
            if a[row, col] < 0:
                multiply_row(a, shadow=transform, row=row, factor=-1)
            # Now finalize the HNF by using the just-finished row to make the elements
            # above the pivot element nonnegative and smaller than the pivot
            for target_row in range(row):
                reduce_row(a, shadow=transform, pivot_col=col, source_row=row,
                           target_row=target_row)
            row += 1
            # It's tempting to also update col here
            # but then we'd have to duplicate the rank-checking code
    return (a, transform, rank)
