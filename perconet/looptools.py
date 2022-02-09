import numpy as np


def swaprows(a: np.ndarray, i1: int, i2: int):
    """
    Swap (in-place) rows i1 and i2 in twodimensional :obj:`ndarray` a.

    Will throw an error if called on an array with fewer than 2 dimensions.
    """
    a[(i1, i2), :] = a[(i2, i1), :]


def integer_gaussian_elimination(input: np.ndarray):
    """
    Construct a row echelon form of the input matrix using only row swapping and
    addition of integer multiples of other rows to rows.

    Returns:
        Tuple[:obj:`np.ndarray`, :obj:`np.ndarray`, int]:
            Tuple (r, u, rank) with r representing the row echelon form of input,
            u representing the orthogonal transformation matrix and rank the
            number of nonzero rows of r.
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
            # we will be adding multiples of row source_row to all other rows
            # to make their elements in the current column as close to 0 as we can
            if target_row == source_row:
                continue
            q = a[target_row, col] // a[source_row, col]  # integer division
            if q != 0:
                # update matrix a and keep track of the corresponding transformation matrix
                a[target_row, :] -= q * a[source_row, :]
                transform[target_row, :] -= q * transform[source_row, :]
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
            row += 1
            # It's tempting to also update col here
            # but then we'd have to duplicate the rank-checking code
    return (a, transform, rank)
