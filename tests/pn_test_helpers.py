
import perconet as pn
import numpy as np


def initialize_test(bondlist):
    assert len(bondlist) > 0
    for bond in bondlist:
        assert len(bond) == 5
    bonds = np.asarray(bondlist).astype(int)
    nodelist_bonds = bonds[:, 0:2]
    n_nodes = nodelist_bonds.max()+1
    assert n_nodes < 100000
    assert n_nodes > 0
    nodenum, coordination = np.unique(nodelist_bonds, return_counts=True)
    max_coord = coordination[np.argmax(coordination)]
    assert max_coord > 0
    network = pn.PeriodicNetwork(n_nodes, max_degree=max_coord, verbose=False)
    for bond in bondlist:
        assert network.add_edge(bond[0], bond[1], bond[2:])
    return network
