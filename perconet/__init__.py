# This file is part of the perconet package
# (c) 2022 Eindhoven University of Technology
# Released under EUPL v1.2
# See LICENSE file for details
# Contributors:
# * Chiara Raffaelli
# * Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet
"""
Perconet is a package for analyzing the percolation properties of periodic
nets.

Simple examples of periodic nets are a cubic lattice in three dimensions with
periodic boundary conditions, or a honeycomb lattice in two dimensions with
periodic boundary conditions.

In physical systems, the nodes (points) in the network will typically have
a position associated with them, and the bonds (edges) can be defined to
connect two nodes. The effect of periodic boundaries is captured by associating
a boundary-crossing vector with each edge.

Mathematically speaking, any graph that is embedded in a 3-torus can be
analyzed with this package. It is not necessary to a associate positions with
the nodes, as long as the boundary-crossing properties of the edges are
prescribed correctly.
"""

from perconet.periodicnetwork import PeriodicNetwork
from perconet.loopfinder import LoopFinder

__all__ = ["PeriodicNetwork", "LoopFinder"]

__version__ = "0.2.2"
