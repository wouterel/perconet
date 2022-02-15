Information for...
==================

Those confused about terminology
--------------------------------
Most of the jargon used comes from the mathematics of graphs, with nodes (points) connected
by edges (lines). In the context of gelation, the nodes will represent molecules or
colloids, and the edges will represent chemical or physical bonds. In many cases, edges
or bonds may also be called *links*.

We use the mathematical term *graph* to denote any graph, and the term *periodic nets*
or *periodic networks* to denote graphs that are embedded in a periodic box.

Chemists
--------
The most obvious use case in chemistry for **perconet** is detecting *gelation*.
Models for a gelation process with periodic boundary conditions will typically lead
to data that specifies positions for the building blocks (monomers) and a list
of bonds that have been generated during the gelation process. There will typically
be one large molecule and many small ones, and **perconet** will determine
for you whether that large molecule connects to itself around the periodic
boundary, signalling the presence of an infinite molecule, the *gel*.

While primarily written for periodic systems, it is also possible to use **perconet**
for percolation analysis of systems with simple boundaries. To this end, denote
a certain subset of the nodes to be one *side* of the system, and another subset
to be the other side, and then ask **perconet** if the two sides are connected.
With this approach, even the output of an experimental image analysis 
process could be used as input.


Mathematicians
--------------
The three-dimensional periodic boxes that inspired this package are a 
topological space known as a 3-torus. The use of the package is, however, not
limited to three dimensions and can be used to analyze graphs embedded
in any cartesian power of the circle :math:`\mathbb{T}^d=S^d`.

A loop in such an embedded graph is characterized by an element of the fundamental group
of the `d`-torus, which is :math:`\mathbb{Z}^d`. The element specifies, for each
periodic boundary, how often the loop in question goes around that boundary.
 
Not all elements of the fundamental group of the `d`-torus are necessarily represented
in every graph: Perhaps it only wraps around one of the boundaries,
or some boundary can only be looped around an even number of times.
The subgroup of :math:`\mathbb{Z}^d` that is actually realized by the periodic
net is a lattice, for which the method
:py:meth:`perconet.LoopFinder.get_independent_loops` provides a basis.
One can also say these are the generators of the subgroup. The basis
is provided in a row echelon form which is not unique: To determine
whether two graphs have the exact same percolation structure one should
compare the Hermite Normal Forms of the matrices representing the bases.
This is a feature we may add to **perconet** in the near future. For some
applications it may be desirable to have a near-orthogonal basis, in which
case improving it via the LLL-algorithm may prove useful.


Physicists and mechanical engineers
-----------------------------------

bla