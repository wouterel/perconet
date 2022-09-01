When and why to use
===================

`Insert graphics and main motivation` bla

.. _Loop independence:

Loop independence
-----------------

If the loop finder identifies a loop that goes around both the :math:`+x` and :math:`+y`
boundaries :math:`\left[\vec{b}_1=(1,1,0)\right]`, and another loop that only goes around the :math:`+x`
boundary :math:`\left[\vec{b}_2=(1,0,0)\right]`, we can construct a loop with :math:`\vec{b}=(0,1,0)` by first
going around the first loop and then going around the second loop in reverse: :math:`\vec{b}=\vec{b}_1-\vec{b}_2`.
Generalizing, any linear combination of loops with integer coefficients is also a loop. Thus it makes
sense to reduce the list of loops to a list of `independent` loops by constructing a basis
of independent loops. Because the basis is to be used only with integer coefficients (one cannot
go around a loop half a time), it is a lattice basis and the space of allowed loops is a lattice.
Writing the list of loops as a matrix (each row representing a loop), the reduction is like gaussian
elimination, but with the constraint that only integer multiples of loops can be added to other
loops and multiplying a row by a constant is not allowed (except for -1 which is just reversing the
direction of a loop).

A way of reducing the list of loops to a list of independent loops that gives a unique result, so
one can compare different loop structures, is to cast it into Hermite normal form. See
`Wikipedia <https://en.wikipedia.org/wiki/Hermite_normal_form>`_ or your favorite linear algebra text
for details. This is the form :py:meth:`perconet.LoopFinder.get_independent_loops` returns. Note that
the exact definition of Hermite normal form varies slightly between authors.