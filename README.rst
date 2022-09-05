========
perconet
========

Overview
========

The **perconet** package provides tools to analyze the percolation properties of
periodic nets or networks. The terminology to describe such systems varies between
fields (see the sections *for mathematicians* and *for chemists* in the documentation).

Periodic networks arise often in simulations of chemical or physical systems with
*periodic boundary conditions*. Since periodic nets have no actual boundaries, the
question whether they percolate must be interpreted as *does the structure wrap around
the periodic boundary*?

**perconet** implements a loop finder algorithm that reports the number of *independent*
ways in which the structure wraps around the boundary.

Documentation
=============
Documentation is generated using Sphinx and hosted on `ReadTheDocs <https://perconet.readthedocs.io/>`_.

Release Status
==============
This is a development release that has been extensively tested in a few contexts but
we cannot guarantee it will work as you will expect. If you get unexpected results
and suspect you found a bug, please open an issue on github.


Credits and License
===================
Perconet was written by Chiara Raffaelli and Wouter G. Ellenbroek.
Issue reports and contributions are welcome through our `GitHub repository <https://github.com/wouterel/perconet>`_

We share **perconet** under the European Union Public License. See LICENSE file for details.