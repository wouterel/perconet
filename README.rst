========
perconet
========

.. image:: https://shields.io/pypi/v/perconet
    :target: https://pypi.org/project/perconet
    :alt: PyPI Version


.. image:: https://readthedocs.org/projects/perconet/badge/?version=stable
    :target: https://perconet.readthedocs.io/en/stable/?badge=stable
    :alt: Documentation Status


.. image:: https://shields.io/pypi/l/perconet
    :target: https://github.com/wouterel/perconet/blob/develop/LICENSE
    :alt: EUPL-1.2

.. image:: https://github.com/wouterel/perconet/actions/workflows/test.yaml/badge.svg?branch=develop
   :target: https://github.com/wouterel/perconet/actions/workflows/test.yaml
   :alt: Test

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

Installation
============
**perconet** can be installed in any Python environment (version 3.8 or higher) using ``pip install perconet``.
See `https://pypi.org/project/perconet/ <https://pypi.org/project/perconet/>`_ for package details.

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
