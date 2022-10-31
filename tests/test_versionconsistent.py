# This file is part of the perconet package
# (c) 2022 Eindhoven University of Technology
# Released under EUPL v1.2
# See LICENSE file for details
# Contributors:
# * Chiara Raffaelli
# * Wouter G. Ellenbroek
# Contributing: https://github.com/wouterel/perconet

import perconet as pn
import toml
import os


def test_version():
    testdir = os.path.dirname(os.path.abspath(__file__))
    tomlfile = os.path.join(testdir, "..", "pyproject.toml")
    release = toml.load(tomlfile)['tool']['poetry']['version']
    assert release == pn.__version__


if __name__ == "__main__":
    test_version()
