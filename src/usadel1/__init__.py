"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

A solver for Keldysh-Usadel 1D circuit equations.
 
To find out how this library can be used, you should peek at

- `solver.Geometry`: How to specify a geometry
- `currents`: Simple interface to calculating currents
- `selfconsistentiteration`: A self-consistent iteration
- `solver`: Low-level solver interface

"""

# Modules to import
from solver import *
from selfconsistentiteration import *
from nonlinearsolver import *
from currents import *
from util import *
from error import *
from version import __version__

__docformat__ = "restructuredtext en"
