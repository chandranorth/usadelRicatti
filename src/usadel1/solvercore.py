"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

Solver core interface.

This module adds error handling to the routines imported from the
Fortran module.

Note that not all logic is contained in the Fortran module `_solvercore`,
but an important part is in `solver`.
"""

import _solvercore as _sc
from error import *
import warnings

__docformat__ = "restructuredtext en"
__all__ = ['set_delta', 'set_kinetic', 'set_params',
           'set_solvers', 'sp_initialize', 'kin_initialize',
           'sp_solve', 'kin_solve']

## Decorate _sc:

def decorate_sc_func(func, docstring=""):
    def _wrapper(*args, **kwargs):
        ret = func(*args, **kwargs)
        errmsg, errflag, warnmsg,warnflag,nwarnings = _sc.spectral2_get_error()
        if errflag != 0:
            raise CoreError(errmsg.strip())
        if warnflag != 0:
            if nwarnings <= 1:
                warnings.warn(warnmsg.strip(), CoreWarning)
            else:
                warnings.warn(warnmsg.strip()
                              + " (+ %d other warnings)" % nwarnings,
                              CoreWarning)
        return ret
    basedoc = "    " + func.__doc__.replace("\n", "\n    ")
    if docstring:
        _wrapper.__doc__ = docstring + "\n\n::\n\n" + basedoc
    else:
        _wrapper.__doc__ = "::\n\n" + basedoc
    return _wrapper

set_delta = decorate_sc_func(_sc.spectral2_set_delta, 
                             """Set order parameter""")
set_kinetic = decorate_sc_func(_sc.spectral2_set_kinetic,
                               """Set kinetic coefficients""")

set_params = decorate_sc_func(_sc.spectral2_set_params,
                              """Set geometry""")
set_solvers = decorate_sc_func(_sc.spectral2_set_solvers,
                               """Set solver parameters""")

sp_initialize = decorate_sc_func(_sc.spectral2_sp_initialize,
                                 """Initialize spectral solver""")
kin_initialize=decorate_sc_func(_sc.spectral2_kin_initialize,
                                  """Initialize kinetic solver""")

sp_solve = decorate_sc_func(_sc.spectral2_sp_solve,
                            """Solve spectral equations""")
kin_solve = decorate_sc_func(_sc.spectral2_kin_solve,
                             """Solve kinetic equations""")

