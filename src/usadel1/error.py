"""
:Author:       Pauli Virtanen <pauli@ltl.tkk.fi>
:Organization: Low Temperature Laboratory, Helsinki University of Technology
:Date:         2005-2006

Errors, exceptions and warnings in usadel1.
"""

import sys

__docformat__ = "restructuredtext en"
__all__ = ['NoConvergenceException',
           'CoreError',
           'CoreWarning',
           'clear_warnings',
           'get_warnings',
           'set_warnings',
           'suppress_warnings']

import warnings

class CoreError(RuntimeError):
    """An error in usadel1 Fortran core."""
    pass

class CoreWarning(Warning):
    pass

class NoConvergenceWarning(CoreWarning):
    pass

class NoConvergenceException(Exception):
    """A numerical algorithm failed to converge.

    :Ivariables:
      - `result`: a partial result or other information on the calculation
    """
    def __init__(self, msg=None, result=None):
        Exception.__init__(self, msg)
        self.result = result



def clear_warnings():
    """Clear all warnings. (DEPRECATED)"""
    warnings.warn("clear_warnings is deprecated", DeprecationWarning)

def get_warnings():
    """Get current warnings. (DEPRECATED)"""
    warnings.warn("get_warnings is deprecated", DeprecationWarning)
    return None

def set_warnings(exception, nwarnings=1):
    """Set a warning. (DEPRECATED)"""
    warnings.warn("get_warnings is deprecated", DeprecationWarning)
    warnings.warn(str(exception), CoreWarning)

def suppress_warnings(supressed):
    """Suppress/enable automatic printing of warning messages. (DEPRECATED)"""
    warnings.warn("get_warnings is deprecated", DeprecationWarning)
    if suppressed:
        warnings.simplefilter("ignore", CoreWarning)
    else:
        warnings.simplefilter("default", CoreWarning)
