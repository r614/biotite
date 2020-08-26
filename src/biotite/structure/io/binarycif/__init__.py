# This source code is part of the Biotite package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

"""
This subpackage provides support for the the modern BinaryCIF file
format. The :class:`BinaryCIFFile` class provides dictionary-like access to
every field in BinaryCIF files.
Additional utility functions allow conversion of these dictionaries to
:class:`AtomArray` and :class:`AtomArrayStack` objects and vice versa.
"""

__name__ = "biotite.structure.io.binarycif"
__author__ = "Roshan Pawar"

from .file import *
