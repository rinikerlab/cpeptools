"""
cPepTools
Various tools to aid the simulational study of cyclic (mainly) peptidic species.[D[D[
"""

# Make Python 2 and 3 imports work the same
# Safe to remove with Python 3-only code
from __future__ import absolute_import

# Add imports here
from .cpeptools import *
# from .confgen.confgen import *
from .confgen.utils import bound_matrix_from_ellipse, get_amide_pairwise_coulomb_interaction
from .mol_ops import *
# from .confgen.postprocessing import *
# from .molgen.peptide import *
# from .molgen.fragment import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
