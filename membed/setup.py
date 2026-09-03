# Cython build helper for the GloVe co-occurrence counter in glove_build/.
# NOT the package build config: the distributable package is built by the
# setup.py at the repository root. Run this from glove_build/ (requires Cython)
# to rebuild the `cooccur` executable from compute_cooccur.pyx.
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("compute_cooccur.pyx")
)
