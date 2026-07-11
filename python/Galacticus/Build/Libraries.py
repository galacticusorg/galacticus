"""Central registry of external-library knowledge for the build system.

This is the single place that records how Galacticus relates to external
libraries. Two build scripts consume it, each using a different half of the
model, and adding support for a new external library usually means touching
several of these tables together:

* ``useDependencies.py`` (dependency discovery) uses EXTERNAL_MODULES,
  MODULE_LIBRARIES, and INCLUDE_LIBRARIES to decide which ``use``/``include``
  statements imply a link-time library requirement (recorded in per-object
  ``.fl`` files).
* ``libraryDependencies.py`` (link-time flag generation) uses DEPENDENCIES,
  STATIC_LINK_DEPENDENCIES, and CONDITIONAL_FLAGS to expand the ``.fl``
  requirements into a correctly ordered ``-l...`` link line.

Historically these lived as separate hardcoded tables in the two scripts
with no cross-reference between them.
"""

__all__ = ['EXTERNAL_MODULES', 'MODULE_LIBRARIES', 'INCLUDE_LIBRARIES',
           'DEPENDENCIES', 'STATIC_LINK_DEPENDENCIES', 'CONDITIONAL_FLAGS']


# Fortran modules provided by the toolchain; the dependency scanner must not
# try to build them from Galacticus source.
EXTERNAL_MODULES = frozenset({
    'omp_lib', 'hdf5', 'h5tb', 'h5lt', 'h5global', 'h5fortran_types',
    'fox_common', 'fox_dom', 'fox_wxml', 'fox_utils', 'mpi', 'mpi_f08',
})

# Fortran modules whose `use` triggers a link-time library dependency.
MODULE_LIBRARIES = {
    'nearest_neighbors':  'ANN',
    'points_convex_hull': 'qhullcpp',
    'fftw3':              'fftw3',
    'fox_common':         'FoX_common',
    'fox_dom':            'FoX_dom',
    'fox_wxml':           'FoX_wxml',
    'fox_utils':          'FoX_utils',
    'hdf5':               'hdf5_fortran',
    'h5tb':               'hdf5_hl_fortran',
    'vectors':            'blas',
    'models_likelihoods': 'matheval',
    'input_parameters':   'matheval',
    'interface_gsl':      'gsl',
    'output_versioning':  'git2',
}

# Include files whose inclusion triggers a link-time library dependency.
INCLUDE_LIBRARIES = {
    'crypt': 'crypt',
}

# Library dependency graph: key depends on values. Used to transitively
# expand the set of libraries to link. (For static builds
# libraryDependencies.py adds an extra `hdf5 -> dl` edge to its own copy.)
DEPENDENCIES = {
    'hdf5_hl_fortran': ['hdf5_hl'],
    'hdf5_hl'        : ['hdf5'],
    'hdf5_fortran'   : ['hdf5'],
    'hdf5'           : ['z'],
    'gsl'            : ['gslcblas'],
    'FoX_dom'        : ['FoX_fsys', 'FoX_utils', 'FoX_sax'],
    'FoX_sax'        : ['FoX_common'],
    'FoX_utils'      : ['FoX_wxml'],
    'qhullcpp'       : ['qhull_r', 'stdc++'],
}

# Static-link ordering: key must appear before values on the link line.
STATIC_LINK_DEPENDENCIES = {
    'hdf5'           : ['z', 'dl'],
    'hdf5_hl'        : ['hdf5'],
    'hdf5_fortran'   : ['hdf5'],
    'hdf5_hl_fortran': ['hdf5_hl'],
    'gsl'            : ['gslcblas'],
    'FoX_dom'        : ['FoX_fsys', 'FoX_utils', 'FoX_sax', 'FoX_wxml'],
    'FoX_sax'        : ['FoX_common'],
    'FoX_wxml'       : ['FoX_utils'],
    'FoX_common'     : ['FoX_fsys'],
    'qhullcpp'       : ['stdc++'],
}

# Conditionally compiled libraries: linked only when the corresponding
# availability flag (from the Makefile_Config_* probes) is present in the
# compile flags.
CONDITIONAL_FLAGS = {
    'fftw3'   : '-DFFTW3AVAIL',
    'ANN'     : '-DANNAVAIL',
    'qhullcpp': '-DQHULLAVAIL',
    'qhull_r' : '-DQHULLAVAIL',
    'matheval': '-DMATHEVALAVAIL',
    'git2'    : '-DGIT2AVAIL',
}
