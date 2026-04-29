# Meta-property type table shared by the components-build pipeline.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Classes/MetaProperties.pm — for
# now we only port the data table (`meta_property_types`).  The two
# `Class_Meta_Property_Get` / `Class_Meta_Property_Set` hooks register on
# the `classIteratedFunctions` phase, which is a sub-iteration owned by
# Classes.pm; they will land alongside that port.

# The order of entries is preserved verbatim from the Perl module so that
# generated identifier suffixes (`Float`, `LongInteger`, `Integer` × ranks
# 0/1) keep matching what the rest of the pipeline expects.
meta_property_types = [
    {'label': 'float',       'intrinsic': 'double precision',                          'rank': 0},
    {'label': 'float',       'intrinsic': 'double precision',                          'rank': 1},
    {'label': 'longInteger', 'intrinsic': 'integer',          'type': 'kind_int8',     'rank': 0},
    {'label': 'longInteger', 'intrinsic': 'integer',          'type': 'kind_int8',     'rank': 1},
    {'label': 'integer',     'intrinsic': 'integer',                                   'rank': 0},
    {'label': 'integer',     'intrinsic': 'integer',                                   'rank': 1},
]
