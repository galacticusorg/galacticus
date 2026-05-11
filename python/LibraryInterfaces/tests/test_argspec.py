"""Tests for `LibraryInterfaces.ArgSpec`.

Tiny dataclass — but it's the IR shared across every emitter and
pipeline stage in the cross-language build, so the field defaults
and the `from_raw` adapter need to behave exactly as documented.
"""

from LibraryInterfaces.ArgSpec import ArgSpec


# ---------------------------------------------------------------------------
# Field defaults
# ---------------------------------------------------------------------------

def test_minimal_construction_sets_documented_defaults():
    """An ArgSpec built with only `name` carries the documented default for
    every other field."""
    a = ArgSpec(name='x')
    assert a.name              == 'x'
    assert a.intrinsic         == ''
    assert a.type_spec         == ''
    assert a.attributes        == []
    assert a.is_optional       is False
    assert a.is_function_class is False
    assert a.ctype             == ''
    assert a.ctype_pointer     is False
    assert a.fort_type         == ''
    assert a.fort_is_present   is True
    assert a.fort_attributes   == []
    assert a.fort_pass_as      == ''
    assert a.py_is_present     is True
    assert a.galacticus_is_present is True
    assert a.pass_by           == ''


def test_attributes_default_is_independent_per_instance():
    """The `attributes` list must NOT be shared between instances (a
    classic mutable-default-argument bug — `field(default_factory=list)`
    prevents it)."""
    a = ArgSpec(name='a')
    b = ArgSpec(name='b')
    a.attributes.append('intent(in)')
    assert b.attributes == []


def test_fort_modules_default_is_independent_per_instance():
    a = ArgSpec(name='a')
    b = ArgSpec(name='b')
    a.fort_modules['mymod'] = {'sym': 1}
    assert b.fort_modules == {}


# ---------------------------------------------------------------------------
# from_raw — adapter from raw Fortran-declaration dict
# ---------------------------------------------------------------------------

def test_from_raw_minimal_dict():
    """A dict with only `name` produces an ArgSpec with empty strings
    everywhere else."""
    a = ArgSpec.from_raw({'name': 'x'})
    assert a.name       == 'x'
    assert a.intrinsic  == ''
    assert a.type_spec  == ''
    assert a.attributes == []


def test_from_raw_populated_dict():
    a = ArgSpec.from_raw({
        'name':       'mass',
        'intrinsic':  'double precision',
        'type':       '',
        'attributes': ['intent(in)', 'optional'],
    })
    assert a.name       == 'mass'
    assert a.intrinsic  == 'double precision'
    assert a.attributes == ['intent(in)', 'optional']


def test_from_raw_None_intrinsic_becomes_empty_string():
    """Defensive: if `intrinsic` is missing or None, use `''` so downstream
    string comparisons work without isinstance() checks."""
    a = ArgSpec.from_raw({'name': 'x', 'intrinsic': None})
    assert a.intrinsic == ''


def test_from_raw_None_attributes_becomes_empty_list():
    a = ArgSpec.from_raw({'name': 'x', 'attributes': None})
    assert a.attributes == []


def test_from_raw_missing_name_raises():
    """Mirrors the Perl behaviour: every argument must have a name."""
    import pytest
    with pytest.raises(KeyError):
        ArgSpec.from_raw({'intrinsic': 'integer'})


def test_from_raw_does_not_alias_input_attributes_list():
    """Mutating the ArgSpec's attributes must not affect the source dict
    (the ctor copies via `list(d.get('attributes') or [])`)."""
    src = {'name': 'x', 'attributes': ['a', 'b']}
    a = ArgSpec.from_raw(src)
    a.attributes.append('c')
    assert src['attributes'] == ['a', 'b']
