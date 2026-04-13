#!/usr/bin/env python3
"""
Add `units()` method overrides to all concrete nodePropertyExtractor implementations.

For each concrete extractor file, this script:
 1. Parses the class name, base-class type, and function-name prefix.
 2. Inspects the existing `unitsInSI` function body to infer description/quantity.
 3. Injects `procedure :: units => <prefix>Units` into the type definition.
 4. Appends the `<prefix>Units` function to the contains section.

Run from the Galacticus root directory:
    python3 scripts/add_units_overrides.py [--dry-run]
"""

import re
import sys
import os
import argparse

# ---------------------------------------------------------------------------
# Unit-expression → (description, quantity) mapping.
# Keys are normalised Fortran expressions (lower-case, spaces stripped).
# ---------------------------------------------------------------------------
UNIT_MAP = {
    "masssolar":                    ("Solar masses",       "solMass"),
    "megaparsec":                   ("Mpc",                "Mpc"),
    "gigayear":                     ("Gyr",                "Gyr"),
    "kilo":                         ("km/s",               "km/s"),
    "masssolar/gigayear":           ("M\u2609/Gyr",        "solMass/Gyr"),
    "masssolar/megaparsec**3":      ("M\u2609/Mpc\u00b3",  "solMass/Mpc**3"),
    "masssolar/megaparsec**2":      ("M\u2609/Mpc\u00b2",  "solMass/Mpc**2"),
    "ergs":                         ("erg",                "erg"),
    "ergs*gigayear":                ("erg Gyr",            "erg*Gyr"),
    "ergs/centi**2":                ("erg/cm\u00b2",       "erg/cm**2"),
    "ergs*centi**3":                ("erg cm\u00b3",       "erg*cm**3"),
    "kilo*electronvolt":            ("keV",                "keV"),
    "luminosityzeroPointab":        ("????",               "????"),
    "luminosityzeroPointab":        ("????",               "????"),
    "luminositysolar":              ("L\u2609",            "solLum"),
    "1.0d0/megaparsec**3":          ("Mpc\u207b\u00b3",    "Mpc**-3"),
    "1.0d0/gigayear**2":            ("Gyr\u207b\u00b2",    "Gyr**-2"),
    "degreestoradians":             ("degrees",            "deg"),
    "megaparsec/gigayear":          ("Mpc/Gyr",            "Mpc/Gyr"),
    # dimensionless
    "0.0d0":                        ("",                   ""),
    "1.0d0":                        ("",                   ""),
    "0.0_double":                   ("",                   ""),
    "1.0_double":                   ("",                   ""),
}

# Normalise a Fortran expression for lookup.
def normalise(expr):
    return re.sub(r"\s+", "", expr.lower().lstrip("+"))

def lookup_units(expr):
    key = normalise(expr)
    return UNIT_MAP.get(key, ("????", "????"))

# ---------------------------------------------------------------------------
# Detect which abstract base class a type extends.
# ---------------------------------------------------------------------------
BASE_SCALAR          = "scalar"
BASE_TUPLE           = "tuple"
BASE_ARRAY           = "array"
BASE_LIST            = "list"
BASE_LIST2D          = "list2d"
BASE_INT_SCALAR      = "integer_scalar"
BASE_INT_TUPLE       = "integer_tuple"
BASE_INT_LIST        = "integer_list"

BASE_MAP = {
    "nodePropertyExtractorScalar":         BASE_SCALAR,
    "nodePropertyExtractorTuple":          BASE_TUPLE,
    "nodePropertyExtractorArray":          BASE_ARRAY,
    "nodePropertyExtractorList":           BASE_LIST,
    "nodePropertyExtractorList2D":         BASE_LIST2D,
    "nodePropertyExtractorIntegerScalar":  BASE_INT_SCALAR,
    "nodePropertyExtractorIntegerTuple":   BASE_INT_TUPLE,
    "nodePropertyExtractorIntegerList":    BASE_INT_LIST,
}

# ---------------------------------------------------------------------------
# Per-base-class: the template pieces we need.
# ---------------------------------------------------------------------------
# return-type declaration, argument list, and body pattern for unitType construction.

def scalar_units_func(class_name, prefix, desc, qty):
    unused = "!$GLC attributes unused :: self" if not desc.startswith("????") or not qty.startswith("????") else ""
    desc_arg = f",description='{desc}'" if desc else ""
    qty_arg  = f",quantity='{qty}'"     if qty  else ""
    return f"""\
  function {prefix}Units(self) result(units)
    !!{{
    Return the units of the {prefix} property.
    !!}}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType    )                :: units
    class({class_name}), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(self%unitsInSI(){desc_arg}{qty_arg})
    return
  end function {prefix}Units
"""

def array_units_func(class_name, prefix, per_element):
    """per_element: list of (desc, qty) one per element, or a single (desc,qty) for all."""
    if isinstance(per_element, tuple):
        # uniform
        desc, qty = per_element
        desc_arg = f",description='{desc}'" if desc else ""
        qty_arg  = f",quantity='{qty}'"     if qty  else ""
        body = f"""\
    siValues=self%unitsInSI(time)
    allocate(units(size(siValues)))
    do i=1,size(siValues)
       units(i)=unitType(siValues(i){desc_arg}{qty_arg})
    end do"""
    else:
        # per-element list
        lines = ["    siValues=self%unitsInSI(time)",
                 f"    allocate(units({len(per_element)}))"]
        for idx, (d, q) in enumerate(per_element, 1):
            d_arg = f",description='{d}'" if d else ""
            q_arg = f",quantity='{q}'"    if q else ""
            lines.append(f"    units({idx})=unitType(siValues({idx}){d_arg}{q_arg})")
        body = "\n".join(lines)
    return f"""\
  function {prefix}Units(self,time) result(units)
    !!{{
    Return the units of the {prefix} properties.
    !!}}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType    ), dimension(:), allocatable :: units
    class           ({class_name}), intent(inout)             :: self
    double precision              , intent(in   )             :: time
    double precision              , dimension(:), allocatable :: siValues
    integer                                                   :: i
    !$GLC attributes unused :: self

{body}
    return
  end function {prefix}Units
"""

def array_units_func_opt_time(class_name, prefix, per_element):
    """Array base-class uses optional time."""
    if isinstance(per_element, tuple):
        desc, qty = per_element
        desc_arg = f",description='{desc}'" if desc else ""
        qty_arg  = f",quantity='{qty}'"     if qty  else ""
        body = f"""\
    siValues=self%unitsInSI(time)
    allocate(units(size(siValues)))
    do i=1,size(siValues)
       units(i)=unitType(siValues(i){desc_arg}{qty_arg})
    end do"""
    else:
        lines = ["    siValues=self%unitsInSI(time)",
                 f"    allocate(units({len(per_element)}))"]
        for idx, (d, q) in enumerate(per_element, 1):
            d_arg = f",description='{d}'" if d else ""
            q_arg = f",quantity='{q}'"    if q else ""
            lines.append(f"    units({idx})=unitType(siValues({idx}){d_arg}{q_arg})")
        body = "\n".join(lines)
    return f"""\
  function {prefix}Units(self,time) result(units)
    !!{{
    Return the units of the {prefix} properties.
    !!}}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType    ), dimension(:), allocatable :: units
    class           ({class_name}), intent(inout)             :: self
    double precision              , intent(in   ), optional   :: time
    double precision              , dimension(:), allocatable :: siValues
    integer                                                   :: i
    !$GLC attributes unused :: self

{body}
    return
  end function {prefix}Units
"""

def list_units_func(class_name, prefix, per_element):
    """List base-class: no time argument."""
    if isinstance(per_element, tuple):
        desc, qty = per_element
        desc_arg = f",description='{desc}'" if desc else ""
        qty_arg  = f",quantity='{qty}'"     if qty  else ""
        body = f"""\
    siValues=self%unitsInSI()
    allocate(units(size(siValues)))
    do i=1,size(siValues)
       units(i)=unitType(siValues(i){desc_arg}{qty_arg})
    end do"""
    else:
        lines = ["    siValues=self%unitsInSI()",
                 f"    allocate(units({len(per_element)}))"]
        for idx, (d, q) in enumerate(per_element, 1):
            d_arg = f",description='{d}'" if d else ""
            q_arg = f",quantity='{q}'"    if q else ""
            lines.append(f"    units({idx})=unitType(siValues({idx}){d_arg}{q_arg})")
        body = "\n".join(lines)
    return f"""\
  function {prefix}Units(self) result(units)
    !!{{
    Return the units of the {prefix} properties.
    !!}}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType    ), dimension(:), allocatable :: units
    class({class_name}), intent(inout)             :: self
    double precision   , dimension(:), allocatable :: siValues
    integer                                        :: i
    !$GLC attributes unused :: self

{body}
    return
  end function {prefix}Units
"""

def int_scalar_units_func(class_name, prefix, desc, qty):
    desc_arg = f",description='{desc}'" if desc else ""
    qty_arg  = f",quantity='{qty}'"     if qty  else ""
    return f"""\
  function {prefix}Units(self) result(units)
    !!{{
    Return the units of the {prefix} property.
    !!}}
    use :: Units_MetaData, only : unitType
    implicit none
    type (unitType    )                :: units
    class({class_name}), intent(inout) :: self
    !$GLC attributes unused :: self

    units=unitType(self%unitsInSI(){desc_arg}{qty_arg})
    return
  end function {prefix}Units
"""

# ---------------------------------------------------------------------------
# Parse a Fortran file and return information needed to generate the override.
# ---------------------------------------------------------------------------

def parse_file(path):
    with open(path, encoding="utf-8") as fh:
        src = fh.read()

    # Already has a units override?
    if re.search(r"procedure\s*::\s*units\s*=>", src, re.IGNORECASE):
        return None  # skip

    # Find "type, extends(<base>) :: <className>" — use first occurrence
    m = re.search(
        r"type\s*,\s*extends\s*\(\s*(\w+)\s*\)\s*::\s*(\w+)",
        src, re.IGNORECASE)
    if not m:
        return None
    base_raw, class_name = m.group(1), m.group(2)
    base_type = BASE_MAP.get(base_raw)
    if base_type is None:
        return None  # not a recognised extractor base

    # Derive the function-name prefix from the class name.
    # Convention: strip leading "nodePropertyExtractor" and lowercase first char.
    prefix = class_name
    if prefix.lower().startswith("nodepropertyextractor"):
        prefix = prefix[len("nodePropertyExtractor"):]
    prefix = prefix[0].lower() + prefix[1:]

    # Find the unitsInSI function body (everything up to the next "end function").
    units_body = ""
    m2 = re.search(
        r"(function\s+\w+UnitsInSI\b.*?end\s+function\s+\w+UnitsInSI)",
        src, re.IGNORECASE | re.DOTALL)
    if m2:
        units_body = m2.group(1)

    return {
        "src":        src,
        "path":       path,
        "base_type":  base_type,
        "class_name": class_name,
        "prefix":     prefix,
        "units_body": units_body,
    }

# ---------------------------------------------------------------------------
# Infer (description, qty) for scalar and simple uniform-array extractors.
# ---------------------------------------------------------------------------

def infer_scalar(units_body):
    """Return (desc, qty) from a unitsInSI function body for a scalar extractor."""
    # Extract the assignment: <name>UnitsInSI = <expr>
    m = re.search(r"\w+UnitsInSI\s*=\s*(.+)", units_body)
    if not m:
        return ("????", "????")
    expr = m.group(1).strip().rstrip("&").strip()
    # Strip line-continuation and collapse whitespace.
    expr = re.sub(r"&\s*\n\s*&?", " ", expr)
    expr = re.sub(r"\s+", "", expr)
    return lookup_units(expr)

# ---------------------------------------------------------------------------
# Infer per-element (desc, qty) list for array/tuple/list extractors.
# ---------------------------------------------------------------------------

def infer_array_elements(units_body, n_hint=None):
    """
    Try to extract per-element units from a unitsInSI function body.
    Handles:
      - Array literal:   name = [a, b, c]  or  name = (/ a, b, c /)
      - Element-wise:    name(1)=a  /  name(2)=b  ...
      - Uniform scalar:  name = expr
    Returns either a uniform (desc,qty) tuple if all elements agree, or a list.
    """
    # Flatten continuation lines (& at end-of-line, optional & at start of next)
    flat = re.sub(r"&\s*\n\s*&?", " ", units_body)

    # --- 1. Array literal assignment -----------------------------------------
    # Match bare "unitsInSI" or any prefix like "fooUnitsInSI".
    # Anchor to start-of-line to avoid false matches on allocate(name(n)).
    m = re.search(r"^\s*\w*[Uu]nitsInSI\s*=\s*\[(.+?)\]",
                  flat, re.MULTILINE | re.DOTALL)
    if not m:
        m = re.search(r"^\s*\w*[Uu]nitsInSI\s*=\s*\(/(.+?)/\)",
                      flat, re.MULTILINE | re.DOTALL)
    if m:
        items_str = m.group(1)
        items = [x.strip().lstrip("+").rstrip("+").strip()
                 for x in re.split(r",", items_str)]
        items = [x for x in items if x]
        result = [lookup_units(x) for x in items]
        if len(set(result)) == 1:
            return result[0]
        return result

    # --- 2. Element-wise indexed assignments: name(i) = expr -----------------
    indexed = re.findall(
        r"^\s*\w*[Uu]nitsInSI\s*\(\d+\)\s*=\s*([^!\n]+)",
        flat, re.MULTILINE)
    if indexed:
        result = [lookup_units(re.sub(r"\s+", "", e.strip().rstrip("&")))
                  for e in indexed]
        if len(set(result)) == 1:
            return result[0]
        return result

    # --- 3. Uniform scalar assignment -----------------------------------------
    m2 = re.search(r"^\s*\w*[Uu]nitsInSI\s*=\s*([^!\n]+)", flat, re.MULTILINE)
    if m2:
        expr = re.sub(r"\s+", "", m2.group(1).strip().rstrip("&"))
        return lookup_units(expr)

    return ("????", "????")

# ---------------------------------------------------------------------------
# Generate the override function text and the updated source.
# ---------------------------------------------------------------------------

def make_override(info):
    base  = info["base_type"]
    cn    = info["class_name"]
    pfx   = info["prefix"]
    body  = info["units_body"]

    if base == BASE_SCALAR or base == BASE_INT_SCALAR:
        desc, qty = infer_scalar(body)
        if base == BASE_SCALAR:
            fn_text = scalar_units_func(cn, pfx, desc, qty)
        else:
            fn_text = int_scalar_units_func(cn, pfx, desc, qty)

    elif base == BASE_TUPLE or base == BASE_INT_TUPLE:
        per_elem = infer_array_elements(body)
        fn_text = array_units_func(cn, pfx, per_elem)

    elif base == BASE_ARRAY:
        per_elem = infer_array_elements(body)
        fn_text = array_units_func_opt_time(cn, pfx, per_elem)

    elif base in (BASE_LIST, BASE_INT_LIST, BASE_LIST2D):
        per_elem = infer_array_elements(body)
        fn_text = list_units_func(cn, pfx, per_elem)

    else:
        return None

    return fn_text

def inject(info, fn_text, dry_run=False):
    src = info["src"]
    pfx = info["prefix"]
    path = info["path"]

    # 1. Insert "procedure :: units => <pfx>Units" after the last existing
    #    "procedure :: " line in the type definition block.
    #    We look for the last occurrence of "procedure ::" before "end type".
    end_type_pat = re.compile(
        r"(end\s+type\s+" + re.escape(info["class_name"]) + r")", re.IGNORECASE)
    m = end_type_pat.search(src)
    if not m:
        print(f"  SKIP (no end type found): {path}")
        return

    # Find the last "procedure ::" before end type
    before_end = src[:m.start()]
    last_proc = list(re.finditer(r"([ \t]*procedure\b[^\n]*\n)", before_end, re.IGNORECASE))
    if not last_proc:
        print(f"  SKIP (no procedure lines): {path}")
        return

    insert_after = last_proc[-1].end()
    new_line = f"     procedure :: units       => {pfx}Units\n"
    src = src[:insert_after] + new_line + src[insert_after:]

    # 2. Append the function before the final "end" line of the file
    #    (just before the last blank-or-end line), or at end of file.
    src = src.rstrip() + "\n\n" + fn_text

    if dry_run:
        print(f"  [dry-run] would update {path}")
        return

    with open(path, "w", encoding="utf-8") as fh:
        fh.write(src)
    print(f"  updated {path}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

ABSTRACT_FILES = {
    "nodes.property_extractor.scalar.F90",
    "nodes.property_extractor.tuple.F90",
    "nodes.property_extractor.array.F90",
    "nodes.property_extractor.list.F90",
    "nodes.property_extractor.list2D.F90",
    "nodes.property_extractor.integer_scalar.F90",
    "nodes.property_extractor.integer_tuple.F90",
    "nodes.property_extractor.integer_list.F90",
    "nodes.property_extractor.multi.F90",
    "nodes.property_extractor.F90",
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be done without modifying files")
    parser.add_argument("source_dir", nargs="?", default="source",
                        help="Directory containing the .F90 files (default: source)")
    args = parser.parse_args()

    source_dir = args.source_dir
    files = sorted(f for f in os.listdir(source_dir)
                   if f.startswith("nodes.property_extractor.") and f.endswith(".F90")
                   and f not in ABSTRACT_FILES)

    updated = skipped = errors = 0
    for fname in files:
        path = os.path.join(source_dir, fname)
        try:
            info = parse_file(path)
            if info is None:
                skipped += 1
                continue
            fn_text = make_override(info)
            if fn_text is None:
                skipped += 1
                continue
            inject(info, fn_text, dry_run=args.dry_run)
            updated += 1
        except Exception as exc:
            print(f"  ERROR {path}: {exc}")
            errors += 1

    print(f"\nDone: {updated} updated, {skipped} skipped, {errors} errors.")

if __name__ == "__main__":
    main()
