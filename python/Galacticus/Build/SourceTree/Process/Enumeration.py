"""Processes `enumeration` directives: for each enum declared under a module
or file, synthesizes the enumeration type, its equality operator, and
(optionally) validator / encode / decode / describe helpers.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/SourceTree/Process/Enumeration.pm
"""

from List.ExtraUtils                                         import as_array
from Galacticus.Build.SourceTree                             import (
    walk_tree, parse_code, children, set_visibility,
    insert_after_node, insert_pre_contains, insert_post_contains,
)
from Galacticus.Build.SourceTree.Process                     import register_process
from Galacticus.Build.SourceTree.Parse.ModuleUses            import add_uses
from Galacticus.Build.SourceTree.Process.SourceIntrospection import location


def _ucfirst(s):
    return s[:1].upper() + s[1:] if s else s


def _insert_code_tree(parent, source_text, inserter):
    """Parse `source_text` into a subtree and insert its direct children into
    `parent` via `inserter` (one of insert_after_node / insert_pre_contains /
    insert_post_contains).

    Mirrors the Perl idiom
        my $sub  = ParseCode($text, '…');
        my @kids = Children($sub);
        Insert…($target, \\@kids);
    used throughout Enumeration.pm.

    The synthetic sub-tree is parsed with the OUTER source's name as its
    `source` attribution.  Keeping the outer file's name on every node
    means that when FunctionClass later threads a child class file's
    enumerations into a parent's `<module>` interfaces list, the parent's
    `classDocumentation` postprocess pass sees `source != outer_source`
    on those nodes and correctly skips them — letting the child file's
    own preprocess.py invocation own its enumeration's documentation.
    A bare `'Enumeration'` placeholder source would always count as
    "not a real .F90 file" and so always get processed in the parent's
    pass — leading to "missing function type for method 'operator(==)'"
    warnings for every enumeration declared in a functionClass child
    file (treeStatistic, wagner2016*, …).
    """
    src_name = parent.get('source') or 'Enumeration'
    sub_tree = parse_code(source_text, name=src_name)
    kids = children(sub_tree)
    for k in kids:
        k['parent'] = None
    inserter(parent, kids)


# ---------------------------------------------------------------------------
# Generated-code builders (one per Perl sub-block)
# ---------------------------------------------------------------------------

def _emit_type_definition(name, entries, indexing, visibility, validator):
    """Emit the `type, extends(enumerationType) :: enumerationXType` block,
    the per-member `parameter` values, and optionally the Min/Max/Count.

    Mirrors Enumeration.pm:36-58.
    """
    lines  = "  ! Auto-generated enumeration\n"
    lines += f"  type, extends(enumerationType) :: enumeration{name}Type\n"
    lines += "  contains\n"
    lines += "    !![\n"
    lines += "    <methods>\n"
    lines += (
        "      <method method=\"operator(==)\" description=\"Test the "
        "equality of two members of the enumeration.\"/>\n"
    )
    lines += "    </methods>\n"
    lines += "    !!]\n"
    lines += f"    procedure ::                  enumeration{name}IsEqual\n"
    lines += f"    generic   :: operator(==) =>  enumeration{name}IsEqual\n"
    lines += f"  end type enumeration{name}Type\n"

    i = indexing - 1
    for entry in entries:
        i += 1
        label = entry.get('label', '')
        lines += (
            f"  type(enumeration{name}Type), parameter, {visibility} :: "
            f"{name}{_ucfirst(label)}=enumeration{name}Type({i})\n"
        )

    count = i + 1 - indexing
    if validator == 'yes':
        lines += f"  integer, parameter, {visibility} :: {name}Min  ={indexing}\n"
        lines += f"  integer, parameter, {visibility} :: {name}Max  ={i}\n"
        lines += f"  integer, parameter, {visibility} :: {name}Count={count}\n"

    lines += "  ! End auto-generated enumeration\n\n"
    return lines


def _emit_equality_function(name):
    """Emit the `enumerationXIsEqual` pure elemental function.  Mirrors
    Enumeration.pm:78-91.
    """
    fn = f"enumeration{_ucfirst(name)}IsEqual"
    out  = "\n"
    out += "  ! Auto-generated enumeration function\n"
    out += (
        f"  pure elemental logical function {fn}(enumerationA,enumerationB) "
        "result(isEqual)\n"
    )
    out += "    !!{\n"
    out += f"    Validate a \\mono{{{name}}} enumeration value.\n"
    out += "    !!}\n"
    out += "    implicit none\n\n"
    out += (
        f"    class(enumeration{name}Type), intent(in   ) :: "
        "enumerationA, enumerationB\n"
    )
    out += "    isEqual=enumerationA%ID == enumerationB%ID\n"
    out += "    return\n"
    out += f"  end function {fn}\n"
    out += "  ! End auto-generated enumeration function\n"
    return out, fn


def _emit_validator_function(name):
    """Emit the `enumerationXIsValid` validator.  Mirrors Enumeration.pm:98-112."""
    fn = f"enumeration{_ucfirst(name)}IsValid"
    out  = "\n"
    out += "  ! Auto-generated enumeration function\n"
    out += f"  logical function {fn}(enumerationValue)\n"
    out += "    !!{\n"
    out += f"    Validate a \\mono{{{name}}} enumeration value.\n"
    out += "    !!}\n"
    out += "    implicit none\n\n"
    out += f"    type(enumeration{name}Type), intent(in   ) :: enumerationValue\n"
    out += (
        f"    {fn}=(enumerationValue%ID >= {name}Min .and. "
        f"enumerationValue%ID <= {name}Max)\n"
    )
    out += "    return\n"
    out += f"  end function {fn}\n"
    out += "  ! End auto-generated enumeration function\n"
    return out, fn


def _emit_encode_function(name, entries, indexing, on_error, node):
    """Emit the four-overload encode function family + its two interfaces.

    Mirrors Enumeration.pm:121-239.  Returns (function_text, interface_text,
    public_function_name).
    """
    fn = f"enumeration{_ucfirst(name)}Encode"
    status_arg = "" if on_error else ",status"

    interface  = f" interface {fn}\n"
    interface += f"  module procedure {fn}Char\n"
    interface += f"  module procedure {fn}VarStr\n"
    interface += f" end interface {fn}\n\n"
    interface += f" interface {fn}ID\n"
    interface += f"  module procedure {fn}IDChar\n"
    interface += f"  module procedure {fn}IDVarStr\n"
    interface += f" end interface {fn}ID\n\n"

    out  = "\n"
    out += "  ! Auto-generated enumeration functions\n"

    # -- IDVarStr --
    out += f"  integer function {fn}IDVarStr(name,includesPrefix)\n"
    out += "    !!{\n"
    out += (
        f"    Encode a \\mono{{{name}}} enumeration from a string, returning "
        "the appropriate identifier ID.\n"
    )
    out += "    !!}\n"
    out += "    use :: ISO_Varying_String\n"
    out += "    implicit none\n\n"
    out += "    type   (varying_string), intent(in   )           :: name\n"
    out += "    logical                , intent(in   ), optional :: includesPrefix\n"
    out += f"    {fn}IDVarStr={fn}ID(char(name),includesPrefix)\n"
    out += "    return\n"
    out += f"  end function {fn}IDVarStr\n\n"

    # -- IDChar --
    out += f"  integer function {fn}IDChar(name,includesPrefix)\n"
    out += "    !!{\n"
    out += (
        f"    Encode a \\mono{{{name}}} enumeration from a string, returning "
        "the appropriate identifier ID.\n"
    )
    out += "    !!}\n"
    out += "    use :: ISO_Varying_String\n"
    out += "    implicit none\n\n"
    out += "    character(len=*), intent(in   )           :: name\n"
    out += "    logical         , intent(in   ), optional :: includesPrefix\n"
    out += f"    type(enumeration{name}Type) :: member\n"
    out += f"    member={fn}(name,includesPrefix)\n"
    out += f"    {fn}IDChar=member%ID\n"
    out += "    return\n"
    out += f"  end function {fn}IDChar\n\n"

    # -- VarStr --
    out += f"  function {fn}VarStr(name,includesPrefix{status_arg})\n"
    out += "    !!{\n"
    out += (
        f"    Encode a \\mono{{{name}}} enumeration from a string, returning "
        "the appropriate identifier.\n"
    )
    out += "    !!}\n"
    out += "    use :: ISO_Varying_String\n"
    out += "    implicit none\n\n"
    out += f"    type   (enumeration{name}Type) :: {fn}VarStr\n"
    out += "    type   (varying_string), intent(in   )           :: name\n"
    out += "    logical                , intent(in   ), optional :: includesPrefix\n"
    if not on_error:
        out += "    integer                , intent(  out), optional :: status\n"
    out += f"    {fn}VarStr={fn}(char(name),includesPrefix{status_arg})\n"
    out += "    return\n"
    out += f"  end function {fn}VarStr\n\n"

    # -- Char (the workhorse) --
    out += f"  function {fn}Char(name,includesPrefix{status_arg})\n"
    out += "    !!{\n"
    out += (
        f"    Encode a \\mono{{{name}}} enumeration from a string, returning "
        "the appropriate identifier.\n"
    )
    out += "    !!}\n"
    if not on_error:
        out += (
            "    use :: Error             , only : Error_Report, "
            "errorStatusSuccess, errorStatusFail\n"
        )
    out += "    use :: ISO_Varying_String, only : var_str     , operator(//)\n"
    out += "    implicit none\n\n"
    out += f"    type   (enumeration{name}Type) :: {fn}Char\n"
    out += "    character(len=*), intent(in   )           :: name\n"
    out += "    logical         , intent(in   ), optional :: includesPrefix\n"
    if not on_error:
        out += "    integer         , intent(  out), optional :: status\n"
    out += "    logical                                   :: includesPrefix_\n\n"
    out += "    includesPrefix_=.true.\n"
    out += "    if (present(includesPrefix)) includesPrefix_=includesPrefix\n"
    if not on_error:
        out += "    if (present(status)) status=errorStatusSuccess\n"
    for j in range(2):
        if j == 0:
            out += "    if (includesPrefix_) then\n"
        else:
            out += "    else\n"
        out += "      select case (trim(name))\n"
        i = indexing - 1
        for entry in entries:
            i += 1
            label = entry.get('label', '')
            case_text = f"{name}{_ucfirst(label)}" if j == 0 else label
            out += f"      case ('{case_text}')\n"
            out += f"        {fn}Char=enumeration{name}Type({i})\n"
        out += "      case default\n"
        if on_error:
            out += f"      {fn}Char=enumeration{name}Type({on_error})\n"
        else:
            out += f"      {fn}Char=enumeration{name}Type(-1)\n"
            out += "      if (present(status)) then\n"
            out += "         status=errorStatusFail\n"
            out += "      else\n"
            loc_expr = location(node, node.get('line', 0))
            out += (
                "         call Error_Report(var_str('unrecognized enumeration "
                "member [')//trim(name)//']'//"
                f"enumeration{_ucfirst(name)}Describe()//"
                f"{loc_expr})\n"
            )
            out += "      end if\n"
        out += "      end select\n"
    out += "    end if\n"
    out += "    return\n"
    out += f"  end function {fn}Char\n\n"
    out += "  ! End auto-generated enumeration functions\n"
    return out, interface, fn


def _emit_decode_function(name, entries, indexing, on_error, node):
    """Emit the decode function family + its interface.

    Mirrors Enumeration.pm:243-315.
    """
    fn = f"enumeration{_ucfirst(name)}Decode"

    interface  = f" interface {fn}\n"
    interface += f"  module procedure {fn}Enumerator\n"
    interface += f"  module procedure {fn}ID\n"
    interface += f" end interface {fn}\n\n"

    out  = "\n"
    out += "  ! Auto-generated enumeration function\n"

    # -- Enumerator wrapper --
    out += f"  function {fn}Enumerator(enumerationValue,includePrefix)\n"
    out += "    !!{\n"
    out += f"    Decode a \\mono{{{name}}} enumeration to a string.\n"
    out += "    !!}\n"
    out += "    use ISO_Varying_String\n"
    out += "    implicit none\n\n"
    out += (
        f"    type   (varying_string)                                        "
        f"                     :: {fn}Enumerator\n"
    )
    out += (
        f"    type   (enumeration{name}Type), intent(in   )           :: "
        "enumerationValue\n"
    )
    out += (
        "    logical                                                   , "
        "intent(in   ), optional :: includePrefix\n\n"
    )
    out += f"    {fn}Enumerator={fn}(enumerationValue%ID,includePrefix)\n"
    out += "    return\n"
    out += f"  end function {fn}Enumerator\n"

    # -- ID worker --
    out += f"  function {fn}ID(enumerationValue,includePrefix)\n"
    out += "    !!{\n"
    out += f"    Decode a \\mono{{{name}}} enumeration to a string.\n"
    out += "    !!}\n"
    out += "    use ISO_Varying_String\n"
    if not on_error:
        out += "    use Error\n"
    out += "    implicit none\n\n"
    out += f"    type   (varying_string)                          :: {fn}ID\n"
    out += "    integer                , intent(in   )           :: enumerationValue\n"
    out += "    logical                , intent(in   ), optional :: includePrefix\n"
    for j in range(2):
        if j == 0:
            out += "    if (present(includePrefix).and.includePrefix) then\n"
            out += f"      {fn}ID='{name}'\n"
        else:
            out += "    else\n"
            out += f"      {fn}ID=''\n"
        i = indexing - 1
        out += "    select case(enumerationValue)\n"
        for entry in entries:
            i += 1
            label = entry.get('label', '')
            label_text = _ucfirst(label) if j == 0 else label
            out += f"    case ({i})\n"
            out += f"       {fn}ID={fn}ID//'{label_text}'\n"
        out += "    case default\n"
        if on_error:
            out += f"      {fn}ID={fn}ID//'Error'\n"
        else:
            loc_expr = location(node, node.get('line', 0))
            out += (
                f"      call Error_Report('invalid enumeration value'//"
                f"{loc_expr})\n"
            )
        out += "    end select\n"
    out += "    end if\n"
    out += "    return\n"
    out += f"  end function {fn}ID\n"
    out += "  ! End auto-generated enumeration function\n"
    return out, interface, fn


def _emit_description_function(name, entries, indexing):
    """Emit `enumerationXDescription{Enumerator,ID}` + its interface.

    Mirrors Enumeration.pm:319-385.  The body of the ID function includes a
    `select case` whose text is built separately in Perl and spliced via
    `InsertAfterNode`; here we concatenate it directly into the source text
    before re-parsing.
    """
    fn = f"enumeration{_ucfirst(name)}Description"

    interface  = f" interface {fn}\n"
    interface += f"  module procedure {fn}Enumerator\n"
    interface += f"  module procedure {fn}ID\n"
    interface += f" end interface {fn}\n\n"

    wrapper  = "\n"
    wrapper += "  ! Auto-generated enumeration function\n"
    wrapper += f"  function {fn}Enumerator(enumerationValue)\n"
    wrapper += "    !!{\n"
    wrapper += (
        f"    Return a description of a \\mono{{{name}}} enumeration member.\n"
    )
    wrapper += "    !!}\n"
    wrapper += "    use ISO_Varying_String\n"
    wrapper += "    implicit none\n\n"
    wrapper += (
        f"    type   (varying_string)                                        "
        f"                     :: {fn}Enumerator\n"
    )
    wrapper += (
        f"    type   (enumeration{name}Type), intent(in   )           :: "
        "enumerationValue\n\n"
    )
    wrapper += f"    {fn}Enumerator={fn}(enumerationValue%ID)\n"
    wrapper += "    return\n"
    wrapper += f"  end function {fn}Enumerator\n"

    body  = f"  function {fn}ID(enumerationValue) result(description)\n"
    body += "    !!{\n"
    body += f"    Return a description of a \\mono{{{name}}} enumeration value.\n"
    body += "    !!}\n"
    body += "    use :: ISO_Varying_String, only : varying_string, assignment(=)\n"
    body += "    implicit none\n"
    body += "    type(varying_string) :: description\n"
    body += "    integer, intent(in) :: enumerationValue\n\n"
    body += "    select case (enumerationValue)\n"
    i = indexing - 1
    for entry in entries:
        i += 1
        desc = entry.get('description', '') or ''
        body += f"   case ({i})\n"
        body += f"    description='{desc}'\n"
    body += "    end select\n"
    body += "    return\n"
    body += f"  end function {fn}ID\n"
    body += "  ! End auto-generated enumeration function\n"
    return wrapper, body, interface, fn


def _emit_describe_function(name, entries):
    """Emit the unconditional `enumerationXDescribe` function.  Mirrors
    Enumeration.pm:387-431.
    """
    fn = f"enumeration{_ucfirst(name)}Describe"

    out  = "\n"
    out += "  ! Auto-generated enumeration function\n"
    out += f"  function {fn}() result(description)\n"
    out += "    !!{\n"
    out += f"    Return a description of the \\mono{{{name}}} enumeration.\n"
    out += "    !!}\n"
    out += (
        "    use :: ISO_Varying_String, only : varying_string, var_str, "
        "operator(//)\n"
    )
    out += "    implicit none\n"
    out += "    type(varying_string) :: description\n\n"
    out += (
        f"    description=var_str(char(10))//\"Enumeration '{name}' has the "
        "following members:\"\n"
    )
    max_len = max((len(e.get('label', '')) for e in entries), default=0)
    last = len(entries) - 1
    for i, entry in enumerate(entries):
        label = entry.get('label', '')
        desc  = entry.get('description')
        sep   = '.' if i == last else ';'
        pad   = ' ' * (max_len - len(label))
        body  = f": {desc}{sep}" if desc else ''
        out += (
            f"    description=description//char(10)//\"   {pad}{label}{body}\"\n"
        )
    out += "    \n"
    out += "    return\n"
    out += f"  end function {fn}\n"
    out += "  ! End auto-generated enumeration function\n"
    return out, fn


# ---------------------------------------------------------------------------
# Main pass
# ---------------------------------------------------------------------------

def process_enumerations(tree, options):
    """Mirrors Process_Enumerations() from Enumeration.pm."""
    for node in walk_tree(tree):
        if node.get('type') != 'enumeration':
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue

        parent = node.get('parent') or {}
        if parent.get('type') not in ('module', 'file'):
            raise RuntimeError(
                "process_enumerations: parent node must be a module or file")

        directive['processed'] = True
        name        = directive['name']
        visibility  = directive.get('visibility', 'public')
        validator   = directive.get('validator',  'no')
        indexing    = int(directive.get('indexing', 0))
        entries     = list(as_array(directive.get('entry')))
        encode_on   = directive.get('encodeFunction') == 'yes'
        decode_on   = directive.get('decodeFunction') == 'yes'
        on_error    = directive.get('errorValue')  # optional

        # --- Type definition (+ optional Min/Max/Count params) ---
        type_text = _emit_type_definition(
            name, entries, indexing, visibility, validator)
        add_uses(parent, {
            'moduleUse':   {'Enumerations': {
                'intrinsic': False, 'only': {'enumerationType': True}}},
            'moduleOrder': ['Enumerations'],
        })
        set_visibility(parent, f"enumeration{name}Type", visibility)
        _insert_code_tree(
            parent, type_text,
            inserter=lambda p, kids: insert_after_node(node, kids))

        # --- Equality function ---
        eq_text, _ = _emit_equality_function(name)
        _insert_code_tree(parent, eq_text, inserter=insert_post_contains)

        # --- Validator (optional) ---
        if validator == 'yes':
            val_text, val_name = _emit_validator_function(name)
            _insert_code_tree(parent, val_text, inserter=insert_post_contains)
            set_visibility(parent, val_name, visibility)

        # --- Encode (optional) ---
        if encode_on:
            enc_text, enc_iface, enc_name = _emit_encode_function(
                name, entries, indexing, on_error, node)
            _insert_code_tree(parent, enc_text, inserter=insert_post_contains)
            _insert_code_tree(parent, enc_iface, inserter=insert_pre_contains)
            set_visibility(parent, enc_name, visibility)

        # --- Decode + Description (paired in Perl, toggled by decodeFunction) ---
        if decode_on:
            dec_text, dec_iface, dec_name = _emit_decode_function(
                name, entries, indexing, on_error, node)
            _insert_code_tree(parent, dec_text, inserter=insert_post_contains)
            _insert_code_tree(parent, dec_iface, inserter=insert_pre_contains)
            set_visibility(parent, dec_name, visibility)

            wrapper, body, desc_iface, desc_name = _emit_description_function(
                name, entries, indexing)
            _insert_code_tree(parent, body,    inserter=insert_post_contains)
            _insert_code_tree(parent, wrapper, inserter=insert_post_contains)
            _insert_code_tree(parent, desc_iface, inserter=insert_pre_contains)
            set_visibility(parent, desc_name, visibility)

        # --- Describe (always) ---
        desc_text, desc_name = _emit_describe_function(name, entries)
        _insert_code_tree(parent, desc_text, inserter=insert_post_contains)
        set_visibility(parent, desc_name, visibility)


register_process('enumerations', process_enumerations)
