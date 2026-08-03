"""Contains a Python module which implements parsing of directives in the Galacticus preprocessor system.

Andrew Benson (ported to Python 2026)
"""

import re
import os
import xml.etree.ElementTree as ET


from XML.Utils import xml_to_dict

try:
    from lxml import etree as _lxml_etree  # type: ignore[import-not-found]
    _HAS_LXML = True
except ImportError:
    _lxml_etree = None
    _HAS_LXML = False

# Cached stateStorables.xml content, loaded lazily from $BUILDPATH.
_state_storables = None
_state_storables_loaded = False

# Define the directory where schemas are stored
EXEC_PATH = os.environ.get('GALACTICUS_EXEC_PATH')
if EXEC_PATH:
    SCHEMAS_DIR = os.path.abspath(os.path.join(EXEC_PATH,"/schemas"))
else:
    SCHEMAS_DIR = None

if _HAS_LXML:
    class SchemaResolver(_lxml_etree.Resolver):
        """Custom resolver to set the absolute path for schema imports.

        This is needed as we have schemas that are constructed internally from
        text, so contain no path information when parsed.
        """
        def resolve(self, system_url, public_id, context):
            # Check if the import path refers to the file we want to redirect
            if 'commonTypes.xsd' in system_url:
                # Construct the absolute path to the file
                full_path = os.path.join(SCHEMAS_DIR, 'commonTypes.xsd')
                # Return the resolved file path to lxml
                return self.resolve_filename(full_path, context)
            return None

def _load_state_storables():
    """Load $BUILDPATH/stateStorables.xml once and cache it.

    Returns None if BUILDPATH is unset or the file is missing.
    """
    global _state_storables, _state_storables_loaded
    if _state_storables_loaded:
        return _state_storables
    _state_storables_loaded = True
    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        return None
    path = os.path.join(build_path, 'stateStorables.xml')
    if not os.path.exists(path):
        return None
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return None
    _state_storables = xml_to_dict(root)
    return _state_storables


# XSD schema templates.  `{name}` is filled in with the directive name.
_FUNCTION_CLASS_SCHEMA = """<?xml version="1.0"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xs:element name="{name}">
    <xs:complexType>
      <xs:sequence>
       <xs:element name="description"       type="xs:string" minOccurs="1" maxOccurs="1"/>
       <xs:element name="descriptorSpecial" type="xs:string" minOccurs="0" maxOccurs="1"/>
       <xs:element name="linkedList"                         minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:attribute name="type"       use="required"/>
         <xs:attribute name="variable"   use="required"/>
         <xs:attribute name="next"       use="required"/>
         <xs:attribute name="object"     use="required"/>
         <xs:attribute name="objectType" use="required"/>
         <xs:attribute name="module"     use="optional"/>
        </xs:complexType>
       </xs:element>
       <xs:element name="deepCopy"                     minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:sequence>
          <xs:element name="ignore"        minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="functionClass" minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="increment"     minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
            <xs:attribute name="atomic"    use="optional" >
             <xs:simpleType>
              <xs:restriction base="xs:string">
               <xs:enumeration value="no" />
               <xs:enumeration value="yes"/>
              </xs:restriction>
             </xs:simpleType>
            </xs:attribute>
           </xs:complexType>
          </xs:element>
          <xs:element name="setTo"         minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
            <xs:attribute name="value"     use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="deallocate"    minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
         </xs:sequence>
        </xs:complexType>
       </xs:element>
       <xs:element name="assignment"                   minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:sequence>
          <xs:element name="functionClass" minOccurs="0" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
         </xs:sequence>
         <xs:attribute name="forceArrayAssign" use="optional"/>
        </xs:complexType>
       </xs:element>
       <xs:element name="stateStorable"                minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:sequence>
          <xs:element name="functionClass" minOccurs="0" maxOccurs="1"         >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="restoreTo"     minOccurs="0" maxOccurs="unbounded" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
            <xs:attribute name="state"     use="required"/>
           </xs:complexType>
          </xs:element>
          <xs:element name="exclude"       minOccurs="0" maxOccurs="1"         >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
           </xs:complexType>
          </xs:element>
         </xs:sequence>
        </xs:complexType>
       </xs:element>
       <xs:element name="stateStore"                   minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:sequence>
          <xs:element name="stateStore" minOccurs="1" maxOccurs="1" >
           <xs:complexType>
            <xs:attribute name="variables" use="required"/>
            <xs:attribute name="store"     use="required"/>
            <xs:attribute name="restore"   use="required"/>
            <xs:attribute name="module"    use="required"/>
           </xs:complexType>
          </xs:element>
         </xs:sequence>
        </xs:complexType>
       </xs:element>
       <xs:element name="runTimeFileDependencies"      minOccurs="0" maxOccurs="1" >
        <xs:complexType>
         <xs:attribute name="paths" use="required"/>
        </xs:complexType>
       </xs:element>
      </xs:sequence>
      <xs:attribute name="name"      use="required"/>
      <xs:attribute name="recursive" use="optional" >
       <xs:simpleType>
        <xs:restriction base="xs:string">
         <xs:enumeration value="no" />
         <xs:enumeration value="yes"/>
        </xs:restriction>
       </xs:simpleType>
      </xs:attribute>
      <xs:attribute name="abstract"  use="optional" >
       <xs:simpleType>
        <xs:restriction base="xs:string">
         <xs:enumeration value="no" />
         <xs:enumeration value="yes"/>
        </xs:restriction>
       </xs:simpleType>
      </xs:attribute>
      <!-- docformat="rst" marks the description as reStructuredText, rendered
           on ReadTheDocs (see scripts/doc/extractDocsRST.py). -->
      <xs:attribute name="docformat" use="optional"/>
    </xs:complexType>
  </xs:element>
</xs:schema>
"""

_EVENT_HOOK_STATIC_SCHEMA = """<?xml version="1.0"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xs:simpleType name="yesNo">
    <xs:restriction base="xs:string">
      <xs:enumeration value="yes"/>
      <xs:enumeration value="no" />
    </xs:restriction>
  </xs:simpleType>
  <xs:element name="{name}">
    <xs:complexType>
      <xs:attribute name="function"  use="required"             />
      <xs:attribute name="after"     use="optional"             />
      <xs:attribute name="before"    use="optional"             />
      <xs:attribute name="useGlobal" use="optional" type="yesNo"/>
    </xs:complexType>
  </xs:element>
</xs:schema>
"""


def _validate_directive(directive_name, xml_text, context_node, line_number):
    """Validate the directive XML against the appropriate XSD schema if possible.

    If lxml is not installed, or no schema applies, validation is skipped
    silently.
    """
    if not _HAS_LXML:
        return
    state_storables = _load_state_storables()
    schema_xml = None

    from Galacticus.Build.StateStorables import (
        function_class_names    as _fcn_names,
        event_hook_static_names as _ehs_names,
    )
    function_class_names = _fcn_names(state_storables)
    event_hook_statics   = _ehs_names(state_storables)

    if state_storables and (directive_name + "Class") in function_class_names:
        schema_xml = _FUNCTION_CLASS_SCHEMA.format(name=directive_name)
    elif state_storables and directive_name in event_hook_statics:
        schema_xml = _EVENT_HOOK_STATIC_SCHEMA.format(name=directive_name)
    else:
        if SCHEMAS_DIR:
            candidate = os.path.join(SCHEMAS_DIR, directive_name + '.xsd')
            if os.path.exists(candidate):
                with open(candidate, 'r') as fh:
                    schema_xml = fh.read()

    if schema_xml is None:
        return

    try:
        parser = _lxml_etree.XMLParser(remove_blank_text=True)
        parser.resolvers.add(SchemaResolver())
        schema_doc = _lxml_etree.fromstring(schema_xml.encode('utf-8'), parser=parser)
        schema = _lxml_etree.XMLSchema(schema_doc)
        document = _lxml_etree.fromstring(xml_text.encode('utf-8'))
    except _lxml_etree.XMLSyntaxError as exc:
        raise RuntimeError(
            f"Parse_Directives: failed parsing directive '{directive_name}': {exc}")
    if not schema.validate(document):
        # Walk up to the enclosing file node for the error message.
        file_node = context_node
        while file_node is not None and file_node.get('type') != 'file':
            file_node = file_node.get('parent')
        file_name = file_node.get('name') if file_node else 'unknown'
        raise RuntimeError(
            f"Parse_Directives: validation of directive '{directive_name}' failed "
            f"in file {file_name} at line {line_number}:\n"
            f"{schema.error_log}")


def _parse_directive_xml(xml_text, context_node, line_number):
    """Parse accumulated XML text into a directive node dict.

    `line_number` is the absolute (1-based) line number, in the original
    source file, of the directive's first line — it becomes the directive
    node's `line`, which downstream process hooks (ObjectBuilder, …) embed
    in `{introspection:location}` expansions.

    Returns None if the text does not parse as XML.
    """
    try:
        elem = ET.fromstring(xml_text)
    except ET.ParseError:
        # Wrap in a root element and try again — handles bare fragments
        # (e.g. leading whitespace before a single directive).  Only accept
        # the wrap when it produces exactly one child; otherwise the input
        # is multiple sibling directives, which the caller should have
        # split into individual directives upstream.
        try:
            elem = ET.fromstring('<root>' + xml_text + '</root>')
        except ET.ParseError:
            return None
        if len(elem) != 1:
            return None
        elem = list(elem)[0]

    directive_name = elem.tag
    _validate_directive(directive_name, xml_text, context_node, line_number)
    directive_dict = xml_to_dict(elem)
    return {
        'type':       directive_name,
        'directive':  directive_dict,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':     context_node.get('source', 'unknown'),
        'line':       line_number,
    }


def parse_directives(tree):
    """Walk the tree replacing XML directive comment blocks with directive nodes.

    Directives are delimited by:
      !![          (opening marker)
      !< <tagname ...>  (XML content lines — '!<' prefix stripped)
      !!]          (closing marker)

    A single `!![ ... !!]` block may contain *multiple* sibling directive
    tags (e.g. six `<constant ... />` in a row).  Each is emitted as its
    own directive node, wrapped in synthetic `!![ ... !!]` markers — this
    prevents the XML parser from collapsing the siblings into a synthetic
    `<root>` wrapper node.
    """
    from Galacticus.Build.SourceTree import walk_tree, replace_node, _make_code_node

    nodes_to_replace = []

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue

        content = node.get('content', '')
        new_nodes      = []
        raw_code_buf   = []
        # Per-directive accumulators (reset each time a directive completes):
        raw_dir_buf    = []   # `!<`-stripped XML lines for current directive
        raw_dir_lines  = []   # original raw lines for current directive
        # Synthetic wrapping for each emitted directive:
        raw_opener     = None  # the `!![\n` line that opened the current XML block
        raw_closer     = None  # the corresponding `!!]\n` line (derived from opener)
        in_xml         = False
        in_directive   = False
        directive_root = None
        # Absolute (1-based) line tracking in the original source file.
        current_line   = node['line']   # line number of the line being read
        code_run_line  = node['line']   # first line of the buffered code run
        dir_run_line   = node['line']   # first line of the current directive

        def _flush_code_buf():
            if raw_code_buf:
                new_nodes.append(_make_code_node(
                    ''.join(raw_code_buf), node['source'], code_run_line))
                raw_code_buf.clear()

        def _emit_directive(pending_dir):
            """Wrap a parsed directive in synthetic `!![ ... !!]` and push it."""
            pending_dir['firstChild'] = {
                'type':       'code',
                'content':    (raw_opener or '') + ''.join(raw_dir_lines)
                              + (raw_closer or ''),
                'parent':     pending_dir,
                'sibling':    None,
                'firstChild': None,
                'source':     node['source'],
                # The serialized content starts with the `!![` marker, one
                # line before the directive tag itself, so the line-map
                # anchor sits on the marker; the directive node's own
                # `line` stays on the tag for {introspection:location}.
                'line':       dir_run_line - 1,
            }
            _flush_code_buf()
            new_nodes.append(pending_dir)

        for raw_line in content.splitlines(keepends=True):
            stripped = re.sub(r'^\s*!<\s*', '', raw_line)

            if re.match(r'^\s*!!\]', raw_line):
                # Block is closing.  Flush any directive that didn't have its
                # own end tag (rare — usually directives self-close before this).
                if raw_dir_buf:
                    xml_text    = ''.join(raw_dir_buf)
                    pending_dir = _parse_directive_xml(xml_text, node, dir_run_line)
                    if pending_dir:
                        _emit_directive(pending_dir)
                    else:
                        # Unparseable trailing fragment — emit as raw code with
                        # the synthetic wrapper restored.
                        _flush_code_buf()
                        new_nodes.append(_make_code_node(
                            (raw_opener or '') + ''.join(raw_dir_lines)
                            + raw_line,
                            node['source'], dir_run_line))
                    raw_dir_buf.clear()
                    raw_dir_lines.clear()
                raw_opener     = None
                raw_closer     = None
                in_directive   = False
                directive_root = None
                in_xml         = False
                current_line  += 1
                continue

            if re.match(r'^\s*!!\[', raw_line):
                # Close the current code run here: the marker line (and any
                # non-directive XML block it opens) is dropped from the
                # output, so code resuming after the block must start a new
                # node — otherwise its `.lmap` anchor would be off by the
                # number of dropped lines.
                _flush_code_buf()
                in_xml     = True
                raw_opener = raw_line
                # Derive `!!]` from `!![` by replacing `[` → `]`.
                raw_closer = raw_opener.replace('[', ']')
                current_line += 1
                continue

            if in_xml:
                m = re.match(r'^\s*<([^\s>/]+)', stripped)
                if m and not in_directive:
                    directive_root = m.group(1)
                    in_directive   = True
                    dir_run_line   = current_line
                if in_directive:
                    stripped = stripped.replace('&nbsp;', ' ')
                    raw_dir_buf.append(stripped)
                    raw_dir_lines.append(raw_line)
                    # Three end-tag forms: a closing `</tag>`, a self-closing
                    # `<tag attr=…/>` (attributes may contain `/` — e.g. URLs —
                    # so use `.*` greedily), or a bare `<tag/>`.
                    end1 = re.search(
                        r'</\s*' + re.escape(directive_root) + r'\s*>', stripped)
                    end2 = re.search(
                        r'\s*<' + re.escape(directive_root) + r'\s.*/>',
                        stripped)
                    end3 = re.search(
                        r'\s*<' + re.escape(directive_root) + r'\s*/>', stripped)
                    if end1 or end2 or end3:
                        xml_text    = ''.join(raw_dir_buf)
                        pending_dir = _parse_directive_xml(xml_text, node, dir_run_line)
                        if pending_dir:
                            _emit_directive(pending_dir)
                        else:
                            # Unparseable directive fragment — emit verbatim.
                            _flush_code_buf()
                            new_nodes.append(_make_code_node(
                                (raw_opener or '') + ''.join(raw_dir_lines)
                                + (raw_closer or ''),
                                node['source'], dir_run_line))
                        raw_dir_buf.clear()
                        raw_dir_lines.clear()
                        in_directive   = False
                        directive_root = None
                current_line += 1
                continue

            if not raw_code_buf:
                code_run_line = current_line
            raw_code_buf.append(raw_line)
            current_line += 1

        if raw_code_buf:
            new_nodes.append(_make_code_node(
                ''.join(raw_code_buf), node['source'], code_run_line))

        if len(new_nodes) != 1 or new_nodes[0] is not node:
            nodes_to_replace.append((node, new_nodes))

    for old_node, new_node_list in nodes_to_replace:
        replace_node(old_node, new_node_list)

def post_process_directives(tree):
    """Verify that every directive node has been processed.

    Called after the Process/* passes complete; raises RuntimeError on the
    first unprocessed directive encountered.

    Directives whose type is in the `NonProcessed` exemption list are
    forgiven even if they have no `processed` flag — code-generating hooks
    (DeepCopyActions, StateStorable, Enumeration, FunctionClass, …) inject
    fresh `<methods>` blocks late in the pipeline, after `nonProcessed`
    has already walked the tree, so those directives never get marked.
    """
    from Galacticus.Build.SourceTree import walk_tree
    from Galacticus.Build.SourceTree.Process.NonProcessed import (
        is_non_processed_type,
    )

    for node in walk_tree(tree):
        if 'directive' not in node:
            continue
        directive = node.setdefault('directive', {})
        if directive.get('processed'):
            continue
        if is_non_processed_type(node.get('type')):
            continue
        file_node = node
        while file_node is not None and file_node.get('type') != 'file':
            file_node = file_node.get('parent')
        file_name = file_node.get('name') if file_node else 'unknown'
        raise RuntimeError(
            f"directive '{node.get('type')}' was not processed at line "
            f"{node.get('line', 0)} in {file_name}")
