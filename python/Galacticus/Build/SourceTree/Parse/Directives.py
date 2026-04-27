# Contains a Python module which implements parsing of directives in the Galacticus preprocessor system.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/SourceTree/Parse/Directives.pm

import re
import os
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', ''), 'python'))

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

    Mirrors the Perl `our $stateStorables` caching pattern at Directives.pm:20-25.
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
# Mirrors the `fill_in_string` templates at Directives.pm:123-250 (functionClass)
# and Directives.pm:255-273 (eventHookStatic).
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


def _validate_directive(directive_name, xml_text, context_node):
    """Validate the directive XML against the appropriate XSD schema if possible.

    Mirrors the validation block at Directives.pm:119-288.  If lxml is not
    installed, or no schema applies, validation is skipped silently (matching
    the Perl code's behavior when $BUILDPATH is unset).
    """
    if not _HAS_LXML:
        return
    state_storables = _load_state_storables()
    schema_xml = None

    function_classes = (state_storables or {}).get('functionClasses') or {}
    event_hook_statics = (state_storables or {}).get('eventHookStatics') or []
    if isinstance(event_hook_statics, dict):
        event_hook_statics = [event_hook_statics]

    if state_storables and (directive_name + "Class") in function_classes:
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
            f"in file {file_name} at line {context_node.get('line', 0)}:\n"
            f"{schema.error_log}")


def _parse_directive_xml(xml_text, context_node):
    """Parse accumulated XML text into a directive node dict.

    Mirrors the `XML::Simple->XMLin(..., keepRoot => 1)` call at Directives.pm:111.
    Returns None if the text does not parse as XML.
    """
    try:
        elem = ET.fromstring(xml_text)
    except ET.ParseError:
        # Wrap in a root element and try again (handles bare fragments).
        try:
            elem = ET.fromstring('<root>' + xml_text + '</root>')
            if len(elem) == 1:
                elem = list(elem)[0]
        except ET.ParseError:
            return None

    directive_name = elem.tag
    _validate_directive(directive_name, xml_text, context_node)
    directive_dict = xml_to_dict(elem)
    return {
        'type':       directive_name,
        'directive':  directive_dict,
        'parent':     None,
        'firstChild': None,
        'sibling':    None,
        'source':     context_node.get('source', 'unknown'),
        'line':       context_node.get('line',   0),
    }


def parse_directives(tree):
    """Walk the tree replacing XML directive comment blocks with directive nodes.

    Mirrors Parse_Directives() from perl/Galacticus/Build/SourceTree/Parse/Directives.pm.

    Directives are delimited by:
      !![          (opening marker)
      !< <tagname ...>  (XML content lines — '!<' prefix stripped)
      !!]          (closing marker)

    A single `!![ ... !!]` block may contain multiple directives back-to-back;
    each is emitted as its own node (the Perl flushes on every closing tag).
    """
    from Galacticus.Build.SourceTree import walk_tree, replace_node, _make_code_node

    nodes_to_replace = []

    for node in walk_tree(tree):
        if node.get('type') != 'code':
            continue

        content = node.get('content', '')
        new_nodes      = []
        raw_code_buf   = []
        raw_dir_buf    = []   # XML text for the directive currently being assembled.
        raw_dir_lines  = []   # raw source lines for the directive currently being assembled.
        in_xml         = False
        in_directive   = False
        directive_root = None
        raw_opener     = None  # the actual `!![` line that opened the current block.

        def _emit(directive_node):
            """Attach a synthesized firstChild raw-text block, flush any
            pending plain-code buffer, and append the directive to new_nodes."""
            raw_closer = re.sub(r'\[', ']', raw_opener) if raw_opener else '!!]\n'
            directive_node['firstChild'] = {
                'type':       'code',
                'content':    (raw_opener or '') + ''.join(raw_dir_lines) + raw_closer,
                'parent':     directive_node,
                'sibling':    None,
                'firstChild': None,
                'source':     node['source'],
                'line':       node['line'],
            }
            if raw_code_buf:
                new_nodes.append(_make_code_node(
                    ''.join(raw_code_buf), node['source'], node['line']))
                raw_code_buf.clear()
            new_nodes.append(directive_node)

        for raw_line in content.splitlines(keepends=True):
            stripped = re.sub(r'^\s*!<\s*', '', raw_line)

            if re.match(r'^\s*!!\]', raw_line):
                # End of XML block.  If a directive is mid-assembly (malformed
                # source missing a closing tag), parse what we have.
                if raw_dir_buf:
                    pending_dir = _parse_directive_xml(''.join(raw_dir_buf), node)
                    raw_dir_buf = []
                    if pending_dir:
                        _emit(pending_dir)
                in_xml         = False
                in_directive   = False
                directive_root = None
                raw_opener     = None
                raw_dir_lines  = []
                continue

            if re.match(r'^\s*!!\[', raw_line):
                in_xml        = True
                raw_opener    = raw_line
                raw_dir_lines = []
                continue

            if in_xml:
                m = re.match(r'^\s*<([^\s>/]+)', stripped)
                if m and not in_directive:
                    directive_root = m.group(1)
                    in_directive   = True
                    raw_dir_lines  = []
                if in_directive:
                    raw_dir_lines.append(raw_line)
                    stripped = stripped.replace('&nbsp;', ' ')
                    raw_dir_buf.append(stripped)
                    end1 = re.search(
                        r'</\s*' + re.escape(directive_root) + r'\s*>', stripped)
                    end2 = re.match(
                        r'^\s*<' + re.escape(directive_root) + r'(\s[^/]*)?\s*/>',
                        stripped)
                    if end1 or end2:
                        pending_dir  = _parse_directive_xml(
                            ''.join(raw_dir_buf), node)
                        raw_dir_buf  = []
                        if pending_dir:
                            _emit(pending_dir)
                        in_directive   = False
                        directive_root = None
                        raw_dir_lines  = []
                continue

            raw_code_buf.append(raw_line)

        if raw_code_buf:
            new_nodes.append(_make_code_node(
                ''.join(raw_code_buf), node['source'], node['line']))

        if len(new_nodes) != 1 or new_nodes[0] is not node:
            nodes_to_replace.append((node, new_nodes))

    for old_node, new_node_list in nodes_to_replace:
        replace_node(old_node, new_node_list)


def post_process_directives(tree):
    """Verify that every directive node has been processed.

    Mirrors PostProcess_Directives() at Directives.pm:347-361.  Called after
    the Process/* passes complete; raises RuntimeError on the first unprocessed
    directive encountered.
    """
    from Galacticus.Build.SourceTree import walk_tree

    for node in walk_tree(tree):
        if 'directive' not in node:
            continue
        directive = node.get('directive') or {}
        if not directive.get('processed'):
            file_node = node
            while file_node is not None and file_node.get('type') != 'file':
                file_node = file_node.get('parent')
            file_name = file_node.get('name') if file_node else 'unknown'
            raise RuntimeError(
                f"directive '{node.get('type')}' was not processed at line "
                f"{node.get('line', 0)} in {file_name}")
