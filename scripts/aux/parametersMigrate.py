#!/usr/bin/env python3
"""Update a Galacticus parameter file from its last modified revision to the current revision."""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from lxml import etree


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Migrate a Galacticus parameter file to the current revision."
    )
    parser.add_argument("inputFile", help="Input parameter file")
    parser.add_argument("outputFile", help="Output parameter file")
    parser.add_argument("--validate", default="yes", help="Validate the parameter file (default: yes)")
    parser.add_argument("--prettyify", default="no", help="Pretty-print the output (default: no)")
    parser.add_argument("--ignoreWhiteSpaceChanges", default="yes", help="Ignore whitespace changes (default: yes)")
    parser.add_argument("--outputFormatVersion", type=int, default=2, help="Output format version (default: 2)")
    parser.add_argument("--lastModifiedRevision", default=None, help="Override last modified revision hash")
    parser.add_argument("--timeStamp", default=None, help="Override timestamp")
    return parser.parse_args()


def preprocess_multiline_attributes(input_filename):
    """Pre-process the file to concatenate any attribute values that are split across multiple lines.

    Replaces newlines inside unclosed quotes with %%NEWLINE%%, skipping XML comments.
    Returns the path to the temporary pre-processed file.
    """
    fout = tempfile.NamedTemporaryFile(mode='w', delete=False)
    tmp_filename = fout.name
    in_multiline = False
    in_comment = False
    with open(input_filename, "r") as fin:
        for line in fin:
            count_quotes = line.count('"')
            if "<!--" in line:
                in_comment = True
            if count_quotes % 2 == 1 and not in_comment:
                in_multiline = not in_multiline
            if in_multiline and not in_comment:
                line = line.replace("\n", "%%NEWLINE%%")
            if "-->" in line:
                in_comment = False
            fout.write(line)
    return tmp_filename


def restore_multiline_attributes(input_filename, output_filename):
    """Restore %%NEWLINE%% placeholders back to actual newlines."""
    with open(input_filename, "r") as fin, open(output_filename, "w") as fout:
        for line in fin:
            fout.write(line.replace("%%NEWLINE%%", "\n"))


def parse_migrations(path):
    """Parse migrations.xml into a structured dict.

    Returns a dict with:
        'migration': list of migration dicts, each with 'commit' and 'translation' (list)
        'default': list of default dicts, each with 'parameter' and 'value'
    """
    tree = etree.parse(path)
    root = tree.getroot()
    migrations = {"migration": [], "default": []}
    for migration_elem in root.findall("migration"):
        migration = {
            "commit": migration_elem.get("commit"),
            "translation": [],
        }
        for trans_elem in migration_elem.findall("translation"):
            translation = {}
            if trans_elem.get("xpath") is not None:
                translation["xpath"] = trans_elem.get("xpath")
            if trans_elem.get("function") is not None:
                translation["function"] = trans_elem.get("function")
            name_elem = trans_elem.find("name")
            if name_elem is not None:
                translation["name"] = {"new": name_elem.get("new")}
            value_elem = trans_elem.find("value")
            if value_elem is not None:
                translation["value"] = {
                    "old": value_elem.get("old"),
                    "new": value_elem.get("new"),
                }
            remove_elem = trans_elem.find("remove")
            if remove_elem is not None:
                translation["remove"] = True
            migration["translation"].append(translation)
        migrations["migration"].append(migration)
    for default_elem in root.findall("default"):
        migrations["default"].append(
            {
                "parameter": default_elem.get("parameter"),
                "value": default_elem.get("value"),
            }
        )
    return migrations


def as_array(item):
    """Normalize a value to a list (like Perl's List::ExtraUtils::as_array)."""
    if item is None:
        return []
    if isinstance(item, list):
        return item
    return [item]


def insert_after(parent, new_node, ref_node):
    """Insert new_node as a child of parent immediately after ref_node."""
    children = list(parent)
    if ref_node is None:
        parent.append(new_node)
    else:
        idx = children.index(ref_node)
        parent.insert(idx + 1, new_node)


def insert_before(parent, new_node, ref_node):
    """Insert new_node as a child of parent immediately before ref_node."""
    if ref_node is None:
        parent.append(new_node)
    else:
        children = list(parent)
        idx = children.index(ref_node)
        parent.insert(idx, new_node)


def git_is_tracked(filename):
    """Check if a file is tracked by git."""
    result = subprocess.run(
        ["git", "ls-files", "--error-unmatch", filename],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def git_head_hash():
    """Get the current HEAD commit hash."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        message = (
            "Failed to determine the current git HEAD revision. "
            "Ensure that git is installed and that this script is run inside a git repository."
        )
        if e.stderr:
            message += f" Git error: {e.stderr.strip()}"
        raise RuntimeError(message) from e
    return result.stdout.strip()


def git_last_modified_hash(filename):
    """Get the hash of the last commit that modified the given file."""
    try:
        result = subprocess.run(
            ["git", "log", "-n", "1", "--pretty=format:%H", "--", filename],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        message = (
            f"Failed to determine the last git commit that modified '{filename}'. "
            "Ensure that git is installed, this script is run inside a git repository, "
            "and that the file is tracked."
        )
        if e.stderr:
            message += f" Git error: {e.stderr.strip()}"
        raise RuntimeError(message) from e
    return result.stdout.strip()


def git_ancestry(hash_from, hash_to):
    """Get the ancestry path between two commits (oldest first)."""
    try:
        result = subprocess.run(
            ["git", "rev-list", "--ancestry-path", f"{hash_from}..{hash_to}"],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        message = (
            f"Failed to determine git ancestry from '{hash_from}' to '{hash_to}'. "
            "Ensure that git is installed, this script is run inside a git repository, "
            "and that both revisions are valid."
        )
        if e.stderr:
            message += f" Git error: {e.stderr.strip()}"
        raise RuntimeError(message) from e
    hashes = [h for h in result.stdout.strip().split("\n") if h]
    hashes.reverse()  # Oldest first, matching Perl's reverse(@ancestry)
    return hashes


# ---------------------------------------------------------------------------
# Migrate function (Perl lines 136-369) -- placeholder for Step 3
# ---------------------------------------------------------------------------

def migrate(input_doc, parameters, root_level, is_grid, input_filename, options, hash_head, is_in_git, migrations):
    """Apply all migrations to a <parameters> element."""

    output_format_version = str(options.outputFormatVersion)

    # Check for a format version element.
    format_version_elems = parameters.findall("formatVersion")
    format_version = format_version_elems[0] if format_version_elems else None
    if format_version is not None:
        format_version.text = output_format_version
    elif root_level:
        format_version_node = etree.Element("formatVersion")
        format_version_node.text = output_format_version
        format_version_node.tail = "\n  "
        parameters.insert(0, format_version_node)

    # Validate the parameter file.
    if options.validate == "yes":
        exec_path = os.environ.get("GALACTICUS_EXEC_PATH", ".")
        result = subprocess.run(
            [os.path.join(exec_path, "scripts", "aux", "validateParameters.py"), input_filename]
        )
        if result.returncode != 0:
            sys.exit(f'input file "{input_filename}" is not a valid Galacticus parameter file')

    # Find last modified hash.
    hash_last_modified = None
    element_last_modified_list = parameters.findall("lastModified")
    element_last_modified = element_last_modified_list[0] if element_last_modified_list else None
    if element_last_modified is None:
        element_last_modified = etree.Element("lastModified")
        element_last_modified.tail = "\n  "
        parameters.insert(0, element_last_modified)

    if options.lastModifiedRevision is not None:
        hash_last_modified = options.lastModifiedRevision
    else:
        # Look for a last modification hash.
        if element_last_modified is not None and element_last_modified.get("revision") is not None:
            hash_last_modified = element_last_modified.get("revision")
        elif is_in_git:
            hash_last_modified = git_last_modified_hash(input_filename)
        else:
            hash_last_modified = "6eab8997cd73cb0a474228ade542d133890ad138^"

    # Update the last modified metadata.
    time_stamp = options.timeStamp if options.timeStamp is not None else datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    element_last_modified.set("revision", hash_head)
    element_last_modified.set("time", time_stamp)

    # Find the ancestry.
    ancestry = git_ancestry(hash_last_modified, hash_head)

    # Iterate over the revision ancestry.
    for hash_ancestor in ancestry:
        matched_migrations = [m for m in as_array(migrations["migration"]) if m["commit"] == hash_ancestor]
        if len(matched_migrations) != 1:
            continue
        migration = matched_migrations[0]
        # Report.
        print(f"Updating to revision {hash_ancestor}")

        # Iterate over translations.
        for translation in as_array(migration["translation"]):

            # Handle special cases (function dispatch).
            if "function" in translation:
                func_name = translation["function"]
                if func_name in SPECIAL_FUNCTIONS:
                    SPECIAL_FUNCTIONS[func_name](input_doc, parameters, is_grid)
                else:
                    sys.exit(f"Unknown special migration function: {func_name}")

            # Handle removals.
            if "remove" in translation:
                for parameter in parameters.xpath(translation["xpath"]):
                    print(f"   remove parameter: {parameter.tag}")
                    parameter.getparent().remove(parameter)

            # Translate names.
            if "name" in translation:
                for parameter in parameters.xpath(translation["xpath"]):
                    print(f"   translate parameter name: {parameter.tag} --> {translation['name']['new']}")
                    leaf_name = translation["name"]["new"]
                    parameter.tag = leaf_name

            # Translate values.
            if "value" in translation:
                for parameter in parameters.xpath(translation["xpath"]):
                    # Find all value sub-elements, or use the parameter itself.
                    value_children = parameter.findall("value")
                    if value_children:
                        all_values = value_children
                    else:
                        all_values = [parameter]
                    for value_elem in all_values:
                        # Get the current value text.
                        if value_elem is parameter:
                            values_text = value_elem.get("value", "")
                        else:
                            values_text = value_elem.text or ""
                        values_text = values_text.strip()
                        values = values_text.split() if values_text else []
                        for i, this_value in enumerate(values):
                            if this_value == translation["value"]["old"]:
                                print(f"   translate parameter value: {translation['xpath']}")
                                print(f"                                 {this_value} --> {translation['value']['new']}")
                                values[i] = translation["value"]["new"]
                        # Update the values in the parameter.
                        if value_elem is parameter:
                            value_elem.set("value", " ".join(values))
                        else:
                            value_elem.text = " ".join(values)

    # Put subparameters into their host parameter.
    # We need to iterate over a snapshot since we're modifying the tree.
    for parameter in list(parameters):
        if not isinstance(parameter.tag, str):
            continue  # Skip comments/PIs
        m = re.match(r"(.+)--(.+)", parameter.tag)
        if m:
            host_name = m.group(1)
            sub_name = m.group(2)
            host_name_trimmed = re.sub(r"\..+\.", "", host_name)
            sibling = parameter.getnext()
            parameter.tag = sub_name
            host_found = False
            for host_parameter in list(parameters):
                if host_parameter.tag == host_name_trimmed:
                    # Move the sub-parameter into the existing host.
                    parameters.remove(parameter)
                    if len(host_parameter) == 0:
                        host_parameter.text = "\n    "
                    parameter.tail = "\n  "
                    host_parameter.append(parameter)
                    host_found = True
                    break
            if not host_found:
                # Create the new host node.
                host_leaf_name = re.sub(r"\..+\.", "", host_name)
                if re.search(r"\..+\.", host_name):
                    default_value = re.sub(r".*\.(.+)\.", r"\1", host_name)
                else:
                    defaults = [d for d in migrations["default"] if d["parameter"] == host_leaf_name]
                    if len(defaults) != 1:
                        sys.exit(
                            f'parametersMigrate.py: attempting to insert a "{host_leaf_name}" element, '
                            "but no default value is known"
                        )
                    default_value = defaults[0]["value"]
                parameter_node = etree.Element(host_leaf_name)
                parameter_node.set("value", default_value)
                parameter_node.text = "\n    "
                parameters.remove(parameter)
                parameter.tail = "\n  "
                parameter_node.append(parameter)
                if sibling is not None:
                    insert_before(parameters, parameter_node, sibling)
                else:
                    parameters.append(parameter_node)

    # Search for duplicated parameters.
    duplicate_wrappers = {"mergerTreeOperator": "sequence"}
    processed_tags = set()
    for parameter in list(parameters):
        if not isinstance(parameter.tag, str):
            continue  # Skip comments/PIs
        node_name = parameter.tag
        if node_name in processed_tags:
            continue
        duplicates = [p for p in parameters if p.tag == node_name]
        if len(duplicates) > 1 and node_name in duplicate_wrappers:
            processed_tags.add(node_name)
            sibling = duplicates[-1].getnext()
            wrapper_node = etree.Element(node_name)
            wrapper_node.set("value", duplicate_wrappers[node_name])
            wrapper_node.text = "\n    "
            for duplicate in duplicates:
                parameters.remove(duplicate)
                duplicate.tail = "\n  "
                wrapper_node.append(duplicate)
            if sibling is not None:
                insert_before(parameters, wrapper_node, sibling)
            else:
                parameters.append(wrapper_node)


# ---------------------------------------------------------------------------
# Special migration functions (Perl lines 371-1110)
# ---------------------------------------------------------------------------


def radiation_field_intergalactic_background_cmb(input_doc, parameters, is_grid):
    """Special handling to add CMB radiation into the intergalactic background radiation."""
    for node in parameters.xpath(".//radiationFieldIntergalacticBackground[@value]"):
        print("   translate special './/radiationFieldIntergalacticBackground[@value]'")
        # Construct new summation and CMB nodes, and a copy of the original node.
        summation_node = etree.Element("radiationFieldIntergalacticBackground")
        cmb_node = etree.Element("radiationField")
        original_node = etree.Element("radiationField")
        summation_node.set("value", "summation")
        cmb_node.set("value", "cosmicMicrowaveBackground")
        original_node.set("value", node.get("value"))
        # Move any children of the original node to the copy.
        for child in list(node):
            node.remove(child)
            original_node.append(child)
        # Also move text content.
        if node.text:
            original_node.text = node.text
            node.text = None
        # Assemble our new nodes.
        summation_node.append(cmb_node)
        summation_node.append(original_node)
        # Insert the new parameters and remove the original.
        parent = node.getparent()
        insert_after(parent, summation_node, node)
        parent.remove(node)


def black_hole_seed_mass(input_doc, parameters, is_grid):
    """Special handling to add black hole seed mass models."""
    do_translate = False
    mass_seed = 100.0  # Default value.
    component_black_hole = None
    for node in parameters.xpath(
        ".//componentBlackHole[@value='simple' or @value='standard' or @value='nonCentral']/massSeed[@value]"
    ):
        print("   translate special './/componentBlackHole[@value]/massSeed[@value]'")
        do_translate = True
        component_black_hole = node.getparent()
        mass_seed = node.get("value")
        node.getparent().remove(node)
    if not do_translate:
        return
    # Find nodeOperators.
    node_operators = parameters.xpath(".//nodeOperator[@value='multi']")
    if len(node_operators) == 0:
        sys.exit("can not find any `nodeOperator[@value='multi']` into which to insert a black hole seed operator")
    if len(node_operators) > 1:
        sys.exit("found multiple `nodeOperator[@value='multi']` nodes - unknown into which to insert a black hole seed operator")
    # Create a new node operator and insert into the list.
    operator_node = etree.Element("nodeOperator")
    seed_node = etree.Element("blackHoleSeeds")
    mass_node = etree.Element("mass")
    spin_node = etree.Element("spin")
    operator_node.set("value", "blackHolesSeed")
    seed_node.set("value", "fixed")
    mass_node.set("value", str(mass_seed))
    spin_node.set("value", "0.0")
    if is_grid:
        operator_node.set("iterable", "no")
        seed_node.set("iterable", "no")
    # Assemble our new nodes.
    seed_node.append(mass_node)
    seed_node.append(spin_node)
    # Insert the new parameters.
    node_operators[0].append(operator_node)
    if component_black_hole is not None:
        insert_after(component_black_hole.getparent(), seed_node, component_black_hole)
    else:
        parameters.append(seed_node)


def black_hole_physics(input_doc, parameters, is_grid):
    """Special handling to move black hole physics from components to operators."""
    defaults = {
        "simple": {
            "heatsHotHalo": "false",
            "efficiencyHeating": "1.0e-3",
            "efficiencyWind": "2.2157e-3",
            "growthRatioToStellarSpheroid": "1.0e-3",
        },
        "standard": {
            "bondiHoyleAccretionEnhancementSpheroid": "5.0e0",
            "bondiHoyleAccretionEnhancementHotHalo": "6.0e0",
            "bondiHoyleAccretionHotModeOnly": "true",
            "bondiHoyleAccretionTemperatureSpheroid": "1.0e2",
            "efficiencyWind": "2.4e-3",
            "efficiencyWindScalesWithEfficiencyRadiative": "false",
            "efficiencyRadioMode": "1.0",
            "heatsHotHalo": "true",
        },
    }
    component_properties = {}
    component_node = None
    do_translate = False
    for node in parameters.xpath(
        ".//componentBlackHole[@value='simple' or @value='standard' or @value='nonCentral']"
    ):
        print("   translate special './/componentBlackHole[@value]'")
        do_translate = True
        component_node = node
        component_type = node.get("value")
        if component_type == "simple":
            component_properties["type"] = "simple"
        elif component_type in ("standard", "nonCentral"):
            component_properties["type"] = "standard"
        # Extract all sub-parameters.
        for node_child in node.xpath("*[@value]"):
            component_properties[node_child.tag] = node_child.get("value")
            node.remove(node_child)
        # Insert default parameters.
        if component_properties["type"] in defaults:
            for param_name, param_value in defaults[component_properties["type"]].items():
                if param_name not in component_properties:
                    component_properties[param_name] = param_value
    if not do_translate:
        return
    # Find nodeOperators.
    node_operators = parameters.xpath(".//nodeOperator[@value='multi']")
    if len(node_operators) == 0:
        sys.exit("can not find any `nodeOperator[@value='multi']` into which to insert a black hole seed operator")
    if len(node_operators) > 1:
        sys.exit("found multiple `nodeOperator[@value='multi']` nodes - unknown into which to insert a black hole seed operator")
    # Build node operators.
    operator_accretion_node = etree.Element("nodeOperator")
    operator_winds_node = etree.Element("nodeOperator")
    operator_cgm_heat_node = etree.Element("nodeOperator")
    operator_accretion_node.set("value", "blackHolesAccretion")
    operator_winds_node.set("value", "blackHolesWinds")
    operator_cgm_heat_node.set("value", "blackHolesCGMHeating")
    if is_grid:
        operator_accretion_node.set("iterable", "no")
        operator_winds_node.set("iterable", "no")
        operator_cgm_heat_node.set("iterable", "no")
    node_operators[0].append(operator_accretion_node)
    node_operators[0].append(operator_winds_node)
    if component_properties.get("heatsHotHalo") == "true":
        node_operators[0].append(operator_cgm_heat_node)
    # Handle the "simple" black hole component.
    if component_properties["type"] == "simple":
        # Accretion.
        accretion_node = etree.Element("blackHoleAccretionRate")
        growth_rate_node = etree.Element("growthRatioToStellarSpheroid")
        accretion_node.set("value", "spheroidTracking")
        growth_rate_node.set("value", component_properties["growthRatioToStellarSpheroid"])
        accretion_node.append(growth_rate_node)
        insert_after(component_node.getparent(), accretion_node, component_node)
        # Winds.
        wind_node = etree.Element("blackHoleWind")
        efficiency_node = etree.Element("efficiencyWind")
        wind_node.set("value", "simple")
        efficiency_node.set("value", component_properties["efficiencyWind"])
        wind_node.append(efficiency_node)
        insert_after(component_node.getparent(), wind_node, component_node)
        if component_properties.get("heatsHotHalo") == "true":
            # CGM heating.
            heating_node = etree.Element("blackHoleCGMHeating")
            eff_node = etree.Element("efficiencyHeating")
            heating_node.set("value", "quasistatic")
            eff_node.set("value", component_properties["efficiencyHeating"])
            heating_node.append(eff_node)
            insert_after(component_node.getparent(), heating_node, component_node)
    # Handle the "standard" black hole component.
    if component_properties["type"] == "standard":
        # Accretion.
        accretion_node = etree.Element("blackHoleAccretionRate")
        spheroid_node = etree.Element("bondiHoyleAccretionEnhancementSpheroid")
        hot_halo_node = etree.Element("bondiHoyleAccretionEnhancementHotHalo")
        hot_mode_node = etree.Element("bondiHoyleAccretionHotModeOnly")
        temperature_node = etree.Element("bondiHoyleAccretionTemperatureSpheroid")
        accretion_node.set("value", "standard")
        spheroid_node.set("value", component_properties["bondiHoyleAccretionEnhancementSpheroid"])
        hot_halo_node.set("value", component_properties["bondiHoyleAccretionEnhancementHotHalo"])
        hot_mode_node.set("value", component_properties["bondiHoyleAccretionHotModeOnly"])
        temperature_node.set("value", component_properties["bondiHoyleAccretionTemperatureSpheroid"])
        accretion_node.append(spheroid_node)
        accretion_node.append(hot_halo_node)
        accretion_node.append(hot_mode_node)
        accretion_node.append(temperature_node)
        insert_after(component_node.getparent(), accretion_node, component_node)
        # Winds.
        wind_node = etree.Element("blackHoleWind")
        efficiency_node = etree.Element("efficiencyWind")
        scale_node = etree.Element("efficiencyWindScalesWithEfficiencyRadiative")
        wind_node.set("value", "ciotti2009")
        efficiency_node.set("value", component_properties["efficiencyWind"])
        scale_node.set("value", component_properties["efficiencyWindScalesWithEfficiencyRadiative"])
        wind_node.append(efficiency_node)
        wind_node.append(scale_node)
        insert_after(component_node.getparent(), wind_node, component_node)
        if component_properties.get("heatsHotHalo") == "true":
            # CGM heating.
            heating_node = etree.Element("blackHoleCGMHeating")
            eff_node = etree.Element("efficiencyRadioMode")
            heating_node.set("value", "jetPower")
            eff_node.set("value", component_properties["efficiencyRadioMode"])
            heating_node.append(eff_node)
            insert_after(component_node.getparent(), heating_node, component_node)


def model_parameter_xpath(input_doc, parameters, is_grid):
    """Special handling to switch `modelParameter` names to use XPath syntax."""
    for model_parameter in parameters.xpath(".//modelParameter[@value='active' or @value='inactive']"):
        print("   translate special './/modelParameter[@value]'")
        for name_node in model_parameter.findall("name"):
            value = name_node.get("value")
            value = value.replace("::", "/")
            # Increment indices to XPath standard 1-indexing.
            value = re.sub(
                r"([\[\{])(\d+)([\]\}])",
                lambda m: m.group(1) + str(int(m.group(2)) + 1) + m.group(3),
                value,
            )
            name_node.set("value", value)


def method_suffix_remove(input_doc, parameters, is_grid):
    """Special handling to remove the 'Method' suffix from parameter names."""
    print("   translate special - remove Method suffixes")
    # Names where the Method suffix persists.
    exceptions = {
        "readSubhaloAngularMomentaMethod",
        "subhaloAngularMomentaMethod",
        "diskRadiusSolverCole2000Method",
        "duttonMaccio2014DensityContrastMethod",
        "duttonMaccio2014DensityProfileMethod",
    }
    for parameter in parameters.xpath(".//*"):
        node_name = parameter.tag
        if node_name in exceptions:
            continue
        m = re.match(r"^treeNodeMethod(.*)$", node_name)
        if m:
            parameter.tag = "component" + m.group(1)
        else:
            m = re.match(r"^(.*)Method$", node_name)
            if m:
                parameter.tag = m.group(1)


def satellite_orphanize(input_doc, parameters, is_grid):
    """Special handling to add a nodeOperator to orphanize satellites."""
    nodes = parameters.xpath(".//componentSatellite[@value='preset']")
    if len(nodes) <= 0:
        return
    print("   translate special './/componentSatellite[@value='preset']'")
    # Find nodeOperators.
    node_operators = parameters.xpath(".//nodeOperator[@value='multi']")
    if len(node_operators) == 0:
        sys.exit("can not find any `nodeOperator[@value='multi']` into which to insert a satellite orphanizer operator")
    if len(node_operators) > 1:
        sys.exit("found multiple `nodeOperator[@value='multi']` nodes - unknown into which to insert a satellite orphanizer operator")
    # Build node operator.
    operator_orphanize = etree.Element("nodeOperator")
    operator_orphanize.set("value", "satelliteOrphanize")
    if is_grid:
        operator_orphanize.set("iterable", "no")
    node_operators[0].append(operator_orphanize)


def black_hole_non_central(input_doc, parameters, is_grid):
    """Special handling to add nodeOperators for non-central black hole evolution."""
    nodes = parameters.xpath(".//componentBlackHole[@value='nonCentral']")
    if len(nodes) <= 0:
        return
    if len(nodes) > 1:
        sys.exit("found multiple `.//componentBlackHole[@value='nonCentral']` nodes - unknown what should be done in this situation")
    print("   translate special './/componentBlackHole[@value='nonCentral']'")
    # Look for any "tripleInteraction" option.
    triple_interactions = nodes[0].findall("tripleInteraction[@value]")
    if len(triple_interactions) > 1:
        sys.exit("found multiple `.//tripleInteraction[@value]` nodes - unknown what should be done in this situation")
    triple_interaction = (
        triple_interactions[0].get("value") == "true" if len(triple_interactions) == 1 else True
    )
    for ti in triple_interactions:
        ti.getparent().remove(ti)
    # Find nodeOperators.
    node_operators = parameters.xpath(".//nodeOperator[@value='multi']")
    if len(node_operators) == 0:
        sys.exit("can not find any `nodeOperator[@value='multi']` into which to insert a black hole operator")
    if len(node_operators) > 1:
        sys.exit("found multiple `nodeOperator[@value='multi']` nodes - unknown into which to insert a black hole operator")
    # Build node operators.
    operator_migration = etree.Element("nodeOperator")
    operator_triple = etree.Element("nodeOperator")
    if is_grid:
        operator_migration.set("iterable", "no")
        operator_triple.set("iterable", "no")
    operator_migration.set("value", "blackHolesRadialMigration")
    node_operators[0].append(operator_migration)
    if triple_interaction:
        operator_triple.set("value", "blackHolesTripleInteraction")
        node_operators[0].append(operator_triple)


def hot_halo_very_simple(input_doc, parameters, is_grid):
    """Special handling to add nodeOperators for 'very simple' hot halo evolution."""
    nodes = parameters.xpath(".//componentHotHalo[@value='verySimple' or @value='verySimpleDelayed']")
    if len(nodes) <= 0:
        return
    if len(nodes) > 1:
        sys.exit(
            "found multiple `.//componentHotHalo[@value='verySimple' or @value='verySimpleDelayed']` nodes"
            " - unknown what should be done in this situation"
        )
    print("   translate special './/componentHotHalo[@value='verySimple' or @value='verySimpleDelayed']'")
    # Find nodeOperators.
    node_operators = parameters.xpath(".//nodeOperator[@value='multi']")
    if len(node_operators) == 0:
        sys.exit("can not find any `nodeOperator[@value='multi']` into which to insert CGM operators")
    if len(node_operators) > 1:
        sys.exit("found multiple `nodeOperator[@value='multi']` nodes - unknown into which to insert CGM operators")
    # Build node operators.
    operator_outer_radius = etree.Element("nodeOperator")
    operator_cooling = etree.Element("nodeOperator")
    operator_accretion = etree.Element("nodeOperator")
    operator_starvation = etree.Element("nodeOperator")
    operator_outflow_reinc = etree.Element("nodeOperator")
    component = etree.Element("component")
    operator_outer_radius.set("value", "CGMOuterRadiusVirialRadius")
    operator_cooling.set("value", "CGMCoolingInflow")
    operator_accretion.set("value", "CGMAccretion")
    operator_starvation.set("value", "CGMStarvation")
    operator_outflow_reinc.set("value", "CGMOutflowReincorporation")
    component.set("value", "disk")
    if is_grid:
        operator_outer_radius.set("iterable", "no")
        operator_cooling.set("iterable", "no")
        operator_accretion.set("iterable", "no")
        operator_starvation.set("iterable", "no")
        operator_outflow_reinc.set("iterable", "no")
    operator_cooling.append(component)
    node_operators[0].append(operator_outer_radius)
    node_operators[0].append(operator_cooling)
    node_operators[0].append(operator_accretion)
    node_operators[0].append(operator_starvation)
    if nodes[0].get("value") == "verySimpleDelayed":
        node_operators[0].append(operator_outflow_reinc)


def collaborative_mpi(input_doc, parameters, is_grid):
    """Special handling to migrate the `collaborativeMPI` parameter."""
    nodes = parameters.xpath(".//posteriorSampleLikelihood[@value='galaxyPopulation']/collaborativeMPI")
    if len(nodes) <= 0:
        return
    print("   translate special './/posteriorSampleLikelihood[@value='galaxyPopulation']/collaborativeMPI'")
    # Replace each node.
    for node in nodes:
        count_groups = 1 if node.get("value") == "true" else -1
        count_collaborative_groups = etree.Element("countCollaborativeGroups")
        first_come_first_served = etree.Element("firstComeFirstServed")
        count_collaborative_groups.set("value", str(count_groups))
        first_come_first_served.set("value", "false")
        parent = node.getparent()
        # Insert the new elements after the old node, then remove the old node.
        insert_after(parent, count_collaborative_groups, node)
        insert_after(parent, first_come_first_served, count_collaborative_groups)
        parent.remove(node)


def hot_halo_standard_accretion(input_doc, parameters, is_grid):
    """Special handling to add and modify nodeOperators for CGM accretion in the standard hot halo component."""
    nodes = parameters.xpath(
        ".//componentHotHalo[@value='verySimple' or @value='verySimpleDelayed' "
        "or @value='standard' or @value='coldMode' or @value='outflowTracking']"
    )
    if len(nodes) <= 0:
        # None found - check if we have any componentHotHalo.
        nodes_any = parameters.xpath(".//componentHotHalo[@value]")
        if len(nodes_any) > 0:
            return
        # Check if we have nodeOperators present.
        node_ops = parameters.xpath(".//nodeOperator[@value='multi']")
        if len(node_ops) == 0:
            return
        # No componentHotHalo - insert one with the default value.
        component_hot_halo = etree.SubElement(parameters, "componentHotHalo")
        component_hot_halo.set("value", "standard")
        nodes = parameters.xpath(
            ".//componentHotHalo[@value='verySimple' or @value='verySimpleDelayed' "
            "or @value='standard' or @value='coldMode' or @value='outflowTracking']"
        )
    if len(nodes) > 1:
        sys.exit(
            "found multiple `.//componentHotHalo[@value='verySimple' or ...]` nodes"
            " - unknown what should be done in this situation"
        )
    print(
        "   translate special './/componentHotHalo[@value='verySimple' or @value='verySimpleDelayed'"
        " or @value='standard' or @value='coldMode' or @value='outflowTracking']'"
    )
    # Determine the type of hot halo component.
    hh_type = nodes[0].get("value")
    # Find nodeOperators.
    node_operators = parameters.xpath(".//nodeOperator[@value='multi']")
    if len(node_operators) == 0:
        sys.exit("can not find any `nodeOperator[@value='multi']` into which to insert CGM operators")
    if len(node_operators) > 1:
        sys.exit("found multiple `nodeOperator[@value='multi']` nodes - unknown into which to insert CGM operators")
    # Handle the verySimple or verySimpleDelayed cases.
    if hh_type in ("verySimple", "verySimpleDelayed"):
        # Find existing nodeOperatorCGMOutflowReincorporation.
        ops_reinc = node_operators[0].xpath(".//nodeOperator[@value='CGMOutflowReincorporation']")
        if len(ops_reinc) > 1:
            sys.exit("found multiple `nodeOperator[@value='CGMOutflowReincorporation']` nodes - unknown which to use")
        if len(ops_reinc) == 1:
            include_satellites = etree.Element("includeSatellites")
            include_satellites.set("value", "false")
            ops_reinc[0].append(include_satellites)
        # Find existing nodeOperatorCGMAccretion operator.
        ops_accretion = node_operators[0].xpath(".//nodeOperator[@value='CGMAccretion']")
        if len(ops_accretion) == 0:
            sys.exit("can not find any `nodeOperator[@value='CGMAccretion']` operators")
        if len(ops_accretion) > 1:
            sys.exit("found multiple `nodeOperator[@value='CGMAccretion']` nodes - unknown which to use")
        allow_negative = etree.Element("allowNegativeCGMMass")
        allow_negative.set("value", "false")
        ops_accretion[0].append(allow_negative)
        # Find existing nodeOperatorStarvation operator.
        ops_starvation = node_operators[0].xpath(".//nodeOperator[@value='CGMStarvation']")
        if len(ops_starvation) == 0:
            sys.exit("can not find any `nodeOperator[@value='CGMStarvation']` operators")
        if len(ops_starvation) > 1:
            sys.exit("found multiple `nodeOperator[@value='CGMStarvation']` nodes - unknown which to use")
        starve_outflows_only = etree.Element("starveOutflowsOnly")
        fraction_baryon_limit = etree.Element("fractionBaryonLimitInNodeMerger")
        starve_outflows_only.set("value", "false")
        fraction_baryon_limit.set("value", "false")
        ops_starvation[0].append(starve_outflows_only)
        ops_starvation[0].append(fraction_baryon_limit)
    # Handle the standard, coldMode, and outflowTracking cases.
    if hh_type in ("standard", "coldMode", "outflowTracking"):
        # Find pre-existing options.
        def _get_prior_value(xpath_expr, default):
            prior = nodes[0].xpath(xpath_expr)
            return prior[0].get("value") if len(prior) == 1 else default

        fraction_baryon_value = _get_prior_value(".//fractionBaryonLimitInNodeMerger[@value]", "false")
        starve_satellites_value = _get_prior_value(".//starveSatellites[@value]", "false")
        starve_satellites_outflowed_value = _get_prior_value(".//starveSatellitesOutflowed[@value]", "false")
        angular_momentum_value = _get_prior_value(".//angularMomentumAlwaysGrows[@value]", "false")
        outflow_to_cold_mode_value = _get_prior_value(".//outflowToColdMode[@value]", "false")
        # Remove obsoleted parameters.
        if hh_type != "coldMode":
            for xpath_expr in (".//fractionBaryonLimitInNodeMerger[@value]", ".//angularMomentumAlwaysGrows[@value]"):
                for elem in nodes[0].xpath(xpath_expr):
                    elem.getparent().remove(elem)
        # Remove starveSatellites and starveSatellitesOutflowed regardless.
        for xpath_expr in (".//starveSatellites[@value]", ".//starveSatellitesOutflowed[@value]"):
            for elem in nodes[0].xpath(xpath_expr):
                elem.getparent().remove(elem)
        # Add the accretion operator.
        operator_accretion = etree.Element("nodeOperator")
        operator_accretion.set("value", "CGMAccretion")
        if is_grid:
            operator_accretion.set("iterable", "no")
        node_operators[0].append(operator_accretion)
        allow_negative = etree.SubElement(operator_accretion, "allowNegativeCGMMass")
        allow_negative.set("value", "true")
        angular_momentum = etree.SubElement(operator_accretion, "angularMomentumAlwaysGrows")
        angular_momentum.set("value", angular_momentum_value)
        # Add a starvation operator if needed.
        if starve_satellites_value == "true" or starve_satellites_outflowed_value == "true":
            operator_starvation = etree.Element("nodeOperator")
            operator_starvation.set("value", "CGMStarvation")
            if is_grid:
                operator_starvation.set("iterable", "no")
            node_operators[0].append(operator_starvation)
            starve_outflows_only_value = (
                "true" if starve_satellites_outflowed_value == "true" and starve_satellites_value == "false" else "false"
            )
            starve_outflows_only = etree.SubElement(operator_starvation, "starveOutflowsOnly")
            starve_outflows_only.set("value", starve_outflows_only_value)
            fraction_baryon_limit = etree.SubElement(operator_starvation, "fractionBaryonLimitInNodeMerger")
            fraction_baryon_limit.set("value", fraction_baryon_value)
        # Add a reincorporation operator if needed.
        if not (hh_type == "coldMode" and outflow_to_cold_mode_value == "true"):
            operator_reinc = etree.Element("nodeOperator")
            operator_reinc.set("value", "CGMOutflowReincorporation")
            if is_grid:
                operator_reinc.set("iterable", "no")
            node_operators[0].append(operator_reinc)
            include_satellites_value = (
                "false" if starve_satellites_outflowed_value == "true" or starve_satellites_value == "true" else "true"
            )
            include_satellites = etree.SubElement(operator_reinc, "includeSatellites")
            include_satellites.set("value", include_satellites_value)


def hot_halo_standard_ram_pressure_stripping(input_doc, parameters, is_grid):
    """Special handling to add and modify nodeOperators for ram pressure stripping of the CGM in the standard hot halo component."""
    nodes = parameters.xpath(
        ".//componentHotHalo[@value='standard' or @value='coldMode']"
    )
    if len(nodes) <= 0:
        # None found - check if we have any componentHotHalo.
        nodes_any = parameters.xpath(".//componentHotHalo[@value]")
        if len(nodes_any) > 0:
            return
        # Check if we have nodeOperators present.
        node_ops = parameters.xpath(".//nodeOperator[@value='multi']")
        if len(node_ops) == 0:
            return
        # No componentHotHalo - insert one with the default value.
        component_hot_halo = etree.SubElement(parameters, "componentHotHalo")
        component_hot_halo.set("value", "standard")
        nodes = parameters.xpath(
            ".//componentHotHalo[@value='standard' or @value='coldMode' or @value='outflowTracking']"
        )
    if len(nodes) > 1:
        sys.exit(
            "found multiple `.//componentHotHalo[@value='standard' or ...]` nodes"
            " - unknown what should be done in this situation"
        )
    print(
        "   translate special './/componentHotHalo[@value='standard' or @value='coldMode']'"
    )
    # Determine the type of hot halo component.
    hh_type = nodes[0].get("value")
    # Find nodeOperators.
    node_operators = parameters.xpath(".//nodeOperator[@value='multi']")
    if len(node_operators) == 0:
        sys.exit("can not find any `nodeOperator[@value='multi']` into which to insert CGM operators")
    if len(node_operators) > 1:
        print(etree.tostring(parameters, encoding="unicode"))
        sys.exit("found multiple `nodeOperator[@value='multi']` nodes - unknown into which to insert CGM operators")
    # Add the ram pressure stripping operator.
    operator_accretion = etree.Element("nodeOperator")
    operator_accretion.set("value", "CGMOuterRadiusRamPressureStripping")
    operator_accretion.tail = "\n"
    if is_grid:
        operator_accretion.set("iterable", "no")
    node_operators[0].append(operator_accretion)
    # Handle the coldMode cases.
    if hh_type == "coldMode":
        # Find pre-existing options.
        def _get_prior_value(xpath_expr, default):
            prior = nodes[0].xpath(xpath_expr)
            return prior[0].get("value") if len(prior) == 1 else default

        outflow_to_cold_mode_value = _get_prior_value(".//outflowToColdMode[@value]", "false")
        # Remove obsoleted parameters.
        for elem in nodes[0].xpath(".//outflowToColdMode[@value]"):
            elem.getparent().remove(elem)
        # Set the option in the CGMOutflowReincorporation nodeOperator.
        operator_reincorporation = node_operators[0].xpath(".//nodeOperator[@value='CGMOutflowReincorporation']")
        if len(operator_reincorporation) == 0:
            sys.exit("can not find any `nodeOperator[@value='CGMOutflowReincorporation']` into which to insert `outflowToColdMode`")
        if len(operator_reincorporation) > 1:
            sys.exit("found multiple `nodeOperator[@value='CGMOutflowReincorporation']` nodes - unknown into which to insert `outflowToColdMode`")
        outflow_to_cold_mode = etree.SubElement(operator_reincorporation, "outflowToColdMode")
        outflow_to_cold_mode.set("value", outflow_to_cold_mode_value)

def hot_halo_standard_inflow_outflow(input_doc, parameters, is_grid):
    """Special handling to add and modify nodeOperators for CGM inflow/outflow in the standard hot halo component."""
    # Rename any `nodeOperatorCGMCoolingInflow` to `CGMCoolingHeating`.
    for op in parameters.xpath(".//nodeOperator[@value='multi']/nodeOperator[@value='CGMCoolingInflow']"):
        op.tag = "CGMCoolingHeating"
    # Look for "componentHotHalo" parameters.
    nodes = parameters.xpath(
        ".//componentHotHalo[@value='standard' or @value='coldMode' or @value='outflowTracking']"
    )
    if len(nodes) <= 0:
        nodes_any = parameters.xpath(".//componentHotHalo[@value]")
        if len(nodes_any) > 0:
            return
        node_ops = parameters.xpath(".//nodeOperator[@value='multi']")
        if len(node_ops) == 0:
            return
        component_hot_halo = etree.SubElement(parameters, "componentHotHalo")
        component_hot_halo.set("value", "standard")
        nodes = parameters.xpath(
            ".//componentHotHalo[@value='standard' or @value='coldMode' or @value='outflowTracking']"
        )
    if len(nodes) > 1:
        sys.exit(
            "found multiple `.//componentHotHalo[@value='standard' or ...]` nodes"
            " - unknown what should be done in this situation"
        )
    print(
        "   translate special './/componentHotHalo[@value='standard' or @value='coldMode'"
        " or @value='outflowTracking']'"
    )
    hh_type = nodes[0].get("value")
    # Find the component to which cooling should be directed.
    nodes_disk = parameters.xpath(".//componentDisk[@value]")
    if len(nodes_disk) == 0:
        component_cooling = "disk"
    else:
        type_disk = nodes_disk[0].get("value")
        component_cooling = "none" if type_disk == "null" else "disk"
    # Find nodeOperators.
    node_operators = parameters.xpath(".//nodeOperator[@value='multi']")
    if len(node_operators) == 0:
        sys.exit("can not find any `nodeOperator[@value='multi']` into which to insert CGM operators")
    if len(node_operators) > 1:
        sys.exit("found multiple `nodeOperator[@value='multi']` nodes - unknown into which to insert CGM operators")

    # Find pre-existing options.
    def _get_prior(xpath_expr, default):
        prior = nodes[0].xpath(xpath_expr)
        return prior[0].get("value") if len(prior) == 1 else default

    cooling_from_node_value = _get_prior(".//coolingFromNode[@value]", "currentNode")
    excess_heat_drives_outflow_value = _get_prior(".//hotHaloExcessHeatDrivesOutflow[@value]", "true")
    rate_maximum_expulsion_value = _get_prior(".//rateMaximumExpulsion[@value]", "1.0")
    fraction_loss_angular_momentum_value = _get_prior(".//fractionLossAngularMomentum[@value]", "0.3")
    efficiency_stripping_outflow_value = _get_prior(".//efficiencyStrippingOutflow[@value]", "0.1")
    track_stripped_gas_value = _get_prior(".//trackStrippedGas[@value]", "true")
    # Remove obsoleted parameters.
    for xpath_expr in (
        ".//coolingFromNode[@value]",
        ".//hotHaloExcessHeatDrivesOutflow[@value]",
        ".//rateMaximumExpulsion[@value]",
        ".//fractionLossAngularMomentum[@value]",
        ".//efficiencyStrippingOutflow[@value]",
        ".//trackStrippedGas[@value]",
    ):
        for elem in nodes[0].xpath(xpath_expr):
            elem.getparent().remove(elem)
    # Remove any blackHolesCGMHeating nodeOperator.
    bh_cgm_heatings = parameters.xpath(".//nodeOperator[@value='multi']/nodeOperator[@value='blackHolesCGMHeating']")
    if len(bh_cgm_heatings) > 0:
        cgm_heating = etree.Element("circumgalacticMediumHeating")
        cgm_heating.set("value", "AGNFeedback")
        if is_grid:
            cgm_heating.set("iterable", "no")
        insert_after(parameters, cgm_heating, nodes[0])
        for bh in bh_cgm_heatings:
            bh.getparent().remove(bh)
    # Add a cooling/heating operator.
    operator_cooling_heating = etree.Element("nodeOperator")
    operator_cooling_heating.set("value", "CGMCoolingHeating")
    if is_grid:
        operator_cooling_heating.set("iterable", "no")
    node_operators[0].append(operator_cooling_heating)
    component_elem = etree.SubElement(operator_cooling_heating, "component")
    cooling_from = etree.SubElement(operator_cooling_heating, "coolingFrom")
    excess_heat = etree.SubElement(operator_cooling_heating, "excessHeatDrivesOutflow")
    rate_max = etree.SubElement(operator_cooling_heating, "rateMaximumExpulsion")
    component_elem.set("value", component_cooling)
    cooling_from.set("value", cooling_from_node_value)
    excess_heat.set("value", excess_heat_drives_outflow_value)
    rate_max.set("value", rate_maximum_expulsion_value)
    # Add a cold mode inflow operator.
    if hh_type == "coldMode":
        operator_cold_mode = etree.Element("nodeOperator")
        operator_cold_mode.set("value", "CGMColdModeInflow")
        if is_grid:
            operator_cold_mode.set("iterable", "no")
        node_operators[0].append(operator_cold_mode)
        comp = etree.SubElement(operator_cold_mode, "component")
        cf = etree.SubElement(operator_cold_mode, "coolingFrom")
        comp.set("value", component_cooling)
        cf.set("value", cooling_from_node_value)
    # Add coolingInfallTorque object.
    cooling_infall_torque = etree.Element("coolingInfallTorque")
    cooling_infall_torque.set("value", "fixed")
    insert_after(parameters, cooling_infall_torque, nodes[0])
    fraction_loss = etree.SubElement(cooling_infall_torque, "fractionLossAngularMomentum")
    fraction_loss.set("value", fraction_loss_angular_momentum_value)
    # Add hotHaloOutflowStripping object.
    if track_stripped_gas_value == "true":
        outflow_stripping = etree.Element("hotHaloOutflowStripping")
        outflow_stripping.set("value", "standard")
        insert_after(parameters, outflow_stripping, nodes[0])
        efficiency = etree.SubElement(outflow_stripping, "efficiency")
        efficiency.set("value", efficiency_stripping_outflow_value)
    else:
        outflow_stripping = etree.Element("hotHaloOutflowStripping")
        outflow_stripping.set("value", "zero")
        insert_after(parameters, outflow_stripping, nodes[0])
    # Replace the "outflowTracking" component with a nodeOperator and nodePropertyExtractor.
    if hh_type == "outflowTracking":
        nodes[0].set("value", "standard")
        # Find nodePropertyExtractors, inserting one if none exists.
        extractors = parameters.xpath(".//nodePropertyExtractor[@value='multi']")
        if len(extractors) > 1:
            sys.exit("found multiple `nodePropertyExtractor[@value='multi']` nodes - unknown into which to insert CGM operators")
        if len(extractors) == 0:
            extractor_multi = etree.SubElement(parameters, "nodePropertyExtractor")
            extractor_indices = etree.Element("nodePropertyExtractor")
            extractor_multi.set("value", "multi")
            extractor_indices.set("value", "nodeIndices")
            extractor_multi.append(extractor_indices)
            extractors = parameters.xpath(".//nodePropertyExtractor[@value='multi']")
        # Insert nodeOperator.
        node_operator = etree.Element("nodeOperator")
        node_operator.set("value", "trackOutflowedMass")
        node_operators[0].append(node_operator)
        # Insert nodePropertyExtractor.
        node_extractor = etree.Element("nodePropertyExtractor")
        node_extractor.set("value", "trackOutflowedMass")
        extractors[0].append(node_extractor)


# ---------------------------------------------------------------------------
# Dispatch table for special migration functions
# ---------------------------------------------------------------------------

SPECIAL_FUNCTIONS = {
    "radiation_field_intergalactic_background_cmb": radiation_field_intergalactic_background_cmb,
    "black_hole_seed_mass": black_hole_seed_mass,
    "black_hole_physics": black_hole_physics,
    "model_parameter_xpath": model_parameter_xpath,
    "method_suffix_remove": method_suffix_remove,
    "satellite_orphanize": satellite_orphanize,
    "black_hole_non_central": black_hole_non_central,
    "hot_halo_very_simple": hot_halo_very_simple,
    "collaborative_mpi": collaborative_mpi,
    "hot_halo_standard_accretion": hot_halo_standard_accretion,
    "hot_halo_standard_inflow_outflow": hot_halo_standard_inflow_outflow,
    "hot_halo_standard_ram_pressure_stripping": hot_halo_standard_ram_pressure_stripping,
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    options = parse_arguments()

    # Write starting message.
    print(f"Translating file: {options.inputFile}")

    # Pre-process the file to concatenate any attribute values that are split across multiple lines.
    tmp_input = preprocess_multiline_attributes(options.inputFile)

    # Parse the input file.
    parser = etree.XMLParser(remove_blank_text=False)
    input_doc = etree.parse(tmp_input, parser)
    is_grid = False
    parameters_list = input_doc.xpath("/parameters")
    if parameters_list:
        root = parameters_list[0]
    else:
        parameter_grid_list = input_doc.xpath("/parameterGrid")
        if parameter_grid_list:
            is_grid = True
            root = parameter_grid_list[0]
            parameters_list = root.xpath("parameters")
            if not parameters_list:
                sys.exit("can not find <parameters>")
        else:
            sys.exit(f"can not find <parameters> in file `{options.inputFile}`")

    # Check if the file is under version control.
    is_in_git = git_is_tracked(options.inputFile)

    # Find the current hash.
    hash_head = git_head_hash()

    # Read migration rules.
    exec_path = os.environ.get("GALACTICUS_EXEC_PATH", ".")
    migrations = parse_migrations(os.path.join(exec_path, "scripts", "aux", "migrations.xml"))

    # Iterate over parameter sets.
    for parameters in parameters_list:
        migrate(input_doc, parameters, True, is_grid, options.inputFile, options, hash_head, is_in_git, migrations)

    # Output the resulting file.
    if options.prettyify == "yes":
        etree.indent(input_doc, space="  ")
    serialized = etree.tostring(input_doc, xml_declaration=True, encoding="UTF-8", pretty_print=False).decode("UTF-8")
    serialized = re.sub(r"><!--", ">\n\n  <!--", serialized)
    serialized = re.sub(r"><(?!/)", ">\n\n  <", serialized)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        tmp_output = f.name
        f.write(serialized)

    # If requested, ignore whitespace changes.
    if options.ignoreWhiteSpaceChanges == "yes":
        # Make a patch from the old to the new file, but ignoring changes in whitespace.
        with tempfile.NamedTemporaryFile(mode='w') as patch_file:
            subprocess.run(
                ["diff", "-w", "-u", tmp_input, tmp_output],
                stdout=patch_file,
            )
            # Apply the patch to the old file.
            subprocess.run(
                ["patch", tmp_input, patch_file.name, f"--output={tmp_output}"],
                stdout=subprocess.DEVNULL,
            )

    # Undo any split line reformatting that we previously applied.
    restore_multiline_attributes(tmp_output, options.outputFile)

    # Clean up.
    for f in [tmp_input, tmp_output]:
        if os.path.exists(f):
            os.unlink(f)


if __name__ == "__main__":
    main()
