#!/usr/bin/env python3
"""Generate Galacticus parameter files for the mass-resolution scaling experiment.

Builds a fully-resolved (XInclude-expanded) parameter document from the
reference dark-matter-only subhalo model, then applies per-run edits. Usage:

  python3 make_params.py <outputParamFile> <massResolution> <treeCount> <outputHDF5> [key=value ...]

Recognized extra keys:
  fractionTimestepSatelliteMinimum, backtrackToSatellites, profileOdeEvolver,
  reuseODEStepSize, timeOffsetMaximumAbsolute, timeOffsetMaximumRelative,
  profiler (set to 'simple' to add an evolve profiler), massTree, seed,
  timestepHostRelative, outputRedshifts (colon-separated list, adds standard
  outputter for census runs), verbosity
"""
import sys
from lxml import etree

REFDIR = "/home/user/galacticus/parameters/reference"

SKELETON = f"""<?xml version='1.0' encoding='UTF-8'?>
<parameters>
  <formatVersion>2</formatVersion>
  <verbosityLevel value="warn"/>
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/darkMatterParticleCDM.xml"       xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/cosmologyDarkMatterOnly.xml"     xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/powerSpectrum.xml"               xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/structureFormation.xml"          xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/darkMatterHalosProfile.xml"      xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/darkMatterHalosTidalHeating.xml" xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/darkMatterHalosStructure.xml"    xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/subhaloOrbits.xml"               xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/mergerTrees.xml"                 xpointer="xpointer(parameters/*)" />
  <xi:include xmlns:xi="http://www.w3.org/2001/XInclude" href="{REFDIR}/evolutionDarkMatterOnly.xml"     xpointer="xpointer(parameters/*)" />

  <randomNumberGenerator value="GSL">
    <seed value="8122"/>
  </randomNumberGenerator>

  <task                   value="evolveForests"/>
  <evolveForestsWorkShare value="cyclic"       />

  <mergerTreeBuildMasses value="fixedMass">
    <massTree  value="1.0e13"/>
    <treeCount value="1"     />
  </mergerTreeBuildMasses>
  <mergerTreeMassResolution value="fixed">
    <massResolution value="1.0e8"/>
  </mergerTreeMassResolution>

  <mergerTreeOperator value="treeProcessingTimer"/>

  <outputFileName value="out.hdf5"/>
  <mergerTreeOutputter value="null"/>
  <outputTimes value="list">
    <redshifts value="0.0"/>
  </outputTimes>
</parameters>
"""


def find1(root, path):
    nodes = root.findall(path)
    if len(nodes) != 1:
        raise RuntimeError(f"expected exactly one node at {path}, found {len(nodes)}")
    return nodes[0]


def main():
    outParam, massRes, treeCount, outHDF5 = sys.argv[1:5]
    extra = dict(kv.split("=", 1) for kv in sys.argv[5:])

    parser = etree.XMLParser(remove_blank_text=False)
    doc = etree.fromstring(SKELETON.encode(), parser).getroottree()
    doc.xinclude()
    root = doc.getroot()

    # Replace the CAMB transfer function with the analytic Eisenstein & Hu form:
    # numerically almost identical for this purpose, and avoids a long CAMB
    # computation (approved by user for this performance study).
    tf = find1(root, "./transferFunction")
    if tf.get("value") == "CAMB":
        for child in list(tf):
            tf.remove(child)
        tf.set("value", "eisensteinHu1999")
        for name, value in [("neutrinoNumberEffective", "3.046"), ("neutrinoMassSummed", "0.000")]:
            el = etree.SubElement(tf, name)
            el.set("value", value)

    find1(root, "./mergerTreeMassResolution/massResolution").set("value", massRes)
    find1(root, "./mergerTreeBuildMasses/treeCount").set("value", treeCount)
    find1(root, "./outputFileName").set("value", outHDF5)
    if "massTree" in extra:
        find1(root, "./mergerTreeBuildMasses/massTree").set("value", extra["massTree"])
    if "seed" in extra:
        find1(root, "./randomNumberGenerator/seed").set("value", extra["seed"])
    if "verbosity" in extra:
        find1(root, "./verbosityLevel").set("value", extra["verbosity"])
    if "fractionTimestepSatelliteMinimum" in extra:
        find1(root, "./mergerTreeEvolver/fractionTimestepSatelliteMinimum").set("value", extra["fractionTimestepSatelliteMinimum"])
    if "backtrackToSatellites" in extra:
        find1(root, "./mergerTreeEvolver/backtrackToSatellites").set("value", extra["backtrackToSatellites"])
    if "timestepHostRelative" in extra:
        find1(root, "./mergerTreeEvolver/timestepHostRelative").set("value", extra["timestepHostRelative"])
    if "reuseODEStepSize" in extra:
        find1(root, "./mergerTreeNodeEvolver/reuseODEStepSize").set("value", extra["reuseODEStepSize"])
    if "profileOdeEvolver" in extra:
        ev = find1(root, "./mergerTreeNodeEvolver")
        el = etree.SubElement(ev, "profileOdeEvolver")
        el.set("value", extra["profileOdeEvolver"])
    if "timeOffsetMaximumAbsolute" in extra or "timeOffsetMaximumRelative" in extra:
        # The satellite timestep lives inside the multi timestep container.
        multi = find1(root, "./mergerTreeEvolveTimestep")
        sat = [c for c in multi.findall("mergerTreeEvolveTimestep") if c.get("value") == "satellite"]
        assert len(sat) == 1
        if "timeOffsetMaximumAbsolute" in extra:
            find1(sat[0], "./timeOffsetMaximumAbsolute").set("value", extra["timeOffsetMaximumAbsolute"])
        if "timeOffsetMaximumRelative" in extra:
            find1(sat[0], "./timeOffsetMaximumRelative").set("value", extra["timeOffsetMaximumRelative"])
    if extra.get("profiler") == "simple":
        el = etree.SubElement(root, "mergerTreeEvolveProfiler")
        el.set("value", "simple")
        sub = etree.SubElement(el, "timeStepMinimum")
        sub.set("value", "1.0e-8")
    if "outputRedshifts" in extra:
        find1(root, "./outputTimes/redshifts").set("value", extra["outputRedshifts"].replace(":", " "))
        outp = find1(root, "./mergerTreeOutputter")
        outp.set("value", "standard")
        sub = etree.SubElement(outp, "outputReferences")
        sub.set("value", "false")

    etree.indent(root, space="  ")
    doc.write(outParam, xml_declaration=True, encoding="UTF-8", pretty_print=True)
    print(f"wrote {outParam}")


if __name__ == "__main__":
    main()
