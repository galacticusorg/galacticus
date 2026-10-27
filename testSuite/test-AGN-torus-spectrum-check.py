#!/usr/bin/env python3
"""Check that an AGN spectrum which already includes the emission of a dusty torus can not be attenuated by a torus.

The default accretion disk spectrum, `accretionDiskSpectraHopkins2007`, is that of AGN as observed, and so includes the
infrared emission of their tori. Attenuating it by a model of the torus would count that emission twice, so property
extractors refuse the combination at construction. Starting from the AGN torus misalignment model, with that spectrum,
this checks that:

* a `dustAttenuation` extractor whose attenuator is a sequence containing a torus is refused;
* a `SEDDustEmission` extractor whose attenuator is a torus is refused;
* a `dustAttenuation` extractor whose attenuator has no torus is accepted, and the model runs.

Andrew Benson, Claude
"""
import os
import subprocess
import sys

from lxml import etree

# Text of the error reported when the combination is refused.
MESSAGE = "already includes the emission of the dusty torus"

# This script runs with the working directory set to testSuite/, while the model is run from the repository root.
parameterFileBase = "parameters/agnTorusMisalignment.xml"
os.makedirs("outputs", exist_ok=True)

# Extractors for each case: the attenuator, and whether the model must be refused.
SCREEN = """
<dustAttenuation value="screenSurfaceDensityMetals">
  <coefficient value="1.0"/>
  <dustExtinctionCurve value="powerLaw">
    <exponent value="0.7"/>
    <wavelengthReference value="5500.0"/>
  </dustExtinctionCurve>
</dustAttenuation>
"""
SED_AGN = """
<nodePropertyExtractor value="SEDAGN">
  <frame value="rest"/>
  <wavelengthMinimum value="8.0e2"/>
  <wavelengthMaximum value="4.0e4"/>
  <resolution value="50.0"/>
</nodePropertyExtractor>
"""
CASES = (
    ("sequence containing a torus, in a dustAttenuation extractor", True, f"""
<nodePropertyExtractor value="dustAttenuation">
  <dustAttenuation value="sequence">
    <dustAttenuation value="AGNTorus">
      <galacticInclination value="fixed">
        <inclination value="60.0"/>
      </galacticInclination>
      <dustExtinctionCurve value="powerLaw">
        <exponent value="0.7"/>
        <wavelengthReference value="5500.0"/>
      </dustExtinctionCurve>
    </dustAttenuation>
    {SCREEN}
  </dustAttenuation>
  {SED_AGN}
</nodePropertyExtractor>
"""),
    ("SKIRTOR torus, in a SEDDustEmission extractor", True, f"""
<nodePropertyExtractor value="SEDDustEmission">
  <frame value="rest"/>
  <wavelengthMinimum value="1.0e3"/>
  <wavelengthMaximum value="1.0e8"/>
  <resolution value="50.0"/>
  <dustAttenuation value="AGNTorusSKIRTOR">
    <galacticInclination value="fixed">
      <inclination value="60.0"/>
    </galacticInclination>
  </dustAttenuation>
  <dustEmissionSpectrum value="AGNTorusSKIRTOR"/>
  {SED_AGN}
</nodePropertyExtractor>
"""),
    ("host-galaxy screen without a torus, in a dustAttenuation extractor", False, f"""
<nodePropertyExtractor value="dustAttenuation">
  {SCREEN}
  {SED_AGN}
</nodePropertyExtractor>
"""),
)

for index, (label, refused, extractor) in enumerate(CASES):
    parameters = etree.parse(parameterFileBase)
    root = parameters.getroot()
    # Use the default spectrum, which includes the emission of the torus.
    spectra = root.find("accretionDiskSpectra")
    if spectra is None:
        print(f"FAILED: [{label}] no accretionDiskSpectra parameter found in {parameterFileBase}")
        sys.exit(0)
    spectra.set("value", "hopkins2007")
    # Replace the extractors with those of this case.
    root.remove(root.find("nodePropertyExtractor"))
    root.append(etree.fromstring(extractor))
    outputName = f"testSuite/outputs/agnTorusSpectrumCheck{index}.hdf5"
    root.find("outputFileName").set("value", outputName)
    parameterFile = f"outputs/agnTorusSpectrumCheck{index}.xml"
    parameters.write(parameterFile)
    outputPath = os.path.join("..", outputName)
    if os.path.exists(outputPath):
        os.remove(outputPath)
    status = subprocess.run(f"cd ..; ./Galacticus.exe testSuite/{parameterFile}", shell=True, capture_output=True)
    log = status.stdout.decode() + status.stderr.decode()
    if refused:
        if status.returncode != 0 and MESSAGE in log:
            print(f"SUCCESS: [{label}] refused, as the spectrum includes the emission of the torus")
        else:
            print(f"FAILED: [{label}] was not refused, although the spectrum includes the emission of the torus")
    else:
        if status.returncode == 0 and os.path.exists(outputPath):
            print(f"SUCCESS: [{label}] accepted, and the model ran")
        else:
            print(f"FAILED: [{label}] did not run, although no torus attenuates the spectrum")
            print(log[-2000:])
sys.exit(0)
