#!/usr/bin/env python3
"""Generate parameter files for the halo mass function validation models.

This wraps `haloMassFunctionGenerateContent.py` to generate base parameter files
for a representative sample of the halo mass function calibration models of
Benson et al. (2026; https://ui.adsabs.harvard.edu/abs/2026arXiv260612137B).
These models are used for validation in the CI/CD pipeline (see
`testSuite/validate-haloMassFunction.py`).

Run from the repository root:

  export GALACTICUS_DATA_PATH=/path/to/datasets
  PYTHONPATH=$PWD/python python3 constraints/pipelines/darkMatter/haloMassFunctionValidationGenerate.py

This regenerates, in `testSuite/parameters/validation/haloMassFunction/`:

  * `haloMassFunctionBase_*.xml`  - base parameter files (one per simulation/redshift);
  * `haloMassFunction_<suite>.xml` - halo mass function model definitions per suite;
  * `manifest.json`               - per-case likelihood metadata (target data file,
                                    mass range, etc.) used by the validation driver.

The best-fit parameter file `haloMassFunctionParameters.xml` in the same directory
is *not* regenerated - it holds the calibrated parameter values and is maintained
by hand.
"""

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import tempfile

import lxml.etree as ET

# The validation groups. Each group corresponds to a single CI/CD job, driver
# invocation, and published validation metric. Selection strings have the format
# `suite::group::resolution::simulation::realization::redshift` with
# comma-separated alternatives at each level and `*` as a wildcard.
validationGroups = {
    "MDPL": [
        "MDPL::*::*::*::*::0.000",
        "MDPL::MDPL2::*::*::*::0.987,3.127",
        "MDPL::HugeMDPL::*::*::*::0.987,3.037",
    ],
    "SymphonyLMC": [
        "Symphony::LMC::resolutionX1::CDM::*::0.000",
    ],
    "SymphonyGroup": [
        "Symphony::Group::resolutionX1::CDM::*::0.000",
    ],
    "SymphonyMilkyWayZ0": [
        "Symphony::MilkyWay::resolutionX1::CDM::*::0.000",
        "Symphony::MilkyWay::resolutionX64::CDM::*::0.000",
    ],
    "SymphonyMilkyWayZ1": [
        "Symphony::MilkyWay::resolutionX1::CDM::*::0.990",
    ],
    "SymphonyMilkyWayZ4": [
        "Symphony::MilkyWay::resolutionX1::CDM::*::3.984",
    ],
    "COZMIC": [
        "COZMIC::MilkyWay::resolutionX8::WDM:3keV,WDM:3keV:f0.6,WDM:3keV:bumpCutoff,IDM:1e-4GeV:envelope,FDM:25.9e-22eV::*::0.000",
    ],
}

pathPipeline = "constraints/pipelines/darkMatter/"
pathTarget   = "testSuite/parameters/validation/haloMassFunction/"
pathOutputs  = "testSuite/outputs/validation/haloMassFunction/"
# Relative path from the target directory back to the pipeline directory - used
# to rewrite XInclude hrefs, which Galacticus resolves relative to the including
# file.
pathPipelineRelative = "../../../../constraints/pipelines/darkMatter/"


def parseLikelihoods(fileNameConfig):
    """Extract per-case likelihood metadata from a generated MCMC config file."""
    cases  = []
    config = ET.parse(fileNameConfig).getroot()
    for likelihood in config.iter("posteriorSampleLikelihood"):
        if likelihood.get("value") != "haloMassFunction":
            continue
        fileNamesTarget  = likelihood.find("fileNames"       ).get("value").split()
        redshifts        = likelihood.find("redshifts"       ).get("value").split()
        massRangeMinimum = float(likelihood.find("massRangeMinimum").get("value"))
        massRangeMaximum = None
        if likelihood.find("massRangeMaximum") is not None:
            massRangeMaximum = float(likelihood.find("massRangeMaximum").get("value"))
        fileNameBase = os.path.basename(likelihood.find("baseParametersFileName").get("value"))
        # The config names only the first base parameter file; files for the other
        # redshifts follow the same naming convention.
        match = re.match(r"^(haloMassFunctionBase_.+_z)\d+\.\d+\.xml$", fileNameBase)
        if match is None:
            raise RuntimeError(f"unable to parse base parameter file name '{fileNameBase}'")
        for redshift, fileNameTarget in zip(redshifts, fileNamesTarget):
            case = {
                "parameterFile"   : match.group(1)+redshift+".xml",
                "targetData"      : fileNameTarget,
                "redshift"        : float(redshift),
                "massRangeMinimum": massRangeMinimum,
            }
            if massRangeMaximum is not None:
                case["massRangeMaximum"] = massRangeMaximum
            cases.append(case)
    return cases


def postProcess(content, pathGenerated):
    """Rewrite paths in a generated base parameter file for its new location."""
    # Model output files go to the test suite outputs directory.
    content = re.sub(
        r'(<outputFileName\s+value=")'+re.escape(pathGenerated),
        r"\g<1>"+pathOutputs,
        content,
    )
    # Suite halo mass function definitions live alongside the base files.
    content = content.replace('href="'+pathGenerated, 'href="')
    # Pipeline includes are found relative to the target directory.
    content = content.replace('href="'+pathPipeline, 'href="'+pathPipelineRelative)
    return content


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--group", action="append", default=None,
        help="validation group to (re)generate (may be repeated; default: all)",
    )
    args   = parser.parse_args()
    groups = args.group if args.group is not None else list(validationGroups.keys())
    for group in groups:
        if group not in validationGroups:
            sys.exit(f"unknown validation group '{group}' - known groups: {', '.join(validationGroups)}")
    if not os.path.isdir(pathPipeline) or not os.path.isdir("testSuite"):
        sys.exit("this script must be run from the Galacticus repository root")
    if "GALACTICUS_DATA_PATH" not in os.environ:
        sys.exit("GALACTICUS_DATA_PATH must be set (the generator reads simulation data files)")

    # Read any existing manifest so that single-group regeneration preserves other groups.
    fileNameManifest = pathTarget+"manifest.json"
    manifest         = {"groups": {}}
    if os.path.exists(fileNameManifest):
        with open(fileNameManifest) as file:
            manifest = json.load(file)

    for group in groups:
        print(f"Generating validation group '{group}'...")
        with tempfile.TemporaryDirectory(prefix="haloMassFunctionValidation") as pathGenerated:
            # Ensure a trailing slash - the generator and the generated content assume one.
            pathGenerated = pathGenerated+"/"
            command       = [
                sys.executable, pathPipeline+"haloMassFunctionGenerateContent.py",
                "--pipelinePath"     , pathPipeline ,
                "--outputDirectory"  , pathGenerated,
                # The accelerator is an interpolation cache used to speed up MCMC
                # sampling - it is not needed for single model evaluations.
                "--removeAccelerator", "true"       ,
            ]
            for select in validationGroups[group]:
                command.extend(["--select", select])
            status = subprocess.run(command, capture_output=True, text=True)
            if status.returncode != 0:
                print(status.stdout)
                print(status.stderr, file=sys.stderr)
                sys.exit(f"generation failed for validation group '{group}'")
            # Extract likelihood metadata for this group.
            cases = parseLikelihoods(pathGenerated+"haloMassFunctionConfig.xml")
            # Remove base parameter files from any prior generation of this group.
            if group in manifest["groups"]:
                for case in manifest["groups"][group]["cases"]:
                    fileNameOld = pathTarget+case["parameterFile"]
                    if os.path.exists(fileNameOld):
                        os.unlink(fileNameOld)
            # Install post-processed base parameter files, and suite halo mass
            # function definitions.
            for fileName in sorted(glob.glob(pathGenerated+"haloMassFunctionBase_*.xml")):
                with open(fileName) as file:
                    content = file.read()
                with open(pathTarget+os.path.basename(fileName), "w") as file:
                    file.write(postProcess(content, pathGenerated))
            for fileName in sorted(glob.glob(pathGenerated+"haloMassFunction_*.xml")):
                if re.search(r"haloMassFunction(Config|ConfigResume|Parameters)\.xml$", fileName):
                    continue
                with open(fileName) as file:
                    content = file.read()
                with open(pathTarget+os.path.basename(fileName), "w") as file:
                    file.write(content)
            manifest["groups"][group] = {"cases": cases}
            print(f"   ...generated {len(cases)} cases")

    # The model discrepancy variance is defined in the best-fit parameters file -
    # record it in the manifest for use by the validation driver.
    parameters = ET.parse(pathTarget+"haloMassFunctionParameters.xml").getroot()
    manifest["varianceFractionalModelDiscrepancy"] = float(
        parameters.find("varianceFractionalModelDiscrepancy").get("value")
    )
    with open(fileNameManifest, "w") as file:
        json.dump(manifest, file, indent=2)
        file.write("\n")
    print("Done.")


if __name__ == "__main__":
    main()
