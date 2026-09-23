#!/usr/bin/env python3
"""Regression tests for failures found by fuzzing the parameters of `parameters/quickTest.xml`.

Each case applies a parameter change file to `parameters/quickTest.xml` and runs the model. Every change here once crashed the
executable (a segmentation fault, or a Fortran runtime error), or ended with a failed task yet an exit status of zero. The
cases are of three kinds:

  * invalid configurations, which must be rejected with a clean `Fatal error:` naming the problem - so the model must exit
    with a non-zero status, must not crash, and must emit the expected diagnostic;
  * valid configurations, which must run to completion, and must then satisfy a check specific to the case;
  * configurations for which evolution fails (an ODE integration failure), which must report that the task failed and exit
    with a non-zero status. The failure here is a numerical one which may one day be fixed - if so, this case needs another
    configuration which causes a failure.

Andrew Benson
"""

import os
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import numpy as np

# Ensure output directory exists.
os.makedirs("outputs", exist_ok=True)

# Each case: (label, changes, expected outcome, OpenMP thread count). The expected outcome is one of:
#   ("rejected"  , [strings which must all appear in the output]),
#   ("runs"      , check) where `check` is a function returning a list of problems (empty if none), or
#   ("taskFailed", None).


def tagFor(label):
    """Return the tag used to name the files of the case with this label."""
    return "".join(c if c.isalnum() else "_" for c in label)


def outputFileNameFor(label):
    """Return the name of the model output file written by the case with this label."""
    return f"outputs/parameterFuzzingRegressions_{tagFor(label)}.hdf5"


# Labels of the cases whose checks must find their own output file.
labelMergerTreeMassUnion = "union of merger tree masses with two members"
labelOutputTimesUnion    = "union of output times with two members"


def checkMergerTreeMassUnion():
    """Check that the union built one tree from each of its two `fixedMass` members."""
    fileName = outputFileNameFor(labelMergerTreeMassUnion)
    if not os.path.exists(fileName):
        return [f"output file '{fileName}' was not written"]
    with h5py.File(fileName, "r") as file:
        trees  = file["Outputs/Output1/mergerTreeIndex"     ][:]
        masses = file["Outputs/Output1/nodeData/basicMass"  ][:]
    problems = []
    if len(trees) != 2:
        problems.append(f"expected 2 trees, one from each member of the union, but found {len(trees)}")
    for massExpected in (1.0e12, 1.0e13):
        if not np.any(np.isclose(masses, massExpected, rtol=1.0e-3)):
            problems.append(f"no halo of mass {massExpected:g} M_Solar - a member of the union contributed no tree")
    return problems


def checkOutputTimesUnion():
    """Check that the union output at the times of both of its `list` members."""
    fileName = outputFileNameFor(labelOutputTimesUnion)
    if not os.path.exists(fileName):
        return [f"output file '{fileName}' was not written"]
    with h5py.File(fileName, "r") as file:
        expansionFactors = sorted(file[f"Outputs/{name}"].attrs["outputExpansionFactor"] for name in file["Outputs"].keys())
    problems = []
    if len(expansionFactors) != 2:
        problems.append(f"expected 2 outputs, one from each member of the union, but found {len(expansionFactors)}")
    # The members ask for redshifts 1 and 0, i.e. expansion factors 0.5 and 1.
    for expansionFactorExpected in (0.5, 1.0):
        if not any(np.isclose(expansionFactors, expansionFactorExpected, rtol=1.0e-3)):
            problems.append(f"no output at expansion factor {expansionFactorExpected} - a member of the union contributed no time")
    return problems


def checkEvolutionOutput():
    """Check that the evolution output file is complete, well-formed XML with at least one node record."""
    fileName = "outputs/parameterFuzzingRegressions_evolution.xml"
    if not os.path.exists(fileName):
        return [f"evolution output file '{fileName}' was not written"]
    try:
        root = ET.parse(fileName).getroot()
    except ET.ParseError as error:
        return [f"evolution output file is not well-formed XML ({error})"]
    if root.tag != "evolution":
        return [f"evolution output file has root element <{root.tag}>, not <evolution>"]
    if len(root) == 0:
        return ["evolution output file contains no node records"]
    return []


cases = [
    (
        "cooling from the formation node without formation nodes (CGMCoolingHeating)",
        """  <change type="update" path="nodeOperator/nodeOperator[@value='CGMCoolingHeating']/coolingFrom" value="formationNode"/>
""",
        ("rejected", ["no formation node exists", "[coolingFrom]=formationNode in the [CGMCoolingHeating] node operator"]),
        1,
    ),
    (
        "cooling rate from the formation node without formation nodes (cole2000)",
        """  <change type="replace" path="coolingRate">
    <coolingRate value="cole2000"/>
  </change>
""",
        ("rejected", ["no formation node exists", "[coolingRate]=cole2000"]),
        1,
    ),
    (
        "union of merger tree masses with no members",
        """  <change type="replace" path="mergerTreeBuildMasses">
    <mergerTreeBuildMasses value="union"/>
  </change>
""",
        ("rejected", ["at least one [mergerTreeBuildMasses] must be specified"]),
        1,
    ),
    (
        labelMergerTreeMassUnion,
        """  <change type="replace" path="mergerTreeBuildMasses">
    <mergerTreeBuildMasses value="union">
      <mergerTreeBuildMasses value="fixedMass">
        <massTree value="1.0e12"/>
        <treeCount value="1"/>
      </mergerTreeBuildMasses>
      <mergerTreeBuildMasses value="fixedMass">
        <massTree value="1.0e13"/>
        <treeCount value="1"/>
      </mergerTreeBuildMasses>
    </mergerTreeBuildMasses>
  </change>
""",
        ("runs", checkMergerTreeMassUnion),
        1,
    ),
    (
        "union of output times with no members",
        """  <change type="replaceOrAppend" path="outputTimes">
    <outputTimes value="union"/>
  </change>
""",
        ("rejected", ["at least one [outputTimes] must be specified"]),
        1,
    ),
    (
        labelOutputTimesUnion,
        """  <change type="replaceOrAppend" path="outputTimes">
    <outputTimes value="union">
      <outputTimes value="list">
        <redshifts value="0.0"/>
      </outputTimes>
      <outputTimes value="list">
        <redshifts value="1.0"/>
      </outputTimes>
    </outputTimes>
  </change>
""",
        ("runs", checkOutputTimesUnion),
        1,
    ),
    (
        "constrained merger tree builder with no builders",
        """  <change type="replace" path="mergerTreeBuilder">
    <mergerTreeBuilder value="constrained"/>
  </change>
""",
        ("rejected", ["at least one [mergerTreeBuilder] must be specified"]),
        1,
    ),
    (
        # Both parameters of the extended ETHOS window function which are raised to a non-integer power must be positive.
        "negative beta1 in the extended ETHOS window function",
        """  <change type="replaceOrAppend" path="powerSpectrumWindowFunction">
    <powerSpectrumWindowFunction value="ETHOSExtended">
      <beta1 value="-1.3"/>
    </powerSpectrumWindowFunction>
  </change>
""",
        ("rejected", ["beta1] has value", "must be greater than 0"]),
        1,
    ),
    (
        # A zero unresolved energy ignores the energy of unresolved accretion, and should be usable. `factorMassResolution` is
        # lowered so that the energy model is applied to halos at `quickTest`'s resolution - otherwise it applies to none.
        "zero unresolved energy in the Johnson et al. (2021) scale radius model",
        """  <change type="replace" path="darkMatterProfileScaleRadius/darkMatterProfileScaleRadius">
    <darkMatterProfileScaleRadius value="johnson2021">
      <unresolvedEnergy     value="0.0"/>
      <factorMassResolution value="1.0"/>
      <darkMatterProfileScaleRadius value="concentration"/>
    </darkMatterProfileScaleRadius>
  </change>
""",
        ("runs", lambda: []),
        1,
    ),
    (
        "array parameter with the wrong number of values",
        """  <change type="replace" path="virialOrbit">
    <virialOrbit value="jiang2014">
      <gammaRatioLow value="0.2"/>
    </virialOrbit>
  </change>
""",
        ("rejected", ["parameter [gammaRatioLow] requires 3 values, but 1 was given"]),
        1,
    ),
    (
        # CAMB itself rejects this cosmology - its reason should be quoted in the error message.
        "CAMB given an invalid cosmology",
        """  <change type="replace" path="transferFunction">
    <transferFunction value="CAMB">
      <cosmologyParameters value="simple">
        <HubbleConstant  value="67.36"/>
        <OmegaMatter     value="13.4" />
        <OmegaDarkEnergy value="0.685"/>
        <OmegaBaryon     value="0.049"/>
        <temperatureCMB  value="2.725"/>
      </cosmologyParameters>
    </transferFunction>
  </change>
""",
        ("rejected", ["CAMB failed - its output ends:", "Your matter densities are strange"]),
        1,
    ),
    (
        # CLASS itself fails for this cosmology - its reason should be quoted in the error message.
        "CLASS given an invalid cosmology",
        """  <change type="replace" path="transferFunction">
    <transferFunction value="CLASSCDM">
      <cosmologyParameters value="simple">
        <HubbleConstant  value="67.36"/>
        <OmegaMatter     value="13.4" />
        <OmegaDarkEnergy value="0.685"/>
        <OmegaBaryon     value="0.049"/>
        <temperatureCMB  value="2.725"/>
      </cosmologyParameters>
    </transferFunction>
  </change>
""",
        ("rejected", ["CLASS failed - its output ends:", "Error in"]),
        1,
    ),
    (
        # Without `darkMatterProfileScaleSet` no node operator gives the halo a dark matter profile, so the enclosed mass within
        # `radiusTree` is zero whatever the halo mass, and no root can be found. The scale radius interpolation and halo spin
        # operators are removed too, as they would otherwise report the missing profile first.
        "tree radius with no dark matter profile",
        """  <change type="replace" path="mergerTreeBuildMasses">
    <mergerTreeBuildMasses value="fixedMass">
      <massTree   value="1.0e12"/>
      <radiusTree value="0.2"   />
      <treeCount  value="1"     />
    </mergerTreeBuildMasses>
  </change>
  <change type="remove" path="nodeOperator/nodeOperator[@value='darkMatterProfileScaleSet']"/>
  <change type="remove" path="nodeOperator/nodeOperator[@value='darkMatterProfileScaleInterpolate']"/>
  <change type="remove" path="nodeOperator/nodeOperator[@value='haloAngularMomentumRandom']"/>
  <change type="remove" path="nodeOperator/nodeOperator[@value='haloAngularMomentumInterpolate']"/>
""",
        ("rejected", ["the halo has no dark matter mass distribution, so [radiusTree]",
                      "[nodeOperator]=darkMatterProfileScaleSet"]),
        1,
    ),
    (
        "halo spin with no dark matter profile",
        """  <change type="remove" path="nodeOperator/nodeOperator[@value='darkMatterProfileScaleSet']"/>
  <change type="remove" path="nodeOperator/nodeOperator[@value='darkMatterProfileScaleInterpolate']"/>
""",
        ("rejected", ["the halo has no dark matter mass distribution, so its angular momentum scale",
                      "[nodeOperator]=darkMatterProfileScaleSet"]),
        1,
    ),
    (
        # Nested inside `filteredMainBranch` the scale radius is set only for halos on the main branch.
        "dark matter profile scale radius never set",
        """  <change type="replace" path="nodeOperator/nodeOperator[@value='darkMatterProfileScaleSet']">
    <nodeOperator value="filteredMainBranch">
      <nodeOperator value="darkMatterProfileScaleSet"/>
    </nodeOperator>
  </change>
  <change type="remove" path="nodeOperator/nodeOperator[@value='darkMatterProfileScaleInterpolate']"/>
  <change type="remove" path="nodeOperator/nodeOperator[@value='haloAngularMomentumRandom']"/>
  <change type="remove" path="nodeOperator/nodeOperator[@value='haloAngularMomentumInterpolate']"/>
""",
        ("rejected", ["requires a positive scale radius", "has not been set"]),
        1,
    ),
    (
        "zero dark matter profile scale radius with an NFW profile",
        """  <change type="replace" path="darkMatterProfileScaleRadius">
    <darkMatterProfileScaleRadius value="zero"/>
  </change>
""",
        ("rejected", ["requires a positive scale radius", "[darkMatterProfileScaleRadius]"]),
        1,
    ),
    (
        # `quickTest` uses the Eisenstein & Hu (1999) transfer function, which has no small-scale cutoff.
        "prompt cusps with no small-scale cutoff in the power spectrum",
        """  <change type="append" path="nodeOperator">
    <nodeOperator value="darkMatterProfilePromptCusps"/>
  </change>
""",
        ("rejected", ["of the power spectrum does not converge", "small-scale cutoff"]),
        1,
    ),
    (
        # `accretionHalo` is changed too, so that something asks the IGM for its state: nothing in `quickTest` does otherwise.
        "internal IGM state without the universe operator which solves for it",
        """  <change type="replaceOrAppend" path="intergalacticMediumState">
    <intergalacticMediumState value="internal"/>
  </change>
  <change type="replace" path="accretionHalo">
    <accretionHalo value="naozBarkana2007"/>
  </change>
""",
        ("rejected", ["the state of the IGM has not been solved for",
                      "[universeOperator]=intergalacticMediumStateEvolve"]),
        1,
    ),
    (
        "evolution output node operator",
        """  <change type="append" path="nodeOperator">
    <nodeOperator value="evolutionOutput">
      <outputFileName value="testSuite/outputs/parameterFuzzingRegressions_evolution.xml"/>
    </nodeOperator>
  </change>
""",
        ("runs", checkEvolutionOutput),
        # Use several threads, so that per-thread copies of the node operator must share the output file.
        4,
    ),
    (
        "mass movement destination which the component does not implement",
        """  <change type="replaceOrAppend" path="mergerMassMovements">
    <mergerMassMovements value="simple">
      <destinationStarsMinorMerger value="unmoved"/>
    </mergerMassMovements>
  </change>
""",
        ("rejected", ["the `unmoved` destination of [destinationStarsSatellite] is not supported"]),
        1,
    ),
    (
        "star formation rate set for a disk which does not exist",
        """  <change type="replaceOrAppend" path="starFormationRateDisks">
    <starFormationRateDisks value="fixed"/>
  </change>
""",
        ("rejected", ["a non-zero rate was set for the `massGas` property of the `disk` component"]),
        1,
    ),
    (
        "operator which needs a component that is not present",
        """  <change type="replaceOrAppend" path="nodeOperator/nodeOperator[@value='CGMOuterRadiusRamPressureStripping']">
    <nodeOperator value="positionToHost"/>
  </change>
""",
        ("rejected", ["attempt to set the `position` property of the null `position` component"]),
        1,
    ),
    (
        "SIDM satellite evaporation without satellite orbits",
        """  <change type="replaceOrAppend" path="nodeOperator/nodeOperator[@value='CGMCoolingHeating']">
    <nodeOperator value="satelliteEvaporationSIDM"/>
  </change>
""",
        ("rejected", ['the "satellite" component must provide gettable "position" and "velocity"']),
        1,
    ),
    (
        "satellite tidal mass loss without satellite orbits",
        """  <change type="replaceOrAppend" path="nodeOperator/nodeOperator[@value='DMOInterpolate']">
    <nodeOperator value="satelliteTidalMassLoss"/>
  </change>
""",
        ("rejected", ['the "satellite" component must provide gettable "position" and "velocity"']),
        1,
    ),
    (
        "primordial power spectrum which diverges",
        """  <change type="update" path="powerSpectrumPrimordial/running" value="0.807252"/>
""",
        ("rejected", ["primordial power spectrum overflows at wavenumber"]),
        1,
    ),
    (
        # The original configuration for this case, a [rateMaximumExpulsion] of 1e30, no longer causes evolution to fail, so
        # tolerances far below what the solver can achieve are used instead.
        "failed evolution",
        """  <change type="update" path="mergerTreeNodeEvolver/odeToleranceAbsolute" value="1.0e-100"/>
  <change type="update" path="mergerTreeNodeEvolver/odeToleranceRelative" value="1.0e-100"/>
""",
        ("taskFailed", None),
        1,
    ),
    (
        "failed evolution with failures tolerated",
        """  <change type="update" path="mergerTreeNodeEvolver/odeToleranceAbsolute" value="1.0e-100"/>
  <change type="update" path="mergerTreeNodeEvolver/odeToleranceRelative" value="1.0e-100"/>
  <change type="replaceOrAppend" path="task">
    <task value="evolveForests">
      <tolerateFailures value="true"/>
    </task>
  </change>
""",
        ("runs", lambda: []),
        1,
    ),
]


def run(label, changes, countThreads):
    """Run quickTest with the given changes and number of OpenMP threads, returning (return code, output)."""
    tag            = tagFor(label)
    changeFileName = f"outputs/parameterFuzzingRegressions_{tag}.changes.xml"
    with open(changeFileName, "w") as changeFile:
        changeFile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<changes>\n")
        changeFile.write(changes)
        changeFile.write(f"""  <change type="replaceOrAppend" path="outputFileName">
    <outputFileName value="testSuite/outputs/parameterFuzzingRegressions_{tag}.hdf5"/>
  </change>
""")
        changeFile.write("</changes>\n")
    environment                    = dict(os.environ)
    environment["OMP_NUM_THREADS"] = str(countThreads)
    try:
        process = subprocess.run(
            f"cd ..; ./Galacticus.exe parameters/quickTest.xml testSuite/{changeFileName}",
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=600,
            env=environment,
        )
    except subprocess.TimeoutExpired:
        return None, ""
    return process.returncode, process.stdout or ""


def crashed(returnCode, output):
    """Return true if the model crashed rather than stopping cleanly."""
    # Galacticus' own signal handlers report the signal ("Galacticus experienced a segfault", "...a floating point
    # exception", etc.) and exit with the signal number as the status. An unhandled signal is reported by the shell as 128+N, or
    # as a negative number if not run through a shell. A Fortran runtime error exits with status 2 but writes no `Fatal error:`
    # message.
    return (
        returnCode < 0
        or returnCode > 128
        or "Galacticus experienced" in output
        or "Fortran runtime error" in output
        or "Segmentation fault" in output
    )


failed = False
for label, changes, (kind, expectation), countThreads in cases:
    returnCode, output = run(label, changes, countThreads)
    if returnCode is None:
        print(f"FAILED: {label}: model did not finish within the time limit")
        failed = True
        continue
    if crashed(returnCode, output):
        print(f"FAILED: {label}: model crashed (return code {returnCode})")
        print(output[-3000:])
        failed = True
        continue
    if kind == "rejected":
        if returnCode == 0:
            print(f"FAILED: {label}: model ran to completion, but the configuration is invalid and should have been rejected")
            failed = True
            continue
        missing = [expected for expected in expectation if expected not in output]
        if "Fatal error:" not in output or missing:
            print(f"FAILED: {label}: model stopped (return code {returnCode}) but without the expected diagnostic {missing}")
            print(output[-3000:])
            failed = True
            continue
        print(f"SUCCESS: {label}: rejected cleanly (return code {returnCode})")
    elif kind == "taskFailed":
        if returnCode == 0:
            print(f"FAILED: {label}: model exited with status 0, but its task should have failed")
            print(output[-3000:])
            failed = True
            continue
        if "task failed" not in output:  # markers: exempt
            print(f"FAILED: {label}: model exited with status {returnCode}, but did not report a failed task")
            print(output[-3000:])
            failed = True
            continue
        print(f"SUCCESS: {label}: task failure reported, with non-zero exit status {returnCode}")
    else:
        if returnCode != 0:
            print(f"FAILED: {label}: model did not run to completion (return code {returnCode})")
            print(output[-3000:])
            failed = True
            continue
        problems = expectation()
        if problems:
            for problem in problems:
                print(f"FAILED: {label}: {problem}")
            failed = True
            continue
        print(f"SUCCESS: {label}: ran to completion")

# Always exit with status 0 - failure is signaled by "FAILED" in the output above.
if failed:
    print("FAILED: parameter fuzzing regression tests")
else:
    print("SUCCESS: parameter fuzzing regression tests")
sys.exit(0)
