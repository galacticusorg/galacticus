#!/usr/bin/env python3
"""
Format a Galacticus parameter file to make it more easily comprehensible.
Port of parametersFormat.pl (Andrew Benson, 08-March-2016)
"""

import re
import sys
import xml.etree.ElementTree as ET

# Groups define the desired ordering of parameters. Each group has a name,
# a human-readable description (written as a section comment), and a list of
# members. A member is either an exact element tag name (string) or a compiled
# regular expression. The first matching group wins, and within a group
# parameters appear in member-list order (alphabetically for ties).
GROUPS = [
    {
        "name": "logging",
        "description": "Logging",
        "members": [
            "verbosityLevel",
        ],
    },
    {
        "name": "cosmology",
        "description": "Cosmological model",
        "members": [
            "cosmologyParameters",
            "cosmologyFunctions",
            re.compile(r"cosmology\d"),
        ],
    },
    {
        "name": "powerSpectrum",
        "description": "Power spectrum",
        "members": [
            "powerSpectrum",
            "powerSpectrumPrimordial",
            "transferFunction",
            "powerSpectrumPrimordialTransferred",
            "cosmologicalMassVariance",
            "powerSpectrumNonlinear",
        ],
    },
    {
        "name": "structureFormation",
        "description": "Structure formation",
        "members": [
            "linearGrowth",
            "criticalOverdensity",
            "virialDensityContrast",
            re.compile(r"virialDensityContrast.*"),
            "haloMassFunction",
            re.compile(r"haloMassFunctionSystematic\d"),
            "gravitationalLensing",
        ],
    },
    {
        "name": "igmReionization",
        "description": "Intergalactic medium and reionization",
        "members": [
            "intergalacticMediumState",
            re.compile(r"intergalaticMediumState.*"),
        ],
    },
    {
        "name": "components",
        "description": "Tree node component selection",
        "members": [
            re.compile(r"^component"),
            "nodeComponentBasicExtendedSphericalCollapseType",
        ],
    },
    {
        "name": "mergerTreeConstruction",
        "description": "Merger tree construction",
        "members": [
            "mergerTreeConstruct",
            "mergerTreeBuilder",
            "mergerTreeMassResolution",
            re.compile(r"^mergerTreeBuild.*"),
            re.compile(r"^mergerTreesBuild.*"),
            "treeBranching",
            re.compile(r"^modifiedPressSchechter.*"),
            "haloMassFunctionSampling",
            re.compile(r"^haloMassFunctionSampling.*"),
            re.compile(r"^mergerTreeRead.*"),
            re.compile(r"mergerTreeImport.*"),
        ],
    },
    {
        "name": "mergerTreeOperators",
        "description": "Merger tree operators",
        "members": [
            "mergerTreePruneBaryons",
        ],
    },
    {
        "name": "haloHierarchy",
        "description": "Merger tree halo hierarchy",
        "members": [
            "nodeMergers",
            "nodePromotionIndexShift",
        ],
    },
    {
        "name": "darkMatterHalo",
        "description": "Dark matter halo structure",
        "members": [
            "darkMatterProfile",
            "darkMatterProfileConcentration",
            "darkMatterProfileMinimumConcentration",
            "darkMatterProfileScaleCorrectForConcentrationDefinition",
        ],
    },
    {
        "name": "haloSpin",
        "description": "Dark matter halo spin",
        "members": [
            "haloSpinDistribution",
            "randomSpinResetMassFactor",
        ],
    },
    {
        "name": "satelliteOrbits",
        "description": "Satellite orbits",
        "members": [
            "virialOrbit",
            re.compile(r"^virialOrbits.*"),
            "satelliteMergingTimescales",
            "mergingTimescaleMultiplier",
            "satellitesTidalField",
            "satelliteOrbitStoreOrbitalParameters",
        ],
    },
    {
        "name": "hotAtmosphere",
        "description": "Hot atmosphere",
        "members": [
            "hotHaloOutflowReincorporation",
            re.compile(r"^hotHaloOutflowReincorporation.*"),
            "hotHaloTemperatureProfile",
            "hotHaloTrackStrippedGas",
            "hotHaloOutflowReturnRate",
            "hotHaloOutflowToColdMode",
            "hotHaloAngularMomentumAlwaysGrows",
            "hotHaloAngularMomentumLossFraction",
            "hotHaloColdModeCoredIsothermalCoreRadii",
            "hotHaloCoreRadiusOverVirialRadius",
            "hotHaloExcessHeatDrivesOutflow",
            "hotHaloMassDistribution",
            "hotHaloNodeMergerLimitBaryonFraction",
        ],
    },
    {
        "name": "coldAtmosphere",
        "description": "Cold atmosphere",
        "members": [
            "coldModeIsothermalCoreRadiusOverVirialRadius",
            "coldModeMassDistribution",
        ],
    },
    {
        "name": "accretionOntoHalos",
        "description": "Accretion onto halos",
        "members": [
            "accretionHaloTotal",
        ],
    },
    {
        "name": "accretionFromIGM",
        "description": "Accretion from the IGM",
        "members": [
            "accretionHalo",
        ],
    },
    {
        "name": "hotHaloRamPressure",
        "description": "Ram pressure stripping of hot atmosphere",
        "members": [
            "starveSatellites",
            "hotHaloRamPressureForce",
            "hotHaloRamPressureStripping",
            "hotHaloRamPressureStrippingTimescale",
            "hotHaloOutflowStrippingEfficiency",
            "ramPressureStrippingFormFactor",
        ],
    },
    {
        "name": "coolingAndInfall",
        "description": "Cooling and infall",
        "members": [
            "coolingRate",
            re.compile(r"^coolingRate.*"),
            "coolingFunction",
            "coolingMeanAngularMomentumFrom",
            "coolingRadius",
            "infallRadius",
            "coolingRotationVelocityFrom",
            "coolingSpecificAngularMomentum",
            "coolingTime",
            re.compile(r"^coolingTime.*"),
            "coldModeInfallRate",
            re.compile(r"^coldModeInfallRate.*"),
            "zeroCoolingRateAboveVelocity",
            "chemicalState",
        ],
    },
    {
        "name": "galacticStructure",
        "description": "Galactic structure",
        "members": [
            "diskMassDistribution",
            "spheroidMassDistribution",
            "galacticStructureRadiusSolver",
            re.compile(r"^galacticStructureRadii.*"),
            "adiabaticContractionGnedinA",
            "adiabaticContractionGnedinOmega",
            "spheroidAngularMomentumAtScaleRadius",
        ],
    },
    {
        "name": "barInstability",
        "description": "Galactic disk bar instability",
        "members": [
            "barInstability",
            "stabilityThresholdGaseous",
            "stabilityThresholdStellar",
        ],
    },
    {
        "name": "diskRamPressure",
        "description": "Ram pressure stripping of galactic disk",
        "members": [
            "ramPressureStrippingMassLossRateDisks",
        ],
    },
    {
        "name": "spheroidRamPressure",
        "description": "Ram pressure stripping of galactic spheroid",
        "members": [
            "ramPressureStrippingMassLossRateSpheroids",
        ],
    },
    {
        "name": "diskTidal",
        "description": "Tidal stripping of galactic disk",
        "members": [
            "tidalStrippingMassLossRateDisks",
        ],
    },
    {
        "name": "spheroidTidal",
        "description": "Tidal stripping of galactic spheroid",
        "members": [
            "tidalStrippingMassLossRateSpheroids",
        ],
    },
    {
        "name": "diskStarFormation",
        "description": "Star formation in disks",
        "members": [
            "starFormationTimescaleDisks",
            re.compile(r"^starFormationTimescaleDisks.*"),
            re.compile(r"^diskStarFormation.*"),
            re.compile(r"^diskVerySimpleSurfaceDensity.*"),
            "starFormationDiskMinimumTimescale",
            "starFormationFrequencyKMT09",
            "molecularComplexClumpingFactorKMT09",
            "molecularFractionFastKMT09",
            "starFormationRateSurfaceDensityDisks",
            "starFormationTimescaleIntegratedSurfaceDensityTolerance",
        ],
    },
    {
        "name": "spheroidStarFormation",
        "description": "Star formation in spheroids",
        "members": [
            "starFormationTimescaleSpheroids",
            re.compile(r"^starFormationTimescaleSpheroids.*"),
            re.compile(r"^spheroidStarFormation.*"),
            "starFormationSpheroidEfficiency",
            "starFormationSpheroidMinimumTimescale",
            "starFormationSpheroidVelocityExponent",
        ],
    },
    {
        "name": "diskStellarFeedback",
        "description": "Stellar feedback in disks",
        "members": [
            "starFormationFeedbackDisks",
            re.compile(r"^diskOutflow.*"),
        ],
    },
    {
        "name": "spheroidStellarFeedback",
        "description": "Stellar feedback in spheroids",
        "members": [
            "starFormationFeedbackSpheroids",
            re.compile(r"^spheroidOutflow.*"),
        ],
    },
    {
        "name": "initialMassFunction",
        "description": "Stellar initial mass function",
        "members": [
            "imfSelection",
            "imfSelectionFixed",
            re.compile(r"^imf.*RecycledInstantaneous$"),
            re.compile(r"^imf.*YieldInstantaneous$"),
        ],
    },
    {
        "name": "stellarEvolution",
        "description": "Stellar evolution",
        "members": [
            "stellarPopulationProperties",
            "stellarPopulationSpectra",
            "stellarPopulationLuminosityIntegrationToleranceRelative",
        ],
    },
    {
        "name": "smbhAccretionDisks",
        "description": "Super-massive black hole accretion disks",
        "members": [
            "accretionDisks",
            "accretionDiskSwitchedScaleAdafRadiativeEfficiency",
            "accretionRateThinDiskMaximum",
            "accretionRateThinDiskMinimum",
            "accretionRateTransitionWidth",
            "adafAdiabaticIndex",
            "adafEnergyOption",
            "adafRadiativeEfficiency",
            "adafRadiativeEfficiencyType",
            "adafViscosityOption",
        ],
    },
    {
        "name": "superMassiveBlackHoles",
        "description": "Super massive black holes",
        "members": [
            "blackHoleBinaryMergers",
            "blackHoleHeatsHotHalo",
            "blackHoleSeedMass",
            "blackHoleWindEfficiency",
            "blackHoleWindEfficiencyScalesWithRadiativeEfficiency",
            "blackHoleRadioModeFeedbackEfficiency",
            "bondiHoyleAccretionEnhancementHotHalo",
            "bondiHoyleAccretionEnhancementSpheroid",
            "bondiHoyleAccretionHotModeOnly",
            "bondiHoyleAccretionTemperatureSpheroid",
            "spheroidEnergeticOutflowMassRate",
        ],
    },
    {
        "name": "galaxyMerging",
        "description": "Galaxy merging",
        "members": [
            "satelliteMergingMassMovements",
            "satelliteMergingRemnantSize",
            "majorMergerMassRatio",
            "mergerRemnantSizeOrbitalEnergy",
            "minorMergerGasMovesTo",
        ],
    },
    {
        "name": "solvers",
        "description": "Solvers and time-stepping",
        "members": [
            "odeAlgorithm",
            "odeToleranceAbsolute",
            "odeToleranceRelative",
            "diskMassToleranceAbsolute",
            "spheroidMassToleranceAbsolute",
            "diskVerySimpleMassScaleAbsolute",
            "hotHaloVerySimpleDelayedMassScaleRelative",
            re.compile(r"^timestep.*"),
            re.compile(r".*AnalyticSolver.*"),
        ],
    },
    {
        "name": "output",
        "description": "Output",
        "members": [
            "outputFileName",
            "hdf5CompressionLevel",
            "mergerTreeOutput",
            "mergerTreeOutputReferences",
            "outputRedshifts",
            "outputSatelliteStatus",
            "outputColdModeInfallRate",
        ],
    },
    {
        "name": "luminosities",
        "description": "Stellar luminosities",
        "members": [
            "luminosityFilter",
            "luminosityRedshift",
            "luminosityType",
        ],
    },
    {
        "name": "numerics",
        "description": "Numerics",
        "members": [
            "randomSeed",
        ],
    },
    {
        "name": "analysis",
        "description": "Analysis",
        "members": [
            "mergerTreeAnalyses",
            re.compile(r"^analysisMassFunction.*"),
            re.compile(r"^analysisMassFunctions.*"),
        ],
    },
    {
        "name": "haloModel",
        "description": "Halo model [clustering]",
        "members": [
            "haloModelPowerSpectrumModifier",
        ],
    },
    {
        "name": "stellarMassSystematics",
        "description": "Stellar mass random and systematic error models [constraints]",
        "members": [
            re.compile(r".*StellarMassFunctionZ[\d.]+MassSystematic\d+$"),
            re.compile(r"^stellarMassSystematics\d+$"),
            re.compile(r".*StellarMassFunctionZ[\d.]+MassRandom\d+$"),
            re.compile(r".*StellarMassFunctionZ[\d.]+MassRandomMinimum$"),
            re.compile(r".*StellarMassFunctionZ[\d.]+MassRandomMaximum$"),
            re.compile(r"^sdssStellarMassFunctionHighMass.*"),
            re.compile(r"^gamaStellarMassFunctionZ0\.03"),
        ],
    },
    {
        "name": "hiMassSystematics",
        "description": "HI mass random and systematic error models [constraints]",
        "members": [
            re.compile(r"alfalfaHiMassFunctionZ0\.00.*$"),
        ],
    },
    {
        "name": "sizeSystematics",
        "description": "Galaxy size random and systematic error models [constraints]",
        "members": [
            re.compile(r"^.*SizeFunctionZ[\d.]+MassSystematic\d+$"),
            re.compile(r"^.*SizeFunctionZ[\d.]+MassRandom\d+$"),
            re.compile(r"^.*SizeFunctionZ[\d.]+MassRandomMinimum$"),
            re.compile(r"^.*SizeFunctionZ[\d.]+MassRandomMaximum$"),
            re.compile(r"^.*SizeFunctionZ[\d.]+RadiusSystematic\d+$"),
        ],
    },
    {
        "name": "clusteringSystematics",
        "description": "Clustering random and systematic error models [constraints]",
        "members": [
            re.compile(r"^.*ClusteringZ[\d.]+MassSystematic\d+$"),
            re.compile(r"^.*ClusteringZ[\d.]+MassRandom\d+$"),
            re.compile(r"^.*ClusteringZ[\d.]+MassRandomMinimum$"),
        ],
    },
]


def _is_comment(node):
    """Return True if the ElementTree node is an XML comment."""
    return callable(node.tag)


def _node_has_value(node):
    """Return True if the element has a 'value' attribute or a 'value' child."""
    if _is_comment(node):
        return False
    return "value" in node.attrib or node.find("value") is not None


def _member_rank(group, tag):
    """Return the index of the first matching member pattern, or None."""
    for i, member in enumerate(group["members"]):
        if isinstance(member, str):
            if tag == member:
                return i
        else:
            if member.search(tag):
                return i
    return None


def _process_parameters(parameters):
    """Reorganise one <parameters> element in-place."""
    new_groups = {}

    # Snapshot the current value-bearing child elements.
    param_nodes = [c for c in list(parameters) if _node_has_value(c)]

    # Assign each parameter to the first group it matches and remove it from
    # the parent temporarily.
    for param in param_nodes:
        for group in GROUPS:
            rank = _member_rank(group, param.tag)
            if rank is not None:
                parameters.remove(param)
                new_groups.setdefault(group["name"], []).append(param)
                break

    # If any ungrouped value parameters remain, add a "Miscellaneous" header
    # at the very start (it will be pushed to the end as groups are prepended).
    if any(_node_has_value(c) for c in parameters):
        parameters.insert(0, ET.Comment(" Miscellaneous parameters "))

    # Insert groups in reverse order so that the first group ends up at the top.
    for group in reversed(GROUPS):
        if group["name"] not in new_groups:
            continue

        # Sort by position in the member list; alphabetical order breaks ties.
        sorted_params = sorted(
            new_groups[group["name"]],
            key=lambda n: (_member_rank(group, n.tag), n.tag),
        )

        # Insert in reverse-sorted order at position 0 so the final order is
        # ascending (mirrors the Perl insertBefore-at-firstChild logic).
        for param in reversed(sorted_params):
            parameters.insert(0, param)

        parameters.insert(0, ET.Comment(" " + group["description"] + " "))


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: parametersFormat.py <inputFile> <outputFile>",
            file=sys.stderr,
        )
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Parse while preserving comment nodes.
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    tree = ET.parse(input_file, parser)
    root = tree.getroot()

    if root.tag == "parameters":
        parameter_sets = [root]
    else:
        parameter_sets = list(root.iter("parameters"))

    if not parameter_sets:
        sys.exit("parametersFormat.py: cannot find parameters")

    for parameters in parameter_sets:
        _process_parameters(parameters)

    # Apply 2-space pretty-print indentation.
    ET.indent(tree, space="  ")

    # Insert a blank line before every section comment that is a direct child
    # of a <parameters> element (matches the Perl "\n  " text-node insertion).
    for parameters in parameter_sets:
        children = list(parameters)
        for i, child in enumerate(children):
            if _is_comment(child):
                if i == 0:
                    parameters.text = "\n\n  "
                else:
                    children[i - 1].tail = "\n\n  "

    tree.write(output_file, xml_declaration=True, encoding="UTF-8")


if __name__ == "__main__":
    main()
