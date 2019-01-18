#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Scalar::Util 'reftype';
use XML::LibXML qw(:libxml);
use XML::LibXML::PrettyPrint;
use Data::Dumper;
use Galacticus::Options;
use List::ExtraUtils;

# Update a Galacticus parameter file from one version to a later version.
# Andrew Benson (13-September-2014)

# Get arguments.
die("Usage: parametersMigrate.pl <inputFile> <outputFile> [options]")
    unless ( scalar(@ARGV) >= 2 );
my $inputFileName  = $ARGV[0];
my $outputFileName = $ARGV[1];
my %options =
    (
     inputVersion        => "0.9.3",
     outputVersion       => "0.9.4",
     validate            => "yes"  ,
     prettyify           => "no"   ,
     inputFormatVersion  => 1      ,
     outputFormatVersion => 2
    );
# Parse options.
my %optionsDefined = &Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Define translations.
my @translations =
    (
      {
	 inputVersion  => "0.9.0",
	 outputVersion => "0.9.1",
 	 names         =>
	 {
 	     Omega_0 => "Omega_Matter",
	 },
	 values        =>
	 {
	     adafEnergyOption                      =>
	     {
		 "pure ADAF"                => "pureADAF"               ,
	     },
	     coolingFunctionMethods                =>
	     {
		 "atomic_CIE_Cloudy"        => "atomicCIECloudy"        ,
		 "CIE_from_file"            => "cieFromFile"            ,
	     },
	     coolingSpecificAngularMomentumMethod  => 
	     {
		 "constant rotation"        => "constantRotation"       ,
	     },
	     coolingTimeAvailableMethod            => 
	     {
		 "White-Frenk"              => "White-Frenk1991"        ,
	     },
	     coolingRateMethod                     => 
	     {
		 "White + Frenk"            => "White-Frenk1991"        ,
	     },
	     cosmologyMethod                       =>
	     {
		 "matter + lambda"          => "matter-lambda"          ,
	     },
	     criticalOverdensityMethod             =>
	     {
		 "spherical top hat"        => "sphericalTopHat"        ,
	     },
	     darkMatterConcentrationMethod         =>
	     {
		 "Gao 2008"                 => "Gao2008"                ,
	     },
	     hotHaloCoolingFromNode                =>
	     {
		 "current node"             => "currentNode"            ,
		 "formation node"           => "formationNode"
	     },
	     hotHaloDensityMethod                  =>
	     {
		 "cored isothermal"         => "coredIsothermal"
	     },
	     infallRadiusMethod                    =>
	     {
		 "cooling and freefall"     => "coolingAndFreefall"     ,
	     },
	     nodeMergersMethod                     => 
	     {
		 "single level hierarchy"   => "singleLevelHierarchy"   ,
	     },
	     powerSpectrumMethod                   => 
	     {
		 "power law"                => "powerLaw"               ,
	     },
	     satelliteMergingMethod                => 
	     {
		 "Lacey-Cole + Tormen"      => "Lacey-Cole+Tormen"      ,
	     },
	     starFormationHistoriesMethod          => 
	     {
		 "metallicity split"        => "metallicitySplit"       ,
	     },
	     starFormationHistoriesMethod          => 
	     {
		 "metallicity split"        => "metallicitySplit"       ,
	     },
	     starFormationFeedbackDisksMethod      => 
	     {
		 "power law"                => "powerLaw"               ,
	     },
	     starFormationFeedbackSpheroidsMethod  => 
	     {
		 "power law"                => "powerLaw"               ,
	     },
	     starFormationTimescaleDisksMethod     => 
	     {
		 "dynamical time"           => "dynamicalTime"          ,
	     },
	     starFormationTimescaleSpheroidsMethod => 
	     {
		 "dynamical time"           => "dynamicalTime"          ,
	     },
	     transferFunctionMethod                => 
	     {
		 "Eisenstein + Hu"          => "Eisenstein-Hu1999"      ,
	     },
	     treeBranchingMethod                   => 
	     {
		 "modified Press-Schechter" => "modifiedPress-Schechter",
	     },
	     virialDensityContrastMethod           =>
	     {
		 "spherical top hat"        => "sphericalTopHat"        ,
	     },
	 },
      },
      {
	 inputVersion  => "0.9.1",
	 outputVersion => "0.9.2",
 	 names         =>
	 {
 	     treeNodeMethodSatelliteOrbit                => "treeNodeMethodSatellite",
	 },
	 values        =>
	 {
 	     treeNodeMethodSatellite              =>
	     {
		 "simple"                    => "standard"
	     },
	     treeNodeMethodSpheroid               =>
	     {
		 "Hernquist"                 => {
		                                 value => "standard",
						 new   => [
						           {
							    name  => "spheroidMassDistribution",
							    value => "hernquist"
						           }
                                                          ]
		                                },
		 "sersic"                    => {
		                                 value => "standard",
						 new   => [
						           {
							    name  => "spheroidMassDistribution",
							    value => "sersic"
						           }
                                                          ]
		                                }
	     },
	 }
     },
     {
	 inputVersion  => "0.9.2",
	 outputVersion => "0.9.3",
	 names         =>
	 {
	     accretionHalosMethod                        => "accretionHaloMethod"                             ,
	     cosmologyMethod                             => "cosmologyFunctionsMethod"                        ,
	     darkMatterConcentrationMethod               => "darkMatterProfileConcentrationMethod"            ,
	     hotHaloCoredIsothermalCoreRadiiMethod       => "hotHaloColdModeCoredIsothermalCoreRadiiMethod"   ,
	     hotHaloDensityMethod                        => "hotHaloMassDistributionMethod"                   ,
	     ionizationStateFile                         => "chemicalStateFile"                               ,
	     isothermalCoreRadiusOverScaleRadius         => "hotHaloCoreRadiusOverScaleRadius"                ,
	     isothermalCoreRadiusOverVirialRadius        => "hotHaloCoreRadiusOverVirialRadius"               ,
	     isothermalCoreRadiusOverVirialRadiusMaximum => "hotHaloCoreRadiusOverVirialRadiusMaximum"        ,
	     mergerTreeBuildCole2000MassResolution       => "mergerTreeHaloMassResolution"                    ,
	     nfw96ConcentrationC                         => "nfw1996ConcentrationC"                           ,
	     satelliteMergingMethod                      => "satelliteMergingTimescalesMethod"                ,
	     treeNodeMethodFormationTimes                => "treeNodeMethodFormationTime"                     ,
	     luminosityFilterAbsoluteMagnitudeThresholds => "luminosityFilterAbsoluteMagnitudeThresholdMaxima"
	 },
	 values        =>
         {
	     cosmologyFunctionsMethod             =>
	     {
		 "matter-lambda"             => "matterLambda"                 ,
		 "matter-darkEnergy"         => "matterDarkEnergy"             ,
	     },
	     darkMatterProfileConcentrationMethod =>
	     {
		 "Gao2008"                   => "gao2008"                      ,
		 "Munoz-Cuartas2011"         => "munozCuartas2011"             ,
		 "Prada2011"                 => "prada2011"                    ,
		 "Zhao2009"                  => "zhao2009"
	     },
	     hotHaloMassDistributionMethod        =>
	     {
		 "coredIsothermal"           => "betaProfile"
	     },
	     satelliteMergingTimescalesMethod     => 
	     {
		 "BoylanKolchin2008"         => "boylanKolchin2008"            ,
		 "Jiang2008"                 => "jiang2008"                    ,
		 "Lacey-Cole"                => "laceyCole1993"                ,
		 "Lacey-Cole+Tormen"         => "laceyCole1993Tormen"          ,
		 "Wetzel-White2010"          => "wetzelWhite2010"
	     },
	     virialDensityContrastMethod          =>
	     {
		 "Bryan-Norman1998"          => "bryanNorman1998"              ,
		 "sphericalTopHatDarkEnergy" => "sphericalCollapseMatterDE"    ,
		 "sphericalTopHat"           => "sphericalCollapseMatterLambda",
		 "Kitayama-Suto1996"         => "kitayamaSuto1996"
	     },
	     virialOrbitsMethod                   =>
	     {
		 "Benson2005"                => "benson2005"                   ,
		 "Wetzel2010"                => "wetzel2010"
 
	     },
	 }
     },
     {
	 inputVersion  => "0.9.3",
	 outputVersion => "0.9.4",
	 names         =>
	 {
	     "hotHaloTemperatureMethod"                                              => "hotHaloTemperatureProfileMethod"                                                         ,
	     "mergerTreeBuildMethod"                                                 => "mergerTreeBuilderMethod"                                                                 ,
	     "darkMatterShapeMethod"                                                 => "darkMatterProfileShapeMethod"                                                            ,
	     "H_0"                                                                   => "cosmologyParametersMethod--HubbleConstant"                                               ,
	     "Omega_Matter"                                                          => "cosmologyParametersMethod--OmegaMatter"                                                  ,
	     "Omega_DE"                                                              => "cosmologyParametersMethod--OmegaDarkEnergy"                                              ,
	     "Omega_b"                                                               => "cosmologyParametersMethod--OmegaBaryon"                                                  ,
	     "T_CMB"                                                                 => "cosmologyParametersMethod--temperatureCMB"                                               ,
	     "effectiveNumberNeutrinos"                                              => "transferFunctionMethod--neutrinoNumberEffective"                                         ,
	     "summedNeutrinoMasses"                                                  => "transferFunctionMethod--neutrinoMassSummed"                                              ,
	     "transferFunctionWDMFreeStreamingLength"                                => "transferFunctionMethod--freeStreamingLength"                                             ,
	     "transferFunctionWdmCutOffScale"                                        => "transferFunctionMethod--scaleCutOff"                                                     ,
	     "transferFunctionWdmEpsilon"                                            => "transferFunctionMethod--epsilon"                                                         ,
	     "transferFunctionWdmEta"                                                => "transferFunctionMethod--eta"                                                             ,
	     "transferFunctionWdmNu"                                                 => "transferFunctionMethod--nu"                                                              ,
	     "stellarPopulationSpectraFileForceZeroMetallicity"                      => "stellarPopulationSpectraMethod--forceZeroMetallicity"                                    ,
	     "stellarPopulationSpectraForChabrierIMF"                                => "stellarPopulationSpectraMethod--fileNameForChabrierIMF"                                  ,
	     "stellarPopulationSpectraForBaugh2005TopHeavyIMF"                       => "stellarPopulationSpectraMethod--fileNameForBaugh2005TopHeavyIMF"                         ,
	     "stellarPopulationSpectraForKroupaIMF"                                  => "stellarPopulationSpectraMethod--fileNameForKroupaIMF"                                    ,
	     "stellarPopulationSpectraForMiller-ScaloIMF"                            => "stellarPopulationSpectraMethod--fileNameForMiller-ScaloIMF"                              ,
	     "stellarPopulationSpectraForSalpeterIMF"                                => "stellarPopulationSpectraMethod--fileNameForSalpeterIMF"                                  ,
	     "stellarPopulationSpectraForScaloIMF"                                   => "stellarPopulationSpectraMethod--fileNameForScaloIMF"                                     ,
	     "stellarPopulationSpectraForKennicuttIMF"                               => "stellarPopulationSpectraMethod--fileNameForKennicuttIMF"                                 ,
	     "accretionDiskSpectraFileName"                                          => "accretionDiskSpectraMethod--fileName"                                                    ,
	     "chemicalStateFile"                                                     => "chemicalStateMethod--fileName"                                                           ,
	     "coolingFunctionMethods"                                                => "coolingFunctionMethod"                                                                   ,
	     "powerSpectrumMethod"                                                   => "powerSpectrumPrimordialMethod"                                                           ,
	     "powerSpectrumIndex"                                                    => "powerSpectrumPrimordialMethod--index"                                                    ,
	     "powerSpectrumRunning"                                                  => "powerSpectrumPrimordialMethod--running"                                                  ,
	     "powerSpectrumReferenceWavenumber"                                      => "powerSpectrumPrimordialMethod--wavenumberReference"                                      ,
	     "mergerTreeBuildCole2000AccretionLimit"                                 => "mergerTreeBuilderMethod--accretionLimit"                                                 ,
	     "mergerTreeBuildCole2000MergeProbability"                               => "mergerTreeBuilderMethod--mergeProbability"                                               ,
	     "mergerTreeBuildCole2000HighestRedshift"                                => "mergerTreeBuilderMethod--redshiftMaximum"                                                ,
	     "mergerTreeBuildCole2000FixedRandomSeeds"                               => "mergerTreeBuilderMethod--randomSeedsFixed"                                               ,
	     "mergerTreeRegridTimes"                                                 => "mergerTreeOperatorMethod.regridTimes."                                                   ,
	     "mergerTreeRegridDumpTrees"                                             => "mergerTreeOperatorMethod.regridTimes.--dumpTrees"                                        ,
	     "mergerTreeRegridCount"                                                 => "mergerTreeOperatorMethod.regridTimes.--regridCount"                                      ,
	     "mergerTreeRegridStartExpansionFactor"                                  => "mergerTreeOperatorMethod.regridTimes.--expansionFactorStart"                             ,
	     "mergerTreeRegridEndExpansionFactor"                                    => "mergerTreeOperatorMethod.regridTimes.--expansionFactorEnd"                               ,
	     "mergerTreeRegridSpacing"                                               => "mergerTreeOperatorMethod.regridTimes.--snapshotSpacing"                                  ,
	     "mergerTreePruneBranches"                                               => "mergerTreeOperatorMethod.pruneByMass."                                                   ,
	     "mergerTreePruningMassThreshold"                                        => "mergerTreeOperatorMethod.pruneByMass.--massThreshold"                                    ,
	     "mergerTreePruneHierarchyAtDepth"                                       => "mergerTreeOperatorMethod.pruneHierarchy.--hierarchyDepth"                                ,
	     "mergerTreePruneNonEssential"                                           => "mergerTreeOperatorMethod.pruneNonEssential"                                              ,
	     "mergerTreePruningNonEssentialID"                                       => "mergerTreeOperatorMethod.pruneNonEssential.--essentialNodeID"                            ,
	     "mergerTreePruningNonEssentialTime"                                     => "mergerTreeOperatorMethod.pruneNonEssential.--essentialNodeTime"                          ,
	     "mergerTreeComputeConditionalMassFunction"                              => "mergerTreeOperatorMethod.conditionalMF."                                                 ,
	     "mergerTreeComputeConditionalMassFunctionParentMassCount"               => "mergerTreeOperatorMethod.conditionalMF.--parentMassCount"                                ,
	     "mergerTreeComputeConditionalMassFunctionParentMassMinimum"             => "mergerTreeOperatorMethod.conditionalMF.--parentMassMinimum"                              ,
	     "mergerTreeComputeConditionalMassFunctionParentMassMaximum"             => "mergerTreeOperatorMethod.conditionalMF.--parentMassMaximum"                              ,
	     "mergerTreeComputeConditionalMassFunctionMassRatioCount"                => "mergerTreeOperatorMethod.conditionalMF.--massRatioCount"                                 ,
	     "mergerTreeComputeConditionalMassFunctionMassRatioMinimum"              => "mergerTreeOperatorMethod.conditionalMF.--massRatioMinimum"                               ,
	     "mergerTreeComputeConditionalMassFunctionMassRatioMaximum"              => "mergerTreeOperatorMethod.conditionalMF.--massRatioMaximum"                               ,
	     "mergerTreeComputeConditionalMassFunctionParentRedshifts"               => "mergerTreeOperatorMethod.conditionalMF.--parentRedshifts"                                ,
	     "mergerTreeComputeConditionalMassFunctionProgenitorRedshifts"           => "mergerTreeOperatorMethod.conditionalMF.--progenitorRedshifts"                            ,
	     "mergerTreeComputeConditionalMassFunctionPrimaryProgenitorDepth"        => "mergerTreeOperatorMethod.conditionalMF.--primaryProgenitorDepth"                         ,
	     "mergerTreeConditionalMassFunctionFormationRateTimeFraction"            => "mergerTreeOperatorMethod.conditionalMF.--formationRateTimeFraction"                      ,
	     "mergerTreeExportFileName"                                              => "mergerTreeOperatorMethod.export.--outputFileName"                                        ,
	     "mergerTreeExportOutputFormat"                                          => "mergerTreeOperatorMethod.export.--exportFormat"                                          ,
	     "warmDarkMatterCriticalOverdensityGX"                                   => "criticalDensityMethod.barkana2001WDM.--gX"                                               ,
	     "warmDarkMatterCriticalOverdensityMX"                                   => "criticalDensityMethod.barkana2001WDM.--mX"                                               ,
	     "warmDarkMatterCriticalOverdensityUseFittingFunction"                   => "criticalDensityMethod.barkana2001WDM.--useFittingFunction"                               ,
	     "sigma_8"                                                               => "cosmologicalMassVarianceMethod.filteredPower.--sigma_8"                                  ,
	     "mergerTreeExportOutputFormat"                                          => "mergerTreeOperatorMethod.write.--exportFormat"                                           ,
	     "darkMatterHaloConcentrationCorrea2015AScaling"                         => "darkMatterProfileConcentrationMethod.correa2015.--A"                                     ,
	     "darkMatterProfileConcentrationDiemerKravtsov2014Kappa"                 => "darkMatterProfileConcentrationMethod.diemerKravtsov2014.--kappa"                         ,
	     "darkMatterProfileConcentrationDiemerKravtsov2014Phi0"                  => "darkMatterProfileConcentrationMethod.diemerKravtsov2014.--phi0"                          ,
	     "darkMatterProfileConcentrationDiemerKravtsov2014Phi1"                  => "darkMatterProfileConcentrationMethod.diemerKravtsov2014.--phi1"                          ,
	     "darkMatterProfileConcentrationDiemerKravtsov2014Eta0"                  => "darkMatterProfileConcentrationMethod.diemerKravtsov2014.--eta0"                          ,
	     "darkMatterProfileConcentrationDiemerKravtsov2014Eta1"                  => "darkMatterProfileConcentrationMethod.diemerKravtsov2014.--eta1"                          ,
	     "darkMatterProfileConcentrationDiemerKravtsov2014Alpha"                 => "darkMatterProfileConcentrationMethod.diemerKravtsov2014.--alpha"                         ,
	     "darkMatterProfileConcentrationDiemerKravtsov2014Beta"                  => "darkMatterProfileConcentrationMethod.diemerKravtsov2014.--beta"                          ,
	     "duttonMaccio2014FitType"                                               => "darkMatterProfileConcentrationMethod.duttonMaccio2014.--fitType"                         ,
	     "klypin2015ConcentrationSample"                                         => "darkMatterProfileConcentrationMethod.klypin2015.--sample"                                ,
	     "nfw1996ConcentrationF"                                                 => "darkMatterProfileConcentrationMethod.nfw1996.--f"                                        ,
	     "nfw1996ConcentrationC"                                                 => "darkMatterProfileConcentrationMethod.nfw1996.--C"                                        ,
	     "prada2011ConcentrationA"                                               => "darkMatterProfileConcentrationMethod.prada2011.--A"                                      ,
	     "prada2011ConcentrationB"                                               => "darkMatterProfileConcentrationMethod.prada2011.--B"                                      ,
	     "prada2011ConcentrationC"                                               => "darkMatterProfileConcentrationMethod.prada2011.--C"                                      ,
	     "prada2011ConcentrationD"                                               => "darkMatterProfileConcentrationMethod.prada2011.--D"                                      ,
	     "prada2011ConcentrationC0"                                              => "darkMatterProfileConcentrationMethod.prada2011.--C0"                                     ,
	     "prada2011ConcentrationC1"                                              => "darkMatterProfileConcentrationMethod.prada2011.--C1"                                     ,
	     "prada2011ConcentrationX0"                                              => "darkMatterProfileConcentrationMethod.prada2011.--X0"                                     ,
	     "prada2011ConcentrationX1"                                              => "darkMatterProfileConcentrationMethod.prada2011.--X1"                                     ,
	     "prada2011ConcentrationInverseSigma0"                                   => "darkMatterProfileConcentrationMethod.prada2011.--inverseSigma0"                          ,
	     "prada2011ConcentrationInverseSigma1"                                   => "darkMatterProfileConcentrationMethod.prada2011.--inverseSigma1"                          ,
	     "prada2011ConcentrationAlpha"                                           => "darkMatterProfileConcentrationMethod.prada2011.--alpha"                                  ,
	     "prada2011ConcentrationBeta"                                            => "darkMatterProfileConcentrationMethod.prada2011.--beta"                                   ,
	     "darkMatterProfileConcentrationCDMMethod"                               => "darkMatterProfileConcentrationMethod.wdm.--darkMatterProfileConcentrationMethod"         ,
	     "mergerTreeBuildMassResolutionFixed"                                    => "mergerTreeMassResolutionMethod.fixed.--massResolution"                                   ,
	     "mergerTreeBuildMassResolutionScaledMinimum"                            => "mergerTreeMassResolutionMethod.scaled.--massResolutionMinimum"                           ,
	     "mergerTreeBuildMassResolutionScaledFraction"                           => "mergerTreeMassResolutionMethod.scaled.--massResolutionFractional"                        ,
	     "powerSpectrumWindowFunctionSharpKSpaceNormalization"                   => "powerSpectrumWindowFunctionMethod.sharpKSpace.--normalization"                           ,
	     "powerSpectrumWindowFunctionSharpKSpaceTopHatRadiiRatio"                => "powerSpectrumWindowFunctionMethod.topHatSharpKHybrid.--radiiRatio"                       ,
	     "haloMassFunctionSimpleSystematicAlpha"                                 => "haloMassFunctionMethod.simpleSystematic.--alpha"                                         ,
	     "haloMassFunctionSimpleSystematicBeta"                                  => "haloMassFunctionMethod.simpleSystematic.--beta"                                          ,
	     "mergerTreeBuildMassResolutionFixed"                                    => "mergerTreeMassResolutionMethod.fixed.--massResolution"                                   ,
	     "spinDistributionBett2007Lambda0"	                                     => "haloSpinDistributionMethod.bett2007.--lambda0"                                           ,
	     "spinDistributionBett2007Alpha"	                                     => "haloSpinDistributionMethod.bett2007.--alpha"                                             ,
	     "lognormalSpinDistributionMedian"                                       => "haloSpinDistributionMethod.logNormal.--median"                                           ,
	     "lognormalSpinDistributionSigma"	                                     => "haloSpinDistributionMethod.logNormal.--sigma"                                            ,
	     "deltaFunctionSpinDistributionSpin"	                             => "haloSpinDistributionMethod.deltaFunction.--spin"                                         ,
	     "stellarMassFilterThreshold"                                            => "mergerTreeOutput--galacticFilterMethod.stellarMass.--massThreshold"                      ,
	     "mergerTreeOutput"                                                      => "mergerTreeOutput--outputMergerTrees"                                                     ,
	     "mergerTreeOutputReferences"                                            => "mergerTreeOutput--outputReferences"                                                      ,
	     "mergerTreeAnalyses"                                                    => "mergerTreeOutput--analyses"                                                              ,
	     "mergerTreeConstructMethod"                                             => "mergerTreeConstructorMethod"                                                             ,
	     "mergerTreesBuildFixedThreadAssignment"                                 => "evolveForestsWorkShareMethod.cyclic."                                                    ,
	     "mergerTreeBuildTreeMassesFile"                                         => "mergerTreeBuildMasses.readHDF5.--fileName"                                               ,
	     "mergerTreeBuildHaloMassMinimum"                                        => "mergerTreeBuildMassesMethod.sampledDistributionUniform.--massTreeMinimum"                ,
	     "mergerTreeBuildHaloMassMaximum"                                        => "mergerTreeBuildMassesMethod.sampledDistributionUniform.--massTreeMaximum"                ,
	     "mergerTreeBuildTreesPerDecade"                                         => "mergerTreeBuildMassesMethod.sampledDistributionUniform.--treesPerDecade"                 ,
	     "haloMassSampleModifierFixedMassHalosMass"                              => "mergerTreeBuildMassesMethod.fixedMass.--massTree"                                        , 
	     "haloMassSampleModifierFixedMassHalosCount"                             => "mergerTreeBuildMassesMethod.fixedMass.--treeCount"                                       ,
	     "haloMassSampleModifierFixedMassHalosRadius"                            => "mergerTreeBuildMassesMethod.fixedMass.--radiusTree"                                      ,
	     "mergerTreeHaloMass"                                                    => "mergerTreeConstructorMethod.smoothAccretion.--massHalo"                                  ,
	     "mergerTreeHaloMassResolution"                                          => "mergerTreeConstructorMethod.smoothAccretion.--massHaloResolution"                        ,
	     "mergerTreeHaloMassDeclineFactor"                                       => "mergerTreeConstructorMethod.smoothAccretion.--massHaloDeclineFactor"                     ,
	     "mergerTreeBaseRedshift"                                                => "mergerTreeConstructorMethod.smoothAccretion.--redshiftBase"                              ,
	     "mergerTreeConstructFullySpecifiedFileName"                             => "mergerTreeConstructorMethod.fullySpecified.--fileName"                                   ,
	     "mergerTreeStateStoreFile"                                              => "mergerTreeConstructorMethod.stateRestore.--fileName"                                     ,
	     "mergerTreeReadFileName"                                                => "mergerTreeConstructorMethod.read.--fileNames"                                            ,
	     "mergerTreeReadForestSizeMaximum"                                       => "mergerTreeConstructorMethod.read.--forestSizeMaximum"                                    ,
	     "mergerTreeReadPresetMergerTimes"                                       => "mergerTreeConstructorMethod.read.--presetMergerTimes"                                    ,
	     "mergerTreeReadPresetMergerNodes"                                       => "mergerTreeConstructorMethod.read.--presetMergerNodes"                                    ,
	     "mergerTreeReadPresetSubhaloMasses"                                     => "mergerTreeConstructorMethod.read.--presetSubhaloMasses"                                  ,
	     "mergerTreeReadSubhaloAngularMomentaMethod"                             => "mergerTreeConstructorMethod.read.--subhaloAngularMomentaMethod"                          ,
	     "mergerTreeReadPresetSubhaloIndices"                                    => "mergerTreeConstructorMethod.read.--presetSubhaloIndices"                                 ,
	     "mergerTreeReadPresetPositions"                                         => "mergerTreeConstructorMethod.read.--presetPositions"                                      ,
	     "mergerTreeReadPresetScaleRadii"                                        => "mergerTreeConstructorMethod.read.--presetScaleRadii"                                     ,
	     "mergerTreeReadPresetScaleRadiiFailureIsFatal"                          => "mergerTreeConstructorMethod.read.--scaleRadiiFailureIsFatal"                             ,
	     "mergerTreeReadPresetScaleRadiiConcentrationMinimum"                    => "mergerTreeConstructorMethod.read.--presetScaleRadiiMinimumMass"                          ,
	     "mergerTreeReadPresetScaleRadiiConcentrationMaximum"                    => "mergerTreeConstructorMethod.read.--presetScaleRadiiConcentrationMaximum"                 ,
	     "mergerTreeReadPresetScaleRadiiMinimumMass"                             => "mergerTreeConstructorMethod.read.--presetScaleRadiiMinimumMass"                          ,
	     "mergerTreeReadPresetUnphysicalSpins"                                   => "mergerTreeConstructorMethod.read.--presetUnphysicalSpins"                                ,
	     "mergerTreeReadPresetSpins"                                             => "mergerTreeConstructorMethod.read.--presetSpins"                                          ,
	     "mergerTreeReadPresetSpins3D"                                           => "mergerTreeConstructorMethod.read.--presetSpins3D"                                        ,
	     "mergerTreeReadPresetOrbits"                                            => "mergerTreeConstructorMethod.read.--presetOrbits"                                         ,
	     "mergerTreeReadPresetOrbitsSetAll"                                      => "mergerTreeConstructorMethod.read.--presetOrbitsSetAll"                                   ,
	     "mergerTreeReadPresetOrbitsAssertAllSet"                                => "mergerTreeConstructorMethod.read.--presetOrbitsAssertAllSet"                             ,
	     "mergerTreeReadPresetOrbitsBoundOnly"                                   => "mergerTreeConstructorMethod.read.--presetOrbitsBoundOnly"                                ,
	     "mergerTreeReadPresetParticleCounts"                                    => "mergerTreeConstructorMethod.read.--presetParticleCounts"                                 ,
	     "mergerTreeReadPresetVelocityMaxima"                                    => "mergerTreeConstructorMethod.read.--presetVelocityMaxima"                                 ,
	     "mergerTreeReadPresetVelocityDispersions"                               => "mergerTreeConstructorMethod.read.--presetVelocityDispersions"                            ,
	     "mergerTreeReadBeginAt"                                                 => "mergerTreeConstructorMethod.read.--beginAt"                                              ,
	     "mergerTreeReadOutputTimeSnapTolerance"                                 => "mergerTreeConstructorMethod.read.--outputTimeSnapTolerance"                              ,
	     "mergerTreeReadMissingHostsAreFatal"                                    => "mergerTreeConstructorMethod.read.--missingHostsAreFatal"                                 ,
	     "mergerTreeReadTreeIndexToRootNodeIndex"                                => "mergerTreeConstructorMethod.read.--treeIndexToRootNodeIndex"                             ,
	     "mergerTreeReadAllowBranchJumps"                                        => "mergerTreeConstructorMethod.read.--allowBranchJumps"                                     ,
	     "mergerTreeReadAllowSubhaloPromotions"                                  => "mergerTreeConstructorMethod.read.--allowSubhaloPromotions"                               ,
	     "mergerTreeReadSubresolutionMergingMethod"                              => "mergerTreeConstructorMethod.read.--satelliteMergingTimescalesMethod"                     ,
	     "mergerTreeHaloMassFunctionSamplingMethod"                              => "mergerTreeBuildMassDistributionMethod"                                                   ,
	     "starFormationSpheroidEfficiency"                                       => "starFormationTimescaleSpheroidsMethod.dynamicalTime.--efficiency"                        ,
	     "starFormationSpheroidVelocityExponent"                                 => "starFormationTimescaleSpheroidsMethod.dynamicalTime.--exponentVelocity"                  ,
	     "starFormationSpheroidMinimumTimescale"                                 => "starFormationTimescaleSpheroidsMethod.dynamicalTime.--timescaleMinimum"                  ,
	     "starFormationTimescaleSpheroidsVelocityMaximumScalingTimescale"        => "starFormationTimescaleSpheroidsMethod.velocityMaxScaling.--timescale"                    ,
	     "starFormationTimescaleSpheroidsVelocityMaximumScalingVelocityExponent" => "starFormationTimescaleSpheroidsMethod.velocityMaxScaling.--exponentVelocity"             ,
	     "starFormationTimescaleSpheroidsVelocityMaximumScalingRedshiftExponent" => "starFormationTimescaleSpheroidsMethod.velocityMaxScaling.--exponentRedshift"             ,
	     "starFormationDiskEfficiency"                                           => "starFormationTimescaleDisksMethod.dynamicalTime.--efficiency"                            ,
	     "starFormationDiskVelocityExponent"                                     => "starFormationTimescaleDisksMethod.dynamicalTime.--exponentVelocity"                      ,
	     "starFormationDiskMinimumTimescale"                                     => "starFormationTimescaleDisksMethod.dynamicalTime.--timescaleMinimum"                      ,
	     "starFormationTimescaleDisksVelocityMaximumScalingTimescale"            => "starFormationTimescaleDisksMethod.velocityMaxScaling.--timescale"                        ,
	     "starFormationTimescaleDisksVelocityMaximumScalingVelocityExponent"     => "starFormationTimescaleDisksMethod.velocityMaxScaling.--exponentVelocity"                 ,
	     "starFormationTimescaleDisksVelocityMaximumScalingRedshiftExponent"     => "starFormationTimescaleDisksMethod.velocityMaxScaling.--exponentRedshift"                 ,
	     "starFormationTimescaleDisksFixedTimescale"                             => "starFormationTimescaleDisksMethod.fixed.--timescale"                                     ,
	     "starFormationTimescaleIntegratedSurfaceDensityTolerance"               => "starFormationTimescaleDisksMethod.intgrtdSurfaceDensity.--tolerance"                     ,
	     "starFormationTimescaleDisksHaloScalingTimescale"                       => "starFormationTimescaleDisksMethod.haloScaling.--timescale"                               ,
	     "starFormationTimescaleDisksHaloScalingVirialVelocityExponent"          => "starFormationTimescaleDisksMethod.haloScaling.--exponentVelocityVirial"                  ,
	     "starFormationTimescaleDisksHaloScalingRedshiftExponent"                => "starFormationTimescaleDisksMethod.haloScaling.--exponentRedshift"                        ,
	     "barInstabilityMethod"                                                  => "galacticDynamicsBarInstabilityMethod"                                                    ,
	     "stabilityThresholdStellar"                                             => "galacticDynamicsBarInstabilityMethod.efstathiou1982.--stabilityThresholdStellar"         ,
	     "stabilityThresholdGaseous"                                             => "galacticDynamicsBarInstabilityMethod.efstathiou1982.--stabilityThresholdGaseous"         ,
	     "stabilityTimescaleMinimum"                                             => "galacticDynamicsBarInstabilityMethod.efstathiou1982.--timescaleMinimum"                  ,
	     "harrassmentMassThreshold"                                              => "galacticDynamicsBarInstabilityMethod.efstathiou1982Tidal.--massThresholdHarrassment"     ,
	     "tidalStrippingMassLossRateDisksMethod"                                 => "tidalStrippingDisksMethod"                                                               ,
	     "tidalStrippingMassLossRateDiskSimpleFractionalRateMaximum"             => "tidalStrippingDisksMethod.simple.--rateFractionalMaximum"                                ,
	     "tidalStrippingMassLossRateSpheroidsMethod"                             => "tidalStrippingSpheroidsMethod"                                                           ,
	     "tidalStrippingMassLossRateSpheroidSimpleFractionalRateMaximum"         => "tidalStrippingSpheroidsMethod.simple.--rateFractionalMaximum"                            ,
	     "ramPressureStrippingMassLossRateDisksMethod"                           => "ramPressureStrippingDisksMethod"                                                         ,
	     "ramPressureStrippingMassLossRateDiskSimpleFractionalRateMaximum"       => "ramPressureStrippingDisksMethod.simple.--rateFractionalMaximum"                          ,
	     "ramPressureStrippingMassLossRateSpheroidsMethod"                       => "ramPressureStrippingSpheroidsMethod"                                                     ,
	     "ramPressureStrippingMassLossRateSpheroidSimpleFractionalRateMaximum"   => "ramPressureStrippingSpheroidsMethod.simple.--rateFractionalMaximum"                      ,
	     "satellitesTidalFieldMethod"                                            => "satelliteTidalFieldMethod"                                                               ,
	     "hotHaloColdModeCoredIsothermalCoreRadiiMethod"                         => "hotHaloColdModeCoreRadiiMethod"                                                          ,
	     "coldModeIsothermalCoreRadiusOverVirialRadius"                          => "hotHaloColdModeCoreRadiiMethod.virialFraction.--coreRadiusOverVirialRadius"              ,
	     "hotHaloRamPressureStrippingTimescaleMethod"                            => "hotHaloRamPressureTimescaleMethod"                                                       ,
	     "ramPressureStrippingFormFactor"                                        => "hotHaloRamPressureStrippingMethod.font2008.--formFactor"                                 ,
	     "starFormationFrequencyKMT09"                                           => "starFormationRateSurfaceDensityDisksMethod.krumholz2009.--frequencyStarFormation"        , 
	     "molecularComplexClumpingFactorKMT09"                                   => "starFormationRateSurfaceDensityDisksMethod.krumholz2009.--clumpingFactorMolecularComplex", 
	     "molecularFractionFastKMT09"                                            => "starFormationRateSurfaceDensityDisksMethod.krumholz2009.--molecularFractionFast"         , 
	     "assumeMonotonicSurfaceDensityKMT09"                                    => "starFormationRateSurfaceDensityDisksMethod.krumholz2009.--assumeMonotonicSurfaceDensity" , 
	     "accretionColdModeShockStabilityThreshold"                              => "accretionHaloMethod.coldMode.--thresholdStabilityShock"                                  ,
	     "accretionColdModeShockStabilityTransitionWidth"                        => "accretionHaloMethod.coldMode.--widthTransitionStabilityShock"                            ,
	     "accretionHalosSimpleNegativeAccretionAllowed"                          => "accretionHaloMethod.simple.--accretionNegativeAllowed"                                   ,
	     "accretionHalosSimpleAccreteNewGrowthOnly"                              => "accretionHaloMethod.simple.--accretionNewGrowthOnly"                                     ,
	     "reionizationSuppressionVelocity"                                       => "accretionHaloMethod.simple.--velocitySuppressionReionization"                            ,
	     "reionizationSuppressionRedshift"                                       => "accretionHaloMethod.simple.--redshiftReionization"                                       ,
	     "virialDensityContrastPercolationLinkingLength"                         => "virialDensityContrastMethod.percolation.--linkingLength"
	 },
	 values        =>
	 {
	     mergerTreeBuilderMethod        =>
	     {
		 "Cole2000"                  => "cole2000"
	     },
		   hotHaloRamPressureStrippingMethod =>
	       {
		   "Font2008" => "font2008"
	       },  
		   hotHaloRamPressureForceMethod =>
	       {
		   "Font2008" => "font2008"
	       },  
	      starFormationTimescaleDisksMethod =>
	     {
	      "integratedSurfaceDensity" => "intgrtdSurfaceDensity"
	     },
	     starFormationRateSurfaceDensityDisksMethod => 
	     {
		 "KMT09" => "krumholz2009"
	     }, 
	     starFormationExpulsiveFeedbackDisksMethod =>
	     {
		 "null" => "zero"
             },
	     starFormationExpulsiveFeedbackSpheroidsMethod =>
	     {
		 "null" => "zero"
	     },
	     hotHaloColdModeCoredIsothermalCoreRadiiMethod =>
	     {
		 "virialRadiusFraction" => "virialFraction"
	     },
	     barInstabilityMethod =>	 
	     {
		 "ELN"       => "efstathiou1982",
		 "ELN+tidal" => "efstathiou1982Tidal",
		 "null"      => "stable"
	     },
	     treeNodeMethodDisk      =>
	     {
		 "exponential"           => {
		                             value => "standard",
		                             new   => [
				    	               {
			                                name  => "diskMassDistribution",
			                                value => "exponentialDisk"
			                               }
			                              ]
		                            }
	     },
	     darkMatterProfileShapeMethod   =>
	     {
		 "Gao2008"                   => "gao2008"
	     },
	     transferFunctionMethod         =>
	     {
		 "null"                      => "identity"        ,
		 "Eisenstein-Hu1999"         => "eisensteinHu1999"
	     },
	     stellarPopulationSpectraMethod =>
	     {
	         "Conroy-White-Gunn2009"     => "FSPS"
             },
	     haloSpinDistributionMethod =>
	     {
	         "Bett2007"                  => "bett2007"
	     },
	     criticalOverdensityMethod =>	     {
		 "Kitayama-Suto1996"         => "kitayamaSuto1996"             ,
		 "sphericalTopHat"           => "sphericalCollapseMatterLambda",
		 "sphericalTopHatDarkEnergy" => "sphericalCollapseMatterDE"
	     },
	     mergerTreeImporterMethod        =>
	     {
	         "sussing"                   => "sussingASCII"
	     },
	     starFormationFeedbackDisksMethod =>    
	     {
		 "Creasey2012"               => "creasey2012"
	     }
	 }
	     "hotHaloOutflowReturnRate"                                  => "hotHaloOutflowReincorporationMethod.haloDynamicalTime.--multiplier"
    }
    );

# Define known defaults.
my %knownDefaults =
    (
     "cosmologyParametersMethod"            => "simple"       ,
     "mergerTreeBuilderMethod"              => "cole2000"     ,
     "cosmologicalMassVarianceMethod"       => "filteredPower",
     "cosmologyParametersMethod"            => "simple"       ,
     "darkMatterProfileConcentrationMethod" => "gao2008"      ,
     "powerSpectrumPrimordialMethod"        => "powerLaw"     ,
     "haloMassFunctionMethod"               => "tinker2008"   ,
     "mergerTreeOutput"                     => "true"
  );

# Parse the input file.
my $parser     = XML::LibXML->new();
my $input      = $parser->parse_file($inputFileName);
my @parameterSets;
if ( $input->findnodes('parameters') ) {
    @parameterSets = $input->findnodes('parameters');
} elsif ( $input->findnodes('parameterGrid') ) {
    my $parameterGrid = $input->findnodes('parameterGrid')->[0];
    if ( $parameterGrid->findnodes('parameters') ) {
	@parameterSets = $parameterGrid->findnodes('parameters');
    } else {
	die('can not find parameters')
    }
} else {
    die('can not find parameters')
}

# Write starting message.
print "Translating file: ".$inputFileName."\n";

# Iterate over parameter sets.
foreach my $parameters ( @parameterSets ) {
    &Translate($parameters,1,$inputFileName);
}

# Output the resulting file.
my $pp = XML::LibXML::PrettyPrint->new(indent_string => "  ");
$pp->pretty_print($input)
    if ( $options{'prettyify'} eq "yes" );
my $serialized = $input->toString();
$serialized =~ s/><!\-\-/>\n\n  <!\-\-/gm;
$serialized =~ s/></>\n\n  </gm;
open(my $outputFile,">",$outputFileName);
print $outputFile $serialized;
close($outputFile);

exit;

sub Translate {
    my $parameters    = shift();
    my $rootLevel     = shift();
    my $inputFileName = shift();

    # Set initial input/output versions.
    my $inputVersion  = $options{'inputVersion' };
    my $outputVersion = $options{'outputVersion'};

    # Check for a version element.
    my $version    = $parameters->findnodes('version')->[0];
    if ( defined($version) ) {
	$inputVersion = $version->textContent()
	    unless ( exists($optionsDefined{'inputVersion'}) );
	$version->firstChild()->setData("0.9.4");
    } elsif ( $rootLevel ) {
	my $versionNode = $input->createElement ("version");
	my $newVersion  = $input->createTextNode("0.9.4"  );
	my $newBreak    = $input->createTextNode("\n  "   );
	$versionNode->addChild($newVersion);
	$parameters->insertBefore($versionNode,$parameters->firstChild());
	$parameters->insertBefore($newBreak   ,$parameters->firstChild());
    }
    
    # Determine output format version.
    unless ( exists($options{'outputFormatVersion'}) ) {
	$options{'outputFormatVersion'} = 1;
	my @outputVersionSubstrings = split(/\./,$outputVersion);
	$options{'outputFormatVersion'} = 2
	    if (
		$outputVersionSubstrings[0] >   0
		||
		(
		 $outputVersionSubstrings[0] == 0
		 &&    
		 $outputVersionSubstrings[1] >  9
		)
		||
		(
		 $outputVersionSubstrings[0] == 0
		 &&    
		 $outputVersionSubstrings[1] == 9
		 &&    
		 $outputVersionSubstrings[2] >  3
		)
	    );
    }
    
    # Check for a format version element.
    my $formatVersion = $parameters->findnodes('formatVersion')->[0];
    if ( defined($formatVersion) ) {
	$options{'inputFormatVersion'} = $formatVersion->textContent();
	$formatVersion->firstChild()->setData($options{'outputFormatVersion'});
    } elsif ( $rootLevel ) {
	my $formatVersionNode = $input->createElement ("formatVersion"                );
	my $newFormatVersion  = $input->createTextNode($options{'outputFormatVersion'});
	my $newBreak          = $input->createTextNode("\n  "                         );
	$formatVersionNode->addChild($newFormatVersion);
	$parameters->insertBefore($formatVersionNode,$parameters->firstChild());
	$parameters->insertBefore($newBreak         ,$parameters->firstChild());
    }

    # Validate the parameter file.
    if ( $options{'validate'} eq "yes" ) {
	system($ENV{'GALACTICUS_EXEC_PATH'}."/scripts/aux/validateParameters.pl ".$inputFileName);
	die('input file "'.$inputFileName.'"is not a valid Galacticus parameter file')
	    unless ( $? == 0 );
    }
    
    # Iterate over translations.
    foreach my $translation ( @translations ) {
	# Skip irrelevant translation.
	next
	    unless ( $translation->{'inputVersion'} eq $inputVersion );
	# Report.
	print "Translating from v".$translation->{'inputVersion'}." to v".$translation->{'outputVersion'}."\n";
	$inputVersion = $translation->{'outputVersion'};
	# Find parameter nodes.
	my @parameterNodes;
	if ( $options{'inputFormatVersion'} <= 1 ) {
	    @parameterNodes = $parameters->findnodes('parameter');
	} else {
	    @parameterNodes = map {( $_->exists('value') || $_->hasAttribute('value') ) ? $_ : ()} $parameters->findnodes('*');
	}
	# Apply translations.
	for my $parameter ( @parameterNodes ) {
	    # Get name and value text elements.
	    my $name;
	    my $nameText;
	    my @allValues;
	    if ( $options{'inputFormatVersion'} <= 1 ) {
		$name      = $parameter->findnodes('name' )->[0]->firstChild();
		$nameText  = $name->textContent();
		@allValues = $parameter->findnodes('value');
	    } else {
		$name     = $parameter;
		$nameText = $name->nodeName();
		if ( $parameter->exists('value') ) {
		    @allValues = $parameter->findnodes('value');
		} else {
		    @allValues = $parameter;
		}
	    }
	    # Ignore parameters with no value - sometimes this is used as a placeholder.
	    next
		if ( ! grep {my $valueText = $_->isSameNode($name) ? $_->getAttribute('value') : $_->firstChild()->textContent(); $valueText ne ""} @allValues );	    
	    # Translate names.
	    if ( exists($translation->{'names'}->{$nameText}) ) {
		print "   translate parameter name: ".$nameText." --> ".$translation->{'names'}->{$nameText}."\n";
		if ( $options{'inputFormatVersion'} <= 1 ) {
		    $name->setData    ($translation->{'names'}->{$nameText});
		} else {
		    my $leafName = $translation->{'names'}->{$nameText};
		    my $valueTo;
		    unless ( $translation->{'names'}->{$nameText} =~ m/\-\-/ ) {
			if ( $leafName =~ m/(.*)\.(.*)\./ ) {
			    $valueTo    = $2;
			}
		    }
		    $name->setNodeName($leafName);
		    $name->setAttribute('value',$valueTo)
			if ( $valueTo );
		}
	    }
	    # Translate values.
	    foreach my $value ( @allValues ) {
		if ( exists($translation->{'values'}->{$nameText}) ) {
		    # Split values.
		    my $valuesText;
		    if ( $value->isSameNode($name) ) {
			$valuesText = $value->getAttribute('value');
		    } else {
			$valuesText = $value->firstChild()->textContent();
		    }
		    $valuesText =~ s/^\s*//;
		    $valuesText =~ s/\s*$//;
		    my @values;
		    if ( $translation->{'inputVersion'} eq "0.9.0" ) {
			# For v0.9.0, we cannot split as method names were permitted to contain spaces.
			push(@values,$valuesText);
		    } else {
			@values = split(/\s+/,$valuesText);
		    }
		    foreach my $thisValue ( @values ) {
			if ( exists($translation->{'values'}->{$nameText}->{$thisValue}) ) {
			    print "   translate parameter value: ".$nameText."\n";
			    if ( ref($translation->{'values'}->{$nameText}->{$thisValue}) ) {
				my $newValue = $translation->{'values'}->{$nameText}->{$thisValue};
				print "                                 ".$thisValue." --> ".$newValue->{'value'}."\n";
				$thisValue = $newValue->{'value'};
				if ( exists($newValue->{'new'}) ) {
				    foreach my $newParameter ( &List::ExtraUtils::as_array($newValue->{'new'}) ) {
					print "      add parameter: ".$newParameter->{'name'}." = ".$newParameter->{'value'}."\n";
					if ( $options{'inputFormatVersion'} <= 1 ) {
					    my $parameterNode = $input->createElement (                "parameter" );
					    my $name          = $input->createElement (                "name"      );
					    my $value         = $input->createElement (                "value"     );
					    my $nameText      = $input->createTextNode($newParameter->{'name'     });
					    my $valueText     = $input->createTextNode($newParameter->{'value'    });
					    $name ->addChild($nameText );
					    $value->addChild($valueText);
					    $parameterNode->addChild($input->createTextNode("\n    "));
					    $parameterNode->addChild($name  );
					    $parameterNode->addChild($input->createTextNode("\n    "));
					    $parameterNode->addChild($value );
					    $parameterNode->addChild($input->createTextNode("\n  "));
					    $parameters->insertAfter($parameterNode,$parameter);
					    $parameters->insertAfter($input->createTextNode("\n  "),$parameter);
					} else {
					    my $parameterNode = $input->createElement($newParameter->{'name'});
					    $parameterNode->setAttribute('value',$newParameter->{'value'});
					    $parameters->insertAfter($parameterNode,$parameter);
					    $parameters->insertAfter($input->createTextNode("\n  "),$parameter);
					}
				    }
				}
			    } else {
				print "                                 ".$thisValue." --> ".$translation->{'values'}->{$nameText}->{$thisValue}."\n";
				$thisValue = $translation->{'values'}->{$nameText}->{$thisValue};
			    }
			}
			if ( $value->isSameNode($name) ) {
			    $value->setAttribute('value',join(" ",@values));
			} else {
			    $value->firstChild()->setData(join(" ",@values));
			}
		    }
		}
	    }
	}
	# Finished if output version is reached.
	last
	    if ( $outputVersion eq $inputVersion );
    }
    
    # Handle transition from old to new.
    if ( $options{'outputFormatVersion'} >= 2 && $options{'inputFormatVersion'} < 2) {
	print "Converting to new format (v2)...\n";
	for my $parameter ( $parameters->findnodes('parameter') ) {
	    # Get name and value text elements.
	    my $name   = $parameter->findnodes('name' )->[0]->firstChild();
	    my @values = $parameter->findnodes('value');
	    # Create the new node.
	    my $parameterNode = $input->createElement($name->textContent());
	    # Determine if we can use a value attribute or not.
	    my $useAttribute = 1;
	    $useAttribute = 0
		if ( scalar(@values) > 1 );
	    foreach my $valueNode ( @values ) {
		$useAttribute = 0
		    if ( $valueNode->findnodes('*') );
	    }
	    # Add values.
	    foreach my $valueNode ( @values ) {
		my $value  = $valueNode->firstChild();
		# Find any subparameters.
		my @subParameters = $valueNode->findnodes('*');
		if ( $useAttribute ) {
		    $parameterNode->setAttribute('value',$value->textContent());
		} else {
		    &Translate($valueNode,0,$inputFileName)
			if ( @subParameters );
		    $parameterNode->addChild($input->createTextNode("\n"));
		    $parameterNode->addChild($valueNode);
		}
	    }
	    $parameterNode->addChild($input->createTextNode("\n"))
		unless ( $useAttribute );
	    # Insert the new node.
	    $parameters->insertAfter($parameterNode,$parameter);
	    # Remove the old parameter.
	    $parameters->removeChild($parameter);
	}
    }
    # Put subparameters into their host parameter.
    for my $parameter ( $parameters->getChildrenByTagName("*") ) {
	if ( $parameter->nodeName() =~ m/(.+)\-\-(.+)/ ) {
	    my $hostName = $1;
	    my $subName  = $2;
	    (my $hostNameTrimmed = $hostName) =~ s/\..+\.//;
	    $parameters->removeChild($parameter->nextSibling());
	    my $sibling = $parameter->nextSibling();
	    my $hostFound;
	    $parameter->setNodeName($subName);
	    for my $hostParameter ( $parameters->getChildrenByTagName("*") ) {
		if ( $hostParameter->nodeName() eq $hostNameTrimmed ) {
		    my $hostChildren = $hostParameter->getChildrenByTagName("*");
		    $hostParameter->addChild($input     ->createTextNode("\n"      ))
			if ( $hostChildren->size() == 0 );
		    $hostParameter->addChild($input     ->createTextNode("  "      ));
		    $hostParameter->addChild($parameters->removeChild   ($parameter));
		    $hostParameter->addChild($input     ->createTextNode("\n  "    ));
		    $hostFound = 1;
		}
	    }
	    unless ( $hostFound ) {
		# Create the new node.
		(my $hostLeafName = $hostName) =~ s/\..+\.//;
		my $defaultValue;
		if ( $hostName =~ m/\..+\./ ) {
		    ($defaultValue = $hostName) =~ s/.*\.(.+)\./$1/;
		} else {
		    die('parametersMigrate.pl: attempting to insert a "'.$hostLeafName.'" element, but no default value is known')
			unless ( exists($knownDefaults{$hostLeafName}) );
		    $defaultValue = $knownDefaults{$hostLeafName};
		}
		my $parameterNode = $input->createElement($hostLeafName);
		$parameterNode->setAttribute('value',$defaultValue);
		$parameterNode->addChild($input     ->createTextNode("\n    "  ));
		$parameterNode->addChild($parameters->removeChild   ($parameter));
		$parameterNode->addChild($input     ->createTextNode("\n  "    ));
		$parameters->insertBefore($parameterNode,$sibling);
	    }
	}
    }
    # Strip out value extensions.
    for my $parameter ( $parameters->getChildrenByTagName("*") ) {
	if ( $parameter->nodeName() =~ m/([^\.]+)\./ ) {
	    my $nodeName = $1;
	    $parameter->setNodeName($nodeName);
	}
    }
    # Search for duplicated parameters.
    my %duplicates =
	(
	 mergerTreeOperatorMethod => "sequence"
	);
   for my $parameter ( $parameters->getChildrenByTagName("*") ) {
	my $nodeName = $parameter->nodeName();
	my $duplicates = $parameters->getChildrenByTagName($nodeName);
	if ( $duplicates->size() > 1 && exists($duplicates{$nodeName}) ) {
	    my $sibling = $parameter->nextSibling();
	    my $wrapperNode = $input->createElement($nodeName);
	    $wrapperNode->setAttribute('value',$duplicates{$nodeName});
	    foreach my $duplicate ( $duplicates->get_nodelist() ) {
		$wrapperNode->addChild($input     ->createTextNode("\n    "  ));
		$wrapperNode->addChild($parameters->removeChild   ($duplicate));
		$wrapperNode->addChild($input     ->createTextNode("\n  "    ));
	    }
	    $parameters->insertBefore($wrapperNode,$sibling);
	}
    }
}
