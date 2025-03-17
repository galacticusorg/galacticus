#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use Text::Template qw(fill_in_string);
use XML::Simple;
use PDL;
use PDL::IO::HDF5;
use Galacticus::Options;
use Galacticus::Constraints::HaloMassFunctions qw(iterate);
use List::Util;
use List::ExtraUtils;
use Data::Dumper;

# Script to generate content for halo mass function constraint pipeline.
# Andrew Benson 8-April-2022

# Generates the haloMassFunctionConfig.xml file.
# Also generates haloMassFunctionBase_*.xml files for all halo mass functions (including inserting the measured environment properties).

# Get command line options.
my %options =
    (
     binAverage                   => "true" , # If true halo mass functions will be averaged over each bin.
     initializeToPosteriorMaximum => "false", # f true, initialize state to the posterior maximumSpecify the state initializor to use.
     removeAccelerator            => "false", # If true any `<haloMassFunction value="accelerator"/>` will be removed.
     removeDetectionEfficiency    => "false", # If true any `<haloMassFunction value="detectionEfficiency"/>` will be removed.
     removeErrorConvolved         => "false", # If true any `<haloMassFunction value="errorConvolved"/>` will be removed.
     removeSimulationVariance     => "false", # If true any `<haloMassFunction value="simulationVariance"/>` will be removed.
     removeMultiplier             => "false"  # If true any `<haloMassFunction value="multiplier"/>` (used for isolation bias models) will be removed.
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);
die("no `--pipelinePath` option given")
    unless ( exists($options{'pipelinePath'}) );
die("no `--outputDirectory` option given")
    unless ( exists($options{'outputDirectory'}) );

# Export options.
$content::binAverage = $options{'binAverage'};

# Ensure paths are correctly suffixed.
foreach my $path ( 'pipelinePath', 'outputDirectory' ) {
    $options{$path} .= "/"
	unless ( $options{$path} =~ m/\/$/ );
}

# Set global parameters.
$content::countParticlesMinimum = 100  ;
$content::fractionMassPrimary   =   0.1;
$content::pipelinePath          = $options{'pipelinePath'   };
$content::outputDirectory       = $options{'outputDirectory'};

# Parse the simulations definition file.
my $xml = new XML::Simple();
my $simulations = $xml->XMLin(
    $options{'pipelinePath'}."haloMassFunctionSimulations.xml",
    ForceArray => [ "suite"         , "group"         , "simulation"          ],
    KeyAttr    => {  suite => "name",  group => "name",  simulation => "name" }
    );

# Initialize hashes of isolationBias and perturbation parameters needed.
my %isolationBiases;
my %perturbations ;

# Iterate over simulation suites/groups building isolation bias and perturbation parameters as needed.
print "Creating parameters...\n";
my $haveGroups = 0;
foreach my $entry ( &iterate($simulations,\%options, stopAfter => "group") ) {
    $haveGroups = 1;
    print
	"  ".
	$entry->{'suite'}->{'name'}."\t".
	$entry->{'group'}->{'name'}."\n";
    # Generate isolation bias parameters.
    if ( $entry->{'suite'}->{'includeIsolationBias'}->{'value'} eq "true" && $options{'removeMultiplier'} eq "false" ) {
	($entry->{'group'}->{'labelIsolationBias'} = $entry->{'suite'}->{'name'}.$entry->{'group'}->{'name'}) =~ s/:/_/g;
	$entry->{'group'}->{'isolationBias'     } =
	    "haloMassFunctionParameters/isolationBias"          .$entry->{'group'}->{'labelIsolationBias'}.
	    "  haloMassFunctionParameters/isolationBiasExponent".$entry->{'group'}->{'labelIsolationBias'};
	++$isolationBiases{$entry->{'group'}->{'labelIsolationBias'}};
    } else {
	$entry->{'group'}->{'labelIsolationBias'} = "";
    }
    # Generate perturbation parameters.
    if ( $entry->{'suite'}->{'includePerturbation'}->{'value'} eq "true" && $options{'removeSimulationVariance'} eq "false" ) {
	$entry->{'group'}->{'labelPerturbation'} = $entry->{'suite'}->{'name'}.$entry->{'group'}->{'name'};
	$entry->{'group'}->{'perturbation'     } = "haloMassFunctionParameters/perturbation".$entry->{'group'}->{'labelPerturbation'};
	++$perturbations{$entry->{'group'}->{'labelPerturbation'}};
    } else {
	$entry->{'group'}->{'perturbation'} = "";
    }
}
die("no groups match this selection")
    unless ( $haveGroups );

# Iterate over suites building custom halo mass function parameter files.
foreach my $entry ( &iterate($simulations,\%options, stopAfter => "suite") ) {
    my $haloMassFunctionParameterFileNameOriginal = $options{'pipelinePath'   }."haloMassFunction_".$entry->{'suite'}->{'name'}.".xml";
    my $haloMassFunctionParameterFileNameNew      = $options{'outputDirectory'}."haloMassFunction_".$entry->{'suite'}->{'name'}.".xml";
    my $haloMassFunctionParameters                = $xml->XMLin($haloMassFunctionParameterFileNameOriginal);
    foreach my $removal ( "accelerator", "detectionEfficiency", "errorConvolved", "simulationVariance", "multiplier" ) {
	next
	    unless ( $options{"remove".ucfirst($removal)} eq "true" );
	my $haloMassFunctionPrior = $haloMassFunctionParameters;
	my $haloMassFunction      = $haloMassFunctionParameters->{'haloMassFunction'};
	while ( $haloMassFunction ) {
	    if ( $haloMassFunction->{'value'} eq $removal ) {
		$haloMassFunctionPrior->{'haloMassFunction'} = $haloMassFunction->{'haloMassFunction'};
	    } else {
		$haloMassFunctionPrior = $haloMassFunction;
	    }
	    $haloMassFunction = $haloMassFunction->{'haloMassFunction'};
	}
    }
    open(my $haloMassFunctionParameterFileNew,">",$haloMassFunctionParameterFileNameNew);
    print $haloMassFunctionParameterFileNew $xml->XMLout($haloMassFunctionParameters,rootName => "parameters");
    close($haloMassFunctionParameterFileNew);
}

# Iterate over all suites/groups/simulations/redshifts/realizations generating parameter files and config entries.
print "Generating base parameter files...\n";
my $configLikelihood;
my %detectionEfficiencyClasses;
my $haveModels = 0;
foreach $content::entry ( &iterate($simulations,\%options) ) {
    $haveModels = 1;
    print
	"  ".
	     $content::entry->{'suite'      }->{'name'}."\t".
     	     $content::entry->{'group'      }->{'name'}."\t".
	     $content::entry->{'simulation' }->{'name'}."\t".
	     $content::entry->{'realization'}          ."\t".
	"z=".$content::entry->{'redshift'   }          ."\n";
    # Determine the minimum and maximum halo masses.
    $content::massHaloMinimum = sprintf("%11.5e",$content::countParticlesMinimum*$content::entry->{'group'}->{'massParticle'});
    $content::massHaloMaximum = sprintf("%11.5e",$content::fractionMassPrimary  *$content::entry           ->{'massPrimary' })
	if ( $content::entry->{'suite'}->{'limitMassMaximum'}->{'value'} eq "primaryFraction" );
    # Generate file names.
    $content::fileNameBase   = $options{'outputDirectory'}  ."haloMassFunctionBase_".$content::entry->{'suite'}->{'name'}."_".$content::entry->{'group'}->{'name'}."_".$content::entry->{'simulation'}->{'name'}."_".$content::entry->{'realization'}."_z".$content::entry->{'redshift'}.".xml" ;
    $content::fileNameTarget = "\%DATASTATICPATH\%/darkMatter/haloMassFunction_"    .$content::entry->{'suite'}->{'name'}."_".$content::entry->{'group'}->{'name'}."_".$content::entry->{'simulation'}->{'name'}."_".$content::entry->{'realization'}."_z".$content::entry->{'redshift'}.".hdf5";
    # Determine detection efficiency class.
    (my $suiteName = $content::entry->{'suite'}->{'name'}) =~ s/://g;
    $content::class = $suiteName.(exists($content::entry->{'group'}->{'detectionEfficiencyClass'}) ? $content::entry->{'group'}->{'detectionEfficiencyClass'} : "");
    if ( $options{'removeDetectionEfficiency'} eq "false" ) {
	++$detectionEfficiencyClasses{$content::class};
	$content::detectionEfficiency1 = fill_in_string(<<'CODE', PACKAGE => 'content');
haloMassFunctionParameters/massMinimumParticleCount{$class} haloMassFunctionParameters/efficiencyAtMassMinimum{$class}
CODE
	$content::detectionEfficiency2 = fill_in_string(<<'CODE', PACKAGE => 'content');
haloMassFunctionParameters/exponentMassDetection{$class}    haloMassFunctionParameters/exponentRedshiftDetection{$class}
CODE
    } else {
	$content::detectionEfficiency1 = "";
	$content::detectionEfficiency2 = "";
    }
    # Determine if the likelihood should be heated.
    $content::likelihoodHeatingOpener = "";
    $content::likelihoodHeatingCloser = "";
    if ( exists($content::entry->{'options'}->{'likelihoodTemperature'}) ) {
	$content::likelihoodHeatingOpener .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <posteriorSampleLikelihood value="heated">
       <temperature value="{$entry->{'options'}->{'likelihoodTemperature'}}"/>
CODE
	$content::likelihoodHeatingCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
</posteriorSampleLikelihood>
CODE
    }
    # Generate XML for the model likelihood entry for this simulation.
    $configLikelihood .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <!-- Suite: {$entry->{'suite'}->{'name'}}; Group: {$entry->{'group'}->{'name'}}; Simulation: {$entry->{'simulation'}->{'name'}} -->
    <parameterMap value="haloMassFunctionParameters/a                                haloMassFunctionParameters/b
                         haloMassFunctionParameters/c                                haloMassFunctionParameters/p
                         haloMassFunctionParameters/q                                haloMassFunctionParameters/normalization 

                         haloMassFunctionParameters/cW0                              haloMassFunctionParameters/beta0
                         haloMassFunctionParameters/cW1                              haloMassFunctionParameters/beta1
                         haloMassFunctionParameters/wavenumberScaledMinimum          haloMassFunctionParameters/powerSpectrumSmoothingWidth

                         haloMassFunctionParameters/artificialNormalization          haloMassFunctionParameters/artificialExponent
                         haloMassFunctionParameters/artificialCountParticles

                         {$detectionEfficiency1}
                         {$detectionEfficiency2}
 
                         varianceFractionalModelDiscrepancy

                         {$entry->{'group'}->{'perturbation' }}
                         {$entry->{'group'}->{'isolationBias'}}
                        "/>
    <parameterInactiveMap value=""     ignoreWarnings="true"/>
{$likelihoodHeatingOpener}
    <posteriorSampleLikelihood value="haloMassFunction">
      <!-- Options matched to those of Benson (2017; https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.3454B) -->
      <baseParametersFileName value="{$fileNameBase}"       />
      <fileName               value="{$fileNameTarget}"     />
      <redshift               value="{$entry->{'redshift'}}"/>
      <massRangeMinimum       value="{$massHaloMinimum}"/> <!-- {$countParticlesMinimum} times zoom-in {$entry->{'suite'}->{'name'}} {$entry->{'group'}->{'name'}} particle mass -->
CODE
    if ( $content::entry->{'suite'}->{'limitMassMaximum'}->{'value'} eq "primaryFraction" ) {
	$configLikelihood .= fill_in_string(<<'CODE', PACKAGE => 'content');
      <massRangeMaximum       value="{$massHaloMaximum}"/> <!-- {$fractionMassPrimary} of the target halo mass                                                                   -->
CODE
    }
    $configLikelihood .= fill_in_string(<<'CODE', PACKAGE => 'content');
      <binCountMinimum        value="0"                 />    
      <likelihoodPoisson      value="true"              />
      <binAverage             value="{$binAverage}"     />
    </posteriorSampleLikelihood>
{$likelihoodHeatingCloser}
CODE

    # Generate the base parameter file for this instance.
    my $base = fill_in_string(<<'CODE', PACKAGE => 'content');
<?xml version="1.0" encoding="UTF-8"?>
<parameters>
  <formatVersion>2</formatVersion>
  <version>0.9.4</version>

  <!-- Output control -->
  <outputFileName value="{$outputDirectory}haloMassFunction_{$entry->{'suite'}->{'name'}}_{$entry->{'group'}->{'name'}}_{$entry->{'simulation'}->{'name'}}_{$entry->{'realization'}}_z{$entry->{'redshift'}}.hdf5"/>
  <outputTimes value="list">
    <redshifts value="{$entry->{'redshift'}}"/>
  </outputTimes>  

  <!-- Include cosmology and mass function parameters -->
  <xi:include href="haloMassFunctionParameters.xml"                                                                       xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="{$pipelinePath}simulation_{$entry->{'suite'}->{'name'}}_{$entry->{'group'}->{'name'}}.xml"            xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="{$pipelinePath}cosmology_{$entry->{'suite'}->{'name'}}.xml"                                           xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="{$outputDirectory}haloMassFunction_{$entry->{'suite'}->{'name'}}.xml"                                 xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>
  <xi:include href="{$pipelinePath}transferFunction_{$entry->{'suite'}->{'name'}}_{$entry->{'simulation'}->{'name'}}.xml" xpointer="xpointer(parameters/*)" xmlns:xi="http://www.w3.org/2001/XInclude"/>

  <!-- Detection effiency -->
  <detectionMassMinimumParticleCount value="=[haloMassFunctionParameters/massMinimumParticleCount{$class}]"  ignoreWarnings="true"/>
  <detectionEfficiencyAtMassMinimum  value="=[haloMassFunctionParameters/efficiencyAtMassMinimum{$class}]"   ignoreWarnings="true"/>
  <detectionExponentMass             value="=[haloMassFunctionParameters/exponentMassDetection{$class}]"     ignoreWarnings="true"/>
  <detectionExponentRedshift         value="=[haloMassFunctionParameters/exponentRedshiftDetection{$class}]" ignoreWarnings="true"/>
CODE
    if ( $content::entry->{'suite'}->{'includePerturbation'}->{'value'} eq "true" && $options{'removeSimulationVariance'} eq "false" ) {
	$base .= fill_in_string(<<'CODE', PACKAGE => 'content');

  <!-- Use the simulation variance perturbation relevant to this simulation -->
  <perturbationFractional value="=[haloMassFunctionParameters/perturbation{$entry->{'group'}->{'labelPerturbation'}}]"/>
CODE
    }
    if ( $content::entry->{'suite'}->{'includeIsolationBias'}->{'value'} eq "true" && $options{'removeMultiplier'} eq "false" ) {
	$base .= fill_in_string(<<'CODE', PACKAGE => 'content');

  <!-- Isolation bias -->
  <isolationBias         value="=[haloMassFunctionParameters/isolationBias{$entry->{'group'}->{'labelIsolationBias'}}]"         ignoreWarnings="true"/>
  <isolationBiasExponent value="=[haloMassFunctionParameters/isolationBiasExponent{$entry->{'group'}->{'labelIsolationBias'}}]" ignoreWarnings="true"/>

CODE
    }
    if ( $content::entry->{'suite'}->{'includeEnvironment'}->{'value'} eq "true" ) {
	$base .= fill_in_string(<<'CODE', PACKAGE => 'content');

  <!-- Halo environments -->
  <haloEnvironment value="fixed">
    <massEnvironment value="{$entry->{'environment'}->{'massEnvironment'       }}"                      />
    <overdensity     value="{$entry->{'environment'}->{'overdensityEnvironment'}}"                      />
    <redshift        value="{$entry                 ->{'redshift'              }}" ignoreWarnings="true"/>
  </haloEnvironment>
CODE
    }
    $base .= fill_in_string(<<'CODE', PACKAGE => 'content');

</parameters>
CODE

    # Generate the base parameter file.
    open(my $baseFile,">".$content::fileNameBase);
    print $baseFile $base;
    close($baseFile);
}
die("no models match this selection")
    unless ( $haveModels );
print "...done\n";

# Generate openers and closers for the config and parameter files.
my $configOpener = fill_in_string(<<'CODE', PACKAGE => 'content');
<?xml version="1.0" encoding="UTF-8"?>
<parameters>
  <!-- Posterior sampling simulation parameter file for constraining to dark matter halo mass functions. -->
  <!-- Andrew Benson (17-September-2020)                                                                 -->  
  <formatVersion>2</formatVersion>
  <version>0.9.4</version>

  <verbosityLevel value="standard"/>

  <task value="posteriorSample">
    <initializeNodeClassHierarchy value="true"/>
  </task>

  <outputFileName value="{$outputDirectory}haloMassFunction.hdf5"/>

  <!-- Likelihood -->
  <posteriorSampleLikelihood value="independentLikelihoods">
    <orderRotation value="byRankOnNode"/> <!-- Rotate likelihoods by the on-node rank to ensure that each process begins with a different power spectrum class. -->
CODE
my $countParameters = 16+scalar(keys(%perturbations))+2*scalar(keys(%isolationBiases))+4*scalar(keys(%detectionEfficiencyClasses));
$content::gammaInitial = 2.35/sqrt($countParameters);
my $configInitializer = fill_in_string(<<'CODE', PACKAGE => 'content');
  </posteriorSampleLikelihood>

CODE
if ( exists($options{'initializeToPosteriorMaximum'}) ) {
    $content::priorLogFileRoot = $options{'initializeToPosteriorMaximum'};
    $configInitializer .= fill_in_string(<<'CODE', PACKAGE => 'content');
  <posteriorSampleStateInitialize value="posteriorMaximumGaussianSphere">
     <logFileRoot      value="{$priorLogFileRoot}"/>
     <radiusSphere     value="1.0e-6"/>
     <radiusIsRelative value="true"  />
  </posteriorSampleStateInitialize>   

CODE
} else {
    $configInitializer .= fill_in_string(<<'CODE', PACKAGE => 'content');
  <posteriorSampleStateInitialize value="gaussianSphere">
     <radiusSphere     value="1.0e-6"/>
     <radiusIsRelative value="true"  />
  </posteriorSampleStateInitialize>   

CODE
}
$configInitializer .= fill_in_string(<<'CODE', PACKAGE => 'content');
  <posteriorSampleDffrntlEvltnProposalSize value="adaptive" >
    <gammaInitial          value="{$gammaInitial}"/>
    <gammaAdjustFactor     value="1.10e+0"/>
    <gammaMinimum          value="1.00e-4"/>
    <gammaMaximum          value="3.00e+0"/>
    <acceptanceRateMinimum value="0.10e+0"/>
    <acceptanceRateMaximum value="0.90e+0"/>
    <updateCount           value="10"     />
    <appendLog             value="false"  />
    <flushLog              value="true"   />
    <restoreFromLog        value="false"  />
    <logFileName           value="{$outputDirectory}haloMassFunctionGamma.log"/>
  </posteriorSampleDffrntlEvltnProposalSize>

  <!-- MCMC -->
  <posteriorSampleSimulation value="differentialEvolution">
    <stepsMaximum           value="100000"                                  />
    <acceptanceAverageCount value="    10"                                  />
    <stateSwapCount         value="    11"                                  /> <!-- Offset swaps from reporting, otherwise we only get reports for swap steps, which gives a biased view of progress. -->
    <slowStepCount          value="    11"                                  />
    <logFileRoot            value="{$outputDirectory}haloMassFunctionChains"/>
    <reportCount            value="    10"                                  />
    <sampleOutliers         value="false"                                   />
    <logFlushCount          value="     1"                                  />
    <appendLogs             value="false"                                   />
    <loadBalance            value="false"                                   />
CODE
my $configResumer = fill_in_string(<<'CODE', PACKAGE => 'content');
  </posteriorSampleLikelihood>

  <posteriorSampleStateInitialize value="resume">
    <logFileRoot  value="{$outputDirectory}haloMassFunctionChains"/>
    <restoreState value="true"                                    />
  </posteriorSampleStateInitialize>   

  <posteriorSampleDffrntlEvltnProposalSize value="adaptive" >
    <gammaInitial          value="2.00e+0"/>
    <gammaAdjustFactor     value="1.10e+0"/>
    <gammaMinimum          value="1.00e-4"/>
    <gammaMaximum          value="3.00e+0"/>
    <acceptanceRateMinimum value="0.10e+0"/>
    <acceptanceRateMaximum value="0.90e+0"/>
    <updateCount           value="10"     />
    <appendLog             value="true"   />
    <flushLog              value="true"   />
    <restoreFromLog        value="true"   />
    <logFileName           value="{$outputDirectory}haloMassFunctionGamma.log"/>
  </posteriorSampleDffrntlEvltnProposalSize>

  <!-- MCMC -->
  <posteriorSampleSimulation value="differentialEvolution">
    <stepsMaximum           value="100000"                                  />
    <acceptanceAverageCount value="    10"                                  />
    <stateSwapCount         value="    20"                                  />
    <slowStepCount          value="    10"                                  />
    <logFileRoot            value="{$outputDirectory}haloMassFunctionChains"/>
    <reportCount            value="    10"                                  />
    <sampleOutliers         value="false"                                   />
    <logFlushCount          value="     1"                                  />
    <appendLogs             value="true"                                    />

CODE
my $configCloser = fill_in_string(<<'CODE', PACKAGE => 'content');
    <posteriorSampleState value="correlation">
      <acceptedStateCount value="100"/>
    </posteriorSampleState>

    <posteriorSampleConvergence value="gelmanRubin">
      <thresholdHatR              value=" 1.20"/>
      <burnCount                  value="200"  />
      <testCount                  value=" 20"  />
      <outlierCountMaximum        value=" 3"   />
      <outlierSignificance        value=" 0.95"/>
      <outlierLogLikelihoodOffset value="60"   />
      <reportCount                value=" 1"   />
      <logFileName                value="{$outputDirectory}haloMassFunctionConvergence.log"/>
    </posteriorSampleConvergence>
    
    <posteriorSampleStoppingCriterion value="correlationLength">
      <stopAfterCount value="1000"/>
    </posteriorSampleStoppingCriterion>

    <posteriorSampleDffrntlEvltnRandomJump   value="adaptive"/>

    <!-- Parameters of the dark matter halo mass function. -->
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/a"/>
      <label value="a" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="uniform">
	<limitLower value=" 0.03"/>
	<limitUpper value="10.00"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
	<median value="0.0"/>
	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/b"/>
      <label value="b" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="uniform">
	<limitLower value="-3.00"/>
	<limitUpper value="+3.00"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
	<median value="0.0"/>
	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
     <modelParameter value="active">
      <name value="haloMassFunctionParameters/c" />
      <label value="c" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="+0.10" />
        <limitUpper value="+5.00" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/p"/>
      <label value="p" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="uniform">
	<limitLower value="-3.0"/>
	<limitUpper value="+3.0"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
	<median value="0.0"/>
	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/q"/>
      <label value="q" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="uniform">
	<limitLower value="-3.00"/>
	<limitUpper value="+3.00"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
	<median value="0.0"/>
	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/normalization"/>
      <label value="A" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="logUniform">
	<limitLower value="1.0e-3"/>
	<limitUpper value="1.0e+3"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="logarithm"/>
      <distributionFunction1DPerturber value="cauchy">
	<median value="0.0"/>
	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>

    <!-- Artificial halo model parameters -->
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/artificialNormalization"/>
      <label value="A" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="logUniform">
	<limitLower value="1.0e-2"/>
	<limitUpper value="1.0e+1"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="logarithm"/>
      <distributionFunction1DPerturber value="cauchy">
	<median value="0.0"/>
	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/artificialExponent"/>
      <label value="\alpha" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="uniform">
	<limitLower value="-3.0"/>
	<limitUpper value="+0.0"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
	<median value="0.0"/>
	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/artificialCountParticles"/>
      <label value="N" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="logUniform">
	<limitLower value="1.0e+0"/>
	<limitUpper value="1.0e+3"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="logarithm"/>
      <distributionFunction1DPerturber value="cauchy">
	<median value="0.0"/>
	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>

    <!-- Window function parameters -->
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/cW0" />
      <label value="c_\mathrm\{W,0\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="0.50" />
        <limitUpper value="6.00" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
      <slow value="false" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/beta0" />
      <label value="\beta_0" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value=" 0.50" />
        <limitUpper value="10.00" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
      <slow value="false" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/cW1" />
      <label value="c_\mathrm\{W,1\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="logUniform">
        <limitLower value="0.1" />
        <limitUpper value="5.0" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="logarithm" />
      <slow value="false" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/beta1" />
      <label value="\beta_1" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="0.1" />
        <limitUpper value="5.0" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="logarithm" />
      <slow value="false" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/wavenumberScaledMinimum" />
      <label value="x_\mathrm\{min\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="0.0" />
        <limitUpper value="5.0" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
      <slow value="false" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/powerSpectrumSmoothingWidth" />
      <label value="\log k_\mathrm\{width\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="logUniform">
        <limitLower value="0.01" />
        <limitUpper value="10.0" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="logarithm" />
      <slow value="false" />
    </modelParameter>

    <modelParameter value="active">
      <name value="varianceFractionalModelDiscrepancy"/>
      <label value="\mathcal\{C\}_\mathrm\{disc\}" ignoreWarnings="true"/>
      <distributionFunction1DPrior value="logUniform">
    	<limitLower value="1.0e-6"/>
    	<limitUpper value="1.0e+0"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="logarithm"/>
      <distributionFunction1DPerturber value="cauchy">
    	<median value="0.0"/>
    	<scale value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
CODE

my $parametersOpener = fill_in_string(<<'CODE', PACKAGE => 'content');
<?xml version="1.0" encoding="UTF-8"?>
<parameters>
  <formatVersion>2</formatVersion>
  <version>0.9.4</version>

  <!-- Controllable parameters of the halo mass function -->
  <haloMassFunctionParameters value="" ignoreWarnings="true">
    <!-- Halo mass function -->
    <a                                value="+0.75"/>
    <b                                value="+0.84"/>
    <c                                value="+2.26"/>
    <p                                value="-1.09"/>
    <q                                value="+0.51"/>
    <normalization                    value="+0.57"/>
    <!-- Artificial halo variance -->
    <artificialNormalization          value="+1.00e+0"/>
    <artificialExponent               value="-0.50e+0"/>
    <artificialCountParticles         value="+1.00e+2"/>
    <!-- ETHOS window function -->
    <cW0                              value="+2.59"/>
    <beta0                            value="+3.51"/>
    <cW1                              value="+0.74"/>
    <beta1                            value="+4.83"/>
    <wavenumberScaledMinimum          value="+0.00"/>
    <powerSpectrumSmoothingWidth      value="+1.00"/>
CODE
my $parametersCloser;

# Add perturbation parameters as needed.
if ( %perturbations ) {
    $configCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');

    <!-- Perturbation model parameters -->
CODE
    $parametersCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <!-- Perturbation model parameters -->
CODE
    foreach $content::perturbationLabel ( sort(keys(%perturbations)) ) {
	$configCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/perturbation{$perturbationLabel}" />
      <label value="\epsilon_\mathrm\{{$perturbationLabel}\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="normal">
        <limitLower value="-10.0" />
        <limitUpper value="+10.0" />
        <mean value="0.0" />
        <variance value="1.0" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
    </modelParameter>
CODE
	$parametersCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <perturbation{$perturbationLabel} value="0.0"/>
CODE
    }
}

# Add isolation bias parameters as needed.
if ( %isolationBiases ) {
    $configCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');

    <!-- Isolation bias model parameters -->
CODE
    $parametersCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <!-- Isolation bias model parameters -->
CODE
    foreach $content::isolationBiasLabel ( sort(keys(%isolationBiases)) ) {
	$configCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/isolationBias{$isolationBiasLabel}" />
      <label value="\mathcal\{I\}_\mathrm\{{$isolationBiasLabel}\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="logNormal">
        <limitLower value="1.0e-3" />
        <limitUpper value="1.0e+3" />
        <mean value="1.0" />
        <variance value="0.25" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
    </modelParameter>    
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/isolationBiasExponent{$isolationBiasLabel}" />
      <label value="\alpha_\mathrm\{{$isolationBiasLabel}\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="-3.0" />
        <limitUpper value="+0.0" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
    </modelParameter>
CODE
	$parametersCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <isolationBias{$isolationBiasLabel}         value="1.0"/>
    <isolationBiasExponent{$isolationBiasLabel} value="0.0"/>
CODE
    }
}

# Add detection efficiency model parameters.
if ( scalar(keys(%detectionEfficiencyClasses)) > 0 ) {
    $configCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');

    <!-- Detection efficiency model parameters -->
CODE
    $parametersCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <!-- Detection efficiency model parameters -->
CODE
}
# Generate parameters for each detection effiency class.
foreach $content::class ( sort(keys(%detectionEfficiencyClasses)) ) {
    $configCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/massMinimumParticleCount{$class}" />
      <label value="N_\mathrm\{min,{$class}\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="logUniform">
        <limitLower value="  1.00" />
        <limitUpper value="100.00" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="logarithm" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/efficiencyAtMassMinimum{$class}" />
      <label value="\epsilon_\mathrm\{min,{$class}\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="+0.00" />
        <limitUpper value="+1.00" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/exponentMassDetection{$class}" />
      <label value="\alpha_\mathrm\{min,{$class}\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="-3.00" />
        <limitUpper value="+0.00" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
    </modelParameter>
    <modelParameter value="active">
      <name value="haloMassFunctionParameters/exponentRedshiftDetection{$class}" />
      <label value="\beta_\mathrm\{min,{$class}\}" ignoreWarnings="true"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0" />
        <scale value="1.0e-4" />
      </distributionFunction1DPerturber>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="-3.00" />
        <limitUpper value="+3.00" />
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity" />
    </modelParameter>
CODE
    $parametersCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
    <massMinimumParticleCount{$class}  value="10.0"/>
    <efficiencyAtMassMinimum{$class}   value=" 1.0"/>
    <exponentMassDetection{$class}     value=" 0.0"/>
    <exponentRedshiftDetection{$class} value=" 0.0"/>
CODE
}

# Finish the closers.
$configCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');

  </posteriorSampleSimulation>

  <!-- Random seed -->
  <randomNumberGenerator value="GSL">
    <seed          value="219" />
    <mpiRankOffset value="true"/>
  </randomNumberGenerator>

</parameters>
CODE
$parametersCloser .= fill_in_string(<<'CODE', PACKAGE => 'content');
  </haloMassFunctionParameters>

  <!-- Parameter controlling the fractional variance due to model discrepancy -->
  <varianceFractionalModelDiscrepancy value="0.0"/>

</parameters>
CODE

# Generate the config file.
open(my $configFile    ,">",$options{'outputDirectory'}."haloMassFunctionConfig.xml"      );
print $configFile     $configOpener     ;
print $configFile     $configLikelihood ;
print $configFile     $configInitializer;
print $configFile     $configCloser     ;
close($configFile);

# Generate the resume file.
open(my $resumeFile    ,">",$options{'outputDirectory'}."haloMassFunctionConfigResume.xml");
print $resumeFile     $configOpener     ;
print $resumeFile     $configLikelihood ;
print $resumeFile     $configResumer    ;
print $resumeFile     $configCloser     ;
close($resumeFile);

# Generate the parameters file.
open(my $parametersFile,">",$options{'outputDirectory'}."haloMassFunctionParameters.xml"  );
print $parametersFile $parametersOpener ;
print $parametersFile $parametersCloser ;
close($parametersFile);

exit;
