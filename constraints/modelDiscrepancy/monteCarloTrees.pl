#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
use XML::Simple;
require Galacticus::Options;
require Galacticus::Constraints::Parameters;
require Galacticus::Constraints::DiscrepancyModels;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Run calculations to determine the model discrepancy arising from the use of Monte Carlo merger trees.
# Andrew Benson (16-November-2012)

# Get arguments.
die("Usage: monteCarloTrees.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make              => "yes",
     plot              => "no" ,
     waitSleepDuration => 10   ,
     monotonizeGrowth  => "no"
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Define constants.
my $massSolar = pdl 1.93392e30;

# Parse the constraint config file.
my $xml    = new XML::Simple;
my $config = $xml->XMLin($configFile, KeyAttr => 0);

# Validate the config file.
die("monteCarloTrees.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("monteCarloTrees.pl: compilation must be specified in config file"   ) unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("monteCarloTrees.pl: baseParameters must be specified in config file") unless ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Validate options.
die("monteCarloTrees.pl: a tree directory must be specified")
    unless ( exists($arguments{'treeDirectory'}) );

# Determine the scratch and work directories.
my $workDirectory    =        $config->{'likelihood'}->{'workDirectory'   } ;
my $scratchDirectory = exists($config->{'likelihood'}->{'scratchDirectory'})
                       ?
                              $config->{'likelihood'}->{'scratchDirectory'}
                       :
                              $config->{'likelihood'}->{'workDirectory'   } ;

# Create the work and scratch directories.
system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'});

# Determine base parameters to use.
my $baseParameters = exists($arguments{'baseParameters'}) 
	             ?
	                    $arguments{'baseParameters'}
                     :
	                    $config->{'likelihood'}->{'baseParameters'};

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("monteCarloTrees.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'likelihood'}->{'compilation'},$baseParameters);
my @constraints = @{$constraintsRef};

# Extract all existing analyses.
my $analyses = $parameters->{'mergerTreeAnalyses'}->{'value'};
   
# Ensure we have redshift zero included in the outputs.
my $redshifts = $parameters->{'outputRedshifts'}->{'value'};
$redshifts .= " 0.0"
    unless ( grep {$_ == 0.0}  split(" ",$redshifts) );

# Define parameters for N-body tree calculation.
my $models =
{
    default =>
    {
	label      => "monteCarlo",
	parameters =>
	    [
	     {
		 name  => "outputRedshifts",
		 value => $redshifts
	     },
	     {
		 name  => "mergerTreeOperatorMethod",
		 value => "sequence"
	     },
	     {
		 name  => "mergerTreeOperatorMethod->\@mergerTreeOperatorMethod",
		 value => "regridTimes"
	     },
	     {
		 name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->dumpTrees",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->snapTolerance",
		 value => 0.0
	     },
	     {
		 name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->snapshotSpacing",
		 value => "list"
	     },
	     {
		 name  => "mergerTreeImportGalacticusReweightTrees",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeConstructMethod",
		 value => "build"
	     },
	     {
		 name  => "mergerTreeBuildTreesHaloMassDistribution",
		 value => "read"
	     },
	     {
		 name  => "mergerTreeBuildTreeMassesFile",
		 value => $workDirectory."/modelDiscrepancy/monteCarloTrees/treeRootMasses.hdf5"
	     },	     
	     {
		 name  => "mergerTreeBuilderMethod->redshiftMaximum",
		 value => 1.0e30
	     }
	    ]
    },
    alternate =>
    {
	label      => "nBody",
	parameters =>
	    [
	     {
		 name  => "outputRedshifts",
		 value => $redshifts
	     },
	     {
		 name  => "mergerTreeConstructMethod",
		 value => "read"
	     },
	     {
		 name  => "mergerTreeReadPresetMergerTimes",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetMergerNodes",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetSubhaloMasses",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetSubhaloIndices",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetPositions",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetScaleRadii",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetSpins",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetSpins",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetOrbits",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetParticleCounts",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetVelocityMaxima",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadPresetVelocityDispersions",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadAllowBranchJumps",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeReadAllowSubhaloPromotions",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeImportGalacticusReweightTrees",
		 value => "false"
	     },
	     {
		 name  => "allTreesExistAtFinalTime",
		 value => "false"
	     },
	     {
		 name  => "mergerTreeOperatorMethod",
		 value => "sequence"
	     },
	     {
		 name  => "mergerTreeOperatorMethod->\@mergerTreeOperatorMethod",
		 value => "monotonizeMassGrowth"
	     },
	     {
		 name  => "mergerTreeOperatorMethod->\@mergerTreeOperatorMethod",
		 value => "outputRootMasses"
	     },
	     {
		 name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->redshift",
		 value => 0.0
	     },
	     {
		 name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->alwaysIsolatedHalosOnly",
		 value => "true"
	     },
	     {
		 name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->fileName",
		 value => $workDirectory."/modelDiscrepancy/monteCarloTrees/treeRootMasses.hdf5"
	     },
	     {
		 name  => "treeNodeMethodMergingStatistics",
		 value => "standard"
	     },
	     {
		 name  => "hierarchyLevelResetFactor",
		 value => 2.0
	     },
	     {
		 name  => "mergerTreeReadFixedThreadAssignment",
		 value => "false"
	     }
	    ]
    }
};

# Iterate over available merger trees.
my @treeFiles;
my $snapshotRedshifts = pdl [];
my $haloMassMinimum   = pdl 1.0e60;
opendir(my $treeDirectory,$arguments{'treeDirectory'});
while ( my $fileName = readdir($treeDirectory) ) {
    next
	unless ( $fileName =~ m/\.hdf5$/ );    
    print "Processing tree file: ".$fileName."\n";
    # Extract cosmological parameters if necessary.
    foreach my $model ( "default", "alternate" ) {
	my $parameters = $models->{$model}->{'parameters'};
	unless ( grep {$_ eq "cosmologyParametersMethod->OmegaMatter"} @{$parameters} ) {
	    my $treeFile      = new PDL::IO::HDF5($arguments{'treeDirectory'}."/".$fileName);
	    my $treeCosmology = $treeFile->group('cosmology');
	    my %mapping =
		(
		 'HubbleParam'        => 'cosmologyParametersMethod->HubbleConstant',
		 'OmegaBaryon'        => 'cosmologyParametersMethod->OmegaBaryon'   ,
		 'OmegaMatter'        => 'cosmologyParametersMethod->OmegaMatter'   ,
		 'powerSpectrumIndex' => 'powerSpectrumPrimordialMethod->index'     ,
		 'sigma_8'            => 'cosmologicalMassVarianceMethod->sigma_8'
		);
	    foreach ( keys(%mapping) ) {
		(my $value) = $treeCosmology->attrGet($_);
		$value *= 100.0
		    if ( $_ eq "HubbleParam" );
		push(
		    @{$parameters},
		    {
			name  => $mapping{$_},
			value => $value->sclr()
		    }
		    );
		push(
		    @{$parameters},
		    {
			name  => "cosmologyParametersMethod->OmegaDarkEnergy",
			value => 1.0-$value->sclr()
		    }		
		    )
		    if ( $_ eq "OmegaMatter" );
	    }
	    (my $transferFunction) = $treeCosmology->attrGet('transferFunction');
	    if ( $transferFunction eq "CAMB" || $transferFunction eq "CMBFast" ) {
		my ($sigma8) = grep {$_->{'name'} eq "cosmologicalMassVarianceMethod->sigma_8"} @{$parameters};
		push(
		    @{$parameters},
		    {
			name  => "transferFunctionMethod",
			value => "CAMB"
		    },	
		    {
			name  => "cosmologicalMassVarianceMethod->tolerance",
			value => 3.0e-4
		    },		
		    {
			name  => "cosmologicalMassVarianceMethod->toleranceTopHat",
			value => 1.0e-4
		    },
		    {
			name  => "powerSpectrumNonlinearMethod",
			value => "peacockDodds1996"
		    },
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod",
			value => "standard"
		    },
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod->transferFunctionMethod",
			value => "eisensteinHu1999"
		    },
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod->transferFunctionMethod->neutrinoMassSummed", 
			value => 0.0
		    },
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod->transferFunctionMethod->neutrinoNumberEffective",
			value => 3.04
		    },		    
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod->powerSpectrumPrimordialTransferredMethod",
			value => "simple"
		    },
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod->cosmologicalMassVarianceMethod", 
			value => "filteredPower"
		    },
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod->cosmologicalMassVarianceMethod->sigma_8", 
			value => $sigma8->{'value'}
		    },
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod->cosmologicalMassVarianceMethod->toleranceTopHat", 
			value => 0.0001
		    },
		    {
			name  => "powerSpectrumNonlinearMethod->powerSpectrumMethod->cosmologicalMassVarianceMethod->tolerance",
			value => 0.0003
		    }		    
		    );
	    } else {
		die('monteCarloTrees.pl: unrecognized transfer function');
	    }
	}
    }	
    # Extract snapshot redshifts.
    {
	my $treeFile       = new PDL::IO::HDF5($arguments{'treeDirectory'}."/".$fileName);
	my $forestHalos    = $treeFile   ->group  ('forestHalos')       ;
	my $haloRedshifts  = $forestHalos->dataset('redshift'   )->get();
	$snapshotRedshifts = $snapshotRedshifts->append($haloRedshifts->uniq());
    }
    # Extract virial overdensity definition information.
    {
	my $treeFile          = new PDL::IO::HDF5($arguments{'treeDirectory'}."/".$fileName);
	my $groupFinder       = $treeFile   ->group  ('groupFinder');
	(my $groupFinderCode) = $groupFinder->attrGet('code'       );
	if ( $groupFinderCode eq "ROCKSTAR" || $groupFinderCode eq "SUBFIND" ) {
	    # RockStar and SubFind use an FoF finder. Find the linking length and set a percolation theory-based virial density
	    # contrast threshold.
	    (my $linkingLength) = $groupFinder->attrGet('linkingLength');
	    foreach ( "default", "alternate" ) {
		my @virialDensityParameters =
		    (
		     {
			 name  => "virialDensityContrastMethod",
			 value => "percolation"
		     },
		     {
			 name  => "virialDensityContrastPercolationLinkingLength",
			 value => $linkingLength->sclr()
		     }
		    );
		push(
		    @{$models->{$_}->{'parameters'}},
		    @virialDensityParameters
		    );
	    }
	}
    }
    # Extract minimum halo mass.
    {
	my  $treeFile         = new PDL::IO::HDF5($arguments{'treeDirectory'}."/".$fileName);
	my  $forestHalos      = $treeFile   ->group  ('forestHalos')       ;
	my  $units            = $treeFile   ->group  ('units'      )       ;
	my  $cosmology        = $treeFile   ->group  ('cosmology'  )       ;
	my  $haloMasses       = $forestHalos->dataset('nodeMass'   )->get();
	my  $haloRedshifts    = $forestHalos->dataset('redshift'   )->get();
	(my $hubbleParameter) = $cosmology  ->attrGet('HubbleParam')       ; 
	(my $unitsSI, my $unitsExpScaleFactor, my $unitsExpHubble)
	    = $units->attrGet('massUnitsInSI','massScaleFactorExponent','massHubbleExponent');
	my $thisHaloMassMinimum = 
	    minimum($haloMasses)
	    *($unitsSI/$massSolar)
	    *($hubbleParameter**$unitsExpHubble)
	    *((1.0/(1.0+$haloRedshifts))**$unitsExpScaleFactor);
	$haloMassMinimum = $thisHaloMassMinimum->((0))
	    if ( $thisHaloMassMinimum->((0)) < $haloMassMinimum );
    }
    # Set the tree file name.
    push(@treeFiles,$fileName);
}
closedir($treeDirectory);
# Set baryon pruning if requested.
if ( exists($arguments{'pruneBaryons'}) ) {
    push(
        @{$models->{$_}->{'parameters'}},
	{
	    name  => "mergerTreePruneBaryons",
	    value => $arguments{'pruneBaryons'}
	}
	)
	foreach ( "default", "alternate" );
}
# Set regrid operator for Monte Carlo trees.
$snapshotRedshifts = $snapshotRedshifts->uniq();
push(
    @{$models->{'default'}->{'parameters'}},
    {
	name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->snapshotRedshifts",
	value => join(" ",$snapshotRedshifts->list())
    },
    {
	name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->regridCount",
	value => nelem($snapshotRedshifts)
    },
    );
# Strip mass growth monotonizer if not required.
@{$models->{'alternate'}->{'parameters'}} = map {$_->{'value'} eq "monotonizeMassGrowth" ? () : $_} @{$models->{'alternate'}->{'parameters'}}
    unless ( $arguments{'monotonizeGrowth'} eq "yes" );
# Set options for growth on non-new accretion.
if ( exists($arguments{'accreteNewGrowthOnly'}) ) {
    push(
	@{$models->{'alternate'}->{'parameters'}},
	{
	    name  => "treeNodeMethodBasic",
	    value => "extendedTracking"
	},
	{
	    name  => "accreteNewGrowthOnly",
	    value => $arguments{'accreteNewGrowthOnly'} eq "yes" ? "true" : "false"
	}
	);
}
# Add halo mass function analyses.
push(
        @{$models->{$_}->{'parameters'}},
    {
	name  => "analysisHaloMassFunctionsMassMinimum",
	value => sclr($haloMassMinimum/2.0)
    }
    )
    foreach ( "default", "alternate" );
my $step = int(nelem($snapshotRedshifts)/5);
for (my $i=0;$i<nelem($snapshotRedshifts);$i+=$step) {
    $analyses .= " haloMassFunctionZ".$snapshotRedshifts->(($i));
}
push(
    @{$models->{'default'}->{'parameters'}},
    {
	name  => "mergerTreeAnalyses",
	value => $analyses
    }
    );
for (my $i=0;$i<nelem($snapshotRedshifts);$i+=$step) {
    $analyses .= " isolatedHaloMassFunctionZ".$snapshotRedshifts->(($i));
}
push(
    @{$models->{'alternate'}->{'parameters'}},
    {
	name  => "mergerTreeAnalyses",
	value => $analyses
    }
    );

# Set mass resolution for Monte Carlo trees case.
my $massResolution = exists($arguments{'massResolution'}) ? $arguments{'massResolution'} : $haloMassMinimum->sclr();
push(
    @{$models->{'default'}->{'parameters'}},
    {
	name  => "mergerTreeMassResolutionMethod",
	value => "fixed"
    },
    {
	name  => "mergerTreeMassResolutionMethod->massResolution",
	value => $massResolution
    },
    );

# Set tree file names.
push(
    @{$models->{'alternate'}->{'parameters'}},
    {
	name  => "mergerTreeReadFileName",
	value => join(" ",map {$arguments{'treeDirectory'}."/".$_} @treeFiles)
    },
    );

# Add a mass resolution if provided.
if ( exists($arguments{'massThreshold'}) ) {
    push(
	@{$models->{$_}->{'parameters'}},
	{
	    name  => "mergerTreeOperatorMethod->\@mergerTreeOperatorMethod",
	    value => "pruneByMass"
	},
	{
	    name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->massThreshold",
	    value => $arguments{'massThreshold'}
	},
	{
	    name  => "mergerTreeOperatorMethod->mergerTreeOperatorMethod->preservePrimaryProgenitor",
	    value => "true"
	}	
	)
	foreach ( "default", "alternate" );
}

# Run the models.
my @modelRuns =
    (
     {
     	 runOnly => "alternate",
     	 analyze => "no"
     },
     {
	 runOnly => "default",
	 analyze => "yes"
     }
    );
foreach my $model ( @modelRuns ) {
    $arguments{$_} = $model->{$_}
        foreach ( keys(%{$model}) );
    &DiscrepancyModels::RunModels(
	    "monteCarloTrees"                              ,
	    "the use of Monte Carlo generated merger trees",
	    $configFile                                    ,
	    \%arguments                                    ,
	    $models
	);
}

# Make plots of the halo mass functions.
for (my $i=0;$i<nelem($snapshotRedshifts);$i+=$step) {
    my ($gnuPlot, $plot);
    my $z = $snapshotRedshifts->(($i));
    (my $pz = $z) =~ s/\./_/g;
    # Open a pipe to GnuPlot.
    my $plotFileTeX = $workDirectory."/modelDiscrepancy/monteCarloTrees/haloMassFunctionZ".$pz.".tex";
    open($gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal cairolatex standalone pdf transparent color lw 2 size 4in,2.8in font 'cmr,m,n'\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.96\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.94\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.4,0.2\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [1.0e8:1.0e15]\n";
    print $gnuPlot "set yrange [1.0e-9:1.0e0]\n";
    print $gnuPlot "set xlabel '\$M_\\mathrm{halo}/\\mathrm{M}_\\odot\$'\n";
    print $gnuPlot "set ylabel '\$\\mathrm{d}n/\\mathrm{d}\\log M_\\mathrm{halo}/\\hbox{Mpc}^{-3}\$'\n";
    my $j = -1;
    foreach my $modelName ( "nBody", "monteCarlo" ) {
	my $model = new PDL::IO::HDF5($workDirectory."/modelDiscrepancy/monteCarloTrees/".$modelName."/galacticus.hdf5");
	my $analysis = $model->group('analysis');
	my @types = $modelName eq "nBody" ? ( "halo", "isolatedHalo" ) : ( "halo" );
	foreach my $type ( @types ) {
	    ++$j;
	    my $hmf      = $analysis->group($type.'MassFunctionZ'.$z);
	    my $mass     = $hmf->dataset('mass')->get();
	    my $massFunction = $hmf->dataset('massFunction')->get();
	    my $massFunctionCovariance = $hmf->dataset('massFunctionCovariance')->get();
	    &PrettyPlots::Prepare_Dataset(
		\$plot,
		$mass,
		$massFunction,
		errorUp   => sqrt($massFunctionCovariance->diagonal(0,1)),
		errorDown => sqrt($massFunctionCovariance->diagonal(0,1)),
		style     => "point",
		symbol    => [6,7], 
		weight    => [5,3],
		pointSize => 0.1,
		color     => $PrettyPlots::colorPairs{$PrettyPlots::colorPairSequences{'contrast'}[$j]},
		title     => $modelName." : ".$type
		);
	}
    }
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($plotFileTeX,margin => 1);
}

exit;
