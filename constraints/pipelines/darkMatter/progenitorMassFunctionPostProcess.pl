#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use XML::Simple;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Galacticus::Options;
use Data::Dumper;

# Generate a progenitor mass function using the optimal parameters.
# Andrew Benson (22-September-2020)

# Get arguments.
die('Usage: progenitorMassFunctionPostProcess.pl <outputDirectory>')
    unless ( scalar(@ARGV) >= 1 );
my $outputDirectory = $ARGV[0];
my %options;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Define simulations to process.
my @simulations =
(
 {
     label               => "VSMDPL",
     description         => "Progenitor halo mass function for non-backsplash z=0 halos from the VSMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/vsmdpl/",
     color               => "cornflowerBlue"
 },
 {
     label               => "SMDPL",
     description         => "Progenitor halo mass function for non-backsplash z=0 halos from the SMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/smdpl/",
     color               => "mediumSeaGreen"
 },
 {
     label               => "MDPL2",
     description         => "Progenitor halo mass function for non-backsplash z=0 halos from the MDPL2 simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/mdpl2/",
     color               => "salmon"
 },
 {
     label               => "BigMDPL",
     description         => "Progenitor halo mass function for non-backsplash z=0 halos from the BigMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/bigmdpl/",
     color               => "indianRed"
 },
 {
     label               => "HugeMDPL",
     description         => "Progenitor halo mass function for non-backsplash z=0 halos from the HugeMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/hugemdpl/",
     color               => "midnightBlue"
 },
 {
     label               => "Caterpillar_LX12",
     description         => "Progenitor halo mass function for non-backsplash z=0 halos from the Caterpillar LX12 simulations.",
     simulationReference => "Griffen et al.; 2016; ApJ; 818; 10",
     simulationURL       => "https://www.caterpillarproject.org/",
     color               => "midnightBlue"
 },
 {
     label               => "Caterpillar_LX13",
     description         => "Progenitor halo mass function for non-backsplash z=0 halos from the Caterpillar LX13 simulation.",
     simulationReference => "Griffen et al.; 2016; ApJ; 818; 10",
     simulationURL       => "https://www.caterpillarproject.org/",
     color               => "midnightBlue"
 },
 {
     label               => "Caterpillar_LX14",
     description         => "Progenitor halo mass function for non-backsplash z=0 halos from the Caterpillar LX14 simulation.",
     simulationReference => "Griffen et al.; 2016; ApJ; 818; 10",
     simulationURL       => "https://www.caterpillarproject.org/",
     color               => "midnightBlue"
 }
);

# Define parameter sets to use.
my @parameterSets =
    (
     {
	 # Best fit model from this pipeline.
	 label  => "default"
     },
     {
	 # The original Parkinson, Cole, & Helly (2008; http://adsabs.harvard.edu/abs/2008MNRAS.383..557P) parameters.
	 label  => "PCH",
	 G0     => +0.57,
	 gamma1 => +0.38,
	 gamma2 => -0.01,
	 gamma3 => +0.00
     }
    );

# Factor by which to increase the number of trees used.
my $treeBoostFactor = 4;

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Get an XML object.
my $xml = new XML::Simple();

# Specify the minimum number of particles used in fitting progenitor mass functions.
my $countParticlesMinimum = 3000;

# Iterate over simulations to get the model mass function.
my @jobs;
foreach my $simulation ( @simulations ) {
    foreach my $parameterSet ( @parameterSets ) {
	# Read the parameter file.
	my $parameters = $xml->XMLin($outputDirectory."/progenitorMassFunctionBase".$simulation->{'label'}.".xml");
	# Modify parameters.
	$parameters->{'galacticusOutputFileName'}->{'value'} =  $outputDirectory."/progenitorMassFunction".$simulation->{'label'}.$parameterSet->{'label'}.".hdf5";
	if      ( $parameters->{'mergerTreeBuildMassesMethod'}->{'value'} eq "sampledDistributionUniform" ) {
	    $parameters->{'mergerTreeBuildMassesMethod'}->{'treesPerDecade'  }->{'value'} *= $treeBoostFactor;
	} elsif ( $parameters->{'mergerTreeBuildMassesMethod'}->{'value'} eq "replicate"                  ) {
	    $parameters->{'mergerTreeBuildMassesMethod'}->{'replicationCount'}->{'value'} *= $treeBoostFactor;
	} else {
	    die('unknown "mergerTreeBuildMassesMethod"');
	}
	foreach my $name ( 'G0', 'gamma1', 'gamma2', 'gamma3' ) {
	    $parameters->{'mergerTreeBranchingProbabilityMethod'}->{$name}->{'value'} = $parameterSet->{$name}
	        if ( exists($parameterSet->{$name}) );
	}
	# Write the parameters.
	open(my $parameterFile,">".$outputDirectory."/progenitorMassFunctionBase".$simulation->{'label'}.$parameterSet->{'label'}.".xml");
	print $parameterFile $xml->XMLout($parameters, RootName => "parameters");
	close($parameterFile);
	# Generate a job.
	my $job;
	$job->{'command'   } =
	    "Galacticus.exe ".$outputDirectory."/progenitorMassFunctionBase".$simulation->{'label'}.$parameterSet->{'label'}.".xml" ;
	$job->{'launchFile'} = $outputDirectory."/progenitorMassFunctions_".$simulation->{'label'}.$parameterSet->{'label'}.".sh" ;
	$job->{'logFile'   } = $outputDirectory."/progenitorMassFunctions_".$simulation->{'label'}.$parameterSet->{'label'}.".log";
	$job->{'label'     } =                   "progenitorMassFunctions_".$simulation->{'label'}.$parameterSet->{'label'}       ;
	$job->{'ppn'       } = 16;
	$job->{'nodes'     } = 1;
	$job->{'mpi'       } = "yes";
	push(@jobs,$job)
	    unless ( -e $outputDirectory."/progenitorMassFunction".$simulation->{'label'}.$parameterSet->{'label'}.":MPI0000.hdf5" );
    }
}
&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs)
    if ( scalar(@jobs) > 0 );

# Read all  data.
foreach my $simulation ( @simulations ) {
    foreach my $parameterSet ( @parameterSets ) {
	# Read the Galacticus model results.
	my $model                                           = new PDL::IO::HDF5($outputDirectory."/progenitorMassFunction".$simulation->{'label'}.$parameterSet->{'label'}.":MPI0000.hdf5" );
	my $analyses                                        = $model   ->group ('analyses');
	my @analysisGroupNames                              = $analyses->groups(          );
	foreach my $analysisGroupName ( @analysisGroupNames ) {
	    my $analysisGroup = $analyses->group($analysisGroupName);
	    $simulation->{$parameterSet->{'label'}}->{'model'}->{$analysisGroupName}->{$_} = $analysisGroup->dataset($_)->get()
		foreach ( 'massRatio', 'progenitorMassFunction', 'progenitorMassFunctionTarget', 'progenitorMassFunctionCovariance', 'progenitorMassFunctionCovarianceTarget' );
	    ($simulation->{$parameterSet->{'label'}}->{'model'}->{$analysisGroupName}->{$_}) = $analysisGroup->attrGet($_)
		foreach ( 'description', 'targetLabel', 'massRatioLikelihoodMinimum', 'massRatioLikelihoodMaximum', 'logLikelihood' );
	}
    }
}

# Begin creating the plots.
{
    # Iterate over simulations.
    foreach my $simulation ( @simulations ) {
	# Iterate over analyses.
	foreach my $analysisName ( keys(%{$simulation->{'default'}->{'model'}}) ) {
	    my $analysisDefault = $simulation->{'default'}->{'model'}->{$analysisName};
	    ## Progenitor mass function.
	    my $plot;
	    my $gnuPlot;
	    my $plotFileTeX = $outputDirectory."/".$analysisName.".tex";
	    open($gnuPlot,"|gnuplot");
	    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
	    print $gnuPlot "set output '".$plotFileTeX."'\n";
	    print $gnuPlot "set xlabel '\$ x=M_\\mathrm{progenitor}/M_\\mathrm{parent} \$ '\n";
	    print $gnuPlot "set ylabel '\$ \\mathrm{d} f / \\mathrm{d} \\log x \$ \n";
	    print $gnuPlot "set lmargin screen 0.15\n";
	    print $gnuPlot "set rmargin screen 0.95\n";
	    print $gnuPlot "set bmargin screen 0.15\n";
	    print $gnuPlot "set tmargin screen 0.95\n";
	    print $gnuPlot "set title offset 0.0,-0.75 '\\scriptsize ".$analysisDefault->{'description'}."'\n";
	    print $gnuPlot "set key spacing 1.2\n";
	    print $gnuPlot "set key at screen 0.7,0.89\n";
	    print $gnuPlot "set logscale xy\n";
	    print $gnuPlot "set mxtics 10\n";
	    print $gnuPlot "set mytics 10\n";
	    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	    print $gnuPlot "set xrange [1.0e-8:1.0e1]\n";
	    print $gnuPlot "set yrange [1.0e-6:2.0e0]\n";
	    print $gnuPlot "set pointsize 1.0\n";
	    # Plot gray boxes indicate non-constrained regions.
	    my $x1 = pdl [ 1.0e-8                                                  , $analysisDefault->{'massRatioLikelihoodMinimum'}->sclr() ];
	    my $x2 = pdl [ $analysisDefault->{'massRatioLikelihoodMaximum'}->sclr(), 1.0e+1                                                   ];
	    my $y1 = pdl [ 1.0e-6                                                  , 1.0e-6                                                   ];
	    my $y2 = pdl [ 2.0e+0                                                  , 2.0e+0                                                   ];
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                        ,
		$x1                                                           ,
		$y1                                                           ,
		y2           => $y2                                           ,
		style        => "filledCurve"                                 ,
		color        => $GnuPlot::PrettyPlots::colorPairs{'blackGray'},
		transparency => 0.8
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                        ,
		$x2                                                           ,
		$y1                                                           ,
		y2           => $y2                                           ,
		style        => "filledCurve"                                 ,
		color        => $GnuPlot::PrettyPlots::colorPairs{'blackGray'},
		transparency => 0.8
		);
	    # Plot target and model results.
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                             ,
		$analysisDefault->{'massRatio'                   }                                                 ,
		$analysisDefault->{'progenitorMassFunctionTarget'}                                                 ,
		errorUp      => $analysisDefault->{'progenitorMassFunctionCovarianceTarget'}->diagonal(0,1)->sqrt(),
		errorDown    => $analysisDefault->{'progenitorMassFunctionCovarianceTarget'}->diagonal(0,1)->sqrt(),
		style        => "point"                                                                            ,
		weight       => [2,1]                                                                              ,
		symbol       => [6,6]                                                                              ,
		pointSize    => 1.0                                                                                ,
		color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}                          ,
		title        => $analysisDefault->{'targetLabel'}
		);
	    foreach my $parameterSet ( @parameterSets ) {
		my $analysis = $simulation->{$parameterSet->{'label'}}->{'model'}->{$analysisName};
		my $title     = "Galacticus";
		my $pointSize = 0.5;
		my $color     = "redYellow";
		if ( $parameterSet->{'label'} ne "default" ) {
		    $title     = $parameterSet->{'label'};
		    $pointSize = 0.1;
		    if      ( $parameterSet->{'label'} eq "PCH"        ) {
			$color = "plum";
		    } elsif ( $parameterSet->{'label'} eq "Newton2020" ) {
			$color = "lightGoldenrod";
		    }
		}
		&GnuPlot::PrettyPlots::Prepare_Dataset(
		     \$plot                                                                                ,
		     $analysis->{'massRatio'             }                                                 ,
		     $analysis->{'progenitorMassFunction'}                                                 ,
		     errorUp      => $analysis->{'progenitorMassFunctionCovariance'}->diagonal(0,1)->sqrt(),
		     errorDown    => $analysis->{'progenitorMassFunctionCovariance'}->diagonal(0,1)->sqrt(),
		     style        => "point"                                                               ,
		     weight       => [2,1]                                                                 ,
		     symbol       => [6,7]                                                                 ,
		     pointSize    => $pointSize                                                            ,
		     color        => $GnuPlot::PrettyPlots::colorPairs{$color}                             ,
		     title        => $title
		    );
	    }
	    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	    close($gnuPlot);
	    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
	}
    }
}

exit 0;
