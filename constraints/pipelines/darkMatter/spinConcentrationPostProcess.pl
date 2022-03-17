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

# Generate spin and concentration distribution functions using the optimal parameters.
# Andrew Benson (27-October-2021)

# Get arguments.
die('Usage: spinConcentrationPostProcess.pl <outputDirectory>')
    unless ( scalar(@ARGV) >= 1 );
my $outputDirectory = $ARGV[0];
my %options;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Define simulations to process.
my @simulations =
(
 {
     label               => "VSMDPL",
     description         => "non-backsplash z=0 halos from the VSMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/vsmdpl/",
     color               => "cornflowerBlue"
 },
 {
     label               => "SMDPL",
     description         => "non-backsplash z=0 halos from the SMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/smdpl/",
     color               => "mediumSeaGreen"
 },
 {
     label               => "MDPL2",
     description         => "non-backsplash z=0 halos from the MDPL2 simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/mdpl2/",
     color               => "salmon"
 },
 {
     label               => "BigMDPL",
     description         => "non-backsplash z=0 halos from the BigMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/bigmdpl/",
     color               => "indianRed"
 },
 {
     label               => "HugeMDPL",
     description         => "non-backsplash z=0 halos from the HugeMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/hugemdpl/",
     color               => "midnightBlue"
 }
);

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Get an XML object.
my $xml = new XML::Simple();

# Iterate over simulations to compute the distributions.
my @jobs;
foreach my $simulation ( @simulations ) {
    # Read the parameter file.
    my $parameters = $xml->XMLin($outputDirectory."/spinConcentrationBase".$simulation->{'label'}.".xml");
    # Modify parameters.
    $parameters->{'outputFileName'}->{'value'} =  $outputDirectory."/spinConcentration".$simulation->{'label'}.".hdf5";
    # Write the parameters.
    open(my $parameterFile,">".$outputDirectory."/spinConcentrationBase".$simulation->{'label'}.".xml");
    print $parameterFile $xml->XMLout($parameters, RootName => "parameters");
    close($parameterFile);
    # Generate a job.
    my $job;
    $job->{'command'   } =
	"Galacticus.exe ".$outputDirectory."/spinConcentrationBase".$simulation->{'label'}.".xml" ;
    $job->{'launchFile'} = $outputDirectory."/spinConcentrations_".$simulation->{'label'}.".sh" ;
    $job->{'logFile'   } = $outputDirectory."/spinConcentrations_".$simulation->{'label'}.".log";
    $job->{'label'     } =                   "spinConcentrations_".$simulation->{'label'}       ;
    $job->{'ppn'       } = 16;
    $job->{'nodes'     } = 1;
    $job->{'mpi'       } = "yes";
    push(@jobs,$job)
	unless ( -e $outputDirectory."/spinConcentration".$simulation->{'label'}.":MPI0000.hdf5" );
}
&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs)
    if ( scalar(@jobs) > 0 );

# Read all  data.
foreach my $simulation ( @simulations ) {
    # Read the Galacticus model results.
    my $model                                           = new PDL::IO::HDF5($outputDirectory."/spinConcentration".$simulation->{'label'}.":MPI0000.hdf5" );
    my $analyses                                        = $model   ->group ('analyses');
    my @analysisGroupNames                              = $analyses->groups(          );
    foreach my $analysisGroupName ( @analysisGroupNames ) {
	my $analysisGroup = $analyses->group($analysisGroupName);
	if ( $analysisGroupName =~ m/^spinDistribution/ ) {
	    $simulation->{'model'}->{$analysisGroupName}->{$_} = $analysisGroup->dataset($_)->get()
		foreach ( 'spin', 'spinDistributionFunction', 'spinDistributionFunctionTarget', 'spinDistributionFunctionCovariance', 'spinDistributionFunctionCovarianceTarget' );
	} else {
	    $simulation->{'model'}->{$analysisGroupName}->{$_} = $analysisGroup->dataset($_)->get()
		foreach ( 'concentration', 'concentrationFunction', 'concentrationFunctionTarget', 'concentrationFunctionCovariance', 'concentrationFunctionCovarianceTarget' );
	}
	($simulation->{'model'}->{$analysisGroupName}->{$_}) = $analysisGroup->attrGet($_)
	    foreach ( 'description', 'targetLabel', 'logLikelihood' );
    }
}

# Begin creating the plots.
foreach my $simulation ( @simulations ) {
    # Iterate over analyses.
    foreach my $analysisName ( keys(%{$simulation->{'model'}}) ) {
	my $analysis = $simulation->{'model'}->{$analysisName};
	## Spin distribution function.
	if ( $analysisName =~ m/^spinDistribution/ ) {
	    my $plot;
	    my $gnuPlot;
	    my $plotFileTeX = $outputDirectory."/".$analysisName.".tex";
	    open($gnuPlot,"|gnuplot");
	    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
	    print $gnuPlot "set output '".$plotFileTeX."'\n";
	    print $gnuPlot "set xlabel '\$ \\lambda \$ '\n";
	    print $gnuPlot "set ylabel '\$ \\mathrm{d} f / \\mathrm{d} \\log \\lambda \$ \n";
	    print $gnuPlot "set lmargin screen 0.15\n";
	    print $gnuPlot "set rmargin screen 0.95\n";
	    print $gnuPlot "set bmargin screen 0.15\n";
	    print $gnuPlot "set tmargin screen 0.95\n";
	    print $gnuPlot "set title offset 0.0,-0.75 '\\scriptsize ".$analysis->{'description'}."'\n";
	    print $gnuPlot "set key spacing 1.2\n";
	    print $gnuPlot "set key at screen 0.7,0.89\n";
	    print $gnuPlot "set logscale xy\n";
	    print $gnuPlot "set mxtics 10\n";
	    print $gnuPlot "set mytics 10\n";
	    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	    print $gnuPlot "set xrange [1.0e-3:6.0e-1]\n";
	    print $gnuPlot "set yrange [1.0e-4:1.0e+0]\n";
	    print $gnuPlot "set pointsize 1.0\n";
	    # Plot target and model results.
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                        ,
		$analysis->{'spin'                          }                                                 ,
		$analysis->{'spinDistributionFunctionTarget'}                                                 ,
		errorUp      => $analysis->{'spinDistributionFunctionCovarianceTarget'}->diagonal(0,1)->sqrt(),
		errorDown    => $analysis->{'spinDistributionFunctionCovarianceTarget'}->diagonal(0,1)->sqrt(),
		style        => "point"                                                                       ,
		weight       => [2,1]                                                                         ,
		symbol       => [6,6]                                                                         ,
		pointSize    => 1.0                                                                           ,
		color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}                     ,
		title        => $analysis->{'targetLabel'}
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                  ,
		$analysis->{'spin'                    }                                                 ,
		$analysis->{'spinDistributionFunction'}                                                 ,
		errorUp      => $analysis->{'spinDistributionFunctionCovariance'}->diagonal(0,1)->sqrt(),
		errorDown    => $analysis->{'spinDistributionFunctionCovariance'}->diagonal(0,1)->sqrt(),
		style        => "point"                                                                 ,
		weight       => [2,1]                                                                   ,
		symbol       => [6,7]                                                                   ,
		pointSize    => 0.5                                                                     ,
		color        => $GnuPlot::PrettyPlots::colorPairs{'redYellow'}                          ,
		title        => "Galacticus"
		);	    
	    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	    close($gnuPlot);
	    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
	} else {
	    ## Concentration distribution functions.
	    my $plot;
	    my $gnuPlot;
	    my $plotFileTeX = $outputDirectory."/".$analysisName.".tex";
	    open($gnuPlot,"|gnuplot");
	    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
	    print $gnuPlot "set output '".$plotFileTeX."'\n";
	    print $gnuPlot "set xlabel '\$ c_\\mathrm{vir} \$ '\n";
	    print $gnuPlot "set ylabel '\$ \\mathrm{d} f / \\mathrm{d} \\log c_\\mathrm{vir} \$ \n";
	    print $gnuPlot "set lmargin screen 0.15\n";
	    print $gnuPlot "set rmargin screen 0.95\n";
	    print $gnuPlot "set bmargin screen 0.15\n";
	    print $gnuPlot "set tmargin screen 0.95\n";
	    print $gnuPlot "set title offset 0.0,-0.75 '\\scriptsize ".$analysis->{'description'}."'\n";
	    print $gnuPlot "set key spacing 1.2\n";
	    print $gnuPlot "set key at screen 0.7,0.89\n";
	    print $gnuPlot "set logscale x\n";
	    print $gnuPlot "set mxtics 10\n";
	    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	    print $gnuPlot "set xrange [1.0e0:1.0e2]\n";
	    print $gnuPlot "set yrange [-0.05:2.0e+0]\n";
	    print $gnuPlot "set pointsize 1.0\n";
	    # Plot target and model results.
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                                     ,
		$analysis->{'concentration'              }                                                 ,
		$analysis->{'concentrationFunctionTarget'}                                                 ,
		errorUp      => $analysis->{'concentrationFunctionCovarianceTarget'}->diagonal(0,1)->sqrt(),
		errorDown    => $analysis->{'concentrationFunctionCovarianceTarget'}->diagonal(0,1)->sqrt(),
		style        => "point"                                                                    ,
		weight       => [2,1]                                                                      ,
		symbol       => [6,6]                                                                      ,
		pointSize    => 1.0                                                                        ,
		color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}                  ,
		title        => $analysis->{'targetLabel'}
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset(
		\$plot                                                                               ,
		$analysis->{'concentration'        }                                                 ,
		$analysis->{'concentrationFunction'}                                                 ,
		errorUp      => $analysis->{'concentrationFunctionCovariance'}->diagonal(0,1)->sqrt(),
		errorDown    => $analysis->{'concentrationFunctionCovariance'}->diagonal(0,1)->sqrt(),
		style        => "point"                                                              ,
		weight       => [2,1]                                                                ,
		symbol       => [6,7]                                                                ,
		pointSize    => 0.5                                                                  ,
		color        => $GnuPlot::PrettyPlots::colorPairs{'redYellow'}                       ,
		title        => "Galacticus"
		);	  
	    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	    close($gnuPlot);
	    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
	}
    }
}

exit 0;
