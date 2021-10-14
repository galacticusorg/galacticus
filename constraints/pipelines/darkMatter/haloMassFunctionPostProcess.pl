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

# Generate a halo mass function using the optimal parameters.
# Andrew Benson (22-September-2020)

# Get arguments.
die('Usage: haloMassFunctionPostProcess.pl <outputDirectory>')
    unless ( scalar(@ARGV) == 1 );
my $outputDirectory = $ARGV[0];
my %options;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Define simulations to process.
my @simulations =
(
 {
     label               => "VSMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the VSMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/vsmdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 6.2e6,
     color               => "cornflowerBlue"
 },
 {
     label               => "SMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the SMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/smdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 9.63e7,
     color               => "mediumSeaGreen"
 },
 {
     label               => "MDPL2",
     description         => "Halo mass function for non-backsplash z=0 halos from the MDPL2 simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/mdpl2/",
     hubbleConstant      => 0.6777,
     massParticle        => 1.51e9,
     color               => "salmon"
 },
 {
     label               => "BigMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the BigMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/bigmdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 2.359e10,
     color               => "indianRed"
 },
 {
     label               => "HugeMDPL",
     description         => "Halo mass function for non-backsplash z=0 halos from the HugeMDPL simulation.",
     simulationReference => "Klypin, Yepes, Gottlober, Hess; 2016; MNRAS; 457; 4340",
     simulationURL       => "https://www.cosmosim.org/cms/simulations/hugemdpl/",
     hubbleConstant      => 0.6777,
     massParticle        => 7.9e10,
     color               => "midnightBlue"
 }
);

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Get an XML object.
my $xml = new XML::Simple();

# Specify the minimum number of particles used in fitting halo mass functions.
my $countParticlesMinimum = 3000;

# Iterate over simulations to get the model mass function.
my @jobs;
foreach my $simulation ( @simulations ) {

    # Convert particle mass to Solar masses.
    $simulation->{'massParticle'} /= $simulation->{'hubbleConstant'};

    # Parse the base parameters.
    my $parameters = $xml->XMLin($outputDirectory."/haloMassFunctionBase.xml");

    # Modify particle mass and output file name.
    $parameters->{'nbodyHaloMassError'      }->{'massParticle'}->{'value'} =                                       $simulation->{'massParticle'}        ;
    $parameters->{'galacticusOutputFileName'}                  ->{'value'} = $outputDirectory."/haloMassFunction_".$simulation->{'label'       }.".hdf5";

    # Write parmeter file.
    my $parameterFileName = $outputDirectory."/haloMassFunctionBase_".$simulation->{'label'}.".xml";
    open(my $outputFile,">",$parameterFileName);
    print $outputFile $xml->XMLout($parameters, RootName => "parameters");
    close($outputFile);
    
    # Generate a job.
    my $job;
    $job->{'command'   } =
	"Galacticus.exe ".$parameterFileName;
    $job->{'launchFile'} = $outputDirectory."/haloMassFunction_".$simulation->{'label'}.".sh" ;
    $job->{'logFile'   } = $outputDirectory."/haloMassFunction_".$simulation->{'label'}.".log";
    $job->{'label'     } =                   "haloMassFunction_".$simulation->{'label'}       ;
    $job->{'ppn'       } = 1;
    $job->{'nodes'     } = 1;
    $job->{'mpi'       } = "yes";
    push(@jobs,$job)
	unless ( -e $outputDirectory."/haloMassFunction_".$simulation->{'label'}.":MPI0000.hdf5" );
}
&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs)
    if ( scalar(@jobs) > 0 );

# Read all  data.
foreach my $simulation ( @simulations ) {
    # Read the resulting halo mass function.
    my $model              = new PDL::IO::HDF5($outputDirectory."/haloMassFunction_".$simulation->{'label'}.":MPI0000.hdf5");
    my $outputs                                    = $model           ->group  ('Outputs'                       )       ;
    my $output                                     = $outputs         ->group  ('Output1'                       )       ;
    $simulation->{'model'}->{'mass'              } = $output          ->dataset('haloMass'                      )->get();
    $simulation->{'model'}->{'massFunction'      } = $output          ->dataset('haloMassFunctionLnM'           )->get();
    $simulation->{'model'}->{'massFunctionBinned'} = $output          ->dataset('haloMassFunctionLnMBinAveraged')->get();
    
    # Read the target dataset.
    my $target             = new PDL::IO::HDF5($ENV{'GALACTICUS_DATA_PATH'}."/static/darkMatter/haloMassFunction_".$simulation->{'label'}."_z0.000.hdf5");
    my $targetSimulation                           = $target          ->group  ('simulation0001'                )       ;
    $simulation->{'target'}->{'mass'             } = $targetSimulation->dataset('mass'                          )->get();
    $simulation->{'target'}->{'massFunction'     } = $targetSimulation->dataset('massFunction'                  )->get();
    $simulation->{'target'}->{'count'            } = $targetSimulation->dataset('count'                         )->get();
    
    # Construct the target dataset errors.
    $simulation->{'target'}->{'massFunctionError'}                    = $simulation->{'target'}->{'massFunction'}->copy();
    my $nonZeroTarget                                                 = which($simulation->{'target'}->{'count'} > 0);
    $simulation->{'target'}->{'massFunctionError'}->($nonZeroTarget) /= $simulation->{'target'}->{'count'}->($nonZeroTarget)->sqrt();    

    # Interpolate mass function.
    ($simulation->{'model'}->{'massFunctionInterpolated'})  
	= interpolate(
	$simulation->{'target'}->{'mass'              }->log(),
	$simulation->{'model' }->{'mass'              }->log(),
	$simulation->{'model' }->{'massFunctionBinned'}->log()
	);
    $simulation->{'model'}->{'massFunctionInterpolated'}   .= exp($simulation->{'model'}->{'massFunctionInterpolated'});
}

# Begin creating the plots.
unless ( -e $outputDirectory."/haloMassFunction.pdf" ) {
    ## Halo mass function.
    my $plot;
    my $gnuPlot;
    my $plotFileTeX = $outputDirectory."/haloMassFunction.tex";
    open($gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
    print $gnuPlot "set ylabel '\$ \\mathrm{d} n / \\mathrm{d} \\log M \$ [Mpc\$^{-3}\$]\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.65,0.49\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [8.0e7:1.0e16]\n";
    print $gnuPlot "set yrange [1.0e-9:3.0]\n";
    print $gnuPlot "set pointsize 1.0\n";
    foreach my $simulation ( @simulations ) {
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	     \$plot                                                                   ,
	     $simulation->{'target'}->{'mass'        }                                ,
	     $simulation->{'target'}->{'massFunction'}                                ,
	     errorUp      => $simulation->{'target'}->{'massFunctionError'}           ,
	     errorDown    => $simulation->{'target'}->{'massFunctionError'}           ,
	     style        => "point"                                                  ,
	     weight       => [2,1]                                                    ,
	     symbol       => [6,6]                                                    ,
	     pointSize    => 0.5                                                      ,
	     color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}},
	     title        => $simulation->{'label'}." N-body"
	    );
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	     \$plot                                                                   ,
	     $simulation->{'model'}->{'mass'        }                                  ,
	     $simulation->{'model'}->{'massFunction'}                                ,
	     style        => "line"                                                   ,
	     weight       => [3,1]                                                    ,
	     color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}
	    );
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	     \$plot                                                                   ,
	     $simulation->{'model'}->{'mass'              }                           ,
	     $simulation->{'model'}->{'massFunctionBinned'}                           ,
	     style        => "point"                                                  ,
	     weight       => [2,1]                                                    ,
	     symbol       => [2,2]                                                    ,
	     pointSize    => 0.3                                                      ,
	     color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}
	    );
	my $massHaloMinimum = $simulation->{'massParticle'}*$countParticlesMinimum;
	my $xLimit          = pdl [ $massHaloMinimum, $massHaloMinimum ];
	my $yLimit          = pdl [ 1.0e-9          , 3.0              ];
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	     \$plot                                                                   ,
	     $xLimit                                                                  ,
	     $yLimit                                                                  ,
	     style        => "line"                                                   ,
	     weight       => [1,1]                                                    ,
	     linePattern  => 3                                                        ,
	     color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}
	    );
    }
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
}

unless ( -e $outputDirectory."/haloMassFunctionResidualsFractional.pdf" ) {
    ## Halo mass function residuals, fractional.
    my $plot;
    my $gnuPlot;
    my $plotFileTeX = $outputDirectory."/haloMassFunctionResidualsFractional.tex";
    open($gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
    print $gnuPlot "set ylabel '\$ [(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{model}-(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{target}]/(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{target} \$\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.65,0.49\n";
    print $gnuPlot "set logscale x\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [8.0e7:1.0e16]\n";
    print $gnuPlot "set yrange [-1.0:+1.0]\n";
    print $gnuPlot "set pointsize 1.0\n";
    foreach my $simulation ( @simulations ) {
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	     \$plot                                                                   ,
	     $simulation->{'target'}->{'mass'}                                        ,
	     (
	      +$simulation->{'model' }->{'massFunctionInterpolated'}
	      -$simulation->{'target'}->{'massFunction'            }
	     )
	     /$simulation->{'target'}->{'massFunction'}                               ,
	     style        => "point"                                                  ,
	     weight       => [2,1]                                                    ,
	     symbol       => [1,1]                                                    ,
	     pointSize    => 0.2                                                      ,
	     color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}
	    );
	my $massHaloMinimum = $simulation->{'massParticle'}*$countParticlesMinimum;
	my $xLimit          = pdl [ $massHaloMinimum, $massHaloMinimum ];
	my $yLimit          = pdl [ -1.0            , 1.0              ];
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	     \$plot                                                                   ,
	     $xLimit                                                                  ,
	     $yLimit                                                                  ,
	     style        => "line"                                                   ,
	     weight       => [1,1]                                                    ,
	     linePattern  => 3                                                        ,
	     color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}
	    );
    }
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
}

unless ( -e $outputDirectory."/haloMassFunctionResidualsNormalized.pdf" ) {
    ## Halo mass function residuals, normalized.
    my $plot;
    my $gnuPlot;
    my $plotFileTeX = $outputDirectory."/haloMassFunctionResidualsNormalized.tex";
    open($gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
    print $gnuPlot "set ylabel '\$ |(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{model}-(\\mathrm{d} n / \\mathrm{d} \\log M)_\\mathrm{target}|/\\sigma_\\mathrm{target} \$\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.65,0.49\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [8.0e7:1.0e16]\n";
    print $gnuPlot "set yrange [0.1:1.0e3]\n";
    print $gnuPlot "set pointsize 1.0\n";
    foreach my $simulation ( @simulations ) {
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	     \$plot                                                                   ,
	     $simulation->{'target'}->{'mass'}                                        ,
	     abs(
		 +$simulation->{'target'}->{'massFunction'            }
		 -$simulation->{'model' }->{'massFunctionInterpolated'}
	     )
	     /$simulation->{'target'}->{'massFunctionError'}                          ,
	     style        => "point"                                                  ,
	     weight       => [2,1]                                                    ,
	     symbol       => [1,1]                                                    ,
	     pointSize    => 0.2                                                      ,
	     color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}
	    );
	my $massHaloMinimum = $simulation->{'massParticle'}*$countParticlesMinimum;
	my $xLimit          = pdl [ $massHaloMinimum, $massHaloMinimum ];
	my $yLimit          = pdl [ 0.1             , 1.0e3            ];
	&GnuPlot::PrettyPlots::Prepare_Dataset(
	     \$plot                                                                   ,
	     $xLimit                                                                  ,
	     $yLimit                                                                  ,
	     style        => "line"                                                   ,
	     weight       => [1,1]                                                    ,
	     linePattern  => 3                                                        ,
	     color        => $GnuPlot::PrettyPlots::colorPairs{$simulation->{'color'}}
	    );
    }
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
}

exit 0;
