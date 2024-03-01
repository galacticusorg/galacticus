#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use PDL::Math;
use XML::Simple;
use Data::Dumper;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use Galacticus::Options;
use Galacticus::Launch::Hooks;
use Galacticus::Launch::PBS;
use Galacticus::Launch::Slurm;
use Galacticus::Launch::Local;
use Galacticus::Constraints::HaloMassFunctions qw(iterate);

# Generate a halo mass function using the optimal parameters.
# Andrew Benson (22-September-2020)

# Get arguments.
my %options;
$options{'force'} = "no";
&Galacticus::Options::Parse_Options(\@ARGV,\%options);
die("no `--pipelinePath` option given")
    unless ( exists($options{'pipelinePath'}) );
die("no `--outputDirectory` option given")
    unless ( exists($options{'outputDirectory'}) );

# Ensure paths are correctly suffixed.
foreach my $path ( 'pipelinePath', 'outputDirectory' ) {
    $options{$path} .= "/"
	unless ( $options{$path} =~ m/\/$/ );
}

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Parse the simulations definition file.
my $xml = new XML::Simple();
my $simulations = $xml->XMLin(
    $options{'pipelinePath'}."haloMassFunctionSimulations.xml",
    ForceArray => [ "suite"         , "group"         , "simulation"          ],
    KeyAttr    => {  suite => "name",  group => "name",  simulation => "name" }
    );

# Specify the minimum number of particles used in fitting halo mass functions.
my $countParticlesMinimum = 100;

# Determine if models/plots should be forcibly remade.
$options{'force'} = "yes"
    if ( exists($options{'select'}) );

# Construct file names for models.
my @entries = &iterate($simulations,\%options);
foreach my $entry ( @entries ) {
    # Construct the parameter and model file names.
    my $identifier = $entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_".$entry->{'realization'}."_z".$entry->{'redshift'};
    $entry->{'fileParameters'} = $options{'outputDirectory'}."haloMassFunctionBase_".$identifier.".xml";
    $entry->{'filePrefix'    } = $options{'outputDirectory'}."haloMassFunction_"    .$identifier       ;
}

# Create sets of models that can be plotted together.
my %setsRealization;
foreach my $entry ( @entries ) {
    my $identifier = $entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_z".$entry->{'redshift'};
    push(
	@{$setsRealization{$identifier}},
	$entry
	);
}
my %identifiersHiRes;
my %setsResolution;
foreach my $entry ( @entries ) {
    next
	unless ( $entry->{'group'}->{'name'} =~ m/:hires$/ );
    (my $identifierLoRes = $entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_z".$entry->{'redshift'}) =~ s/:hires//;
    ++$identifiersHiRes{$identifierLoRes};
}
foreach my $entry ( @entries ) {
    (my $identifierLoRes = $entry->{'suite'}->{'name'}."_".$entry->{'group'}->{'name'}."_".$entry->{'simulation'}->{'name'}."_z".$entry->{'redshift'}) =~ s/:hires//;
    next
	unless ( exists($identifiersHiRes{$identifierLoRes}) );
    if ( $entry->{'group'}->{'name'} =~ m/:hires$/ ) {
	push(@{$setsResolution{$identifierLoRes}->{'hiRes'}},$entry);
    } else {
    	push(@{$setsResolution{$identifierLoRes}->{'loRes'}},$entry);
    }
}

# Generate jobs to create the models.
my @jobs =
    map
    {
	-e $_->{'filePrefix'} && $options{'force'} eq "no"
	?
        ()
	:
	{
	    command    => $ENV{'GALACTICUS_EXEC_PATH'}."/Galacticus.exe ".$_->{'fileParameters'},
	    launchFile => $_->{'filePrefix'}.".sh" ,
	    logFile    => $_->{'filePrefix'}.".log",
	    label      => "haloMassFunction_".$_->{'suite'}->{'name'}."_".$_->{'group'}->{'name'}."_".$_->{'simulation'}->{'name'}."_".$_->{'realization'}."_z".$_->{'redshift'},
	    ppn        => 1,
	    nodes      => 1,
	    mpi        => "no"
	}
    }
    @entries;
&{$Galacticus::Launch::Hooks::moduleHooks{$queueManager->{'manager'}}->{'jobArrayLaunch'}}(\%options,@jobs)
    if ( scalar(@jobs) > 0 );

# Create plots of each individual model.
foreach my $entry ( @entries ) {
    next
	if ( -e $entry->{'filePrefix'}.".pdf" && $options{'force'} eq "no" );
    # Read the model halo mass function.
    my $dataModel;
    my $fileModel = new PDL::IO::HDF5($entry->{'filePrefix'}.":MPI0000.hdf5");
    my $output    = $fileModel->group('Outputs/Output1');
    $dataModel->{$_} = $output->dataset($_)->get()
	foreach ( 'haloMass', 'haloMassFunctionLnMBinAveraged' );
    $dataModel->{'nonZero'} = which($dataModel->{'haloMassFunctionLnMBinAveraged'} > 0);
    # Determine the minimum halo mass that was used in the fit.
    my $massHaloMinimum = $countParticlesMinimum*$entry->{'group'}->{'massParticle'};
    # Read the target data.
    my $dataTarget;
    my $fileTarget       = new PDL::IO::HDF5($entry->{'fileTargetData'});
    my $targetSimulation = $fileTarget->group('simulation0001');
    $dataTarget->{$_} = $targetSimulation->dataset($_)->get()
	foreach ( 'mass', 'massFunction', 'count' );
    # Construct the target dataset errors.
    $dataTarget->{'massFunctionError'} = $dataTarget->{'massFunction'}->copy();
    $dataTarget->{'nonZero'          } = which(($dataTarget->{'count'} > 0) & ($dataTarget->{'mass'} >= $massHaloMinimum));
    $dataTarget->{'nonFit'           } = which(($dataTarget->{'count'} > 0) & ($dataTarget->{'mass'} <  $massHaloMinimum));
    $dataTarget->{'massFunctionError'}->($dataTarget->{'nonZero'}) /= $dataTarget->{'count'}->double()->($dataTarget->{'nonZero'})->sqrt();
    # Create the plot.
    my $plot;
    my $gnuPlot;
    my $plotFileTeX = $entry->{'filePrefix'}.".tex";
    my $xMinimum = sprintf("%7.1e",10.0**floor($dataTarget->{'mass'        }->($dataTarget->{'nonZero'})->minimum()->log10()));
    my $xMaximum = sprintf("%7.1e",10.0**ceil ($dataTarget->{'mass'        }->($dataTarget->{'nonZero'})->maximum()->log10()));
    my $yMinimum = sprintf("%7.1e",10.0**floor($dataTarget->{'massFunction'}->($dataTarget->{'nonZero'})->minimum()->log10()));
    my $yMaximum = sprintf("%7.1e",10.0**ceil ($dataTarget->{'massFunction'}->($dataTarget->{'nonZero'})->maximum()->log10()));
    open($gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set title offset 0.0,-0.5 '".$entry->{'suite'}->{'name'}." ".$entry->{'group'}->{'name'}." ".$entry->{'simulation'}->{'name'}." ".$entry->{'realization'}." \$z=".$entry->{'redshift'}."\$'\n";
    print $gnuPlot "set xlabel '\$ M \$ [\$\\mathrm{M}_\\odot\$]'\n";
    print $gnuPlot "set ylabel '\$ \\mathrm{d} n / \\mathrm{d} \\log M \$ [Mpc\$^{-3}\$]\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
    print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
    print $gnuPlot "set pointsize 1.0\n";
    my $lineX = pdl [ $massHaloMinimum, $massHaloMinimum ];
    my $lineY = pdl [ $yMinimum       , $yMaximum        ];
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot                                                        ,
	$lineX                                                        ,
	$lineY                                                        ,
	style        => "line"                                        ,
	weight       => [2,1]                                         ,
	linePattern  => 2                                             ,
	color        => $GnuPlot::PrettyPlots::colorPairs{'blackGray'}
	);   
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot                                                                   ,
	$dataModel->{'haloMass'                      }->($dataModel->{'nonZero'}),
	$dataModel->{'haloMassFunctionLnMBinAveraged'}->($dataModel->{'nonZero'}),
	style        => "line"                                                   ,
	weight       => [2,1]                                                    ,
	color        => $GnuPlot::PrettyPlots::colorPairs{'redYellow'}
	);   
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot                                                        ,
	$dataTarget->{'mass'        }->($dataTarget->{'nonFit'}),
	$dataTarget->{'massFunction'}->($dataTarget->{'nonFit'}),
	errorDown    => $dataTarget->{'massFunctionError'}->($dataTarget->{'nonZero'}),
	errorUp      => $dataTarget->{'massFunctionError'}->($dataTarget->{'nonZero'}),
	style        => "point"                                        ,
	symbol       => [6,7]                                         ,
	weight       => [2,1]                                         ,
	pointSize    => 0.25                                             ,
	transparency => 0.9,
	color        => $GnuPlot::PrettyPlots::colorPairs{'mediumSeaGreen'}
	);   
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	\$plot                                                        ,
	$dataTarget->{'mass'        }->($dataTarget->{'nonZero'}),
	$dataTarget->{'massFunction'}->($dataTarget->{'nonZero'}),
	errorDown    => $dataTarget->{'massFunctionError'}->($dataTarget->{'nonZero'}),
	errorUp      => $dataTarget->{'massFunctionError'}->($dataTarget->{'nonZero'}),
	style        => "point"                                        ,
	symbol       => [6,7]                                         ,
	weight       => [2,1]                                         ,
	pointSize    => 0.25                                             ,
	color        => $GnuPlot::PrettyPlots::colorPairs{'mediumSeaGreen'}
	);   
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
}

exit;

exit 0;
