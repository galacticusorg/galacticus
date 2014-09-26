#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use Clone qw(clone);
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
require Galacticus::Constraints::Parameters;
require Galacticus::Constraints::DiscrepancySystematics;
require Galacticus::Options;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require Galacticus::Launch::PBS;

# Run calculations to determine the model discrepancy arising from the use of randomized halo spins.
# Andrew Benson (23-March-2014)

# Get arguments.
die("Usage: randomSpins.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make         => "yes" ,
     plot         => "no"
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Parameters::Parse_Config($configFile);

# Validate the config file.
die("randomSpins.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'workDirectory' }) );
die("randomSpins.pl: compilation must be specified in config file"   ) unless ( exists($config->{'compilation'   }) );
die("randomSpins.pl: baseParameters must be specified in config file") unless ( exists($config->{'baseParameters'}) );

# Determine the scratch and work directories.
my $workDirectory    = $config->{'workDirectory'};
my $scratchDirectory = $config->{'workDirectory'};
$scratchDirectory    = $config->{'scratchDirectory'}
    if ( exists($config->{'scratchDirectory'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'workDirectory'});

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("randomSpins.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'compilation'},$config->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Switch off thread locking.
$parameters->{'parameter'}->{'treeEvolveThreadLock'}->{'value'} = "false";

# Initialize a stack for PBS models.
my @pbsStack;

# Extract the standard reset factor for spins.
my $resetFactor = $parameters->{'parameter'}->{'randomSpinResetMassFactor'};

# Specify models to run.
my @models = 
    (
     {
	 label      => "standard",
	 parameters =>
	     [
	      # Use the standard factor for spin reset.
	      {
		  name  => "randomSpinResetMassFactor",
		  value => $resetFactor
	      }
	     ]
     },
     {
	 label      => "short",
	 parameters => 
	     [
	      # Use a short factor for spin reset.
	      {
		  name  => "randomSpinResetMassFactor",
		  value => $resetFactor/1.5
	      },
	     ]
     }
    );

# Iterate over models.
foreach my $model ( @models ) {
# Specify the output name.
    my $modelDirectory = $workDirectory."/modelDiscrepancy/randomSpins/".$model->{'label'};
    system("mkdir -p ".$modelDirectory);
    my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
    push(
	@{$model->{'parameters'}},
	{
	    name  => "galacticusOutputFileName",
	    value => $galacticusFileName
	}
	);
    my $newParameters = clone($parameters);
    $newParameters->{'parameter'}->{$_->{'name'}}->{'value'} = $_->{'value'}
        foreach ( @{$model->{'parameters'}} );
    # Adjust the number of trees to run if specified.
    $newParameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $arguments{'treesPerDecade'}
        if ( exists($arguments{'treesPerDecade'}) );
    # Set the fixed mass resolution.
    $newParameters->{'parameter'}->{'mergerTreeBuildMassResolutionScaledMinimum'}->{'value'} =
	$newParameters->{'parameter'}->{'mergerTreeBuildMassResolutionScaledMinimum'}->{'value'};
    $newParameters->{'parameter'}->{'mergerTreeBuildMassResolutionScaledMinimum'}->{'value'} =
	$arguments{'massResolutionFixed'}
           if ( exists($arguments{'massResolutionFixed'}) );
    # Run the model.
    unless ( -e $galacticusFileName ) {
	# Generate the parameter file.
	my $parameterFileName = $modelDirectory."/parameters.xml";
	&Parameters::Output($newParameters,$parameterFileName);
	# Create a job for PBS.
	my $command = "mpirun --bynode -np 1 Galacticus.exe ".$parameterFileName."\n";
	foreach my $constraint ( @constraints ) {
	    # Parse the definition file.
	    my $xml                  = new XML::Simple;
	    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
	    # Insert code to run the analysis code.
	    my $analysisCode = $constraintDefinition->{'analysis'};
	    $command .= $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5\n";
	}
	# Queue the calculation.
	my %job =
	    (
	     launchFile => $modelDirectory."/launch.pbs",
	     label      => "randomSpins".$model->{'label'},
	     logFile    => $modelDirectory."/launch.log",
	     command    => $command
	    );
	foreach ( 'ppn', 'walltime', 'memory' ) {
	    $job{$_} = $arguments{$_}
	    if ( exists($arguments{$_}) );
	}
	# Queue the calculation.
	push(
	    @pbsStack,
	    \%job
	    );
    }
}
# Send jobs to PBS.
&PBS::SubmitJobs(\%arguments,@pbsStack)
    if ( scalar(@pbsStack) > 0 );

# Iterate over constraints.
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $xml                  = new XML::Simple;
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Locate the model results.
    my $standardFileName   = $workDirectory."/modelDiscrepancy/randomSpins/standard/".$constraintDefinition->{'label'}.".hdf5";
    my $shortFileName      = $workDirectory."/modelDiscrepancy/randomSpins/short/"   .$constraintDefinition->{'label'}.".hdf5";
    # Read the results.
    my $standard           = new PDL::IO::HDF5($standardFileName);
    my $short              = new PDL::IO::HDF5($shortFileName   );
    # Extract the results.
    my $shortX             = $short   ->dataset('x'         )->get();
    my $shortY             = $short   ->dataset('y'         )->get();
    my $shortCovariance    = $standard->dataset('covariance')->get();
    my $standardY          = $standard->dataset('y'         )->get();
    my $standardCovariance = $standard->dataset('covariance')->get();
    # Apply any systematics models.
    my %systematicResults;
    foreach my $argument ( keys(%arguments) ) {
	if ( $argument =~ m/^systematic(.*)/ ) {
	    my $model = $1;
	    if ( exists($DiscrepancySystematics::models{$model}) ) {
		%{$systematicResults{$model}} =
		    &{$DiscrepancySystematics::models{$model}}(
		    \%arguments          ,
		    $constraintDefinition,
		    $shortX              ,
		    $shortY              ,
		    $shortCovariance     ,
		    $standardY           ,
		    $standardCovariance
		);
	    }
	}
    }
    # Compute the covariance.
    my $modelDiscrepancyCovarianceMultiplicative = 
	+$standardCovariance      
	*outer($shortY/$standardY**2,$shortY/$standardY**2);
    # Output the model discrepancy to file.
    my $outputFile = new PDL::IO::HDF5(">".$workDirectory."/modelDiscrepancy/randomSpins/discrepancy".ucfirst($constraintDefinition->{'label'}).".hdf5");
    $outputFile->dataset('multiplicativeCovariance')->set($modelDiscrepancyCovarianceMultiplicative);
    $outputFile->attrSet(
     	description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to use of random halo spins."
     	);
    # Add results of systematics models.
    my $systematicGroup = $outputFile->group("systematicModels");
    foreach my $model ( keys(%systematicResults) ) {
	my $modelGroup = $systematicGroup->group($model);
	my %modelResults = %{$systematicResults{$model}};
	foreach my $parameter ( keys(%modelResults) ) {
	    $modelGroup->attrSet($parameter => $modelResults{$parameter});
	}
    }
    # Create a plot if requested.
    if ( $arguments{'plot'} eq "yes" ) {
	# Determine suitable x and y ranges.
	my $x    = 10.0**$shortX;
	my $xMin = 0.67*min($x     );
	my $xMax = 1.50*max($x     );
	my $yMin = 0.67*min($shortY->append($standardY));
	my $yMax = 1.50*max($shortY->append($standardY));
	# Declare standards for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	(my $label = ucfirst($constraintDefinition->{'label'})) =~ s/\./_/;
	my $plotFile = $workDirectory."/modelDiscrepancy/randomSpins/discrepancy".$label.".pdf";
	($plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
	open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$plotFileEPS."'\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.4,0.2\n";
	print $gnuPlot "set key left\n";
	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [".$xMin.":".$xMax."]\n";
	print $gnuPlot "set yrange [".$yMin.":".$yMax."]\n";
	print $gnuPlot "set xlabel '\$x\$'\n";
	print $gnuPlot "set ylabel '\$y\$'\n";
	&PrettyPlots::Prepare_Dataset(\$plot,
				      $x,$shortY,
				      errorUp   => sqrt($shortCovariance->diagonal(0,1)),
				      errorDown => sqrt($shortCovariance->diagonal(0,1)),
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
				      title     => "Short"
	    );
	&PrettyPlots::Prepare_Dataset(\$plot,
				      $x,$standardY,
				      errorUp   => sqrt($standardCovariance->diagonal(0,1)),
				      errorDown => sqrt($standardCovariance->diagonal(0,1)),
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'mediumSeaGreen'},
				      title     => "Standard"
	    );
	&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
    }
}

exit;
