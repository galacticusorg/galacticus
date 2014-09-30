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
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require Galacticus::Launch::PBS;

# Run calculations to determine the model discrepancy arising from the use of variable mass resolution when building trees.
# Andrew Benson (28-January-2013)

# Get arguments.
die("Usage: variableTreeMassResolution.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make         => "yes" ,
     plot         => "no"  ,
     massFraction => 1.0e-3
    );
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Parse the constraint config file.
my $xml    = new XML::Simple;
my $config = $xml->XMLin($configFile, KeyAttr => 0);

# Validate the config file.
die("variableTreeMassResolution.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("variableTreeMassResolution.pl: compilation must be specified in config file"   ) unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("variableTreeMassResolution.pl: baseParameters must be specified in config file") unless ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Determine the scratch and work directories.
my $workDirectory    = $config->{'likelihood'}->{'workDirectory'};
my $scratchDirectory = $config->{'likelihood'}->{'workDirectory'};
$scratchDirectory    = $config->{'likelihood'}->{'scratchDirectory'} if ( exists($config->{'likelihood'}->{'scratchDirectory'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'});

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("variableTreeMassResolution.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'likelihood'}->{'compilation'},$config->{'likelihood'}->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Switch off thread locking.
$parameters->{'parameter'}->{'treeEvolveThreadLock'}->{'value'} = "false";

# Set the mass resolution.
if ( exists($arguments{'massResolution'}) ) {
    $parameters->{'parameter'}->{'mergerTreeBuildMassResolutionFixed'        }->{'value'} = $arguments{'massResolution'};
    $parameters->{'parameter'}->{'mergerTreeBuildMassResolutionScaledMinimum'}->{'value'} = $arguments{'massResolution'};
}

# Initialize a stack for PBS models.
my @pbsStack;

# Specify models to run.
my @models = 
    (
     {
	 label      => "fixedResolution",
	 parameters =>
	     [
	      # Switch to using fixed mass resolution.
	      {
		  name  => "mergerTreesBuildMassResolutionMethod",
		  value => "fixed"
	      }
	     ]
     },
     {
	 label      => "variableResolution",
	 parameters => 
	     [
	      # Switch to using variable mass resolution.
	      {
		  name  => "mergerTreesBuildMassResolutionMethod",
		  value => "scaled"
	      },
	     ]
     }
    );

# Iterate over models.
foreach my $model ( @models ) {
# Specify the output name.
    my $modelDirectory = $workDirectory."/modelDiscrepancy/variableTreeMassResolution/".$model->{'label'};
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
	    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
	    # Insert code to run the analysis code.
	    my $analysisCode = $constraintDefinition->{'analysis'};
	    $command .= $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5\n";
	}
	# Queue the calculation.
	my %job =
	    (
	     launchFile => $modelDirectory."/launch.pbs",
	     label      => "variableTreeMassResolution".$model->{'label'},
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
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Locate the model results.
    my $variableResolutionResultFileName = $workDirectory."/modelDiscrepancy/variableTreeMassResolution/variableResolution/".$constraintDefinition->{'label'}.".hdf5";
    my $fixedResolutionResultFileName    = $workDirectory."/modelDiscrepancy/variableTreeMassResolution/fixedResolution/"   .$constraintDefinition->{'label'}.".hdf5";
    # Read the results.
    my $variableResolutionResult = new PDL::IO::HDF5($variableResolutionResultFileName);
    my $fixedResolutionResult    = new PDL::IO::HDF5($fixedResolutionResultFileName   );
    # Extract the results.
    my $fixedX             = $fixedResolutionResult   ->dataset('x'         )->get();
    my $fixedY             = $fixedResolutionResult   ->dataset('y'         )->get();
    my $fixedCovariance    = $variableResolutionResult->dataset('covariance')->get();
    my $variableY          = $variableResolutionResult->dataset('y'         )->get();
    my $variableCovariance = $variableResolutionResult->dataset('covariance')->get();
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
		    $fixedX              ,
		    $fixedY              ,
		    $fixedCovariance     ,
		    $variableY           ,
		    $variableCovariance
		);
	    }
	}
    }
    # Find the multiplicative discrepancy between these two models.
    (my $nonZero, my $zero)            = which_both($fixedY > 0.0);
    my $modelDiscrepancyMultiplicative = $variableY->copy();
    $modelDiscrepancyMultiplicative->($nonZero) /= $fixedY->($nonZero);
    $modelDiscrepancyMultiplicative->($zero   ) .= 1.0;
    # Compute the covariance.
    my $modelDiscrepancyCovarianceMultiplicative = 
	 $variableCovariance*outer(       1.0/$fixedY   ,       1.0/$fixedY   )
	+$fixedCovariance   *outer($variableY/$fixedY**2,$variableY/$fixedY**2);
    # Compute additional covariance arising from the fact that the mass
    # systematic model is imperfect. We assume here that the shift in
    # abundance should be by a constant factor. Any deviation from this
    # is therefor taken to be a failing of the mass systematic model and
    # is added to the covariance.
    my $multiplierNormalized                   = $modelDiscrepancyMultiplicative->copy();
    $multiplierNormalized                     /= $modelDiscrepancyMultiplicative->((0));
    my $covarianceExtra                        = outer(log($multiplierNormalized),log($multiplierNormalized));
    $modelDiscrepancyMultiplicative           .= $modelDiscrepancyMultiplicative->((0));
    $modelDiscrepancyCovarianceMultiplicative += $covarianceExtra;
    # Output the model discrepancy to file.
    my $outputFile = new PDL::IO::HDF5(">".$workDirectory."/modelDiscrepancy/variableTreeMassResolution/discrepancy".ucfirst($constraintDefinition->{'label'}).".hdf5");
    $outputFile->dataset('multiplicative'          )->set($modelDiscrepancyMultiplicative          );
    $outputFile->dataset('multiplicativeCovariance')->set($modelDiscrepancyCovarianceMultiplicative);
    $outputFile->attrSet(
	description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to use of variable mass resolution in merger trees."
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
	my $x    = 10.0**$fixedX;
	my $xMin = 0.67*min($x     );
	my $xMax = 1.50*max($x     );
	my $yMin = 0.67*min($fixedY->append($variableY));
	my $yMax = 1.50*max($fixedY->append($variableY));
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	(my $label = ucfirst($constraintDefinition->{'label'})) =~ s/\./_/;
	my $plotFile = $workDirectory."/modelDiscrepancy/variableTreeMassResolution/discrepancy".$label.".pdf";
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
				      $x,$fixedY,
				      errorUp   => sqrt($fixedCovariance->diagonal(0,1)),
				      errorDown => sqrt($fixedCovariance->diagonal(0,1)),
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
				      title     => "Fixed"
	    );
	&PrettyPlots::Prepare_Dataset(\$plot,
				      $x,$variableY,
				      errorUp   => sqrt($variableCovariance->diagonal(0,1)),
				      errorDown => sqrt($variableCovariance->diagonal(0,1)),
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'mediumSeaGreen'},
				      title     => "Variable"
	    );
	&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
    }
}

exit;
