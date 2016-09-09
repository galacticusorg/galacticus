#!/usr/bin/env perl
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use strict;
use warnings;
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::MatrixOps;
use PDL::LinearAlgebra;
use PDL::IO::HDF5;
use UNIVERSAL;
use Data::Dumper;
use Clone qw(clone);
use Text::Table;
use Math::SigFigs;
use Galacticus::Constraints::Parameters;
use Galacticus::Constraints::Covariances;
use Galacticus::Launch::PBS;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Run calculations to test the convergence of a Galacticus model for the given constraints compilation.
# Andrew Benson (14-November-2012)

# Get arguments.
die("Usage: testConvergence.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make      => "yes"        ,
     directory => "convergence"
    );
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::Parse_Config($configFile);

# Validate the config file.
die("testConvergence.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("testConvergence.pl: compilation must be specified in config file"   ) unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("testConvergence.pl: baseParameters must be specified in config file")
    unless ( exists($config->{'likelihood'}->{'baseParameters'}) || exists($arguments{'baseParameters'}) );

# Determine base parameters.
my $baseParameters;
if ( exists($arguments{'baseParameters'}) ) {
    $baseParameters = $arguments{'baseParameters'};
} else {
    $baseParameters = $config->{'likelihood'}->{'baseParameters'};
}

# Determine the scratch and work directories.
my $workDirectory    = $config->{'likelihood'}->{'workDirectory'};
my $scratchDirectory = $config->{'likelihood'}->{'workDirectory'};
$scratchDirectory    = $config->{'likelihood'}->{'scratchDirectory'}
    if ( exists($config->{'likelihood'}->{'scratchDirectory'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'});

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("testConvergence.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Galacticus::Constraints::Parameters::Compilation($config->{'likelihood'}->{'compilation'},$baseParameters);
my @constraints = @{$constraintsRef};

# Set an initial random number seed.
$parameters->{'randomSeed'}->{'value'} = 824;

# Set a dummy baseline variable.
$parameters->{'baseline'  }->{'value'} = 1.0;

# Ensure that trees will be built, not read.
$parameters->{'mergerTreeConstructMethod'}->{'value'} = "build";

# Set the number of trees per decade of halo mass if specified on the command line.
$parameters->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $arguments{'treesPerDecade'}
   if ( exists($arguments{'treesPerDecade'}) );

# Define parameters to test for convergence.
my @convergences =
 (
 {
     parameter => "baseline",
     factor    => 1.1,
     steps     => 21,
     ideal     => "smallest"
 },
 {
     parameter => "mergerTreeBuildHaloMassMinimum",
     factor    => 1.259,
     steps     => 21,
     ideal     => "smallest"
 },
 {
     parameter => "mergerTreeBuildHaloMassMaximum",
     factor    => 0.794,
     steps     => 21,
     ideal     => "largest"
 },
 {
     parameter => "mergerTreeBuilderMethod:mergeProbability",
     factor    => 1.259,
     steps     => 21,
     ideal     => "smallest"
 },
 {
     parameter => "mergerTreeBuilderMethod:accretionLimit",
     factor    => 1.259,
     steps     => 21,
     ideal     => "smallest",
     maximum   => 1.0
 },
 {
     parameter => "modifiedPressSchechterFirstOrderAccuracy",
     factor    => 1.259,
     steps     => 21,
     ideal     => "smallest"
 },
 {
     parameter => "timestepHostAbsolute",
     factor    => 0.794,
     steps     => 21,
     ideal     => "smallest"
 },
 {
     parameter => "timestepHostRelative",
     factor    => 0.794,
     steps     => 21,
     ideal     => "smallest"
 },
 {
     parameter => "timestepSimpleAbsolute",
     factor    => 0.794,
     steps     => 21,
     ideal     => "smallest"
 },
 {
     parameter => "timestepSimpleRelative",
     factor    => 0.794,
     steps     => 31,
     ideal     => "smallest"
 },
 {
     parameter => "odeToleranceAbsolute",
     factor    => 0.794,
     steps     => 21,
     ideal     => "smallest"
 },
 {
     parameter => "odeToleranceRelative",
     factor    => 0.794,
     steps     => 21,
     ideal     => "smallest"
 }
);

# Add convergence checks for mass resolution.
if ( $parameters->{'mergerTreeMassResolutionMethod'}->{'value'} eq "fixed" ) {
    push(
	@convergences,
	{
	    parameter => "mergerTreeMassResolutionMethod:massResolution",
	    factor    => 0.794,
	    steps     => 11,
	    ideal     => "smallest"
	}
	);
} elsif ( $parameters->{'mergerTreeMassResolutionMethod'}->{'value'} eq "scaled" ) {
    push(
	@convergences,
	{
	    parameter => "mergerTreeMassResolutionMethod:massResolutionMinimum",
	    factor    => 0.794,
	    steps     => 11,
	    ideal     => "smallest"
	},
	{
	    parameter => "mergerTreeMassResolutionMethod:massResolutionFractional",
	    factor    => 0.794,
	    steps     => 11,
	    ideal     => "smallest"
	}
	);
}

# Iterate over parameters to run models that will test for convergence.
my @pbsStack;
foreach my $convergence ( @convergences ) {
    # Report on activity.
    print "Running convergence models for: ".$convergence->{'parameter'}."\n";    
    # Make a copy of the parameters.
    my $currentParameters = clone($parameters);
    # Find the active parameter.
    my $activeParameter = $currentParameters;
    $activeParameter = $activeParameter->{$_}
        foreach ( split(":",$convergence->{'parameter'}) );
    # Override parameters.
    foreach ( keys(%arguments) ) {
	if ( $_ =~ m/^parameterOverride:(.+)/ ) {
	    my $parameterName = $1;
	    if ( $arguments{$_} =~ m/^\%\%(.*)\%\%$/ ) {
		my $copyParameter = $1;
		$currentParameters->{$parameterName}->{'value'} = $currentParameters->{$copyParameter}->{'value'};
	    } else {
		$currentParameters->{$parameterName}->{'value'} = $arguments{$_};
	    }
	}
    }
    # Step through values of this parameter.
    for(my $i=0;$i<$convergence->{'steps'};++$i) {
	# Adjust the parameter.
	$activeParameter->{'value'} *= $convergence->{'factor'}
	   if ( $i > 0 );
	# Check if parameter is in range.
	next
	    if 
	    (
	     exists($convergence->{'maximum'})
	     &&
	     $activeParameter->{'value'} >= $convergence->{'maximum'}
	    );
	# Create a directory for output.
	(my $parameterSafe = $convergence->{'parameter'}) =~ s/://g;
	my $modelDirectory = $workDirectory."/".$arguments{'directory'}."/".$parameterSafe."/".$i;
	system("mkdir -p ".$modelDirectory);
	# Specify the output file name.
	my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
	$currentParameters->{'galacticusOutputFileName'                }->{'value'} = $galacticusFileName;
	# Increment the random number seed.
        $currentParameters->{'randomSeed'                              }->{'value'} += 12;
	# Switch off resource sharing.
        $currentParameters->{'treeEvolveThreadLock'                    }->{'value'} = "false";
	# Switch off fixed random seeds.
        $currentParameters->{'mergerTreeBuildCole2000FixedRandomSeeds' }->{'value'} = "false";
	# Ensure we randomly sample from the halo mass function.
        $currentParameters->{'mergerTreeBuildTreesHaloMassDistribution'}->{'value'} = "random";
	# Check if the model has already been run.
	unless ( -e $galacticusFileName ) {
	    # Generate the parameter file.
	    my $parameterFileName = $modelDirectory."/parameters.xml";
	    &Galacticus::Constraints::Parameters::Output($currentParameters,$parameterFileName);
	    # Create a PBS job.
	    my $command = "mpirun --bynode -np 1 Galacticus.exe ".$parameterFileName."\n";
	    foreach my $constraint ( @constraints ) {
		# Parse the definition file.
		my $xml = new XML::Simple;
		my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
		# Insert code to run the analysis code.
		my $analysisCode = $constraintDefinition->{'analysis'};
		$command .= $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5\n";
	    }
	    my %job =
		(
		 launchFile  => $modelDirectory."/launch.pbs",
		 label       => "convergence".ucfirst($parameterSafe).$i,
		 logFile     => $modelDirectory."/launch.log",
		 ppn         => 16,
		 command     => $command,
		 environment =>
		 [
		  "LD_LIBRARY_PATH=/home/abenson/Galacticus/Tools/lib:/home/abenson/Galacticus/Tools/lib64:\$LD_LIBRARY_PATH",
		  "PATH=/home/abenson/Galacticus/Tools/bin:\$PATH",
		  "GFORTRAN_ERROR_DUMPCORE=NO"
		 ]
		);
	    # Queue the calculation.
	    push(
		@pbsStack,
		\%job
		);
	}
    }
}
# Send jobs to PBS.
&Galacticus::Launch::PBS::SubmitJobs(\%arguments,@pbsStack)
    if ( scalar(@pbsStack) > 0 );
# Iterate over constraints.
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $xml = new XML::Simple;
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Create a table for our output report.
    my $reportTable = Text::Table->new(
	{
	    title  => "Parameter",
	    align  => "left"
	},
	{
	    is_sep => 1,
	    body   => " "
	},
	{
	    title  => "Convergence measure", 
	    align  => "num"
	},
	{
	    is_sep => 1,
	    body   => " +/- "
	},
	{
	    title  => "Error",
	    align  => "num"
	},
	{
	    is_sep => 1,
	    body   => " ("
	},
	{
	    title  => "Normalized measure",
	    align  => "num"
	},
	{
	    is_sep => 1,
	    body   => ")"
	}
	);
    # Iterate over convergence criteria, creating a plot of the convergence and reporting on its success or failure.
    my $baselineTestStatistic;
    my $baselineTestStatisticVariance;
    foreach my $convergence ( @convergences ) {
	# Report on activity.
	print "Analysing convergence in ".$constraintDefinition->{'name'}." for: ".$convergence->{'parameter'}."\n";
	# Initialize array of results.
	my @results;
	# Make a copy of the parameters.
	my $currentParameters = clone($parameters);
	# Find the active parameter.
	my $activeParameter = $currentParameters;
	$activeParameter = $activeParameter->{$_}
           foreach ( split(":",$convergence->{'parameter'}) );
	# Step through values of this parameter.
	for(my $i=0;$i<$convergence->{'steps'};++$i) {
	    # Adjust the parameter.
	    $activeParameter->{'value'} *= $convergence->{'factor'}
	       if ( $i > 0 );
	    # Check if parameter is in range.
	    next
		if 
		(
		 exists($convergence->{'maximum'})
		 &&
		 $activeParameter->{'value'} >= $convergence->{'maximum'}
		);
	    # Locate the model directory.
	    (my $parameterSafe = $convergence->{'parameter'}) =~ s/://g;
	    my $modelDirectory = $workDirectory."/".$arguments{'directory'}."/".$parameterSafe."/".$i;
	    # Read the results.
	    my $result = new PDL::IO::HDF5($modelDirectory."/".$constraintDefinition->{'label'}.".hdf5");
	    # Store the results for later use.
	    my $x          = $result->dataset('x'         )->get();
	    my $y          = $result->dataset('y'         )->get();
	    my $covariance = $result->dataset('covariance')->get();
	    my $activeParameterCurrent = clone($activeParameter);
	    push
		(
		 @results,
		 {
		     parameter  => $activeParameterCurrent,
		     x          => $x,
		     y          => $y,
		     covariance => $covariance
		 }
		);
	}
	# Construct the convergence measure.
	my $parameter          = pdl [];
	my $convergenceMeasure = pdl [];
	my $convergenceError   = pdl [];
	my $optimalEntry;
	if      ( $convergence->{'ideal'} eq "largest"  ) {
	    $optimalEntry = scalar(@results)-1;
	    $optimalEntry = 0
		if ( $convergence->{'factor'} < 1.0 );
	} elsif ( $convergence->{'ideal'} eq "smallest" ) {
	    $optimalEntry = 0;
	    $optimalEntry = scalar(@results)-1
		if ( $convergence->{'factor'} < 1.0 );
	} else {
	    die("testConvergence.pl: unrecognized ideal for convergence");
	}
	my $optimal            = $results[$optimalEntry];
	foreach ( @results ) {
	    # Construct measure.
	    my $measure;
	    my $measureError;
	    # Find non-empty entries.
	    my $nonEmpty = 
		which
		(
		 ($_      ->{'y'} > 0.0)
		 |
		 ($optimal->{'y'} > 0.0)
		);
	    # Find difference between models and its covariance.
	    my $covariance = $_->{'covariance'}->($nonEmpty,$nonEmpty)+$optimal->{'covariance'}->($nonEmpty,$nonEmpty);
	    # Catch empty covariance matrices.
	    if ( all($covariance == 0.0) ) {
		$measure      = pdl 0.0;
		$measureError = pdl 0.0;
	    } else {
		# Evaluate the convergence measure.
		$measure =
		    &Galacticus::Constraints::Covariances::ComputeLikelihood
		    (
		     $_         ->{'y'}->($nonEmpty),
		     $optimal   ->{'y'}->($nonEmpty),
		     $covariance                    ,
		     quiet                  =>  $arguments{'quiet'} ,
		     inversionMethod        => "eigendecomposition" ,
		     productMethod          => "linearSolver"       ,
		     normalized             => 0                    ,
		     assumePositiveDefinite => 0
		    );
		if ( $measure > 0.0 ) {
		    $measureError = sqrt(2.0*$measure);
		} else {
		    $measure      .=     0.0;
		    $measureError  = pdl 0.0;
		}
	    }
	    $parameter            = $parameter         ->append($_->{'parameter'}->{'value'});
	    $convergenceMeasure   = $convergenceMeasure->append($measure                    );
	    $convergenceError     = $convergenceError  ->append($measureError               );
	}
    	if ( $convergence->{'parameter'} eq "baseline" ) {
	    $baselineTestStatistic         = average( $convergenceMeasure                                                          );
	    $baselineTestStatisticVariance = average(($convergenceMeasure-$baselineTestStatistic)**2)/(nelem($convergenceMeasure)-1);
	} else {
	    $convergenceError =
		($convergenceMeasure/$baselineTestStatistic)
		*sqrt(
		     ($convergenceError/$convergenceMeasure)**2
		    +$baselineTestStatisticVariance/$baselineTestStatistic**2
		);
	    $convergenceMeasure /= $baselineTestStatistic;
	    die('testConvergence.pl: the baseline test statistic is not defined')
		if ( ! defined($baselineTestStatistic) );
	}
	# Report on convergence in this parameter.
	if ( nelem($convergenceMeasure) > 0 ) {
	    my $current       = 0;
	    $current = 1
		if ( $optimal == 0 );
	    my $currentStatus = ($convergenceMeasure->($current)-1.0)/$convergenceError->($current);
	    print "  --> Normalized convergence measure: ".$currentStatus."\n";
	    $reportTable->add(
		$convergence->{'parameter'},
		FormatSigFigs($convergenceMeasure->($current)->sclr(),4),
		FormatSigFigs($convergenceError  ->($current)->sclr(),4),
		FormatSigFigs($currentStatus                 ->sclr(),4)
	    );
	    # Construct a plot of the convergence.
	    my $plotFileName = $constraintDefinition->{'label'}."_".$convergence->{'parameter'};
	    $plotFileName =~ s/\./_/g;
	    $plotFileName = $workDirectory."/".$arguments{'directory'}."/".$plotFileName.".pdf";
	    my $activeUsedParameter = $parameters;
	    $activeUsedParameter = $activeUsedParameter->{$_}
	       foreach ( split(":",$convergence->{'parameter'}) );
	    my $usedValueX = pdl ( $activeUsedParameter->{'value'} );
	    my $usedValueY = pdl ( 1.0 );
	    my $parameterFull = $parameter->append($usedValueX); 
	    my $xMinimum = minimum($parameterFull     )/1.05;
	    my $xMaximum = maximum($parameterFull     )*1.05;
	    my $yMinimum = 0.3;
	    my $yMaximum = maximum($convergenceMeasure)*1.05;
	    $yMaximum = 1.5
		if ( $yMaximum < 1.5 );
	    my ($gnuPlot, $outputFile, $outputFileEPS, $plot);
	    ($outputFileEPS = $plotFileName) =~ s/\.pdf$/.eps/;
	    open($gnuPlot,"|gnuplot");
	    print $gnuPlot "set terminal epslatex color colortext lw 2 7\n";
	    print $gnuPlot "set output '".$outputFileEPS."'\n";
	    print $gnuPlot "set lmargin screen 0.15\n";
	    print $gnuPlot "set rmargin screen 0.95\n";
	    print $gnuPlot "set bmargin screen 0.15\n";
	    print $gnuPlot "set tmargin screen 0.95\n";
	    print $gnuPlot "set key spacing 1.2\n";
	    print $gnuPlot "set key at screen 0.2,0.2\n";
	    print $gnuPlot "set key left\n";
	    print $gnuPlot "set key bottom\n";
	    print $gnuPlot "set logscale xy\n";
	    print $gnuPlot "set mxtics 10\n";
	    print $gnuPlot "set mytics 10\n";
	    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	    print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
	    print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
	    print $gnuPlot "set title 'Convergence of ".$constraintDefinition->{'name'}." with ".$convergence->{'parameter'}."'\n";
	    print $gnuPlot "set xlabel '{\\tt ".$convergence->{'parameter'}."}'\n";
	    print $gnuPlot "set ylabel 'Convergence measure; \$\\chi^2\$'\n";
	    my $convergedX = pdl ( $xMinimum, $xMaximum );
	    my $convergedY = pdl (       1.00,     1.00 );
	    # Add an arrow to show optimal direction.
	    my $xTip = 0.35;
	    $xTip = 0.65
		if ( $convergence->{'ideal'} eq "largest" );
	    print $gnuPlot "set arrow from graph 0.5, first 0.5477 to graph ".$xTip.", first 0.5477 ls 1 lw 5 filled lc rgbcolor \"#3CB371\"\n";
	    &GnuPlot::PrettyPlots::Prepare_Dataset
		(
		 \$plot,
		 $convergedX,
		 $convergedY,
		 style       => "line",
		 weight      => [5,3],
		 color       => $GnuPlot::PrettyPlots::colorPairs{'redYellow'}
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset
		(
		 \$plot,
		 $parameter,
		 $convergenceMeasure,
		 errorUp     => $convergenceError,
		 errorDown   => $convergenceError,
		 style       => "point",
		 symbol      => [6,7],
		 weight      => [5,3],
		 color       => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'}
		);
	    &GnuPlot::PrettyPlots::Prepare_Dataset
		(
		 \$plot,
		 $usedValueX,
		 $usedValueY,
		 style       => "point",
		 symbol      => [6,7],
		 weight      => [5,3],
		 color       => $GnuPlot::PrettyPlots::colorPairs{'indianRed'},
		);
	    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	    close($gnuPlot);
	    &GnuPlot::LaTeX::GnuPlot2PDF($outputFileEPS);
	    # Make a plot of the results.
	    $plotFileName = $constraintDefinition->{'label'}."_".$convergence->{'parameter'}."_results";
	    $plotFileName =~ s/\./_/g;
	    $plotFileName = $workDirectory."/".$arguments{'directory'}."/".$plotFileName.".pdf";
	    ($outputFileEPS = $plotFileName) =~ s/\.pdf$/.eps/;
	    undef($plot);
	    open($gnuPlot,"|gnuplot");
	    print $gnuPlot "set terminal epslatex color colortext lw 2 7\n";
	    print $gnuPlot "set output '".$outputFileEPS."'\n";
	    print $gnuPlot "set lmargin screen 0.15\n";
	    print $gnuPlot "set rmargin screen 0.95\n";
	    print $gnuPlot "set bmargin screen 0.15\n";
	    print $gnuPlot "set tmargin screen 0.95\n";
	    print $gnuPlot "set key spacing 1.2\n";
	    print $gnuPlot "set key at screen 0.2,0.2\n";
	    print $gnuPlot "set key left\n";
	    print $gnuPlot "set key bottom\n";
	    print $gnuPlot "set logscale xy\n";
	    print $gnuPlot "set mxtics 10\n";
	    print $gnuPlot "set mytics 10\n";
	    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	    print $gnuPlot "set title 'Convergence of ".$constraintDefinition->{'name'}." with ".$convergence->{'parameter'}."'\n";
	    print $gnuPlot "set xlabel '\$x\$'\n";
	    print $gnuPlot "set ylabel '\$y\$'\n";
	    $xMinimum = pdl +1.0e30;
	    $xMaximum = pdl -1.0e30;
	    $yMinimum = pdl +1.0e30;
	    $yMaximum = pdl -1.0e30;
	    foreach ( @results ) {
		my $nonZeroY = which($_->{'y'} > 0.0);
		$xMinimum = minimum($_->{'x'})
		    if ( minimum($_->{'x'}) < $xMinimum );
		$xMaximum = maximum($_->{'x'})
		    if ( maximum($_->{'x'}) > $xMaximum );
		$yMinimum = minimum($_->{'y'}->($nonZeroY))
		    if ( minimum($_->{'y'}->($nonZeroY)) < $yMinimum );
		$yMaximum = maximum($_->{'y'}->($nonZeroY))
		    if ( maximum($_->{'y'}->($nonZeroY)) > $yMaximum );
	    }
	    $xMinimum /= 2.0;
	    $xMaximum *= 2.0;
	    $yMinimum /= 2.0;
	    $yMaximum *= 2.0;
	    print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
	    print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
	    my $i = -1;
	    foreach ( @results ) {
		++$i;
		my $fraction = $i/(scalar(@results)-1);
		&GnuPlot::PrettyPlots::Prepare_Dataset
		    (
		     \$plot,
		     $_->{'x'},
		     $_->{'y'},
		     style  => "point",
		     symbol => [6,7],
		     weight => [5,3],
		     color  => [
			 &GnuPlot::PrettyPlots::Color_Gradient($fraction,[0.0,1.0,0.5],[240.0,1.0,0.5]),
			 &GnuPlot::PrettyPlots::Color_Gradient($fraction,[0.0,1.0,0.5],[240.0,1.0,0.5])
		     ]
		    );
	    }
	    print $gnuPlot "set xrange [".$xMinimum.":".$xMaximum."]\n";
	    print $gnuPlot "set yrange [".$yMinimum.":".$yMaximum."]\n";
	    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	    close($gnuPlot);
	    &GnuPlot::LaTeX::GnuPlot2PDF($outputFileEPS);
	}
    }
    # Write the convergence report.
    open(my $report,">".$workDirectory."/".$arguments{'directory'}."/".$constraintDefinition->{'label'}."Report.txt");
    print $report "Report on convergence in Galacticus numerical parameters\n\n";
    print $report "For constraint: ".$constraintDefinition->{'name'}."\n\n";
    print $report "The convergence measure, computed for the current default value of the\n";
    print $report "parameter in question, is the root mean square difference between the\n";
    print $report "relevant statistic (as computed for this constraint) and that for the\n";
    print $report "same statistic from the most \"ideal\" model (i.e. the model in which\n";
    print $report "the parameter takes the value which should minimize numerical\n";
    print $report "artefacts), divided by the error in that statistic. number should be\n";
    print $report "consistent with unity if convergence is achieved. The normalized\n";
    print $report "convergence measure is (C-1)/E where C is the convergence measure and E\n";
    print $report "is its error. It approximately represents the number of sigma deviation\n";
    print $report "from convergence.\n\n";
    print $report $reportTable->table()."\n";
    close($report);
}

exit;
