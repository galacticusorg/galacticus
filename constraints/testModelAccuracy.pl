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
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use UNIVERSAL;
use Data::Dumper;
use Clone qw(clone);
use Text::Table;
use Math::SigFigs;
use PDL::IO::HDF5;
use PDL::LinearAlgebra;
use PDL::Matrix;
use PDL::MatrixOps;
require Galacticus::Constraints::Parameters;
require Galacticus::Constraints::Covariances;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Run calculations to test the accuracy a Galacticus model for the given constraints compilation.
# Andrew Benson (15-November-2012)

# Get arguments.
die("Usage: testModelAccuracy.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make => "yes"
    );
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Parse the constraint config file.
my $config = &Parameters::Parse_Config($configFile);

# Validate the config file.
die("testModelAccuracy.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("testModelAccuracy.pl: compilation must be specified in config file"   ) unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("testModelAccuracy.pl: baseParameters must be specified in config file") unless ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Determine the scratch and work directories.
my $workDirectory    = $config->{'likelihood'}->{'workDirectory'   };
my $scratchDirectory = $config->{'likelihood'}->{'workDirectory'   };
$scratchDirectory    = $config->{'likelihood'}->{'scratchDirectory'}
    if ( exists($config->{'likelihood'}->{'scratchDirectory'}) );

# Create the work and scratch directories.
system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'});

# Ensure that Galacticus is built.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("testModelAccuracy.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'likelihood'}->{'compilation'},$config->{'likelihood'}->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Set an initial random number seed.
$parameters->{'parameter'}->{'randomSeed'}->{'value'} = 824;

# Switch off thread locking.
$parameters->{'parameter'}->{'treeEvolveThreadLock'}->{'value'} = "false";

# Set the default number of trees per decade.
$parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $arguments{'treesPerDecade'}
	if ( exists($arguments{'treesPerDecade'}) );

# Ensure no abundance limits are applied for halo mass function sampling.
$parameters->{'parameter'}->{'haloMassFunctionSamplingAbundanceMinimum'}->{'value'} = -1.0;
$parameters->{'parameter'}->{'haloMassFunctionSamplingAbundanceMaximum'}->{'value'} = -1.0;

# Define parameters to test for accuracy.
my @accuracies =
(
 {
     parameter => "mergerTreeBuildTreesPerDecade",
     factor    => 0.5,
     steps     => 8
 },
);

# Define sampling methods.
my @samplingMethods = 
    ( 
      {
	  name  => "p1_0.0_p2_0.0"            ,
	  label => "\$(p_1,p_2)=( 0.0, 0.0)\$",
	  p1    =>  0.0                       ,
	  p2    =>  0.0
      },
      {
	  name  => "p1_m1.0_p2_0.0"           ,
	  label => "\$(p_1,p_2)=(-1.0, 0.0)\$",
	  p1    => -1.0                       ,
	  p2    =>  0.0
      },
      {
	  name  => "p1_p1.0_p2_0.0"           ,
	  label => "\$(p_1,p_2)=(+1.0, 0.0)\$",
	  p1    => +1.0                       ,
	  p2    =>  0.0
      }
    );

# Iterate over parameters to run models that will test for accuracy.
my @pbsStack;
foreach my $accuracy ( @accuracies ) {
    # Report on activity.
    print "Running accuracy models for: ".$accuracy->{'parameter'}."\n";
    # Make a copy of the parameters.
    my $currentParameters = clone($parameters);
    # Step through values of this parameter.
    for(my $i=0;$i<$accuracy->{'steps'};++$i) {
	# Adjust the parameter.
	$currentParameters->{'parameter'}->{$accuracy->{'parameter'}}->{'value'} *= $accuracy->{'factor'}
	    if ( $i > 0 );
	# Iterate over sampling methods.
	foreach my $samplingMethod ( @samplingMethods ) {
	    # Create a directory for output.
	    my $modelDirectory = $workDirectory."/accuracy/".$accuracy->{'parameter'}."/".$samplingMethod->{'name'}.$i;
	    system("mkdir -p ".$modelDirectory);
	    # Specify the sampling method.
	    $currentParameters->{'parameter'}->{'haloMassFunctionSamplingMethod'          }->{'value'} = "haloMassFunction";
	    $currentParameters->{'parameter'}->{'haloMassFunctionSamplingAbundanceMinimum'}->{'value'} = -1.0;
	    $currentParameters->{'parameter'}->{'haloMassFunctionSamplingAbundanceMaximum'}->{'value'} = -1.0;
	    $currentParameters->{'parameter'}->{'haloMassFunctionSamplingModifier1'       }->{'value'} = $samplingMethod->{'p1'};
	    $currentParameters->{'parameter'}->{'haloMassFunctionSamplingModifier2'       }->{'value'} = $samplingMethod->{'p2'};
	    # Specify the output file name.
	    my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
	    $currentParameters->{'parameter'}->{'galacticusOutputFileName'}->{'value'} = $galacticusFileName;
	    # Increment the random number seed.
	    $currentParameters->{'parameter'}->{'randomSeed'}->{'value'} += 1;
	    # Check if the model has already been run.
	    unless ( -e $galacticusFileName ) {
		# Generate the parameter file.
		my $parameterFileName = $modelDirectory."/parameters.xml";
		&Parameters::Output($currentParameters,$parameterFileName);
		# Create a batch script for PBS.
		my $batchScriptFileName = $modelDirectory."/launch.pbs";
		open(oHndl,">".$batchScriptFileName);
		print oHndl "#!/bin/bash\n";
		print oHndl "#PBS -N accuracy".ucfirst($accuracy->{'parameter'}).$samplingMethod->{'name'}.$i."\n";
		print oHndl "#PBS -l nodes=1:ppn=12\n";
		print oHndl "#PBS -j oe\n";
		print oHndl "#PBS -o ".$modelDirectory."/launch.log\n";
		print oHndl "#PBS -V\n";
		print oHndl "cd \$PBS_O_WORKDIR\n";
		print oHndl "export LD_LIBRARY_PATH=/home/abenson/Galacticus/Tools/lib:/home/abenson/Galacticus/Tools/lib64:\$LD_LIBRARY_PATH\n";
		print oHndl "export PATH=/home/abenson/Galacticus/Tools/bin:\$PATH\n";
		print oHndl "export GFORTRAN_ERROR_DUMPCORE=YES\n";
		print oHndl "ulimit -t unlimited\n";
		print oHndl "ulimit -c unlimited\n";
		print oHndl "export OMP_NUM_THREADS=12\n";
		print oHndl "mpirun --bynode -np 1 /usr/bin/time --format='Time: %S %U' Galacticus.exe ".$parameterFileName."\n";
		foreach my $constraint ( @constraints ) {
		    # Parse the definition file.
		    my $xml = new XML::Simple;
		    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
		    # Insert code to run the analysis code.
		    my $analysisCode = $constraintDefinition->{'analysis'};
		    print oHndl $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5\n";
		}
		close(oHndl);
		# Queue the calculation.
		push(
		    @pbsStack,
		    $batchScriptFileName
		    );
	    }
	}
    }
    # Send jobs to PBS.
    &PBS_Submit(@pbsStack)
	if ( scalar(@pbsStack) > 0 );
}

# Iterative over accuracy criteria, creating a plot of the accuracy.
foreach my $accuracy ( @accuracies ) {
    # Report on activity.
    print "Analysing accuracy in ".$accuracy->{'parameter'}."\n";
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
	    title  => "Sampling p1", 
	    align  => "num"
	},
	{
	    is_sep => 1,
	    body   => " "
	},
	{
	    title  => "Sampling p2", 
	    align  => "num"
	},
	{
	    is_sep => 1,
	    body   => " "
	},
	{
	    title  => "Accuracy measure", 
	    align  => "num"
	},
	{
	    is_sep => 1,
	    body   => " "
	},
	{
	    title  => "CPU time [s]", 
	    align  => "num"
	}
	);
    # Begin constructing a plot of the accuracy.
    my $plotFileName = $workDirectory."/accuracy/".$accuracy->{'parameter'}.".pdf";
    my ($gnuPlot, $outputFile, $outputFileEPS, $plot);
    ($outputFileEPS = $plotFileName) =~ s/\.pdf$/.eps/;
    open($gnuPlot,"|gnuplot");
    print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
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
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [1.0e2:1.0e8]\n";
    print $gnuPlot "set yrange [1.0e-1:1.0e4]\n";
    print $gnuPlot "set title 'Convergence with tree processing time'\n";
    print $gnuPlot "set xlabel 'Tree processing time [s]'\n";
    print $gnuPlot "set ylabel 'Convergence measured []'\n";
    # Iterate over sampling methods.
    my $i = -1;
    foreach my $samplingMethod ( @samplingMethods ) {
	++$i;
	# Initialize array of results.
	my $parameter       = pdl [];
	my $timing          = pdl [];
	my $accuracyMeasure = pdl [];
	# Make a copy of the parameters.
	my $currentParameters = clone($parameters);
	# Step through values of this parameter.
	for(my $i=0;$i<$accuracy->{'steps'};++$i) {
	    # Adjust the parameter.
	    $currentParameters->{'parameter'}->{$accuracy->{'parameter'}}->{'value'} *= $accuracy->{'factor'}
	        if ( $i > 0 );
	    # Locate the model directory.
	    my $modelDirectory = $workDirectory."/accuracy/".$accuracy->{'parameter'}."/".$samplingMethod->{'name'}.$i."/";
	    # Get the model timing.
	    my $modelTiming;
	    open(my $logFile,$modelDirectory."launch.log");	    
	    while ( my $line = <$logFile> ) {
		$modelTiming = $1+$2
		    if ( $line =~ m/^Time:\s+([0-9\.]+)\s+([0-9\.]+)/ );
	    }
	    close($logFile);
	    # Iterate over constraints.
	    my $thisAccuracyMeasure = pdl 0.0;
	    foreach my $constraint ( @constraints ) {
		# Parse the definition file.
		my $xml = new XML::Simple;
		my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
		# Check if the accuracy data was reported. If not, assume this model just didn't complete
		# in the available time.
		my $resultFileName = $modelDirectory."/".$constraintDefinition->{'label'}.".hdf5";
		if ( -e $resultFileName ) {
		    # Open the file.
		    my $resultFile = new PDL::IO::HDF5($resultFileName);
		    # Read covariance matrices.
		    my $covarianceModel = $resultFile->dataset('covariance'    )->get();
		    my $covarianceData  = $resultFile->dataset('covarianceData')->get();
		    # Normalize the matrices.
		    my $normalization = $covarianceData->((0),(0))->copy();
		    $covarianceModel /= $normalization;
		    $covarianceData  /= $normalization;
		    # For any empty rows in the model covariance, set the diagonal to a small value.
		    my $epsilon                                = pdl 1.0e-3;		    
		    my $empty                                  = which($covarianceModel->diagonal(0,1) < $epsilon*$covarianceData->diagonal(0,1));
		    $covarianceModel->diagonal(0,1)->($empty) .= $epsilon*$covarianceData->diagonal(0,1)->($empty);
		    # Do Cholesky decompositions.
		    my $choleskyData  = mchol($covarianceData );
		    my $choleskyModel = eval{mchol($covarianceModel)};
		    if ( defined($choleskyModel) ) {
			# Construct P and Q matrices.
			my $P = minv($choleskyData) x transpose($choleskyModel);
			my $Q = minv($choleskyData) x           $choleskyModel;
			# Construct accuracy measure matrix.
			my $PTP = transpose($P) x $P;
			$PTP->diagonal(0,1) += 1.0;
			my $PTPinv = minv($PTP);
			my $accuracyMatrix = $Q x $PTPinv x transpose($Q);
			# Construct accuracy measure.
			$thisAccuracyMeasure += sum($accuracyMatrix);
		    }
		}
	    }
	    # Append this to the results arrays.
	    $parameter       = $parameter      ->append($currentParameters->{'parameter'}->{$accuracy->{'parameter'}}->{'value'});
	    $accuracyMeasure = $accuracyMeasure->append($thisAccuracyMeasure                                                    );
	    $timing          = $timing         ->append($modelTiming                                                            );
	    # Append to the table of results.
	    $reportTable->add(
	    	$accuracy->{'parameter'}             ,
	    	$samplingMethod->{'p1'}              ,
	    	$samplingMethod->{'p2'}              ,
	    	FormatSigFigs($thisAccuracyMeasure,5),
	    	FormatSigFigs($modelTiming        ,5)
	    	);
	}
	# Add to the plot.
	&PrettyPlots::Prepare_Dataset
	    (
	     \$plot,
	     $timing,
	     $accuracyMeasure,
	     style  => "point",
	     symbol => [6,7],
	     weight => [5,3],
	     color  => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[$i]},
	     title  => $samplingMethod->{'name'}
	    );
    }
    # Finalize the plot.
    &PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &LaTeX::GnuPlot2PDF($outputFileEPS);    
    # Write the accuracy report.
    open(my $report,">".$workDirectory."/accuracy/".$accuracy->{'parameter'}."Report.txt");
    print $report "Report on accuracy in Galacticus model\n\n";
    print $report "For parameter: ".$accuracy->{'parameter'}."\n\n";
    print $report $reportTable->table()."\n";
    close($report);
}

exit;

sub PBS_Submit {
    # Submit jobs to PBS and wait for them to finish.
    my @pbsStack = @_;
    my %pbsJobs;
    # Determine maximum number allowed in queue at once.
    my $jobMaximum = 10;
    $jobMaximum = $arguments{'pbsJobMaximum'}
        if ( exists($arguments{'pbsJobMaximum'}) );
    # Submit jobs and wait.
    print "Waiting for PBS jobs to finish...\n";
    while ( scalar(keys %pbsJobs) > 0 || scalar(@pbsStack) > 0 ) {
	# Find all PBS jobs that are running.
	my %runningPBSJobs;
	undef(%runningPBSJobs);
	open(pHndl,"qstat -f|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^Job\sId:\s+(\S+)/ ) {$runningPBSJobs{$1} = 1};
	}
	close(pHndl);
	foreach my $jobID ( keys(%pbsJobs) ) {
	    unless ( exists($runningPBSJobs{$jobID}) ) {
		print "PBS job ".$jobID." has finished.\n";
		# Remove the job ID from the list of active PBS jobs.
		delete($pbsJobs{$jobID});
	    }
	}
	# If fewer than ten jobs are in the queue, pop one off the stack.
	if ( scalar(@pbsStack) > 0 && scalar(keys %pbsJobs) < 20 ) {
	    my $batchScript = pop(@pbsStack);
	    # Submit the PBS job.
	    open(pHndl,"qsub ".$batchScript."|");
	    my $jobID = "";
	    while ( my $line = <pHndl> ) {
	    	if ( $line =~ m/^(\d+\S+)/ ) {$jobID = $1};
	    }
	    close(pHndl);	    
	    # Add the job number to the active job hash.
	    unless ( $jobID eq "" ) {
	    	$pbsJobs{$jobID} = 1;
	    }
	    sleep 5;
	} else {
	    # Wait.
	    sleep 60;
	}
    }
}
