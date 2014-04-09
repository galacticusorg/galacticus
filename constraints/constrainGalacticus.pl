#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use Data::Dumper;
use Fcntl qw(:DEFAULT :flock);
use MIME::Lite;
use List::Util qw(min max);
use Storable;
require File::Which;
require File::NFSLock;
require System::Redirect;
require Galacticus::Constraints::Parameters;
require List::ExtraUtils;
# Script timing
# use Time::HiRes qw(gettimeofday tv_interval);
# my $self      = $0;
# my $timeStart = [gettimeofday];

# Finds constraints on Galacticus parameters.
# Andrew Benson (02-September-2011)

# Get command-line parameters.
die("Usage: constrainGalacticus.pl <configFile> <mpiRank> <likelihoodFile> <temperature> <param1> [<param2>......]") 
    unless ( scalar(@ARGV) > 4 );
my $configFile = $ARGV[0];
my $config     = &Parameters::Parse_Config($configFile,useStored => 1);
my @parameters;
if ( UNIVERSAL::isa($config->{'parameters'}->{'parameter'},"ARRAY") ) {
    @parameters = @{$config->{'parameters'}->{'parameter'}};
} else {
    push(@parameters,$config->{'parameters'}->{'parameter'});
}

# Count active parameters.
my $parameterCount = 0;
for(my $i=0;$i<scalar(@parameters);++$i) {
    ++$parameterCount 
	if ( exists($parameters[$i]->{'prior'}) );
}

# Get the MPI rank.
my $mpiRank        = $ARGV[1];

# Get the name for the likelihood file.
my $likelihoodFile = $ARGV[2];

# Get the temperature.
my $temperature    = $ARGV[3];

# Convert command line arguments to a parameter structure.
die("constrainGalacticus.pl: number of supplied arguments does not match number of parameters") 
    unless ( scalar(@ARGV) == $parameterCount+4 );
my $j = -1;
my %parameterValues;
for(my $i=0;$i<scalar(@parameters);++$i) {
    if ( exists($parameters[$i]->{'prior'}) ) {
	++$j;
	$parameterValues{$parameters[$i]->{'name'}} = $ARGV[$j+4];
    }
}

# Set the values of any parameters that are defined in terms of other parameters.
my $failCount = 1;
while ( $failCount > 0 ) {
    $failCount = 0;
    for(my $i=0;$i<scalar(@parameters);++$i) {
	if ( exists($parameters[$i]->{'define'}) ) {
	    die ("constrainGalacticusWrapper.pl: cannot specify a prior for a defined parameter")
		if ( exists($parameters[$i]->{'prior'}) );
	    # Attempt to replace named parameters in the definition with their values.
	    while ( $parameters[$i]->{'define'} =~ m/\%([a-zA-Z0-9_]+)/ ) {
		my $parameterName = $1;
		if ( exists($parameterValues{$parameterName}) ) {
		    $parameters[$i]->{'define'} =~ s/\%$parameterName/$parameterValues{$parameterName}/g;
		} else {
		    ++$failCount;
		    last;
		}
		$parameterValues{$parameters[$i]->{'name'}} = eval($parameters[$i]->{'define'})
		    unless ( $parameters[$i]->{'define'} =~ m/\%([a-zA-Z0-9_]+)/ );
	    }
	}
    }
}

# Create an array of new parameters.
my $newParameters;
for(my $i=0;$i<scalar(@parameters);++$i) {
    push(
	 @{$newParameters->{'parameter'}},
	 {
	     name  =>                  $parameters[$i]->{'name'} ,
	     value => $parameterValues{$parameters[$i]->{'name'}}
	 }
	 );
}

# Find the scratch directory.
my $scratchDirectory = $config->{'likelihood'}->{'workDirectory'}."/mcmc";
$scratchDirectory = $config->{'likelihood'}->{'scratchDirectory'}
    if ( exists($config->{'likelihood'}->{'scratchDirectory'}) );

# Expand any environment variable names in the scratch directory.
while ( $scratchDirectory =~ m/\$([_A-Z]+)/ ) {
    my $environmentVariableName  = $1;
    my $environmentVariableValue = $ENV{$environmentVariableName};
    $scratchDirectory =~ s/\$$environmentVariableName/$environmentVariableValue/g;
}

# Ensure scratch and work directories exist.
system("mkdir -p ".$config->{'likelihood'}->{'workDirectory'}."/mcmc")
    unless ( -e $config->{'likelihood'}->{'workDirectory'}."/mcmc" );
system("mkdir -p ".$scratchDirectory)
    unless ( -e $scratchDirectory );

# Report.
if ( exists($config->{'likelihood'}->{'report'}) ) {
    if ( $config->{'likelihood'}->{'report'} eq "yes" ) {
	print "Report from constrainGalacticus.pl:\n";
	print "  MPI rank is : ".$mpiRank."\n";
	print "  Output in   : ".$scratchDirectory."/newParameters_".$mpiRank.".xml\n";
	print "  Parameters  : \n";
	print Dumper($newParameters);
    }
}

# Find the set of base parameters to use.
my $baseParameters = "parameters.xml";
$baseParameters = $config->{'likelihood'}->{'baseParameters'}
   if ( exists($config->{'likelihood'}->{'baseParameters'}) );

# Find the CPU limit.
my $cpuLimit;
if ( exists($config->{'likelihood'}->{'cpulimit'}) ) {
    $cpuLimit = $config->{'likelihood'}->{'cpulimit'};
    if ( exists($config->{'likelihood'}->{'threads'}) ) { 
	$cpuLimit *= $config->{'likelihood'}->{'threads'};
    } else {
	$cpuLimit *= Sys::CPU::cpu_count();
    }
}

# Extract compilation file and project directories.
my $compilationFile  = $config->{'likelihood'}->{'compilation'  };
my $projectDirectory = $config->{'likelihood'}->{'workDirectory'};
my $galacticusFile   = "constrainGalacticus_".$mpiRank.".hdf5";

# Remove any old semaphore file.
unlink("/dev/shm/sem.galacticus")
    if ( -e "/dev/shm/sem.galacticus" );

# Bad log likelihood (highly improbable) which we will return in failure conditions.
my $badLogLikelihood = -1.0e30;

# Initialize a list of temporary files to remove after we're finished.
my @temporaryFiles;

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($compilationFile,$baseParameters);
my @constraints = @{$constraintsRef};

# Set an output file name.
$parameters->{'parameter'}->{'galacticusOutputFileName'}->{'value'} = $scratchDirectory."/".$galacticusFile;
push(@temporaryFiles,$parameters->{'parameter'}->{'galacticusOutputFileName'}->{'value'});

# Set state file names.
my $stateFileRoot;
if ( exists($config->{'likelihood'}->{'saveState'}) && $config->{'likelihood'}->{'saveState'} eq "yes" ) {
    $stateFileRoot = $scratchDirectory."/".$galacticusFile;
    $stateFileRoot =~ s/\.hdf5//;
    $parameters->{'parameter'}->{'stateFileRoot'}->{'value'} = $stateFileRoot;
    push(@temporaryFiles,$stateFileRoot."*");
}

# Set a random number seed.
$parameters->{'parameter'}->{'randomSeed'}->{'value'} = int(rand(10000))+1
    unless ( exists($config->{'likelihood'}->{'randomize'}) && $config->{'likelihood'}->{'randomize'} eq "no" );

# If running at a high temperature, modify the number of merger trees per decade.
my $temperatureEffective = 1.0;
if ( 
    exists($parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}) 
    &&
    $parameters->{'parameter'}->{'mergerTreeConstructMethod'}->{'value'} eq "build" 
    ) {
    my $treesPerDecadeEffective =
	max(
	    int(
		$parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'}
		/$temperature
	    )
	    ,$config->{'likelihood'}->{'treesPerDecadeMinimum'}
	);
    $temperatureEffective = $parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'}/$treesPerDecadeEffective;
    $parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'} = $treesPerDecadeEffective;
}

# If fixed sets of trees are to be used, create them as necessary, and store to a file.
if ( exists($config->{'likelihood'}->{'useFixedTrees'}) && $config->{'likelihood'}->{'useFixedTrees'} eq "yes" ) {
    # Record the required set of output redshifts.
    my $outputRedshifts = $parameters->{'parameter'}->{'outputRedshifts'   }->{'value'};
    # Record and remove any analyses.
    my $savedAnalyses   = $parameters->{'parameter'}->{'mergerTreeAnalyses'}->{'value'};
    delete($parameters->{'parameter'}->{'mergerTreeAnalyses'});
    # Get a lock on the tree file.
    my $fixedTreeDirectory;
    if ( exists($config->{'likelihood'}->{'fixedTreesInScratch'}) && $config->{'likelihood'}->{'fixedTreesInScratch'} eq "yes" ) {
	$fixedTreeDirectory = $config->{'likelihood'}->{'scratchDirectory'}."/";
    } else {
	$fixedTreeDirectory = $config->{'likelihood'}->{'workDirectory'   }."/";
    }
    my $fixedTreeFile      = $fixedTreeDirectory                           ."fixedTrees".$parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".hdf5";
    my $buildFixedTreeFile = $config->{'likelihood'}->{'workDirectory'}."/"."fixedTrees".$parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".hdf5";
    if ( 
	my $lock = new File::NFSLock {
	    file               => $fixedTreeFile,
	    lock_type          => LOCK_EX
	}
	)
    {
	unless ( -e $fixedTreeFile ) {
	    if ( 
		my $buildLock = new File::NFSLock {
		    file               => $buildFixedTreeFile,
		    lock_type          => LOCK_EX
		}
		)
	    {
		unless ( -e $buildFixedTreeFile ) {
		    # Create the tree file if necessary. (Set output redshift to a very large value to avoid any galaxy formation
		    # calculation being carried out - we only want to build the trees.)
		    $parameters->{'parameter'}->{'outputRedshifts'             }->{'value'} = "10000.0";
		    $parameters->{'parameter'}->{'mergerTreesWrite'            }->{'value'} = "true";
		    $parameters->{'parameter'}->{'mergerTreeExportFileName'    }->{'value'} = $buildFixedTreeFile;
		    $parameters->{'parameter'}->{'mergerTreeExportOutputFormat'}->{'value'} = "galacticus";
		    my $treeParameters;
		    push(@{$treeParameters->{'parameter'}},{name => $_, value => $parameters->{'parameter'}->{$_}->{'value'}})
			foreach ( keys(%{$parameters->{'parameter'}}) );
		    my $treeXML = new XML::Simple (RootName=>"parameters", NoAttr => 1);
		    open(pHndl,">".$config->{'likelihood'}->{'workDirectory'}."/treeBuildParameters".$parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".xml");
		    print pHndl $treeXML->XMLout($treeParameters);
		    close pHndl;		  
		    my $treeCommand;
		    $treeCommand .= "ulimit -t ".$cpuLimit."; "
			if ( defined($cpuLimit) );
		    $treeCommand .= "ulimit -c unlimited; GFORTRAN_ERROR_DUMPCORE=YES; ./Galacticus.exe ".$config->{'likelihood'}->{'workDirectory'}."/treeBuildParameters".$parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".xml";
		    my $treeLog = $config->{'likelihood'}->{'workDirectory'}."/treeBuildParameters".$parameters->{'parameter'}->{'mergerTreeBuildTreesPerDecade'}->{'value'}.".log";
		    SystemRedirect::tofile($treeCommand,$treeLog);
		    unless ( $? == 0 ) {
			system("mv ".$treeLog." ".$treeLog.".failed.".$$);
			die("constrainGalacticus.pl: Galacticus model failed");		
		    }
		    sleep(1)
			until ( -e $buildFixedTreeFile );
		}
		$buildLock->unlock();
	    }
	    system("cp -f ".$buildFixedTreeFile." ".$fixedTreeFile);
	}
	$lock->unlock();
    }
    # Modify parameters to use the tree file.
    $parameters->{'parameter'}->{'mergerTreeConstructMethod'}->{'value'} = "read";
    $parameters->{'parameter'}->{'mergerTreeReadFileName'   }->{'value'} = $fixedTreeFile;
    $parameters->{'parameter'}->{'outputRedshifts'          }->{'value'} = $outputRedshifts;
    $parameters->{'parameter'}->{'mergerTreeAnalyses'       }->{'value'} = $savedAnalyses;
    $parameters->{'parameter'}->{'mergerTreesWrite'         }->{'value'} = "false";
}

# Extract new parameters.
if ( defined($newParameters->{'parameter'}) ) {
    my @newParameterList;
    if ( ref($newParameters->{'parameter'}) eq "ARRAY" ) {
	@newParameterList = @{$newParameters->{'parameter'}};
    } else {
	push(@newParameterList,$newParameters->{'parameter'});
    }
    for my $newParameter ( @newParameterList ) {
	$parameters->{'parameter'}->{$newParameter->{'name'}}->{'value'} = $newParameter->{'value'};
    }
}

# Write the modified parameters to file.
&Parameters::Output($parameters,$scratchDirectory."/constrainGalacticusParameters".$mpiRank.".xml");
push(@temporaryFiles,$scratchDirectory."/constrainGalacticusParameters".$mpiRank.".xml");

# Run the Galacticus model.
my $glcCommand;
$glcCommand .= "ulimit -t ".$cpuLimit."; "
    if ( defined($cpuLimit) );
$glcCommand .= "export OMP_NUM_THREADS=".$config->{'likelihood'}->{'threads'}."; "
    if ( exists($config->{'likelihood'}->{'threads'}) );
$glcCommand .= "export ".$_."; "
    foreach ( &ExtraUtils::as_array($config->{'likelihood'}->{'environment'}) );
$glcCommand .= "ulimit -c unlimited; ./Galacticus.exe ".$scratchDirectory."/constrainGalacticusParameters".$mpiRank.".xml";
my $logFile = $scratchDirectory."/constrainGalacticusParameters".$mpiRank.".log";
push(@temporaryFiles,$logFile);
# my $timeGalacticusStart = [gettimeofday];
SystemRedirect::tofile($glcCommand,$logFile);
unless ( $? == 0 ) {
    # Issue a failure.
    print "ERROR: Galacticus model failed to complete\n";
    &reportFailure($config,$scratchDirectory,$logFile,$stateFileRoot);
    # Try running the model again - in case this was a random error.
    SystemRedirect::tofile($glcCommand,$logFile);
    unless ( $? == 0 ) {
	# Display the final likelihood.
	&outputLikelihood($config,$badLogLikelihood);
	print "constrainGalacticus.pl: Galacticus model failed";
	system("rm ".join(" ",@temporaryFiles))
	    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
	exit;
    }
}
# my $timeGalacticusElapsed = tv_interval($timeGalacticusStart,[gettimeofday]);
# print "%% timing : Galacticus : ".$timeGalacticusElapsed."\n";

# Perform processing of the model, accumulating likelihood as we go.
my $logLikelihood = 0.0;
my $xml           = new XML::Simple;
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $constraintDefinition;
    if ( -e $constraint->{'definition'}.".store" ) {
	$constraintDefinition = retrieve($constraint->{'definition'}.".store");
    } else {
	$constraintDefinition = $xml->XMLin($constraint->{'definition'});
    }
    # Run the analysis code.
    my $analysisCommand = $constraintDefinition->{'analysis'};
    $analysisCommand   .= " ".$scratchDirectory."/".$galacticusFile." --outputFile ".$scratchDirectory."/likelihood".$mpiRank.".xml";
    $analysisCommand .= " --temperature ".$temperatureEffective;
    $analysisCommand .= " --modelDiscrepancies ".$projectDirectory."/modelDiscrepancy"
	if ( -e $projectDirectory."/modelDiscrepancy" );
    $analysisCommand .= " --resultFile ".$scratchDirectory."/results".$mpiRank.".xml"
	if ( exists($config->{'likelihood'}->{'storeResults'}) && $config->{'likelihood'}->{'storeResults'} eq "yes" );
    system($analysisCommand);
    unless ( $? == 0 ) {
	# Issue a failure.
	print "ERROR: Analysis script failed to complete\n";
	&reportFailure($config,$scratchDirectory,$logFile,$stateFileRoot);
	# Display the final likelihood.
	&outputLikelihood($config,$badLogLikelihood);
	print "constrainGalacticus.pl: analysis code failed";
	system("rm ".join(" ",@temporaryFiles))
	    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
	exit;
    }
    # Store the results.
    if ( exists($config->{'likelihood'}->{'storeResults'}) && $config->{'likelihood'}->{'storeResults'} eq "yes" ) {
	my $results = $xml->XMLin($scratchDirectory."/results".$mpiRank.".xml");
	open(oHndl,">>".$config->{'likelihood'}->{'workDirectory'}."/mcmc/results".ucfirst($constraintDefinition->{'label'})."_".$mpiRank.".txt");
	for(my $i=0;$i<scalar(@{$results->{'y'}});++$i) {
	    print oHndl "\t"
		unless ( $i == 0 );
	    print oHndl ${$results->{'y'}}[$i]."\t".${$results->{'error'}}[$i];
	}
	print oHndl "\n";
	close(oHndl);
	unlink($scratchDirectory."/results".$mpiRank.".xml");
    }
    # Read the likelihood.
    my $likelihood = $xml->XMLin($scratchDirectory."/likelihood".$mpiRank.".xml");
    if ( $likelihood->{'logLikelihood'} eq "nan" ) {
	# Issue a failure.
        print "ERROR: Likelihood is NaN\n";
	&reportFailure($config,$scratchDirectory,$logFile,$stateFileRoot);
	# Display the final likelihood.
	&outputLikelihood($config,$badLogLikelihood);
	print "constrainGalacticus.pl: likelihood calculation failed";
	system("rm ".join(" ",@temporaryFiles))
	    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" && scalar(@temporaryFiles) > 0 );
	exit;
    }
    # Extract the likelihood and weight it.
    my $thisLogLikelihood = $likelihood->{'logLikelihood'};
    $thisLogLikelihood *= $constraintDefinition->{'weight'}
        if ( defined($constraintDefinition->{'weight'}) );
    # Accumulate the likelihood.
    $logLikelihood += $thisLogLikelihood;
    # Clean up.
    unlink($scratchDirectory."/likelihood".$mpiRank.".xml")
	if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} eq "yes" );
}

# Remove the model.
system("rm ".join(" ",@temporaryFiles))
    if ( exists($config->{'likelihood'}->{'cleanUp'}) && $config->{'likelihood'}->{'cleanUp'} && scalar(@temporaryFiles) > 0 );

# Display the final likelihood.
&outputLikelihood($config,$logLikelihood);

# Script timing.
# my $timeElapsed = tv_interval($timeStart,[gettimeofday]);
# print "%% timing : ".$self." : ".$timeElapsed."\n";
exit;

sub outputLikelihood {
    # Output the log likelihood.
    my $config        = shift;
    my $logLikelihood = shift;
    open(oHndl,">".$likelihoodFile);
    print oHndl $logLikelihood."\n";
    close(oHndl);
}

sub reportFailure {
    # Handle failures of model or analysis.
    my $config           = shift;
    my $scratchDirectory = shift;
    my $logFile       = shift;
    my $stateFileRoot = shift;
    if ( exists($config->{'likelihood'}->{'failArchive'}) ) {
	my $failArchiveName = $config->{'likelihood'}->{'failArchive'}.".tar.bz2";
	if ( ! -e $failArchiveName && -e "galacticusConfig.xml" ) {
	    # Send an email if possible.
	    my $xml     = new XML::Simple;
	    my $config  = $xml->XMLin("galacticusConfig.xml");
	    if ( exists($config->{'contact'}->{'email'}) && &File::Which::which('sendmail') ) {
		my $message  = "A Galacticus model failed while seeking constraints.\n";
		$message    .= "Failed model is in: ".$failArchiveName."\n";
		my $msg = MIME::Lite->new(
		    From    => '',
		    To      => $config->{'contact'}->{'email'},
		    Subject => 'Galacticus model failed while seeking constraints',
		    Type    => 'TEXT',
		    Data    => $message
		    );
		$msg->send();
	    } else {
		open(my $eHndl,">constrainGalacticusError.txt");
		print $eHndl "A Galacticus model failed while seeking constraints.\n";
		print $eHndl "Failed model is in: ".$failArchiveName."\n";
		close($eHndl);
	    }
	}
	if ( ! -e $failArchiveName ) {
	    system("touch ".$failArchiveName);
	    my $tarCommand = "tar cvfj ".$failArchiveName." ".$scratchDirectory."/constrainGalacticusParameters".$mpiRank.".xml ".$logFile." ".$scratchDirectory."/".$galacticusFile;
	    opendir(dHndl,".");
	    while ( my $fileName = readdir(dHndl) ) {
		$tarCommand .= " ".$fileName
		    if ( $fileName eq "core.".$$ );
	    }
	    closedir(dHndl);
	    $tarCommand .= " ".$stateFileRoot.".*state*"
		if ( $config->{'likelihood'}->{'saveState'} eq "yes" );
	    system($tarCommand);
	}
    }
    if ( exists($config->{'likelihood'}->{'failCount'}) ) {
	my $count = 0;
	if ( -e $config->{'likelihood'}->{'failCount'} ) {
	    open(iHndl,$config->{'likelihood'}->{'failCount'});
	    $count = <iHndl>;
	    close(iHndl);
	    if (defined($count) ) {
		chomp($count);
	    } else {
		$count = 0;
	    }
	}
	++$count;
	open(oHndl,">".$config->{'likelihood'}->{'failCount'});
	print oHndl $count."\n";
	close(oHndl);
    }
}
