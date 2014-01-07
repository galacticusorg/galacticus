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
require File::Which;
require File::NFSLock;
require System::Redirect;
require Galacticus::Constraints::Parameters;
require List::ExtraUtils;

# Finds constraints on Galacticus parameters.
# Andrew Benson (02-September-2011)

# Get command line arguments.
die("Usage: constrainGalacticus.pl <configFile> <workDirectory> <compilationFile> <parameterFile> <projectDirectory> [options]")
    unless ( scalar(@ARGV) >= 5 );
my $configFile       = $ARGV[0];
my $workDirectory    = $ARGV[1];
my $compilationFile  = $ARGV[2];
my $parameterFile    = $ARGV[3];
my $projectDirectory = $ARGV[4];

# Parse the config file.
my $config     = &Parameters::Parse_Config($configFile);

# Create a hash of named arguments.
my $iArg = -1;
my %arguments = (
    make                => "no",
    timing              => "no",
    output              => "stdout",
    galacticusFile      => "constrainGalacticus.hdf5",
    saveState           => "no",
    threads             => 1,
    reuseTrees          => "",
    suffix              => "",
    cleanUp             => "yes",
    randomize           => "yes",
    temperature         => 1.0
    );
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Remove any old semaphore file.
unlink("/dev/shm/sem.galacticus")
    if ( -e "/dev/shm/sem.galacticus" );

# Bad log likelihood (highly improbable) which we will return in failure conditions.
my $badLogLikelihood = -1.0e30;

# Get a hash of the parameter values.
my $baseParametersFile = "constraints/baseParameters.xml";
$baseParametersFile    = $arguments{'baseParameters'}
    if ( defined($arguments{'baseParameters'}) );
(my $constraintsRef, my $parameters) = &Parameters::Compilation($compilationFile,$baseParametersFile);
my @constraints = @{$constraintsRef};

# Set an output file name.
$parameters->{'parameter'}->{'galacticusOutputFileName'}->{'value'} = $workDirectory."/".$arguments{'galacticusFile'};

# Set state file names.
my $stateFileRoot;
if ( $arguments{'saveState'} eq "yes" ) {
    $stateFileRoot = $workDirectory."/".$arguments{'galacticusFile'};
    $stateFileRoot =~ s/\.hdf5//;
    $parameters->{'parameter'}->{'stateFileRoot'}->{'value'} = $stateFileRoot;
}

# Set a random number seed.
$parameters->{'parameter'}->{'randomSeed'}->{'value'} = int(rand(10000))+1
    unless ( $arguments{'randomize'} eq "no" );

# Ensure that timing data is collected.
$parameters->{'parameter'}->{'metaCollectTimingData'}->{'value'} = "true"
    if ( $arguments{'timing'} eq "yes" );

# Parse the modifications to the parameters.
my $xml              = new XML::Simple;
my $newParameterData = $xml->XMLin($parameterFile, KeyAttr => "");
if ( defined($newParameterData->{'parameter'}) ) {
    my @newParameters;
    if ( ref($newParameterData->{'parameter'}) eq "ARRAY" ) {
	@newParameters = @{$newParameterData->{'parameter'}};
    } else {
	push(@newParameters,$newParameterData->{'parameter'});
    }
    for my $newParameter ( @newParameters ) {
	$parameters->{'parameter'}->{$newParameter->{'name'}}->{'value'} = $newParameter->{'value'};
    }
}

# If trees are to be reused, create them, and store to a file.
unless ( $arguments{'reuseTrees'} eq "" ) {
    # Record the required set of output redshifts.
    my $outputRedshifts = $parameters->{'parameter'}->{'outputRedshifts'}->{'value'};
    # Record and remove any analyses.
    my $savedAnalyses = $parameters->{'parameter'}->{'mergerTreeAnalyses'}->{'value'};
    delete($parameters->{'parameter'}->{'mergerTreeAnalyses'});
    # Get a lock on the tree file.
    if ( 
	my $lock = new File::NFSLock {
	    file               => $arguments{'reuseTrees'},
	    lock_type          => LOCK_EX
	}
	)
    {
	unless ( -e $arguments{'reuseTrees'} ) {
	    # Create the tree file if necessary. (Set output redshift to a very large value to avoid any galaxy formation
	    # calculation being carried out - we only want to build the trees.)
	    $parameters->{'parameter'}->{'outputRedshifts'             }->{'value'} = "10000.0";
	    $parameters->{'parameter'}->{'mergerTreesWrite'            }->{'value'} = "true";
	    $parameters->{'parameter'}->{'mergerTreeExportFileName'    }->{'value'} = $arguments{'reuseTrees'};
	    $parameters->{'parameter'}->{'mergerTreeExportOutputFormat'}->{'value'} = "galacticus";
	    $parameters->{'parameter'}->{'mergerTreeExportOutputFormat'}->{'value'} = "galacticus";
	    my $treeParameters;
	    push(@{$treeParameters->{'parameter'}},{name => $_, value => $parameters->{'parameter'}->{$_}->{'value'}})
		foreach ( keys(%{$parameters->{'parameter'}}) );
	    my $treeXML = new XML::Simple (RootName=>"parameters", NoAttr => 1);
	    open(pHndl,">".$workDirectory."/treeBuildParameters.xml");
	    print pHndl $treeXML->XMLout($treeParameters);
	    close pHndl;
	    if ( $arguments{'make'} eq "yes" ) {
		system("make Galacticus.exe");
		die("constrainGalacticus.pl: failed to build Galacticus.exe") unless ( $? == 0 );
	    }
	    my $treeCommand;
	    $treeCommand .= "ulimit -t ".$arguments{'cpulimit'}."; " if ( exists($arguments{'cpulimit'}) );
	    $treeCommand .= "ulimit -c unlimited; GFORTRAN_ERROR_DUMPCORE=YES; ./Galacticus.exe ".$workDirectory."/treeBuildParameters.xml";
	    my $treeLog = $workDirectory."/treeBuildParameters.log";
	    SystemRedirect::tofile($treeCommand,$treeLog);
	    unless ( $? == 0 ) {
		system("mv ".$treeLog." ".$treeLog.".failed.".$$);
		die("constrainGalacticus.pl: Galacticus model failed");
	    }
	}
	$lock->unlock();
    }
    # Modify parameters to use the tree file.
    $parameters->{'parameter'}->{'mergerTreeConstructMethod'}->{'value'} = "read";
    $parameters->{'parameter'}->{'mergerTreeReadFileName'   }->{'value'} = $arguments{'reuseTrees'};
    $parameters->{'parameter'}->{'outputRedshifts'          }->{'value'} = $outputRedshifts;
    $parameters->{'parameter'}->{'mergerTreeAnalyses'       }->{'value'} = $savedAnalyses;
}

# Write the modified parameters to file.
&Parameters::Output($parameters,$workDirectory."/constrainGalacticusParameters".$arguments{'suffix'}.".xml");

# Run the Galacticus model.
if ( $arguments{'make'} eq "yes" ) {
    system("make Galacticus.exe");
    die("constrainGalacticus.pl: failed to build Galacticus.exe") unless ( $? == 0 );
}
my $glcCommand;
$glcCommand .= "ulimit -t ".$arguments{'cpulimit'}."; "
    if ( exists($arguments{'cpulimit'}) );
$glcCommand .= "export OMP_NUM_THREADS=".$arguments{'threads'}."; "
    if ( exists($arguments{'threads'}) );
$glcCommand .= "export ".$_."; "
    foreach ( &ExtraUtils::as_array($config->{'likelihood'}->{'environment'}) );
$glcCommand .= "ulimit -c unlimited; ./Galacticus.exe ".$workDirectory."/constrainGalacticusParameters".$arguments{'suffix'}.".xml";
my $logFile = $workDirectory."/constrainGalacticusParameters".$arguments{'suffix'}.".log";
SystemRedirect::tofile($glcCommand,$logFile);
unless ( $? == 0 ) {
    # Issue a failure.
    print "ERROR: Galacticus model failed to complete\n";
    &reportFailure(\%arguments,$workDirectory,$logFile,$stateFileRoot);
    # Try running the model again - in case this was a random error.
    SystemRedirect::tofile($glcCommand,$logFile);
    unless ( $? == 0 ) {
	# Display the final likelihood.
	&outputLikelihood(\%arguments,$badLogLikelihood);
	print "constrainGalacticus.pl: Galacticus model failed";
	exit;
    }
}

# Perform processing of the model, accumulating likelihood as we go.
my $logLikelihood = 0.0;
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Run the analysis code.
    my $analysisCommand = $constraintDefinition->{'analysis'};
    $analysisCommand   .= " ".$workDirectory."/".$arguments{'galacticusFile'}." --outputFile ".$workDirectory."/likelihood".$arguments{'suffix'}.".xml";
    $analysisCommand .= " --modelDiscrepancies ".$projectDirectory."/modelDiscrepancy"
	if ( -e $projectDirectory."/modelDiscrepancy" );
    $analysisCommand .= " --resultFile ".$workDirectory."/results".$arguments{'suffix'}.".xml"
	if ( exists($arguments{'storeResults'}) );
    system($analysisCommand);
    unless ( $? == 0 ) {
	# Issue a failure.
	print "ERROR: Analysis script failed to complete\n";
	&reportFailure(\%arguments,$workDirectory,$logFile,$stateFileRoot);
	# Display the final likelihood.
	&outputLikelihood(\%arguments,$badLogLikelihood);
	print "constrainGalacticus.pl: analysis code failed";
	exit;
    }
    # Store the results.
    if ( exists($arguments{'storeResults'}) ) {
	my $results = $xml->XMLin($workDirectory."/results".$arguments{'suffix'}.".xml");
	open(oHndl,">>".$arguments{'storeResults'}."/results".ucfirst($constraintDefinition->{'label'})."_".$arguments{'suffix'}.".txt");
	for(my $i=0;$i<scalar(@{$results->{'y'}});++$i) {
	    print oHndl "\t"
		unless ( $i == 0 );
	    print oHndl ${$results->{'y'}}[$i]."\t".${$results->{'error'}}[$i];
	}
	print oHndl "\n";
	close(oHndl);
	unlink($workDirectory."/results".$arguments{'suffix'}.".xml");
    }
    # Read the likelihood.
    my $likelihood = $xml->XMLin($workDirectory."/likelihood".$arguments{'suffix'}.".xml");
    if ( $likelihood->{'logLikelihood'} eq "nan" ) {
	# Issue a failure.
        print "ERROR: Likelihood is NaN\n";
	&reportFailure(\%arguments,$workDirectory,$logFile,$stateFileRoot);
	# Display the final likelihood.
	&outputLikelihood(\%arguments,$badLogLikelihood);
	print "constrainGalacticus.pl: likelihood calculation failed";
	exit;
    }
    # Extract the likelihood and weight it.
    my $thisLogLikelihood = $likelihood->{'logLikelihood'};
    $thisLogLikelihood *= $constraintDefinition->{'weight'}
        if ( defined($constraintDefinition->{'weight'}) );
    # Accumulate the likelihood.
    $logLikelihood += $thisLogLikelihood;
    # Clean up.
    unlink($workDirectory."/likelihood".$arguments{'suffix'}.".xml")
	if ( $arguments{'cleanUp'} eq "yes" );
}

# Extract tree timing information.
system("scripts/analysis/treeTiming.pl ".$workDirectory."/".$arguments{'galacticusFile'}." --maxPoints 10000 --outputFile ".$workDirectory."/constrainGalacticusTiming.xml --accumulate")
    if ( $arguments{'timing'} eq "yes" );

# Remove the model.
unlink($workDirectory."/".$arguments{'galacticusFile'})
    if ( $arguments{'cleanUp'} eq "yes" );

# Adjust likelihood for temperature.
$logLikelihood /= $arguments{'temperature'};

# Display the final likelihood.
&outputLikelihood(\%arguments,$logLikelihood);

exit;

sub outputLikelihood {
    # Output the log likelihood.
    my %arguments     = %{shift()};
    my $logLikelihood =   shift   ;
    if ( $arguments{'output'} eq "stdout" ) {
	print $logLikelihood."\n";
    } else {
	open(oHndl,">".$arguments{'output'});
	print oHndl $logLikelihood."\n";
	close(oHndl);
    }
}

sub reportFailure {
    # Handle failures of model or analysis.
    my %arguments     = %{shift()};
    my $workDirectory = shift;
    my $logFile       = shift;
    my $stateFileRoot = shift;
    if ( exists($arguments{'failArchive'}) ) {
	my $failArchiveName = $arguments{'failArchive'}.".tar.bz2";
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
	    my $tarCommand = "tar cvfj ".$failArchiveName." ".$workDirectory."/constrainGalacticusParameters".$arguments{'suffix'}.".xml ".$logFile." ".$workDirectory."/".$arguments{'galacticusFile'};
	    opendir(dHndl,".");
	    while ( my $fileName = readdir(dHndl) ) {
		$tarCommand .= " ".$fileName
		    if ( $fileName eq "core.".$$ );
	    }
	    closedir(dHndl);
	    $tarCommand .= " ".$stateFileRoot.".*state*"
		if ( $arguments{'saveState'} eq "yes" );
	    system($tarCommand);
	}
    }
    if ( exists($arguments{'failCount'}) ) {
	my $count = 0;
	if ( -e $arguments{'failCount'} ) {
	    open(iHndl,$arguments{'failCount'});
	    $count = <iHndl>;
	    close(iHndl);
	    if (defined($count) ) {
		chomp($count);
	    } else {
		$count = 0;
	    }
	}
	++$count;
	open(oHndl,">".$arguments{'failCount'});
	print oHndl $count."\n";
	close(oHndl);
    }
}
