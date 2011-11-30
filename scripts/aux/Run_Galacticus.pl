#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use Thread;
use XML::Simple;
use Data::Dumper;
use Data::Compare;
use File::Copy;
use File::Slurp qw( slurp );
require System::Redirect;
use Sys::CPU;
use MIME::Lite;
require IO::Compress::Simple;
use Switch;

# Script to run sets of Galacticus models, looping through sets of parameters and performing analysis
# on the results. Supports launching of multiple threads (each with a different model) on a single
# machine and submission of jobs to a Condor cluster. Contains error reporting functionality.
# Andrew Benson (11-June-2010)

# Get command line arguments.
if ( $#ARGV < 0 ) {die("Usage: Run_Galacticus.pl <runFile>")};
my $runFile = $ARGV[0];

# Create a hash of named arguments.
my $iArg = -1;
my %arguments;
while ( $iArg < $#ARGV ) {
    ++$iArg;
    if ( $ARGV[$iArg] =~ m/^\-\-(.*)/ ) {
	$arguments{$1} = $ARGV[$iArg+1];
	++$iArg;
    }
}

# Check for an instance number for this launch.
my $thisInstance  = 1;
my $instanceCount = 1;
if ( exists($arguments{"instance"}) && $arguments{"instance"} =~ m/(\d+):(\d+)/ ) {
    $thisInstance  = $1;
    $instanceCount = $2;
}

# Read in the file of models to be run.
my $xml         = new XML::Simple;
my $modelsToRun = $xml->XMLin($runFile, KeyAttr => "", ForceArray => [ "value", "parameter", "parameters", "requirement" ]);

# Read in any configuration options.
my $config;
if ( -e "galacticusConfig.xml" ) {
    my $xml = new XML::Simple;
    $config = $xml->XMLin("galacticusConfig.xml");
}

# Determine root directory for models.
my $rootDirectory = "models";
if ( exists($modelsToRun->{'modelRootDirectory'}) ) {$rootDirectory  = $modelsToRun->{'modelRootDirectory'}};

# Determine the base set of parameters to use.
my $baseParameters = "";
if ( exists($modelsToRun->{'baseParameters'    }) ) {$baseParameters = $modelsToRun->{'baseParameters'    }};

# Determine how many threads to launch.
my $threadCount = 1;
if ( exists($modelsToRun->{'threadCount'       }) ) {
    $threadCount = $modelsToRun->{'threadCount'};
    $threadCount = Sys::CPU::cpu_count() if ( $threadCount eq "maximum" );
}

# Determine how to run the models.
my $launchMethod = "local";
$launchMethod    = "condor"
    if (
	exists($modelsToRun->{'condor'})
	&& exists($modelsToRun->{'condor'}->{'useCondor'})
	&&        $modelsToRun->{'condor'}->{'useCondor'} eq "true"
	);

# Record where we are running.
my $pwd = `pwd`;
chomp($pwd);

# Define a hash for Condor job IDs.
my %condorJobs;

# Launch threads.
if ( $threadCount <= 1 || $launchMethod eq "condor" ) {
    &Launch_Models($modelsToRun);
} else {
    my @threads = ();
    for(my $iThread=0;$iThread<$threadCount;++$iThread) {
	$threads[++$#threads] = new Thread \&Launch_Models, $modelsToRun;
	sleep 1; # Adds a pause to minimize file race conditions.
    }
    foreach my $thread ( @threads ) {
	my @returnData = $thread->join;
    }
}

# Wait for Condor models to finish.
if ( $launchMethod eq "condor" ) {
    print "Waiting for Condor jobs to finish...\n";
    while ( scalar(keys %condorJobs) > 0 ) {
	# Find all Condor jobs that are running.
	my %runningCondorJobs;
	undef(%runningCondorJobs);
	my $condorQuerySuccess = 0;
	open(pHndl,"condor_q|");
	while ( my $line = <pHndl> ) {
	    if ( $line =~ m/^\s*(\d+\.\d+)\s/ ) {$runningCondorJobs{$1} = 1};
	    if ( $line =~ m/^\d+\s+jobs;/ ) {$condorQuerySuccess = 1};
	}
	close(pHndl);
	if ( $condorQuerySuccess == 1 ) {
	    foreach my $jobID ( keys(%condorJobs) ) {
		unless ( exists($runningCondorJobs{$jobID}) ) {
		    print "Condor job ".$jobID." has finished. Post-processing....\n";
		    &Model_Finalize(
				    $condorJobs{$jobID}->{'directory'},
				    $condorJobs{$jobID}->{'directory'}."/galacticus.hdf5",
				    0
				    );
		    # Remove any temporary files associated with this job.
		    unlink(@{$condorJobs{$jobID}->{'temporaryFiles'}});
		    # Remove the job ID from the list of active Condor jobs.
		    delete($condorJobs{$jobID});
		}
	    }
	    sleep 5;
	} else {
	    sleep 10;
	}
    }
}

exit;

sub Launch_Models {
    # Grab input arguments.
    my $modelsToRun = shift;

    # Set initial value of random seed.
    my $randomSeed = 219;

    # Initialize a model counter.
    my $modelCounter = -1;

    # Loop through all model sets.
    my $iModelSet = -1;
    foreach my $parameterSet ( @{$modelsToRun->{'parameters'}} ) {
	# Increment model set counter.
	++$iModelSet;
	
	# Set base model name.
	my $modelBaseName = "galacticus";
	$modelBaseName = ${$parameterSet->{'label'}}[0] if ( exists($parameterSet->{'label'}) );
	
	# Create an array of hashes giving the parameters for this parameter set.
	my @parameterHashes = &Create_Parameter_Hashes($parameterSet);
	
	# Loop over all models and run them.
	my $iModel = 0;
	foreach my $parameterData ( @parameterHashes ) {
	    ++$iModel;
	    ++$modelCounter;
	    
	    # Increment random seed.
	    ++$randomSeed;

	    # Evaluate if this model should be launched.
	    my $runOnInstance = ($modelCounter % $instanceCount) + 1;
	    
	    # Specify the output directory.
	    my $galacticusOutputDirectory = $rootDirectory."/".$modelBaseName."_".$iModelSet.":".$iModel;
	    $galacticusOutputDirectory .= "_".$parameterData->{'label'} if ( exists($parameterData->{'label'}) );
	    # If the output directory does not exist, then create it.
	    unless ( -e $galacticusOutputDirectory || $runOnInstance != $thisInstance ) {		
		system("mkdir -p ".$galacticusOutputDirectory);
		
		# Specify the output file.
		my $galacticusOutputFile = $galacticusOutputDirectory."/galacticus.hdf5";
		
		# Read the default set of parameters.
		my %parameters;
		unless ( $baseParameters eq "" ) {
		    my $xml = new XML::Simple;
		    my $data = $xml->XMLin($baseParameters,ForceArray => 1);
		    my @parameterArray = @{$data->{'parameter'}};
		    for(my $i=0;$i<=$#parameterArray;++$i) {
			$parameters{${$parameterArray[$i]->{'name'}}[0]} = ${$parameterArray[$i]->{'value'}}[0];
		    }
		}
		
		# Set the output file name.
		switch ( $launchMethod ) {
		    case ( "local"  ) {
			$parameters{'galacticusOutputFileName'} = $galacticusOutputFile;
		    }
		    case ( "condor" ) {
			$parameters{'galacticusOutputFileName'} = "galacticus.hdf5";
		    }
		}

		# Set the random seed.
		$parameters{'randomSeed'} = $randomSeed unless ( exists($parameters{'randomSeed'}) );
		
		# Set a state restore file.
		(my $stateFile = $parameters{'galacticusOutputFileName'}) =~ s/\.hdf5//;
		$parameters{'stateFileRoot'}      = $stateFile;
		
		# Transfer parameters for this model from the array of model parameter hashes to the active hash.
		foreach my $parameter ( keys(%{$parameterData}) ) {
		    $parameters{$parameter} = ${$parameterData}{$parameter};
		}
		
		# Transfer values from the active hash to an array suitable for XML output.
		my $data;
		my @parameterArray;
		undef($data);
		undef(@parameterArray);
		foreach my $name ( sort(keys(%parameters)) ) {
		    unless ( $name eq "label" ) {
			$parameterArray[++$#parameterArray]->{'name'}  = $name;
			$parameterArray[  $#parameterArray]->{'value'} = $parameters{$name};
		    }
		}
		$data->{'parameter'} = \@parameterArray;
		
		# Output the parameters as an XML file.
		my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
		open(outHndl,">".$galacticusOutputDirectory."/newParameters.xml");
		print outHndl $xmlOutput->XMLout($data);
		close(outHndl);
		undef($data);
		undef(%parameters);

		# Launch the model.
		switch ( $launchMethod ) {
		    case ( "local" ) {
			&SystemRedirect::tofile("ulimit -t unlimited; ulimit -c unlimited; GFORTRAN_ERROR_DUMPCORE=YES; time Galacticus.exe "
						.$galacticusOutputDirectory."/newParameters.xml",$galacticusOutputDirectory."/galacticus.log");
			my $exitStatus = $?;
			&Model_Finalize($galacticusOutputDirectory,$galacticusOutputFile,$exitStatus);
		    }
		    case ( "condor" ) {
			# Extract Condor options.
			my $condorGalacticusDirectory = "/home/condor/Galacticus/v0.9.0";
			$condorGalacticusDirectory    =   $modelsToRun->{'condor'}->{'galacticusDirectory'}  if ( exists($modelsToRun->{'condor'}->{'galacticusDirectory'}) );
			my $condorUniverse            = "vanilla";
			$condorUniverse               =   $modelsToRun->{'condor'}->{'universe'           }  if ( exists($modelsToRun->{'condor'}->{'universe'           }) );
			my $condorEnvironment         = "";
			$condorEnvironment            =   $modelsToRun->{'condor'}->{'environment'        }  if ( exists($modelsToRun->{'condor'}->{'environment'        }) );
			my @condorRequirements;
			@condorRequirements           = @{$modelsToRun->{'condor'}->{'requirement'        }} if ( exists($modelsToRun->{'condor'}->{'requirement'        }) );
			my @condorTransferFiles       = ( $pwd."/Galacticus.exe", "newParameters.xml" );
			if ( exists($modelsToRun->{'condor'}->{'transferFile'}) ) {
			    if ( UNIVERSAL::isa($modelsToRun->{'condor'}->{'transferFile'},"ARRAY") ) {
				push(@condorTransferFiles,@{$modelsToRun->{'condor'}->{'transferFile'}});
			    } else {
				push(@condorTransferFiles,$modelsToRun->{'condor'}->{'transferFile'});             
			    }
			    foreach ( @condorTransferFiles ) {
				$_ =~ s/{PWD}/$pwd/;
			    }
			}

			# Create the script that Condor will execute.
			my $condorScript = "condor_run_".$modelCounter."_".$$.".csh";
			open(oHndl,">".$condorScript);
			print oHndl "#!/bin/csh\n";
			print oHndl "ln -sf ".$condorGalacticusDirectory."/aux\n";
			print oHndl "ln -sf ".$condorGalacticusDirectory."/data\n";
			print oHndl "ln -sf ".$condorGalacticusDirectory."/perl\n";
			print oHndl "ln -sf ".$condorGalacticusDirectory."/scripts\n";
			print oHndl "ln -sf ".$condorGalacticusDirectory."/work\n";
			print oHndl "exec ./Galacticus.exe newParameters.xml\n";
			print oHndl "rm -f aux data perl scripts work\n";
			print oHndl "exit\n";
			close(oHndl);
			# Creat a submit file.
			my $condorSubmit = "condor_submit_".$modelCounter."_".$$.".txt";
			open(oHndl,">".$condorSubmit);
			print oHndl "Executable              = ".$condorScript."\n";
			print oHndl "Error                   = condor.err\n";
			print oHndl "Output                  = condor.out\n";
			print oHndl "Log                     = condor.log\n";
			print oHndl "InitialDir              = ".$galacticusOutputDirectory."\n";
			print oHndl "Should_Transfer_Files   = YES\n";
			print oHndl "When_To_Transfer_Output = ON_EXIT\n";
			print oHndl "Transfer_Input_Files    = ".join(",",@condorTransferFiles)."\n";
			print oHndl "Universe                = ".$condorUniverse."\n";
			print oHndl "Allow_Startup_Script    = True\n" if ( $condorUniverse eq "standard" );
			print oHndl "Environment             = ".$condorEnvironment."\n" unless ( $condorEnvironment eq "" );
			print oHndl "Requirements            = (".join(") && (",@condorRequirements).")\n" if ( $#condorRequirements >= 0 );
			print oHndl "+RequiresWholeMachine   = True\n" if ( $modelsToRun->{'condor'}->{'wholeMachine'} eq "true" );
			print oHndl "Queue\n";
			close(oHndl);

			# Submit the job - capture the job number?
			open(pHndl,"condor_submit -verbose ".$condorSubmit."|");
			my $jobID = "";
			while ( my $line = <pHndl> ) {
			    if ( $line =~ m/^\*\*\s+Proc\s+([\d\.]+):/ ) {$jobID = $1};
			}
			close(pHndl);

                        # Add job number to active job hash
			unless ( $jobID eq "" ) {
			    $condorJobs{$jobID}->{'directory'} = $galacticusOutputDirectory;
			    @{$condorJobs{$jobID}->{'temporaryFiles'}} = [$condorScript,$condorSubmit];
			}
			sleep 10;

		    }
		}

	    }

	}

    }

}

sub Model_Finalize {
    # Finalize a model by performing analysis and compressing.
    my ($galacticusOutputDirectory,$galacticusOutputFile,$exitStatus) = @_;
    
    if ( $exitStatus == 0 ) {
	# Model finished successfully.
	# Generate plots.
	unless ( $modelsToRun->{'doAnalysis'} eq "no" ) {
	    my $analysisScript = $galacticusPath."data/Galacticus_Compute_Fit_Analyses.xml";
	    $analysisScript = $modelsToRun->{'analysisScript'} if ( exists($modelsToRun->{'analysisScript'}) );
	    system("./scripts/analysis/Galacticus_Compute_Fit.pl ".$galacticusOutputFile." ".$galacticusOutputDirectory." ".$analysisScript);
	}
    } else {
	# The run failed for some reason.
	# Move the core file to the output directory.
	opendir(gDir,".");
	while ( my $file = readdir(gDir) ) {
	    if ( $file =~ m/core\.\d+/ ) {move($file,$galacticusOutputDirectory."/core")};
	}
	closedir(gDir);
	# Report the model failure (by e-mail if we have an e-mail address to send a report to and if so requested).
	my $message  = "FAILED: A Galacticus model failed to finish:\n\n";
	$message .= "  Host:\t".$ENV{"HOSTNAME"}."\n";
	$message .= "  User:\t".$ENV{"USER"}."\n\n";
	$message .= "Model output is in: ".$pwd."/".$galacticusOutputDirectory."\n\n";
	if ( exists($config->{'contact'}->{'email'}) && $config->{'contact'}->{'email'} =~ m/\@/ && $modelsToRun->{'emailReport'} eq "yes" ) {
	    $message .= "Log file is attached.\n";
	    my $msg = MIME::Lite->new(
		From    => '',
		To      => $config->{'contact'}->{'email'},
		Subject => 'Galacticus model failed',
		Type    => 'TEXT',
		Data    => $message
		);
	    $msg->attach(
		Type     => "text/plain",
		Path     => $galacticusOutputDirectory."/galacticus.log",
		Filename => "galacticus.log"
		);
	    $msg->send;
	} else {
	    print $message;
	    print "Log follows:\n";
	    print slurp($galacticusOutputDirectory."/galacticus.log");
	}
    }
    
    # Compress all files in the output directory.
    &Simple::Compress_Directory($galacticusOutputDirectory);
    
}

sub Create_Parameter_Hashes {
    # Create an array of hashes which give the parameter values for each model.

    # Get the input parameters structure.
    my $parameterSet = shift;

    # Initalizes hash which record which parameter combinations have already been processed.
    my %doneCases;

    # Initalize array of hashes that we will return.
    my @parameterHashes;

    # Create a list of all parameters in the structure.
    my @parameterPointer;
    my $parameterID        = 0;
    my $processedParameter = -1;

    # Set the first parameter pointer to point to the input object.
    $parameterPointer[0] = $parameterSet;
    while ( $processedParameter < $#parameterPointer ) {
	# Move to the next parameter to process.
	++$processedParameter;
	# Add the array of parameters to the list.
	my @thisParameterArrays;
	if ( exists($parameterPointer[$processedParameter]->{'parameter'}) ) {
	    $thisParameterArrays[0] = $parameterPointer[$processedParameter]->{'parameter'};
	} elsif ( exists($parameterPointer[$processedParameter]->{'value'}) ) {
	    foreach my $valueElement ( @{$parameterPointer[$processedParameter]->{'value'}} ) {
		$thisParameterArrays[++$#thisParameterArrays] = $valueElement->{'parameter'}
		if ( UNIVERSAL::isa($valueElement,"HASH") && exists($valueElement->{'parameter'}) );
	    }
	}
	if ( @thisParameterArrays ) {
	    # Add any detected parameter arrays to the list.
	    foreach my $thisParameterArray ( @thisParameterArrays ) {
		# Loop over each parameter in the array.
		foreach my $parameter ( @{$thisParameterArray} ) {
		    # Store a pointer to the parameter.
		    $parameterPointer[++$#parameterPointer] = $parameter;
		    # Set an index counter (for its <value> elements) to zero.
		    $parameter->{'index'} = 0;
		    # Set a unique ID for this parameter.
		    ++$parameterID;
		    $parameter->{'ID'}    = $parameterID;
		}
	    }
	}
    }

    # Step through parameter combinations until all have been done.
    my $done = 0;
    until ( $done == 1 ) {
	# Build parameter hash.
	my @currentParameterPointer;
	my $currentProcessedParameter = -1;

	# Build a list of currently active parameters (i.e. those which are defined at the base level or within a currently used <value> element).
	$currentParameterPointer[0] = $parameterSet;
	while ( $currentProcessedParameter < $#currentParameterPointer ) {
	    ++$currentProcessedParameter;
	    # Add the array of parameters to the list.
	    my $thisParameterArray;
	    if ( exists($currentParameterPointer[$currentProcessedParameter]->{'parameter'}) ) {
		$thisParameterArray = $currentParameterPointer[$currentProcessedParameter]->{'parameter'};
	    } elsif ( exists($currentParameterPointer[$currentProcessedParameter]->{'value'}) ) {
		my $valueElement = ${$currentParameterPointer[$currentProcessedParameter]->{'value'}}[$currentParameterPointer[$currentProcessedParameter]->{'index'}];
		$thisParameterArray = $valueElement->{'parameter'} if ( UNIVERSAL::isa( $valueElement, "HASH" ) && exists($valueElement->{'parameter'}) );
	    }
	    if ( defined($thisParameterArray) ) {
		my $iParameter = -1;
		foreach my $parameter ( @{$thisParameterArray} ) {
		    ++$iParameter;
		    $currentParameterPointer[++$#currentParameterPointer] = $parameter;
		}
	    }
	}

	# Build a parameter hash including only active parameters.
	my %parameters;
	my $label  = "";
	my $joiner = "";
	foreach my $parameter ( @currentParameterPointer ) {
	    if ( exists($parameter->{'name'}) ) {
		my $value;
		if ( UNIVERSAL::isa(${$parameter->{'value'}}[$parameter->{'index'}], "HASH" ) && exists(${$parameter->{'value'}}[$parameter->{'index'}]->{'content'}) ) {
		    $value = ${$parameter->{'value'}}[$parameter->{'index'}]->{'content'};
		    $value =~ s/^\s*//;
		    $value =~ s/\s*$//;
		} else {
		    $value = ${$parameter->{'value'}}[$parameter->{'index'}];
		}
		$parameters{$parameter->{'name'}} = $value;
		$label .= ":".$parameter->{'ID'}.".".$parameter->{'index'};
		if ( exists($parameter->{'label'}) && ${$parameter->{'label'}}[0] eq "yes" ) {
		    $parameters{'label'} .= $joiner.$parameter->{'name'}."=".$value;
		    $joiner = ":";
		}
	    }
	}

	# If the parameter set with this label has not yet been added, add it now.
	unless ( exists($doneCases{$label}) ) {	   
	    %{$parameterHashes[++$#parameterHashes]} = %parameters;
	    $doneCases{$label} = 1;
	}

	# See if we find a parameter to increment.
	$done = 1;
	for(my $iParameter=$#parameterPointer;$iParameter>0;--$iParameter) {
	    if ( $parameterPointer[$iParameter]->{'index'} < $#{$parameterPointer[$iParameter]->{'value'}} ) {
		++$parameterPointer[$iParameter]->{'index'};
		for(my $jParameter=$iParameter+1;$jParameter<=$#parameterPointer;++$jParameter) {
		    $parameterPointer[$jParameter]->{'index'} = 0;
		}
		$done = 0;
		last;
	    }
	}
    }

    # Return the array of parameter hashes.
    return @parameterHashes;
}
