#!/usr/bin/env perl
use lib "./perl";
use XML::Simple;
use Data::Dumper;
use Data::Compare;
use File::Copy;
use System::Redirect;
use MIME::Lite;
use IO::Compress::Simple;

# Script to run sets of Galacticus models, looping through sets of parameters and performing analysis on the results.
# Contains error reporting functionality.
# Andrew Benson (11-June-2010)

# Get command line arguments.
if ( $#ARGV != 0 ) {die("Usage: Run_Galacticus.pl <runFile>")};
$runFile = $ARGV[0];

# Read in the file of models to be run.
$xml = new XML::Simple;
$modelsToRun = $xml->XMLin($runFile, KeyAttr => "", ForceArray => 1);

# Read in any configuration options.
if ( -e "galacticusConfig.xml" ) {
    $xml = new XML::Simple;
    $config = $xml->XMLin("galacticusConfig.xml");
}

# Determine root directory for models.
$rootDirectory = "models";
if ( exists($modelsToRun->{'modelRootDirectory'}) ) {$rootDirectory  = ${$modelsToRun->{'modelRootDirectory'}}[0]};

# Determine the base set of parameters to use.
if ( exists($modelsToRun->{'baseParameters'}    ) ) {$baseParameters = ${$modelsToRun->{'baseParameters'}    }[0]};

# Record where we are running.
$pwd = `pwd`;
chomp($pwd);

# Set initial value of random seed.
$randomSeed = 219;

# Loop through all model sets.
$iModelSet = -1;
foreach $parameterSet ( @{$modelsToRun->{'parameters'}} ) {
    # Increment model set counter.
    ++$iModelSet;

    # Create an array of hashes giving the parameters for this parameter set.
    @parameterHashes = &Create_Parameter_Hashes($parameterSet);

    # Loop over all models and run them.
    $iModel = 0;
    foreach $parameterHash ( @parameterHashes ) {
	++$iModel;

	# Increment random seed.
	++$randomSeed;

	# Specify the output directory.
	$galacticusOutputDirectory = $rootDirectory."/galacticus_".$iModelSet.":".$iModel;
	system("mkdir -p ".$galacticusOutputDirectory);

	# Specify the output file.
	$galacticusOutputFile = $galacticusOutputDirectory."/galacticus.hdf5";

	# If the output file does not exist, then generate it.
	unless ( -e $galacticusOutputFile || -e $galacticusOutputFile.".bz2" ) {

	    # Read the default set of parameters.
	    unless ( $baseParameters eq "" ) {
		$xml = new XML::Simple;
		$data = $xml->XMLin($baseParameters,ForceArray => 1);
		@parameterArray = @{$data->{'parameter'}};
		for($i=0;$i<=$#parameterArray;++$i) {
		    $parameterHash{${$parameterArray[$i]->{'name'}}[0]} = ${$parameterArray[$i]->{'value'}}[0];
		}
	    }

	    # Set the output file name.
	    $parameterHash{'galacticusOutputFileName'} = $galacticusOutputFile;

	    # Set the random seed.
	    $parameterHash{'randomSeed'} = $randomSeed unless ( exists($parameterHash{'randomSeed'}) );

	    # Set a state restore file.
	    ($stateFile = $galacticusOutputFile) =~ s/\.hdf5//;
	    $parameterHash{'stateFileRoot'}      = $stateFile;

	    # Transfer parameters for this model from the array of model parameter hashes to the active hash.
	    foreach $parameter ( keys(%{$parameterHash}) ) {
		$parameterHash{$parameter} = ${$parameterHash}{$parameter};
	    }
	    
	    # Transfer values from the active hash to an array suitable for XML output.
	    delete($data->{'parameter'});
	    undef(@parameterArray);
	    foreach $name ( keys(%parameterHash) ) {
		$parameterArray[++$#parameterArray]->{'name'}  = $name;
		$parameterArray[  $#parameterArray]->{'value'} = $parameterHash{$name};
	    }
	    $data->{'parameter'} = \@parameterArray;

	    # Output the parameters as an XML file.
	    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
	    open(outHndl,">".$galacticusOutputDirectory."/newParameters.xml");
	    print outHndl $xmlOutput->XMLout($data);
	    close(outHndl);
	    undef($data);
	    undef(%parameterHash);

 	    # Run Galacticus.
	    &SystemRedirect::tofile("ulimit -t unlimited; ulimit -c unlimited; GFORTRAN_ERROR_DUMPCORE=YES; time Galacticus.exe "
				    .$galacticusOutputDirectory."/newParameters.xml",$galacticusOutputDirectory."/galacticus.log");
 	    if ( $? == 0 ) {
 		# Model finished successfully.
 		# Generate plots.
 		system("./scripts/analysis/Galacticus_Compute_Fit.pl ".$galacticusOutputFile." ".$galacticusOutputDirectory) unless ( ${$modelsToRun->{'doAnalysis'}}[0] eq "no" );
 	    } else {
 		# The run failed for some reason.
 		# Move the core file to the output directory.
 		opendir(gDir,".");
 		while ( $file = readdir(gDir) ) {
 		    if ( $file =~ m/core\.\d+/ ) {move($file,$galacticusOutputDirectory."/core")};
 		}
 		closedir(gDir);
 		# Report the model failure (by e-mail if we have an e-mail address to send a report to and if so requested).
		$message  = "FAILED: A Galacticus model failed to finish:\n\n";
		$message .= "  Host:\t".$ENV{"HOSTNAME"}."\n";
		$message .= "  User:\t".$ENV{"USER"}."\n\n";
		$message .= "Model output is in: ".$pwd."/".$galacticusOutputDirectory."\n\n";
 		if ( $config->{'contact'}->{'email'} =~ m/\@/ && ${$modelsToRun->{'emailReport'}}[0] eq "yes" ) {
 		    $message .= "Log file is attached.\n";
 		    $msg = MIME::Lite->new(
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
    }
}

exit;

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
		if ( exists($valueElement->{'parameter'}) ) {
		    $thisParameterArrays[++$#thisParameterArrays] = $valueElement->{'parameter'};
		}
	    }
	}
	if ( defined(@thisParameterArrays) ) {
	    # Add any detected parameter arrays to the list.
	    foreach my $thisParameterArray ( @thisParameterArrays ) {
		# Loop over each parameter in the array.
		foreach $parameter ( @{$thisParameterArray} ) {
		    # Store a pointer to the parameter.
		    $parameterPointer[++$#parameterPointer] = $parameter;
		    # Set an index counter (for its <value> elements) to zero.
		    $parameter->{'index'}  = 0;
		    # Set a unique ID for this parameter.
		    ++$parameterID;
		    $parameter->{'ID'}     = $parameterID;
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
		if ( exists($valueElement->{'parameter'}) ) {
		    $thisParameterArray = $valueElement->{'parameter'};
		}
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
	my $label = "";
	foreach $parameter ( @currentParameterPointer ) {
	    if ( exists($parameter->{'name'}) ) {
		if ( exists(${$parameter->{'value'}}[$parameter->{'index'}]->{'content'}) ) {
		    $value = ${$parameter->{'value'}}[$parameter->{'index'}]->{'content'};
		    $value =~ s/^\s*//;
		    $value =~ s/\s*$//;
		} else {
		    $value = ${$parameter->{'value'}}[$parameter->{'index'}];
		}
		$parameters{${$parameter->{'name'}}[0]} = $value;
		$label .= ":".$parameter->{'ID'}.".".$parameter->{'index'};
	    }
	}

	# If the parameter set with this label has not yet been added, add it now.
	unless ( exists($doneCases{$label}) ) {	   
	    %{$parameterHashes[++$#parameterHashes]} = %parameters;
	    $doneCases{$label} = 1;
	}

	# See if we find a parameter to increment.
	$done = 1;
	for($iParameter=$#parameterPointer;$iParameter>0;--$iParameter) {
	    if ( $parameterPointer[$iParameter]->{'index'} < $#{$parameterPointer[$iParameter]->{'value'}} ) {
		++$parameterPointer[$iParameter]->{'index'};
		for($jParameter=$iParameter+1;$jParameter<=$#parameterPointer;++$jParameter) {
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
