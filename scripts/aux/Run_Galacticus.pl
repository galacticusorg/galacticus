#!/usr/bin/env perl
use lib "./perl";
use XML::Simple;
use Data::Dumper;
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
if ( exists($modelsToRun->{'modelRootDirectory'}) ) {$rootDirectory = ${$modelsToRun->{'modelRootDirectory'}}[0]};

# Determine the base set of parameters to use.
if ( exists($modelsToRun->{'baseParameters'}) ) {$baseParameters = ${$modelsToRun->{'baseParameters'}}[0]};

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

    # Count up the number of models to run and the periodicity in each parameter.
    undef(%modelsPeriod);
    $modelsCount = 1;
    foreach $parameter ( @{$parameterSet->{'parameter'}} ) {
	$modelsPeriod{${$parameter->{'name'}}[0]} = $modelsCount;
	$modelsCount *= $#{$parameter->{'value'}}+1;
    }
    
    # Generate an array of parameter hashes which specifies the parameter values for each model.
    for($iModel=0;$iModel<$modelsCount;++$iModel) {
	undef(%{$models[$iModel]});
	foreach $parameter ( @{$parameterSet->{'parameter'}} ) {
	    $index = int($iModel/$modelsPeriod{${$parameter->{'name'}}[0]}) % ($#{$parameter->{'value'}}+1);
	    ${$models[$iModel]}{${$parameter->{'name'}}[0]} = ${$parameter->{'value'}}[$index];
	}
    }

    # Loop over all models and run them.
    for($iModel=0;$iModel<$modelsCount;++$iModel) {

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
	    $parameterHash{'galacticusOutputFile'} = $galacticusOutputFile;

	    # Set the random seed.
	    $parameterHash{'randomSeed'} = $randomSeed unless ( exists($parameterHash{'randomSeed'}) );

	    # Set a state restore file.
	    ($stateFile = $galacticusOutputFile) =~ s/\.hdf5//;
	    $parameterHash{'stateFileRoot'}      = $stateFile;

	    # Transfer parameters for this model from the array of model parameter hashes to the active hash.
	    foreach $parameter ( keys(%{$models[$iModel]}) ) {
		$parameterHash{$parameter} = ${$models[$iModel]}{$parameter};
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
 		system("./scripts/analysis/Galacticus_Compute_Fit.pl ".$galacticusOutputFile." ".$galacticusOutputDirectory);
 	    } else {
 		# The run failed for some reason.
 		# Move the core file to the output directory.
 		opendir(gDir,".");
 		while ( $file = readdir(gDir) ) {
 		    if ( $file =~ m/core\.\d+/ ) {move($file,$galacticusOutputDirectory."/core")};
 		}
 		closedir(gDir);
 		# If we have an e-mail address to send a report to, then do so.
 		if ( $config->{'contact'}->{'email'} =~ m/\@/ ) {
 		    $message  = "A Galacticus model failed to finish:\n\n";
 		    $message .= "  Host:\t".$ENV{"HOSTNAME"}."\n";
 		    $message .= "  User:\t".$ENV{"USER"}."\n\n";
 		    $message .= "Model output is in: ".$pwd."/".$galacticusOutputDirectory."\n\n";
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
 		}
 	    }

 	    # Compress all files in the output directory.
	    &Simple::Compress_Directory($galacticusOutputDirectory);

	}
    }
}


exit;
