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
use List::Util;
require Galacticus::Constraints::Parameters;
require Galacticus::Constraints::DiscrepancySystematics;
require Galacticus::Launch::PBS;

# Run calculations to determine the model discrepancy arising from the use of a cosmology-independent set of merger trees.
# Andrew Benson (09-January-2014)

# Get arguments.
die("Usage: fixedTreeCosmology.pl <configFile> [options]") unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments = 
    (
     make                              => "yes",
     sampleCount                       =>  100 ,
     sampleCosmology                   => "yes",
     sampleHaloMassFunctionSystematics => "yes"
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
my $config = &Parameters::Parse_Config($configFile);
my @simulationParameters;
if ( UNIVERSAL::isa($config->{'parameters'}->{'parameter'},"ARRAY") ) {
    @simulationParameters = @{$config->{'likelihood'}->{'parameters'}->{'parameter'}};
} else {
    push(@simulationParameters,$config->{'likelihood'}->{'parameters'}->{'parameter'});
}

# Validate the config file.
die("fixedTreeCosmology.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'likelihood'}->{'workDirectory' }) );
die("fixedTreeCosmology.pl: compilation must be specified in config file"   ) unless ( exists($config->{'likelihood'}->{'compilation'   }) );
die("fixedTreeCosmology.pl: baseParameters must be specified in config file") unless ( exists($config->{'likelihood'}->{'baseParameters'}) );

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
    die("fixedTreeCosmology.pl: failed to build Galacticus.exe")
	unless ( $? == 0 );
}

# Get a hash of the parameter values.
(my $constraintsRef, my $parameters) = &Parameters::Compilation($config->{'likelihood'}->{'compilation'},$config->{'likelihood'}->{'baseParameters'});
my @constraints = @{$constraintsRef};

# Switch off thread locking.
$parameters->{'parameter'}->{'treeEvolveThreadLock'}->{'value'} = "false";

# Specify the directory and file name for the merger tree file.
my $treeDirectory = $workDirectory."/data/";
system("mkdir -p ".$treeDirectory);
my $treeFile         = $treeDirectory."fixedTrees.hdf5";
my $treeFileInternal = $treeDirectory."fixedTreesInternal.hdf5";

# Specify models to run.
my @parameters =
    (
     # Simply output the merger trees.
     {
	 name  => "mergerTreeConstructMethod",
	 value => "build"
     },
     {
	 name  => "mergerTreesWrite",
	 value => "true"
     },
     {
	 name  => "mergerTreeExportOutputFormat",
	 value => "galacticus"
     },
     {
	 name  => "mergerTreeExportFileName",
	 value => $treeFile
     },
     {
	 name  => "cosmology0",
	 value => 0.0
     },
     {
	 name  => "cosmology1",
	 value => 0.0
     },
     {
	 name  => "cosmology2",
	 value => 0.0
     },
     {
	 name  => "cosmology3",
	 value => 0.0
     },
     {
	 name  => "cosmology4",
	 value => 0.0
     },
     {
	 name  => "cosmology5",
	 value => 0.0
     },
     {
	 name  => "haloMassFunctionSystematic1",
	 value => 0.0
     },
     {
	 name  => "haloMassFunctionSystematic2",
	 value => 0.0
     }
    );
my @parametersInternal =
    (
     # Simply output the merger trees.
     {
	 name  => "mergerTreeConstructMethod",
	 value => "build"
     },
     {
	 name  => "mergerTreesWrite",
	 value => "true"
     },
     {
	 name  => "mergerTreeExportOutputFormat",
	 value => "galacticus"
     },
     {
	 name  => "mergerTreeExportFileName",
	 value => $treeFileInternal
     },
     {
	 name  => "cosmology0",
	 value => 0.0
     },
     {
	 name  => "cosmology1",
	 value => 0.0
     },
     {
	 name  => "cosmology2",
	 value => 0.0
     },
     {
	 name  => "cosmology3",
	 value => 0.0
     },
     {
	 name  => "cosmology4",
	 value => 0.0
     },
     {
	 name  => "cosmology5",
	 value => 0.0
     },
     {
	 name  => "haloMassFunctionSystematic1",
	 value => 0.0
     },
     {
	 name  => "haloMassFunctionSystematic2",
	 value => 0.0
     }
    );
# Adjust the number of trees to run if specified.
push(
    @parametersInternal,
    {
	name  => "mergerTreeBuildTreesPerDecade",
	value => $arguments{'treesPerDecade'}
    }
    )
    if ( exists($arguments{'treesPerDecade'}) );
my @generator = 
    (
     {
	 label      => "fixedTreeInternalGeneration",
	 parameters => \@parametersInternal
     },
     {
	 label      => "fixedTreeGeneration",
	 parameters => \@parameters
     }
    );
my @sample;
for(my $i=0;$i<$arguments{'sampleCount'};++$i) {
    # Specify base parameters.
    my @parameters = 
	(
	 {
	     name  => "mergerTreeConstructMethod",
	     value => "read"
	 },
	 {
	     name  => "mergerTreeReadFileName",
	     value => $treeFileInternal
	 },
	 {
	     name  => "mergerTreeReadPresetMergerTimes",
	     value => "false"
	 },
	 {
	     name  => "mergerTreeReadPresetMergerNodes",
	     value => "false"
	 },
	 {
	     name  => "mergerTreeReadPresetSubhaloMasses",
	     value => "false"
	 },
	 {
	     name  => "mergerTreeReadPresetPositions",
	     value => "false"
	 },
	 {
	     name  => "mergerTreeReadPresetScaleRadii",
	     value => "false"
	 },
	 {
	     name  => "mergerTreeReadPresetSpins",
	     value => "false"
	 },
	 {
	     name  => "mergerTreeReadPresetOrbits",
	     value => "false"
	 },
	 {
	     name  => "mergerTreeImportGalacticusMismatchIsFatal",
	     value => "false"
	 },
	 {
	     name  => "mergerTreeImportGalacticusReweightTrees",
	     value => "true"
	 }
	);
    # Generate random deviates.
    my $cosmology;
    if ( $arguments{'sampleCosmology'} eq "yes" ) {
	$cosmology = pdl grandom(6);
    } else {
	$cosmology = pdl zeroes(6);
    }
    for(my $j=0;$j<6;++$j) {
	push(
	    @parameters,
	    {
		name => "cosmology".$j,
		value => $cosmology->(($j))->sclr()
	    }
	    );
    }
    my $haloMassFunction;
    if ( $arguments{'sampleHaloMassFunctionSystematics'} eq "yes" ) {
	$haloMassFunction = pdl grandom(2);
    } else {
	$haloMassFunction = pdl zeroes(2);
    }
    for(my $j=0;$j<2;++$j) {
	my $k = $j+1;
	push(
	    @parameters,
	    {
		name => "haloMassFunctionSystematic".$k,
		value => $haloMassFunction->(($j))->sclr()
	    }
	    );
    }
    # Adjust the number of trees to run if specified.
    push(
	@parameters,
	{
	    name  => "mergerTreeBuildTreesPerDecade",
	    value => $arguments{'treesPerDecade'}
	}
	)
	if ( exists($arguments{'treesPerDecade'}) );
    # Push the model onto the stack.
    push(
	@sample,
	{
	    label      => "fixedTreeSample".$i,
	    parameters => \@parameters
	}
	);
}
for(my $i=0;$i<$arguments{'sampleCount'};++$i) {
    # Specify base parameters.
    my @parameters = 
	(
	 {
	     name  => "mergerTreeConstructMethod",
	     value => "build"
	 }
	);
    # Generate random deviates.
    my $cosmology;
    if ( $arguments{'sampleCosmology'} eq "yes" ) {
	$cosmology = pdl grandom(6);
    } else {
	$cosmology = pdl zeroes(6);
    }
    for(my $j=0;$j<6;++$j) {
	push(
	    @parameters,
	    {
		name => "cosmology".$j,
		value => $cosmology->(($j))->sclr()
	    }
	    );
    }
    my $haloMassFunction;
    if ( $arguments{'sampleHaloMassFunctionSystematics'} eq "yes" ) {
	$haloMassFunction = pdl grandom(2);
    } else {
	$haloMassFunction = pdl zeroes(2);
    }
    for(my $j=0;$j<2;++$j) {
	my $k = $j+1;
	push(
	    @parameters,
	    {
		name => "haloMassFunctionSystematic".$k,
		value => $haloMassFunction->(($j))->sclr()
	    }
	    );
    }
    # Adjust the number of trees to run if specified.
    push(
	@parameters,
	{
	    name  => "mergerTreeBuildTreesPerDecade",
	    value => $arguments{'treesPerDecade'}
	}
	)
	if ( exists($arguments{'treesPerDecade'}) );
    # Push the model onto the stack.
    push(
	@sample,
	{
	    label      => "newTreeSample".$i,
	    parameters => \@parameters
	}
	);
}
&Generate_Models(@generator);
&Generate_Models(@sample   );

# Iterate over constraints.
foreach my $constraint ( @constraints ) {
    # Parse the definition file.
    my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
    # Initialize data.
    my $treeMass;
    my $treeMassFunction;
    my $treeMassFunctionMean;
    my $treeMassFunctionCovariance;
    # Iterate over tree type.
    foreach my $type ( "fixed", "new" ) {
	# Iterate over tree models.
	for(my $i=0;$i<$arguments{'sampleCount'};++$i) {
	    # Locate the tree model results.
	    my $treeResultFileName   = $workDirectory."/modelDiscrepancy/fixedTreeCosmology/".$type."TreeSample".$i."/".$constraintDefinition->{'label'}.".hdf5";
	    # Read the results.
	    my $treeResult           = new PDL::IO::HDF5($treeResultFileName);
	    # Extract the results.
	    $treeMass                = $treeResult->dataset('x')->get()
	       unless ( defined($treeMass) );
	    my $thisTreeMassFunction = $treeResult->dataset('y')->get();
	    push(
		@{$treeMassFunction->{$type}},
		$thisTreeMassFunction
		);
	}
	# Find the mean mass function,
	foreach ( @{$treeMassFunction->{$type}} ) {
	    $treeMassFunctionMean->{$type} += $_;
	}
	$treeMassFunctionMean->{$type} /= $arguments{'sampleCount'};
	# Find the covariance of the mass function.
	foreach ( @{$treeMassFunction->{$type}} ) {
	    my $delta = $_-$treeMassFunctionMean->{$type};
	    $treeMassFunctionCovariance->{$type} += outer($delta,$delta);
	}
	$treeMassFunctionCovariance->{$type} /= $arguments{'sampleCount'}-1;
    }
    # Apply any systematics models.
    my %systematicResults;
    foreach my $argument ( keys(%arguments) ) {
	if ( $argument =~ m/^systematic(.*)/ ) {
	    my $model = $1;
	    if ( exists($DiscrepancySystematics::models{$model}) ) {
		%{$systematicResults{$model}} =
		    &{$DiscrepancySystematics::models{$model}}(
		    \%arguments                           ,
		    $treeMass                             ,
		    $treeMassFunctionMean      ->{'fixed'},
		    $treeMassFunctionCovariance->{'fixed'},
		    $treeMassFunctionMean      ->{'new'  },
		    $treeMassFunctionCovariance->{'new'  }
		);
	    }
	}
    }
    # Find the covariance due to tree variation.
    my $treeVariationCovariance = $treeMassFunctionCovariance->{'new'}-$treeMassFunctionCovariance->{'fixed'};

    # Ensure diagonal elements are non-negative.
    my $negativeDiagonal = which($treeVariationCovariance->diagonal(0,1) < 0.0);
    $treeVariationCovariance->diagonal(0,1)->($negativeDiagonal) .= 0.0;
    # Convert to a multiplicative covariance.
    my $multiplicativeCovariance = $treeVariationCovariance*outer(1.0/$treeMassFunctionMean->{'new'},1.0/$treeMassFunctionMean->{'new'});
    # Output the model discrepancy to file.
    my $outputFile = new PDL::IO::HDF5(">".$workDirectory."/modelDiscrepancy/fixedTreeCosmology/discrepancy".ucfirst($constraintDefinition->{'label'}).".hdf5");
    $outputFile->dataset('multiplicativeCovariance')->set($multiplicativeCovariance);
    $outputFile->attrSet(
    	description => "Model discrepancy for ".$constraintDefinition->{'name'}." due to use of merger trees that are cosmology-independent."
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
}

exit;

sub Generate_Models {
    my @models = @_;
    # Initialize a stack for PBS models.
    my @pbsStack;
    # Iterate over models.
    foreach my $model ( @models ) {
	# Specify the output name.
	my $modelDirectory = $workDirectory."/modelDiscrepancy/fixedTreeCosmology/".$model->{'label'};
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
	foreach ( @simulationParameters ) {
	    if ( exists($_->{'define'}) ) {
		die ("fixedTreeCosmology.pl: cannot specify a prior for a defined parameter")
		    if ( exists($_->{'prior'}) );
		$newParameters->{'parameter'}->{$_->{'name'}}->{'value'} = $_->{'define'}
	    }
	}
	# Set the values of any parameters that are defined in terms of other parameters.
	my $failCount = 1;
	while ( $failCount > 0 ) {
	    $failCount = 0;
	    foreach ( keys(%{$newParameters->{'parameter'}}) ) {
		# Attempt to replace named parameters in the definition with their values.
		while ( $newParameters->{'parameter'}->{$_}->{'value'} =~ m/\%([a-zA-Z0-9_]+)/ ) {
		    my $parameterName = $1;
		    if ( exists($newParameters->{'parameter'}->{$parameterName}) && $newParameters->{'parameter'}->{$parameterName}->{'value'} !~ m/\%([a-zA-Z0-9_]+)/ ) {
			$newParameters->{'parameter'}->{$_}->{'value'} =~ s/\%$parameterName/$newParameters->{'parameter'}->{$parameterName}->{'value'}/g;			
		    } else {
			++$failCount;
			last;
		    }
		    $newParameters->{'parameter'}->{$_}->{'value'} = eval($newParameters->{'parameter'}->{$_}->{'value'})
			unless ( $newParameters->{'parameter'}->{$_}->{'value'} =~ m/\%([a-zA-Z0-9_]+)/ );
		}
	    }
	}

	# Run the model.
	unless ( -e $galacticusFileName ) {
	    # Generate the parameter file.
	    my $parameterFileName = $modelDirectory."/parameters.xml";
	    &Parameters::Output($newParameters,$parameterFileName);
	    # Create a batch script for PBS.
	    my $command = "time mpirun --bynode -np 1 Galacticus.exe ".$parameterFileName."\n";
	    foreach my $constraint ( @constraints ) {
		# Parse the definition file.
		my $constraintDefinition = $xml->XMLin($constraint->{'definition'});
		# Insert code to run the analysis code.
		my $analysisCode = $constraintDefinition->{'analysis'};
		$command .= $analysisCode." ".$galacticusFileName." --resultFile ".$modelDirectory."/".$constraintDefinition->{'label'}.".hdf5";
	    }
	    my %job =
		(
		 launchFile => $modelDirectory."/launch.pbs",
		 label      => "fixedTreeCosmology".$model->{'label'},
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
}
