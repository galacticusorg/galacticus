#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;

# Run Galacticus models with and without parallelism in tree building.
# Andrew Benson (05-July-2023)

# Ensure output directory exists.
system("mkdir -p outputs");

# Run models.
system("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/parallelTreeBuildSerial.xml");
die("FAILED: failed to run serial model")
    unless ( $? == 0 );
system("export OMP_NUM_THREADS=4; cd ..; ./Galacticus.exe testSuite/parameters/parallelTreeBuildParallel.xml");
die("FAILED: failed to run parallel model")
    unless ( $? == 0 );

# Read model data.
my $data;
foreach my $modelName ( "parallelTreeBuildSerial", "parallelTreeBuildParallel" ) {
    my $model    = new PDL::IO::HDF5("outputs/".$modelName.".hdf5");
    my $analysis = $model->group("conditionalMassFunction");
    $data->{$modelName}->{$_} = $analysis->dataset($_)->get()
	foreach ( "conditionalMassFunction", "conditionalMassFunctionError" );
}

# Compare results.
my $selection       = which($data->{'parallelTreeBuildSerial'}->{'conditionalMassFunction'}->flat() > 0.0);
my $errorNormalized = abs($data->{'parallelTreeBuildParallel'}->{'conditionalMassFunction'}->flat()->($selection)-$data->{'parallelTreeBuildSerial'}->{'conditionalMassFunction'}->flat()->($selection))/$data->{'parallelTreeBuildSerial'}->{'conditionalMassFunctionError'}->flat()->($selection);
my $status          = any($errorNormalized > 4.0) ? "FAIL" : "SUCCESS";
print $status.": parallel tree build\n";
if ( $status eq "FAIL" ) {
    print "Conditional mass function:\n";
    print "\t  Serial build: ".$data->{'parallelTreeBuildSerial'  }->{'conditionalMassFunction'     }->flat()->($selection)."\n";
    print "\tParallel build: ".$data->{'parallelTreeBuildParallel'}->{'conditionalMassFunction'     }->flat()->($selection)."\n";
    print "\t   Uncertainty: ".$data->{'parallelTreeBuildSerial'  }->{'conditionalMassFunctionError'}->flat()->($selection)."\n";
}

exit;
