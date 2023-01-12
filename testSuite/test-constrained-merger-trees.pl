#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use File::Slurp qw(slurp);
use System::Redirect;

# Check construction of constrained merger trees.
# Andrew Benson (01-December-2022)

# Run the model.
&System::Redirect::tofile("cd ..; ./Galacticus.exe testSuite/parameters/constrainedMergerTrees.xml","outputs/constrainedMergerTrees.log");
unless ( $? == 0 ) {
    print "FAILED:  model run:\n";
    system("cat outputs/constrainedMergerTrees.log");
} else {
    print "SUCCESS: model run\n";
}

# Define the constraint.
my $redshiftConstraint = 8.0;
my $massConstraint     = 1.0e12;

# Read the model data and check constraint is met.
my $model = new PDL::IO::HDF5("outputs/constrainedMergerTrees.hdf5");
my $trees = $model->group('mergerTreeStructures');
my $massesClosest = pdl [];
foreach my $treeName ( $trees->groups() ) {
    my $tree                  = $trees->group  ($treeName          )       ;
    my $isOnConstrainedBranch = $tree ->dataset('nodeIsConstrained')->get();
    my $redshift              = $tree ->dataset('redshift'         )->get();
    my $mass                  = $tree ->dataset('massBasic'        )->get();
    my $constrainedBranch     = which(($isOnConstrainedBranch == 1) & ($redshift >= $redshiftConstraint));
    my $redshiftDelta         = abs($redshift->($constrainedBranch)-$redshiftConstraint);
    my $redshiftClosest       = minimum_ind($redshiftDelta);
    $massesClosest            = $massesClosest->append($mass->($constrainedBranch)->(($redshiftClosest)));
}   
# Assert that constrained branch has the expected mass.
my $tolerance      = 1.0e-3;
my $status         = all(abs($massesClosest-$massConstraint) < $massConstraint*$tolerance) ? "SUCCESS" : "FAILED";
print $status.": constrained merger trees\n";

exit 0;
