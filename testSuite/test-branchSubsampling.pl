#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test that the merger tree branch subsampling algorithm.
# Andrew Benson (14-July-2021)

# Make output directory.
system("mkdir -p outputs/");

# Run the subsampled model.
system("cd ..; Galacticus.exe testSuite/parameters/mergerTreeBranchSubsampled.xml");
unless ( $? == 0 ) {
    print "FAIL: merger tree branch subsampling model failed to run\n";
    exit;
}

# Run the not subsampled model.
system("cd ..; Galacticus.exe testSuite/parameters/mergerTreeBranchNotSubsampled.xml");
unless ( $? == 0 ) {
    print "FAIL: merger tree branch no subsampling model failed to run\n";
    exit;
}

# Read data and construct counts of subhalos.
my @models =
    (
     {
	 label    => "subsampled",
	 fileName => "outputs/mergerTreeBranchSubsampled.hdf5"
     },
     {
	 label    => "notSubsampled",
	 fileName => "outputs/mergerTreeBranchNotSubsampled.hdf5"
     }
    );
foreach my $model ( @models ) {
    my $modelFile = new PDL::IO::HDF5($model->{'fileName'});
    my $outputs   = $modelFile->group('Outputs' );
    my $output    = $outputs  ->group('Output1' );
    my $nodeData  = $output   ->group('nodeData');
    $model->{$_} = $nodeData->dataset($_)->get()
	foreach ( "basicMass", "nodeIsIsolated", "nodeSubsamplingWeight" );
}
my @massesLogarithmic = ( 11, 10, 9, 8 );
foreach my $model ( @models ) {
    $model->{'countSubhalos'     } = pdl zeros(scalar(@massesLogarithmic));
    $model->{'countSubhalosError'} = pdl zeros(scalar(@massesLogarithmic));
    for(my $i=0;$i<scalar(@massesLogarithmic);++$i) {
	my $selectionSubhalos                   = which(($model->{'nodeIsIsolated'} == 0) & ($model->{'basicMass'}->log10() >= $massesLogarithmic[$i]));
	my $selectionHalos                      = which(($model->{'nodeIsIsolated'} == 1)                                                       );
	$model->{'countSubhalos'     }->(($i)) .= $model->{'nodeSubsamplingWeight'}->($selectionSubhalos)        ->sumover()        /$model->{'nodeSubsamplingWeight'}->($selectionHalos)->sumover();
	$model->{'countSubhalosError'}->(($i)) .= $model->{'nodeSubsamplingWeight'}->($selectionSubhalos)->pow(2)->sumover()->sqrt()/$model->{'nodeSubsamplingWeight'}->($selectionHalos)->sumover();
    }
}
# Compare counts.
my $offsetScaled = abs(+$models[0]->{'countSubhalos'}-$models[1]->{'countSubhalos'})/sqrt(+$models[0]->{'countSubhalosError'}**2+$models[1]->{'countSubhalosError'}**2);
if ( any($offsetScaled > 3.0) ) {
    print "FAIL: merger tree branch subsampling changes subhalo mass function at > 3σ\n";
} else {
    print "SUCCESS: merger tree branch subsampling does not change subhalo mass function at > 3σ\n";
}
exit;
