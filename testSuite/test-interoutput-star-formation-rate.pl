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

# Check calculations of mean interoutput star formation rates.
# Andrew Benson (23-August-2021)

# Run the model.
system("cd ..; ./Galacticus.exe testSuite/parameters/interoutputStarFormationRate.xml");
unless ( $? == 0 ) {
    print "FAILED: model run\n";
} else {
    print "SUCCESS: model run\n";
}

# Read the model data and check for consistency.
my $model              = new PDL::IO::HDF5("outputs/interoutputStarFormationRate.hdf5");
my $parameters         = $model             ->group ('Parameters'       );
my $outputs            = $model             ->group ('Outputs'          );
my $stellarPopulation  = $parameters        ->group ('stellarPopulation');
(my $recycledFraction) = $stellarPopulation->attrGet('recycledFraction' ); 
my $timePrevious       = pdl 0.0;
my $massPrevious       = pdl zeros(6);
my $status             = 1;
for(my $i=1;$i<4;++$i) {
    my  $output                        = $outputs  ->group  ('Output'                              .$i)          ;
    my  $nodeData                      = $output   ->group  ('nodeData'                               )          ;
    (my $time                        ) = $output   ->attrGet('outputTime'                             )          ;
    my  $rateStarFormationInterOutput  = $nodeData ->dataset('diskStarFormationRateInterOutputMean'   )->get   ();
    my  $massStellar                   = $nodeData ->dataset('diskMassStellar'                        )->get   ();
    my  $treeIndex                     = $nodeData ->dataset('mergerTreeIndex'                        )->get   ();
    my $order                          = $treeIndex                                                    ->qsorti();
    # Compute increase in mass diretly and from star formation rate.
    my $massIncrease                   = $massStellar                 ->($order)                                      -$massPrevious;
    my $timeIncrease                   = $time                                                                        -$timePrevious;
    my $massIncreaseFromRate           = $rateStarFormationInterOutput->($order)*$timeIncrease*(1.0-$recycledFraction)              ;
    $status = 0
	if ( any((abs($massIncreaseFromRate-$massIncrease) > 1.0e-2*$massIncrease) & ($massIncrease > 1.0)) );
    $timePrevious .= $time                 ;
    $massPrevious .= $massStellar->($order);
}
my $success = $status ? "success" : "FAIL";
print $success.": mean interoutput star formation rate\n";

exit 0;
