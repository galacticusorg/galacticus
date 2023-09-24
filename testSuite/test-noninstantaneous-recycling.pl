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

# Check calculations of noninstantaneous recycling.
# Andrew Benson (23-September-2023)

# Run the model.
&System::Redirect::tofile("mkdir -p outputs; cd ..; ./Galacticus.exe testSuite/parameters/noninstantaneous_recycling.xml","outputs/noninstantaneous_recycling.log");
unless ( $? == 0 ) {
    print "FAILED:  model run:\n";
    system("cat outputs/noninstantaneous_recycling.log");
} else {
    print "SUCCESS: model run\n";
}
# Read the model data and check for consistency.
my $model   = new PDL::IO::HDF5("outputs/noninstantaneous_recycling.hdf5");
my $outputs  = $model  ->group('Outputs' );
my $output   = $outputs->group('Output10');
my $nodeData = $output ->group('nodeData');
my $data;
$data->{$_} = $nodeData->dataset($_)->get()
    foreach ( 
	"spheroidAbundancesStellarMetals",
	"spheroidAbundancesStellarFe"    ,
	"diskAbundancesStellarMetals"    ,
	"diskAbundancesStellarFe"
    );
my $massMetals = $data->{'spheroidAbundancesStellarMetals'}+$data->{'diskAbundancesStellarMetals'};
my $massFe     = $data->{'spheroidAbundancesStellarFe'    }+$data->{'diskAbundancesStellarFe'    };
my $nonZero    = which($massMetals > 1.0);
my $ratio      = +$massFe    ->($nonZero)
                 /$massMetals->($nonZero); 
my $status = all($ratio < 0.16) ? "SUCCESS" : "FAILED";
print $status.": Fe/Z ratio\n";
exit 0;
