#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::IO::Misc;
use PDL::NiceSlice;
use XML::Simple;

# Test that the "half-mode slope" transfer function agrees with results originated from Daniel Gilman.
# Andrew Benson (11-March-2024)

# Read the data from Daniel Gilman's calculation.
my $transferFunctionGilman;
(my $wavenumberGilman, $transferFunctionGilman->{'-1.0'}, $transferFunctionGilman->{'-2.5'}, $transferFunctionGilman->{'-4.0'} ) = rcols("data/transferFunctionHalfModeSlopeGilman.txt", 0,1,2,3);

# Parse the base parameter file.
my $xml        = new XML::Simple();
my $parameters = $xml->XMLin("parameters/transferFunctionHalfModeSlope.xml");

# Iterate over slopes.
foreach my $slope ( sort(keys(%{$transferFunctionGilman})) ) {
    # Set Hubble constant - values in Daniel's file are in "h" units.
    my $h            = 0.674;
    # Set the half-mode mass.
    my $massHalfMode = 10.0**7.7/$h;
    # Write an updated parameter file.
    $parameters->{'task'            }->{'wavenumberMinimum'}->{'value'} = $wavenumberGilman->(( 0))->sclr()*$h;
    $parameters->{'task'            }->{'wavenumberMaximum'}->{'value'} = $wavenumberGilman->((-1))->sclr()*$h;
    $parameters->{'task'            }->{'pointsPerDecade'  }->{'value'} = sclr((nelem($wavenumberGilman)-1)/log10($wavenumberGilman->((-1))/$wavenumberGilman->((0))));
    $parameters->{'transferFunction'}->{'massHalfMode'     }->{'value'} = $massHalfMode;
    $parameters->{'transferFunction'}->{'slopeHalfMode'    }->{'value'} = $slope; 
    $parameters->{'outputFileName'  }                       ->{'value'} = "testSuite/outputs/transferFunctionHalfModeSlope.hdf5";
    open(my $outputFile,">","outputs/transferFunctionHalfModeSlope.xml");
    print $outputFile $xml->XMLout($parameters, RootName => "parameters");
    close($outputFile);
    # Run the model.
    system("mkdir -p outputs; cd ..; ./Galacticus.exe testSuite/outputs/transferFunctionHalfModeSlope.xml");
    unless ( $? == 0 ) {
	print "FAILED: [d㏒T/d㏒k=".$slope."] model did not run\n";
	next;
    }
    # Read the results.
    my $model            = new PDL::IO::HDF5("outputs/transferFunctionHalfModeSlope.hdf5");
    my $output           = $model ->group  ('Outputs/Output1' )       ;
    my $wavenumber       = $output->dataset('wavenumber'      )->get();
    my $transferFunction = $output->dataset('transferFunction')->get();
    $wavenumber    /= $h;
    if (nelem($transferFunction) != nelem($wavenumberGilman)) {
	print "FAILED: [d㏒T/d㏒k=".$slope."] count of wavenumbers does not match\n";
    } else {
	print "SUCCESS: [d㏒T/d㏒k=".$slope."] count wavenumbers matches\n";
    }
    if (any(abs($wavenumber-$wavenumberGilman) > 1.0e-3*$wavenumberGilman)) {
	print "FAILED: [d㏒T/d㏒k=".$slope."] wavenumbers do not match\n";
    } else {
	print "SUCCESS: [d㏒T/d㏒k=".$slope."] wavenumbers match\n";
    }
    if (any(abs($transferFunction-$transferFunctionGilman->{$slope}) > 1.0e-3*$transferFunctionGilman->{$slope})) {
	print "FAILED: [d㏒T/d㏒k=".$slope."] transfer functions do not match\n";
    } else {
	print "SUCCESS: [d㏒T/d㏒k=".$slope."] transfer functions match\n";
    }
}

exit;
