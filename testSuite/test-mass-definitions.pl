#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);

# Test mass definition conversions.
# Andrew Benson (08-February-2023)

# Make output directory.
system("mkdir -p outputs/");

# Iterate over time options.
foreach my $optionTime ( "current", "infall" ) {
    # Run the model.
    system("cd ..; ./Galacticus.exe testSuite/parameters/massDefinitionsTime".ucfirst($optionTime).".xml");
    if ( $? == 0 ) {
	# Read all data.
	my $model                             = new PDL::IO::HDF5("outputs/massDefinitionsTime".ucfirst($optionTime).".hdf5");
	my $cosmology                         = $model    ->group  ('Parameters'                  )->group('cosmologyParameters');
	my $outputs                           = $model    ->group  ('Outputs'                     )                              ;
	my $output                            = $outputs  ->group  ('Output1'                     )                              ;
	my $nodes                             = $output   ->group  ('nodeData'                    )                              ;
	(my $OmegaMatter, my $HubbleConstant) = $cosmology->attrGet("OmegaMatter","HubbleConstant")                              ;
	my $halos;
	$halos->{$_} = $nodes->dataset($_)->get()
	    foreach ( "redshift", "redshiftLastIsolated", "basicMass", "massHaloEnclosedCurrent", "darkMatterOnlyRadiusVirial" );
	# Find mean density of the universe at the present time.
	my $gravitationalConstant = pdl 4.301e-9;
	my $densityMean           = $OmegaMatter*3.0*$HubbleConstant**2/8.0/PI/$gravitationalConstant;
	if ( $optionTime eq "current" ) {
	    # Density definition is at the current time - so just use the mean density.
	    $halos->{'densityMean'} = $densityMean;
	} elsif ( $optionTime eq "infall" ) {
	    # Density definition is at the infall time - so scale the mean density back to the infall redshift of each halo.
	    $halos->{'densityMean'} = $densityMean*((1.0+$halos->{'redshiftLastIsolated'})/(1.0+$halos->{'redshift'}))**3;
	} else {
	    die('unrecognized time option');
	}
	# Compute the target density, which is 200 times the mean density.
	my $densityTarget = 200.0*$halos->{'densityMean'};
	# Find the density of halos under the internal mass definition.
	$halos->{'density'} = $halos->{'basicMass'}/(4.0*PI*$halos->{'darkMatterOnlyRadiusVirial'}**3/3.0);
	# For isothermal halos, the enclosed mass is inversely proportional to the square root of the enclosed density. So, find
	# the mass under the alternative mass definition using this scaling.
	$halos->{'massConverted'} = $halos->{'basicMass'}*sqrt($halos->{'density'}/$densityTarget);
	# Compute the difference between what the mass should be, and what Galacticus computed.
	my $error = abs($halos->{'massHaloEnclosedCurrent'}-$halos->{'massConverted'})/$halos->{'massConverted'};
	if ( all($error < 1.0e-4) ) {
	    print "SUCCESS: mass definitions for '".$optionTime."' time\n";
	} else {
	    print "FAIL: mass definitions for '"   .$optionTime."' time\n";
	}
    } else {
	print "FAIL: failed to run mass definitions model\n";
    }
}

exit;
