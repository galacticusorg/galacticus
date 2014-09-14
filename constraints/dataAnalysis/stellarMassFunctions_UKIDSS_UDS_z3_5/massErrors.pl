#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::Fit::Polynomial;
use Astro::Cosmology;

# Estimate random errors on stellar masses in the UKISS UDS survey of Caputi et al. (2013).
# Andrew Benson (07-May-2014)

# Central redshifts of the three redshift bins considered.
my $redshift = pdl [ 3.250, 3.875, 4.625 ];

# Fractional error in dz/(1+z) reported by Caputi et al. (2013).
my $sigmaLogZ = pdl 0.05;

# Speed of light in km/s.
my $speedLight = pdl 2.99792458e5;

# Construct cosmological model.
my $cosmology = Astro::Cosmology->new(omega_matter => 0.283812448723631, omega_lambda => 0.716187551276369, h0 => 69.5723630486537);

# Find comoving distances to each redshift.
my $distanceComoving = $cosmology->comoving_distance($redshift);

# Find Hubble parameter at each redshift.
my $hubbleParameter = $cosmology->h0()*sqrt($cosmology->omega_matter()*(1.0+$redshift)**3+$cosmology->omega_lambda());

# Compute the error in log10(M*) from photometric redshift errors.
my $sigmaLog10MassRedshift = (2.0+2.0*(1.0+$redshift)**2*$speedLight/$hubbleParameter/$distanceComoving)*$sigmaLogZ/log(10.0);

# Additional error arising from SED fitting (judged from Figure 8 of Caputi et al. 2013).
my $sigmaLog10MassSED      = pdl 0.2/log(10.0);

# Compute final error.
my $sigmaLog10Mass         =sqrt($sigmaLog10MassRedshift**2+$sigmaLog10MassSED**2);

# Display the results.
print "Redshift\tDispersion in log10 stellar mass\n";
for(my $i=0;$i<nelem($redshift);++$i) {
    print $redshift->(($i))."\t".$sigmaLog10Mass->(($i))."\n";
}

exit;
