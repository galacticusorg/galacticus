#!/usr/bin/env perl
use strict;
use warnings;
use PDL;

# Compute the mean and variance of the subhalo halo occupation
# distribution slope from the simulation results reported by Kravtsov et
# al. (2004, ApJ, 609, 35; http://adsabs.harvard.edu/abs/2004ApJ...609...35K).
# Andrew Benson (23-July-2012)

# Set arrays of the reported slopes, alpha, and the errors on these values.
my $alpha = pdl ( 0.99, 0.92, 0.96, 1.04, 0.61 );
my $error = pdl ( 0.01, 0.03, 0.08, 0.08, 0.21 );

# Compute a weight equal to the inverse variance on each measurement.
my $weight = 1.0/$error**2;

# Compute summary statistics of the slope.
(my $mean, my $rms, my $median, my $min, my $max, my $adev, my $meanrms) = statsover($alpha,$weight);
my $variance = $rms**2;

# Report.
print "Statistics o of the slope of the subhalo halo occupation distribut:on.\n";
print "                Mean: ".$mean."\n";
print "  Standard deviation: ".$rms."\n";
print "            Variance: ".$variance."\n";

exit;
