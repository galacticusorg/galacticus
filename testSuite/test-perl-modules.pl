#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
 $ENV{"GALACTICUS_ROOT_V092"} = "/";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
require Stats::Means;

# Test Perl modules.
# Andrew Benson (7-July-2013)

# Check individual compilation.
system("find . -type f -name \"*.pl\" | xargs -n 1 perl -c 1>perl.tmp 2>&1; sed -r s/failed/FAILED/g perl.tmp; rm perl.tmp");
system("find perl -type f -name \"*.pm\" | xargs -n 1 perl -c 1>perl.tmp 2>&1; sed -r s/failed/FAILED/g perl.tmp; rm perl.tmp");

# Statistics
#  Create bins.
my $bins = pdl sequence(3);

## Means
my $x = pdl ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2 );
my $y = pdl ( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 );
my $w = pdl ( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1,  1,  1,  1,   1,  1, 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1 );
my $eps      = pdl 1.0e-3;
my $mExpect  = pdl ( 4.5       , 14.5       , 24.5        );
my $meExpect = pdl ( 0.90829511,  0.90829511,  0.90829511 );
my $dExpect  = pdl ( 3.0276503 ,  3.0276503 ,  3.0276503  );
my $deExpect = pdl ( 1.6371658 ,  1.6371658 ,  1.6371658  );
(my $m, my $me, my $d, my $de) = &Means::BinnedMean($bins,$x,$y,$w);
print "FAILED: Stats::Means::BinnedMean fails to compute correct means\n"
    unless ( all(abs($m - $mExpect) < $eps*$mExpect ) );
print "FAILED: Stats::Means::BinnedMean fails to compute correct errors on means\n"
    unless ( all(abs($me-$meExpect) < $eps*$meExpect) );
print "FAILED: Stats::Means::BinnedMean fails to compute correct standard deviation\n"
    unless ( all(abs($d - $dExpect) < $eps*$dExpect ) );
print "FAILED: Stats::Means::BinnedMean fails to compute correct errors on standard deviation\n"
    unless ( all(abs($de-$deExpect) < $eps*$deExpect) );

exit;
