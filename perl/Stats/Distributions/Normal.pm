package Stats::Distributions::Normal;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw(PI);

# Functions related to the normal distribution.
# Andrew Benson (10-October-2016)

sub TruncatedMean {
    # Compute the mean of a normal distribution truncated below some value. Integrals evaluated using Mathematica.
    my $mean          = shift();
    my $rootVariance  = shift();
    my $truncate      = shift();
    my $meanTruncated = pdl zeroes(nelem($mean));
    (my $special, my $general) = which_both($mean == $truncate);
    # Handle the special case where the distribution is truncated below the mean.
    $meanTruncated->($special) .= $truncate->($special)+sqrt(2.0/PI)*$rootVariance->($special);
    # General case.
    my $argument                = ($truncate->($general)-$mean->($general))/sqrt(2.0)/$rootVariance->($general);
    my $erfc                    = erfc($argument);
    $meanTruncated->($general) .= 2.0*(exp(-$argument**2)*$rootVariance->($general)/sqrt(2.0*PI)+0.5*$mean->($general)*$erfc)/$erfc;
    return $meanTruncated;
}

sub TruncatedRootVariance {
    # Compute the root-variance of a normal distribution truncated below some value. Integrals evaluated using Mathematica.
    my $mean                  = shift();
    my $rootVariance          = shift();
    my $truncate              = shift();
    my $rootVarianceTruncated = pdl zeroes(nelem($mean));
    (my $special, my $general) = which_both($mean == $truncate);
    # Handle the special case where the distribution is truncated below the mean.
    $rootVarianceTruncated->($special) .= sqrt(PI-2.0)*$rootVariance->($special)*$mean->($special)/(sqrt(PI)*$truncate->($special)+sqrt(2.0)*$rootVariance->($special));
    # General case.
    my $argument                        = ($truncate->($general)-$mean->($general))/sqrt(2.0)/$rootVariance->($general);
    my $erfc                            = erfc($argument);
    $rootVarianceTruncated->($general) .= $mean->($general)*sqrt($erfc*(exp(-$argument**2)*$rootVariance->($general)*($truncate->($general)+$mean->($general))/sqrt(2.0*PI)+0.5*($rootVariance->($general)**2+$mean->($general)**2)*$erfc)/2.0/(exp(-$argument**2)*$rootVariance->($general)/sqrt(2.0*PI)+0.5*$mean->($general)*$erfc)**2-1.0);
    return $rootVarianceTruncated;
}

1;
