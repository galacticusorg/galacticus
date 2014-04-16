#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Find the maximum likelihood estimate of the covariance matrix for the Caputi et al. (2011) UKIDSS UDS stellar mass functions.
# Andrew Benson (05-July-2012)

# Get the argument specifying which redshift bin to use.
die("Usage: covarianceMatrix.pl <redshiftBin>")
    unless ( scalar(@ARGV) == 1 );
my $redshiftBin = $ARGV[0];
die("covarianceMatrix.pl: redshiftBin must be 0, 1, or 2")
    if ( $redshiftBin < 0 || $redshiftBin > 2 );

# Simply run the generic script with our config file as argument.
system($galacticusPath."constraints/dataAnalysis/scripts/covarianceMatrix.pl ".$galacticusPath."constraints/dataAnalysis/stellarMassFunctions_UKIDSS_UDS_z3_5/covarianceMatrixControl.xml ".$redshiftBin);

exit;
