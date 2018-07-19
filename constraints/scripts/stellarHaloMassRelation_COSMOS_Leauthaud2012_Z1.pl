#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Constraints::StellarHaloMassRelation;

# Compute likelihood (and make a plot) for a Galacticus model given the stellar mass-halo mass relation of Leauthaud et
# al. (2012).

# Get name of input and output files.
die("stellarHaloMassRelation_Leauthaud2012_Z1.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named options.
my $iArg = -1;
my %options =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Evaluate the model.
&Galacticus::Constraints::StellarHaloMassRelation::COSMOS2012($galacticusFileName,1,%options);

exit;
