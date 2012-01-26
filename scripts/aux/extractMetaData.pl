#!/usr/bin/env perl
use lib "./perl";
use strict;
use warnings;
use XMP::MetaData;

# A simple script to extract the metadata embedded into plots created by Galacticus.
die("Usage: extractMetaData.pl <plotFile> <outputPrefix>") unless ( scalar(@ARGV) == 2 );
my $plotFile     = $ARGV[0];
my $outputPrefix = $ARGV[1];

# Specify output files.
my $parameterFile = $outputPrefix."Parameters.xml";
my $scriptFile    = $outputPrefix."Script.pl";
my $mergeFile     = $outputPrefix.".merge";
my $patchFile     = $outputPrefix.".patch";

# Extract the metadata.
&MetaData::Read($plotFile,$parameterFile,$scriptFile,$mergeFile,$patchFile);

exit;
