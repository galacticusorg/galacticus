#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V094'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V094'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");
use DateTime;
use File::Copy;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
require Galacticus::AGNLuminosities;

# Driver script for Hopkins (2007) AGN SED model.
# Andrew Benson (04-August-2014)

# Get arguments.
die "Usage: Hopkins2007_AGN_SED_Driver.pl <fileFormatVersion>"
    unless ( scalar(@ARGV) == 1 );
my $fileFormat = $ARGV[0];

# Ensure the requested file format version is compatible.
my $fileFormatCurrent = pdl long(1);
die('Hopkins2007_AGN_SED_Driver.pl: this script supports file format version '.$fileFormatCurrent.' but version '.$fileFormat.' was requested')
    unless ( $fileFormat == $fileFormatCurrent );

# Determine if we need to make the file.
my $fileName = $galacticusPath."data/blackHoles/AGN_SEDs_Hopkins2007.hdf5";
my $makeFile = 0;
if ( -e $fileName ) {
    my $hdfFile = new PDL::IO::HDF5($fileName);
    my @attributes = $hdfFile->attrs();
    if ( grep {$_ eq "fileFormat"} @attributes  ) {
	my @fileFormatCurrentFile = $hdfFile->attrGet('fileFormat');
	$makeFile = 1
	    unless ( $fileFormatCurrentFile[0] == $fileFormatCurrent );
    } else {
	$makeFile = 1;
    }
} else {
    $makeFile = 1;
}

# Check if we need to make the file.
&AGNLuminosities::Build_AGN_Spectra($fileFormat)
    if ( $makeFile == 1 );

exit;

