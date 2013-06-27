#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{'GALACTICUS_ROOT_V092'}) ) {
    $galacticusPath = $ENV{'GALACTICUS_ROOT_V092'};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");
use PDL;
use Data::Dumper;
use XML::Simple;
require Galacticus::HDF5;

# Extracts parameter values from a Galacticus output file and writes them to an XML file in a format suitable to re-use by Galacticus.
# Andrew Benson (10-Mar-2010)

if ( $#ARGV != 1 ) {die("Usage: Extract_Parameter_File.pl <inputGalacticusFile> <outputParameterFile>")};
my $galacticusFile = $ARGV[0];
my $parametersFile = $ARGV[1];

my $dataSet;
$dataSet->{'file'} = $galacticusFile;
&HDF5::Get_Parameters($dataSet);

my %data;
my $iParameter = -1;
foreach my $parameter ( keys(%{$dataSet->{'parameters'}}) ) {
    ++$iParameter;
    my $value;
    if ( ref($dataSet->{'parameters'}->{$parameter}) eq "PDL" ) {
	$value = join(" ",list($dataSet->{'parameters'}->{$parameter}));
    } elsif ( ref($dataSet->{'parameters'}->{$parameter}) eq "PDL::Char" ) {
	my @dims = $dataSet->{'parameters'}->{$parameter}->dims();
	$value = "";
	my $join  = "";
	for (my $i=0;$i<$dims[1];++$i) {
	    $value .= $join.$dataSet->{'parameters'}->{$parameter}->atstr($i);
	    $join   = " ";
	}
    } else {
	$value = $dataSet->{'parameters'}->{$parameter};
    }
    ${$data{'parameter'}}[$iParameter] = {
	"name"  => $parameter,
	"value" => $value
    };
}
my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
open(oHndl,">".$parametersFile);
print oHndl $xmlOutput->XMLout(\%data);
close(oHndl);

exit;
