#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
require Galacticus::HDF5;
use Data::Dumper;
use XML::Simple;

# Extracts parameter values from a Galacticus output file and writes them to an XML file in a format suitable to re-use by Galacticus.
# Andrew Benson (10-Mar-2010)

if ( $#ARGV != 1 ) {die("Usage: Extract_Parameter_File.pl <inputGalacticusFile> <outputParameterFile>")};
$galacticusFile = $ARGV[0];
$parametersFile = $ARGV[1];

$dataSet{'file'} =$galacticusFile;
&HDF5::Get_Parameters(\%dataSet);

$iParameter = -1;
foreach $parameter ( keys(%{$dataSet{'parameters'}}) ) {
    ++$iParameter;
    if ( ref(${$dataSet{'parameters'}}{$parameter}) eq "PDL" ) {
	$value = join(" ",list(${$dataSet{'parameters'}}{$parameter}));
    } elsif ( ref(${$dataSet{'parameters'}}{$parameter}) eq "PDL::Char" ) {
	@dims = ${$dataSet{'parameters'}}{$parameter}->dims;
	$value = "";
	$join  = "";
	for ($i=0;$i<$dims[1];++$i) {
	    $value .= $join.${$dataSet{'parameters'}}{$parameter}->atstr($i);
	    $join   = " ";
	}
    } else {
	$value = ${$dataSet{'parameters'}}{$parameter};
    }
    ${$data{'parameter'}}[$iParameter] = {
	"name"  => $parameter,
	"value" => $value
    };
}
$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"parameters");
open(oHndl,">".$parametersFile);
print oHndl $xmlOutput->XMLout(\%data);
close(oHndl);

exit;
