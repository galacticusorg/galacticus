#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use Data::Dumper;
use XML::Simple;
use Galacticus::HDF5;

# Extracts parameter values from a Galacticus output file and writes them to an XML file in a format suitable to re-use by Galacticus.
# Andrew Benson (10-Mar-2010)

die("Usage: Extract_Parameter_File.pl <inputGalacticusFile> <outputParameterFile>")
    unless ( scalar(@ARGV) == 2 );
my $galacticusFile = $ARGV[0];
my $parametersFile = $ARGV[1];

my $dataSet;
$dataSet->{'file'} = $galacticusFile;
&Galacticus::HDF5::Get_Parameters($dataSet);

my $outputData;
my @stack = map {{name => $_, node => $dataSet->{'parameters'}->{$_}, to => \%{$outputData}}} keys(%{$dataSet->{'parameters'}});
while ( scalar(@stack) > 0 ) {
    my $parameter = pop(@stack);
    my $value;
    if ( $parameter->{'name'} =~ m/^sub:(.*)/ ) {
	my $parameterName = $1;
	push
	    (
	     @stack,
	     map {{name => $_, node => $parameter->{'node'}->{$_}, to => \%{$parameter->{'to'}->{$parameterName}}}} keys(%{$parameter->{'node'}})
	     );

    } else {
	if ( ref($parameter->{'node'}->{'value'}) eq "PDL" ) {
	    $value = join(" ",list($parameter->{'node'}->{'value'}));
	} elsif ( ref($parameter->{'node'}->{'value'}) eq "PDL::Char" ) {
	    my @dims = $parameter->{'node'}->{'value'}->dims();
	    $value = "";
	    my $join  = "";
	    for (my $i=0;$i<$dims[1];++$i) {
		(my $text = $parameter->{'node'}->{'value'}->atstr($i)) =~ s/\s+$//;
		$value .= $join.$text;
		$join   = " ";
	    }
	} else {
	    ($value = $parameter->{'node'}->{'value'}) =~ s/\s+$//;
	}
	$parameter->{'to'}->{$parameter->{'name'}}->{'value'} = $value;
    }
}
my $xmlOutput = new XML::Simple (RootName=>"parameters");
open(oHndl,">".$parametersFile);
print oHndl $xmlOutput->XMLout($outputData);
close(oHndl);

exit;
