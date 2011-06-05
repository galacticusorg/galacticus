#!/usr/bin/env perl

# Extracts parameter values from a Galacticus output file and writes them to an XML file in a format suitable to re-use by Galacticus.
# Andrew Benson (10-Mar-2010)

if ( $#ARGV != 1 ) {die("Usage: Extract_Parameter_File.pl <inputGalacticusFile> <outputParameterFile>")};
$galacticusFile = $ARGV[0];
$parametersFile = $ARGV[1];

open(oHndl,">".$parametersFile);
print oHndl "<parameters>\n";
open(pHndl,"h5ls -d ".$galacticusFile."/Parameters|");
while ( $line = <pHndl> ) {
    @columns = split(/\s+/,$line);
    $name = $columns[0];
    $line = <pHndl>;
    $line = <pHndl>;
    if ( $line =~ m/^\s*\(\d+\)\s*(.*)/ ) {
	$value = $1;
	$value =~ s/\s*\"\s*,\s*\"\s*/ /g;
	$value =~ s/\s*\,\s*/ /g;
	$value =~ s/\s*\"\s*//g;
	$value =~ s/&/&amp;/g;
	print oHndl "  <parameter>\n";
	print oHndl "    <name>".$name."</name>\n";
	print oHndl "    <value>".$value."</value>\n";
	print oHndl "  </parameter>\n";
    }
}
close(pHndl);
print oHndl "</parameters>\n";
close(oHndl);

exit;
