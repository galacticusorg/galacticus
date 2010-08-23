#!/usr/bin/env perl
use XML::Simple;
use Data::Dumper;

# Compile goodness of fits for Galacticus model.
# Andrew Benson (01-Mar-2010)

# Get the name of the Galacticus file to analyze.
if ( $#ARGV != 1 ) {die("Galacticus_Compute_Fit.pl <galacticusFile> <outputDirectory>")};
$galacticusFile  = $ARGV[0];
$outputDirectory = $ARGV[1];
system("mkdir -p $outputDirectory");

# Open the descriptor file that explains what analysis files to run.
$xml = new XML::Simple;
$data = $xml->XMLin("data/Galacticus_Compute_Fit_Analyses.xml");

# Initialize net chi^2 variables.
$chiSquaredNet = 0.0;
$degreesOfFreedomNet = 0.0;

# Loop through all analyses and accumulate fit results.
foreach $analysis ( @{$data->{'analysis'}} ) {
    print "Running analysis script: ".$analysis->{'script'}."\n";
    $fitXML = "";
    $inXML = 0;
    open(pipeHndl,$analysis->{'script'}." ".$galacticusFile." ".$outputDirectory." showFit |");
    while ( $line = <pipeHndl> ) {
	if ( $line =~ m/^\s*<galacticusFit>/   ) {$inXML = 1};
	if ( $inXML == 1 ) {$fitXML .= $line};
	if ( $line =~ m/<\/galacticusFit>\s*$/ ) {$inXML = 1};
    }
    close(pipeHndl);
    undef($fitData);
    if ( $fitXML =~ m/<galacticusFit>/ ) {
	$fitData = $xml->XMLin($fitXML);
	$fitData->{'weight'} = $analysis->{'weight'};
	${$fitsData->{'galacticusFit'}}[++$#{$fitsData->{'galacticusFit'}}] = $fitData;
	$chiSquaredNet += $fitData->{'chiSquared'}*$fitData->{'weight'};
	$degreesOfFreedomNet += $fitData->{'degreesOfFreedom'}*$fitData->{'weight'};
    }
}

# Store the net results.
$reducedChiSquaredNet = $chiSquaredNet/$degreesOfFreedomNet;
undef($fitData);
$fitData->{'chiSquared'} = $chiSquaredNet;
$fitData->{'reducedChiSquared'} = $reducedChiSquaredNet;
$fitData->{'degreesOfFreedom'} = $degreesOfFreedomNet;
$fitData->{'name'} = "net";
${$fitsData->{'galacticusFit'}}[++$#{$fitsData->{'galacticusFit'}}] = $fitData;

# Output the accumulated results.
$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFits");
open(outHndl,">".$outputDirectory."/galacticusFits.xml");
print outHndl $xmlOutput->XMLout($fitsData);
close(outHndl);

exit;
