#!/usr/bin/env perl
use XML::Simple;
use Data::Dumper;

# Compile goodness of fits for Galacticus model.
# Andrew Benson (01-Mar-2010)

# Get the name of the Galacticus file to analyze.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Galacticus_Compute_Fit.pl <galacticusFile> <outputDirectory> [<analysisScript>]")};
$galacticusFile  = $ARGV[0];
$outputDirectory = $ARGV[1];
$analysisScript  = $galacticusPath."data/analyses/Galacticus_Compute_Fit_Analyses.xml";
$analysisScript  = $ARGV[2] if ( $#ARGV == 2 );
system("mkdir -p $outputDirectory");

# Open the descriptor file that explains what analysis files to run.
$xml = new XML::Simple;
$data = $xml->XMLin($analysisScript);

# Initialize net chi^2 variables.
$chiSquaredNet       = 0.0;
$degreesOfFreedomNet = 0.0;

# Loop through all analyses and accumulate fit results.
my @analyses;
if ( UNIVERSAL::isa($data->{'analysis'},"ARRAY") ) {
    push(@analyses,@{$data->{'analysis'}});
} else {
    push(@analyses,$data->{'analysis'});
}
foreach $analysis ( @analyses ) {
    print "Running analysis script: ".$analysis->{'script'}."\n";
    $fitXML = "";
    $inXML = 0;
    system($analysis->{'script'}." ".$galacticusFile." ".$outputDirectory."/".$analysis->{'name'}.".xml ".$outputDirectory."/".$analysis->{'name'}.".pdf");
    open(iHndl,$outputDirectory."/".$analysis->{'name'}.".xml");
    while ( $line = <iHndl> ) {
	if ( $line =~ m/^\s*<constraint>/   ) {$inXML = 1};
	if ( $inXML == 1 ) {$fitXML .= $line};
	if ( $line =~ m/<\/constraint>\s*$/ ) {$inXML = 1};
    }
    close(iHndl);
    undef($fitData);
    if ( $fitXML =~ m/<constraint>/ ) {
	$fitData = $xml->XMLin($fitXML);
	$fitData->{'weight'} = $analysis->{'weight'};
	${$fitsData->{'galacticusFit'}}[++$#{$fitsData->{'galacticusFit'}}] = $fitData;
	if      ( exists($fitData->{'chiSquared'  }) ) {
	    $chiSquaredNet       +=      $fitData->{'chiSquared'      }*$fitData->{'weight'};
	    $degreesOfFreedomNet +=      $fitData->{'degreesOfFreedom'}*$fitData->{'weight'};
	} elsif ( exists($fitData->{'logLikelihood'}) ) {
	    $chiSquaredNet       += -2.0*$fitData->{'logLikelihood'   }*$fitData->{'weight'};
	    $degreesOfFreedomNet +=  1.0                               *$fitData->{'weight'};;    
	}
    }
}

# Store the net results.
$reducedChiSquaredNet = $chiSquaredNet/$degreesOfFreedomNet;
undef($fitData);
$fitData->{'chiSquared'       } = $chiSquaredNet;
$fitData->{'reducedChiSquared'} = $reducedChiSquaredNet;
$fitData->{'degreesOfFreedom' } = $degreesOfFreedomNet;
$fitData->{'name'             } = "net";
${$fitsData->{'galacticusFit'}}[++$#{$fitsData->{'galacticusFit'}}] = $fitData;

# Output the accumulated results.
$xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFits");
open(outHndl,">".$outputDirectory."/galacticusFits.xml");
print outHndl $xmlOutput->XMLout($fitsData);
close(outHndl);

exit;
