#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use Data::Dumper;

# Compile goodness of fits for Galacticus model.
# Andrew Benson (01-Mar-2010)

# Get the name of the Galacticus file to analyze.
die("Galacticus_Compute_Fit.pl <galacticusFile> <outputDirectory> [<analysisScript>]")
    unless ( scalar(@ARGV) == 2 || scalar(@ARGV) == 3 );
my $galacticusFile  = $ARGV[0];
my $outputDirectory = $ARGV[1];
my $analysisScript  = $galacticusPath."data/analyses/Galacticus_Compute_Fit_Analyses.xml";
$analysisScript  = $ARGV[2] 
    if ( scalar(@ARGV) == 3 );
system("mkdir -p $outputDirectory");

# Open the descriptor file that explains what analysis files to run.
my $xml = new XML::Simple;
my $data = $xml->XMLin($analysisScript, KeyAttr => 0);

# Initialize net chi^2 variables.
my $chiSquaredNet       = 0.0;
my $degreesOfFreedomNet = 0.0;

# Loop through all analyses and accumulate fit results.
my @analyses;
my $fitsData;
if ( UNIVERSAL::isa($data->{'analysis'},"ARRAY") ) {
    push(@analyses,@{$data->{'analysis'}});
} else {
    push(@analyses,$data->{'analysis'});
}
foreach my $analysis ( @analyses ) {
    print "Running analysis script: ".$analysis->{'script'}."\n";
    my $fitXML = "";
    my $inXML = 0;
    system($analysis->{'script'}." ".$galacticusFile." ".$outputDirectory."/".$analysis->{'name'}.".xml ".$outputDirectory."/".$analysis->{'name'}.".pdf");
    open(iHndl,$outputDirectory."/".$analysis->{'name'}.".xml");
    while ( my $line = <iHndl> ) {
	if ( $line =~ m/^\s*<constraint>/   ) {$inXML = 1};
	if ( $inXML == 1 ) {$fitXML .= $line};
	if ( $line =~ m/<\/constraint>\s*$/ ) {$inXML = 1};
    }
    close(iHndl);
    if ( $fitXML =~ m/<constraint>/ ) {
	my $fitData = $xml->XMLin($fitXML);
	$fitData->{'weight'} = $analysis->{'weight'};
	push(@{$fitsData->{'galacticusFit'}},$fitData);
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
my $reducedChiSquaredNet = $chiSquaredNet/$degreesOfFreedomNet;
my $fitData;
$fitData->{'chiSquared'       } = $chiSquaredNet;
$fitData->{'reducedChiSquared'} = $reducedChiSquaredNet;
$fitData->{'degreesOfFreedom' } = $degreesOfFreedomNet;
$fitData->{'name'             } = "net";
push(@{$fitsData->{'galacticusFit'}},$fitData);

# Output the accumulated results.
my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFits");
open(my $outHndl,">".$outputDirectory."/galacticusFits.xml");
print $outHndl $xmlOutput->XMLout($fitsData);
close($outHndl);

exit;
