#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Math::SigFigs;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require Stats::Percentiles;
require XMP::MetaData;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_SDSS_Gas_Metallicity.pl <galacticusFile> <outputDir/File> [<showFit>]")};
$self           = $0;
$galacticusFile = $ARGV[0];
$outputTo       = $ARGV[1];
if ( $#ARGV == 2 ) {
    $showFit    = $ARGV[2];
    if ( lc($showFit) eq "showfit"   ) {$showFit = 1};
    if ( lc($showFit) eq "noshowfit" ) {$showFit = 0};
} else {
    $showFit = 0;
}

# Check if output location is file or directory.
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/SDSS_Gas_Metallicity.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Specify Solar metallicity and oxygen abundance.
$solarMetallicity     = pdl 0.0189;
$solarOxygenAbundance = pdl 4.8977e-4;

# Minimum gas fraction for galaxies to be considered in this plot.
$gasFractionMinimum   = pdl 0.1;

# Initialize chi^2 accumulator.
$chiSquared      = 0.0;
$degreesOfFreedom = 0;

# Create data structure to read the results.
$dataBlock->{'file'}  = $galacticusFile;
$dataBlock->{'store'} = 0;
&HDF5::Get_Parameters($dataBlock    );
&HDF5::Count_Trees   ($dataBlock    );
&HDF5::Select_Output ($dataBlock,0.1);
$dataBlock->{'tree'} = "all";
&HDF5::Get_Dataset($dataBlock,['volumeWeight'
			      ,'magnitudeTotal:SDSS_g:observed:z0.1000:dustAtlas[faceOn]:AB'
			      ,'magnitudeTotal:SDSS_z:observed:z0.1000:AB'
			      ,'diskMassStellar'
			      ,'spheroidMassStellar'
			      ,'diskMassGas'
			      ,'spheroidMassGas'
			      ,'diskAbundancesGasMetals'
			      ,'spheroidAbundancesGasMetals'
		   ]);
$dataSets = $dataBlock->{'dataSets'};
$gasFraction    = ($dataSets->{'diskMassGas'}+$dataSets->{'spheroidMassGas'})/($dataSets->{'diskMassGas'}+$dataSets->{'spheroidMassGas'}+$dataSets->{'diskStellarMass'}+$dataSets->{'spheroidMassStellar'});
$gasMetallicity = where(12.0+log10(($dataSets->{'diskAbundancesGasMetals'}+$dataSets->{'spheroidAbundancesGasMetals'})/($dataSets->{'diskMassGas'}+$dataSets->{'spheroidMassGas'}))-log10($solarMetallicity)+log10($solarOxygenAbundance),$gasFraction > $gasFractionMinimum);

# Open a pipe to GnuPlot.
open(gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
print gnuPlot "set output \"tmp.ps\"\n";

# Read the XML data file.
undef(@tmpFiles);
$xml = new XML::Simple;
$data = $xml->XMLin($galacticusPath."data/observations/abundances/Gas_Phase_Metallicities_SDSS_Tremonti_2004.xml");
$iDataset = 0;
foreach $dataSet ( @{$data->{'gasMetallicity'}} ) {
    ++$iDataset;
    $columns = $dataSet->{'columns'};
    $x = pdl @{$columns->{'magnitude'}->{'data'}};
    $x = $x-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataBlock->{'parameters'}->{'H_0'});

    # Compute the distribution of Galacticus galaxies.
    $filter = $columns->{'magnitude'}->{'filter'};
    $dust   = $columns->{'magnitude'}->{'dust'};
    if ( $dust eq "corrected" ) {$dustLabel = ""};
    if ( $dust eq "face-on" )   {$dustLabel = ":dustAtlas[faceOn]"};

    $property    = "magnitudeTotal:".$filter.":observed:z0.1000".$dustLabel.":AB";
    $magnitude   = where($dataSets->{$property}     ,$gasFraction > $gasFractionMinimum);
    $weight      = where($dataSets->{'volumeWeight'},$gasFraction > $gasFractionMinimum);
    $percentiles = pdl [2.5,16.0,50.0,84.0,97.5];
    $results     = &Percentiles::BinnedPercentiles(
	$x,
	$magnitude,
	$gasMetallicity,
	$weight,
	$percentiles,
	);

    $plotCommand = "";
    $joiner = "";
    $iPercentile = 0;
    $chiSquaredRange = 0.0;
    $degreesOfFreedomRange = 0;
    $observationalError = pdl $dataSet->{'distributionError'};
    foreach $percentile ( @{$columns->{'distributionPercentile'}} ) {
	++$iPercentile;

	# Compute chi^2.
	$yData       = pdl @{$percentile->{'data'}};
	$yGalacticus = $results(:,($iPercentile-1));
	$chiSquaredRange += sum((($yData-$yGalacticus)/$observationalError)**2);
	$degreesOfFreedomRange += nelem($yData);

	$plotCommand .= $joiner."'gnuplot".$iDataset.":".$iPercentile.".tmp' lt ".$iPercentile." pt 6 title \"".$dataSet->{'label'}." [".$percentile->{'percentile'}."%]\"";
	$joiner = ", ";
	push(@tmpFiles,"gnuplot".$iDataset.":".$iPercentile.".tmp");
	open(tmpHndl,">gnuplot".$iDataset.":".$iPercentile.".tmp");
	for ($i=0;$i<nelem($x);++$i) {
	    print tmpHndl $x->index($i)." ".$yData->index($i)."\n";
	}
	close(tmpHndl);
	
	if ( any($yGalacticus > 0.0) ) {
	    $plotCommand .= $joiner."'gnuplot".$iDataset.":".$iPercentile."_glc.tmp' lt ".$iPercentile." pt 4 title \"Galacticus [".$percentile->{'percentile'}."%]\"";
	    $joiner = ", ";
	    push(@tmpFiles,"gnuplot".$iDataset.":".$iPercentile."_glc.tmp");
	    open(tmpHndl,">gnuplot".$iDataset.":".$iPercentile."_glc.tmp");
	    for ($i=0;$i<nelem($x);++$i) {
		if ( $yGalacticus->index($i) > 0.0 ) {print tmpHndl $x->index($i)." ".$yGalacticus->index($i)."\n"};
	    }
	    close(tmpHndl);
	}	
    }
    $plotCommand = "plot ".$plotCommand."\n";
    $plotCommand = "set label \"{/Symbol c}^2=".FormatSigFigs($chiSquaredRange,4)." [".$degreesOfFreedomRange."]\" at screen 0.6, screen 0.2\n".$plotCommand;
    $plotCommand = "unset label\n".$plotCommand;

    # Accumulate chi^2.
    $chiSquared += $chiSquaredRange;
    $degreesOfFreedom += $degreesOfFreedomRange;

    # Make the plot.
    print gnuPlot "set xlabel \"".$columns->{'magnitude'}->{'label'}."\"\n";
    print gnuPlot "set ylabel \"12+[O/H]\"\n";
    print gnuPlot "set title \"".$dataSet->{'title'}."\"\n";
    print gnuPlot "set key left\n";
    print gnuPlot "set key bottom\n";
    print gnuPlot "set mxtics 2\n";
    print gnuPlot "set mytics 2\n";
    print gnuPlot "set yrange [7:10]\n";
    print gnuPlot "set pointsize 1.0\n";
    print gnuPlot $plotCommand;
    
}

# Close the pipe to GnuPlot.
close(gnuPlot);

# Convert to PDF.
system("ps2pdf tmp.ps ".$outputFile);
&MetaData::Write($outputFile,$galacticusFile,$self);

# Clean up files.
unlink("tmp.ps",@tmpFiles);

# Display chi^2 information
if ( $showFit == 1 ) {
    $fitData{'name'} = "Tremonti et al. (2004) SDSS gas-phase metallicity distributions";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

exit;
