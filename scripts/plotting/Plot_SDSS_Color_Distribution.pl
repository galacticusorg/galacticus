#!/usr/bin/env perl
use strict;
use warnings;
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
use Data::Dumper;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require System::Redirect;
require XMP::MetaData;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_SDSS_Colors_Distribution.pl <galacticusFile> <outputDir/File> [<showFit>]")};
my $self           = $0;
my $galacticusFile = $ARGV[0];
my $outputTo       = $ARGV[1];
my $showFit;
if ( $#ARGV == 2 ) {
    $showFit    = $ARGV[2];
    if ( lc($showFit) eq "showfit"   ) {$showFit = 1};
    if ( lc($showFit) eq "noshowfit" ) {$showFit = 0};
} else {
    $showFit = 0;
}

# Check if output location is file or directory.
my $outputFile;
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/SDSS_Color_Distribution.pdf";
}
(my $fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/$1/;

# Create data structure to read the results.
my $dataSet;
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.1);

# Read the XML data file.
my $xml = new XML::Simple;
my $data = $xml->XMLin("data/observations/galaxyColors/Galaxy_Colors_SDSS_Weinmann_2006.xml");
my $columns = $data->{'galaxyColors'}->{'columns'};
my $magnitude = pdl @{$columns->{'magnitude'}->{'data'}};
my $color = pdl @{$columns->{'color'}->{'data'}};
$magnitude = $magnitude-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'});

# Bin data into grid.
my $magnitudePoints = pdl 20;
my $magnitudeMin    = pdl -24;
my $magnitudeMax    = pdl -17;
my $magnitudeBin    = pdl ($magnitudeMax-$magnitudeMin)/$magnitudePoints;
my $magnitudeBins   = pdl (0..$magnitudePoints-1)*$magnitudeBin+$magnitudeMin+0.5*$magnitudeBin;
my $colorPoints     = pdl 20;
my $colorMin        = pdl 0;
my $colorMax        = pdl 1.3;
my $colorBin        = pdl ($colorMax-$colorMin)/$colorPoints;
my $colorBins       = pdl (0..$colorPoints-1)*$colorBin+$colorMin;

# Bin the SDSS data.
my $countSDSS = PDL->zeroes($magnitudePoints,$colorPoints);
my $errorSDSS = PDL->zeroes($magnitudePoints,$colorPoints);
for(my $iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    my $minimumMagnitude = $magnitudeBins->index($iMagnitude)-0.5*$magnitudeBin;
    my $maximumMagnitude = $minimumMagnitude+$magnitudeBin;
    for(my $iColor=0;$iColor<$colorPoints;++$iColor) {
	my $minimumColor = $colorBins->index($iColor)-0.5*$colorBin;
	my $maximumColor = $minimumColor+$colorBin;
	$countSDSS($iMagnitude,$iColor) += nelem(where($magnitude,($magnitude >= $minimumMagnitude) & ($magnitude < $maximumMagnitude) & ($color >= $minimumColor) & ($color < $maximumColor)))
    }
}
for(my $iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    $errorSDSS($iMagnitude,:) .= sqrt($countSDSS($iMagnitude,:))/(sum($countSDSS($iMagnitude,:))*$magnitudeBin*$colorBin);
    $countSDSS($iMagnitude,:) .= $countSDSS($iMagnitude,:)/(sum($countSDSS($iMagnitude,:))*$magnitudeBin*$colorBin);
}

# Bin the Galacticus data.
my $countGalacticus = PDL->zeroes($magnitudePoints,$colorPoints);
my $errorGalacticus = PDL->zeroes($magnitudePoints,$colorPoints);

$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['mergerTreeWeight','magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB','magnitudeTotal:SDSS_g:observed:z0.1000:dustAtlas:AB']);
my $dataSets  = $dataSet->{'dataSets'};
my $magnitudeGalacticus = $dataSets->{'magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB'};
my $colorGalacticus     = $dataSets->{'magnitudeTotal:SDSS_g:observed:z0.1000:dustAtlas:AB'}-$dataSets->{'magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB'};
my $weight    = $dataSets->{'mergerTreeWeight'};

for(my $iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    my $minimumMagnitude = $magnitudeBins->index($iMagnitude)-0.5*$magnitudeBin;
    my $maximumMagnitude = $minimumMagnitude+$magnitudeBin;
    for(my $iColor=0;$iColor<$colorPoints;++$iColor) {
	my $minimumColor = $colorBins->index($iColor)-0.5*$colorBin;
	my $maximumColor = $minimumColor+$colorBin;
	$countGalacticus($iMagnitude,$iColor) .= sum(where($weight,($magnitudeGalacticus >= $minimumMagnitude) & ($magnitudeGalacticus < $maximumMagnitude) & ($colorGalacticus >= $minimumColor) & ($colorGalacticus < $maximumColor)));
	$errorGalacticus($iMagnitude,$iColor) .= sum(where($weight**2,($magnitudeGalacticus >= $minimumMagnitude) & ($magnitudeGalacticus < $maximumMagnitude) & ($colorGalacticus >= $minimumColor) & ($colorGalacticus < $maximumColor)));
    }
}

for(my $iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    if ( sum($countGalacticus($iMagnitude,:)) > 0.0 ) {
	$errorGalacticus($iMagnitude,:) .= sqrt($errorGalacticus($iMagnitude,:))/(sum($countGalacticus($iMagnitude,:))*$magnitudeBin*$colorBin);
	$countGalacticus($iMagnitude,:) .= $countGalacticus($iMagnitude,:)/(sum($countGalacticus($iMagnitude,:))*$magnitudeBin*$colorBin);
    } else {
	$errorGalacticus($iMagnitude,:) .= 0.0;
	$countGalacticus($iMagnitude,:) .= 0.0;
    }
}

# Compute chi^2.
my $chiSquaredList = where((($countSDSS-$countGalacticus)**2)/($errorSDSS**2+$errorGalacticus**2),($errorSDSS > 0.0) & ($errorGalacticus > 0.0));
my $chiSquared = sum($chiSquaredList);
my $degreesOfFreedom = nelem($chiSquaredList);
if ( $showFit == 1 ) {
    my %fitData;
    $fitData{'name'} = "Weinmann et al. (2006) SDSS galaxy colors";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make the plot.
open(pHndl,"|gnuplot 1>/dev/null 2>&1");

# Check which version of GnuPlot we're using.
my $gnuPlotVersionString = `gnuplot -V`;
unless ( $gnuPlotVersionString =~ m/^gnuplot\s+(\d)+\.(\d)+\s+patchlevel\s+(\d+)/ ) {die ("Plot_SDSS_Colors_Distribution.pl: unable to determine GnuPlot version")};
my @gnuPlotVersion = ( $1, $2, $3 );
my $gnuPlotNew;
if ( $gnuPlotVersion[0] < 4 || $gnuPlotVersion[0] == 4 && $gnuPlotVersion[1] <= 2 ) {
    $gnuPlotNew = 0;
    print pHndl "set terminal table\n";
    print pHndl "set output 'contour.dat'\n";
} else {
    $gnuPlotNew = 1;
    print pHndl "set table 'contour.dat'\n";
}
print pHndl "unset surface\n";
print pHndl "set contour base; set cntrparam level 10\n";
print pHndl "splot '-'\n";
for(my $iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    for(my $iColor=0;$iColor<$colorPoints;++$iColor) {
 	print pHndl $magnitudeBins->index($iMagnitude)." ".$colorBins->index($iColor)." ".$countGalacticus(($iMagnitude),($iColor))."\n";
    }
    print pHndl "\n" unless ( $iMagnitude == $magnitudePoints-1 );
}
print pHndl "e\n";
print pHndl "unset table\n" if ( $gnuPlotNew == 1);
close(pHndl);
system("awk \"NF<2{printf\\\"\\n\\\"}{print}\" <contour.dat >contour1.dat");

open(my $gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output 'tmp.eps'\n";
print $gnuPlot "set title 'SDSS Galaxy Color Distribution'\n";
print $gnuPlot "set xlabel '\$^{0.1}\$r'\n";
print $gnuPlot "set ylabel '\$^{0.1}\$g\$-^{0.1}\$r'\n";
print $gnuPlot "set xrange [".$magnitudeMin.":".$magnitudeMax."]\n";
print $gnuPlot "set yrange [".$colorMin.":".$colorMax."]\n";
print $gnuPlot "set pm3d map\n";
print $gnuPlot "set pm3d explicit\n";
print $gnuPlot "set logscale z\n";
print $gnuPlot "set palette rgbformulae 34,35,36\n";
print $gnuPlot "splot '-' with pm3d title \"".$data->{'galaxyColors'}->{'label'}."\", 'contour1.dat' with line lt -1 lc rgbcolor \"#FFFF00\" title \"Galacticus\"\n";
for(my $iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    for(my $iColor=0;$iColor<$colorPoints;++$iColor) {
 	print $gnuPlot $magnitudeBins->index($iMagnitude)." ".$colorBins->index($iColor)." ".$countSDSS(($iMagnitude),($iColor))."\n";
    }
    print $gnuPlot "\n" unless ( $iMagnitude == $magnitudePoints-1 );
}
print $gnuPlot "e\n";
close($gnuPlot);
&LaTeX::GnuPlot2PDF("tmp.eps");
unlink("tmp.eps");
unlink("contour.dat");
unlink("contour1.dat");

# Plot slices through
open(pHndl,"|gnuplot");
print pHndl "set terminal postscript enhanced color lw 3 solid\n";
print pHndl "set output 'tmp.ps'\n";
print pHndl "set title 'SDSS Galaxy Color Distribution'\n";
for(my $iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    my $minimumMagnitude = $magnitudeBins->index($iMagnitude)-0.5*$magnitudeBin;
    my $maximumMagnitude = $minimumMagnitude+$magnitudeBin;
    print pHndl "set xlabel \"^{0.1}g-^{0.1}r\"\n";
    print pHndl "set ylabel \"dn/dlog(^{0.1}g-^{0.1}r)/n\"\n";
    print pHndl "set title \"SDSS Galaxy Color Distribution for ".$minimumMagnitude." < ^{0.1}r < ".$maximumMagnitude."\"\n";
    print pHndl "set xrange [".$colorMin.":".$colorMax."]\n";
    print pHndl "set yrange [0.0:10.0]\n";
    print pHndl "set pointsize 1.0\n";
    print pHndl "plot '-' with errorbars lt 1 pt 6 title \"SDSS\", '-' with errorbars lt 2 pt 6 title \"Galacticus\"";
    print pHndl "\n";
    for(my $iColor=0;$iColor<$colorPoints;++$iColor) {
	print pHndl $colorBins->index($iColor)." ".$countSDSS(($iMagnitude),($iColor))." ".$errorSDSS(($iMagnitude),($iColor))."\n";
    }
    print pHndl "e\n";
    for(my $iColor=0;$iColor<$colorPoints;++$iColor) {
	print pHndl $colorBins->index($iColor)." ".$countGalacticus(($iMagnitude),($iColor))." ".$errorGalacticus(($iMagnitude),($iColor))."\n";
    }
    print pHndl "e\n";
}

# Convert to PDF.
close(pHndl);
system("ps2pdf tmp.ps tmp2.pdf");
unlink($outputFile);
&SystemRedirect::tofile("pdfmerge tmp.pdf tmp2.pdf ".$outputFile,"/dev/null");
unlink("tmp.ps","tmp.pdf","tmp2.pdf");
&MetaData::Write($outputFile,$galacticusFile,$self);

exit;
