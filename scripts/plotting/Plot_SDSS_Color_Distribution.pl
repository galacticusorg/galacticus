#!/usr/bin/env perl
use lib "./perl";
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
    $outputFile = $outputTo."/SDSS_Color_Distribution.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.1);

# Read the XML data file.
$xml = new XML::Simple;
$data = $xml->XMLin("data/observations/galaxyColors/Galaxy_Colors_SDSS_Weinmann_2006.xml");
$columns = $data->{'galaxyColors'}->{'columns'};
$magnitude = pdl @{$columns->{'magnitude'}->{'data'}};
$color = pdl @{$columns->{'color'}->{'data'}};
$magnitude = $magnitude-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'});

# Bin data into grid.
$magnitudePoints = pdl 20;
$magnitudeMin    = pdl -24;
$magnitudeMax    = pdl -17;
$magnitudeBin    = pdl ($magnitudeMax-$magnitudeMin)/$magnitudePoints;
$magnitudeBins   = pdl (0..$magnitudePoints-1)*$magnitudeBin+$magnitudeMin+0.5*$magnitudeBin;
$colorPoints     = pdl 20;
$colorMin        = pdl 0;
$colorMax        = pdl 1.3;
$colorBin        = pdl ($colorMax-$colorMin)/$colorPoints+0.5*$colorBin;
$colorBins       = pdl (0..$colorPoints-1)*$colorBin+$colorMin;

# Bin the SDSS data.
$countSDSS = PDL->zeroes($magnitudePoints,$colorPoints);
$errorSDSS = PDL->zeroes($magnitudePoints,$colorPoints);
for($iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    $minimumMagnitude = $magnitudeBins->index($iMagnitude)-0.5*$magnitudeBin;
    $maximumMagnitude = $minimumMagnitude+$magnitudeBin;
    for($iColor=0;$iColor<$colorPoints;++$iColor) {
	$minimumColor = $colorBins->index($iColor)-0.5*$colorBin;
	$maximumColor = $minimumColor+$colorBin;
	$countSDSS($iMagnitude,$iColor) += nelem(where($magnitude,$magnitude >= $minimumMagnitude & $magnitude < $maximumMagnitude & $color >= $minimumColor & $color < $maximumColor))
    }
}
for($iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    $errorSDSS($iMagnitude,:) .= sqrt($countSDSS($iMagnitude,:))/(sum($countSDSS($iMagnitude,:))*$magnitudeBin*$colorBin);
    $countSDSS($iMagnitude,:) .= $countSDSS($iMagnitude,:)/(sum($countSDSS($iMagnitude,:))*$magnitudeBin*$colorBin);
}

# Bin the Galacticus data.
$countGalacticus = PDL->zeroes($magnitudePoints,$colorPoints);
$errorGalacticus = PDL->zeroes($magnitudePoints,$colorPoints);

$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['volumeWeight','magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB','magnitudeTotal:SDSS_g:observed:z0.1000:dustAtlas:AB']);
$dataSets  = $dataSet->{'dataSets'};
$magnitude = $dataSets->{'magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB'};
$color     = $dataSets->{'magnitudeTotal:SDSS_g:observed:z0.1000:dustAtlas:AB'}-$dataSets->{'magnitudeTotal:SDSS_r:observed:z0.1000:dustAtlas:AB'};
$weight    = $dataSets->{'volumeWeight'};

for($iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    $minimumMagnitude = $magnitudeBins->index($iMagnitude)-0.5*$magnitudeBin;
    $maximumMagnitude = $minimumMagnitude+$magnitudeBin;
    for($iColor=0;$iColor<$colorPoints;++$iColor) {
	$minimumColor = $colorBins->index($iColor)-0.5*$colorBin;
	$maximumColor = $minimumColor+$colorBin;
	$countGalacticus($iMagnitude,$iColor) .= sum(where($weight,$magnitude >= $minimumMagnitude & $magnitude < $maximumMagnitude & $color >= $minimumColor & $color < $maximumColor));
	$errorGalacticus($iMagnitude,$iColor) .= sum(where($weight**2,$magnitude >= $minimumMagnitude & $magnitude < $maximumMagnitude & $color >= $minimumColor & $color < $maximumColor));
    }
}

for($iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    if ( sum($countGalacticus($iMagnitude,:)) > 0.0 ) {
	$errorGalacticus($iMagnitude,:) .= sqrt($errorGalacticus($iMagnitude,:))/(sum($countGalacticus($iMagnitude,:))*$magnitudeBin*$colorBin);
	$countGalacticus($iMagnitude,:) .= $countGalacticus($iMagnitude,:)/(sum($countGalacticus($iMagnitude,:))*$magnitudeBin*$colorBin);
    } else {
	$errorGalacticus($iMagnitude,:) .= 0.0;
	$countGalacticus($iMagnitude,:) .= 0.0;
    }
}

# Compute chi^2.
$chiSquaredList = where((($countSDSS-$countGalacticus)**2)/($errorSDSS**2+$errorGalacticus**2),$errorSDSS > 0.0 & $errorGalacticus > 0.0);
$chiSquared = sum($chiSquaredList);
$degreesOfFreedom = nelem($chiSquaredList);
if ( $showFit == 1 ) {
    $fitData{'name'} = "Weinmann et al. (2006) SDSS galaxy colors";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make the plot.
open(pHndl,"|gnuplot 1>/dev/null 2>&1");

# Check which version of GnuPlot we're using.
$gnuPlotVersionString = `gnuplot -V`;
unless ( $gnuPlotVersionString =~ m/^gnuplot\s+(\d)+\.(\d)+\s+patchlevel\s+(\d+)/ ) {die ("Plot_SDSS_Colors_Distribution.pl: unable to determine GnuPlot version")};
@gnuPlotVersion = ( $1, $2, $3 );
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
for($iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    for($iColor=0;$iColor<$colorPoints;++$iColor) {
 	print pHndl $magnitudeBins->index($iMagnitude)." ".$colorBins->index($iColor)." ".$countGalacticus(($iMagnitude),($iColor))."\n";
    }
    print pHndl "\n" unless ( $iMagnitude == $magnitudePoints-1 );
}
print pHndl "e\n";
print pHndl "unset table\n" if ( $gnuPlotNew == 1);
close(pHndl);
system("awk \"NF<2{printf\\\"\\n\\\"}{print}\" <contour.dat >contour1.dat");

open($gnuPlot,"|gnuplot");
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
for($iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    for($iColor=0;$iColor<$colorPoints;++$iColor) {
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
for($iMagnitude=0;$iMagnitude<$magnitudePoints;++$iMagnitude) {
    $minimumMagnitude = $magnitudeBins->index($iMagnitude)-0.5*$magnitudeBin;
    $maximumMagnitude = $minimumMagnitude+$magnitudeBin;
    print pHndl "set xlabel \"^{0.1}g-^{0.1}r\"\n";
    print pHndl "set ylabel \"dn/dlog(^{0.1}g-^{0.1}r)/n\"\n";
    print pHndl "set title \"SDSS Galaxy Color Distribution for ".$minimumMagnitude." < ^{0.1}r < ".$maximumMagnitude."\"\n";
    print pHndl "set xrange [".$colorMin.":".$colorMax."]\n";
    print pHndl "set yrange [0.0:10.0]\n";
    print pHndl "set pointsize 1.0\n";
    print pHndl "plot '-' with errorbars lt 1 pt 6 title \"SDSS\", '-' with errorbars lt 2 pt 6 title \"Galacticus\"";
    print pHndl "\n";
    for($iColor=0;$iColor<$colorPoints;++$iColor) {
	print pHndl $colorBins->index($iColor)." ".$countSDSS(($iMagnitude),($iColor))." ".$errorSDSS(($iMagnitude),($iColor))."\n";
    }
    print pHndl "e\n";
    for($iColor=0;$iColor<$colorPoints;++$iColor) {
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
