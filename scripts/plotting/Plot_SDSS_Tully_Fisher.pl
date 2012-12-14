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
use PDL::NiceSlice;
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
use Carp 'verbose';
require Stats::Means;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require Galacticus::Luminosities;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_SDSS_Tully_Fisher.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
    $outputFile = $outputTo."/SDSS_Tully_Fisher.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Define magnitude bins.
$magnitudePoints = pdl 10;
$magnitudeMin    = pdl -24;
$magnitudeMax    = pdl -19;
$magnitudeBin    = pdl ($magnitudeMax-$magnitudeMin)/$magnitudePoints;
$magnitudeBins   = pdl (0..$magnitudePoints-1)*$magnitudeBin+$magnitudeMin+0.5*$magnitudeBin;

# Create data structure to read the results.
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.1);
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,['volumeWeight','magnitudeTotal:SDSS_i:observed:z0.1000:dustAtlas[faceOn]:AB','bulgeToTotalLuminosities:SDSS_i:observed:z0.1000:dustAtlas','diskVelocity']);
$dataSets     = $dataSet->{'dataSets'};
$magnitude    = $dataSets->{'magnitudeTotal:SDSS_i:observed:z0.1000:dustAtlas[faceOn]:AB'};
$bulgeToTotal = $dataSets->{'bulgeToTotalLuminosities:SDSS_i:observed:z0.1000:dustAtlas'};
$velocity     = $dataSets->{'diskVelocity'};
$weight       = $dataSets->{'volumeWeight'};
delete($dataSet->{'dataSets'});
# Select galaxies which are disk-dominated.
$selection         = which ($bulgeToTotal < 0.3);
# Create subsets of the galaxy properties including only the disk-dominated galaxies.
$magnitudeSelected = $magnitude->index($selection);
$velocitySelected  = $velocity ->index($selection);
$weightSelected    = $weight   ->index($selection);
($velocityMeanGalacticus,$velocityMeanErrorGalacticus,$velocitySigmaGalacticus,$velocitySigmaErrorGalacticus)
    = &Means::BinnedMean($magnitudeBins,$magnitudeSelected,$velocitySelected,$weightSelected);

# Read the XML data file.
$xml     = new XML::Simple;
$data    = $xml->XMLin($galacticusPath."data/observations/tullyFisherRelation/Tully_Fisher_SDSS_Pizagno_2007.xml");
$columns = $data->{'tullyFisher'}->{'columns'};
$x       = pdl @{$columns->{'magnitude'}->{'data'}};
$x       = $x-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'});
$y       = pdl @{$columns->{'diskVelocity'}->{'data'}};
$xError  = pdl @{$columns->{'magnitudeError'}->{'data'}};
$yError  = pdl @{$columns->{'diskVelocityError'}->{'data'}};
$yWeight = 1.0/$yError**2;
($velocityMean,$velocityMeanError,$velocitySigma,$velocitySigmaError)
    = &Means::BinnedMean($magnitudeBins,$x,$y,$yWeight);

# Compute chi^2.
$nonZeroMeanError  = which ($velocityMeanError **2+$velocityMeanErrorGalacticus **2 > 0.0);
$nonZeroSigmaError = which ($velocitySigmaError**2+$velocitySigmaErrorGalacticus**2 > 0.0);
$degreesOfFreedom = nelem($nonZeroMeanError)+nelem($nonZeroSigmaError);
$chiSquared = sum((($velocityMean->index($nonZeroMeanError)-$velocityMeanGalacticus->index($nonZeroMeanError))**2)
		  /($velocityMeanError->index($nonZeroMeanError)**2+$velocityMeanErrorGalacticus->index($nonZeroMeanError)**2))
    +sum((($velocitySigma->index($nonZeroSigmaError)-$velocitySigmaGalacticus->index($nonZeroSigmaError))**2)
	 /($velocitySigmaError->index($nonZeroSigmaError)**2+$velocitySigmaErrorGalacticus->index($nonZeroSigmaError)**2));
if ( $showFit == 1 ) {
    $fitData{'name'} = "Pizagno et al. (2007) SDSS Tully-Fisher relation";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make plot of stellar mass function.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'SDSS Tully-Fisher Relation'\n";
print $gnuPlot "set xlabel 'SDSS i-band absolute magnitude; \$^{0.1}M_{\\rm i}\$'\n";
print $gnuPlot "set ylabel 'Disk rotation speed; \$V_{\\rm disk}\$ [km/s]'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale y\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [-24:-18]\n";
print $gnuPlot "set yrange [30:400]\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$y,
    errorLeft  => $xError,
    errorRight => $xError,
    errorUp    => $yError,
    errorDown  => $yError,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[0]},
    title      => $data->{'tullyFisher'}->{'label'}.' [observed]'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    $magnitudeBins,$velocityMeanGalacticus,
    errorUp    => $velocitySigmaGalacticus,
    errorDown  => $velocitySigmaGalacticus,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{'redYellow'},
    title      => 'Galacticus (mean+dispersion)'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);
&MetaData::Write($plotFile,$galacticusFile,$self);

exit;
