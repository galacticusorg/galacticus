#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
use Carp 'verbose';
use Stats::Means;
use Galacticus::HDF5;
use Galacticus::Magnitudes;
use Galacticus::Path;
use Galacticus::Luminosities;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use XMP::MetaData;
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };

# Get name of input and output files.
die("Plot_SDSS_Tully_Fisher.pl <galacticusFile> <outputDir/File> [<showFit>]")
    unless ( scalar(@ARGV) == 2 || scalar(@ARGV) == 3 );
my $self           = $0;
my $galacticusFile = $ARGV[0];
my $outputTo       = $ARGV[1];
my $showFit;
if ( scalar(@ARGV) == 3 ) {
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
    $outputFile = $outputTo."/SDSS_Tully_Fisher.pdf";
}
(my $fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/$1/;

# Define magnitude bins.
my $magnitudePoints = pdl 10;
my $magnitudeMin    = pdl -24;
my $magnitudeMax    = pdl -19;
my $magnitudeBin    = pdl ($magnitudeMax-$magnitudeMin)/$magnitudePoints;
my $magnitudeBins   = pdl (0..$magnitudePoints-1)*$magnitudeBin+$magnitudeMin+0.5*$magnitudeBin;

# Create data structure to read the results.
my $dataSet;
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&Galacticus::HDF5::Get_Parameters($dataSet);
&Galacticus::HDF5::Count_Trees($dataSet);
&Galacticus::HDF5::Select_Output($dataSet,0.1);
$dataSet->{'tree'} = "all";
&Galacticus::HDF5::Get_Dataset($dataSet,['mergerTreeWeight','magnitudeTotal:SDSS_i:observed:z0.1000:dustAtlas[faceOn]:AB','bulgeToTotalLuminosities:SDSS_i:observed:z0.1000:dustAtlas','diskVelocity']);
my $dataSets     = $dataSet->{'dataSets'};
my $magnitude    = $dataSets->{'magnitudeTotal:SDSS_i:observed:z0.1000:dustAtlas[faceOn]:AB'};
my $bulgeToTotal = $dataSets->{'bulgeToTotalLuminosities:SDSS_i:observed:z0.1000:dustAtlas'};
my $velocity     = $dataSets->{'diskVelocity'};
my $weight       = $dataSets->{'mergerTreeWeight'};
delete($dataSet->{'dataSets'});
# Select galaxies which are disk-dominated.
my $selection         = which ($bulgeToTotal < 0.3);
# Create subsets of the galaxy properties including only the disk-dominated galaxies.
my $magnitudeSelected = $magnitude->index($selection);
my $velocitySelected  = $velocity ->index($selection);
my $weightSelected    = $weight   ->index($selection);
(my $velocityMeanGalacticus, my $velocityMeanErrorGalacticus, my $velocitySigmaGalacticus, my $velocitySigmaErrorGalacticus)
    = &Stats::Means::BinnedMean($magnitudeBins,$magnitudeSelected,$velocitySelected,$weightSelected);

# Read the XML data file.
my $xml     = new XML::Simple;
my $data    = $xml->XMLin(&galacticusPath()."data/observations/tullyFisherRelation/Tully_Fisher_SDSS_Pizagno_2007.xml");
my $columns = $data->{'tullyFisher'}->{'columns'};
my $x       = pdl @{$columns->{'magnitude'}->{'data'}};
$x          = $x-5.0*log10($columns->{'magnitude'}->{'hubble'}/$dataSet->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant'}->{'value'});
my $y       = pdl @{$columns->{'diskVelocity'}->{'data'}};
my $xError  = pdl @{$columns->{'magnitudeError'}->{'data'}};
my $yError  = pdl @{$columns->{'diskVelocityError'}->{'data'}};
my $yWeight = 1.0/$yError**2;
(my $velocityMean, my $velocityMeanError, my $velocitySigma, my $velocitySigmaError)
    = &Stats::Means::BinnedMean($magnitudeBins,$x,$y,$yWeight);

# Compute chi^2.
my $nonZeroMeanError  = which ($velocityMeanError **2+$velocityMeanErrorGalacticus **2 > 0.0);
my $nonZeroSigmaError = which ($velocitySigmaError**2+$velocitySigmaErrorGalacticus**2 > 0.0);
my $degreesOfFreedom = nelem($nonZeroMeanError)+nelem($nonZeroSigmaError);
my $chiSquared = sum((($velocityMean->index($nonZeroMeanError)-$velocityMeanGalacticus->index($nonZeroMeanError))**2)
		  /($velocityMeanError->index($nonZeroMeanError)**2+$velocityMeanErrorGalacticus->index($nonZeroMeanError)**2))
    +sum((($velocitySigma->index($nonZeroSigmaError)-$velocitySigmaGalacticus->index($nonZeroSigmaError))**2)
	 /($velocitySigmaError->index($nonZeroSigmaError)**2+$velocitySigmaErrorGalacticus->index($nonZeroSigmaError)**2));
if ( $showFit == 1 ) {
    my %fitData;
    $fitData{'name'} = "Pizagno et al. (2007) SDSS Tully-Fisher relation";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

# Make plot of stellar mass function.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;
(my $plotFileTeX = $plotFile) =~ s/\.pdf$/.tex/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
print $gnuPlot "set output '".$plotFileTeX."'\n";
print $gnuPlot "set title offset 0,-1 'SDSS Tully-Fisher Relation'\n";
print $gnuPlot "set xlabel '\$^{0.1}M_\\mathrm{i}\$'\n";
print $gnuPlot "set ylabel '\$V_\\mathrm{disk}\$ [km/s]'\n";
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
&GnuPlot::PrettyPlots::Prepare_Dataset(
    \$plot,
    $x,$y,
    errorLeft  => $xError,
    errorRight => $xError,
    errorUp    => $yError,
    errorDown  => $yError,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    pointSize  => 0.5,
    color      => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
    title      => $data->{'tullyFisher'}->{'label'}
    );
&GnuPlot::PrettyPlots::Prepare_Dataset(
    \$plot,
    $magnitudeBins,$velocityMeanGalacticus,
    errorUp    => $velocitySigmaGalacticus,
    errorDown  => $velocitySigmaGalacticus,
    style      => "point",
    symbol     => [6,7],
    weight     => [5,3],
    pointSize  => 0.5,
    color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
    title      => 'Galacticus (mean+dispersion)'
    );
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
&XMP::MetaData::Write($plotFile,$galacticusFile,$self);

exit;
