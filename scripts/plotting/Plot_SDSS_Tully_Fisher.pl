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
use Graphics::GnuplotIF;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require Galacticus::Luminosities;
use Math::SigFigs;
require Stats::Means;
use Data::Dumper;
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_SDSS_Tully_Fisher.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
&HDF5::Get_Dataset($dataSet,['volumeWeight','magnitudeTotal:SDSS_i:observed:z0.1000:dustAtlas[faceOn]:AB','bulgeToTotalLuminosity:SDSS_i:observed:z0.1000:dustAtlas','diskCircularVelocity']);
$dataSets     = $dataSet->{'dataSets'};
$magnitude    = $dataSets->{'magnitudeTotal:SDSS_i:observed:z0.1000:dustAtlas[faceOn]:AB'};
$bulgeToTotal = $dataSets->{'bulgeToTotalLuminosity:SDSS_i:observed:z0.1000:dustAtlas'};
$velocity     = $dataSets->{'diskCircularVelocity'};
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
$data    = $xml->XMLin($galacticusPath."data/SDSS_Tully_Fisher.xml");
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

# Make the plot.
$plot1  = Graphics::GnuplotIF->new();
$plot1->gnuplot_hardcopy( '| ps2pdf - '.$outputFile, 
			  'postscript enhanced', 
			  'color lw 3 solid' );
$plot1->gnuplot_set_xlabel("^{0.1}i");
$plot1->gnuplot_set_ylabel("V_{disk, 2.2 scale lengths} [km/s]");
$plot1->gnuplot_set_title("SDSS Tully-Fisher Relation");
$plot1->gnuplot_cmd("set label \"{/Symbol c}^2=".FormatSigFigs($chiSquared,4)." [".$degreesOfFreedom."]\" at screen 0.6, screen 0.2");
$plot1->gnuplot_cmd("set xrange [-24:-18]");
$plot1->gnuplot_cmd("set yrange [30:400]");
$plot1->gnuplot_cmd("set logscale y");
$plot1->gnuplot_cmd("set mxtics 2");
$plot1->gnuplot_cmd("set mytics 10");
$plot1->gnuplot_cmd("set format y \"10^{\%L}\"");
$plot1->gnuplot_cmd("set pointsize 1.0");
$plot1->gnuplot_cmd("plot '-' with xyerrorbars pt 6 title \"".$data->{'tullyFisher'}->{'label'}."\", '-' with errorbars pt 4 title \"Galacticus (mean+dispersion)\"");
for ($i=0;$i<nelem($x);++$i) {
   $plot1->gnuplot_cmd($x->index($i)." ".$y->index($i)." ".$xError->index($i)." ".$yError->index($i));
}
$plot1->gnuplot_cmd("e");
for ($i=0;$i<nelem($magnitudeBins);++$i) {
   $plot1->gnuplot_cmd($magnitudeBins->index($i)." ".$velocityMeanGalacticus->index($i)." ".$velocitySigmaGalacticus->index($i));
}
$plot1->gnuplot_cmd("e");

exit;
