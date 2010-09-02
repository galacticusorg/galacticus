#!/usr/bin/env perl
use lib "./perl";
use PDL;
use PDL::IO::HDF5;
use PDL::IO::HDF5::Dataset;
use XML::Simple;
use Galacticus::HDF5;
use Galacticus::Magnitudes;
use Galacticus::HaloModel;
use Math::SigFigs;
use Data::Dumper;

# Compute and plot a 2-point correlation functions and compare to 2dFGRS observations.
# Andrew Benson (30-Aug-2010)

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_2dFGRS_Correlation_Functions.pl <galacticusFile> <outputDir/File> [<showFit>]")};
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
    $outputFile = $outputTo."/2dFGRS_Correlation_Functions.pdf";
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataSet{'file'}  = $galacticusFile;
$dataSet{'store'} = 0;
$dataSet{'tree'}  = "all";
&HDF5::Get_Parameters(\%dataSet    );
&HDF5::Count_Trees   (\%dataSet    );
&HDF5::Get_Times     (\%dataSet    );
&HDF5::Select_Output (\%dataSet,0.1);
&HDF5::Get_Dataset   (\%dataSet,['magnitudeTotal:bJ:observed:z0.1000:dustAtlas:vega','nodeBias']);
$dataSets         = \%{$dataSet{'dataSets'}};

# Read the file of observational data.
$xml     = new XML::Simple;
$data    = $xml->XMLin("data/Correlation_Functions_2dFGRS_Norberg_2002.xml");

# Open a pipe to GnuPlot.
open(gnuPlot,"|gnuplot");
print gnuPlot "set terminal postscript enhanced color lw 3 solid\n";
print gnuPlot "set output \"tmp.ps\"\n";

# Loop over correlation functions.
foreach $correlationFunction ( @{$data->{'correlationFunction'}} ) {
    # Skip correlation functions with a color selection and those that are not redshift space.
    unless ( exists($correlationFunction->{'colorRange'}) || $correlationFunction->{'space'} ne "redshift" ) {
	# Get magnitude ranges for this sample.
	$magnitudeMinimum = $correlationFunction->{'magnitudeRange'}->{'minimum'}
	-5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
	$magnitudeMaximum = $correlationFunction->{'magnitudeRange'}->{'maximum'}
	-5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet{'parameters'}->{'H_0'});
	# Get separation, correlation function and errors.
	$separationData   = pdl @{$correlationFunction->{'separation'                          }->{'datum'}};
	$xiData           = pdl @{$correlationFunction->{'correlationFunction'                 }->{'datum'}};
	$xiBootData       = pdl @{$correlationFunction->{'correlationFunctionBootstrapped'     }->{'datum'}};
	$xiBootErrorData  = pdl @{$correlationFunction->{'correlationFunctionBootstrappedError'}->{'datum'}};
	# Convert separation for Hubble constant.
	$separationData  *= ($dataSet{'parameters'}->{'H_0'}/$correlationFunction->{'separation'}->{'hubble'})**$correlationFunction->{'separation'}->{'hubbleExponent'};
	# Get error on actual correlation function.
	$xiErrorData      = $xiData*$xiBootErrorData/$xiBootData;

	# Select a matching subset of model galaxies.
	$selected         = which(${$dataSets->{'magnitudeTotal:bJ:observed:z0.1000:dustAtlas:vega'}} >= $magnitudeMinimum
				  & ${$dataSets->{'magnitudeTotal:bJ:observed:z0.1000:dustAtlas:vega'}} < $magnitudeMaximum);

	# Skip empty selections.
	unless ( nelem($selected) == 0 ) {

	    # Get the power spectrum for these galaxies.
	    ($waveNumber,$linearPowerSpectrum,$galaxyPowerSpectrum) = &HaloModel::Compute_Power_Spectrum(\%dataSet
													 ,$selected
													 ,space => "redshift");
	    # Compute the two-point correlation functions.
	    $separationMinimum         = $separationData->index(0);
	    $separationMaximum         = $separationData->index(nelem($separationData)-1);
	    $separationPointsPerDecade = (nelem($separationData)-1)/log10($separationMaximum/$separationMinimum);
	    ($separations,$galaxyCorrelationFunction) = &HaloModel::Compute_Correlation_Function( $waveNumber
												  ,$galaxyPowerSpectrum
												  ,$separationMinimum
												  ,$separationMaximum
												  ,$separationPointsPerDecade);
	    
	    # Compute chi^2.
	    $thisChiSquared        = sum(($galaxyCorrelationFunction-$xiData)**2/($xiErrorData**2));
	    $thisDegreesOfFreedom  = nelem($xiData);
	    $chiSquared           += $thisChiSquared;
	    $degreesOfFreedom     += $thisDegreesOfFreedom;
	    
	    # Make the plot.
	    print gnuPlot "set xlabel \"s [Mpc]\"\n";
	    print gnuPlot "set ylabel \"{/Symbol x}(s)\"\n";
	    print gnuPlot "set title \"Galaxy redshift space correlation function (".FormatSigFigs($magnitudeMinimum,4)." < M_{b_J} < ".FormatSigFigs($magnitudeMaximum,4).")\"\n";
	    print gnuPlot "unset label\n";
	    print gnuPlot "set label \"{/Symbol c}^2=".FormatSigFigs($thisChiSquared,4)." [".$thisDegreesOfFreedom."]\" at screen 0.6, screen 0.2\n";
	    print gnuPlot "set key right\n";
	    print gnuPlot "set logscale xy\n";
	    print gnuPlot "set mxtics 10\n";
	    print gnuPlot "set mytics 10\n";
	    print gnuPlot "set format x \"10^{\%L}\"\n";
	    print gnuPlot "set format y \ \"10^{\%L}\"\n";
	    print gnuPlot "set pointsize 1.0\n";
	    $plotCommand  = "plot '-' with errorbars pt 6 title \"".$data->{'label'}."\"";
	    $plotCommand .= ", '-' title \"Galacticus\"";
	    print gnuPlot $plotCommand."\n";
	    for ($i=0;$i<nelem($separationData);++$i) {
		print gnuPlot $separationData->index($i)." ".$xiData->index($i)." ".$xiErrorData->index($i)."\n";
	    }
	    print gnuPlot "e\n";
	    for ($i=0;$i<nelem($separations);++$i) {
		print gnuPlot $separations->index($i)." ".$galaxyCorrelationFunction->index($i)."\n";
	    }
	    print gnuPlot "e\n";
	}
    }
}

# Close the pipe to GnuPlot.
close(gnuPlot);

# Convert to PDF.
system("ps2pdf tmp.ps ".$outputFile);

# Clean up files.
unlink("tmp.ps",@tmpFiles);

# Output fit measure if required.
if ( $showFit == 1 ) {
    $fitData{'name'            } = "Norberg et al. (2002) galaxy correlation functions";
    $fitData{'chiSquared'      } = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'        } = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}


exit;
