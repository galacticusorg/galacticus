#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Galacticus::Path;
use PDL;
use XML::Simple;
use Astro::Cosmology;
use Math::SigFigs;
use Stats::Means;
use Galacticus::HDF5;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use XMP::MetaData;
use Galacticus::Options;

# Get name of input and output files.
die("starFormationHistory.pl <galacticusFile> <plotFile>")
    unless ( scalar(@ARGV) >= 2 );
my $self          = $0;
my $modelFileName = $ARGV[0];
my $plotFileName  = $ARGV[1];
my %options;
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Get parameters for the Galacticus model.
my $model;
$model->{'file'} = $modelFileName;
&Galacticus::HDF5::Get_Parameters($model);

# Extract global data
&Galacticus::HDF5::Get_History($model,['historyExpansion','historyStarFormationRate']);
my $history        = $model->{'history'};
my $time           =     $history->{'historyTime'             }      ;
my $redshift       = 1.0/$history->{'historyExpansion'        }-1.0e0;
my $SFR            =     $history->{'historyStarFormationRate'}/1.0e9;
my $stellarDensity =     $history->{'historyStellarDensity'   }      ;

# Determine IMF correction factor.
my $imfCorrection;
if ( $model->{'parameters'}->{'imfSelectionFixed'}->{'value'} eq "Chabrier" ) {
    $imfCorrection = 0.6;
} else {
    die('unknown IMF');
}

# Define redshift bins.
my $redshiftPoints = pdl 16;
my $redshiftMin    = pdl 0.0;
my $redshiftMax    = pdl 8.0;
my $redshiftBin    = pdl ($redshiftMax-$redshiftMin)/$redshiftPoints;
my $redshiftBins   = pdl (0..$redshiftPoints-1)*$redshiftBin+$redshiftMin+0.5*$redshiftBin;

# Read the XML data file.
my $xml      = new XML::Simple;
my $data     = $xml->XMLin(&galacticusPath()."data/observations/starFormationRate/Star_Formation_Rate_Data.xml");
my @dataSets;
foreach my $dataSet ( @{$data->{'starFormationRate'}} ) {
    # Extract data.
    my $columns     = $dataSet->{'columns'};
    my $x           = pdl @{$columns->{'redshift'         }->{'data'}};
    my $xLowerError = pdl @{$columns->{'redshiftErrorDown'}->{'data'}};
    my $xUpperError = pdl @{$columns->{'redshiftErrorUp'  }->{'data'}};
    $xLowerError    = $x-$xLowerError;
    $xUpperError    = $x+$xUpperError;
    my $y           = pdl @{$columns->{'sfr'              }->{'data'}};
    my $yLowerError = pdl @{$columns->{'sfrErrorDown'     }->{'data'}};
    my $yUpperError = pdl @{$columns->{'sfrErrorUp'       }->{'data'}};
    $yUpperError    = (10.0**($y+$yUpperError));
    $yLowerError    = (10.0**($y-$yLowerError));
    $y              = (10.0** $y              );
    # Compute cosmology corrections.
    my $cosmologyData = Astro::Cosmology->new(
	omega_matter => $columns->{'sfr'}->{'omega' },
	omega_lambda => $columns->{'sfr'}->{'lambda'},
	H0           => $columns->{'sfr'}->{'hubble'}
	);
    my $cosmologyGalacticus = Astro::Cosmology->new(
	omega_matter => $model->{'parameters'}->{'cosmologyParametersMethod'}->{'OmegaMatter'    }->{'value'},
	omega_lambda => $model->{'parameters'}->{'cosmologyParametersMethod'}->{'OmegaDarkEnergy'}->{'value'},
	H0           => $model->{'parameters'}->{'cosmologyParametersMethod'}->{'HubbleConstant' }->{'value'}
	);
    my $volumeElementData            = $cosmologyData      ->differential_comoving_volume($x);
    my $volumeElementGalacticus      = $cosmologyGalacticus->differential_comoving_volume($x);
    my $luminosityDistanceData       = $cosmologyData      ->luminosity_distance         ($x);
    my $luminosityDistanceGalacticus = $cosmologyGalacticus->luminosity_distance         ($x);
    my $cosmologyCorrection          = 
	+(     $volumeElementData/     $volumeElementGalacticus)
	/($luminosityDistanceData/$luminosityDistanceGalacticus)**2;
    # Apply cosmology corrections.
    $y           = $y          *$cosmologyCorrection;
    $yUpperError = $yUpperError*$cosmologyCorrection;
    $yLowerError = $yLowerError*$cosmologyCorrection;
    # Apply IMF correction.
    $y           = $y          *$imfCorrection;
    $yUpperError = $yUpperError*$imfCorrection;
    $yLowerError = $yLowerError*$imfCorrection;
    # Store the dataset.
    push(
	@dataSets,
	{
	    x           =>               $x,
	    xLowerError => -$xLowerError+$x,
	    xUpperError => +$xUpperError-$x,
	    y           =>               $y,
	    yLowerError => -$yLowerError+$y,
	    yUpperError => +$yUpperError-$y,
	    label       => $dataSet->{'label'}
	}
	);
}

# Make plot of redshift evolution.
my $plot;
my $gnuPlot;
(my $plotFileTeX = $plotFileName) =~ s/\.pdf$/.tex/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
print $gnuPlot "set output '".$plotFileTeX."'\n";
print $gnuPlot "set title offset 0,-1 'Star Formation Rate History'\n";
print $gnuPlot "set xlabel '\$z\$'\n";
print $gnuPlot "set ylabel '\$\\dot{\\rho}(z)\\,\\,[\\mathrm{M}_\\odot/\\hbox{yr}/\\hbox{Mpc}^{-3}]\$'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.325,0.16\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale y\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [-0.25:9.0]\n";
print $gnuPlot "set yrange [0.005:1.0]\n";
print $gnuPlot "set pointsize 2.0\n";
for(my $i=0;$i<scalar(@dataSets);++$i) {
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	 \$plot,
	               $dataSets[$i]->{'x'          },
	               $dataSets[$i]->{'y'          },
	 errorRight => $dataSets[$i]->{'xUpperError'},
	 errorLeft  => $dataSets[$i]->{'xLowerError'},
	 errorUp    => $dataSets[$i]->{'yUpperError'},
	 errorDown  => $dataSets[$i]->{'yLowerError'},
	 title      => $dataSets[$i]->{'label'      },
	 style      => "point",
	 symbol     => [6,7],
	 weight     => [5,3],
	 pointSize  => 0.5,
	 color      => $GnuPlot::PrettyPlots::colorPairs{($i == 0 ? 'cornflowerBlue' : 'lightSkyBlue')}
	);
}
my $nonZeroPoints = which($SFR > 0.0);
&GnuPlot::PrettyPlots::Prepare_Dataset(
    \$plot,
    $redshift->index($nonZeroPoints),$SFR->index($nonZeroPoints),
    style      => "line",
    weight     => [5,3],
    color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
    title      => 'Galacticus'
    );
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);
&XMP::MetaData::Write($plotFileName,$modelFileName,$self);

exit;
