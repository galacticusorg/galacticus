#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath  = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
require Galacticus::Options;
require Galacticus::Constraints::Covariances;

# Compute likelihood (and make a plot) for a Galacticus model given the projected correlation function data from Hearin et
# al. (2013; http://adsabs.harvard.edu/abs/2013arXiv1310.6747H).

# Get name of input and output files.
die("projectedCorrelationFunction_SDSS_z0.07.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $self               = $0;
my $galacticusFileName = $ARGV[0];
# Create a hash of named arguments.
my $iArg = -1;
my %arguments =
    (
     quiet => 0
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Object to hold all data.
my $data;

# Load observational data.
my $observations = new PDL::IO::HDF5("data/observations/correlationFunctions/Projected_Correlation_Functions_Hearin_2013.hdf5");
my $massMinimum = $observations->dataset('massMinimum')->get();
my $separation  = $observations->dataset('separation' )->get();
my $correlationFunction = $observations->dataset('projectedCorrelationFunctionObserved')->get();
my $covarianceMatrix    = $observations->dataset('covariance')->get();
for(my $entry=0;$entry<3;++$entry) {
    $data->{'observed'}->{$entry}->{'massMinimum'        } = log10($massMinimum->(($entry)));
    $data->{'observed'}->{$entry}->{'separation'         } = $separation;
    $data->{'observed'}->{$entry}->{'correlationFunction'} = $correlationFunction->(:,($entry));
    $data->{'observed'}->{$entry}->{'error'              } = $covarianceMatrix->diagonal(0,1)->($entry*nelem($separation):($entry+1)*nelem($separation)-1)->sqrt();
}
$data->{'observed'}->{'combined'}->{'correlationFunctionCovariance'} = $covarianceMatrix;

# Load model data.
my $galacticusFile   = new PDL::IO::HDF5($galacticusFileName);
my $correlationGroup = $galacticusFile->group('analysis')->group('sdssClusteringZ0.07');
foreach ( "separation", "correlationFunction", "correlationFunctionCovariance" ) {
    $data->{'model'}->{$_} = $correlationGroup->dataset($_)->get();
}

# Create a plot of the correlation function.
if ( exists($arguments{'plotFile'}) ) {
    require GnuPlot::PrettyPlots;
    require GnuPlot::LaTeX;
    require XMP::MetaData;
    # Iterate over masses.
    for(my $entry=0;$entry<3;++$entry) {
	# Declare variables for GnuPlot;
	my ($gnuPlot, $plotFileEPS, $plot);
	# Open a pipe to GnuPlot.
	(my $mass = $data->{'observed'}->{$entry}->{'massMinimum'}) =~ s/\./p/;
	my $suffix = "_lgM".$mass.".eps";
	($plotFileEPS = $arguments{'plotFile'}) =~ s/\.pdf$/$suffix/;
	open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$plotFileEPS."'\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.4,0.2\n";
	print $gnuPlot "set key left\n";
	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set logscale xy\n";
	print $gnuPlot "set mxtics 10\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format x '\$10^{\%L}\$'\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [0.15:25.0]\n";
	print $gnuPlot "set yrange [1.0:1000.0]\n";
	print $gnuPlot "set title 'Projected correlation function at \$z\\approx 0.07\$ for \$\\log_{10}(M_\\star/M_\\odot) > ".$data->{'observed'}->{$entry}->{'massMinimum'}."\$'\n";
	print $gnuPlot "set xlabel 'Separation; \$r_{\\rm p}\$ [Mpc]'\n";
	print $gnuPlot "set ylabel 'Projected correlation; \$w_{\\rm p}(r_{\\rm p})\$ [Mpc]'\n";
	&PrettyPlots::Prepare_Dataset(\$plot,
				      $data->{'observed'}->{$entry}->{'separation'         },
				      $data->{'observed'}->{$entry}->{'correlationFunction'},
				      errorUp   => $data->{'observed'}->{$entry}->{'error'},
				      errorDown => $data->{'observed'}->{$entry}->{'error'},
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
				      title     => "Hearin et al. (2013)"
	    );
	my $separationCount = nelem($data->{'model'}->{'separation'});
	my $modelError = 
	    sqrt(
		$data
		->{'model'}
		->{'correlationFunctionCovariance'}
		->diagonal(0,1)
		->($entry*$separationCount:($entry+1)*$separationCount-1)
	    );
	&PrettyPlots::Prepare_Dataset(\$plot,
				      $data->{'model'}->{'separation'         },
				      $data->{'model'}->{'correlationFunction'}->(:,($entry)),
				      errorUp   => $modelError,
				      errorDown => $modelError,
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{'redYellow'},
				      title     => "Galacticus"
	    );
	&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);
    }
}

# Construct combined datasets.
my $separationCount = nelem(	$data->{'observed'}->{'0'}->{'separation'});
$data->{'observed'}->{'combined'}->{'separation'                   } = pdl zeroes(3*$separationCount                   );
$data->{'observed'}->{'combined'}->{'correlationFunction'          } = pdl zeroes(3*$separationCount                   );
$data->{'observed'}->{'combined'}->{'correlationFunctionError'     } = pdl zeroes(3*$separationCount                   );
for(my $entry=0;$entry<3;++$entry) {
    $data->{'observed'}->{'combined'}->{'separation'              }->($entry*$separationCount:($entry+1)*$separationCount-1)
	.= $data->{'observed'}->{$entry}->{'separation'         };
    $data->{'observed'}->{'combined'}->{'correlationFunction'     }->($entry*$separationCount:($entry+1)*$separationCount-1)
	.= $data->{'observed'}->{$entry}->{'correlationFunction'};
    $data->{'observed'}->{'combined'}->{'correlationFunctionError'}->($entry*$separationCount:($entry+1)*$separationCount-1)
	.= $data->{'observed'}->{$entry}->{'error'};
}
$data->{'model'}->{'combined'}->{'correlationFunction'          } =      $data->{'model'}->{'correlationFunction'          }->flat();
$data->{'model'}->{'combined'}->{'correlationFunctionCovariance'} =      $data->{'model'}->{'correlationFunctionCovariance'};
$data->{'model'}->{'combined'}->{'correlationFunctionError'     } = sqrt($data->{'model'}->{'correlationFunctionCovariance'}->diagonal(0,1));

# Output the results to file if requested.
if ( exists($arguments{'resultFile'}) ) {
    my $resultsFile = new PDL::IO::HDF5(">".$arguments{'resultFile'});
    $resultsFile->dataset('x'             )->set($data->{'observed'}->{'combined'}->{'separation'                   });
    $resultsFile->dataset('y'             )->set($data->{'model'   }->{'combined'}->{'correlationFunction'          });
    $resultsFile->dataset('error'         )->set($data->{'model'   }->{'combined'}->{'correlationFunctionError'     });
    $resultsFile->dataset('covariance'    )->set($data->{'model'   }->{'combined'}->{'correlationFunctionCovariance'});
    $resultsFile->dataset('yData'         )->set($data->{'observed'}->{'combined'}->{'correlationFunction'          });
    $resultsFile->dataset('covarianceData')->set($data->{'observed'}->{'combined'}->{'correlationFunctionCovariance'});
}

# Compute the likelihood:
if ( exists($arguments{'outputFile'}) ) {
    # Construct the full covariance matrix, which is the covariance matrix of the observations plus that of the model.
    my $fullCovariance =
	$data->{'model'   }->{'combined'}->{'correlationFunctionCovariance'}+
	$data->{'observed'}->{'combined'}->{'correlationFunctionCovariance'};
    # Compute the likelihood.
    my $constraint;
    my $logDeterminant;
    my $offsets;
    my $inverseCovariance;
    my $logLikelihood = 
	&Covariances::ComputeLikelihood(
	$data->{'model'   }->{'combined'}->{'correlationFunction'},
	$data->{'observed'}->{'combined'}->{'correlationFunction'},
	$fullCovariance                         ,
	determinant       => \$logDeterminant   ,
	inverseCovariance => \$inverseCovariance,
	offsets           => \$offsets          ,
	quiet             => 0
	);
    $constraint->{'logLikelihood'} = $logLikelihood;
    # Find the Jacobian of the log-likelihood with respect to the model mass function.
    my $jacobian = pdl zeroes(nelem($data->{'model'}->{'combined'}->{'correlationFunction'}),1);
    for(my $i=0;$i<nelem($data->{'model'}->{'combined'}->{'correlationFunction'});++$i) {
	$jacobian->(($i),(0)) .= sum($inverseCovariance->(($i),:)*$offsets);
    }
    # Compute the variance in the log-likelihood due to errors in the model.
    my $logLikelihoodVariance = $jacobian x $data->{'model'}->{'combined'}->{'correlationFunctionCovariance'} x transpose($jacobian);
    $constraint->{'logLikelihoodVariance'} = $logLikelihoodVariance->sclr();
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$arguments{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}

exit;
