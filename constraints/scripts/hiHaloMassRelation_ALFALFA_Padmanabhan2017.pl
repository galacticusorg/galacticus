#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::Math;
use PDL::IO::HDF5;
use PDL::LinearAlgebra;
use Galacticus::Options;
use Galacticus::HDF5;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use PDL::GSL::INTERP;

# Compute likelihood (and make a plot) for a Galacticus model given the HI mass-halo mass relation of Padmanabhan et
# al. (2017).
# Andrew Benson (05-April-2017)

# Get name of input and output files.
die("hiHaloMassRelation_ALFALFA_Padmanabhan2017.pl <galacticusFile> [options]")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];
# Create a hash of named options.
my $iArg = -1;
my %options =
    (
     quiet => 0
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);


# Define parameters and errors in the measured relation fitting function (from Padmanabhan et al. 2017).
my $N10      = pdl 9.89e-03;
my $N10error = pdl 4.89e-03;
my $M10      = pdl 4.58e+11;
my $M10error = pdl 0.19e+11;
my $b10      = pdl 0.90e+00;
my $b10error = pdl 0.39e+00;
my $y10      = pdl 0.74e+00;
my $y10error = pdl 0.03e+00;
# Construct measured mass relation.
my $dataTableCount                = 1000;
my $massHaloMeanDataLogarithmic10 = pdl sequence($dataTableCount)*5.0/($dataTableCount-1)+10.0;
my $massHaloMeanData              = +10.0**$massHaloMeanDataLogarithmic10;
my $x                             = +$massHaloMeanData/$M10;
my $massHIData                    = +2.0*$N10*$massHaloMeanData/($x**(-$b10)+$x**($y10));
my $jacobianN10                   = +$massHIData/$N10;
my $jacobianM10                   = -$massHIData*($b10*$x**(-$b10+1.0)-$y10*$x**($y10-1.0))/($x**(-$b10)+$x**($y10))**2/$massHaloMeanData;
my $jacobianb10                   = +$massHIData*log($x)/$x**$b10/($x**(-$b10)+$x**$y10);
my $jacobiany10                   = -$massHIData*log($x)*$x**$y10/($x**(-$b10)+$x**$y10);
my $massHIErrorData               = +sqrt(
                                          +($jacobianN10*$N10error)**2
                                          +($jacobianM10*$M10error)**2
                                          +($jacobianb10*$b10error)**2
                                          +($jacobiany10*$y10error)**2
                                         );
my $massHaloMeanDataLogarithmic   = log($massHaloMeanData);
my $massHIDataLogarithmic         = log($massHIData      );
my $massHIErrorDataLogarithmic    = $massHIErrorData/$massHIData;
my $spline                        = PDL::GSL::INTERP->init('cspline',$massHaloMeanDataLogarithmic,$massHIDataLogarithmic);
my $splineError                   = PDL::GSL::INTERP->init('cspline',$massHaloMeanDataLogarithmic,$massHIErrorDataLogarithmic);
# Read model data.
my $model                  = new PDL::IO::HDF5($galacticusFileName);
my $analysis              = $model   ->group  ('analyses'                         )
                                     ->group  ('hiHaloMassRelationPadmanabhan2017')       ;
my $massHaloModel         = $analysis->dataset('massHalo'                         )->get();
my $massHIModel           = $analysis->dataset('massHI'                           )->get();
my $massHICovarianceModel = $analysis->dataset('massHICovariance'                 )->get();
my $modelEntries          = 
    which(
	($massHIModel                          >  0.0                      ) 
	&
	($massHICovarianceModel->diagonal(0,1) >  0.0                      )
	&
	($massHaloModel                        >= $massHaloMeanData->(( 0)))
        &
        ($massHaloModel                        <= $massHaloMeanData->((-1)))
    );
my $massHaloModelLogarithmic    = log($massHaloModel   ->($modelEntries));
my $massHIModelLogarithmic      = log($massHIModel->($modelEntries));
my $massHIErrorModelLogarithmic = $massHICovarianceModel->diagonal(0,1)->($modelEntries)->sqrt()/$massHIModel->($modelEntries);
# Interpolate observational data to model points.
my $massHIDataLogarithmicInterpolated      = $spline     ->eval($massHaloModelLogarithmic);
my $massHIErrorDataLogarithmicInterpolated = $splineError->eval($massHaloModelLogarithmic);
# Output the results to file if requested.
if ( exists($options{'resultFile'}) ) {
    my $resultsFile = new PDL::IO::HDF5(">".$options{'resultFile'});
    $resultsFile->dataset('x'             )->set($massHaloModelLogarithmic                   );
    $resultsFile->dataset('y'             )->set($massHIModelLogarithmic                );
    $resultsFile->dataset('yData'         )->set($massHIErrorDataLogarithmicInterpolated);
}
# Compute the likelihood:
if ( exists($options{'outputFile'}) ) {
    my $constraint;
    if ( nelem($massHIModelLogarithmic) > 0 ) {
	my $logLikelihood = -0.5*($massHIModelLogarithmic-$massHIDataLogarithmicInterpolated)**2/($massHIErrorModelLogarithmic**2+$massHIErrorDataLogarithmicInterpolated**2);
	$constraint->{'logLikelihood'} = $logLikelihood->sclr();
    } else {
	$constraint->{'logLikelihood'} = -1.0e30;
    }
    $constraint->{'logLikelihoodVariance'} = 0.0;
    $constraint->{'label'                } = "alfalfaHIHaloMassRelationZ0.00";
    # Output the constraint.
    my $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"constraint");
    open(oHndl,">".$options{'outputFile'});
    print oHndl $xmlOutput->XMLout($constraint);
    close(oHndl);
}    
# Make a plot if requested.
if ( exists($options{'plotFile'}) ) {
    # Declare variables for GnuPlot;
    my ($gnuPlot, $plotFileTeX, $plot);
    # Open a pipe to GnuPlot.
    ($plotFileTeX = $options{'plotFile'}) =~ s/\.pdf$/.tex/;
    open($gnuPlot,"|gnuplot ");
    print $gnuPlot "set terminal cairolatex pdf standalone color lw 2\n";
    print $gnuPlot "set output '".$plotFileTeX."'\n";
    print $gnuPlot "set lmargin screen 0.15\n";
    print $gnuPlot "set rmargin screen 0.95\n";
    print $gnuPlot "set bmargin screen 0.15\n";
    print $gnuPlot "set tmargin screen 0.95\n";
    print $gnuPlot "set key spacing 1.2\n";
    print $gnuPlot "set key at screen 0.45,0.2\n";
    print $gnuPlot "set key left\n";
    print $gnuPlot "set key bottom\n";
    print $gnuPlot "set logscale xy\n";
    print $gnuPlot "set mxtics 10\n";
    print $gnuPlot "set mytics 10\n";
    print $gnuPlot "set format x '\$10^{\%L}\$'\n";
    print $gnuPlot "set format y '\$10^{\%L}\$'\n";
    print $gnuPlot "set xrange [3.0e10:1.0e15]\n";
    print $gnuPlot "set yrange [1.0e07:1.0e11]\n";
    print $gnuPlot "set title offset 0,-1 'HI mass-halo mass relation at \$z=0\$'\n";
    print $gnuPlot "set xlabel 'Halo mass; \$M_{\\mathrm{FoF},b=0.2}\\, [\\mathrm{M}_\\odot]\$'\n";
    print $gnuPlot "set ylabel 'HI mass; \$M_\\mathrm{HI}\\, [\\mathrm{M}_\\odot]\$'\n";
    print $gnuPlot "set pointsize 1.0\n";
    # Plot measured relation.
    &GnuPlot::PrettyPlots::Prepare_Dataset
	(
	 \$plot,
	 $massHaloMeanData,
	 $massHIData-$massHIErrorData,
	 y2           => $massHIData+$massHIErrorData,
	 style        => "filledCurve",
	 weight       => [3,1],
	 color        => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	 transparency => 0.75
	);
    &GnuPlot::PrettyPlots::Prepare_Dataset
	(
	 \$plot,
	 $massHaloMeanData,
	 $massHIData,
	 style     => "line",
	 weight    => [3,1],
	 color     => $GnuPlot::PrettyPlots::colorPairs{'cornflowerBlue'},
	 title     => "Padmanabhan \\\\& Refregier (2017)"
	);
    # Plot model points.
    &GnuPlot::PrettyPlots::Prepare_Dataset
	(
	 \$plot,
	 $massHaloModel,
	 $massHIModel,
	 errorDown  => $massHICovarianceModel->diagonal(0,1)->sqrt(),
	 errorUp    => $massHICovarianceModel->diagonal(0,1)->sqrt(),
	 style      => "point",
	 symbol     => [6,7], 
	 weight     => [3,1],
	 pointSize  => 0.5,
	 color      => $GnuPlot::PrettyPlots::colorPairs{'redYellow'},
	 title      => "Galacticus"
	);
    # Finalize plot.
    &GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
    close($gnuPlot);
    &GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX,margin => 1);
}

exit;
