#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath  = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use PDL::GSLSF::GAMMA;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Compute mass systematic model coefficients to describe the mass systematic in the SDSS stellar
# mass functon arising from the choice of profile fitting. Based on the results of Bernardi et
# al. (2013; http://adsabs.harvard.edu/abs/2013MNRAS.436..697B).
# Andrew Benson (07-April-2014)

# Specify the various models and mass function fits used by Bernardi et al. (2013).
my @models = 
    (
     {
	 name     => "cmodel",
	 phiStar  => 0.766e-2,
	 mStar    => 0.4103e9,
	 alpha    => 1.764   ,
	 beta     => 0.384   ,
	 phiGamma => 0.557e-2,
	 mGamma   => 4.7802e9,
	 gamma    => 0.053
     },
     {
	 name     => "S\\\\'ersic",
	 phiStar  => 1.040e-2,
	 mStar    => 0.0094e9,
	 alpha    => 1.665   ,
	 beta     => 0.255   ,
	 phiGamma => 0.675e-2,
	 mGamma   => 2.7031e9,
	 gamma    => 0.296
     },
     {
	 name     => "S\\\\'erExp",
	 phiStar  => 0.892e-2,
	 mStar    => 0.0014e9,
	 alpha    => 2.330   ,
	 beta     => 0.239   ,
	 phiGamma => 0.738e-2,
	 mGamma   => 3.2324e9,
	 gamma    => 0.305
     },
     {
	 name     => "S\\\\'ersic (Simard)",
	 phiStar  => 0.820e-2,
	 mStar    => 0.0847e9,
	 alpha    => 1.755   ,
	 beta     => 0.310   ,
	 phiGamma => 0.539e-2,
	 mGamma   => 5.2204e9,
	 gamma    => 0.072
     }
    );

# Initialize suitable ranges of masses.
my $logMassLimited    = pdl sequence(36)*3.5/35.0+9.0;
my $massLimited       = 10.0**$logMassLimited;
my $logMassNormalized = $logMassLimited-10.8;
my $logMass           = pdl sequence(56)*5.5/55.0+8.0;
my $mass              = 10.0**$logMass;

# Declare variables for GnuPlot;
my ($gnuPlot, $plotFileEPS, $plot);
# Open a pipe to GnuPlot.
my $plotFile = "constraints/dataAnalysis/stellarMassFunction_SDSS_z0.07/profileSystematic.pdf";
($plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.4\n";
print $gnuPlot "set key at screen 0.4,0.2\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set xrange [9.0:12.5]\n";
print $gnuPlot "set yrange [-0.75:+0.75]\n";
print $gnuPlot "set xlabel '\$\\log_{10} M_{\\star,{\\tt cmodel}}\$'\n";
print $gnuPlot "set ylabel '\$\\log_{10} M_{\\star,{\\tt cmodel}}-\\log_{10} M_{\\star}\$'\n";

# Compute mass functions for each model.
my $color = -1;
foreach my $model ( @models ) {
    # Evaluate the Gamma function factor.
    my ($gamma,$error) = gsl_sf_gamma($model->{'alpha'}/$model->{'beta'});
    # Evaluate the mass function.
    $model->{'phi'} = ($model->{'phiStar'}*$model->{'beta'}*($mass/$model->{'mStar'})**$model->{'alpha'}*exp(-($mass/$model->{'mStar'})**$model->{'beta'})/$gamma+$model->{'phiGamma'}*($mass/$model->{'mGamma'})**$model->{'gamma'}*exp(-$mass/$model->{'mGamma'}));
    if ( $model->{'name'} eq "cmodel" ) {
	# For cmodel, evaluate the mass function at a limited range of masses onto which we will interpolate.
	$model->{'phiLimited'} = ($model->{'phiStar'}*$model->{'beta'}*($massLimited/$model->{'mStar'})**$model->{'alpha'}*exp(-($massLimited/$model->{'mStar'})**$model->{'beta'})/$gamma+$model->{'phiGamma'}*($massLimited/$model->{'mGamma'})**$model->{'gamma'}*exp(-$massLimited/$model->{'mGamma'}));
	$model->{'logMassLimited'} = $logMassLimited;
    } else {
	# For other models, compute the corresponding masses by abundance matching with the cmodel mass function.
	$model->{'logMassLimited'} = interpol($models[0]->{'phiLimited'},$model->{'phi'},$logMass);
    }
    # Compute and plot the mass offset.
    $model->{'offset'} = $models[0]->{'logMassLimited'}-$model->{'logMassLimited'};
    unless ( $model->{'name'} eq "cmodel" ) {
	++$color;
	&PrettyPlots::Prepare_Dataset(\$plot,
				      $logMassLimited,
				      $model->{'offset'},
				      style     => "point",
				      symbol    => [6,7], 
				      weight    => [5,3],
				      color     => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'sequence1'}}[$color]},
				      title     => $model->{'name'}
	    );
    }
}
# Construct a simple fit to the offsets.
my @fits = 
    (
     {
	 mu0    => -0.10,
	 kappa0 => -0.00,
	 mu1    => -0.00,
	 kappa1 => -0.33,
	 beta   => +0.50
     },
     {
	 mu0    => -0.10,
	 kappa0 => -0.00,
	 mu1    => -0.00,
	 kappa1 => -0.25,
	 beta   => +0.50
     }
    );
$color = -1;
foreach my $fit ( @fits ) {
    my $systematic =
	+($fit->{'mu0'}+$fit->{'kappa0'}*$logMassNormalized)*(1.0-1.0/(1.0+exp(-$logMassNormalized/$fit->{'beta'})))
	+($fit->{'mu1'}+$fit->{'kappa1'}*$logMassNormalized)*(    1.0/(1.0+exp(-$logMassNormalized/$fit->{'beta'})));
    ++$color;
    &PrettyPlots::Prepare_Dataset(\$plot,
				  $logMassLimited,
				  $systematic,
				  style     => "line",
				  weight    => [5,3],
				  color     => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[$color]},
				  title     => "model"
	);
}
# Finish plotting.
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS,margin => 1);

exit;
