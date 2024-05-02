#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Compute the fraction of particles, and kick energy retained in a halo for decaying dark matter models. Uses a Monte Carlo
# simulation to evaluate these fractions, which is used as a target dataset to test our numerical calculations.
# Andrew Benson (01-May-2024)

# Specify the velocity dispersion and the escape velocity.
my $velocityDispersion = pdl 20.0;
my $velocityEscape     = pdl 60.0;

# Specify the width of the truncated Maxwell-Boltzmann distribution (truncated at the escape velocity) required to give an mean
# squared velocity of 3. The numerical factor here was found by solving for this mean squared velocity using Mathematica.
my $velocityWidth      = pdl 1.064752990694562*$velocityDispersion;

# Construct the truncated Maxwell-Boltzmann distribution and find the cumulative distribution function.
my $countVelocity = 100000;
my $velocity      = pdl sequence($countVelocity)/($countVelocity-1)*$velocityEscape;
my $pdf           = $velocity**2*exp(-0.5*($velocity/$velocityWidth)**2);
my $cdf           = $pdf->cumusumover()/$pdf->sum();

# Generate a sequence of velocity kicks for which to compute the retained fractions.
my $velocityKicks = pdl sequence(130)+1.0;

# Specify the number of Monte Carlo trials.
my $countTrials = 30;

# Specify the number of particles to sample for each trial.
my $countSample = 1000000;

# Create arrays to store results.
my $fractionsRetained = pdl zeroes($countTrials,nelem($velocityKicks));
my $energiesRetained  = pdl zeroes($countTrials,nelem($velocityKicks));

# Iterate over velocity kicks.
for(my $trial=0;$trial<$countTrials;++$trial) {
    for(my $i=0;$i<nelem($velocityKicks);++$i) {
	my $velocityKick   = $velocityKicks->(($i));
	# Generate random velocities and kick angles.
	my $x1 = pdl random($countSample);
	my $x2 = pdl random($countSample);
	(my $velocities) = interpolate($x1,$cdf,$velocity);
	my $cosTheta     = 2.0*$x2-1.0;
	# Compute the final velocity.
	my $velocityFinal = sqrt($velocities**2+$velocityKick**2+2.0*$velocities*$velocityKick*$cosTheta);
	# Find those particles that are retained.
	my $retained = which($velocityFinal < $velocityEscape);
	# Compute the energy gain of retain particles.
	my $energyGain = sum(0.5*$velocityFinal->($retained)**2-0.5*$velocities->($retained)**2)/$countSample;
	# Store the retained fractions.
	$energiesRetained ->(($trial),($i)) .= $energyGain/(0.5*$velocityKick**2);
	$fractionsRetained->(($trial),($i)) .= nelem($retained)/$countSample;
    }
}

# Compute the mean fractions and their uncertainties.
my $fractionRetained            = $fractionsRetained->average();
my $fractionRetainedUncertainty = sqrt(($fractionsRetained->pow(2)->average()-$fractionRetained->pow(2))/$countTrials);
my $energyRetained              = $energiesRetained ->average();
my $energyRetainedUncertainty   = sqrt(($energiesRetained ->pow(2)->average()-$energyRetained  ->pow(2))/$countTrials);

# Store results to file.
my $output = new PDL::IO::HDF5(">testSuite/data/decayingDarkMatterRetention.hdf5");
$output->dataset('velocityKick'               )->set($velocityKicks              );
$output->dataset('fractionRetained'           )->set($fractionRetained           );
$output->dataset('fractionRetainedUncertainty')->set($fractionRetainedUncertainty);
$output->dataset('energyRetained'             )->set($energyRetained             );
$output->dataset('energyRetainedUncertainty'  )->set($energyRetainedUncertainty  );
$output->attrSet(velocityDispersion => $velocityDispersion);
$output->attrSet(velocityEscape     => $velocityEscape    );

exit;
