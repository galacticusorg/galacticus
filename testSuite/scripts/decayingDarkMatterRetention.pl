#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::GSL::INTEG;
use PDL::Constants qw(PI);

# Compute the fraction of particles, and kick energy retained in a halo for decaying dark matter models. Uses a Monte Carlo
# simulation to evaluate these fractions, which is used as a target dataset to test our numerical calculations.
# Andrew Benson (01-May-2024)

# Specify the velocity dispersion and the escape velocity.
my $velocityDispersion = pdl 20.0;
my $velocityEscape     = pdl 60.0;

# Generate a sequence of velocity kicks for which to compute the retained fractions.
my $velocityKicks      = pdl sequence(130)+1.0;

# Generate a sequence of escape velocities for which to compute the retained fractions.
my $velocityEscapes    = pdl sequence( 10)+60.0;

# Specify the number of Monte Carlo trials.
my $countTrials = 30;

# Specify the number of particles to sample for each trial.
my $countSample = 1000000;

# Create arrays to store results.
my $fractionsRetained = pdl zeroes($countTrials,nelem($velocityKicks),nelem($velocityEscapes));
my $energiesRetained  = pdl zeroes($countTrials,nelem($velocityKicks),nelem($velocityEscapes));

# Iterate over escape velocities.
my $velocityWidth = $velocityDispersion->copy();
for(my $j=0;$j<nelem($velocityEscapes);++$j) {
    my $velocityEscape   = $velocityEscapes->(($j));
    print "Velocity escape = ".$velocityEscape." [".($j+1)." of ".nelem($velocityEscapes)."]\n";
    # Compute the width of the truncated Maxwell-Boltzmann distribution (truncated at the escape velocity) required to give a mean
    # squared velocity of 3 times the velocity dispersion. We use a simple iterative scheme here.
    my $velocitySquaredMeanTarget = 3.0*$velocityDispersion**2;
    my $boostFactor               = pdl 1.0;
    my $velocityMeanSquared;
    for(my $i=0;$i<10;++$i) {
	$velocityWidth         .= $boostFactor*$velocityDispersion;
	my $limit               = 1000     ;
	my $key                 =    5     ;
	my $toleranceRelative   =    1.0e-4;
	my $toleranceAbsolute   =    0.0e+0;
	($velocityMeanSquared)  = gslinteg_qag(\&velocityMeanSquaredIntegrand ,0.0,$velocityEscape,$toleranceRelative,$toleranceAbsolute,$limit,$key);
	my ($normalization   )  = gslinteg_qag(\&velocitNormalizationIntegrand,0.0,$velocityEscape,$toleranceRelative,$toleranceAbsolute,$limit,$key);
	$velocityMeanSquared   /= $normalization;
	$boostFactor           *= sqrt($velocitySquaredMeanTarget/$velocityMeanSquared);

    }
    print "\tComputed boost factor of: ".$boostFactor."\n";
    print "\t\tMean squared velocity is ".$velocityMeanSquared." (target was ".$velocitySquaredMeanTarget.")\n";
    # Construct the truncated Maxwell-Boltzmann distribution and find the cumulative distribution function.
    my $countVelocity = 100000;
    my $velocity      = pdl sequence($countVelocity)/($countVelocity-1)*$velocityEscape;
    my $pdf           = $velocity**2*exp(-0.5*($velocity/$velocityWidth)**2);
    my $cdf           = $pdf->cumusumover()/$pdf->sum();
    # Iterate over velocity kicks.
    for(my $i=0;$i<nelem($velocityKicks);++$i) {
	# Iterate over trials.
	my $velocityKick   = $velocityKicks->(($i));
	print "\tVelocity kick = ".$velocityKick." [".($i+1)." of ".nelem($velocityKicks)."]\n";
	for(my $trial=0;$trial<$countTrials;++$trial) {
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
	    $energiesRetained ->(($trial),($i),($j)) .= $energyGain/(0.5*$velocityKick**2);
	    $fractionsRetained->(($trial),($i),($j)) .= nelem($retained)/$countSample;
	}
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
$output->dataset('velocityEscape'             )->set($velocityEscapes             );
$output->dataset('fractionRetained'           )->set($fractionRetained           );
$output->dataset('fractionRetainedUncertainty')->set($fractionRetainedUncertainty);
$output->dataset('energyRetained'             )->set($energyRetained             );
$output->dataset('energyRetainedUncertainty'  )->set($energyRetainedUncertainty  );
$output->attrSet(velocityDispersion => $velocityDispersion);

exit;

sub velocityMeanSquaredIntegrand {
    # Velocity distribution function.
    my ($velocity) = @_;
    return sqrt(2.0/PI)*$velocityWidth*($velocity/$velocityWidth)**4*exp(-0.5*($velocity/$velocityWidth)**2);
}

sub velocitNormalizationIntegrand {
    # Velocity distribution function.
    my ($velocity) = @_;
    return sqrt(2.0/PI)*($velocity/$velocityWidth)**2/$velocityWidth*exp(-0.5*($velocity/$velocityWidth)**2);
}
