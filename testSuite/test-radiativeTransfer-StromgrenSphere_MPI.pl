#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use Galacticus::Options;

# Test the radiative transfer code by attempting to reproduce a Strömgren sphere solution.
# Andrew Benson (04-December-2019)

# Read in any configuration options.
my $config = &Galacticus::Options::LoadConfig();

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     })
    if ( defined($queueManager) );

# Get any command line options.
my %options =
    (
     'processesPerNode'  => (defined($queueConfig) && exists($queueConfig->{'ppn'})) ? $queueConfig->{'ppn'} : 1,
     'oversubscribe'     => "no",
     'allow-run-as-root' => "no"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# We need at least 2 processes to run this test.
if ( $options{'processesPerNode'} < 2 ) {
    print "SKIPPED: at least 2 processes per node are required for this test\n";
    exit;
}

# Run the calculation.
system("cd ..; mpirun --oversubscribe -np ".$options{'processesPerNode'}.($options{'allow-run-as-root'} eq "yes" ? " --allow-run-as-root" : "").($options{'oversubscribe'} eq "yes" ? " --oversubscribe" : "")." Galacticus.exe testSuite/parameters/test-radiativeTransfer-StromgrenSphere.xml");
die("FAILED: failed to run calculation")
    unless ( $? == 0 );
# Read model output and parameters.
my $outputFile                                                                    = new PDL::IO::HDF5('outputs/radiativeTransfer-StromgrenSphere:MPI0000.hdf5');
my $parameters                                                                    = $outputFile         ->group  ('Parameters'                                                        )       ;
my $recombination                                                                 = $parameters         ->group  ('atomicRecombinationRateRadiative'                                  )       ;
my $computationalDomain                                                           = $parameters         ->group  ('computationalDomain'                                               )       ;
my $model                                                                         = $outputFile         ->group  ('radiativeTransferModel'                                            )       ;
my $densityNumber                                                                 = $model              ->dataset('densityNumberH'                                                    )->get();
my $fractionHydrogenII                                                            = $model              ->dataset('fractionHII'                                                       )->get();
(my $recombinationCoefficient                                                   ) = $recombination      ->attrGet('rateCoefficient'                                                   )       ;
(my $xBoundaries              , my $yBoundaries, my $zBoundaries, my $countCells) = $computationalDomain->attrGet('xBoundaries'              ,'yBoundaries','zBoundaries','countCells')       ;
(my $rateLymanContinuumEmitted                                                  ) = $model              ->attrGet('rateLymanContinuumEmitted'                                         )       ;
# Assert that the grid is the same in each direction.
if ( 
    $countCells ->((1)) != $countCells ->((0)) || $countCells ->((2)) != $countCells ->((0))
    ||
    $yBoundaries->((0)) != $xBoundaries->((0)) || $zBoundaries->((0)) != $xBoundaries->((0))
    ||
    $yBoundaries->((1)) != $xBoundaries->((1)) || $zBoundaries->((1)) != $xBoundaries->((1))
    ) {
    print "FAILED: grid is not the same along each axis\n";
    exit 0;
}
# Compute computational domain cell volumes.
my $megaParsec = pdl 3.086e+22;
my $centi      = pdl 1.000e-02;
my $cellVolume = 
    (
     +($xBoundaries->((1))-$xBoundaries->((0)))
     *($yBoundaries->((1))-$yBoundaries->((0)))
     *($zBoundaries->((1))-$zBoundaries->((0)))
    )
    /$countCells->prodover()
    *($megaParsec/$centi)**3;
my $cellSize  = $xBoundaries->((1))-$xBoundaries->((0));
# Compute the total recombination rate.
my $recombinationRate = sum(($densityNumber*$fractionHydrogenII)**2)*$recombinationCoefficient*$cellVolume;
# Compute the Strömgren radius.
my $radiusStromgren = (3.0*$rateLymanContinuumEmitted/4.0/PI/$densityNumber->max()**2/$recombinationCoefficient)**(1.0/3.0)*$centi/$megaParsec;
# Compute boundary layer volume relative to Strömgren sphere volume.
my $boundaryVolumeFraction = 3.0*$cellSize/$radiusStromgren;
# Test for success in the recombination rate. The tolerance is twice the boundary layer volume fraction as this should be
# approximately the level of uncertainty we expect in the integrated recombination rate due to the finite resolution of the grid.
my $successRate = abs($recombinationRate-$rateLymanContinuumEmitted) < 2.0*$boundaryVolumeFraction*$rateLymanContinuumEmitted ? "success" : "FAIL";
print $successRate.": ionization balance: recombination / ionization rate = ".$recombinationRate." / ".$rateLymanContinuumEmitted." s⁻¹\n";
# Check each cell ionization state.
my $insideIsIonized  = 1;
my $outsideIsNeutral = 1;
for(my $i=0;$i<$countCells->((0));++$i) {
    my $xLower   = $xBoundaries->((0))+($xBoundaries->((1))-$xBoundaries->((0)))/$countCells->((0))*($i+0);
    my $xUpper   = $xBoundaries->((0))+($xBoundaries->((1))-$xBoundaries->((0)))/$countCells->((0))*($i+1);
    my $xMinimum = abs($xLower) <= abs($xUpper) ? abs($xLower) : abs($xUpper);
    my $xMaximum = abs($xLower) >  abs($xUpper) ? abs($xLower) : abs($xUpper);
    for(my $j=0;$j<$countCells->((1));++$j) {
	my $yLower   = $yBoundaries->((0))+($yBoundaries->((1))-$yBoundaries->((0)))/$countCells->((1))*($j+0);
	my $yUpper   = $yBoundaries->((0))+($yBoundaries->((1))-$yBoundaries->((0)))/$countCells->((1))*($j+1);
	my $yMinimum = abs($yLower) <= abs($yUpper) ? abs($yLower) : abs($yUpper);
	my $yMaximum = abs($yLower) >  abs($yUpper) ? abs($yLower) : abs($yUpper);
	for(my $k=0;$k<$countCells->((2));++$k) {
	    my $zLower   = $zBoundaries->((0))+($zBoundaries->((1))-$zBoundaries->((0)))/$countCells->((2))*($k+0);
	    my $zUpper   = $zBoundaries->((0))+($zBoundaries->((1))-$zBoundaries->((0)))/$countCells->((2))*($k+1);
	    my $zMinimum = abs($zLower) <= abs($zUpper) ? abs($zLower) : abs($zUpper);
	    my $zMaximum = abs($zLower) >  abs($zUpper) ? abs($zLower) : abs($zUpper);
	    my $radiusMinimum = sqrt($xMinimum**2+$yMinimum**2+$zMinimum**2);
	    my $radiusMaximum = sqrt($xMaximum**2+$yMaximum**2+$zMaximum**2);
	    if      ( $radiusMaximum < $radiusStromgren-$cellSize ) {
		$insideIsIonized  = 0
		    unless ( $fractionHydrogenII->(($i),($j),($k)) > 0.99 );
	    } elsif ( $radiusMinimum > $radiusStromgren+$cellSize ) {
		$outsideIsNeutral = 0
		    unless ( $fractionHydrogenII->(($i),($j),($k)) < 0.02 );
	    }
	}
    }
}
my $insideIsIonizedStatus  = $insideIsIonized  ? "success" : "FAIL";
my $outsideIsNeutralStatus = $outsideIsNeutral ? "success" : "FAIL";
print $insideIsIonizedStatus .": inside Strömgren radius is ionized\n";
print $outsideIsNeutralStatus.": outside Strömgren radius is neutral\n";
exit 0;
