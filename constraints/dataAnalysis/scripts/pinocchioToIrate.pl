#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::IO::Misc;

# Convert Pinocchio output catalogs to IRATE format.
# Andrew Benson (28-August-2014)

# Get arguments.
die("pinocchioToIrate.pl <pinocchioDirectoryName> <pinocchioRealization> <irateFileName>")
    unless ( scalar(@ARGV) == 3 );
my $simulationDirectoryName = $ARGV[0];
my $realization             = $ARGV[1];
my $irateFileName           = $ARGV[2];
$simulationDirectoryName .= "/"
    unless ( $simulationDirectoryName =~ m/\/$/ );
# Read the parameter file.
my %parameters;
open(my $parameterFile,$simulationDirectoryName."parametersPinocchio".$realization.".txt");
while ( my $line = <$parameterFile> ) {
    if ( $line =~ m/^([^\#\s]+)\s+([^\s]+)/ ) {
	my $parameterName  = $1;
	my $parameterValue = $2;
	$parameters{$parameterName} = $parameterValue;
    }
}
close($parameterFile);
# Read the list of outputs.
my $outputs = pdl [];
open(my $outputFile,$simulationDirectoryName.$parameters{'OutputList'});
while ( my $line = <$outputFile> ) {
    chomp($line);
    $line =~ s/^\s*//;
    $line =~ s/\s*$//;
    $outputs = $outputs->append($line);
}
close($outputFile);
$outputs = $outputs->qsort();
# Open IRATE file.
my $irateFile = new PDL::IO::HDF5(">".$irateFileName);
# Iterate over outputs.
my $snapshot = 0;
for(my $i=0;$i<nelem($outputs);++$i) {
    # Read catalog.
    ++$snapshot;
    my $redshiftLabel   = sprintf("%6.4f",$outputs->(($i))->sclr());
    my $catalogFileName = $simulationDirectoryName."pinocchio.".$redshiftLabel.".".$parameters{'RunFlag'}.".catalog.out";
    print "Reading ".$catalogFileName."\n";
    (my $mass, my $x, my $y, my $z, my $vx, my $vy, my $vz) = rcols($catalogFileName,1,5,6,7,8,9,10);
    # Construct 3D datasets.
    my $center   = pdl zeroes(3,nelem($x));
    my $velocity = pdl zeroes(3,nelem($x));
    $center  ->((0),:) .= $x;
    $center  ->((1),:) .= $y;
    $center  ->((2),:) .= $z;
    $velocity->((0),:) .= $vx;
    $velocity->((1),:) .= $vy;
    $velocity->((2),:) .= $vz;
    # Snapshot group.
    my $snapshotLabel = sprintf("Snapshot%5.5i",$snapshot);
    my $snapshot = $irateFile->group($snapshotLabel);
    $snapshot->attrSet(Redshift => pdl $outputs->(($i))->sclr());
    # Halo catalog.
    my $haloCatalog = $snapshot->group('HaloCatalog');
    $haloCatalog->dataset('Center'  )->set($center  );
    $haloCatalog->dataset('Velocity')->set($velocity);
    $haloCatalog->dataset('Mass'    )->set($mass    );
    $haloCatalog->dataset('Center'  )->attrSet(unitname => "Mpc"   , unitscgs => pdl [3.08568e+24,  0, -1]);
    $haloCatalog->dataset('Velocity')->attrSet(unitname => "km/s"  , unitscgs => pdl [1.00000e+05,  0,  0]);
    $haloCatalog->dataset('Mass'    )->attrSet(unitname => "Msolar", unitscgs => pdl [1.98892e+33,  0,  0]);
}
# Cosmology.
my $cosmology = $irateFile->group('Cosmology');
$cosmology->attrSet(HubbleParam        => pdl $parameters{'Hubble100'    });
$cosmology->attrSet(OmegaBaryon        => pdl $parameters{'OmegaBaryon'  }); 
$cosmology->attrSet(OmegaLambda        => pdl $parameters{'OmegaLambda'  });   
$cosmology->attrSet(OmegaMatter        => pdl $parameters{'Omega0'       });    
$cosmology->attrSet(PowerSpectrumIndex => pdl $parameters{'PowerSpectrum'});     
$cosmology->attrSet(sigma_8            => pdl $parameters{'Sigma8'       });
# Simulation properties.
my $simulation = $irateFile->group('SimulationProperties');
$simulation->attrSet(
    boxSize => pdl $parameters{'BoxSize'}
    );

exit;
