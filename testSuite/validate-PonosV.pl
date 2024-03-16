#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use JSON::PP;
use Git;

# Run models to validate subhalo projected density against the PonosV simulation of Fiacconi et al. (2016; https://ui.adsabs.harvard.edu/abs/2016ApJ...824..144F).
# Andrew Benson (13-February-2024)

# Make output directory.
system("mkdir -p outputs/");

# Run the validate model.
system("cd ..; export OMP_NUM_THREADS=2; ./Galacticus.exe testSuite/parameters/validate_PonosV.xml");
unless ( $? == 0 ) {
   print "FAIL: PonosV validation model failed to run\n";
   exit;
}

# Read data.
my $model       = new PDL::IO::HDF5("outputs/validate_PonosV.hdf5");
my $outputs     = $model  ->group('Outputs'         );
my $redshift0p0 = $outputs->group('Output2/nodeData');
my $redshift0p7 = $outputs->group('Output1/nodeData');
my $data;
$data->{'0.0'}->{$_} = $redshift0p0->dataset($_)->get()
    foreach ( 'nodeIsIsolated', 'darkMatterOnlyRadiusVirial' );
$data->{'0.7'}->{$_} = $redshift0p7->dataset($_)->get()
    foreach ( 'nodeIsIsolated', 'hostDarkMatterOnlyRadiusVirial', 'positionOrbitalX', 'positionOrbitalY', 'positionOrbitalZ', 'satelliteBoundMass', 'mergerTreeIndex' );

# Validate z=0.0 host halo virial radii.
my $hostsFinal         = which($data->{'0.0'}->{'nodeIsIsolated'} == 1);
my $radiusVirialTarget = pdl 0.6005; # Virial radius of PonosV from Table 1 of Fiacconi et al. (2016; https://ui.adsabs.harvard.edu/abs/2016ApJ...824..144F).
my $offsetFractional   = ($data->{'0.0'}->{'darkMatterOnlyRadiusVirial'}->($hostsFinal)-$radiusVirialTarget)/$radiusVirialTarget;
if ( any($offsetFractional > 0.01) ) {
    print "FAIL: PonosV z=0.0 host virial radii\n";
    print "   Expected: ".$radiusVirialTarget."\n";
    print "   Found: ".$data->{'0.0'}->{'darkMatterOnlyRadiusVirial'}->($hostsFinal)."\n";
} else {
    print "SUCCESS: PonosV z=0.0 host virial radii\n";
}

# Select z=0.7 subhalos.
my $radiusFractionalMinimum = pdl 0.00e0;
my $radiusFractionalMaximum = pdl 0.04e0;
my $massBoundMinimum8       = pdl 0.50e8;
my $massBoundMaximum8       = pdl 2.00e8;
my $massBoundMinimum9       = pdl 0.50e9;
my $massBoundMaximum9       = pdl 2.00e9;
my $kilo                    = pdl 1.00e3;
$data->{'0.7'}->{'radiusOrbital2D'} = sqrt(+$data->{'0.7'}->{'positionOrbitalX'}**2+$data->{'0.7'}->{'positionOrbitalY'}**2                                        );
$data->{'0.7'}->{'radiusOrbital3D'} = sqrt(+$data->{'0.7'}->{'positionOrbitalX'}**2+$data->{'0.7'}->{'positionOrbitalY'}**2+$data->{'0.7'}->{'positionOrbitalZ'}**2);
my $selection
    = which
    (
     ($data->{'0.7'}->{'nodeIsIsolated'    }               == 0                                                                          ) # } Select subhalos.
     &              
     ($data->{'0.7'}->{'radiusOrbital3D'   }               <=                          $data->{'0.7'}->{'hostDarkMatterOnlyRadiusVirial'}) # } Select subhalos within the host virial radius.
     &              
     ($data->{'0.7'}->{'radiusOrbital2D'   }               >= $radiusFractionalMinimum*$data->{'0.7'}->{'hostDarkMatterOnlyRadiusVirial'}) # ⎫ 
     &                                                                                                                                     # ⎬ Select subhalos close to projected radius of 0.02 of host virial radius.
     ($data->{'0.7'}->{'radiusOrbital2D'   }               <= $radiusFractionalMaximum*$data->{'0.7'}->{'hostDarkMatterOnlyRadiusVirial'}) # ⎭
    );
my $selection8
    = which
    (
     ($data->{'0.7'}->{'satelliteBoundMass'}->($selection) >= $massBoundMinimum8                                                         ) # ⎫
     &                                                                                                                                     # ⎬ Select subhalos close to a bound mass of 10⁸M☉.
     ($data->{'0.7'}->{'satelliteBoundMass'}->($selection) <= $massBoundMaximum8                                                         ) # ⎭
    );
my $selection9
    = which
    (
     ($data->{'0.7'}->{'satelliteBoundMass'}->($selection) >= $massBoundMinimum9                                                         ) # ⎫
     &                                                                                                                                     # ⎬ Select subhalos close to a bound mass of 10⁹M☉.
     ($data->{'0.7'}->{'satelliteBoundMass'}->($selection) <= $massBoundMaximum9                                                         ) # ⎭
    );

# Compute the number density of selected subhalos in each tree.
(my $treeCount) = $model->group('Parameters')->group('mergerTreeBuildMasses')->attrGet('treeCount');
my $subhaloSurfaceDensity8 = pdl double zeroes($treeCount->((0))->sclr());
my $subhaloSurfaceDensity9 = pdl double zeroes($treeCount->((0))->sclr());
for(my $i=0;$i<$treeCount;++$i) {
    my $selectTree8    = which($data->{'0.7'}->{'mergerTreeIndex'}->($selection)->($selection8) == $i+1);
    my $selectTree9    = which($data->{'0.7'}->{'mergerTreeIndex'}->($selection)->($selection9) == $i+1);
    if ( nelem($selectTree8) > 0 ) {
	my $countSubhalos8 = pdl double($selectTree8->dim(0));
	$subhaloSurfaceDensity8->(($i)) .=
	    +$countSubhalos8
	    /2.0
	    /PI
	    /$data->{'0.7'}->{'hostDarkMatterOnlyRadiusVirial'}->($selection)->($selectTree8)->((0))**2
	    /$kilo                                                                                  **2
	    /(
		+$radiusFractionalMaximum**2
		-$radiusFractionalMinimum**2
	    )
	    /log10(
		   +$massBoundMaximum8
		   /$massBoundMinimum8
	          );
    } else {
	$subhaloSurfaceDensity8->(($i)) .= 0.0;	
    }
    if ( nelem($selectTree9) > 0 ) {
	my $countSubhalos9 = pdl double($selectTree9->dim(0));
	$subhaloSurfaceDensity9->(($i)) .=
	    +$countSubhalos9
	    /2.0
	    /PI
	    /$data->{'0.7'}->{'hostDarkMatterOnlyRadiusVirial'}->($selection)->($selectTree9)->((0))**2
	    /$kilo                                                                                  **2
	    /(
		+$radiusFractionalMaximum**2
		-$radiusFractionalMinimum**2
	    )
	    /log10(
		   +$massBoundMaximum9
		   /$massBoundMinimum9
	          );
    } else {
	$subhaloSurfaceDensity9->(($i)) .= 0.0;	
    }
}

# Set target values from the PonosV simulation of Fiacconi et al. (2016; https://ui.adsabs.harvard.edu/abs/2016ApJ...824..144F).
my $alphaPonosV                  = pdl 0.850;
my $subhaloSurfaceDensity8PonosV = pdl 0.006;

# Compute percentage of realizations above/below the PonosV subhalo surface density, and report.
(my $above, my $below) = which_both($subhaloSurfaceDensity8 > $subhaloSurfaceDensity8PonosV);
my $percentageAbove = 100.0*double(nelem($above))/$treeCount->((0));
my $percentageBelow = 100.0*double(nelem($below))/$treeCount->((0));
my $statusSurfaceDensity = ($percentageAbove > 5.0 || $percentageBelow > 5.0) ? "SUCCESS" : "FAIL";
print $statusSurfaceDensity.": Percentage of realizations above/below the PonosV subhalo surface density: ".sprintf("%5.1f",$percentageAbove)."/".sprintf("%5.1f",$percentageBelow)."\n";

# Compute the mean slope of the subhalo mass function, and report.
# We exclude models for which there are no subhalos present in one of the mass cuts. This probably introduces some bias.
my $nonZero     = which(($subhaloSurfaceDensity8 > 0.0) & ($subhaloSurfaceDensity9 > 0.0));
my $alphas      = -log($subhaloSurfaceDensity8->($nonZero)        /$subhaloSurfaceDensity9->($nonZero)        )
                  /log(sqrt($massBoundMinimum8*$massBoundMaximum8)/sqrt($massBoundMinimum9*$massBoundMaximum9));
(my $aboveSlope, my $belowSlope) = which_both($alphas > $alphaPonosV);
my $percentageAboveSlope = 100.0*double(nelem($aboveSlope))/nelem($alphas);
my $percentageBelowSlope = 100.0*double(nelem($belowSlope))/nelem($alphas);
my $statusSlope = ($percentageAboveSlope > 5.0 || $percentageBelowSlope > 5.0) ? "SUCCESS" : "FAIL";
print $statusSlope.": Percentage of realizations above/below the PonosV subhalo mass function slope: ".sprintf("%5.1f",$percentageAboveSlope)."/".sprintf("%5.1f",$percentageBelowSlope)."\n";

# Interface with git.
my $repo         = Git->repository(Directory => $ENV{'GALACTICUS_EXEC_PATH'});
my $lastRevision = $repo->command_oneline( [ 'rev-list', '--all' ], STDERR => 0 );
(my $authorName  = $repo->command_oneline( [ 'show', '-s', '--format="%an"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
(my $authorEmail = $repo->command_oneline( [ 'show', '-s', '--format="%ae"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
(my $authorDate  = $repo->command_oneline( [ 'show', '-s', '--format="%aD"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
(my $message     = $repo->command_oneline( [ 'show', '-s', '--format="%s"' , $lastRevision ], STDERR => 0 )) =~ s/"//g;

# Generate content for the validation metrics page.
my $output;
$output =
{
    repoUrl       => "https://github.com/galacticusorg/galacticus",
    parameterFile => "testSuite/parameters/validate_PonosV.xml",
    commit        =>
    {
	author =>
	{
	    name => $authorName,
	    email => $authorEmail
	},
	id        => $lastRevision,
	message   => $message,
	timestamp => $authorDate,
	url       => "https://github.com/galacticusorg/galacticus/commit/".$lastRevision
    },
    surfaceDensity =>
    {
	percentageBelow => $percentageBelow             ->sclr(),
	percentageAbove => $percentageAbove             ->sclr(),
	target          => $subhaloSurfaceDensity8PonosV->sclr()
    },
    slope =>
    {
	percentageBelow => $percentageBelowSlope->sclr(),
	percentageAbove => $percentageAboveSlope->sclr(),
	target          => $alphaPonosV         ->sclr()
    }
};
my $json = JSON::PP->new()->pretty()->encode($output);
open(my $reportFile,">","outputs/results_PonosV.json");
print $reportFile "window.PONOSV_DATA = ";
print $reportFile $json;
close($reportFile);

exit;
