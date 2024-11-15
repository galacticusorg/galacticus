#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::Constants qw(PI);
use JSON::PP;
use Git;
use Stats::Histograms;

# Run models to validate suppression of the halo mass function by baryons.
# Andrew Benson (01-April-2024)

# Define the target data for the mass function suppression. This was read from Figure 2 of Zheng et al. (2024;
# https://ui.adsabs.harvard.edu/abs/2024arXiv240317044Z).
#
# * 'withBaryons' corresponds to the "RI" model of Zheng et al.
# * 'withBaryons_noReionization' corresponds to the "NR" model of Zheng et al.
#
# Array indices correspond to redshift:
#
# * 4 ==> z = 9.27
# * 3 ==> z = 5.72
# * 2 ==> z = 3.06
# * 1 ==> z = 0.00
my @redshifts = ( "", "9.27", "5.72", "3.06", "0.00" );
my $target;
$target->{'withBaryons'               }->[4] =
    pdl [
	0.7115384615384613,
	0.6961538461538460,
	0.7576923076923074,
	0.6961538461538460,
	0.7730769230769230,
	0.7346153846153844,
	0.6423076923076921,
	0.9500000000000000
    ];
$target->{'withBaryons_noReionization'}->[4] =
    pdl [
	0.7346153846153844,
	0.7423076923076921,
	0.8269230769230768,
	0.9038461538461537,
	0.9423076923076922,
	1.0346153846153845,
	0.9038461538461537,
	0.9807692307692306
    ];
$target->{'withBaryons'               }->[3] =
    pdl [
	0.6750000000000000,
	0.7321428571428572,
	0.6821428571428573,
	0.7392857142857143,
	0.5607142857142857,
	0.7464285714285716,
	0.9964285714285716,
	1.0000000000000000
    ];
$target->{'withBaryons_noReionization'}->[3] =
    pdl [
	0.7321428571428572,
	0.7678571428571429,
	0.8250000000000001,
	0.9464285714285715,
	0.8964285714285715,
	0.9892857142857143,
	0.9964285714285716,
	1.0000000000000000
    ];
$target->{'withBaryons'               }->[2] =
    pdl [
	0.6685714285714284,
	0.7028571428571427,
	0.6876190476190476,
	0.7371428571428571,
	0.6609523809523808,
	0.9885714285714284,
	0.9885714285714284,
	1.0000000000000000
    ];
$target->{'withBaryons_noReionization'}->[2] =
    pdl [
	0.6723809523809523,
	0.7104761904761904,
	0.8247619047619046,
	0.9314285714285714,
	0.9466666666666665,
	0.9999999999999999,
	0.9885714285714284,
	1.0000000000000000
    ];
$target->{'withBaryons'               }->[1] =
    pdl [
	0.6756756756756757,
	0.6540540540540540,
	0.7261261261261260,
	1.0000000000000000,
	1.0072072072072071,
	1.0000000000000000,
	1.0072072072072071,
	1.0000000000000000
    ];
$target->{'withBaryons_noReionization'}->[1] =
    pdl [
	0.6756756756756757,
	0.6612612612612612,
	0.7333333333333334,
	1.0000000000000000,
	0.9927927927927928,
	1.0000000000000000,
	1.0000000000000000,
	1.0072072072072071
    ];

# Define χ² targets for each dataset.
my $chiSquaredTarget;
$chiSquaredTarget->{'withBaryons'               }->[1] = 6.0;
$chiSquaredTarget->{'withBaryons_noReionization'}->[1] = 7.0;
$chiSquaredTarget->{'withBaryons'               }->[2] = 4.0;
$chiSquaredTarget->{'withBaryons_noReionization'}->[2] = 3.0;
$chiSquaredTarget->{'withBaryons'               }->[3] = 5.0;
$chiSquaredTarget->{'withBaryons_noReionization'}->[3] = 5.0;
$chiSquaredTarget->{'withBaryons'               }->[4] = 3.0;
$chiSquaredTarget->{'withBaryons_noReionization'}->[4] = 4.0;

# Make output directory.
system("mkdir -p outputs/");

# Run the validate model.
system("cd ..; ./Galacticus.exe testSuite/parameters/validate_baryonicSuppression_IGM_evolution.xml");
unless ( $? == 0 ) {
  print "FAIL: baryonic suppression (IGM evolution) validation model failed to run\n";
  exit;
}

# Read data and repackage into a file suitable for re-reading by other models.
{
    my $model         = new PDL::IO::HDF5( "outputs/validate_baryonicSuppression_IGM_evolution.hdf5");
    my $igmFile       = new PDL::IO::HDF5(">outputs/validate_baryonicSuppression_IGM.hdf5"          );
    my $data;
    my $igmProperties = $model        ->group  ('igmProperties')       ;
    $data->{$_}       = $igmProperties->dataset($_             )->get()
	foreach ( 'redshift', 'temperature', 'densityHydrogen1', 'densityHydrogen2', 'densityHelium1', 'densityHelium2', 'densityHelium3' );
    $data->{'hIonizedFraction' } =  $data->{'densityHydrogen2'}                                                         /($data->{'densityHydrogen1'}+$data->{'densityHydrogen2'}                          );
    $data->{'heIonizedFraction'} = ($data->{'densityHelium2'  }+$data->{'densityHelium3'}                              )/($data->{'densityHelium1'  }+$data->{'densityHelium2'  }+$data->{'densityHelium3'});
    $data->{'electronFraction' } = ($data->{'densityHydrogen2'}+$data->{'densityHelium2'}+2.0*$data->{'densityHelium3'})/($data->{'densityHydrogen1'}+$data->{'densityHydrogen2'}                          );
    $igmFile->dataset('redshift'         )->set($data->{'redshift'         });
    $igmFile->dataset('matterTemperature')->set($data->{'temperature'      });
    $igmFile->dataset('hIonizedFraction' )->set($data->{'hIonizedFraction' });
    $igmFile->dataset('heIonizedFraction')->set($data->{'heIonizedFraction'});
    $igmFile->dataset('electronFraction' )->set($data->{'electronFraction' });
    $igmFile->attrSet(extrapolationAllowed => pdl long 1);
    $igmFile->attrSet(fileFormat           => pdl long 1);
}

# Establish bins in halo mass.
my $massHaloLogarithmicBins = pdl sequence(8)/7.0*3.5+4.0;
my $haloMassFunction;

# Run models with and without baryonic suppression.
foreach my $suffix ( "withoutBaryons", "withBaryons", "withBaryons_noReionization" ) {
    print "Running model: '".$suffix."'\n";
    system("cd ..; ./Galacticus.exe testSuite/parameters/validate_baryonicSuppression_evolve_".$suffix.".xml");
    unless ( $? == 0 ) {
    	print "FAIL: baryonic suppression (evolve: '".$suffix."') validation model failed to run\n";
    	exit;
    }
    my $model                          = new PDL::IO::HDF5("outputs/validate_baryonicSuppression_evolve_".$suffix.".hdf5");
    my $cosmology                      = $model    ->group  ('Parameters/cosmologyParameters');
    (my $OmegaMatter, my $OmegaBaryon) = $cosmology->attrGet('OmegaMatter','OmegaBaryon');
    my $fractionDarkMatter             = ($OmegaMatter-$OmegaBaryon)/$OmegaMatter;
    # Iterate over outputs.
    for(my $outputIndex=1;$outputIndex<=4;++$outputIndex) {
	print "\tProcessing output: ".$outputIndex."\n";
	# Read all required data.
	my $output = $model ->group('Outputs/Output'.$outputIndex);
	my $nodes  = $output->group('nodeData'                   );
	my $data;
	$data->{$_} = $nodes->dataset($_)->get()
	    foreach ( 'nodeIsIsolated', 'mergerTreeIndex', 'nodeIndex', 'parentIndex', 'hotHaloMass', 'diskMassGas', 'massHaloEnclosedCurrent', 'basicMass', 'nodeIsIsolated', 'mergerTreeWeight' );
	# Identify isolated and subhalos.
	(my $isolated, my $subhalo) = which_both($data->{'nodeIsIsolated'} == 1);
	# Build an index mapping each halo to its host halo. We take advantage here of the depth-first ordering of the outputs.
	my $index = pdl zeros(nelem($data->{'nodeIsIsolated'}));
	for(my $i=0;$i<nelem($isolated);++$i) {
	    my $indexStart                   = $i == 0 ? 0: $isolated->(($i-1))+1;
	    my $indexEnd                     =              $isolated->(($i  ))  ;
	    $index->($indexStart:$indexEnd) .=              $isolated->(($i  ))  ;
	}
	# Accumulate masses of isolated halos.
	my $massHalo;
	if ( $suffix eq "withoutBaryons" ) {
	    # In models without baryons, the halo mass is just the dark matter mass.
	    $massHalo            = +$data->{'massHaloEnclosedCurrent'};
	} else {
	    # In models with baryons we assume that the hot gas of the host is distributed as the dark matter, so compute a
	    # correction factor to account for the differing virial density contrast definitions in Galacticus and Zheng et
	    # al. Gas in the disk is assumed to always be included within the virial radius.
	    my $correctionFactor = +$data->{'massHaloEnclosedCurrent'}
	                           /$data->{'basicMass'              };
	    $massHalo            = +$data->{'massHaloEnclosedCurrent'}*$fractionDarkMatter
		                   +$data->{'hotHaloMass'            }*$correctionFactor
		                   +$data->{'diskMassGas'            };
	    # Accumulate masses of subhalos halos. We assume that these are also distributed as the dark matter of the host and so
	    # apply the same correction factor.
	    $massHalo->($index->($subhalo)) += $data->{'hotHaloMass'}->($subhalo)*$correctionFactor->($index->($subhalo));
	    $massHalo->($index->($subhalo)) += $data->{'diskMassGas'}->($subhalo)*$correctionFactor->($index->($subhalo));
	}
	# Construct final quantities needed for the mass function.
	my $weight              = $data    ->{'mergerTreeWeight'}->($isolated)         ;
	my $massHaloLogarithmic = $massHalo                      ->($isolated)->log10();
	# Construct the mass function.
	($haloMassFunction->{$suffix}->[$outputIndex]->{'massFunction'}, $haloMassFunction->{$suffix}->[$outputIndex]->{'massFunctionError'})
	    = &Stats::Histograms::Histogram($massHaloLogarithmicBins,$massHaloLogarithmic,$weight,differential => 1);
    }
}

# Compute ratios of mass functions with the dark matter only model mass function.
my $output;
my $chiSquared;
my $failed = 0;
for(my $outputIndex=1;$outputIndex<=4;++$outputIndex) {
    foreach my $suffix ( "withBaryons", "withBaryons_noReionization" ) {
	$haloMassFunction->{$suffix}->[$outputIndex]->{'ratio'     } = $haloMassFunction->{$suffix}->[$outputIndex]->{'massFunction'}/$haloMassFunction->{'withoutBaryons'}->[$outputIndex]->{'massFunction'};
	$haloMassFunction->{$suffix}->[$outputIndex]->{'ratioError'} =
	    +sqrt(
	          +(
		    +$haloMassFunction->{$suffix         }->[$outputIndex]->{'massFunctionError'}
		    /$haloMassFunction->{$suffix         }->[$outputIndex]->{'massFunction'     }
		   )**2
		  +(
		    +$haloMassFunction->{'withoutBaryons'}->[$outputIndex]->{'massFunctionError'}
		    /$haloMassFunction->{'withoutBaryons'}->[$outputIndex]->{'massFunction'     }
		   )**2
	         )
	    *        $haloMassFunction->{$suffix         }->[$outputIndex]->{'ratio'            };
	$chiSquared->{$suffix}->[$outputIndex] =
	    +sum  (
	           +(
		     +$haloMassFunction->{$suffix}->[$outputIndex]->{'ratio'     }
		     -$target          ->{$suffix}->[$outputIndex]
		    )		                                                  **2
		   /  $haloMassFunction->{$suffix}->[$outputIndex]->{'ratioError'}**2
	          )
	    /nelem(   $target          ->{$suffix}->[$outputIndex]);
	delete($haloMassFunction->{$suffix         }->[$outputIndex]->{'massFunction'     });
	delete($haloMassFunction->{$suffix         }->[$outputIndex]->{'massFunctionError'});
	@{$output->{'model'}->{$suffix         }->[$outputIndex]->{'ratio'     }} = $haloMassFunction->{$suffix         }->[$outputIndex]->{'ratio'     }->list();
	@{$output->{'model'}->{$suffix         }->[$outputIndex]->{'ratioError'}} = $haloMassFunction->{$suffix         }->[$outputIndex]->{'ratioError'}->list();	
	# Report.
	my $status     = $chiSquared->{$suffix}->[$outputIndex] < $chiSquaredTarget->{$suffix}->[$outputIndex] ? "SUCCESS" : "FAILED";
	$failed = 1
	    if ( $status eq "FAILED" );
	my $inequality = $chiSquared->{$suffix}->[$outputIndex] < $chiSquaredTarget->{$suffix}->[$outputIndex] ? "<"       : "≥"     ;
	print $status.": model '".$suffix.(" " x (length("withBaryons_noReionization")-length($suffix)))."' at z=".$redshifts[$outputIndex]." validation (χ² = ".sprintf("%5.3f",$chiSquared->{$suffix}->[$outputIndex])." ".$inequality." ".sprintf("%5.3f",$chiSquaredTarget->{$suffix}->[$outputIndex]).")\n";
    }
    delete($haloMassFunction->{'withoutBaryons'}->[$outputIndex]->{'massFunction'     });
    delete($haloMassFunction->{'withoutBaryons'}->[$outputIndex]->{'massFunctionError'});
}

# Interface with git.
my $repo         = Git->repository(Directory => $ENV{'GALACTICUS_EXEC_PATH'});
my $lastRevision = $repo->command_oneline( [ 'rev-list', '--all' ], STDERR => 0 );
(my $authorName  = $repo->command_oneline( [ 'show', '-s', '--format="%an"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
(my $authorEmail = $repo->command_oneline( [ 'show', '-s', '--format="%ae"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
(my $authorDate  = $repo->command_oneline( [ 'show', '-s', '--format="%aD"', $lastRevision ], STDERR => 0 )) =~ s/"//g;
(my $message     = $repo->command_oneline( [ 'show', '-s', '--format="%s"' , $lastRevision ], STDERR => 0 )) =~ s/"//g;

# Generate content for the validation metrics page.
$output->{'repoUrl'      } = "https://github.com/galacticusorg/galacticus";
$output->{'parameterFile'} = "testSuite/parameters/validate_baryonicSuppression_evolve_withBaryons.xml";
$output->{'commit'       } =
{
    author =>
    {
    	    name  => $authorName,
    	    email => $authorEmail
    },
    id        => $lastRevision,
    message   => $message,
    timestamp => $authorDate,
    url       => "https://github.com/galacticusorg/galacticus/commit/".$lastRevision
};
foreach my $suffix ( "withBaryons", "withBaryons_noReionization" ) {
    for(my $outputIndex=1;$outputIndex<=4;++$outputIndex) {
	@{$output->{'target'}->{$suffix}->[$outputIndex]} = $target->{$suffix}->[$outputIndex]->list();
    }
}
@{$output->{'redshift'}} = @redshifts;
@{$output->{'massHalo'}} = $massHaloLogarithmicBins->list();
my $json = JSON::PP->new()->pretty()->encode($output);
open(my $reportFile,">","outputs/results_baryonicSuppression.json");
print $reportFile "window.BARYONICSUPPRESSION_DATA = ";
print $reportFile $json;
close($reportFile);
if ( $failed ) {
    print "model failed - results were:\n\n";
    system("cat outputs/results_baryonicSuppression.json");
}

exit;
