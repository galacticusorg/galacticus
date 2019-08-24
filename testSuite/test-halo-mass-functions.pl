#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use PDL::IO::Misc;
use XML::Simple;
use Data::Dumper;

# Test dark matter halo mass functions against HMFcalc results.
# Andrew Benson (07-August-2019)

# Specify mass functions to test.
my @massFunctionTypes =
    (
     {
     	 label  => "Press-Schechter",
     	 method => "pressSchechter"
     },
     {
     	 label  => "Sheth-Tormen",
     	 method => "shethTormen"
     },
     {
     	 label  => "Tinker2008",
     	 method => "tinker2008"
     },
     # Disabled due to bug in HMFcalc implementation of Bhattacharya mass function.
     # {
     # 	 label  => "Bhattacharya",
     # 	 method => "bhattacharya2011"
     # }
    );

# Iterate over mass functions.
foreach my $massFunctionType ( @massFunctionTypes ) {

    # Specify HMFcalc directory.
    my $hmfCalcPath = "data/HMFcalc/".$massFunctionType->{'label'}."/";
    system("mkdir -p outputs/HMFcalc");

    # Read the HMFcalc parameter file.
    my %parametersHMFCalc;
    open(my $parameterFile,$hmfCalcPath."parameters.txt");
    while ( my $line = <$parameterFile> ) {
	next
	    unless ( $line =~ m/^([a-zA-Z0-9_]+):\s+([\d\.\+\-]+)\s*$/ );
	$parametersHMFCalc{$1} = $2;
    }
    close($parameterFile);

    # Read the HMFcalc transfer function file.
    (my $wavenumber, my $transferFunction) = rcols($hmfCalcPath."kVector_PLANCK-SMT .txt",0,2);
    $wavenumber *= $parametersHMFCalc{'h'};

    # Read the HMFcalc mass function file.
    (my $mass, my $massFunction) = rcols($hmfCalcPath."mVector_PLANCK-SMT .txt",0,6);
    $mass         /= $parametersHMFCalc{'h'}   ;
    $massFunction *= $parametersHMFCalc{'h'}**3;
    
    # Construct a parameter file for Galacticus.
    my $xml        = new XML::Simple();
    my $parameters = $xml->XMLin('parameters/haloMassFunctionsBase.xml');
    $parameters->{'taskMethod'                    }->{'haloMassMinimum'     }->{'value'} =      $mass->(( 0))->sclr();
    $parameters->{'taskMethod'                    }->{'haloMassMaximum'     }->{'value'} =      $mass->((-1))->sclr();
    $parameters->{'taskMethod'                    }->{'pointsPerDecade'     }->{'value'} =  1.0/$parametersHMFCalc{'dlog10m'};
    $parameters->{'cosmologyParametersMethod'     }->{'temperatureCMB'      }->{'value'} =      $parametersHMFCalc{'t_cmb'  };
    $parameters->{'cosmologyParametersMethod'     }->{'OmegaMatter'         }->{'value'} =      $parametersHMFCalc{'omegam' };
    $parameters->{'cosmologyParametersMethod'     }->{'OmegaDarkEnergy'     }->{'value'} =      $parametersHMFCalc{'omegav' };
    $parameters->{'cosmologyParametersMethod'     }->{'OmegaBaryon'         }->{'value'} =      $parametersHMFCalc{'omegab' };
    $parameters->{'cosmologyParametersMethod'     }->{'HubbleConstant'      }->{'value'} =      $parametersHMFCalc{'H0'     };
    $parameters->{'cosmologicalMassVarianceMethod'}->{'sigma_8'             }->{'value'} =      $parametersHMFCalc{'sigma_8'};
    $parameters->{'powerSpectrumPrimordialMethod' }->{'index'               }->{'value'} =      $parametersHMFCalc{'n'      };
    $parameters->{'criticalOverdensityMethod'     }->{'criticalOverdensity' }->{'value'} =      $parametersHMFCalc{'delta_c'};
    $parameters->{'virialDensityContrastMethod'   }->{'densityContrastValue'}->{'value'} =      $parametersHMFCalc{'delta_h'};
    $parameters->{'transferFunctionMethod'        }->{'fileName'            }->{'value'} = "testSuite/outputs/HMFcalc/".$massFunctionType->{'label'}."_Tk.hdf5" ;
    $parameters->{'haloMassFunctionMethod'        }                          ->{'value'} = $massFunctionType->{'method'};
    $parameters->{'galacticusOutputFileName'      }                          ->{'value'} = "testSuite/outputs/HMFcalc/".$massFunctionType->{'label'}."_HMF.hdf5";
    open(my $parameterOutputFile,">","outputs/HMFcalc/".$massFunctionType->{'label'}.".xml");
    print $parameterOutputFile $xml->XMLout($parameters, RootName => "parameters");
    close($parameterOutputFile);

    # Construct a transfer function file for Galacticus.
    my $transferFunctionFile = new PDL::IO::HDF5(">outputs/HMFcalc/".$massFunctionType->{'label'}."_Tk.hdf5");
    $transferFunctionFile                                             ->dataset('wavenumber'             )->set    (                     $wavenumber                     );
    $transferFunctionFile->group('darkMatter'   )                     ->dataset('transferFunctionZ0.0000')->set    (                     $transferFunction               );
    $transferFunctionFile                                                                                 ->attrSet("fileFormat"      => pdl long(2)                     );
    $transferFunctionFile->group('parameters'   )                                                         ->attrSet('OmegaMatter'     => pdl $parametersHMFCalc{'omegam'});
    $transferFunctionFile->group('parameters'   )                                                         ->attrSet('OmegaDarkEnergy' => pdl $parametersHMFCalc{'omegav'});
    $transferFunctionFile->group('parameters'   )                                                         ->attrSet('OmegaBaryon'     => pdl $parametersHMFCalc{'omegab'});
    $transferFunctionFile->group('parameters'   )                                                         ->attrSet('HubbleConstant'  => pdl $parametersHMFCalc{'H0'    });
    $transferFunctionFile->group('parameters'   )                                                         ->attrSet('temperatureCMB'  => pdl $parametersHMFCalc{'t_cmb' });
    $transferFunctionFile->group('extrapolation')->group('wavenumber')                                    ->attrSet('low'             => 'fix'                           );
    $transferFunctionFile->group('extrapolation')->group('wavenumber')                                    ->attrSet('high'            => 'extrapolate'                   );
    undef($transferFunctionFile);

    # Run Galacticus to generate the mass function.
    system("cd ..; Galacticus.exe testSuite/outputs/HMFcalc/".$massFunctionType->{'label'}.".xml");
    unless ( $? == 0 ) {
    	print "FAILED: Galacticus failed\n";
    	exit 0;
    }

    # Read the mass function generated by Galacticus.
    my $massFunctionFile       = new PDL::IO::HDF5("outputs/HMFcalc/".$massFunctionType->{'label'}."_HMF.hdf5");
    my $massGalacticus         = $massFunctionFile->group('Outputs')->group('Output1')->dataset('haloMass'           )->get();
    my $massFunctionGalacticus = $massFunctionFile->group('Outputs')->group('Output1')->dataset('haloMassFunctionLnM')->get();

    # Check that masses are consistent.
    if ( nelem($mass) != nelem($massGalacticus) ) {
	print "FAILED: [".$massFunctionType->{'label'}."] number of masses differ\n";
	exit 0;
    }
    if ( any(abs($mass-$massGalacticus)/$mass > 1.0e-6) ) {
	print "FAILED: [".$massFunctionType->{'label'}."]masses are inconsistent\n";
	exit 0;
    }

    # Compare mass function with that from HMFcalc.
    my $differenceFractional       = abs($massFunction-$massFunctionGalacticus)/$massFunction;
    my $differenceFractionlMaximum = $differenceFractional->maximum();    
    if ( $differenceFractionlMaximum > 1.0e-3 ) {
	print "FAILED: [" .$massFunctionType->{'label'}."] mass function differs\n";
    } else {
	print "success: [".$massFunctionType->{'label'}."] mass functions agree\n";
    }

}

exit 0;
