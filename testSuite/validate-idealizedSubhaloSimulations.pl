#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Galacticus::Validation;

# Run models to validate a idealized subhalo simulations models.
# Andrew Benson (21-March-2024)

# WORKAROUND: Attempting to run all models from this script results in
# a failure in PDL::IO::HDF5. To workaround this we accept a command
# line argument that specifies which model to run, and run only a
# single model. This script is then driven by a wrapper script that
# runs it for each individual model.
die("Usage: validate-idealizedSubhaloSimulations.pl <modelNumber>")
    unless ( scalar(@ARGV) == 1 );
my $i = $ARGV[0];

# Make output directory.
system("mkdir -p outputs/idealizedSubhaloSimulations");

# Define the models to validate.
my @models =
    (
     {
	 fileName          => "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.005_alpha_1.0_beta_3.0_gamma_1.5.hdf5"                   ,
	 name              => "Idealized Subhalo Simulation (rₚ/rₐ=0.005; γ=1.5)"                                                                     ,
	 suffix            => "idealizedSubhaloSimulation_rpra0.005_gamma1.5"                                                                         ,
	 parameterFileName => "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.005_alpha_1.0_beta_3.0_gamma_1.5.xml"
     },
     {
	 fileName          => "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.01_alpha_1.0_beta_3.0_gamma_1.5.hdf5"                    ,
	 name              => "Idealized Subhalo Simulation (rₚ/rₐ=0.01; γ=1.5)"                                                                      ,
	 suffix            => "idealizedSubhaloSimulation_rpra0.01_gamma1.5"                                                                          ,
	 parameterFileName => "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.01_alpha_1.0_beta_3.0_gamma_1.5.xml"
     },
     {
	 fileName          => "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.05_alpha_1.0_beta_3.0_gamma_0.5.hdf5"                    ,
	 name              => "Idealized Subhalo Simulation (rₚ/rₐ=0.05; γ=0.5)"                                                                      ,
	 suffix            => "idealizedSubhaloSimulation_rpra0.05_gamma0.5"                                                                          ,
	 parameterFileName => "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.05_alpha_1.0_beta_3.0_gamma_0.5.xml"
     },
     {
	 fileName          => "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.05_alpha_1.0_beta_3.0_gamma_1.0.hdf5"                    ,
	 name              => "Idealized Subhalo Simulation (rₚ/rₐ=0.05; γ=1.0)"                                                                      ,
	 suffix            => "idealizedSubhaloSimulation_rpra0.05_gamma1.0"                                                                          ,
	 parameterFileName => "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.05_alpha_1.0_beta_3.0_gamma_1.0.xml"
     },
     {
	 fileName          => "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_0.0.hdf5"                     ,
	 name              => "Idealized Subhalo Simulation (rₚ/rₐ=0.2; γ=0.0)"                                                                       ,
	 suffix            => "idealizedSubhaloSimulation_rpra0.2_gamma0.0"                                                                           ,
	 parameterFileName => "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_0.0.xml"
     },
     {
	 fileName          => "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_0.5.hdf5"                     ,
	 name              => "Idealized Subhalo Simulation (rₚ/rₐ=0.2; γ=0.5)"                                                                       ,
	 suffix            => "idealizedSubhaloSimulation_rpra0.2_gamma0.5"                                                                           ,
	 parameterFileName => "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_0.5.xml"
     },
     {
	 fileName          => "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_1.0.hdf5"                     ,
	 name              => "Idealized Subhalo Simulation (rₚ/rₐ=0.2; γ=1.0)"                                                                       ,
	 suffix            => "idealizedSubhaloSimulation_rpra0.2_gamma1.0"                                                                           ,
	 parameterFileName => "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_1.0.xml"
     },
     {
	 fileName          => "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.4_alpha_1.0_beta_3.0_gamma_0.0.hdf5"                     ,
	 name              => "Idealized Subhalo Simulation (rₚ/rₐ=0.4; γ=0.0)"                                                                       ,
	 suffix            => "idealizedSubhaloSimulation_rpra0.4_gamma0.0"                                                                           ,
	 parameterFileName => "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.4_alpha_1.0_beta_3.0_gamma_0.0.xml"
     }
    );
die("model number is out of range")
    unless ( $i >= 0 && $i < scalar(@models) );

# Run the validation model.
system("cd ..; ./Galacticus.exe ".$models[$i]->{'parameterFileName'});
unless ( $? == 0 ) {
    print "FAIL: idealized subhalo validation model '".$models[$i]->{'suffix'}."' failed to run\n";
    exit;
}    
# Extract and validate the likelihoods.
&Galacticus::Validation::extract
    (
     $models[$i]->{'fileName'         },
     $models[$i]->{'name'             },
     $models[$i]->{'suffix'           },
     $models[$i]->{'parameterFileName'}
    );

print "SUCCESS: idealized subhalo simulations validation model\n";

exit;
