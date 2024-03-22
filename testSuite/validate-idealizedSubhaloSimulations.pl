#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Galacticus::Validation;

# Run models to validate a idealized subhalo simulations models.
# Andrew Benson (21-March-2024)

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

foreach my $model ( @models ) {
    # Run the validation model.
    system("cd ..; ./Galacticus.exe ".$model->{'parameterFileName'});
    unless ( $? == 0 ) {
	print "FAIL: idealized subhalo validation model '".$model->{'suffix'}."' failed to run\n";
	exit;
    }    
    # Pause to ensure file is ready.
    sleep(10);
    # Extract and validate the likelihoods.
    &Galacticus::Validation::extract
	(
	 $model->{'fileName'         },
	 $model->{'name'             },
	 $model->{'suffix'           },
	 $model->{'parameterFileName'}
	);
}

print "SUCCESS: dark matter-only subhalos validation model\n";

exit;
