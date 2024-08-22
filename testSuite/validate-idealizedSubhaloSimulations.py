#!/usr/bin/env python3
import sys
import os
import subprocess
import numpy as np
import h5py
import validate

# Run models to validate a idealized subhalo simulations models.
# Andrew Benson (21-March-2024)

# Define the models to validate.
models = [
     {
	 "fileName"         : "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.005_alpha_1.0_beta_3.0_gamma_1.5.hdf5"                   ,
	 "name"             : "Idealized Subhalo Simulation (rₚ/rₐ=0.005; γ=1.5)"                                                                     ,
	 "suffix"           : "idealizedSubhaloSimulation_rpra0.005_gamma1.5"                                                                         ,
	 "parameterFileName": "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.005_alpha_1.0_beta_3.0_gamma_1.5.xml"
     },
     {
         "fileName"         : "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.01_alpha_1.0_beta_3.0_gamma_1.5.hdf5"                    ,
         "name"             : "Idealized Subhalo Simulation (rₚ/rₐ=0.01; γ=1.5)"                                                                      ,
         "suffix"           : "idealizedSubhaloSimulation_rpra0.01_gamma1.5"                                                                          ,
         "parameterFileName": "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.01_alpha_1.0_beta_3.0_gamma_1.5.xml"
     },
     {
         "fileName"         : "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.05_alpha_1.0_beta_3.0_gamma_0.5.hdf5"                    ,
         "name"             : "Idealized Subhalo Simulation (rₚ/rₐ=0.05; γ=0.5)"                                                                      ,
         "suffix"           : "idealizedSubhaloSimulation_rpra0.05_gamma0.5"                                                                          ,
         "parameterFileName": "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.05_alpha_1.0_beta_3.0_gamma_0.5.xml"
     },
     {
         "fileName"         : "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.05_alpha_1.0_beta_3.0_gamma_1.0.hdf5"                    ,
         "name"             : "Idealized Subhalo Simulation (rₚ/rₐ=0.05; γ=1.0)"                                                                      ,
         "suffix"           : "idealizedSubhaloSimulation_rpra0.05_gamma1.0"                                                                          ,
         "parameterFileName": "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.05_alpha_1.0_beta_3.0_gamma_1.0.xml"
     },
     {
         "fileName"         : "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_0.0.hdf5"                     ,
         "name"             : "Idealized Subhalo Simulation (rₚ/rₐ=0.2; γ=0.0)"                                                                       ,
         "suffix"           : "idealizedSubhaloSimulation_rpra0.2_gamma0.0"                                                                           ,
         "parameterFileName": "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_0.0.xml"
     },
     {
         "fileName"         : "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_0.5.hdf5"                     ,
         "name"             : "Idealized Subhalo Simulation (rₚ/rₐ=0.2; γ=0.5)"                                                                       ,
         "suffix"           : "idealizedSubhaloSimulation_rpra0.2_gamma0.5"                                                                           ,
         "parameterFileName": "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_0.5.xml"
     },
     {
         "fileName"         : "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_1.0.hdf5"                     ,
         "name"             : "Idealized Subhalo Simulation (rₚ/rₐ=0.2; γ=1.0)"                                                                       ,
         "suffix"           : "idealizedSubhaloSimulation_rpra0.2_gamma1.0"                                                                           ,
         "parameterFileName": "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.2_alpha_1.0_beta_3.0_gamma_1.0.xml"
     },
     {
         "fileName"         : "outputs/idealizedSubhaloSimulations/tidalTrack_xc_0.7_ratio_0.4_alpha_1.0_beta_3.0_gamma_0.0.hdf5"                     ,
         "name"             : "Idealized Subhalo Simulation (rₚ/rₐ=0.4; γ=0.0)"                                                                       ,
         "suffix"           : "idealizedSubhaloSimulation_rpra0.4_gamma0.0"                                                                           ,
         "parameterFileName": "testSuite/parameters/idealizedSubhaloSimulations/tidalTrackBestFit_xc_0.7_ratio_0.4_alpha_1.0_beta_3.0_gamma_0.0.xml"
     }
    ]

# Create output path.
try:
    os.mkdir("outputs/idealizedSubhaloSimulations")
except FileExistsError:
    pass
# Iterate over models.
for model in models:
    # Run the validation model.
    status = subprocess.run("cd ..; ./Galacticus.exe "+model['parameterFileName'],shell=True)
    if status.returncode != 0:
        print("FAILED: idealized subhalo validation model '"+model['suffix']+"' failed to run")
        sys.exit()
    # Extract and validate the likelihoods.
    validate.extract(
        model['fileName'         ],
        model['name'             ],
        model['suffix'           ],
        model['parameterFileName']
    )

print("SUCCESS: idealized subhalo simulations validation model")
