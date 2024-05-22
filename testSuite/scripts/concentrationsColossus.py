# Compute concentrations of halos under various models using Colossus.
# Andrew Benson (13-August-2019)
import numpy as np
from colossus.cosmology import cosmology
from colossus.halo import concentration

# Specify a set of concentration models for comparison.
models = [ ( 'prada12', '200c' ), ( 'dutton14', '200c' ), ( 'dutton14', 'vir' ), ( 'diemer15_orig', '200c' ), ( 'diemer19', '200c'), ( 'diemer19', 'vir'), ( 'ludlow16', '200c' ) ]

# Specify a set of redshifts at which to compute concentrations.
redshifts = [ 0.0, 1.0 ]

# Build a range of halo masses at which to compute concentrations.
massLogarithmic = np.linspace(6.0,15.0,19)
mass            = 10.0**massLogarithmic

# Specify cosmology.
cosmology.setCosmology('planck15')
cosmo = cosmology.getCurrent()

# Iterate over models.
for model in models:

    for redshift in redshifts:

        outputFile = open("testSuite/data/concentrationsColossus/"+model[0]+"_"+model[1]+"_z"+str(redshift)+".txt", "w")
    
        # Write cosmological parameters to file.
        outputFile.write("# OmegaMatter = "+str(cosmo.Om0)+"\n")
        outputFile.write("# OmegaDarkEnergy = "+str(cosmo.Ode0)+"\n")
        outputFile.write("# OmegaBaryon = "+str(cosmo.Ob0)+"\n")
        outputFile.write("# HubbleConstant = "+str(cosmo.H0)+"\n")
        outputFile.write("# sigma_8 = "+str(cosmo.sigma8)+"\n")
        outputFile.write("# index = "+str(cosmo.ns)+"\n")
        outputFile.write("# temperatureCMB = "+str(cosmo.Tcmb0)+"\n")
        outputFile.write("# effectiveNumberNeutrinos = "+str(cosmo.Neff)+"\n")

        # Write concentration parameters to file.
        outputFile.write("# concentrationModel = "+str(model[0])+"\n")
        outputFile.write("# virialDensityContrast = "+str(model[1])+"\n")

        # Construct the list of concentrations.
        cvir = concentration.concentration(mass, model[1], redshift, model = model[0])

        # Output concentrations.
        outputFile.write("#\n")
        outputFile.write("# mass\t\tconcentration\n")
        for x in zip(mass,cvir):
            outputFile.write("%12.6e\t%9.6f\n" % x)

        outputFile.close()
