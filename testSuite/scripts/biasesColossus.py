# Compute concentrations of halos under various models using Colossus.
# Andrew Benson (13-August-2019)
import numpy as np
from colossus.cosmology import cosmology
from colossus.lss import bias

# Specify a set of concentration models for comparison.
models = [ ( 'cole89', 'NA' ), ( 'sheth01', 'NA' ), ( 'tinker10', '200c' ) ]

# Build a range of halo masses at which to compute concentrations.
massLogarithmic = np.linspace(10.0,15.0,10)
mass            = 10.0**massLogarithmic

# Specify redshifts.
redshifts = ( 0.0, 0.5, 1.0 )

# Specify cosmology.
cosmology.setCosmology('planck15')
cosmo = cosmology.getCurrent()

# Iterate over models.
for model in models:

    outputFile = open("testSuite/data/haloBiasesColossus/"+model[0]+"_"+model[1]+".txt", "w")
    
    # Write cosmological parameters to file.
    outputFile.write("# OmegaMatter = "+str(cosmo.Om0)+"\n")
    outputFile.write("# OmegaDarkEnergy = "+str(cosmo.Ode0)+"\n")
    outputFile.write("# OmegaBaryon = "+str(cosmo.Ob0)+"\n")
    outputFile.write("# HubbleConstant = "+str(cosmo.H0)+"\n")
    outputFile.write("# sigma_8 = "+str(cosmo.sigma8)+"\n")
    outputFile.write("# index = "+str(cosmo.ns)+"\n")
    outputFile.write("# temperatureCMB = "+str(cosmo.Tcmb0)+"\n")
    outputFile.write("# effectiveNumberNeutrinos = "+str(cosmo.Neff)+"\n")

    # Write bias parameters to file.
    outputFile.write("# biasModel = "+str(model[0])+"\n")
    outputFile.write("# virialDensityContrast = "+str(model[1])+"\n")

    # Output bias.
    outputFile.write("#\n")
    outputFile.write("# mass\t\tredshift\tbias\n")

    # Iterate over redshifts.
    for redshift in redshifts:
    
        # Construct the list of concentrations.
        b = bias.haloBias(mass, model = model[0], z = redshift, mdef = model[1])

        for x in zip(mass,b):
            outputFile.write("%12.6e\t%3.1f\t\t%9.6f\n" % (x[0], redshift, x[1]))

    outputFile.close()
