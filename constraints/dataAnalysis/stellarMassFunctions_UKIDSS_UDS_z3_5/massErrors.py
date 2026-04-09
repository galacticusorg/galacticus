#!/usr/bin/env python3
# Estimate random errors on stellar masses in the UKIDSS UDS survey of Caputi et al. (2013).
# Andrew Benson (07-May-2014)

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

# Central redshifts of the three redshift bins considered.
redshift = np.array([3.250, 3.875, 4.625])

# Fractional error in dz/(1+z) reported by Caputi et al. (2013).
sigmaLogZ = 0.05

# Speed of light in km/s.
speedLight = 2.99792458e5

# Construct cosmological model.
cosmology = FlatLambdaCDM(H0=69.5723630486537, Om0=0.283812448723631)

# Find comoving distances to each redshift (in Mpc).
distanceComoving = cosmology.comoving_distance(redshift).to(u.Mpc).value

# Find Hubble parameter at each redshift (in km/s/Mpc).
hubbleParameter = (cosmology.H(redshift)).to(u.km / u.s / u.Mpc).value

# Compute the error in log10(M*) from photometric redshift errors.
sigmaLog10MassRedshift = (
    2.0 + 2.0 * (1.0 + redshift)**2 * speedLight / hubbleParameter / distanceComoving
) * sigmaLogZ / np.log(10.0)

# Additional error arising from SED fitting (judged from Figure 8 of Caputi et al. 2013).
sigmaLog10MassSED = 0.2 / np.log(10.0)

# Compute final error.
sigmaLog10Mass = np.sqrt(sigmaLog10MassRedshift**2 + sigmaLog10MassSED**2)

# Display the results.
print('Redshift\tDispersion in log10 stellar mass')
for i in range(len(redshift)):
    print(f'{redshift[i]}\t{sigmaLog10Mass[i]}')
