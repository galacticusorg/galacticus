#!/usr/bin/env python3
# Compute the mean and variance of the subhalo halo occupation distribution slope from the
# simulation results reported by Kravtsov et al. (2004, ApJ, 609, 35;
# http://adsabs.harvard.edu/abs/2004ApJ...609...35K).
# Andrew Benson (23-July-2012)

import numpy as np

# Set arrays of the reported slopes, alpha, and the errors on these values.
alpha = np.array([0.99, 0.92, 0.96, 1.04, 0.61])
error = np.array([0.01, 0.03, 0.08, 0.08, 0.21])

# Compute a weight equal to the inverse variance on each measurement.
weight = 1.0 / error**2

# Compute weighted mean and standard deviation.
totalWeight = weight.sum()
mean        = np.sum(weight * alpha) / totalWeight
variance    = np.sum(weight * (alpha - mean)**2) / totalWeight
rms         = np.sqrt(variance)

# Report.
print('Statistics of the slope of the subhalo halo occupation distribution.')
print(f'                Mean: {mean}')
print(f'  Standard deviation: {rms}')
print(f'            Variance: {variance}')
