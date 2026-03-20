#!/usr/bin/env python3
import numpy as np
import h5py
from scipy import integrate

# Compute the fraction of particles, and kick energy retained in a halo for decaying dark matter models. Uses a Monte Carlo
# simulation to evaluate these fractions, which is used as a target dataset to test our numerical calculations.
# Andrew Benson (01-May-2024; ported to Python)

# Specify the velocity dispersion and the escape velocity.
velocity_dispersion = 20.0
velocity_escape     = 60.0

# Generate a sequence of velocity kicks for which to compute the retained fractions.
velocity_kicks   = np.arange(1, 131, dtype=float)  # 1..130

# Generate a sequence of escape velocities for which to compute the retained fractions.
velocity_escapes = np.arange(60, 70,  dtype=float)  # 60..69

# Specify the number of Monte Carlo trials.
count_trials = 30

# Specify the number of particles to sample for each trial.
count_sample = 1_000_000

# Create arrays to store results.
fractions_retained = np.zeros((count_trials, len(velocity_kicks), len(velocity_escapes)))
energies_retained  = np.zeros((count_trials, len(velocity_kicks), len(velocity_escapes)))

# Iterate over escape velocities.
velocity_width = velocity_dispersion
for j, v_escape in enumerate(velocity_escapes):
    print(f"Velocity escape = {v_escape} [{j+1} of {len(velocity_escapes)}]")

    # Compute the width of the truncated Maxwell-Boltzmann distribution (truncated at the escape
    # velocity) required to give a mean squared velocity of 3 times the velocity dispersion. We
    # use a simple iterative scheme here.
    velocity_squared_mean_target = 3.0 * velocity_dispersion**2
    boost_factor = 1.0
    velocity_mean_squared = None
    for _ in range(10):
        velocity_width = boost_factor * velocity_dispersion

        def velocity_mean_squared_integrand(v):
            return (np.sqrt(2.0 / np.pi) * velocity_width
                    * (v / velocity_width)**4
                    * np.exp(-0.5 * (v / velocity_width)**2))

        def velocity_normalization_integrand(v):
            return (np.sqrt(2.0 / np.pi)
                    * (v / velocity_width)**2 / velocity_width
                    * np.exp(-0.5 * (v / velocity_width)**2))

        velocity_mean_squared, _ = integrate.quad(
            velocity_mean_squared_integrand, 0.0, v_escape,
            limit=1000, epsrel=1.0e-4, epsabs=0.0
        )
        normalization, _ = integrate.quad(
            velocity_normalization_integrand, 0.0, v_escape,
            limit=1000, epsrel=1.0e-4, epsabs=0.0
        )
        velocity_mean_squared /= normalization
        boost_factor *= np.sqrt(velocity_squared_mean_target / velocity_mean_squared)

    print(f"\tComputed boost factor of: {boost_factor}")
    print(f"\t\tMean squared velocity is {velocity_mean_squared} (target was {velocity_squared_mean_target})")

    # Construct the truncated Maxwell-Boltzmann distribution and find the cumulative distribution function.
    count_velocity = 100_000
    velocity       = np.arange(count_velocity) / (count_velocity - 1) * v_escape
    pdf            = velocity**2 * np.exp(-0.5 * (velocity / velocity_width)**2)
    cdf            = np.cumsum(pdf) / pdf.sum()

    # Iterate over velocity kicks.
    for i, v_kick in enumerate(velocity_kicks):
        print(f"\tVelocity kick = {v_kick} [{i+1} of {len(velocity_kicks)}]")
        # Iterate over trials.
        for trial in range(count_trials):
            # Generate random velocities and kick angles.
            x1 = np.random.uniform(size=count_sample)
            x2 = np.random.uniform(size=count_sample)
            velocities = np.interp(x1, cdf, velocity)
            cos_theta  = 2.0 * x2 - 1.0
            # Compute the final velocity.
            velocity_final = np.sqrt(velocities**2 + v_kick**2 + 2.0 * velocities * v_kick * cos_theta)
            # Find those particles that are retained.
            retained = velocity_final < v_escape
            # Compute the energy gain of retained particles.
            energy_gain = np.sum(
                0.5 * velocity_final[retained]**2 - 0.5 * velocities[retained]**2
            ) / count_sample
            # Store the retained fractions.
            energies_retained [trial, i, j] = energy_gain / (0.5 * v_kick**2)
            fractions_retained[trial, i, j] = retained.sum() / count_sample

# Compute the mean fractions and their uncertainties.
fraction_retained             = fractions_retained.mean(axis=0)
fraction_retained_uncertainty = np.sqrt(
    ((fractions_retained**2).mean(axis=0) - fraction_retained**2) / count_trials
)
energy_retained               = energies_retained.mean(axis=0)
energy_retained_uncertainty   = np.sqrt(
    ((energies_retained**2).mean(axis=0) - energy_retained**2) / count_trials
)

# Store results to file.
with h5py.File("testSuite/data/decayingDarkMatterRetention.hdf5", "w") as f:
    f.create_dataset("velocityKick",                data=velocity_kicks)
    f.create_dataset("velocityEscape",              data=velocity_escapes)
    f.create_dataset("fractionRetained",            data=fraction_retained)
    f.create_dataset("fractionRetainedUncertainty", data=fraction_retained_uncertainty)
    f.create_dataset("energyRetained",              data=energy_retained)
    f.create_dataset("energyRetainedUncertainty",   data=energy_retained_uncertainty)
    f.attrs["velocityDispersion"] = velocity_dispersion
