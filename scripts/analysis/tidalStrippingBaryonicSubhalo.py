#!/usr/bin/env python3
"""Plot the evolution of a subhalo, containing a stellar disk and spheroid, which is tidally stripped as it orbits in a fixed host
potential. Intended for use with the output of parameters/examples/tidalStrippingBaryonicSubhalo.xml:

  ./Galacticus.exe parameters/examples/tidalStrippingBaryonicSubhalo.xml
  ./scripts/analysis/tidalStrippingBaryonicSubhalo.py tidalStrippingBaryonicSubhalo.hdf5 tidalStrippingBaryonicSubhalo.pdf

Three panels are shown as a function of time:

  1. Radii: the orbital radius of the subhalo, its tidal radius, the scale radii of its disk and spheroid, and the radius at which
     the mass distribution of the spheroid is truncated.
  2. Masses: the bound dark matter mass of the subhalo, and the stellar masses of its disk and spheroid, each relative to its
     initial value.
  3. Degree of truncation: the fraction, f = exp(-x), of the mass of the untruncated mass distribution retained by each of the disk
     and spheroid (from the "tidalTruncation" meta-properties). In this model only tidal stripping changes the masses of the disk and
     spheroid, so f should equal the fraction of the initial stellar mass retained, which is shown for comparison.

The spheroid truncation radius is not output by Galacticus, so is reconstructed here from the degree of truncation, x, for the case
of a Hernquist spheroid truncated by the "exponential" galacticComponentTruncation class (the model used by the example parameter
file). It is omitted for any other model.

Andrew Benson, Claude (2026)
"""

import argparse
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import gamma, gammaincc

parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
parser.add_argument('fileName'  , help='Galacticus output file to read.')
parser.add_argument('outputFile', help='Name of the plot file to write (format inferred from the extension).')
parser.add_argument('--nodeIndex', type=int, default=3, help='Index of the subhalo node to plot (default: %(default)s).')
arguments = parser.parse_args()


def attribute(group, name):
    """Return a parameter value (stored as an HDF5 attribute) as a string."""
    value = group.attrs[name]
    return value.decode() if isinstance(value, bytes) else str(value)


def hernquistTruncationRadius(truncation, ratioDecay):
    """Return the truncation radius (in units of the scale radius) of a dimensionless Hernquist profile exponentially truncated
    (as in the massDistributionSphericalTruncatedExponential class, with r_decay = ratioDecay × r_t) such that it retains a
    fraction exp(-truncation) of its mass. Computed by tabulating the retained mass as a function of truncation radius and
    interpolating."""
    radius = np.logspace(-7.0, 6.0, 2000)
    radiusDecay = ratioDecay*radius
    # Logarithmic slope of the Hernquist density profile, and the exponent κ which makes the logarithmic slope continuous.
    slope = -1.0-3.0*radius/(1.0+radius)
    kappa = radius/radiusDecay+slope
    density = 1.0/(2.0*np.pi*radius*(1.0+radius)**3)
    massInner = radius**2/(1.0+radius)**2
    x = radius/radiusDecay
    # Mass beyond the truncation radius: 4π ρ(r_t) r_d³ exp(x_t) x_t^-κ Γ(3+κ,x_t), with Γ the upper incomplete gamma function.
    massOuter = 4.0*np.pi*density*radiusDecay**3*np.exp(x-kappa*np.log(x))*gammaincc(3.0+kappa, x)*gamma(3.0+kappa)
    truncationTabulated = -np.log(massInner+massOuter)
    # Truncation decreases monotonically with radius - reverse for interpolation.
    order = np.argsort(truncationTabulated)
    radiusTruncation = np.full_like(truncation, np.nan)
    # Truncations below 10⁻⁴ are treated as negligible by Galacticus, so no truncation radius is shown for them.
    truncated = truncation >= 1.0e-4
    radiusTruncation[truncated] = np.exp(np.interp(np.log(truncation[truncated]),
                                                   np.log(truncationTabulated[order]), np.log(radius[order])))
    return radiusTruncation


with h5py.File(arguments.fileName, 'r') as model:
    # Extract the history of the subhalo from the outputs.
    outputs = sorted(model['Outputs'], key=lambda name: int(name.removeprefix('Output')))
    names = ('time', 'radiusOrbitalCurrent', 'satelliteRadiusTidal', 'satelliteBoundMass', 'diskMassStellar', 'diskRadius',
             'spheroidMassStellar', 'spheroidRadius', 'diskTidalTruncation', 'spheroidTidalTruncation')
    history = {name: [] for name in names}
    for output in outputs:
        nodeData = model['Outputs'][output]['nodeData']
        if 'nodeIndex' not in nodeData:
            continue
        selection = np.flatnonzero(nodeData['nodeIndex'][:] == arguments.nodeIndex)
        if selection.size == 0:
            continue
        for name in names:
            history[name].append(nodeData[name][selection[0]])
    if len(history['time']) == 0:
        raise SystemExit(f"node {arguments.nodeIndex} was not found in any output")
    history = {name: np.array(values) for name, values in history.items()}
    # Determine whether the spheroid truncation radius can be reconstructed.
    parameters = model['Parameters']
    radiusTruncationSpheroid = None
    if attribute(parameters, 'galacticComponentTruncation') == 'exponential' \
       and attribute(parameters['componentSpheroid'], 'massDistributionSpheroid') == 'hernquist':
        ratioDecay = float(attribute(parameters['galacticComponentTruncation'], 'radiusDecayOverRadiusTruncate'))
        radiusTruncationSpheroid = hernquistTruncationRadius(history['spheroidTidalTruncation'], ratioDecay) \
            * history['spheroidRadius']

kiloparsec = 1.0e3  # Galacticus radii are in Mpc.
# Component sizes are zero until first found by the galactic structure solver - exclude these from the plot.
for name in ('diskRadius', 'spheroidRadius'):
    history[name] = np.where(history[name] > 0.0, history[name], np.nan)
if radiusTruncationSpheroid is not None:
    radiusTruncationSpheroid = np.where(radiusTruncationSpheroid > 0.0, radiusTruncationSpheroid, np.nan)
time = history['time']
figure, axes = plt.subplots(3, 1, sharex=True, figsize=(7.0, 10.0))

# Radii.
axes[0].plot(time, history['radiusOrbitalCurrent']*kiloparsec, color='black', label='Orbital radius')
axes[0].plot(time, history['satelliteRadiusTidal']*kiloparsec, color='grey' , label='Tidal radius')
axes[0].plot(time, history['diskRadius'          ]*kiloparsec, color='C0'   , label='Disk scale radius')
axes[0].plot(time, history['spheroidRadius'      ]*kiloparsec, color='C1'   , label='Spheroid scale radius')
if radiusTruncationSpheroid is not None:
    axes[0].plot(time, radiusTruncationSpheroid*kiloparsec, color='C1', linestyle='--', label='Spheroid truncation radius')
axes[0].set_yscale('log')
axes[0].set_ylabel('Radius [kpc]')
axes[0].legend(fontsize='small')

# Masses relative to their initial values.
for name, color, label in (('satelliteBoundMass' , 'black', 'Dark matter (bound)'),
                           ('diskMassStellar'    , 'C0'   , 'Disk stars'         ),
                           ('spheroidMassStellar', 'C1'   , 'Spheroid stars'     )):
    axes[1].plot(time, history[name]/history[name][0], color=color, label=label)
axes[1].set_yscale('log')
axes[1].set_ylabel(r'$M/M_\mathrm{initial}$')
axes[1].legend(fontsize='small')

# Degree of truncation, compared with the fraction of stellar mass retained.
for component, color in (('disk', 'C0'), ('spheroid', 'C1')):
    axes[2].plot(time, np.exp(-history[component+'TidalTruncation']), color=color,
                 label=f"{component.capitalize()}: $f=\\exp(-x)$")
    axes[2].plot(time, history[component+'MassStellar']/history[component+'MassStellar'][0], color=color, linestyle=':',
                 label=f"{component.capitalize()}: $M_\\star/M_{{\\star,\\mathrm{{initial}}}}$")
axes[2].set_yscale('log')
axes[2].set_ylabel('Fraction retained')
axes[2].set_xlabel('Time [Gyr]')
axes[2].legend(fontsize='small')

figure.tight_layout()
figure.savefig(arguments.outputFile)
