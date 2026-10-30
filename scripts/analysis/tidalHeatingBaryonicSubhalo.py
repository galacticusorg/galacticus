#!/usr/bin/env python3
"""Plot the effects of tidal heating on the spheroid of a subhalo orbiting in a fixed host potential, by comparing runs of
parameters/examples/tidalStrippingBaryonicSubhalo.xml with and without heating of the spheroid:

  ./Galacticus.exe parameters/examples/tidalStrippingBaryonicSubhalo.xml
  ./Galacticus.exe parameters/examples/tidalStrippingBaryonicSubhalo.xml parameters/examples/tidalStrippingBaryonicSubhaloUnheated.xml
  ./scripts/analysis/tidalHeatingBaryonicSubhalo.py tidalStrippingBaryonicSubhalo.hdf5 tidalStrippingBaryonicSubhaloUnheated.hdf5 \
      tidalHeatingBaryonicSubhalo.pdf

In both runs the dark matter of the subhalo is tidally heated, and its spheroid is tidally stripped and truncated. In the first
run the spheroid is also heated along with the dark matter (the "tracer" galacticComponentHeating class). Four panels are shown:

  1. The half-mass radius of the spheroid as a function of time.
  2. The stellar mass of the spheroid, relative to its initial value, as a function of time. Heating moves stars outward, where
     they are more easily stripped, so the heated spheroid loses more mass.
  3. The tidal track of the spheroid: its half-mass radius against its stellar mass, each relative to their initial values. At a
     given mass the heated spheroid is more extended.
  4. The density profile of the spheroid stars at several times, along with that of the dark matter at the final time.

Andrew Benson, Claude (2026)
"""

import argparse
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
parser.add_argument('fileNameHeated'  , help='Galacticus output file for the run with heating of the spheroid.'   )
parser.add_argument('fileNameUnheated', help='Galacticus output file for the run without heating of the spheroid.')
parser.add_argument('outputFile'      , help='Name of the plot file to write (format inferred from the extension).')
parser.add_argument('--nodeIndex', type=int, default=3, help='Index of the subhalo node to plot (default: %(default)s).')
parser.add_argument('--times', type=float, nargs='+', default=[3.0, 6.0, 13.0],
                    help='Times (in Gyr) at which to show density profiles (default: %(default)s).')
arguments = parser.parse_args()

kiloparsec = 1.0e3  # Galacticus radii are in Mpc.


def readHistory(fileName):
    """Read the history of the subhalo, and its density profiles, from a Galacticus output file."""
    names = ('time', 'spheroidMassStellar', 'satelliteBoundMass')
    history = {name: [] for name in names+('radiusHalfMass', 'densityStellar', 'densityDarkMatter')}
    radii = None
    with h5py.File(fileName, 'r') as model:
        outputs = sorted(model['Outputs'], key=lambda name: int(name.removeprefix('Output')))
        for output in outputs:
            nodeData = model['Outputs'][output]['nodeData']
            if 'nodeIndex' not in nodeData:
                continue
            selection = np.flatnonzero(nodeData['nodeIndex'][:] == arguments.nodeIndex)
            if selection.size == 0:
                continue
            i = selection[0]
            for name in names:
                history[name].append(nodeData[name][i])
            history['radiusHalfMass'].append(nodeData['massProfileRadius'][i][0])
            # Density profiles are of the spheroid stars and the dark matter, at radii given in the column labels.
            columns = [column.decode() for column in nodeData['densityProfileColumns'][:]]
            stellar = [j for j, column in enumerate(columns) if ':spheroid:stellar:' in column]
            dark    = [j for j, column in enumerate(columns) if ':darkHalo:dark:'    in column]
            if radii is None:
                radii = np.array([float(columns[j].split(':')[-1]) for j in stellar])
            history['densityStellar'   ].append(nodeData['densityProfile'][i][stellar])
            history['densityDarkMatter'].append(nodeData['densityProfile'][i][dark   ])
    if len(history['time']) == 0:
        raise SystemExit(f"node {arguments.nodeIndex} was not found in any output of '{fileName}'")
    history = {name: np.array(values) for name, values in history.items()}
    history['radii'] = radii
    # The spheroid has no size until it is first found by the galactic structure solver. Measure its initial half-mass radius at
    # the first output at which it has a size.
    sized = history['radiusHalfMass'] > 0.0
    history['radiusHalfMass'] = np.where(sized, history['radiusHalfMass'], np.nan)
    history['radiusHalfMassInitial'] = history['radiusHalfMass'][np.argmax(sized)]
    return history


runs = {'heated'  : {'history': readHistory(arguments.fileNameHeated  ), 'color': 'C1'  , 'linestyle': '-' ,
                     'label': 'Spheroid heated'    },
        'unheated': {'history': readHistory(arguments.fileNameUnheated), 'color': 'grey', 'linestyle': '--',
                     'label': 'Spheroid not heated'}}

figure, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))
axesRadius, axesMass, axesTrack, axesDensity = axes.flatten()

for run in runs.values():
    history = run['history']
    style   = {'color': run['color'], 'linestyle': run['linestyle'], 'label': run['label']}
    # Half-mass radius.
    axesRadius.plot(history['time'], history['radiusHalfMass']*kiloparsec, **style)
    # Stellar mass.
    axesMass.plot(history['time'], history['spheroidMassStellar']/history['spheroidMassStellar'][0], **style)
    # Tidal track.
    axesTrack.plot(history['spheroidMassStellar']/history['spheroidMassStellar'][0],
                   history['radiusHalfMass']/history['radiusHalfMassInitial'], **style)
# The bound dark matter mass is shown for reference.
history = runs['heated']['history']
axesMass.plot(history['time'], history['satelliteBoundMass']/history['satelliteBoundMass'][0], color='black', linestyle=':',
              label='Dark matter (bound)')

axesRadius.set_yscale('log')
axesRadius.set_xlabel('Time [Gyr]')
axesRadius.set_ylabel(r'Spheroid half-mass radius, $r_{1/2}$ [kpc]')
axesRadius.legend(fontsize='small')

axesMass.set_yscale('log')
axesMass.set_xlabel('Time [Gyr]')
axesMass.set_ylabel(r'$M/M_\mathrm{initial}$')
axesMass.legend(fontsize='small')

axesTrack.set_xscale('log')
axesTrack.set_yscale('log')
axesTrack.invert_xaxis()
axesTrack.set_xlabel(r'Spheroid stellar mass, $M_\star/M_{\star,\mathrm{initial}}$')
axesTrack.set_ylabel(r'$r_{1/2}/r_{1/2,\mathrm{initial}}$')
axesTrack.set_title('Tidal track (time runs right to left)', fontsize='small')
axesTrack.legend(fontsize='small')

# Density profiles at the requested times: one shade per time, darkening with time, solid for the heated run and dashed for the
# unheated run.
shades = plt.get_cmap('Oranges')(np.linspace(0.45, 0.95, len(arguments.times)))
for time, shade in zip(arguments.times, shades):
    for name, run in runs.items():
        history = run['history']
        i = np.argmin(np.abs(history['time']-time))
        density = np.where(history['densityStellar'][i] > 0.0, history['densityStellar'][i], np.nan)
        axesDensity.plot(history['radii']*kiloparsec, density/kiloparsec**3, color=shade, linestyle=run['linestyle'],
                         label=f"Stars, t={history['time'][i]:.1f} Gyr ({'heated' if name == 'heated' else 'not heated'})")
history = runs['heated']['history']
axesDensity.plot(history['radii']*kiloparsec, history['densityDarkMatter'][-1]/kiloparsec**3, color='black', linestyle=':',
                 label=f"Dark matter, t={history['time'][-1]:.1f} Gyr")
axesDensity.set_xscale('log')
axesDensity.set_yscale('log')
densityStellar = np.concatenate([run['history']['densityStellar'][run['history']['densityStellar'] > 0.0] for run in runs.values()])
axesDensity.set_ylim(1.0e-6*densityStellar.max()/kiloparsec**3, 3.0*densityStellar.max()/kiloparsec**3)
axesDensity.set_xlabel('Radius [kpc]')
axesDensity.set_ylabel(r'Density [M$_\odot$ kpc$^{-3}$]')
axesDensity.legend(fontsize='x-small')

figure.tight_layout()
figure.savefig(arguments.outputFile)
