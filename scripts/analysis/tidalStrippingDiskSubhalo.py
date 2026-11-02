#!/usr/bin/env python3
"""Plot the tidal stripping, truncation, and heating of the disk of a subhalo orbiting in a fixed host potential, comparing runs of
parameters/examples/tidalStrippingBaryonicSubhalo.xml with different treatments of the disk, for example:

  ./Galacticus.exe parameters/examples/tidalStrippingBaryonicSubhalo.xml
  ./Galacticus.exe parameters/examples/tidalStrippingBaryonicSubhalo.xml parameters/examples/tidalStrippingBaryonicSubhaloUnheated.xml
  ./scripts/analysis/tidalStrippingDiskSubhalo.py tidalStrippingDiskSubhalo.pdf \\
      tidalStrippingBaryonicSubhalo.hdf5:Heated tidalStrippingBaryonicSubhaloUnheated.hdf5:"Not heated"

Each run is given as a Galacticus output file, optionally followed by a colon and a label. Four panels are shown:

  1. The stellar mass of the disk, relative to its initial value, as a function of time.
  2. The half-mass radius of the disk as a function of time.
  3. The tidal track of the disk: its half-mass radius against its stellar mass, each relative to their initial values.
  4. The degree of truncation of the disk, x=-log(f), where f is the fraction of the mass of its untruncated mass distribution
     which it retains, and its heating parameter, s, as functions of time.

Andrew Benson, Claude (2026)
"""

import argparse
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
parser.add_argument('outputFile', help='Name of the plot file to write (format inferred from the extension).')
parser.add_argument('runs'      , nargs='+', help='Galacticus output files, each optionally followed by ":label".')
parser.add_argument('--nodeIndex', type=int, default=3, help='Index of the subhalo node to plot (default: %(default)s).')
arguments = parser.parse_args()

kiloparsec = 1.0e3  # Galacticus radii are in Mpc.


def readHistory(fileName):
    """Read the history of the disk of the subhalo from a Galacticus output file."""
    names   = ('time', 'diskMassStellar', 'diskTidalTruncation', 'diskTidalHeating')
    history = {name: [] for name in names+('radiusHalfMass',)}
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
            # The disk half-mass radius is the second of the radii output by the massProfile extractor.
            history['radiusHalfMass'].append(nodeData['massProfileRadius'][i][1])
    if len(history['time']) == 0:
        raise SystemExit(f"node {arguments.nodeIndex} was not found in any output of '{fileName}'")
    history = {name: np.array(values) for name, values in history.items()}
    # The disk has no size until it is first found by the galactic structure solver. Measure its initial half-mass radius at the
    # first output at which it has a size.
    sized = history['radiusHalfMass'] > 0.0
    history['radiusHalfMass'] = np.where(sized, history['radiusHalfMass'], np.nan)
    history['radiusHalfMassInitial'] = history['radiusHalfMass'][np.argmax(sized)]
    return history


runs = []
for i, run in enumerate(arguments.runs):
    fileName, _, label = run.partition(':')
    runs.append({'history': readHistory(fileName), 'label': label if label else fileName, 'color': f'C{i}'})

figure, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))
axesMass, axesRadius, axesTrack, axesTidal = axes.flatten()

for run in runs:
    history = run['history']
    style   = {'color': run['color'], 'label': run['label']}
    massRelative = history['diskMassStellar']/history['diskMassStellar'][0]
    axesMass  .plot(history['time'], massRelative                                          , **style)
    axesRadius.plot(history['time'], history['radiusHalfMass']*kiloparsec                  , **style)
    axesTrack .plot(massRelative   , history['radiusHalfMass']/history['radiusHalfMassInitial'], **style)
    axesTidal .plot(history['time'], history['diskTidalTruncation'], linestyle='-' , color=run['color'], label=f"x ({run['label']})")
    axesTidal .plot(history['time'], history['diskTidalHeating'   ], linestyle='--', color=run['color'], label=f"s ({run['label']})")

axesMass.set_yscale('log')
axesMass.set_xlabel('Time [Gyr]')
axesMass.set_ylabel(r'$M_\star/M_{\star,\mathrm{initial}}$')
axesMass.legend(fontsize='small')

axesRadius.set_yscale('log')
axesRadius.set_xlabel('Time [Gyr]')
axesRadius.set_ylabel(r'Disk half-mass radius, $R_{1/2}$ [kpc]')
axesRadius.legend(fontsize='small')

axesTrack.set_xscale('log')
axesTrack.set_yscale('log')
axesTrack.invert_xaxis()
axesTrack.set_xlabel(r'Disk stellar mass, $M_\star/M_{\star,\mathrm{initial}}$')
axesTrack.set_ylabel(r'$R_{1/2}/R_{1/2,\mathrm{initial}}$')
axesTrack.set_title('Tidal track (time runs right to left)', fontsize='small')
axesTrack.legend(fontsize='small')

axesTidal.set_xlabel('Time [Gyr]')
axesTidal.set_ylabel(r'Truncation, $x$ (solid); heating, $s$ (dashed)')
axesTidal.legend(fontsize='x-small')

figure.tight_layout()
figure.savefig(arguments.outputFile)
