#!/usr/bin/env python3
import sys
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Create plots of collected metadata on merger tree ODE evolver.
# Andrew Benson (05-June-2011) [Perl]; ported to Python

if len(sys.argv) != 3:
    print("Usage: metaAnalysis.py <modelFile> <outputFolder>", file=sys.stderr)
    sys.exit(1)

modelFile    = sys.argv[1]
outputFolder = sys.argv[2]

# Create the directory.
os.mkdirs(outputFolder, exist_ok=True)

# Open the file containing the meta-data.
with h5py.File(modelFile, 'r') as HDFfile:
    if 'metaData' not in HDFfile:
        print("metaAnalysis.py: metaData group does not exist", file=sys.stderr)
        sys.exit(1)
    if 'evolverProfiler' not in HDFfile['metaData']:
        print("metaAnalysis.py: metaData/evolverProfiler group does not exist", file=sys.stderr)
        sys.exit(1)

    grp = HDFfile['metaData/evolverProfiler']

    timeSteps                  = grp['timeStep'                  ][()]
    timeStepCount              = grp['timeStepCount'             ][()]
    evaluationCount            = grp['evaluationCount'           ][()]
    timeCPU                    = grp['timeCPU'                   ][()]
    timeStepCountInterrupted   = grp['timeStepCountInterrupted'  ][()]
    evaluationCountInterrupted = grp['evaluationCountInterrupted'][()]
    timeCPUInterrupted         = grp['timeCPUInterrupted'        ][()]
    propertyNames_raw          = grp['propertyNames'             ][()]
    propertyHitCount           = grp['propertyHitCount'          ][()]

# Decode property names.
if propertyNames_raw.dtype.kind in ('S', 'O'):
    propertyNames = [n.decode() if isinstance(n, bytes) else str(n) for n in propertyNames_raw]
else:
    propertyNames = [str(n) for n in propertyNames_raw]

propertyHitRate  = 100.0 * propertyHitCount / propertyHitCount.sum()
propertyIndex    = np.arange(1, len(propertyHitRate) + 1)

# Convert counts into frequencies.
timeStepFrequency                    = timeStepCount.astype(float)            / timeStepCount.sum()
timeStepWeightedFrequency            = evaluationCount.astype(float)          / evaluationCount.sum()
timeCPUFrequency                     = timeCPU.astype(float)                  / timeCPU.sum()
timeStepInterruptedFrequency         = timeStepCountInterrupted.astype(float) / timeStepCountInterrupted.sum()
timeStepInterruptedWeightedFrequency = evaluationCountInterrupted.astype(float) / evaluationCountInterrupted.sum()
timeCPUInterruptedFrequency          = timeCPUInterrupted.astype(float)       / timeCPUInterrupted.sum()

# Find cumulative frequencies.
timeStepFrequencyCumulative                    = np.cumsum(timeStepFrequency)                    / timeStepFrequency.sum()
timeStepWeightedFrequencyCumulative            = np.cumsum(timeStepWeightedFrequency)            / timeStepWeightedFrequency.sum()
timeCPUFrequencyCumulative                     = np.cumsum(timeCPUFrequency)                     / timeCPUFrequency.sum()
timeStepInterruptedFrequencyCumulative         = np.cumsum(timeStepInterruptedFrequency)         / timeStepInterruptedFrequency.sum()
timeStepInterruptedWeightedFrequencyCumulative = np.cumsum(timeStepInterruptedWeightedFrequency) / timeStepInterruptedWeightedFrequency.sum()
timeCPUInterruptedFrequencyCumulative          = np.cumsum(timeCPUInterruptedFrequency)          / timeCPUInterruptedFrequency.sum()

# Find ranges for the timestep plot.
allFrequencies  = np.concatenate([timeStepFrequency, timeStepWeightedFrequency,
                                  timeCPUFrequency, timeStepInterruptedFrequency, 
                                  timeStepInterruptedWeightedFrequency, timeCPUInterruptedFrequency])
nonZero         = allFrequencies > 0.0
timestepMinimum = 10.0 ** np.floor(np.log10(timeSteps[ 0]))
timestepMaximum = 10.0 ** np.ceil (np.log10(timeSteps[-1]))
frequencyMinimum = 10.0 ** np.floor(np.log10(allFrequencies[nonZero].min()))
frequencyMaximum = 10.0 ** np.ceil (np.log10(allFrequencies[nonZero].max()))

cornflowerBlue  = '#6495ED'
mediumSeaGreen  = '#3CB371'
indianRed       = '#CD5C5C'

# --- Plot 1: Timestep frequency distribution ---
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(timeSteps, timeStepFrequency,                    'o', color=cornflowerBlue, markersize=4, label='Steps')
ax.plot(timeSteps, timeStepWeightedFrequency,            'o', color=mediumSeaGreen, markersize=4, label='Evaluations')
ax.plot(timeSteps, timeCPUFrequency,                     'o', color=indianRed,      markersize=4, label='CPU')
ax.plot(timeSteps, timeStepInterruptedFrequency,         's', color=cornflowerBlue, markersize=4, label='Steps (interrupted)')
ax.plot(timeSteps, timeStepInterruptedWeightedFrequency, 's', color=mediumSeaGreen, markersize=4, label='Evaluations (interrupted)')
ax.plot(timeSteps, timeCPUInterruptedFrequency,          's', color=indianRed,      markersize=4, label='CPU (interrupted)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(timestepMinimum, timestepMaximum)
ax.set_ylim(frequencyMinimum, frequencyMaximum)
ax.set_title('Timestep distribution')
ax.set_xlabel('Timestep [Gyr]')
ax.set_ylabel('Frequency')
ax.legend(loc='lower left', fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(outputFolder, 'timesteps.pdf'))
plt.close()

# --- Plot 2: Cumulative timestep frequency distribution ---
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(timeSteps, timeStepFrequencyCumulative,                    'o', color=cornflowerBlue, markersize=4, label='Steps')
ax.plot(timeSteps, timeStepWeightedFrequencyCumulative,            'o', color=mediumSeaGreen, markersize=4, label='Evaluations')
ax.plot(timeSteps, timeCPUFrequencyCumulative,                     'o', color=indianRed,      markersize=4, label='CPU')
ax.plot(timeSteps, timeStepInterruptedFrequencyCumulative,         's', color=cornflowerBlue, markersize=4, label='Steps (interrupted)')
ax.plot(timeSteps, timeStepInterruptedWeightedFrequencyCumulative, 's', color=mediumSeaGreen, markersize=4, label='Evaluations (interrupted)')
ax.plot(timeSteps, timeCPUInterruptedFrequencyCumulative,          's', color=indianRed,      markersize=4, label='CPU (interrupted)')
ax.set_xscale('log')
ax.set_xlim(timestepMinimum, timestepMaximum)
ax.set_ylim(-0.02, 1.02)
ax.set_title('Cumulative timestep distribution')
ax.set_xlabel('Timestep [Gyr]')
ax.set_ylabel('Cumulative fraction')
ax.legend(loc='lower left', fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(outputFolder, 'timestepsCumulative.pdf'))
plt.close()

# --- Plot 3: Property hit rates ---
fig, ax = plt.subplots(figsize=(6, 4))
ax.bar(propertyIndex, propertyHitRate, color=cornflowerBlue)
for iProperty, name in enumerate(propertyNames):
    ax.text(propertyIndex[iProperty], 0, name, rotation=90, va='bottom', ha='center', fontsize=7)
ax.set_xlim(0, len(propertyHitCount) + 1)
ax.set_ylim(0, propertyHitRate.max() * 1.05)
ax.set_title('Hit Frequency')
ax.set_xlabel('Property')
ax.set_ylabel('Hit frequency [%]')
plt.tight_layout()
plt.savefig(os.path.join(outputFolder, 'hitRates.pdf'))
plt.close()
