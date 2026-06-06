#!/usr/bin/env python3
import sys
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Make a plot of the specified cooling function file.
# Andrew Benson (27-Jan-2009) [Perl]; ported to Python

if len(sys.argv) < 2 or len(sys.argv) > 3:
    print("Usage: cooling_function.py <coolingFunctionFile> [<pdfFile>]", file=sys.stderr)
    sys.exit(1)

coolingFunctionFile = sys.argv[1]
pdfFile             = sys.argv[2] if len(sys.argv) == 3 else "cooling_function.pdf"

# Read the data file.
with h5py.File(coolingFunctionFile, 'r') as data:
    description     = data.attrs['description']
    if isinstance(description, bytes):
        description = description.decode()
    coolingFunctions = data['coolingRate'][()]
    metallicities    = data['metallicity'][()]
    temperatures     = data['temperature'][()]

# Build plot titles and find metallicity range (excluding primordial).
lZMin = +1.0e100
lZMax = -1.0e100
for iMetallicity, metallicity in enumerate(metallicities):
    if metallicity > -999.0:
        if metallicity > lZMax:
            lZMax = metallicity
        if metallicity < lZMin:
            lZMin = metallicity
useColorBar = lZMin < lZMax
            
# Make the plot.
fig, ax = plt.subplots(figsize=(5.0, 3.3))

nMetallicities = len(metallicities)
if useColorBar:
    cmap = plt.get_cmap('rainbow')
    norm = mcolors.Normalize(vmin=lZMin, vmax=lZMax)

for iMetallicity in range(nMetallicities):
    metallicity = metallicities[iMetallicity]
    if useColorBar:
        if metallicity <= -999.0:
            color = cmap(0.0)
        else:
            color = cmap(norm(metallicity))
    else:
        color = 'red'
    ax.plot(temperatures, coolingFunctions[:, iMetallicity], color=color, linewidth=1.5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(316.0, 1.0e9)
ax.set_ylim(3.0e-30, 1.0e-19)
ax.set_xlabel('Temperature [K]')
ax.set_ylabel(r'$\Lambda(T)$ [erg cm$^3$ s$^{-1}$]')
ax.set_title(description)

if useColorBar:
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label=r'$\log_{10}(Z/Z_\odot)$')

plt.tight_layout()
plt.savefig(pdfFile)
plt.close()
