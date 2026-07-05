# Galacticus Python library tutorials

Jupyter notebooks introducing the Galacticus Python library interface: the
physics engine of the [Galacticus](https://github.com/galacticusorg/galacticus)
semi-analytic galaxy formation model, exposed to Python through the
`libgalacticus` shared library. Each notebook is committed with executed
output, so you can read the results without running anything.

| Notebook | What it covers |
|---|---|
| [01 — Cosmology calculator](01-cosmology-calculator.ipynb) | Parameters, expansion history, ages, densities, and distances; the object-construction pattern used throughout. |
| [02 — Matter power spectrum](02-matter-power-spectrum.ipynb) | Primordial spectrum → transfer function → growth factor → $P(k,z)$ and $\sigma(M)$; warm dark matter free-streaming and half-mode masses. |
| [03 — Halo mass function](03-halo-mass-function.ipynb) | Spherical-collapse thresholds $\delta_\mathrm{c}(z)$ and $\Delta_\mathrm{vir}(z)$; Sheth–Tormen vs. Tinker et al. mass functions across cosmic time; collapsed mass fractions. |
| [04 — Python callbacks](04-python-callbacks.ipynb) | Galacticus calling *your* Python functions: quadrature of Python integrands over computational domains, and Python cross-sections in radiation fields. |
| [05 — Merger-tree building blocks](05-merger-tree-building-blocks.ipynb) | The characteristic collapsing mass $M_*(z)$, Parkinson–Cole–Helly branching rates, tree-build costs, and sampling root halo masses for tree suites. |
| [06 — Stellar populations](06-stellar-populations.ipynb) | Initial mass functions, stellar lifetimes/remnants/yields from the standard compilation, and Type Ia supernova delay-time distributions. Requires the datasets repository (`GALACTICUS_DATA_PATH`). |

## Getting the library

Two options:

1. **Build from source** (Linux/macOS, requires the Galacticus build
   dependencies):

   ```sh
   make GALACTICUS_BUILD_OPTION=lib libgalacticus.so
   mkdir -p galacticus/lib && cp libgalacticus.so galacticus/lib/
   ```

   which also generates the `galacticus.py` interface module at the
   repository root. Run the notebooks from this `tutorials/` directory and
   the setup cell finds everything automatically.

2. **Download a binary distribution** — CI publishes
   `libgalacticus.tar.bz2` with each
   [bleeding-edge release](https://github.com/galacticusorg/galacticus/releases/tag/bleeding-edge),
   containing a `galacticus/` folder with `lib/libgalacticus.so` and
   `python/galacticus.py`. Unpack it anywhere and set

   ```sh
   export GALACTICUS_LIBRARY_PATH=/path/to/directory/containing/galacticus
   ```

   before starting Jupyter.

Python requirements: `numpy` and `matplotlib` (plus `jupyter` to run the
notebooks).

## Conventions

Galacticus works in masses of $M_\odot$, times in Gyr, distances in Mpc,
wavenumbers in comoving Mpc$^{-1}$, $H_0$ in km/s/Mpc, and wavelengths in
Angstroms. Objects are constructed from parameter values and other objects
(dependency injection), mirroring the object graph Galacticus builds
internally from a parameter file.
