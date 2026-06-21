# Welcome to Galacticus

[![CI/CD](https://github.com/galacticusorg/galacticus/actions/workflows/cicd.yml/badge.svg)](https://github.com/galacticusorg/galacticus/actions/workflows/cicd.yml)

Welcome to the Galacticus project. Galacticus is a semi-analytic model of galaxy formation - a powerful toolkit for modeling the physics of how galaxies form.
For more information please see the [wiki](https://github.com/galacticusorg/galacticus/wiki), and for a description of the physics see the [description paper](https://arxiv.org/abs/1008.1786) on Galacticus.

Have questions? Ask them in the [discussion forum](https://github.com/galacticusorg/galacticus/discussions).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=master&repo=204087682)

---

## Quickstart

This section walks you through building and running a minimal Galacticus model for the first time.

### Prerequisites

Before building, make sure you have the following available:

- A modern Fortran compiler (`gfortran` ≥ 16; earlier versions will not compile Galacticus)
- `make`
- GSL, zlib, HDF5, FoX, and BLAS libraries (HDF5 and FoX must be built with the same compiler used for Galacticus)
- Python 3 (≥ 3.9; used by the build system for code generation and various supporting scripts)

FFTW3, ANN, and `libmatheval` are optional (needed only for specialized calculations). For the full, platform-specific prerequisite list and step-by-step instructions, see the [Installation guide](https://galacticus.readthedocs.io/en/latest/manuals/user-guide/installation/index.html).

> **Tip:** The easiest way to get a fully configured environment is to use [GitHub Codespaces](#open-in-github-codespaces) (click the badge above). All dependencies are pre-installed.

### Installing Python dependencies

The build system imports modules from `python/`. Make those modules importable, and install the third-party packages they depend on (numpy, scipy, h5py, lxml, matplotlib, …), with a single editable install from the repo root:

```bash
pip install -e .
```

This is sufficient for building Galacticus and running its supporting scripts. Heavy emulation pipelines have additional dependencies (TensorFlow, lenstronomy, etc.); install those on demand with `pip install -e '.[emulation]'`. To run the Python test suite, install with `pip install -e '.[test]'`.

### Building Galacticus

From the root of the repository, run:

```bash
make -jN Galacticus.exe
```

Replace `N` with the number of CPU cores you want to use in parallel (e.g., `make -j4 Galacticus.exe` to use 4 cores). Using multiple cores significantly speeds up the build. If you are unsure how many cores are available, you can run `nproc` on Linux to find out.

### Running a minimal model

Once the build completes successfully, run the included quick-test parameter file:

```bash
./Galacticus.exe parameters/quickTest.xml
```

This runs a small, pre-configured galaxy formation model designed to complete quickly.

### What to expect

A successful run will:

- Exit with code `0` (no error message printed to the terminal).
- Write output to an HDF5 file in the working directory. `parameters/quickTest.xml` does not set an output filename, so Galacticus uses the default, `galacticus.hdf5`. Open that file to inspect the results (e.g., with `h5ls` or any HDF5 viewer).
- Print progress information to standard output during the run. The final line should indicate that the run completed without errors.

The output is written to `galacticus.hdf5` in the directory you ran from. To change this, set an `outputFileName` parameter in the parameter file.

### Troubleshooting

| Problem | What to try |
|---|---|
| `make: command not found` | Install `make` (e.g., `sudo apt install make` on Debian/Ubuntu). |
| Compiler errors during build | Check that a supported Fortran compiler is installed and on your `PATH`. |
| Build hangs or fails with parallel jobs | Try removing the `-jN` flag and running `make Galacticus.exe` (single-threaded) to get clearer error output. |
| `./Galacticus.exe: No such file or directory` | The build did not complete successfully. Re-run `make` and check for errors. |
| Parameter file not found | Ensure you are running the command from the repository root directory so that `parameters/quickTest.xml` resolves correctly. |
| Missing library errors at link time | Verify that HDF5, FFTW3, and GSL development packages are installed and that their locations are on the relevant library paths. |

For further help, visit the [wiki](https://github.com/galacticusorg/galacticus/wiki) or ask in the [discussion forum](https://github.com/galacticusorg/galacticus/discussions).

## Analyzing Galacticus output

Analysis and visualization of Galacticus outputs - including plotting of on-the-fly analyses, MCMC chain diagnostics, and posterior corner plots - is provided by the [Dendros](https://github.com/galacticusorg/dendros) package, available on [PyPI](https://pypi.org/project/dendros/):

```bash
python3 -m pip install dendros
```
