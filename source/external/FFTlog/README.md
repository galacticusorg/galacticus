# FFTlog

This directory contains a vendored copy of the FFTLog library, used by Galacticus
for fast Hankel (and related) transforms (see `source/numerical/FFTlog.F90`).

## Provenance

* Upstream repository: https://github.com/emsig/fftlog
* Version: v0.2.2 (the `main` branch as of 2024-12-31; "Update home and build backend")
* The original FFTLog algorithm and Fortran code are by Andrew J. S. Hamilton
  (https://jila.colorado.edu/~ajsh/FFTLog/), repackaged by the emsig project.

Previously these source files were downloaded from GitHub at build time. They are
now shipped directly with Galacticus so that builds do not depend on GitHub being
reachable. This is permitted by the upstream license.

## License

FFTLog is released under the Creative Commons CC0 1.0 Universal public-domain
dedication. The full text is in the `LICENSE` file in this directory.

## Files

The following Fortran source files are taken verbatim from the upstream `src`
directory, except where noted below:

* `cdgamma.f`  — complex gamma function (unmodified).
* `fftlog.f`   — core FFTLog routines (unmodified).
* `drfftb.f`   — real FFT backward transform (PATCHED — see below).
* `drfftf.f`   — real FFT forward transform (PATCHED — see below).
* `drffti.f`   — real FFT initialization (PATCHED — see below).

## Patches

The files `drfftb.f`, `drfftf.f`, and `drffti.f` have been modified relative to
upstream. The patches as applied are kept alongside the source in this directory:

* `drfftb.f.patch`
* `drfftf.f.patch`
* `drffti.f.patch`

These patches change the `ifac` integer work array passed between the FFT
routines into a `real*8` array (`rfac`), using `transfer()` to convert between the
integer and real representations. This is required because FFTLog's `wsave` work
array is declared `real*8`, and storing the integer factorization table directly
in it triggers type-mismatch errors with modern Fortran compilers. The vendored
`.f` files already have these patches applied.

If upstream is updated, re-fetch the pristine `src/*.f` files and re-apply these
patches (e.g. `patch < drfftb.f.patch`), resolving any conflicts, then refresh the
vendored copies here.
