#!/bin/sh

# Build, one at a time, the objects whose compiles need the most memory, so that none of them runs beside another.
#
# Usage: ./scripts/build/buildLargeObjects.sh <buildPath> [make options...]
#
# <buildPath> is the build directory selected by the make options (e.g. ./work/build, or ./work/buildMPI for MPI=yes). The
# make options are passed to each make, and must be those of the build which follows, since they select its configuration.
#
# Used by the CI builds, which run on machines with ~16 GB of memory. There, `make -j2` alone is not enough: the largest
# compile peaks at ~14 GiB, and the next six at 6-10 GiB, so two of them running together can exhaust memory and get the
# compiler killed (`Killed signal terminated program f951`). Each is built here as a make target of its own: make builds
# its prerequisites two at a time, then compiles it with nothing running beside it. Every other compile peaks below
# ~5.5 GiB, so any two of the rest fit together. Peak memories were measured with gfortran 16.0.1, with and without
# link-time optimization (the build profile, from `make BUILDPROFILE=yes`, reports them).
#
# The objects are listed with prerequisites first: `merger_trees/construct/read/_class` depends on all of the others, and
# `hot_halo/standard/_class` on `hydrogen`. A large object built as a prerequisite of a later one would be compiled
# beside whatever else make was running, so any object added here must come after every listed object that it depends on.

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <buildPath> [make options...]" >&2
    exit 1
fi
buildPath=${1%/}
shift

objects="
    chemical/reaction_rates/hydrogen
    objects/nodes/components/spheroid/standard/_class
    objects/nodes/components/disk/standard/_class
    objects/nodes/components/hot_halo/standard/_class
    objects/nodes/components/nuclear_star_cluster/standard/_class
    accretion/halo/simple
    merger_trees/construct/read/_class
"
# Wrappers generated only for the library build.
case " $* " in
    *" LIBRARY=yes "*)
        objects="$objects
    libgalacticus/mergerTreeImporter
    libgalacticus/nodeOperator
    libgalacticus/starFormationHistory
"
        ;;
esac

for object in $objects; do
    make -j2 "$@" "$buildPath/$object.o"
done
