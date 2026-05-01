# Provides functions for working with Cloudy
# Andrew Benson (22-April-2025)
"""Functions for downloading, building, and managing the Cloudy photoionization code.

This module provides utilities for obtaining and compiling the
`Cloudy <https://trac.nublado.org/>`_ spectral synthesis and photoionization
code, which is used by Galacticus to compute cooling rates and other
thermodynamic quantities.
"""
import os
import re
import subprocess
import urllib.request
from urllib.error import HTTPError

def initialize(options):
    """Download, unpack, and compile the Cloudy photoionization code.

    Fetches the requested version of Cloudy from the official distribution
    server, unpacks the tarball, and builds the ``cloudy.exe`` binary.  If any
    of these steps have already been completed, they are skipped so the
    function is safe to call multiple times.

    The Cloudy version is determined in the following order of precedence:

    1. ``options['version']`` if already set by the caller.
    2. The version listed in ``$GALACTICUS_EXEC_PATH/aux/dependencies.yml``,
       when the ``GALACTICUS_EXEC_PATH`` environment variable is set.
    3. The hard-coded fallback ``"23.01"``.

    Parameters
    ----------
    options : dict
        Configuration dictionary.  The following keys are recognised:

        ``'version'`` : str, optional
            Cloudy release string, e.g. ``"23.01"``.  Populated automatically
            from ``dependencies.yml`` when not provided.

    Returns
    -------
    tuple[str, str]
        A 2-tuple ``(cloudy_path, cloudy_version)`` where *cloudy_path* is the
        absolute path to the directory containing ``cloudy.exe`` (with a
        trailing slash) and *cloudy_version* is the version string that was
        used.

    Notes
    -----
    The following environment variables are used:

    ``GALACTICUS_DATA_PATH``
        Root path under which Cloudy is installed (in a ``dynamic/`` sub-
        directory).  Must be set.
    ``GALACTICUS_EXEC_PATH``
        Path to the Galacticus executable tree.  Optional; used to locate
        ``aux/dependencies.yml``.
    ``CLOUDY_COMPILER_PATH``
        If set, prepended to ``PATH`` before invoking the Cloudy build system.
    ``CLOUDY_STATIC_BUILD``
        Set to ``"yes"`` to request a statically linked Cloudy binary.

    Raises
    ------
    RuntimeError
        If any of the download, unpack, or build steps fail.
    """
    # Initialize Cloudy by downloading and compiling.
    # Determine Cloudy version.
    if os.environ.get('GALACTICUS_EXEC_PATH') is not None and 'version' not in options:
        with open(os.environ.get('GALACTICUS_EXEC_PATH')+"/aux/dependencies.yml") as dependencies:
            for line in dependencies:
                match = re.match(r'^cloudy:\s+([\d\.]+)',line)
                if match:
                    options['version'] = match.group(1)
    cloudyVersion      = options.get('version', "23.01")
    cloudyVersionMajor = cloudyVersion[0:cloudyVersion.index(".")]
    # Specify Cloudy path.
    cloudyPath = os.environ['GALACTICUS_DATA_PATH']+"/dynamic/c"+cloudyVersion
    os.makedirs(cloudyPath, exist_ok=True)
    # Download the code.
    if not os.path.isfile(f"{cloudyPath}.tar.gz"):
        print("downloading Cloudy code...")
        try:
            urllib.request.urlretrieve(f"http://data.nublado.org/cloudy_releases/c{cloudyVersionMajor}/c{cloudyVersion}.tar.gz", f"{cloudyPath}.tar.gz")
        except HTTPError as e:
            raise RuntimeError(f"failed to download Cloudy: HTTP Error: {e.code} - {e.reason}") from e
        except Exception as e:
            raise RuntimeError(f"failed to download Cloudy: {e}") from e
    # Unpack the code.
    if not os.path.isfile(f"{cloudyPath}/source/Makefile"):
        print("unpacking Cloudy code...")
        unpack = subprocess.run(['tar', '-x', '-v', '-z', '-C', os.environ['GALACTICUS_DATA_PATH']+'/dynamic', '-f', f"{cloudyPath}.tar.gz"])
        if unpack.returncode != 0:
            raise RuntimeError("failed to unpack Cloudy")
    if not os.path.isfile(f"{cloudyPath}/source/Makefile"):
        raise RuntimeError("failed to unpack Cloudy code.")
    # Build the code.
    if not os.path.isfile(f"{cloudyPath}/source/cloudy.exe"):
        print("compiling Cloudy code...")
        buildCommand      = f"cd {cloudyPath}/source; chmod u=wrx configure.sh capabilities.pl; "
        buildCommand     += "" if os.environ.get('CLOUDY_COMPILER_PATH') is None else f"export PATH={os.environ.get('CLOUDY_COMPILER_PATH')}:$PATH; "
        buildCommand     += " make -f Makefile_modified"
        extra             = "EXTRA = -static\n" if os.environ.get('CLOUDY_STATIC_BUILD') == "yes" else "EXTRA = \n"
        with open(f"{cloudyPath}/source/Makefile", 'r') as makeFile, open(f"{cloudyPath}/source/Makefile_modified", 'w') as makeFileModified:
            for line in makeFile:
                lineModified = re.sub(r'^EXTRA\s*=.*', extra, line)
                makeFileModified.write(lineModified)
        build = subprocess.run(buildCommand, shell=True)
        if build.returncode != 0:
            raise RuntimeError("failed to build Cloudy code")
    # Return path and version.
    return (f"{cloudyPath}/", cloudyVersion)
