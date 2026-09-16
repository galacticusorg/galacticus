#!/usr/bin/env python3
"""Run the dust absorption and emission tutorial, and execute its notebook.

The tutorial parameter file, parameters/tutorials/dustEmission.xml, is run, and the notebook which analyzes its output,
tutorials/models/01-dust-absorption-and-emission.ipynb, is then executed end to end against that output, so that the
parameter file keeps running and the notebook keeps working as Galacticus changes. The notebook is executed in a
temporary copy, so the committed notebook, with its committed output, is not rewritten.

The notebook reads the output with the dendros package, which is installed with pip if it is not already available.

Andrew Benson, with assistance from Claude.
"""
import importlib
import os
import shutil
import subprocess
import sys
import tempfile

# This script runs with the working directory set to testSuite/, while the model is run from the repository root.
root          = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
parameterFile = "parameters/tutorials/dustEmission.xml"
outputFile    = os.path.join(root, "dustEmissionTutorial.hdf5")
notebookPath  = os.path.join(root, "tutorials", "models", "01-dust-absorption-and-emission.ipynb")

# Run the model.
if os.path.exists(outputFile):
    os.remove(outputFile)
status = subprocess.run(f"cd {root}; ./Galacticus.exe {parameterFile}", shell=True, capture_output=True)
log    = status.stdout.decode() + status.stderr.decode()
if status.returncode != 0 or not os.path.exists(outputFile):
    print("FAILED: the dust emission tutorial model did not run")
    print(log[-4000:])
    sys.exit(0)
print("SUCCESS: the dust emission tutorial model ran")

# Make sure the packages needed to execute the notebook are available. The notebook needs dendros 0.8.0 or later (for
# `ndim` on datasets, and compound attributes decoded to dictionaries), so an older installation is upgraded.
def versionAtLeast(version, minimum):
    """Return true if the dotted release number `version` is at least `minimum`."""
    parts = [int(part) for part in version.split("+")[0].split(".")[:3] if part.isdigit()]
    return tuple(parts + [0] * (3 - len(parts))) >= minimum

for package, installName, minimum in (("nbformat", "nbformat", None), ("nbclient", "nbclient", None), ("ipykernel", "ipykernel", None), ("dendros", "dendros[plot]>=0.8.0", (0, 8, 0))):
    try:
        module       = importlib.import_module(package)
        needsInstall = minimum is not None and not versionAtLeast(getattr(module, "__version__", "0"), minimum)
    except ImportError:
        needsInstall = True
    if needsInstall:
        install = subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", "--upgrade", "--break-system-packages", installName], capture_output=True)
        if install.returncode != 0:
            print(f"SKIPPED: the '{package}' package, needed to execute the tutorial notebook, is not available and could not be installed")
            print((install.stdout.decode() + install.stderr.decode())[-2000:])
            sys.exit(0)
import nbformat                         # noqa: E402 -- imported only once known to be available
from nbclient import NotebookClient     # noqa: E402

# Execute a copy of the notebook, from its own directory, so that it finds the model output exactly as a user would.
notebook = nbformat.read(notebookPath, as_version=4)
with tempfile.TemporaryDirectory() as scratch:
    shutil.copy2(notebookPath, scratch)
    client = NotebookClient(notebook, timeout=900, kernel_name="python3", resources={"metadata": {"path": os.path.dirname(notebookPath)}})
    try:
        client.execute()
    except Exception as exception:  # noqa: BLE001 -- report any failure of the notebook
        print(f"FAILED: the dust emission tutorial notebook did not execute: {type(exception).__name__}")
        # Report the end of the message as well as its start: the cause of a failure is at the end of a traceback, and a
        # notebook traceback is easily long enough to be truncated away.
        message = str(exception)
        if len(message) > 6000:
            message = message[:2000] + "\n[...]\n" + message[-4000:]
        print(message)
        sys.exit(0)
print("SUCCESS: the dust emission tutorial notebook executed")
sys.exit(0)
