#!/usr/bin/env python3
"""Execute every tutorial notebook in place, failing on any cell error.

Used by CI (the Python-Interface job) to keep the committed tutorials green
against each freshly built library, and usable locally after a rebuild:

    python3 tutorials/executeNotebooks.py

Each notebook is executed with the working directory set to `tutorials/`
(exactly how a user runs them — the setup cell resolves the library
relative to that) and rewritten in place with its fresh outputs.
"""

import glob
import os
import sys
import time

import nbformat
from nbclient import NotebookClient

TUTORIALS = os.path.dirname(os.path.abspath(__file__))


def main():
    paths = sorted(glob.glob(os.path.join(TUTORIALS, '*.ipynb')))
    if not paths:
        sys.exit("no notebooks found under " + TUTORIALS)
    failures = 0
    for path in paths:
        name     = os.path.basename(path)
        notebook = nbformat.read(path, as_version=4)
        started  = time.time()
        client   = NotebookClient(notebook, timeout=900, kernel_name='python3',
                                  resources={'metadata': {'path': TUTORIALS}})
        try:
            client.execute()
            nbformat.write(notebook, path)
            print(f"PASS {name} ({time.time()-started:.1f}s)")
        except Exception as exc:                       # noqa: BLE001 — report and continue
            failures += 1
            print(f"FAIL {name}: {type(exc).__name__}")
            print(str(exc)[:4000])
    if failures:
        sys.exit(f"{failures} notebook(s) failed")
    print("all notebooks executed successfully")


if __name__ == '__main__':
    main()
