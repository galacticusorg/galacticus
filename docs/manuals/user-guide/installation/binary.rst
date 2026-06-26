Installing a Pre-compiled Binary (Linux)
========================================

By far the easiest way to install and use Galacticus on Linux is to use a pre-compiled binary. For macOS, see `the macOS binary instructions <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/installation/binary-macos.html>`_.

To do this:

* download and unpack the `source <https://github.com/galacticusorg/galacticus/archive/master.zip>`_ and `datasets <https://github.com/galacticusorg/datasets>`_ that are needed at run-time:

.. code-block:: bash

   wget https://github.com/galacticusorg/galacticus/archive/master.zip -O galacticus.zip
   wget https://github.com/galacticusorg/datasets/archive/master.zip -O datasets.zip
   unzip galacticus.zip
   unzip datasets.zip
   mv galacticus-master galacticus
   mv datasets-master datasets

* download and unpack the `tools <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/tools.tar.bz2>`_ that are needed at run-time:

.. code-block:: bash

   cd datasets
   wget https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/tools.tar.bz2
   tar xvfj tools.tar.bz2

* download the `pre-compiled binary <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus.exe>`_, put it inside the ``galacticus`` directory and give it executable permissions:

.. code-block:: bash

   cd ../galacticus
   wget https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus.exe
   chmod u=wrx Galacticus.exe

* set environment variables to indicate the locations at which you downloaded the source and data:

.. code-block:: bash

   export GALACTICUS_EXEC_PATH=/path/to/galacticus
   export GALACTICUS_DATA_PATH=/path/to/datasets

.. note::

   Galacticus needs to write some data files to disk at run time. Usually these are written to ``$GALACTICUS_DATA_PATH/dynamic/``. If you do not have write permission to that location, you should set the environment variable ``GALACTICUS_DYNAMIC_DATA_PATH`` to a path where dynamically-generated files can be written.

.. note::

   The run-time tools (CAMB, CLASS, Cloudy, FSPS, RecFast, AxionCAMB, and mangle) are, by default, stored under the dynamic data path alongside regenerable data. You can relocate them by setting the environment variable ``GALACTICUS_TOOLS_PATH`` to a separate directory. This is useful if you want to keep the (immutable) pre-built tools apart from regenerable cache data — for example so that the cache can be cleared without removing the tools, or to share a read-only tools directory between installs. When unset, ``GALACTICUS_TOOLS_PATH`` defaults to the dynamic data path, so existing setups are unaffected.

You can then run a quick test model using:

.. code-block:: bash

   ./Galacticus.exe parameters/quickTest.xml

Notes for WSL
-------------

We have successfully run Galacticus on Windows using WSL (Windows Subsystem for Linux), version 2, but found that version 1 did not allow Galacticus to run successfully. Therefore, to run Galacticus on WSL you should check that you are running version 2 of WSL and, if necessary, upgrade to version 2. Details on how to do this are available `here <https://learn.microsoft.com/en-us/windows/wsl/install#upgrade-version-from-wsl-1-to-wsl-2>`_.
