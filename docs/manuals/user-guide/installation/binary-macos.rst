Installing a Pre-compiled Binary (macOS)
========================================

.. warning::

   macOS support is in beta-testing - this may or may not work for you. Please report success or failure in the `discussion forum <https://github.com/galacticusorg/galacticus/discussions>`_.

.. note::

   Pre-compiled binaries are provided for Apple Silicon (M1 and newer) only. Builds for Intel (x86-64) Macs have been retired: GitHub Actions is retiring its Intel macOS runners, and Homebrew no longer publishes pre-built packages for that platform. Releases up to and including `v0.9.12 <https://github.com/galacticusorg/galacticus/releases/tag/v0.9.12>`_ still carry a ``Galacticus_MacOS.exe`` built for Intel; for anything newer on an Intel Mac you must `build from source <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/installation/source-macos.html>`_.

This assumes that you want to install a pre-compiled Galacticus in a folder called ``Galacticus`` from your home directory. (You can also attempt to `install Galacticus from source <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/installation/source-macos.html>`_ - usually necessary only if you want to modify the code.)

#. Download and unpack the `source <https://github.com/galacticusorg/galacticus/archive/master.zip>`_ and `datasets <https://github.com/galacticusorg/datasets>`_ that are needed at run-time:

   .. code-block:: bash

      mkdir ~/Galacticus
      cd ~/Galacticus
      curl -L https://github.com/galacticusorg/galacticus/archive/master.zip --output galacticus.zip
      curl -L https://github.com/galacticusorg/datasets/archive/master.zip --output datasets.zip
      unzip galacticus.zip
      unzip datasets.zip
      mv galacticus-master galacticus
      mv datasets-master datasets

#. Download the pre-compiled binary and the tools package, move them to the ``~/Galacticus`` folder, and unpack it:

   .. code-block:: bash

      cd ~/Galacticus
      curl -L https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_MacOS-M1.exe --output Galacticus_MacOS-M1.exe
      mv Galacticus_MacOS-M1.exe galacticus/Galacticus.exe
      chmod u=wrx galacticus/Galacticus.exe
      curl -L https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/toolsMacOSM1.tar.zst --output toolsMacOSM1.tar.zst
      mkdir -p datasets
      mv toolsMacOSM1.tar.zst datasets/
      cd datasets
      tar xf toolsMacOSM1.tar.zst
      cd ..

   .. note::

      The tools archive is compressed with `zstd <https://facebook.github.io/zstd/>`_. The ``tar`` shipped with macOS understands it; if yours does not, install ``zstd`` (for example ``brew install zstd`` or ``sudo port install zstd``) and unpack with ``zstd -d toolsMacOSM1.tar.zst`` followed by ``tar xf toolsMacOSM1.tar``.

#. Set environment variables to indicate the locations at which you downloaded the source and data:

   .. code-block:: bash

      export GALACTICUS_EXEC_PATH=~/Galacticus/galacticus
      export GALACTICUS_DATA_PATH=~/Galacticus/datasets

   .. note::

      Galacticus needs to write some data files to disk at run time. Usually these are written to ``$GALACTICUS_DATA_PATH/dynamic/``. If you do not have write permission to that location, you should set the environment variable ``GALACTICUS_DYNAMIC_DATA_PATH`` to a path where dynamically-generated files can be written.

#. macOS typically won't allow you to run arbitrary unsigned executables that you download. To get around this, find the ``Galacticus.exe`` executable in finder, ctrl-click it, and select "Open". You see a message warning that the executable can't be verified - click "Open" anyway. Then ignore any messages or apps that open as a result. On newer versions of macOS (e.g. Ventura) you may also need to explicitly allow each executable in the "Privacy & Security" settings. If the executable fails to run after doing this steps described above with a message saying that the developer cannot be verified, open "System Settings", go to "Privacy & Security" - you should see a message such as '"Galacticus.exe" was blocked from use because it is not from an identified developer' - click the "Allow Anyway" button next to it and you should now be able to run the executable. You should then be able to run the Galacticus executable. Note that you will need to repeat this procedure for the following executables also:

   * ``~/Galacticus/datasets/dynamic/RecFast/recfast.exe``
   * ``~/Galacticus/datasets/dynamic/CAMB-1.3.2/fortran/camb``
   * ``~/Galacticus/datasets/dynamic/class_public-3.0.2/class``
   * ``~/Galacticus/datasets/dynamic/fsps-3.2/src/autosps.exe``
   * ``~/Galacticus/datasets/dynamic/c17.02/source/cloudy.exe``

#. You can then run a quick test model using:

   .. code-block:: bash

      cd ~/Galacticus/galacticus
      ./Galacticus.exe parameters/quickTest.xml

Debugging
---------

If you run into problems using Galacticus under macOS it can be useful to download the `debug symbols <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/debugSymbolsMacOS-M1.zip>`_ and unpack them into the same folder as your ``Galacticus.exe`` executable to allow backtrace information to be generated.
