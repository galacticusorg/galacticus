Installing a Pre-compiled Binary (macOS)
========================================

.. warning::

   macOS support is in beta-testing - this may or may not work for you. Please report success or failure in the `discussion forum <https://github.com/galacticusorg/galacticus/discussions>`_.

This assumes that you want to install a pre-compiled Galacticus in a folder called ``Galacticus`` from your home directory. (You can also attempt to `install Galacticus from source <https://github.com/galacticusorg/galacticus/wiki/Installation-from-source-on-MacOS>`_ - usually necessary only if you want to modify the code.)

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

   * For x86 chips (older Macs):

     .. code-block:: bash

        cd ~/Galacticus
        curl -L https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_MacOS.exe --output Galacticus_MacOS.exe
        mv Galacticus_MacOS.exe galacticus/Galacticus.exe
        chmod u=wrx galacticus/Galacticus.exe
        curl -L https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/toolsMacOS.zip --output toolsMacOS.zip
        mkdir datasets
        mv toolsMacOS.zip datasets/
        cd datasets
        unzip toolsMacOS.zip
        cd ..

   * For Apple Silicon (M1 and newer chips):

     .. code-block:: bash

        cd ~/Galacticus
        curl -L https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_MacOS.exe --output Galacticus_MacOS-M1.exe
        mv Galacticus_MacOS-M1.exe galacticus/Galacticus.exe
        chmod u=wrx galacticus/Galacticus.exe
        curl -L https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/toolsMacOS.zip --output toolsMacOS-M1.zip
        mkdir datasets
        mv toolsMacOS-M1.zip datasets/
        cd datasets
        unzip toolsMacOS-M1.zip
        cd ..

#. Set environment variables to indicate the locations at which you downloaded the source and data:

   .. code-block:: bash

      export GALACTICUS_EXEC_PATH=~/Galacticus/galacticus
      export GALACTICUS_DATA_PATH=~/Galacticus/datasets

   .. note::

      Galacticus needs to write some data files to disk at run time. Usually these are written to ``$GALACTICUS_DATA_PATH/dynamic/``. If you do not have write permission to that location, you should set the environment variable ``GALACTICUS_DYNAMIC_DATA_PATH`` to a path where dynamically-generated files can be written.

#. macOS typically won't allow you to run arbitrary unsigned executables that you download. To get around this, find the ``Galacticus.exe`` executable in finder, ctrl-click it, and select "Open". You see a message warning that the executable can't be verified - click "Open" anyway. Then ignore any messages or apps that open as a result. On newer versions of macOS (e.g. Ventura) you may also need to explicitly allow each executable in the "Privacy & Security" settings. If the executable fails to run after doing this steps described above with a message saying that the developer cannot be verified, open "System Settings", go to "Privacy & Security" - you should see a message such as '"Galacticus.exe" was blocked from use because it is not from an identified developer' - click the "Allow Anyway" button next to it and you should now be able to run the executable. You should then be able to run the Galacticus executable. Note that you will need to repeat this procedure for the following executables also:

   * ``~/datasets/dynamic/RecFast/recfast.exe``
   * ``~/datasets/dynamic/CAMB-1.3.2/fortran/camb``
   * ``~/datasets/dynamic/class_public-3.0.2/class``
   * ``~/datasets/dynamic/fsps-3.2/src/autosps.exe``
   * ``~/datasets/dynamic/c17.02/source/cloudy.exe``

#. You can then run a quick test model using:

   .. code-block:: bash

      cd ~/Galacticus/galacticus
      ./Galacticus.exe parameters/quickTest.xml

Debugging
---------

If you run into problems using Galacticus under macOS it can be useful to download the `debug symbols <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/debugSymbolsMacOS.zip>`_ (or these `debug symbols <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/debugSymbolsMacOS-M1.zip>`_ for Apple Silicon chips) and unpack them into the same folder as your ``Galacticus.exe`` executable to allow backtrace information to be generated.
