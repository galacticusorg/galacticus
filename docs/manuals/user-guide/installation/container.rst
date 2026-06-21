Running Galacticus in a Container
=================================

We provide containerized versions of Galacticus using both the GitHub
Container Registry (GHCR) and SingularityHub. These provide a complete,
containerized environment with Galacticus pre-installed and ready to run. It's
the easiest way to get started using Galacticus.

(For Docker image build instructions see `here <https://github.com/galacticusorg/galacticus/wiki/Building-Docker-Images>`_.)

Docker
------

The following instructions show how to run Galacticus using a Docker container.

* Pre-requisite: you must have a `Docker <https://www.docker.com/>`_ engine installed and running on your system

* Download the Galacticus image from the GitHub Container Registry (GHCR) repository:

  .. code-block:: bash

     docker pull ghcr.io/galacticusorg/galacticus:latest

* Start a container from the image:

  .. code-block:: bash

     docker run --rm --name galacticus -it ghcr.io/galacticusorg/galacticus:latest bash

* Once inside the container, you can find Galacticus in ``/opt/galacticus``

* Try running a simple model to check everything works:

  .. code-block:: bash

     cd /opt/galacticus
     ./Galacticus.exe parameters/quickTest.xml

* You can exit the container, but note that this will lose all of your work:

  .. code-block:: bash

     exit
