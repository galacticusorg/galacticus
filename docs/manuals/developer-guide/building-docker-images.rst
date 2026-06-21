Building Docker Images
======================

Docker images are typically built and released using GitHub Actions. However, they can also be built locally. This requires a machine with ``docker`` installed and requires root access.

Build Environment
-----------------

To build the Galacticus build environment, follow these steps:

#. Enter the ``galacticusDockerBuildEnv`` repo;
#. ``sudo docker build .`` - to build the image - keep note of the ``imageID`` output at the end of the successful build;
#. ``sudo docker tag <imageID> galacticusorg/buildenv:latest`` - to assign a tag to the image;
#. ``sudo docker push ghcr.io/galacticusorg/buildenv:latest`` - to push the new image to GitHub.

Galacticus itself
-----------------

To build the Galacticus itself, follow these steps:

#. Enter the ``galacticus`` repo;
#. ``sudo docker build --build-arg tag=latest .`` - to build the image - keep note of the ``imageID`` output at the end of the successful build;
#. ``sudo docker tag <imageID> galacticusorg/galacticus:latest`` - to assign a tag to the image;
#. ``sudo docker push ghcr.io/galacticusorg/galacticus:latest`` - to push the new image to GitHub.

Note that the ``tag`` argument in the above controls which tagged version of ``galacticusDockerBuildEnv`` is used as the base image for this build. If you wanted to build using the ``v1.0.0`` tag of ``galacticusDockerBuildEnv``, for example, you would use:

.. code-block:: bash

   sudo docker build --build-arg tag=v1.0.0 .
