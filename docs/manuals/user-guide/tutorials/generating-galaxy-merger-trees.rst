Generating Galaxy Merger Trees
==============================

By "merger tree" we typically mean a merger tree of dark matter halos. However, we can also look at a galaxy merger tree. Since galaxies do not merge at the same time as their halos, the galaxy merger tree will be different from the associated halo merger tree.

To track and output galaxy merger trees requires adding a ``nodeOperator`` and ``nodePropertyExtractor``:

.. code-block:: xml

   <nodeOperator value="galaxyMergerTree">
     <timeStep value="0.025"/>
     <nodePropertyExtractor value="massStellar"/>
     <nodePropertyExtractor value="massISM"/>
     <nodePropertyExtractor value="massBasic"/>
     <nodePropertyExtractor value="starFormationRate">
       <component value="total"/>
     </nodePropertyExtractor>
   </nodeOperator>

This ``nodeOperator`` records information needed to describe the merger tree of each galaxy. The ``timestep`` determines how frequently the properties of each galaxy are sampled in the tree. The ``nodePropertyExtractor``\ s here determine which properties are recorded for each galaxy. (Currently any member of the ``nodePropertyExtractorScalar`` class can be used here.)

.. code-block:: xml

   <nodePropertyExtractor value="galaxyMergerTree"/>

This ``nodePropertyExtractor`` causes the information recorded for the merger tree to be output to the output file.

As an example, you can generate a simple example by running this model:

.. code-block:: bash

   ./Galacticus.exe parameters/tutorials/galaxyMergerTree.xml

and then generate a plot of the tree from it using:

.. code-block:: bash

   ./scripts/analysis/galaxyMergerTree.py

which should look like this:

.. figure:: https://github.com/user-attachments/assets/54b728a5-bd5a-4a32-9a92-e81f95e3db78

   image
