Extracting Star Formation Histories
===================================

Star formation histories (SFHs) are stored in the output ``nodeData`` group but have a more complicated structure than most other quantities. The reason for this is that the shape of each star formation history (e.g., the number of times at which it is tabulated) may vary from one node to another, and some nodes may not have a star formation history at all (e.g., if no disk ever formed in a node, it likely will have no disk star formation history).

If you read a star formation history dataset using ``h5py``, what you get is not a ``numpy`` array directly, but actually a list with one entry per galaxy. Each entry in that list is another list (or the entry might be empty if no star formation history was recorded for that galaxy). That second-level list typically has one entry per metallicity, and each of those entries is a ``numpy`` array containing the SFH for that galaxy in that metallicity bin.

Often, what is of interest is the total SFH summed over all metallicity bins. To get this, first read the dataset, e.g. assuming you have a Galacticus output file called ``myModelFile.hdf5`` and you want the star formation histories from the first output group:

.. code-block:: python

   model   = h5py.File('myModelFile.hdf5','r')
   diskSFH = model['Outputs/Output1/nodeData/diskStarFormationHistoryMass'][:]

Then, the SFH (summed over metallicity) for galaxy ``i`` can be retrieved using:

.. code-block:: python

   np.sum(diskSFH[0],axis=0)

which will give you a ``numpy`` array of the SFH for that galaxy.

To obtain a 2-D ``numpy`` array with SFHs for all galaxies use:

.. code-block:: python

   SFHs = np.array(list(map(lambda x: x if isinstance(x, np.ndarray) else np.zeros(max(map(lambda x: x[0].size if len(x) > 0  else 0,diskSFH))),map(lambda x: np.sum(x,axis=0),diskSFH))))

The ``map`` operations here sum the SFH over metallicity bins for each galaxy, and, for galaxies with no measured SFH, inserts a zero ``numpy`` array of the correct length, and then merges these back into a single, 2-D ``numpy`` array.
