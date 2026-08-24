Merger Tree Structure in the Output
===================================

Each row of the datasets in the ``Outputs/Output``\ *n*\ ``/nodeData`` group of a Galacticus output file describes a single node---a dark matter halo, or a subhalo, together with the galaxy that it contains---at the time of that output. The datasets whose names end in "``Index``" describe how those rows are related to each other: which node hosts which, which node a halo will merge into, and which nodes belong to the same merger tree.

The one thing to know before using them is that the meaning of ``parentIndex``, ``siblingIndex``, and ``satelliteIndex`` **depends on whether the node is an isolated halo or a subhalo**. Galacticus uses a single set of links to represent both the progenitor structure of the tree (which halo descends into which) and the substructure hierarchy (which subhalo lives inside which halo), so the same dataset carries two different meanings depending on the row. The ``nodeIsIsolated`` dataset tells you which case you are looking at. This page explains the conventions, and gives recipes for the operations that most analyses need.

The datasets described here are produced by the :galacticus-class:`nodePropertyExtractorNodeIndices` property extractor, which is included in the output by default. The additional, unambiguously-named indices described in :ref:`sec-output-tree-structure-reference` must be requested explicitly.

Isolated halos and subhalos
---------------------------

A node is an **isolated halo** if it is not contained within any other halo. It is a **subhalo** (equivalently, a satellite) if it has fallen into, and now orbits within, another halo. The distinction is recorded in ``nodeIsIsolated``, which is :math:`1` for an isolated halo and :math:`0` for a subhalo. Galaxies in isolated halos are usually referred to as *central* galaxies, and those in subhalos as *satellite* galaxies.

Subhalos may themselves contain subhalos: a node whose ``nodeIsIsolated`` is :math:`0` can still have a non-negative ``satelliteIndex``. Such sub-subhalos arise when a halo which already contained substructure falls into a larger halo.

Links between isolated halos
----------------------------

For isolated halos the links describe the *progenitor structure* of the merger tree---which halo grows into which, as time advances.

.. mermaid::
   :config: {"block": {"padding": 40}}
   :name: fig-outputTreeStructureIsolated
   :align: center
   :caption: Links between isolated halos. Time increases upward: halos :math:`4` and :math:`5` at time :math:`t_0` merge to form halo :math:`2` at :math:`t_1`, which merges with halo :math:`3` to form halo :math:`1` at :math:`t_2`. Each arrow points in the direction of the index which labels it. Every link which exists here is drawn: halos :math:`5` and :math:`3` have ``siblingIndex`` :math:`=-1`, since each is the last in its list, and halo :math:`1` has ``parentIndex`` :math:`=-1`, since it is at the base of the tree. Note that this is the structure of the tree as it was built, spanning several timesteps---the nodes present in any one output all exist at a single time (see :ref:`sec-output-tree-structure-caveats`).

   block-beta
     columns 8
     t2["t₂"] space:4 n1["1"] space:2
     space:8
     t1["t₁"] space:2 n2["2"] space:3 n3["3"]
     space:8
     t0["t₀"] n4["4"] space:3 n5["5"] space:2
     n4 -- "parentIndex" --> n2
     n5 -- "parentIndex" --> n2
     n4 -- "siblingIndex" --> n5
     n2 -- "parentIndex" --> n1
     n3 -- "parentIndex" --> n1
     n2 -- "siblingIndex" --> n3
     style t0 fill:none,stroke:none
     style t1 fill:none,stroke:none
     style t2 fill:none,stroke:none

* ``parentIndex`` is the halo into which this halo grows at the next timestep of the merger tree. Where several halos share the same parent (halos :math:`4` and :math:`5` in :numref:`fig-outputTreeStructureIsolated`) they are merging with each other, and the parent is the halo which results. A halo at the base of the tree has no parent, and so has ``parentIndex`` :math:`=-1`.
* ``siblingIndex`` is the next halo in the list of halos sharing that same parent, or :math:`-1` for the last halo in the list. Following ``siblingIndex`` from any halo therefore enumerates the remaining halos with which it is merging. (Note that the output contains no dataset naming the *first* halo in that list, so the list can be entered only from a halo already known to be in it---unlike the list of subhalos below, which is entered via ``satelliteIndex``.)

Links between a halo and its subhalos
-------------------------------------

For subhalos the *same three datasets* describe the substructure hierarchy within a single output instead.

.. mermaid::
   :config: {"block": {"padding": 40}}
   :name: fig-outputTreeStructureSubhalos
   :align: center
   :caption: Links between a halo and its subhalos, all at a single output time. Halo :math:`10` hosts subhalos :math:`11` and :math:`12` (drawn with rounded corners), and subhalo :math:`11` in turn hosts the sub-subhalo :math:`13`. As in :numref:`fig-outputTreeStructureIsolated`, each arrow points in the direction of the index which labels it and is drawn only once---subhalo :math:`11` also has ``parentIndex`` :math:`=10`, and sub-subhalo :math:`13` has ``parentIndex`` :math:`=11`.

   block-beta
     columns 6
     space:2 h10["10"] space:3
     space:6
     space s11("11") space:2 s12("12") space
     space:6
     space s13("13") space:4
     h10 -- "satelliteIndex" --> s11
     s12 -- "parentIndex" --> h10
     s11 -- "siblingIndex" --> s12
     s11 -- "satelliteIndex" --> s13

* ``parentIndex`` is the halo which *hosts* this subhalo---not a descendant.
* ``siblingIndex`` is the next subhalo in the same host's list of subhalos, or :math:`-1` for the last subhalo in that list.
* ``satelliteIndex`` is the first subhalo of this node, or :math:`-1` if it has none. It has the same meaning for isolated halos and for subhalos, and is the entry point to the list which ``siblingIndex`` then walks.

The nodes of :numref:`fig-outputTreeStructureSubhalos` appear in the output as:

.. list-table::
   :header-rows: 1
   :widths: auto

   * - ``nodeIndex``
     - ``nodeIsIsolated``
     - ``parentIndex``
     - ``siblingIndex``
     - ``satelliteIndex``
   * - 10
     - 1
     - -1
     - -1
     - 11
   * - 11
     - 0
     - 10
     - 12
     - 13
   * - 12
     - 0
     - 10
     - -1
     - -1
   * - 13
     - 0
     - 11
     - -1
     - -1

(the ``parentIndex`` of :math:`-1` for halo :math:`10` assumes that it is at the base of its merger tree; otherwise it would be the index of the halo into which halo :math:`10` will grow.)

The two pictures are connected: a subhalo is simply a halo whose branch has already joined another branch. In :numref:`fig-outputTreeStructureIsolated`, when halos :math:`4` and :math:`5` merge at :math:`t_1`, halo :math:`4`---the primary progenitor---is promoted to become halo :math:`2`, while halo :math:`5` becomes a subhalo of it, and from that point on halo :math:`5` is linked into the satellite list of halo :math:`2`, exactly as in :numref:`fig-outputTreeStructureSubhalos`.

How to read a row
-----------------

.. list-table::
   :header-rows: 1
   :widths: auto

   * -
     - ``nodeIsIsolated`` :math:`=1` (isolated halo)
     - ``nodeIsIsolated`` :math:`=0` (subhalo)
   * - ``parentIndex``
     - The halo into which this halo grows at the next timestep of the merger tree (:math:`-1` at the base of the tree).
     - The halo which hosts this subhalo.
   * - ``siblingIndex``
     - The next halo sharing the same parent, i.e. the next halo with which this halo is merging (:math:`-1` if there is none).
     - The next subhalo in the same host (:math:`-1` if there is none).
   * - ``satelliteIndex``
     - The first subhalo of this halo (:math:`-1` if it has none).
     - The first sub-subhalo of this subhalo (:math:`-1` if it has none).

.. _sec-output-tree-structure-reference:

Reference: all index datasets
-----------------------------

Only ``nodeIndex``, ``parentIndex``, ``siblingIndex``, ``satelliteIndex``, and ``nodeIsIsolated`` are output by default. The remaining datasets below are produced by property extractors which must be added to the ``nodePropertyExtractor`` parameter to appear in the output---their names are unambiguous, so they are a good choice if you would rather not have to test ``nodeIsIsolated`` in your analysis.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Dataset
     - Extractor
     - Meaning
   * - ``nodeIndex``
     - :galacticus-class:`nodePropertyExtractorNodeIndices`
     - Identifies the node. Unique within a merger tree, but *not* between trees---see :ref:`sec-output-tree-structure-caveats`.
   * - ``parentIndex``
     - :galacticus-class:`nodePropertyExtractorNodeIndices`
     - Descendant or host, depending on ``nodeIsIsolated`` (:math:`-1` if there is neither).
   * - ``siblingIndex``
     - :galacticus-class:`nodePropertyExtractorNodeIndices`
     - Next progenitor of the same parent, or next subhalo of the same host, depending on ``nodeIsIsolated`` (:math:`-1` if there is neither).
   * - ``satelliteIndex``
     - :galacticus-class:`nodePropertyExtractorNodeIndices`
     - First subhalo of this node (:math:`-1` if it has none).
   * - ``nodeIsIsolated``
     - :galacticus-class:`nodePropertyExtractorNodeIndices`
     - :math:`1` for an isolated halo, :math:`0` for a subhalo.
   * - ``hostIndex``
     - :galacticus-class:`nodePropertyExtractorIndicesHost`
     - The host of this subhalo. For an isolated halo the node's *own* ``nodeIndex`` is returned---this extractor never returns :math:`-1`. With ``topLevel`` set to ``true`` the index of the host at the top of the hierarchy is returned instead of that of the direct host, so that a sub-subhalo reports the isolated halo which ultimately contains it.
   * - ``descendantIndex``
     - :galacticus-class:`nodePropertyExtractorDescendants`
     - The node containing this galaxy at the *next output time* (its :term:`forward descendant`). At the final output, or if the galaxy survives to the next output, the node's own index is returned---this extractor never returns :math:`-1`. This is the dataset to use to trace a galaxy from one output to the next.
   * - ``finalDescendantIndex``
     - :galacticus-class:`nodePropertyExtractorFinalDescendant`
     - The node at the base of the merger tree into which this node ultimately merges.
   * - ``mergerTreeIndex``
     - :galacticus-class:`nodePropertyExtractorIndicesTree`
     - The index of the merger tree containing this node. Note that the enclosing ``Outputs/Output``\ *n* group also contains a dataset of this name, which lists the index of each *tree* (one entry per tree, rather than one per node).
   * - ``nodeIndexLastHost``
     - :galacticus-class:`nodePropertyExtractorIndexLastHost`
     - The host halo at the time this node most recently became a subhalo. Requires the :galacticus-class:`nodeOperatorIndexLastHost` node operator to record it.
   * - ``nodeIndexBranchTip``
     - :galacticus-class:`nodePropertyExtractorIndexBranchTip`
     - The earliest progenitor along this node's branch. Requires the :galacticus-class:`nodeOperatorIndexBranchTip` node operator to record it.

Recipes
-------

The examples below use `h5py <https://www.h5py.org/>`_, and all begin by selecting the nodes belonging to a single merger tree. Because ``nodeIndex`` is unique only within a tree, index look-ups must always be done tree-by-tree: the ``mergerTreeIndex``, ``mergerTreeStartIndex``, and ``mergerTreeCount`` datasets in each output group give, for each tree, its index and the range of rows in the ``nodeData`` datasets which belong to it (see :ref:`manual-sec-mergerTreeDatasets` for the full structure of the output file).

.. code-block:: python

   import h5py
   import numpy as np

   model  = h5py.File('galacticus.hdf5','r')
   output = model['Outputs/Output1']
   nodes  = output['nodeData']

   # Select the rows belonging to the first merger tree in this output.
   treeFirst = 0
   start     = output['mergerTreeStartIndex'][treeFirst]
   count     = output['mergerTreeCount'     ][treeFirst]
   rows      = slice(start,start+count)

   nodeIndex      = nodes['nodeIndex'     ][rows]
   parentIndex    = nodes['parentIndex'   ][rows]
   siblingIndex   = nodes['siblingIndex'  ][rows]
   satelliteIndex = nodes['satelliteIndex'][rows]
   isolated       = nodes['nodeIsIsolated'][rows]

   # Build a look-up from node index to row number within this tree.
   row = { index: i for i, index in enumerate(nodeIndex) }

Find the central galaxy and all of its subhalos
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Isolated halos are those with ``nodeIsIsolated`` :math:`=1`. To enumerate the subhalos of one of them, start from its ``satelliteIndex`` and follow ``siblingIndex`` until reaching :math:`-1`:

.. code-block:: python

   def subhalosOf(index):
       """Return the indices of all direct subhalos of the node with the given index."""
       subhalos = []
       subhalo  = satelliteIndex[row[index]]
       while subhalo != -1:
           subhalos.append(subhalo)
           subhalo = siblingIndex[row[subhalo]]
       return subhalos

   central = nodeIndex[np.nonzero(isolated == 1)[0][0]]
   print(central, subhalosOf(central))

This finds only the *direct* subhalos. To include sub-subhalos, apply the same function recursively to each subhalo found---or add the :galacticus-class:`nodePropertyExtractorIndicesHost` extractor with ``topLevel`` set to ``true``, and simply select all nodes whose ``hostIndex`` equals the index of the central (which selects the central itself too, since an isolated halo is its own host).

Note that a merger tree may contain more than one isolated halo at a given output time---any branch which has not yet merged with another is isolated. If you want the most massive one, select on ``basicMass`` rather than assuming the first row.

Enumerate the progenitors of a halo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same ``siblingIndex`` walk applies to isolated halos, where it enumerates the halos which are merging together. Given any halo, the halos which follow it in its parent's list of progenitors are:

.. code-block:: python

   def siblingsOf(index):
       """Return the indices of the halos which follow the given halo in its parent's list of progenitors."""
       siblings = []
       sibling  = siblingIndex[row[index]]
       while sibling != -1:
           siblings.append(sibling)
           sibling = siblingIndex[row[sibling]]
       return siblings

Trace a galaxy from one output to the next
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``parentIndex`` is *not* the right tool for this: for an isolated halo it names a node at the next timestep of the merger tree, which will generally not correspond to any output time. Instead, add the :galacticus-class:`nodePropertyExtractorDescendants` extractor to the output:

.. code-block:: xml

   <nodePropertyExtractor value="multi">
     <nodePropertyExtractor value="nodeIndices"/>
     <nodePropertyExtractor value="descendants"/>
   </nodePropertyExtractor>

and match the ``descendantIndex`` of each node in one output to the ``nodeIndex`` of a node in the next output (within the same merger tree). Note that for merger trees read from file this requires that the ``presetMergerNodes`` and ``presetMergerTimes`` parameters both be set to ``true``, so that Galacticus knows which node each galaxy will merge with, and when. Alternatively, include the ``indexShift`` node operator:

.. code-block:: xml

   <nodeOperator value="indexShift"/>

which changes the convention used for ``nodeIndex`` itself so that a galaxy keeps the same index throughout its evolution---see the description of :galacticus-class:`nodePropertyExtractorNodeIndices` for the two conventions.

.. _sec-output-tree-structure-caveats:

Caveats
-------

* **Node indices are unique only within a merger tree.** Two nodes in different trees may share the same ``nodeIndex``. Always restrict index look-ups to a single tree, using ``mergerTreeIndex``, ``mergerTreeStartIndex``, and ``mergerTreeCount``, or pair ``nodeIndex`` with ``mergerTreeIndex`` to form a globally unique identifier.

* **An isolated halo's** ``parentIndex`` **points into the merger tree, not into the next output.** It is the halo into which this halo grows at the next timestep *of the tree*, and that timestep will generally fall between two output times, so no row with that ``nodeIndex`` need exist in any output. Use ``descendantIndex`` to move between outputs. By contrast, an isolated halo's ``siblingIndex``, and all of the links between a subhalo and its host, generally do refer to nodes which are present in the same output: those nodes are being evolved concurrently.

* **A node's** ``nodeIndex`` **changes when it is promoted.** By default ``nodeIndex`` is the index of the node in the original merger tree, so when a galaxy's halo is promoted into its parent at a tree timestep the index changes, even though the galaxy is the same object. The ``indexShift`` node operator changes this behavior; both conventions are illustrated in the description of :galacticus-class:`nodePropertyExtractorNodeIndices`.

* **The** :math:`-1` **sentinel is not used consistently between extractors.** ``parentIndex``, ``siblingIndex``, and ``satelliteIndex`` are :math:`-1` when the corresponding link does not exist, but ``hostIndex`` reports the node's own index for an isolated halo, and ``descendantIndex`` reports the node's own index when the galaxy survives to the next output. This is deliberate: ``hostIndex`` follows the convention used for the same-named dataset in merger tree files (see :doc:`data/merger-tree-file-format`), in which a halo which is not a subhalo is its own host.

* **The name** ``descendantIndex`` **means something different in merger tree input files.** In a merger tree file read by Galacticus it is the descendant at the next timestep *of that file*; in the output it is the descendant at the next *output time*. See :doc:`data/merger-tree-file-format`.

* If you are working with the ``treeNode`` objects themselves---in Fortran, or through the Python interface---rather than with the output datasets, the same structure is exposed through the ``parent``, ``firstChild``, ``sibling``, and ``firstSatellite`` pointers. Their traversal is described in :doc:`../developer-guide/traversing-a-merger-tree`.
