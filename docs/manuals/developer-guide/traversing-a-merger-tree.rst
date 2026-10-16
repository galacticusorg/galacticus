Traversing a Merger Tree
========================

In Galacticus, merger trees are represented as a collection of interconnected "nodes" - each of which represents a dark matter halo and all of the content (stars, gas, etc.) contained within it. Nodes are connected to other nodes by a set of pointers which represent *two* different structures: the progenitor structure of the tree (which halo grows into which as time advances), and the substructure hierarchy within a halo (which subhalo is hosted by which halo). The same pointers are used for both, so their meaning depends on whether the node is an isolated halo or a satellite. The corresponding conventions in the model *output* are described in :doc:`../user-guide/output-tree-structure`.

Often we need to traverse the tree, moving from one node to the next in a certain order to perform some calculation. This tutorial gives some examples of how to do this traversal.

The ``treeNode`` class
----------------------

A node in the tree is represented by an object of the ``treeNode`` class.

Pointers
~~~~~~~~

Importantly, for tree traversal purposes, each node contains pointers to other nodes to define the tree structure. These pointers are:

* ``firstChild`` - a pointer to the first child node, i.e. to the :term:`primary progenitor` of this node;
* ``parent`` - a pointer to the parent node. For an isolated node this is the node into which it will merge at the next timestep, while for a satellite node it is the node which hosts it;
* ``sibling`` - a pointer to the next node in the list to which this node belongs. For an isolated node that is the list of children of its parent, so ``sibling`` points to the next child of the same parent, while for a satellite node it is the list of satellites of its host, so ``sibling`` points to the next satellite in the same host;
* ``firstSatellite`` - a pointer to the first satellite node hosted by this node. Satellites can themselves host satellites, so this pointer is meaningful for satellite nodes too;
* ``mergeTarget`` - a pointer to the node with which this node will merge (if one has been assigned - see :galacticus-class:`mergerTreeConstructorClass`);
* ``firstMergee`` - a pointer to the first node in the list of nodes which will merge with this node;
* ``siblingMergee`` - a pointer to the next node in that list of mergees.

Not all of these pointers will actually point to anything - some nodes do not have a sibling for example. In those cases, the relevant pointer will be ``null``.

Note that each node belongs to exactly one ``sibling`` list: either the list of children of its parent, or the list of satellites of its host. Which of the two is precisely what ``isSatellite()`` (see below) tests. The list of mergees is a separate, third list, linked by ``siblingMergee`` rather than by ``sibling``.

An example simple tree structure, showing the pointers which link progenitors to their descendants:

.. mermaid::
   :config: {"block": {"padding": 40}}
   :name: fig-traversingProgenitorPointers
   :align: center
   :caption: Pointers linking isolated nodes across timesteps. Time increases upward: nodes :math:`4` and :math:`5` at time :math:`t_0` merge to form node :math:`2` at :math:`t_1`, which merges with node :math:`3` to form node :math:`1` at :math:`t_2`. Each pointer is drawn only once, in the direction in which it points---so, for example, node :math:`4` also has ``parent`` pointing to node :math:`2` (``firstChild`` and ``parent`` are inverses of each other along the primary branch), and node :math:`5` has ``sibling`` pointing to ``null``.

   block-beta
     columns 8
     t2["t₂"] space:4 n1["1"] space:2
     space:8
     t1["t₁"] space:2 n2["2"] space:3 n3["3"]
     space:8
     t0["t₀"] n4["4"] space:3 n5["5"] space:2
     n2 -- "firstChild" --> n4
     n5 -- "parent" --> n2
     n4 -- "sibling" --> n5
     n1 -- "firstChild" --> n2
     n3 -- "parent" --> n1
     n2 -- "sibling" --> n3
     style t0 fill:none,stroke:none
     style t1 fill:none,stroke:none
     style t2 fill:none,stroke:none

Satellite nodes are linked into their host in exactly the same way, but using ``firstSatellite`` in place of ``firstChild``:

.. mermaid::
   :config: {"block": {"padding": 40}}
   :name: fig-traversingSatellitePointers
   :align: center
   :caption: Pointers linking a node to the satellites which it hosts, all at a single time. Node :math:`10` hosts satellites :math:`11` and :math:`12` (drawn with rounded corners), and satellite :math:`11` in turn hosts the satellite :math:`13`. As in :numref:`fig-traversingProgenitorPointers`, each pointer is drawn only once---node :math:`11` also has ``parent`` pointing to node :math:`10`, and node :math:`13` has ``parent`` pointing to node :math:`11`.

   block-beta
     columns 6
     space:2 h10["10"] space:3
     space:6
     space s11("11") space:2 s12("12") space
     space:6
     space s13("13") space:4
     h10 -- "firstSatellite" --> s11
     s12 -- "parent" --> h10
     s11 -- "sibling" --> s12
     s11 -- "firstSatellite" --> s13

Methods
~~~~~~~

The ``treeNode`` class also has some useful methods for tree traversal. For now we'll consider only two:

* ``isPrimaryProgenitor()`` - this returns ``true`` if a node is the primary (most massive) progenitor of its parent - that is, does ``node``-->``parent``-->``firstChild`` point back to ``node``;
* ``isSatellite()`` - this returns ``true`` if a node is a satellite, i.e. if it has a ``parent``, but does not appear in the list of children of that parent (and so must instead appear in its list of satellites). Since ``parent`` and ``sibling`` are shared between the two structures, this is the test to use whenever the meaning of those pointers matters.

Traversing examples
-------------------

Visit all primary progenitors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we are given a ``treeNode`` object ``node``, and we want to visit each of its primary progenitors. To do this we create a worker ``treeNode`` pointer that we can use to traverse the tree. Here's the example:

.. code-block:: fortran

   type(treeNode), pointer :: nodeWorker
   ! Begin by pointer our worker node to the node we have been given.
   nodeWorker => node
   ! Begin a loop to visit each primary progenitor. We check that our worker node is associated, and exit the loop if it is not.
   do while (associated(nodeWorker))
      ! Do whatever calculation we want here on the current primary progenitor node currently pointed to by `nodeWorker`.
      .....
      ! Move to the next primary progenitor of the current worker node - this is found by just following the `firstChild` pointer.
      nodeWorker => nodeWorker%firstChild
      ! Go back to the top of the loop.
   end do

We begin by pointing our worker node at the original ``node`` that we were given. We then enter a loop - in the loop condition we check if our node worker is ``associated()`` - i.e. that it is not ``null``. This will allow us to exit the loop once no more primary progenitors are available. We next do whatever calculation we want to do on ``nodeWorker``, knowing that it points to a primary progenitor of the original ``node``. Then we simply move to the next primary progenitor by following the ``firstChild`` pointer attached to ``nodeWorker``.

Visit all satellites of a node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose instead that we want to visit each satellite hosted by a ``treeNode`` object ``node``. The satellites form a linked list, entered via the ``firstSatellite`` pointer of the host and continued via the ``sibling`` pointer of each satellite:

.. code-block:: fortran

   type(treeNode), pointer :: nodeWorker
   ! Begin by pointing our worker node to the first satellite of the node we have been given.
   nodeWorker => node%firstSatellite
   ! Begin a loop to visit each satellite. We check that our worker node is associated, and exit the loop if it is not.
   do while (associated(nodeWorker))
      ! Do whatever calculation we want here on the current satellite node currently pointed to by `nodeWorker`.
      .....
      ! Move to the next satellite in the list of satellites of `node`.
      nodeWorker => nodeWorker%sibling
      ! Go back to the top of the loop.
   end do

Note that this visits only the *direct* satellites of ``node`` - any satellites of those satellites are not visited. To visit those too, apply the same loop recursively to each satellite found. Note also that, because ``sibling`` is used for both children and satellites, this loop works only because we entered the list through ``firstSatellite``: starting instead from a node found via ``firstChild`` would walk the list of children.

Visit all descendants until our branch merges
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we are given a ``treeNode`` object ``node``, and we want to visit each of its descendant halos along its branch, until that branch merges with another branch. To do this we create a worker ``treeNode`` pointer that we can use to traverse the tree. Here's the example:

.. code-block:: fortran

   type(treeNode), pointer :: nodeWorker
   ! Begin by pointer our worker node to the node we have been given.
   nodeWorker => node
   ! Begin a loop to visit each descendant. We check that our worker node is associated, and exit the loop if it is not.
   do while (associated(nodeWorker))
      ! Do whatever calculation we want here on the current descendant node currently pointed to by `nodeWorker`.
      .....
      ! Check if the current worker node is the primary progenitor of its parent.
      if (nodeWorker%isPrimaryProgenitor()) then
         ! It is the primary progenitor, so simply move to the parent node.
         nodeWorker => nodeWorker%parent
      else
         ! It is not the primary progenitor, so this branch is about to merge into a larger branch. We don't want to follow this larger branch so we are finished. Nullify our worker node so that we will exit the loop.
         nodeWorker => null()
      end if
      ! Go back to the top of the loop.
   end do

We begin by pointing our worker node at the original ``node`` that we were given. We then enter a loop - in the loop condition we check if our node worker is ``associated()`` - i.e. that it is not ``null``. This will allow us to exit the loop once no more descendants are available. We next do whatever calculation we want to do on ``nodeWorker``, knowing that it points to a descendant of the original ``node``. Then we check if the current worker node is the primary progenitor of its parent. If it is, we move to that parent by following the ``firstChild`` pointer attached to ``nodeWorker``, otherwise we make ``nodeWorker`` point to ``null`` so that the loop will exit.
