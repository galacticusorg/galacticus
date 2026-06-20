Traversing a Merger Tree
========================

In Galacticus, merger trees are represented as a collection of interconnected "nodes" - each of which represents a dark matter halo and all of the content (stars, gas, etc.) contained within it. Nodes are connected to other nodes, having one or "child" nodes (those progenitors which exist at the prior timestep), a "parent" node (the node into which the current node will merge at the next timestep) and, possibly, a "sibling" node (the next "child" with the same "parent").

Often we need to traverse the tree, moving from one node to the next in a certain order to perform some calculation. This tutorial gives some examples of how to do this traversal.

The ``treeNode`` class
----------------------

A node in the tree is represented by an object of the ``treeNode`` class.

Pointers
~~~~~~~~

Importantly, for tree traversal purposes, each node contains pointers to other nodes to define the tree structure. These pointers are:

* ``firstChild`` - a pointer to the first child node;
* ``parent`` - a pointer to the parent node;
* ``sibling`` - a pointer to the next child of the same parent.

Not all all these pointers will actually point to anything - some nodes do not have a sibling for example. In those cases, the relevant pointer will be ``null``.

An example simple tree structure:

.. mermaid::

   graph TD
     1--firstChild-->2
     2--parent-->1
     3--parent-->1
     2--sibling-->3
     2--firstChild-->4
     4--parent-->2
     5--parent-->2
     4--sibling-->5
     3--firstChild-->6
     6--parent-->3
     7--parent-->3
     6--sibling-->7

Methods
~~~~~~~

The ``treeNode`` class also has some useful methods for tree traversal. For now we'll consider only one:

* ``isPrimaryProgenitor()`` - this returns ``true`` if a node is the primary (most massive) progenitor of its parent - that is, does ``node``-->``parent``-->``firstChild`` point back to ``node``.

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

We begin by pointing our worker node at the original ``node`` that we were given. We the enter a loop - in the loop condition we check if our node worker is ``associated()`` - i.e. that it is not ``null``. This will allow us to exit the loop once no more primary progenitors are available. We next do whatever calculation we want to do on ``nodeWorker``, knowing that it points to a primary progenitor of the original ``node``. Then we simply move to the next primary progenitor by following the ``firstChild`` pointer attached to ``nodeWorker``.

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

We begin by pointing our worker node at the original ``node`` that we were given. We the enter a loop - in the loop condition we check if our node worker is ``associated()`` - i.e. that it is not ``null``. This will allow us to exit the loop once no more descendants are available. We next do whatever calculation we want to do on ``nodeWorker``, knowing that it points to a descendant of the original ``node``. Then we check if the current worker node is the primary progenitor of its parent. If it is, we move to that parent by following the ``firstChild`` pointer attached to ``nodeWorker``, otherwise we make ``nodeWorker`` point to ``null`` so that the loop will exit.
