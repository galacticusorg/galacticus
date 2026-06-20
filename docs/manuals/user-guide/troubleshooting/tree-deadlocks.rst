Debugging Tree Deadlocks
========================

A tree deadlock occurs when Galacticus is unable to find any node in a merger tree that can be advanced forward in time. This can occur if, for example, a node has an event attached to it that never occurs - thereby preventing it from evolving beyond the time of that event.

Deadlocks can be difficult to diagnose, particularly in very large trees. To aid in this, if Galacticus detects a deadlock it will write an extensive report to the output log, and will also write tree data to a `GraphViz <https://www.graphviz.org/>`_ format file (including information on evolution interdependencies between nodes) for visualization with `dot <https://www.graphviz.org/docs/layouts/dot/>`_ for example, named ``galacticusDeadlockTree_1.gv`` (if more than one deadlocked tree is found then multiple trees are written, with the numerical suffix increasing for 1 for each).

This `GraphViz <https://www.graphviz.org/>`_ file can also be used as a means to detect cycles in the node interdependencies, which can help to isolate the cause of the deadlock. Cycle detection can be performed using `scripts/aux/mergerTreeDeadlockCycleDetector.py <https://github.com/galacticusorg/galacticus/blob/master/scripts/aux/mergerTreeDeadlockCycleDetector.py>`_ script

If cycles are detected, they will be output as follows:

.. code-block:: text

   Cycle in tree 32464292 consisting of 3 nodes:
      (63575) -> (64194) [nodeID: 1272608; time: 0.8684125404; reason: hosted satellite]
      (64195) -> (63575) [nodeID: 880134; time: 0.9984304302; reason: satellite in host]
      (64194) -> (64195) [nodeID: 880133; time: 0.8684125404; reason: mergee ( 0.8684)]

which shows a cycle consisting of three nodes. The numbers displayed in parentheses show indicate a dependency between two nodes, identified by their ``uniqueIDs``, *not* the node indices - the node index is shown as ``nodeID``. Also shown is the time at which the node exists, and the reason why it can not evolve. As can be seen in the above example, these dependencies result in a circular dependence which prevents the tree from evolving.

The `scripts/aux/mergerTreeDeadlockCycleDetector.py <https://github.com/galacticusorg/galacticus/blob/master/scripts/aux/mergerTreeDeadlockCycleDetector.py>`_ script will attempt to diagnose the cause of any cycles. For example in the above case it will show:

.. code-block:: text

   Diagnosis:
      Nodes with mergees present, but these exist after their merging time.
      Check that you have either
         <mergerTreeTimestep value="satellite"/>, or
         <mergerTreeTimestep value="standard"/>
       set in your parameter file.
