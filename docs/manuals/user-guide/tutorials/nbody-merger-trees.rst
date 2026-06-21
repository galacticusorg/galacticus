Using N-body Merger Trees
=========================

See `the documentation <https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/index.html>`_ for details of how to build merger tree
files suitable for input into Galacticus. There are many options which control
precisely how merger trees read from file should be handled. The
following section provides guidance on the best choice of parameters.

Setting Input Parameters
------------------------

To utilize merger trees from the file [#footnote1]_ that you created in a Galacticus run
it's necessary to set two blocks of parameters in the input parameter file that
you will use for the run:

.. code-block:: xml

   <!-- Specify that merger trees are to be read from file, give the name of the file to read, and the type of importer to use. -->
   <mergerTreeConstructor value="read">
     <fileNames                value="treeFile.hdf5"/>
     <treeIndexToRootNodeIndex value="false"        />
     <outputTimeSnapTolerance  value="0.001"        />
     <presetMergerTimes        value="true"         />
     <presetOrbits             value="true"         />
     <presetSpins              value="false"        />
     <presetPositions          value="true"         />
     <presetSubhaloMasses      value="false"        />
     <presetScaleRadii         value="false"        />
     <presetMergerNodes        value="true"         />
   </mergerTreeConstructor>
   <mergerTreeImporter value="galacticus">
     <fatalMismatches          value="true" />
     <reweightTrees            value="false"/>
     <validate                 value="false"/>
   </mergerTreeImporter>

The first of these ``mergerTreeConstruct=read`` tells Galacticus that merger trees will be
constructed by reading them from a file. The second, ``mergerTreeImporter``, gives the type of importer to use - in this case we use the ``galacticus`` importer which knows how to import data from the standard Galacticus merger tree file format.

There are many sub-parameters which can be set to change the behavior of merger tree reading. In order to choose sensible settings for the various parameters that
control merger trees read from file, it is recommended that you read
through each of the items below and follow the guidance given.

**Cosmology:** In addition to specifying that trees should
be read from a file, it's also important to ensure that the values of
cosmological parameters in Galacticus match those in the merger tree file. (If
they don't match, Galacticus will stop with an error message unless you set
``fatalMismatches=false`` in which case you'll just be warned about any
mismatch.) In our case of using merger trees from the Millennium
Simulation, the correct cosmological parameter values can be set as
follows:

.. code-block:: xml

   <!-- Use Millennium Simulation cosmology. -->
   <cosmologyParameters value="simple"/>
    <HubbleConstant  value="73.0"  />
    <OmegaMatter     value="0.25"  />
    <OmegaDarkEnergy value="0.75"  />
    <OmegaBaryon     value="0.0455"/>
    </cosmologyParameters>
   <cosmologicalMassVariance value="filteredPower">
    <sigma_8 value="0.900"/>
    </cosmologicalMassVariance>
   <powerSpectrumPrimordial value="powerLaw">
     <index               value="1.000"/>
     <wavenumberReference value="1.000"/>
     <running             value="0.000"/>
   </powerSpectrumPrimordial>

**Existence at Final Time:** Normally, Galacticus assumes that all
merger trees will exist (i.e. have at least one node present) at the
final output time. This may not be true of trees extracted from an
N-body simulation—in this case Galacticus can be informed of this fact by setting the subparameter ``allTreesExistAtFinalTime=false`` in the ``mergerTreeEvolver`` block:

.. code-block:: xml

   <mergerTreeEvolver value="standard">
     <allTreesExistAtFinalTime value="false"/>
   </mergerTreeEvolver>

**Snapping Nodes to Snapshots:** N-body merger trees are
often built from "snapshots" of the simulation, i.e. all of the nodes
exist at a set of discrete times. Often we want to output nodes at
precisely these output times. In such cases it is useful to set the subparameter in the ``mergerTreeConstructor`` block:

.. code-block:: xml

   <outputTimeSnapTolerance  value="0.001"/>

which ensures that the times of nodes are adjusted to lie at precisely
the output time if that time is within the specified relative tolerance
(this avoids any small differences between node times and output times
that can arise due to rounding errors when converting from redshifts to
times and vice-versa).

**Missing Hosts:** Galacticus expects to find each node's host present in
a merger tree file. If a node's host is not found this is cause for a fatal
error to be issued, since it is impossible to correctly construct and
evolve the corresponding tree. If you absolutely want to run a tree for which one
or more hosts are missing, you can allow this by setting the subparameter in the ``mergerTreeConstructor`` block
``missingHostsAreFatal=false``—in this case missing hosts trigger a warning only
and nodes without a host are forced to become isolated nodes. This will lead to
incorrect tree evolution however, so the recommended setting is:

.. code-block:: xml

   <missingHostsAreFatal value="true"/>

**Branch Jumps and Subhalo Promotions:** If your merger
trees contain subhalos they will most likely exhibit two specific
behaviors [#footnote2]_: i) halos which are subhalos in one timestep may become
non-subhalos (isolated halos) in a subsequent timestep ("subhalo
promotion"), and ii) halos which are subhalos in one branch of the tree may
"jump" to another branch [#footnote3]_ of the tree becoming a subhalo there
("branch jumping"). These behaviors are fully supported by Galacticus and so we
recommend the following settings for subparameters in the ``mergerTreeConstructor`` block:

.. code-block:: xml

   <allowSubhaloPromotions value="true"/>
   <allowBranchJumps       value="true"/>

You may choose to disallow these behaviors by setting either of the
above parameters to ``false``—for example if you wish to
explore how your results would differ if subhalos were forced to remain
subhalos forever in their original branch. Note that allowing subhalo
promotion while not allowing branch jumping can lead to in merger tree
evolution, so change these settings with caution.

Note that for trees which do not contain subhalos these two parameters
are irrelevant.

**Subhalo Masses:** If your trees contain subhalos, the
mass evolution of those subhalos can be preset in the satellite
component of each node. In this way, the subhalo mass in Galacticus will track that
specified by the merger tree file. This requires the use of a satellite
component which allows presetting of subhalo masses. Recommended
settings are therefore:

.. code-block:: xml

   <componentSatellite     value="preset"/>
   <mergerTreeConstructor value="read">
     <presetSubhaloMasses value="true"  />
   </mergerTreeConstructor>

If your trees do not contain subhalos, recommended settings are instead:

.. code-block:: xml

   <componentSatellite     value="standard"/>
   <mergerTreeConstructor value="read">
     <presetSubhaloMasses value="false" />
   </mergerTreeConstructor>

**Halo Positions/Velocities:** If your trees contain
position and velocity information for halos, those positions and
velocities can be preset in the position component of each node. This
requires the use of a position component which allows presetting of
positions and velocities. Recommended settings are therefore:

.. code-block:: xml

   <componentPosition      value="preset"/>
   <mergerTreeConstructor value="read">
     <presetPositions value="true"  />
   </mergerTreeConstructor>

If your trees do not contain position information recommended settings
are:

.. code-block:: xml

   <componentPosition      value="null"/>
   <mergerTreeConstructor value="read">
     <presetPositions value="false" />
   </mergerTreeConstructor>

**Subhalo Orbits:** If your trees contain position and
velocity information they can be used to preset initial orbit
information for subhalos. Note that it is not required that your trees
contain subhalos for this orbit presetting to be performed—Galacticus can follow
subhalo orbits even if subhalos are not included in the trees
themselves. The following settings are recommended:

.. code-block:: xml

   <mergerTreeConstructor value="read">
     <presetOrbits             value="true"/>
     <presetOrbitsSetAll       value="true"/>
     <presetOrbitsAssertAllSet value="true"/>
     <presetOrbitsBoundOnly    value="true"/>
   </mergerTreeConstructor>

These options will cause an orbit to be preset for each subhalo based on
the relative position and velocity of merging halos and assuming that
the orbital energy and angular momentum are conserved between the time
immediately prior to the merger and the time of virial radius crossing.
If the computed orbit does not cross the virial radius of the larger
halo or if the computed orbit is unbound, the above options cause an
orbit to be preset by drawing orbital parameters at random from the
chosen cosmological distribution (see `the documentation <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_).

**Subhalo Merging:** If your merger trees contain subhalo
information, that information can be used to specify when, and with
which other node, each subhalo merges. Specifically, a subhalo is
assumed to merge at the time at which it is not the primary progenitor
of its descendent halo—possibly with some other delay to be described
below. Recommended settings are:

.. code-block:: xml

   <mergerTreeConstructor value="read">
     <presetMergerTimes                             value="true"             />
     <presetMergerNodes                             value="true"             />
     <satelliteMergingTimescalesSubresolution value="boylanKolchin2008"/>
   </mergerTreeConstructor>

The first two options cause subhalos to merge at the time described
above, and with their descendent node. The final option accounts for the
possibility that the subhalo should not actually merge immediately at
this time. For example, in N-body simulations, the subhalo may have
simply been lost due to limitations of resolution or halo finder
algorithms. The final option specifies that some additional time until
merging be added based on the subhalo merging timescale algorithm of `Boylan-Kolchin et al. (2008 <https://ui.adsabs.harvard.edu/abs/2008MNRAS.383...93B>`_; `boylanKolchin2008 <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_),
and computed using the last known orbital properties of the subhalo.

**Halo Scale Radii:** If your merger trees contain
information on halo scale radii or half-mass radii, these can be used to
preset the scale radius of each node. This requires the use of a dark matter
profile component which allows presetting of scale length. Recommended
settings are therefore:

.. code-block:: xml

   <mergerTreeConstructor value="read">
     <presetScaleRadii                     value="true"     />
     <presetScaleRadiiFailureIsFatal       value="true"     />
     <presetScaleRadiiConcentrationMinimum value="3"        />
     <presetScaleRadiiConcentrationMaximum value="60"       />
     <presetScaleRadiiMinimumMass          value="see below"/>
   </mergerTreeConstructor>

Minimum and maximum concentrations are specified—these are used to
restrict the range of scale radii that are allowed for a given halo. If
scale radii are to be determined based on half-mass radii given in the
merger tree file, and if the computed scale radius does not result in a
concentration between these limits, then a fatal error is issued.

Finally, you can set a minimum halo mass via the parameter below which
the scale radii or half-mass radii in your file should be considered not
reliable. For halos below this mass, scale radii will instead be
assigned via the selected dark matter halo concentration method (see
`the documentation <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_).

**Halo Angular Momenta:** If your merger trees contain spin
or angular momentum information these can be preset for each node.
Recommended settings are:

.. code-block:: xml

   <componentSpin          value="preset"/>
   <mergerTreeConstructor value="read">
     <presetSpins           value="true"  />
     <presetUnphysicalSpins value="true"  />
   </mergerTreeConstructor>

The last of these options causes any halos for which the spin given in
the merger tree file is non-positive to be assigned a spin at random
instead, drawn from the specified cosmological distribution (see
`the documentation <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_).

If subhalo masses are not included in their host halo masses in your
merger tree file, you should specify how the angular momenta of subhalos
should be accounted for when adding their mass to their host halo. If
positions and 3D angular momenta are available in your merger tree file,
the recommended setting is:

.. code-block:: xml

   <mergerTreeConstructor value="read">
     <subhaloAngularMomenta value="summation"/>
   </mergerTreeConstructor>

If this information is not present

.. code-block:: xml

   <mergerTreeConstructor value="read">
     <subhaloAngularMomenta value="scale"    />
   </mergerTreeConstructor>

should be used instead.

If your merger tree file contains 3D spin or angular momentum
information, you can choose to make that information available within
Galacticus by using the settings:

.. code-block:: xml

   <componentSpin          value="preset3D"/>
   <mergerTreeConstructor value="read">
     <presetSpins3D value="true"    />
   </mergerTreeConstructor>

**Subhalo Indices:** If your merger trees contain subhalos,
you can tell Galacticus to keep track of the indices of subhalos by setting:

.. code-block:: xml

   <componentSatellite     value="preset"/>
   <mergerTreeConstructor value="read">
     <presetSubhaloIndices value="true"  />
   </mergerTreeConstructor>

The Galacticus output file will then contain ``satelliteNodeIndex``
datasets which list the index (as given in the merger tree file) for all
subhalos and halos. Without specifying this presetting, the index of
subhalos is frozen at the index of the halo immediately prior to it
becoming a subhalo.

The remainder of this section gives more detail about many of the
parameters described above and how they affect handling of merger trees
read from file. Further parameters can be set to control what
information from the stored trees will be used in Galacticus' calculations. Examples are given
below.

Further Details
---------------

Further details of the effects of the many parameters controlling merger
trees read from file are given below.

Node Positions
~~~~~~~~~~~~~~

If position and velocity information for tree nodes is available within
the merger tree file then Galacticus can be instructed to use this information by
using the "preset" method for tree node positions and telling the merger
tree construction method to preset node positions as follows:

.. code-block:: xml

   <!-- Use merger tree node positions -->
   <componentPosition      value="preset"/>
   <mergerTreeConstructor value="read">
     <presetPositions value="true"  />
   </mergerTreeConstructor>

If position information is unavailable, the "null" position method can
be selected and the merger tree construction method instructed not to
preset positions as follows:

.. code-block:: xml

   <!-- Do not use merger tree node positions -->
   <componentPosition      value="null"/>
   <mergerTreeConstructor value="read">
     <presetPositions value="false"  />
   </mergerTreeConstructor>

Virial Orbits
~~~~~~~~~~~~~

If position and velocity information for tree nodes is available within
the merger tree file then Galacticus can be instructed to use this information to
estimate the orbit of each subhalo at the point at which it crosses the
virial radius of its host halo. This "virial orbit" may then be used by,
for example, calculations of merging timescales.

.. code-block:: xml

   <!-- Use merger tree node positions to compute orbits at the virial radius -->
   <mergerTreeConstructor value="read">
     <presetOrbits             value="true"/>
     <presetOrbitsBoundOnly    value="true"/>
     <presetOrbitsSetAll       value="true"/>
     <presetOrbitsAssertAllSet value="true"/>
   </mergerTreeConstructor>

Typically, a merging halo is not seen at precisely the time at which it
crosses the virial radius of its host (due to the fact that N-body
simulations are output at discretely spaced timesteps). Therefore,
Galacticus computes the orbit at the time just prior to merging and assumes that
the orbital parameters (energy and angular momentum) remain fixed to
propagate the orbit to the virial radius of the host. The second
parameter in the above example, ``presetOrbitsBoundOnly``, specifies whether or not only bound
orbits should be set. Some calculations (e.g. of subhalo merging times)
assume bound orbits and may fail if given an unbound orbit. Setting this
option to ``true`` causes only bound orbits to be
preset—unbound orbits are ignored. Note that some orbits cannot be
propagated to the virial radius (i.e. their pericenter is larger than
the virial radius). The option, if true, will cause such orbits to be
assigned randomly using the selected virial orbits distribution function, such that all orbits are
assigned. The option requires that all orbits be set—if
``presetOrbitsSetAll=false`` and ``presetOrbitsAssertAllSet=true`` then Galacticus will exit with an
error message if any orbit cannot be set.

If the satellite component additionally permits setting of the satellite
position and velocity, these properties will also be assigned based on
the relative position and velocity of the satellite and host halos.

Merging Times and Targets
~~~~~~~~~~~~~~~~~~~~~~~~~~

The times at which subhalos merge with their host halo can be determined
directly from the merger tree file if subhalo information is included in
that file. Merging is assumed to occur when the subhalo no longer has a
distinct descendent (i.e. it descends into a non-subhalo). If merging
times are to be computed in this way set

.. code-block:: xml

   <componentSatellite     value="preset"/>
   <mergerTreeConstructor value="read"   >
    <presetMergerTimes value="true"  />
   </mergerTreeConstructor>

which select a satellite orbit method that allows merger times to be
preset and tell the merger tree construction method to preset those
merger times respectively. If merger times are not to be computed in
this way then instead set, for example,

.. code-block:: xml

   <componentSatellite         value="standard" />
   <satelliteMerging          value="jiang2008"/>
   <mergerTreeConstructor value="read">
    <presetMergerNodes value="false"/>
   </mergerTreeConstructor>

which selects a standard satellite orbit method, prevents attempts to
preset the merger times and selects the ``jiang2008`` method
for computing merger times instead.

In addition to setting the times of merger events, it is possible to set
the target node with which a merging node should merge. By default,
Galacticus will assume that all merging occurs with the non-subhalo host node in
which a subhalo is located. This may not be the desired behavior when
using N-body merger trees. For example, such trees may indicate that a
subhalo merges with another subhalo. Setting

.. code-block:: xml

   <mergerTreeConstructor value="read">
    <presetMergerNodes value="true"/>
   </mergerTreeConstructor>

will cause the target node with which each merger should occur to be
determined from the merger tree structure and preset for use in Galacticus.

It is possible to add a delay between the last time at which a subhalo
was seen in a simulation and the time at which it is considered to
merge. This functionality is motivated by the consideration that a
subhalo vanishing from a simulation may be simply due to it dropping
below resolution rather than it actually having undergone a merger. The
parameter can be used to select a satellite merging timescale method
(see `the documentation <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_) to use in this case. (It is
set by default to "``null``" such that no delay before merging
occurs.) The orbit of the subhalo around its parent at the last time it
is present in the merger tree is passed to this method and used to
estimate a time until merging. This delay is added to the time at which
the subhalo merges and, if merge target nodes are being set, the target
node is updated accordingly.

Subhalo Indices
~~~~~~~~~~~~~~~

The indices of subhalos are usually frozen at the index of the halo just
prior to becoming a subhalo. The index of the corresponding halo in the
original tree (as read from file) can be tracked as follows:

.. code-block:: xml

   <componentSatellite     value="preset"/>
   <mergerTreeConstructor value="read">
       <presetSubhaloIndices value="true"  />
   </mergerTreeConstructor>

to first select the "preset" satellite orbit method (which allows
subhalo indices to be preset) and, second, to instruct the merger tree
construction algorithm to preset those indices. The index will then be
available in output as ``satelliteNodeIndex``.

Subhalo Masses
~~~~~~~~~~~~~~

The masses of subhalos (specifically their time evolution after they
become subhalos) can be set using the values stored in the merger tree
file (if available). To set subhalo masses in this way use

.. code-block:: xml

   <componentSatellite      value="preset"/>
   <mergerTreeConstructor value="read">
    <presetSubhaloMasses value="true"  />
   </mergerTreeConstructor>

to first select the "preset" satellite orbit method (which allows
subhalo masses to be preset) and, second, to instruct the merger tree
construction algorithm to preset those masses.

Node Spins
~~~~~~~~~~

If information on the angular momenta of nodes is available in the
merger tree file, this can be used to preset the value of the spin
parameter in each node [#footnote4]_ by setting:

.. code-block:: xml

   <mergerTreeConstructor value="read">
    <presetSpins value="true"/>
   </mergerTreeConstructor>

The spin parameter is set using the spin of each node if available, or
otherwise using the angular momentum of each node stored in the merger
tree file using:
:math:`\lambda = |\mathbf{J}| |E|^{1/2} / G M^{5/2}`,
where :math:`|\mathbf{J}|` is the magnitude of the node's angular momentum, :math:`M` is
the node's mass and :math:`E` is its energy. Additionally, by setting:

.. code-block:: xml

   <mergerTreeReadPresetSpins3D value="true"/>

the spin vector of each node will be set (assuming that the vector spin
or angular momenta of nodes are available in the merger tree file)
using:
:math:`\mathbf{\lambda} = \mathbf{J} |E|^{1/2} / G M^{5/2}`.

If spins could not be determined for some halos the spin (or angular
momentum) should be set to zero in the merger tree file, and the
parameter set to ``true``. Galacticus will then assign a spin to such
halos by sampling from the selected spin distribution (see
`the documentation <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_).

Node Scale Radii
~~~~~~~~~~~~~~~~

If information on the half-mass or scale radii of nodes is available in
the merger tree file, it can be used to preset the value of the dark
matter halo scale radius in each node by setting:

.. code-block:: xml

   <mergerTreeConstructor value="read">
    <presetScaleRadii value="true"/>
   </mergerTreeConstructor>

Before doing this, it is important to be sure that the half-mass or
scale radii of the nodes are reliable. For example, in low mass nodes
extracted from an N-body simulation resolution effect may limit the
accuracy of the measured half-mass or scale radius. In such cases, use
the parameter to specify the lowest mass halos for which the scale radii
should be preset—lower mass halos will be assigned a scale radius using
the method specified by the parameter (which will default to the value
of ; see `the documentation <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_). It is also possible
to specify minimum and maximum allowed concentrations when computing the
scale radius from the half mass radius using the and parameters. If
matching the half mass radius would require a concentration outside of
this range, Galacticus will abort unless ``presetScaleRadiiFailureIsFatal=false``, in which case it
will instead silently use the fallback concentration method described
above.

If only half-mass radii are available, the scale radius is set by using
a root finding algorithm to ensure that half of the total halo mass is
enclosed within the specified half-mass radius.

Miscellaneous N-body Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several miscellaneous properties often available from N-body merger
trees can also be preset by setting the following subparameters of the ``mergerTreeConstructor`` block to
``true``:

``presetParticleCounts``:   Sets the number of particles in each halo (requires the
``particleCount`` dataset to be present in the merger tree
file);

``presetVelocityMaxima``:   Sets the maxima of halo rotation curves (requires the
``velocityMaximum`` dataset to be present in the merger
tree file);

``presetVelocityDispersions``:   Sets the velocity dispersion of halos (requires the
``velocityDispersion`` dataset to be present in the merger
tree file).

Subhalo Promotion
~~~~~~~~~~~~~~~~~

A subhalo may, at a later time, become an isolated halo once again.
Galacticus allows you to control whether such behavior is allowed, or should be
prohibited. To allow such "subhalo promotion", set:

.. code-block:: xml

   <mergerTreeConstructor value="read">
    <allowSubhaloPromotions value="true"/>
   </mergerTreeConstructor>

If you choose to inhibit this behavior by setting the above parameter to
false, a halo that becomes a subhalo will remain a subhalo forever
thereafter. Note that the isolated halo to which it would have been
promoted will still exist, and may therefore form its own galaxy. This
can result in double counting of mass, and so inhibiting subhalo
promotion is not recommended.

"Fly-by" Halos
~~~~~~~~~~~~~~

In some cases, a halo that is part of one tree can later become part of
another tree. This can happen in so-called "fly-by" encounters where a
halo may briefly become a subhalo in a halo in tree A then leave that
halo and become a subhalo in tree B.

The correct way to handle this issue is to combine trees A and B into a
single tree (which will now have multiple base nodes). Galacticus will then
process these two trees simultaneously, correctly handling the fly-by,
and outputting the trees as two separate trees.

If for some reason this is not possible or desired, the fly-by problem
will normally cause Galacticus to complain that the host halo of a node cannot be
found (since it exists in a different tree). This problem can be avoided
by setting:

.. code-block:: xml

   <mergerTreeConstructor value="read">
    <missingHostsAreFatal value="false"/>
   </mergerTreeConstructor>

In this case, nodes with missing hosts are simply treated as being
isolated halos. This will avoid an error condition, but is not a
physically correct way to handle such cases, so use with caution.

Using Particles to Track Unresolved Subhalos
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In N-body simulations it is possible that a subhalo can become "lost"
from the simulation (i.e. can no longer be identified by a halo finder
due to resolution issues) before it has actually merged with the central
galaxy or been completely tidally destroyed. In such cases it is useful
to be able to assign a position to the subhalo at later times. A common
approach to assigning a position (and velocity) is to use that of the
most bound particle in the subhalo at the last time it was identified.
Galacticus allows for particle tracking in this way through the addition of
particle information to the merger tree file.

To add particle tracking data to a merger tree file, follow these steps:

#. Identify all subhalos which are lost from the simulation prior to
   the final timestep;

#. Determine the index of the most bound particle in each such subhalo
   in the last timestep in which it was identified;

#. For each subhalo, extract the redshift, position, and velocity of
   that particle (which is usually trivial to do once its index
   is known) at each subsequent timestep in the simulation;

#. Write these data (along with the particle index) to the
   ``particles`` group in the merger tree file as described
   `here <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/data/merger-tree-file-format.html#particles-group>`_;

#. Add two datasets to the ``forestHalos`` group:

   #. ``particleIndexStart`` which should indicate the index
      in the datasets in the ``particles`` group at which the
      data for each halo begins (or -1 if no particle data is
      included for the halo);

   #. ``particleIndexCount`` which should indicate the number
      of entries in the datasets in the ``particles`` group
      for each halo (or -0 if no particle data is included for
      the halo).

Handling of Extremely Large Merger Tree Forests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Halos can move between merger trees (see
"`Subhalo Promotion <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/tutorials/nbody-merger-trees.html#subhalo-promotion>`_" and
"`Fly-by Halos <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/tutorials/nbody-merger-trees.html#fly-by-halos>`_"), leading to the necessity of
merger tree forests—interconnected groups of merger trees that
Galacticus typically processes as a whole. These forests can become very large—in
some cases so large that they do not fit within the available memory.
Galacticus can handle such forests by splitting them into individual trees. Each
tree is processed separately, and nodes are moved between trees as
needed. If a tree needs a node from another tree before its evolution
can continue, its state can be suspended to disk, and later re-read once
the node it requires becomes available. In this way, very large forests
can be processed without running out of memory (as trees are stored to
disk while they are not being processed).

To cause forests to be split in this way, the following parameters
should be set:

.. code-block:: xml

   <task value="evolveForests">
    <suspendToRAM               value="false"            />
    <suspendPath                value="/my/scratch/path/"/>
   </task>
   <mergerTreeConstructor value="read">
    <forestSizeMaximum          value="10000000"         />
    <subresolutionMerging value="infinite"         />
   </mergerTreeConstructor>

Here, ``suspendToRAM`` specifies that merger trees
should be suspended to disk (i.e. not to RAM which is the default), and
``suspendPath`` gives a path where the suspended
trees can be written—typically scratch space local to the compute node
where Galacticus is being run is a good option.

``forestSizeMaximum`` specifies the maximum
number of nodes allowed in a forest before it will be split. A suitable
number for this depends on the details of the available RAM, the number
of threads sharing that RAM, and the characteristics of the Galacticus model being
used (which will affect the memory required per node).

Finally, ``subresolutionMerging`` is set to
``infinite`` to prevent any merging (which is not supported for
split forests at present, although it should be soon).

Analyzing the Output
--------------------

Positions and Velocities
~~~~~~~~~~~~~~~~~~~~~~~~~~

Components of the position of each node are output as
``positionX``, ``positionY`` and
``positionZ`` and can be accessed in the same way as other
output properties from Galacticus (see `the documentation <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/index.html>`_).

Subhalo Masses
~~~~~~~~~~~~~~

The current mass of subhalos is available via the
``nodeBoundMass`` output dataset and can be accessed in the
same way as other output properties from Galacticus (see `the documentation <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/index.html>`_).
For non-subhalos this property is equal to the usual
``nodeMass`` property.

.. rubric:: Footnotes

.. [#footnote1] The following assumes that merger trees will be read from a file
   following Galacticus's standard HDF5 format which is described
   `here <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/data/merger-tree-file-format.html>`_.

.. [#footnote2] These two behaviors are called out as they specifically *do not*
   occur in merger trees created using Press-Schechter-based algorithms
   for example.

.. [#footnote3] That is, the subhalo's descendent is hosted by a halo other than the
   descendent of the subhalo's host.

.. [#footnote4] Before doing this, it is important to be sure that the angular
   momenta of the nodes are reliable. For example, in low mass nodes
   extracted from an N-body simulation resolution effect may limit the
   accuracy of the measured angular momentum.
