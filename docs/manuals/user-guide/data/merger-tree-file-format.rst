Merger Tree File Format
=======================

Basic File Format
-----------------

Merger trees are stored in HDF5 files for portability and convenience.
Additionally, the format is intended to be sufficiently flexible to
allow it to describe merger trees obtained in a wide variety of ways,
including Monte Carlo algorithms (e.g. extended Press-Schechter
algorithms) and from N-body simulations.

Flexibility and Extensibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All of the groups/datasets in the file except for the
``forestIndex`` and ``forestHalos`` groups are, in
principle, optional. This does not mean that a file created without some
of these optional groups/datasets will be useable by a given code. It is
the responsibility of a given code to check that all data that it
requires is present in the file. You are therefore encouraged to include
as much information as possible when constructing merger tree files.

Additionally, the file format is intended to be extensible. It is
permissible to add additional datasets, for example to describe some
other properties of nodes in each tree. Additional datasets should
follow the structure of currently defined datasets, i.e. they should be
stored as a single dataset combining all trees with nodes listed in the
same order as for other datasets. For additional datasets which might be
of general use you are encouraged to `contact
us <mailto:abenson@carnegiescience.edu>`_ and recommend them for inclusion
in the standard—this allows their name to be standardized.

A Note on Scalar Attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many of the HDF5 attributes discussed in this document are indicated to
be scalar (rank 0) attributes. It is allowable within the standard that
these be pseudo-scalars (rank 1 arrays containing a single element). This
allows such attributes to be created using the
`h5lt <https://support.hdfgroup.org/documentation/hdf5/latest/group___h5_l_t.html>`_ API for
example.

Example File Structure
~~~~~~~~~~~~~~~~~~~~~~~

An example of the structure of such a file, called
“``example.hdf5``” is shown below using the format of
``h5dump``. Each of the groups is described in detail in the
following sections.

.. code-block:: text

   HDF5 "example.hdf5" {
   GROUP "/" {
      ATTRIBUTE "formatVersion"
      GROUP "cosmology" {
      }
      GROUP "groupFinder" {
      }
      GROUP "forestHalos" {
      }
      GROUP "forests" {
      }
      GROUP "particles" {
      }
      GROUP "provenance" {
      }
      GROUP "simulation" {
      }
      GROUP "forestIndex" {
      }
      GROUP "units" {
      }
   }

Format Version Attribute
------------------------

The ``formatVersion`` attribute specifies which version of the
merger tree file format this file follows. The current standard is
version 2.

Cosmology Group
---------------

The ``cosmology`` group describes the cosmological model within
which the merger trees contained in the file were constructed. An
example of this group, showing standard attributes, is given below.

.. code-block:: text

   GROUP "cosmology" {
      ATTRIBUTE "HubbleParam" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0.73
         }
      }
      ATTRIBUTE "OmegaMatter" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0.25
         }
      }
      ATTRIBUTE "OmegaLambda" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0.75
         }
      }
      ATTRIBUTE "OmegaBaryon" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0.045
         }
      }
      ATTRIBUTE "powerSpectrumIndex" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 1
         }
      }
      ATTRIBUTE "sigma_8" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0.9
         }
      }
      ATTRIBUTE "transferFunction" {
         DATATYPE  H5T_STRING {
               STRSIZE 7;
               STRPAD H5T_STR_SPACEPAD;
               CSET H5T_CSET_ASCII;
               CTYPE H5T_C_S1;
            }
         DATASPACE  SCALAR
         DATA {
         (0): "CAMB"
         }
      }
   }

Standard Attributes
~~~~~~~~~~~~~~~~~~~~

The following are standard attributes in the ``cosmology``
group (others may be added as desired).

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Attribute
     - Description
   * - ``HubbleParam``
     - The Hubble parameter in units of 100 km/s/Mpc at :math:`z`\ =0, :math:`h_0`;
   * - ``OmegaMatter``
     - The density of matter (both dark and baryonic matter) in units of
       the critical density at :math:`z`\ =0, :math:`\Omega_\mathrm{M}`;
   * - ``OmegaLambda``
     - The density of dark energy in units of the critical density at
       :math:`z`\ =0, :math:`\Omega_\Lambda`;
   * - ``OmegaBaryon``
     - The density of matter (both dark and baryonic matter) in units of
       the critical density at :math:`z`\ =0, :math:`\Omega_\mathrm{b}`;
   * - ``powerSpectrumIndex``
     - The index of the primordial power spectrum of matter
       fluctuations, i.e. :math:`n_\mathrm{s}` for power spectrum :math:`P(k) \propto k^{n_\mathrm{s}}`;
   * - ``sigma_8``
     - The root-variance of mass fluctuations in real space top-hat spheres
       of radius 8 :math:`h^{-1}` Mpc computed from the :math:`z`\ =0 linear theory power spectrum, :math:`\sigma_8`;
   * - ``transferFunction``
     - A descriptor of the transfer function used to compute the power spectrum.

Group Finder Group
------------------

This group, typically relevant only for merger trees derived from N-body
simulations, describes the characteristics of the group finding
algorithm that was used to find halos in the simulation. An example of
this group, showing standard attributes, is given below.

.. code-block:: text

   GROUP "groupFinder" {
      COMMENT "Group finder parameters."
      ATTRIBUTE "code" {
         DATATYPE  H5T_STRING {
               STRSIZE 7;
               STRPAD H5T_STR_SPACEPAD;
               CSET H5T_CSET_ASCII;
               CTYPE H5T_C_S1;
            }
         DATASPACE  SCALAR
         DATA {
         (0): "SUBFIND"
         }
      }
      ATTRIBUTE "linkingLength" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0.2
         }
      }
      ATTRIBUTE "minimumParticleNumber" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 20
         }
      }
   }

Standard Attributes
~~~~~~~~~~~~~~~~~~~~

The following are standard attributes in the ``groupFinder``
group (others may be added as desired).

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Attribute
     - Description
   * - ``code``
     - The name of the group finding code used in the construction of these
       merger trees;
   * - ``linkingLength``
     - For friends-of-friends group finding algorithms the
       dimensionless (i.e. in units of the mean interparticle spacing)
       linking length used;
   * - ``minimumParticleNumber``
     - The minimum number of particles that a group was required to have in
       order to be included in a merger tree.

Simulation Group
----------------

This group, typically relevant only for merger trees derived from N-body
simulations, describes the characteristics of the simulation from which
the trees were derived. An example of this group, showing standard
attributes, is given below.

.. code-block:: text

    GROUP "simulation" {
      COMMENT "Simulation parameters."
      ATTRIBUTE "ErrTolIntAccuracy" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0.02
         }
      }
      ATTRIBUTE "TypeOfTimestepCriterion" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "boxSize" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 500
         }
      }
      ATTRIBUTE "code" {
         DATATYPE  H5T_STRING {
               STRSIZE 8;
               STRPAD H5T_STR_SPACEPAD;
               CSET H5T_CSET_ASCII;
               CTYPE H5T_C_S1;
            }
         DATASPACE  SCALAR
         DATA {
         (0): "GADGET-2"
         }
      }
      ATTRIBUTE "initialConditions" {
         DATATYPE  H5T_STRING {
               STRSIZE 5;
               STRPAD H5T_STR_SPACEPAD;
               CSET H5T_CSET_ASCII;
               CTYPE H5T_C_S1;
            }
         DATASPACE  SCALAR
         DATA {
         (0): "glass"
         }
      }
      ATTRIBUTE "softeningKernel" {
         DATATYPE  H5T_STRING {
               STRSIZE 6;
               STRPAD H5T_STR_SPACEPAD;
               CSET H5T_CSET_ASCII;
               CTYPE H5T_C_S1;
            }
         DATASPACE  SCALAR
         DATA {
         (0): "spline"
         }
      }
      ATTRIBUTE "softeningPlummerEquivalent" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 0.005
         }
      }
      ATTRIBUTE "startRedshift" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 127
         }
      }
   }

Standard Attributes
~~~~~~~~~~~~~~~~~~~~

The following are standard attributes in the ``simulation``
group (others may be added as desired).

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Attribute
     - Description
   * - ``boxSize``
     - Relevant for cubic volumes typical of cosmological simulations, this
       attributes gives the length of the box in whatever unit system the
       file used (see `Units Group`_);
   * - ``code``
     - The name of the code used to run the simulation;
   * - ``initialConditions``
     - A description of the initial conditions;
   * - ``softeningKernel``
     - A description of the softening kernel used;
   * - ``softeningPlummerEquivalent``
     - The equivalent Plummer softening length;
   * - ``startRedshift``
     - The redshift at which the simulation was begun.

Gadget-specific Standard Attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following are standard attributes in the ``simulation``
group specifically relevant to simulations run with the `Gadget <http://www.mpa-garching.mpg.de/gadget/>`_
code. They typically reflect the values of parameters used by that code.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Attribute
     - Description
   * - ``ErrTolIntAccuracy``
     - The integration accuracy used by Gadget;
   * - ``TypeOfTimestepCriterion``
     - The type of timestepping criterion used by Gadget;
   * - ``SofteningGas``
     - Specifies the (comoving) softening of the first particle group in
       Gadget;
   * - ``SofteningHalo``
     - Specifies the (comoving) softening of the second particle group in
       Gadget;
   * - ``SofteningDisk``
     - Specifies the (comoving) softening of the third particle group in
       Gadget;
   * - ``SofteningBulge``
     - Specifies the (comoving) softening of the fourth particle group in
       Gadget;
   * - ``SofteningStars``
     - Specifies the (comoving) softening of the fifth particle group in
       Gadget;
   * - ``SofteningBndry``
     - Specifies the (comoving) softening of the sixth particle group in
       Gadget;
   * - ``SofteningGasMaxPhys``
     - Specifies the maximum physical softening of the first particle group
       in Gadget;
   * - ``SofteningHaloMaxPhys``
     - Specifies the maximum physical softening of the second particle
       group in Gadget;
   * - ``SofteningDiskMaxPhys``
     - Specifies the maximum physical softening of the third particle group
       in Gadget;
   * - ``SofteningBulgeMaxPhys``
     - Specifies the maximum physical softening of the fourth particle
       group in Gadget;
   * - ``SofteningStarsMaxPhys``
     - Specifies the maximum physical softening of the fifth particle group
       in Gadget;
   * - ``SofteningBndryMaxPhys``
     - Specifies the maximum physical softening of the sixth particle group
       in Gadget.

Units Group
-----------

This group describes the unit system used throughout the file.
Attributes should be included for length, mass and velocity units. In
each case, three attributes are required to describe the units used (in
the following ``quantity`` refers to ``length``,
``mass``, ``time`` or ``velocity``):

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Attribute
     - Description
   * - ``quantityUnitsInSI``
     - The units of this quantity expressed in the SI system;
   * - ``quantityHubbleExponent``
     - The exponent of the reduced Hubble constant, :math:`h`, appearing in the
       units for this quantity;
   * - ``quantityScaleFactorExponent``
     - The exponent, :math:`n`, of the expansion factor, :math:`a`, required to convert
       this quantity into physical units. That is, multiplying this quantity by :math:`a^n` will give the quantity in physical units.

For example, if lengths in the file are expressed in units of comoving
:math:`h^{-1}` Mpc, then we would have

.. code-block:: text

        lengthUnitsInSI           =  3.08568e+22
        lengthHubbleExponent      = -1
        lengthScaleFactorExponent =  1

This allows a code reading the data from a merger tree file to
automatically convert it into whatever unit/coordinate system it
chooses.

An example of this group, showing standard attributes, is given below.

.. code-block:: text

      ATTRIBUTE "lengthHubbleExponent" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): -1
         }
      }
      ATTRIBUTE "lengthScaleFactorExponent" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 1
         }
      }
      ATTRIBUTE "lengthUnitsInSI" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 3.08568e+22
         }
      }
      ATTRIBUTE "massHubbleExponent" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): -1
         }
      }
      ATTRIBUTE "massScaleFactorExponent" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "massUnitsInSI" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 1.98892e+40
         }
      }
      ATTRIBUTE "timeHubbleExponent" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "timeScaleFactorExponent" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "timeUnitsInSI" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 3.1556926e+16
         }
      }
   }   ATTRIBUTE "velocityHubbleExponent" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "velocityScaleFactorExponent" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "velocityUnitsInSI" {
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SCALAR
         DATA {
         (0): 1000
         }
      }
   }

Forest Halos Group
------------------

The ``forestHalos`` group contains the data describing the
actual merger tree forests. Nodes from each forest must be stored
contiguously. An example of this group is given below. In this example,
``<nodeCount>`` is the total number of nodes in all merger
tree forests.

.. code-block:: text

   GROUP "forestHalos" {
      ATTRIBUTE "haloMassesIncludeSubhalos" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "haloAngularMomentaIncludeSubhalos" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "forestsAreSelfContained" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 1
         }
      }
      ATTRIBUTE "treesHaveSubhalos" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 1
         }
      }
      ATTRIBUTE "velocitiesIncludeHubbleFlow" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      ATTRIBUTE "positionsArePeriodic" {
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SCALAR
         DATA {
         (0): 0
         }
      }
      DATASET "descendantIndex" {
      COMMENT "The index of each descendant node."
         DATATYPE  H5T_STD_I64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "expansionFactor" {
      COMMENT "The expansion factor of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "halfMassRadius" {
      COMMENT "The half mass radius of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "hostIndex" {
      COMMENT "The index of each host node."
         DATATYPE  H5T_STD_I64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "nodeIndex" {
      COMMENT "The index of each node."
         DATATYPE  H5T_STD_I64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "nodeMass" {
      COMMENT "The mass of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "particleCount" {
      COMMENT "The number of entries within the particles group for this node."
         DATATYPE  H5T_STD_I64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "particleStart" {
      COMMENT "The index within the particles group at which the particle data for this node is stored."
         DATATYPE  H5T_STD_I64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "position" {
      COMMENT "The position of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount>, 3 ) / ( <nodeCount>, 3 ) }
      }
      DATASET "scaleRadius" {
      COMMENT "The scale radius of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "redshift" {
      COMMENT "The redshift of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "spin" {
      COMMENT "The spin of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount>, 3 ) / ( <nodeCount>, 3 ) }
      }
      DATASET "time" {
      COMMENT "The time of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount> ) / ( <nodeCount> ) }
      }
      DATASET "velocity" {
      COMMENT "The velocity of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <nodeCount>, 3 ) / ( <nodeCount>, 3 ) }
      }
   }

Standard Attributes
~~~~~~~~~~~~~~~~~~~~

The following are standard attributes in the ``forestHalos``
group (others may be added as desired).

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Attribute
     - Description
   * - ``haloMassesIncludeSubhalos``
     - Indicates whether or not the masses of halos include the masses of
       any subhalos that they may contain. A value of 0 implies that halo
       masses do not include masses of subhalos, while a value of 1
       indicates that they do;
   * - ``haloAngularMomentaIncludeSubhalos``
     - Indicates whether or not the angular momenta of halos include the
       angular momenta contributed by any subhalos that they may contain. A
       value of 0 implies that halo masses do not include the contribution
       of subhalos, while a value of 1 indicates that they do. If this
       attribute is not present it is assumed to be equal to the
       ``haloMassesIncludeSubhalos`` attribute;
   * - ``forestsAreSelfContained``
     - Indicates whether or not forests are self-contained, in the sense
       that nodes never transfer from one forest to another. A value of 0
       implies that nodes *can* move from one tree to another, while a
       value of 1 implies that they can not;
   * - ``treesHaveSubhalos``
     - Indicates whether or not trees contain information on subhalos. A
       value of 0 implies that they do not, while a value of 1 implies that
       they do. This attribute is a convenience, as subhalo presence can be
       determined from the node data directly;
   * - ``velocitiesIncludeHubbleFlow``
     - Indicates whether or not velocities include the Hubble flow. A value
       of 0 indicates that they do not, while a value of 1 indicates that
       they do. See `here <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ for important
       notes on velocity definitions in Galacticus.
   * - ``positionsArePeriodic``
     - Indicates whether or not positions are periodic (as in a
       cosmological cube simulation). A value of 0 indicates that they are
       not, while a value of 1 indicates that they are periodic, with a
       period of ``boxSize``.

Standard Datasets
~~~~~~~~~~~~~~~~~

The following are standard datasets in the ``forestHalos``
group.

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Dataset
     - Description
   * - ``angularMomentum``
     - The angular momentum of the halo. This can be either the magnitude
       of the angular momentum or a 3-D vector;
   * - ``expansionFactor``
     - The expansion factor (normalized to unity at the present day) at
       which this node is identified (note that only one of the
       ``expansionFactor``, ``redshift`` and
       ``time`` datasets is required, since they are simply
       related, but multiple *can* be present);
   * - ``descendantIndex``
     - The ``nodeIndex`` of the descendant of this node in the
       merger tree, or -1 if there is no descendant;
   * - ``halfMassRadius``
     - The radius containing half the mass of the node;
   * - ``hostIndex``
     - The ``nodeIndex`` of the node which hosts this node. For
       nodes that are self-hosting (i.e. that are not subhalos inside
       another halo), the value of ``hostIndex`` should be set
       equal to the node’s own ``nodeIndex``;
   * - ``nodeIndex``
     - An ID number for the node, unique at least within each tree. If
       nodes are able to move from one tree to another, the ID must be
       unique within all trees. No other constraints are placed on
       ``nodeIndex`` (e.g. it *does not* have to be monotonically
       increasing within the file for example);
   * - ``nodeMass``
     - The mass of the node;
   * - ``position``
     - The three dimensional position of this node;
   * - ``redshift``
     - The redshift at which this node is identified (note that only one of
       the ``expansionFactor``, ``redshift`` and
       ``time`` datasets is required, since they are simply
       related, but multiple *can* be present);
   * - ``scaleRadius``
     - The characteristic scale radius in the node (typically, but not
       necessarily, the NFW scale radius);
   * - ``spin``
     - The spin parameter, :math:`\lambda`, of the halo. This can be either the
       spin magnitude or a 3-D vector;
   * - ``time``
     - The time at which this node is identified (note that only one of the
       ``expansionFactor``, ``redshift`` and
       ``time`` datasets is required, since they are simply
       related, but multiple *can* be present);
   * - ``velocity``
     - The three dimensional velocity of the node. See
       `here <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ for important notes on
       velocity definitions in Galacticus;
   * - ``particleIndexStart``
     - If the ``particles`` group is included, this dataset should
       give the index of the first entry in the particle datasets that
       corresponds to the particle associated with this node;
   * - ``particleIndexCount``
     - If the ``particles`` group is included, this dataset should
       give the number of entries in particle datasets that correspond to
       the particle associated with this node;

Note that, internally, Galacticus uses a different naming convention for the links
between nodes. Specifically, an isolated node’s descendant is known as
its parent, while a satellite node’s host is also referred to as the
parent. By default, Galacticus will output a ``parentIndex`` dataset,
which therefore specifies either the descendant or host, depending on
whether the node in question is isolated or not. To additionally output
information which matches the use of “descendant” and “host” in these
merger tree files, set both of the input parameters and to
``true``. This will result in the output of two additional
datasets, ``descendantIndex`` and ``hostIndex`` which
correspond to the definitions used in the merger tree file.

Forest Index Group
------------------

The ``forestIndex`` group contains indexing information which
describes which sections of the datasets in the ``forestHalos``
group belong to each forest. An example of this group is given below.

.. code-block:: text

   GROUP "forestIndex" {
      DATASET "firstNode" {
      COMMENT "Position of the first node in each forest in the forestHalos datasets."
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SIMPLE { ( <forestCount> ) / ( <forestCount> ) }
      }
      DATASET "numberOfNodes" {
      COMMENT "Number of nodes in each forst."
         DATATYPE  H5T_STD_I32LE
         DATASPACE  SIMPLE { ( <forestCount> ) / ( <forestCount> ) }
      }
      DATASET "forestIndex" {
      COMMENT "Unique index of forst."
         DATATYPE  H5T_STD_I64LE
         DATASPACE  SIMPLE { ( <forestCount> ) / ( <forestCount> ) }
      }
      DATASET "forestWeight" {
      COMMENT "The number of such forests required per unit volume to create a representative sample."
         DATATYPE  H5T_STD_F64LE
         DATASPACE  SIMPLE { ( <forestCount> ) / ( <forestCount> ) }
      }
   }

Standard Datasets
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Dataset
     - Description
   * - ``firstNode``
     - For each forest, gives the position [1]_ in the
       ``forestHalos`` datasets of the first node in the forest
       (note that dataset indexing begins at 0). This *does not*
       necessarily have to be a root node of the forest—nodes in a single
       forest can be stored in any order in the ``forestHalos``
       datasets, providing that they are contiguous;
   * - ``numberOfNodes``
     - For each forest, gives the number of nodes in the forest;
   * - ``forestIndex``
     - A unique ID number for each forest;
   * - ``forestWeight``
     - A weight factor specifying the number density of each forest
       required to construct a representative sample. If not present, it is
       acceptable to assume that the weight is 1/``boxSize``\ :sup:`3`
       if that attribute is present in the ``simulation`` group. In cases where a
       simulation has been split into multiple subvolumes, the ``forestWeight``
       should still be the weight required to construct a representative sample
       from the *entire* volume. That is, for a cubic simulation of side length
       :math:`L` split into :math:`N` equal-volume subvolumes, ``forestWeight`` should be
       :math:`1/L^3` in each subvolume's merger tree file, *not* :math:`N/L^3`.

Forests Group
-------------

The ``forests`` group is optional and provides a convenience
method for accessing the properties of individual forests. If present,
it contains one group for each tree in the file, named
``forest<forestID>`` where ``<forestID>`` is
the ID of the forest. Each of these groups should contain a set of
scalar references to the sections of the datasets in the
``forestHalos`` group to which this forest corresponds. For
example, the ``descendantIndex`` reference for the forest with
ID number 89 would be as follows:

.. code-block:: text

   GROUP "forests/forest89" {
      DATASET "descendantIndex" {
         DATATYPE  H5T_REFERENCE
         DATASPACE  SCALAR
         DATA {
            DATASET /forestHalos/descendantIndex {(<indexBegin>)-(<indexEnd>)}
         }
      }
   }

Particles Group
---------------

The ``particles`` group is optional and contains information on
particle trajectories. It is intended for use with merger trees derived
from N-body simulations for which it is often useful to track the
location of, for example, the most bound particle associated with a
subhalo even after that subhalo can no longer be tracked in the
simulation. An example of this group is given below. In this example,
``<particleCount>`` is the total number of particles included
in the group.

.. code-block:: text

   GROUP "particles" {
      DATASET "particleID" {
      COMMENT "The ID of each particle."
         DATATYPE  H5T_STD_I64LE
         DATASPACE  SIMPLE { ( <particleCount> ) / ( <particleCount> ) }
      }
      DATASET "redshift" {
      COMMENT "The redshift of each particle."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <particleCount> ) / ( <particleCount> ) }
      }
      DATASET "time" {
      COMMENT "The time of each particle."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <particleCount> ) / ( <particleCount> ) }
      }
      DATASET "expansionFactor" {
      COMMENT "The expansion factor of each particle."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <particleCount> ) / ( <particleCount> ) }
      }
      DATASET "position" {
      COMMENT "The position of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <particleCount>, 3 ) / ( <particleCount>, 3 ) }
      }
      DATASET "velocity" {
      COMMENT "The velocity of each node."
         DATATYPE  H5T_IEEE_F64LE
         DATASPACE  SIMPLE { ( <particleCount>, 3 ) / ( <particleCount>, 3 ) }
      }
   }

Each particle should be stored contiguously (i.e. entries with the same
``particleID`` should be consecutive) and it is frequently
convenient (although not required) that entries for each particle be
arranged in order of increasing cosmic time.

Standard Datasets
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: auto

   * - Dataset
     - Description
   * - ``particleID``
     - An ID, unique within the entire simulation, for each particle;
   * - ``redshift``
     - The redshift at which the particle is recorded (a single particle
       can appear in these datasets at multiple times);
   * - ``time``
     - The time at which the particle is recorded (a single particle can
       appear in these datasets multiple times at different redshifts);
   * - ``expansionFactor``
     - The expansion factor at which the particle is recorded (a single
       particle can appear in these datasets multiple times at different
       expansion factors);
   * - ``position``
     - The spatial position of the particle;
   * - ``velocity``
     - The velocity of the particle. See
       `here <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ for important notes on
       velocity definitions in Galacticus.

Note that only one of the ``expansionFactor``,
``time`` and ``redshift`` datasets is required, as
they are simply related, but multiple of them *can* be present.

.. [1] That is, it gives the array index of the first node of the forest
   in the ``forestHalos/nodeIndex`` dataset for example. It
   does *not* give the ``nodeIndex`` of the first node in the
   forest.
