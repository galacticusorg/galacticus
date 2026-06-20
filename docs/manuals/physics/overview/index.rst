How Galacticus Works
====================

Galacticus evolves halos and galaxies using its "*evolver engine*", which works by applying various "*operators*" to the "*components*" which make up a halo/galaxy system. Here's a schematic - click (you may have to middle-click to make this work) on any item for more details.

.. mermaid::

   flowchart LR
     Components([Components])
     Operators
     Functions
     Engine[[Engine]]
     Engine --> Components
     Components --> Engine
     Operators --> Engine
     Functions --> Operators
     click Components href "https://github.com/galacticusorg/galacticus/wiki/How-Galacticus-Evolves-Halos-and-Galaxies#components" "Things that make up a galaxy (disk, spheroid, etc.)"
     click Engine href "https://github.com/galacticusorg/galacticus/wiki/How-Galacticus-Evolves-Halos-and-Galaxies#engine" "Applies operators to components to evolve them forward in time"
     click Operators href "https://github.com/galacticusorg/galacticus/wiki/How-Galacticus-Evolves-Halos-and-Galaxies#operators" "Processes that act on components"
     click Functions href "https://github.com/galacticusorg/galacticus/wiki/How-Galacticus-Evolves-Halos-and-Galaxies#functions" "Define the functional forms of processes"
    style Operators fill:#ff02ff
    style Components fill:#ff1212
    style Functions fill:#1212ff,color:#ffffff

Consider the following, highly simplified, description of how the process of star formation causes the mass of stars, :math:`M_\star`, to increase with time, while simultaneously decreasing the mass of the ISM, :math:`M_\mathrm{ISM}`:

.. math::

   \frac{\color{magenta}\mathrm{d}\color{red}M_\star}{\color{magenta}\mathrm{d}t} = \color{magenta}+\color{blue}\frac{M_\mathrm{ISM}}{\tau_\star}

.. math::

   \frac{\color{magenta}\mathrm{d}\color{red}M_\mathrm{ISM}}{\color{magenta}\mathrm{d}t} = \color{magenta}-\color{blue}\frac{M_\mathrm{ISM}}{\tau_\star}

Components
----------

The variables highlighted in red represent physical quantities - in this case the masses of stars and ISM in the galaxy. These are provided by a "`component <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_" (e.g. a galactic disk) in Galacticus.

Operators
---------

The operators highlighted in magenta represent a physical process - in this case the process of star formation, which moves mass from the ISM to the stars. Physical processes in Galacticus are implemented by the `nodeOperatorClass <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_.

Functions
---------

The function highlighted in blue represents the actual physics of that process - in this case it describes the rate of star formation. Such functions in Galacticus are provided by numerous different `functionClass <https://galacticus.readthedocs.io/en/latest/manuals/developer-guide/index.html>`_ objects - e.g. `starFormationRateDisksClass <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ in the case of star formation rates in galaxy disks.

Engine
------

Galacticus' evolver engine works by applying a set of `nodeOperatorClass <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ objects to each node in a merger tree in turn - gradually evolving the components in that node forward in time in accordance with the physical processes described by those `nodeOperatorClass <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ objects. The engine consists of the `mergerTreeEvolverClass <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ and `mergerTreeNodeEvolverClass <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ classes.

.. toctree::
   :maxdepth: 1

   structure-formation
   merger-tree-building
   galaxy
   cgm
   cgm-cooling
   star-formation
   galactic-structure
   outflows
   black-holes
   subhalo-evolution
