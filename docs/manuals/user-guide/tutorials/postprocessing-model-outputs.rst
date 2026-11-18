Postprocessing Model Outputs
============================

Once a model has been run, you may want properties of its galaxies which it did not output: spectra over a different range of
wavelengths, emission lines, or luminosities attenuated by a different dust model. Running the model again to get them is
wasteful, and for a large lightcone it may be impractical. Instead, Galacticus can rebuild each galaxy from what was written to
the output file, and compute new properties for it. This is done by the
:galacticus-class:`taskPostprocessOutputs` task, and most easily through the ``galacticus postprocess`` command of the Python
launcher.

What is rebuilt
---------------

Only what was output can be rebuilt. For each galaxy, the following are restored, where datasets for them are present in the
output file:

* the time (which must be output for a lightcone, since each galaxy in a lightcone has its own time), and the mass of the
  basic component;
* the stellar and gas masses, masses of metals in stars and gas, radius, angular momentum, and velocity of the disk and
  spheroid components (only the total mass of metals, not that of individual elements);
* the star formation histories of the disk and spheroid components;
* the inclination of the galaxy.

This is enough for any property computed from the star formation history and from the sizes and gas content of the galaxy:
spectra and luminosities (:galacticus-class:`nodePropertyExtractorSED`), emission lines
(:galacticus-class:`nodePropertyExtractorLuminosityEmissionLine`), and their attenuation by dust
(:galacticus-class:`nodePropertyExtractorDustAttenuation`) using, for example,
:galacticus-class:`dustAttenuationStellarMassRedshift`, or the inclination-dependent tabulations of the dust compendium
(:galacticus-class:`dustAttenuationAtlasCompendium`).

Galaxies are rebuilt in isolation: no host halo, satellite status, or position in a merger tree is restored. A property which
depends on any of these can not be computed (for example, the time at which a galaxy was last isolated is reported as its
current time). Nor can a property which depends on something not listed above: the emission lines of
:galacticus-class:`nodePropertyExtractorLuminosityEmissionLinePanuzzo2003`, for example, require stellar luminosities which are
not restored.

Preparing a model for postprocessing
------------------------------------

The model whose output is to be postprocessed must use the :galacticus-class:`starFormationHistoryFixedAges` star formation
history class (as is usual for lightcones), and output the time of each galaxy and the star formation history of each
component:

.. code-block:: xml

   <starFormationHistory value="fixedAges">
     <ageMinimum         value="0.001"/>
     <countAges          value="10"   />
     <countMetallicities value=" 5"   />
   </starFormationHistory>

   <nodePropertyExtractor value="multi">
     <nodePropertyExtractor value="time"/>
     <nodePropertyExtractor value="starFormationHistory">
       <component value="disk"/>
     </nodePropertyExtractor>
     <nodePropertyExtractor value="starFormationHistory">
       <component value="spheroid"/>
     </nodePropertyExtractor>
     <!-- ...any other properties... -->
   </nodePropertyExtractor>

The properties of the disk and spheroid components are output by default. If inclination-dependent properties are wanted,
give each galaxy a random inclination by adding

.. code-block:: xml

   <nodeOperator value="inclinationRandom"/>

to the model's node operators, and output it with

.. code-block:: xml

   <nodePropertyExtractor value="inclination">
     <galacticInclination value="random"/>
   </nodePropertyExtractor>

The inclination is written in degrees, and is read back as such. (Inclinations generated some other way and added to the file
must therefore also be in degrees, with zero being face-on.) If the file has no ``inclination`` dataset, every galaxy is
postprocessed face-on.

Postprocessing from Python
--------------------------

The ``galacticus postprocess`` command reconstructs the parameters of the original model from the ``Parameters`` group of its
output file, so that the star formation history, stellar populations, and node components are configured exactly as they were.
It then replaces the task, the outputter, the galactic filter, and the property extractors, and runs Galacticus. All that need
be given are the output file, the properties wanted, and the name of the file to write. For example, given a file
``lines.xml`` containing

.. code-block:: xml

   <nodePropertyExtractor value="dustAttenuation">
     <outputUnattenuated value="true"/>
     <dustAttenuation value="atlasCompendium">
       <fileName            value="compendium_exp_sech_Hernquist_hd0.137_hz0.137_dustD03Rv3.1_Attenuations_v2.hdf5"/>
       <galacticInclination value="random"/>
     </dustAttenuation>
     <nodePropertyExtractor value="luminosityEmissionLine">
       <component value="disk"/>
       <lineNames value="balmerAlpha6565 oxygenIII5008"/>
     </nodePropertyExtractor>
   </nodePropertyExtractor>

the disk emission lines of every galaxy in a lightcone, with and without attenuation by dust, are computed with:

.. code-block:: none

   galacticus postprocess lightcone.hdf5 lines.xml -o lines.hdf5 --outputs-group Lightcone

(:galacticus-class:`dustAttenuationAtlasCompendium` checks that the opacity of the dust, set by its ``dustProperties``, matches
that with which the tabulation was computed, and if not reports the value to give as ``opacityExtinctionV``.)

Change files given after the name of the extractors file are applied to the parameters of the original model, so that, for
example, the stellar population library can be changed. ``galacticus postprocess --help`` lists the other options.

The same is available from Python, where the galaxies to postprocess can be selected directly:

.. code-block:: python

   import numpy as np
   from galacticus_launcher.postprocess import postprocess

   # Postprocess galaxies 10, 3, and 7 of the first (and only) lightcone output.
   postprocess("lightcone.hdf5", "lines.xml", "lines.hdf5", outputs_group="Lightcone", selection={1: np.array([10, 3, 7])},
               threads=8)

The selection maps the index of each output (the ``i`` of ``Output[i]``) to the indices, counting from zero, in its
``nodeData`` datasets of the galaxies to postprocess. Indices may be in any order, and may repeat---a repeated galaxy is
postprocessed, and output, each time. Outputs which are not listed are skipped. The extractors may be given as the name of a
file, a string of XML, an ``lxml`` element, or a list of these. ``galacticus_launcher.postprocess.build_parameters()`` and
``galacticus_launcher.postprocess.write_selection()`` give the parameter file and the selection file without running
Galacticus.

Parameters which took their default values in the original model are not carried over: the executable which postprocesses
applies its own defaults afresh. (Some defaults are computed from other parameters, and can not be given explicitly alongside
them.) If that executable is of a different version than the one which ran the model, check that the defaults of interest have
not changed.

Reading the results
-------------------

The output file has the same structure as any Galacticus output, with one row per galaxy postprocessed in each
``[outputsGroupName]/Output[i]/nodeData`` group. Galaxies are *not* written in the order in which they were selected, since
they are processed in parallel. The dataset ``nodeIndexPostprocessSource`` gives, for each row, the index of the galaxy in the
``nodeData`` datasets of the original file, and should be used to match rows to the original galaxies:

.. code-block:: python

   import h5py

   original      = h5py.File("lightcone.hdf5", "r")["Lightcone/Output1/nodeData"]
   postprocessed = h5py.File("lines.hdf5"    , "r")["Lightcone/Output1/nodeData"]
   indices       = postprocessed["nodeIndexPostprocessSource"][:]
   redshift      = original["lightconeRedshiftObserved"][:][indices]
   halpha        = postprocessed["luminosityEmissionLineDisk:balmerAlpha6565:dustAttenuated:atlasCompendium"][:]

Only the properties requested are output: the properties of the rebuilt components would merely duplicate (or, for properties
which are not restored, misrepresent) those in the original file.

Using the task directly
-----------------------

The ``galacticus postprocess`` command writes an ordinary parameter file (``--parameter-file`` keeps a copy). To write one by
hand, start from the parameter file of the original model (or extract it from the output file with
``scripts/parameters/parametersExtract.py``), and replace the task, outputter, and extractors:

.. code-block:: xml

   <task value="postprocessOutputs">
     <fileName          value="lightcone.hdf5"/>
     <outputsGroupName  value="Lightcone"     />
     <!-- Optional: an HDF5 file containing a `Lightcone/Output1/selection` dataset of galaxy indices. -->
     <selectionFileName value="selection.hdf5"/>
     <countBatch        value="1024"          />
   </task>

   <mergerTreeOutputter value="standard">
     <outputsGroupName          value="Lightcone"/>
     <outputComponentProperties value="false"    />
     <countNodesBuffered        value="1000"     />
   </mergerTreeOutputter>

   <galacticFilter value="always"/>

   <nodePropertyExtractor value="multi">
     <nodePropertyExtractor value="indexPostprocessSource"/>
     <!-- ...the properties to compute... -->
   </nodePropertyExtractor>

   <outputFileName value="postprocessed.hdf5"/>

The galactic filter of the original model should usually be replaced, since it was applied when the file was written and may
need properties (positions within a lightcone, for example) which are not restored. Buffering the output of 1000 or so galaxies
at a time (``countNodesBuffered``) greatly reduces the time spent writing it, and allows postprocessing to scale well across
OpenMP threads; without it, writing can dominate the run time.
