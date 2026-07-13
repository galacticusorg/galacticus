Advanced Usage
==============

In this chapter we cover in more detail several aspects of Galacticus, including the structure of parameter files, the structure of output files, and how to perform various types of analysis.

.. _manual-sec-ParameterFiles:

Parameter Files
---------------

As described above, Galacticus requires a file of parameters to be given as a command line argument. The parameter file is an :term:`XML` file (which makes it easy to manipulate and construct these files from within many languages, e.g. Python) with the following structure:

.. code-block:: none

    <parameters>
      <version>0.9.4</version>
      <formatVersion>2</formatVersion>
      <parameter1Name value= "parameter1Value" />
      <parameter2Name>
        <value>parameter2Value</value>
      </parameter2Name>
      <parameter3Name value= "parameter3Value" >
         <subParameter1Name value= "subParameter1Value" />
         <subParameter2Name value= "subParameter2Value" />
         .
         .
         .
      </parameter3Name>
      .
      .
      .
      <parameter4Name value= "parameter4Value" id="myRefParam" >
         <subParameter1Name value= "subParameter1Value" />
         <subParameter2Name value= "subParameter2Value" />
         .
         .
         .
      </parameter3Name>
      .
      .
      .
      <parameter5Name value= "parameter5Value" >
         <parameter4Name idRef="myRefParam"/>
         .
         .
         .
      </parameter5Name>
      .
      .
      .
      <parameter6Name value="=10.0*[parameter4Name/subParameter1Name]"/>
    </parameters>

Each named element must have a ``value`` attribute (preferred), or else contains a value element, which contains the desired value. The value can be a number, word or an array of space-separated numbers or words. Parameters are used to control the values of numerical parameters and also to select methods and other options. In many cases, if a parameter is not specified in the file a default value (hard coded into Galacticus) will be used instead. The default values have been chosen to produce a realistic model of galaxy formation, but may change as Galacticus evolves. Parameters may have sub-parameters embedded within them, as in the example above.

Sub-parameters are used in object composition within Galacticus. For example, the following would specify that linear growth of cosmological large scale structure should be modeled using the ``collisionlessMatter`` method:

.. code-block:: none

     <linearGrowth value="collisionlessMatter">
       <cosmologyParameters value="simple">
         <HubbleConstant   value="70.0" />
         <OmegaMatter      value="0.31"  />
         <OmegaDarkEnergy  value="0.69"  />
         <OmegaBaryon      value="0.045"/>
         <temperatureCMB   value="2.725"/>
       </cosmologyParameters>
       <cosmologyFunctions value="matterLambda"/>
     </linearGrowth>

The linear growth function object requires knowledge about the cosmological parameters and model. In the above, we specify this explicitly by including a definition of the cosmological parameter object and cosmological functions object that our linear growth function object should use. Note that the cosmological functions object also requires knowledge of the cosmological parameters. When the ``cosmologyFunctions`` object is built from the above definition it will first check for a cosmological parameters object defined in its own subparameters. Since it does not find one in this instance it will check for a cosmological parameters definition in its parent object (the ``linearGrowth`` element) and, in this case, will use that definition. If no definition were to be found in any parent element, a default set of cosmological parameters would be used instead\ [#]_.

Parameters can also be defined with an "``id``" attribute. Such parameters are targets for pointers elsewhere in the file (but also act as regular parameters and will be utilized as such by Galacticus). A pointer to a target is created by specifying an element with the same parameter name and an "``idRef``" element with a value equal to that of the ``id`` element of the target. The pointer then acts as a regular parameter, adopting the value and subparameters of the target. Targets can be targeted by multiple pointers. This allows for a single object to be shared between multiple other objects.

It is possible to make a parameter conditional upon the value of another parameter. This can be useful, for example, if you want to include a specific :galacticus-class:`nodeOperatorClass` only if merger trees are being read from file, and not if trees are being constructed internally. This is achieved using the ``active`` attribute of a parameter. For example:

.. code-block:: none

       <nodeOperator value="assemblyHistoryHeuristics" active="[mergerTreeConstructor] == read">

which will result in the ``assemblyHistoryHeuristics`` ``nodeOperator`` being active only if the ``mergerTreeConstructor`` is set to ``read``. The content of the ``active`` element must be a parameter name enclosed in ``[]``, followed by an operator (``==`` or ``!=``), followed by the text to match. The parameter name can include a full path as described below. Note that, while numerical parameters can be used in ``active`` elements, the comparison is always textual, not numerical.

A parameter value can be given by a math expression\ [#]_. In the above example ``parameter5Name`` has a value of "``=10.0*[parameter4Name:subParameter1Name]``". The initial "``=``" indicates that this is a math expression which should be evaluated. In this case, the expression is ten times value value of "``[parameter4Name:subParameter1Name]``", which refers to the value of the sub-parameter ``subParameter1Name`` of parameter ``parameter4Name``. In this way the values of parameters can be derived from those of other parameters.

Lastly, a text parameter can also be evaluated by starting it with an initial "``=``". For text parameters, other parameters can be referenced and their values are simply inserted. So, for example:

.. code-block:: none

   <outputFileName value="=myFile_m[%5.2f|darkMatterParticle/mass]keV.hdf5"/>

would set the output file name by replacing "``[%5.2f|darkMatterParticle/mass]``" in the above with the value of the dark matter particle mass. Note that, for these text parameter substitutions a format specifier must be given---in this case "``%5.2f``"---which specifies how to format the parameter before inserting it into the value. The `format specifier <https://www.w3resource.com/c-programming/stdio/c_library_method_sprintf.php>`_ should follow standard C conventions. Currently the ``s``, ``d``, ``f``, and ``e`` specifiers are supported.

The optional ``version`` element specifies which version of Galacticus this parameter file is intended for. The optional ``formatVersion`` element specifies the parameter file version number (the current standard for parameter files is version 2). While optional, these elements can be useful when migrating parameter files between versions of Galacticus.

All parameter values (both those specified in this file and those set to default) used during a Galacticus run are output to the ``Parameters`` group within the Galacticus output file. If parameters are present in the parameter file which do not match any known parameter in Galacticus then a warning message, listing all unknown parameters, will be given when Galacticus is run. Note that this will *not* prevent Galacticus from running---sometimes it is convenient to include parameters which are not used by Galacticus, but which might be used by some other code.

Extracting Parameter Files From Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The values of all parameters (including those set to defaults) are stored to the Galacticus output file in the ``Parameters`` group (see Section :galacticus-ref:`outputFile:parametersGroup`). These parameter settings can be extracted back to an XML file using the ``parametersExtract.py`` script:

.. code-block:: none

   ./scripts/parameters/parametersExtract.py galacticus.hdf5 extractedParameters.xml

In this example, all parameters that were used to run the ``galacticus.hdf5`` model and that were stored in that file will be extracted and output to the ``extractedParameters.xml`` file.

Differencing Parameter Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The differences between two parameter files can be shown using the ``parametersDiff.py`` script. For example:

.. code-block:: none

   ./scripts/parameters/parametersDiff.py parameters1.xml parameters2.xml

will show the differences between ``parameters1.xml`` and ``parameters2.xml``.

By default, the *order* of differently-named parameters is ignored when looking for differences---the order of differently-named parameters makes no difference to Galacticus. However, this involves re-ordering parameters alphabetically to allow differences to be seen, which can make it more difficult for the user to identify where in the files the differences occur. By adding the option ``--respectOrder`` the order of parameters is preserved. This may result in more differences being shown, but with more useful context for finding them in the parameter files.

Differences are detected via a textual comparison. Consequently, parameters:

.. code-block:: none

   <myParameter value="0.001"/>

and

.. code-block:: none

   <myParameter value="1.0e-3"/>

will be identified as a difference, even though they are numerically identical. To avoid such false differences, numerical values in parameters can be put into a canonical form using the ``--canonicalizeValues`` option with a standard Python `format string <https://docs.python.org/3/library/string.html#formatstrings>`_. For example:

.. code-block:: none

   ./scripts/parameters/parametersDiff.py --canonicalizeValues .4f parameters1.xml parameters2.xml

will convert all numerical values into floating point numbers with 4 digits of provision. So, in the above example the parameters would be rewritten as:/

.. code-block:: none

   <myParameter value="0.0010"/>

and

.. code-block:: none

   <myParameter value="0.0010"/>

and so would be seen as identical.

.. _validating-parameter-files:

Validating Parameter Files
~~~~~~~~~~~~~~~~~~~~~~~~~~

A script, ``scripts/build/parameterValidate.py``, is provided to validate parameter files and thereby ensure that they are consistent with Galacticus's expectations and requirements. For a structural check of a single parameter file (no duplicate parameters, every parameter has a value, etc.) simply execute:

.. code-block:: none

    scripts/build/parameterValidate.py --structural myParameters.xml

No output for a finding (and an exit value of 0) indicates a valid parameter file. Problems are reported with an exit value other than 0 and an error message to help track down the issue.

The validator can additionally perform deeper, catalog-aware checks—verifying that each ``functionClass`` selector names a real implementation, that parameter names are accepted by the selected implementation, and that values respect the parameter's type, enumeration, and any range constraints. These require the parameter catalog (a machine-readable description of every implementation's parameters), generated by ``scripts/build/parameterCatalog.py``:

.. code-block:: none

    scripts/build/parameterCatalog.py `pwd` parameters.catalog.json
    scripts/build/parameterValidate.py --catalog parameters.catalog.json --structural myParameters.xml

A directory may be given in place of a file, in which case all ``*.xml`` files under it are validated and a summary is printed. This same tool is run in continuous integration to validate the bundled parameter files.

.. _editor-parameter-validation:

Editor validation and autocompletion
"""""""""""""""""""""""""""""""""""""

For real-time validation and autocompletion while editing parameter files, an XML Schema is provided at ``schema/parameters.xsd``. It is generated from the parameter catalog and constrains each ``functionClass`` selector's ``value`` to the valid implementations, and each enumeration-valued parameter's ``value`` to the allowed labels. (Content is otherwise permissive; for the precise, per-implementation check use ``parameterValidate.py`` above.)

**VS Code (whole repository).** A workspace configuration is shipped with Galacticus in ``.vscode/settings.json`` that associates the schema with every file under a ``parameters/`` directory. After installing the `XML extension <https://marketplace.visualstudio.com/items?itemName=redhat.vscode-xml>`_, validation and autocompletion therefore work out of the box when you open the Galacticus repository in VS Code. The shipped association is:

.. code-block:: json

    {
      "xml.fileAssociations": [
        { "pattern": "**/parameters/**/*.xml",
          "systemId": "${workspaceFolder}/schema/parameters.xsd" }
      ]
    }

The ``${workspaceFolder}``-anchored path is important: parameter files live at many directory depths, and a bare relative path would be looked for next to each file and not found. If ``${workspaceFolder}`` is not expanded by your version of the extension, use the absolute path to ``schema/parameters.xsd`` in your checkout instead. To cover parameter files kept outside any ``parameters/`` directory, broaden the ``pattern`` (e.g. to ``**/*.xml``).

**A single file (any XSD-aware editor).** Alternatively, reference the schema from the parameter file itself, on its root ``<parameters>`` element:

.. code-block:: xml

    <parameters xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                xsi:noNamespaceSchemaLocation="schema/parameters.xsd">

The location in ``xsi:noNamespaceSchemaLocation`` is resolved **relative to the parameter file's own directory** (not the repository root), so a bare ``schema/parameters.xsd`` only works for a file sitting directly above the ``schema/`` directory. For a file deeper in the tree, use a correctly-rooted relative path (e.g. ``../../schema/parameters.xsd`` for a file in ``parameters/tutorials/``) or an absolute path.

The schema is generated from the sources, so it is regenerated (``make parameters-schema``) and verified in continuous integration. The Galacticus :ref:`git hooks <git-hooks>` additionally check it before a commit is made — regenerating it (only when a staged change could affect it) and blocking the commit if ``schema/parameters.xsd`` is out of date.

Generating Parameter Files
~~~~~~~~~~~~~~~~~~~~~~~~~~

Some scripts are provided which assist in the generation of parameter files. These are located in the ``scripts/parameters/`` folder and are detailed below:

``cosmologicalParametersMonteCarlo.py``
   This script will generate a set of cosmological parameters drawn at random from the WMAP-9 constraints :cite:p:`hinshaw_nine-year_2012`. It uses the covariance matrix (currently defined in ``data/Cosmological_Parameters_WMAP-9.xml``) to produce correlated random variables\ [#]_. The generated parameters are printed to standard output as Galacticus-compatible XML.

Changing Parameter Files
~~~~~~~~~~~~~~~~~~~~~~~~

Galacticus allows parameter files to be modified before they are run, using a set of changes specified in one or more "change files". This can be useful to allow, for example, one base parameter file to be defined, and then to make small changes to allow different models/calculations to be performed (e.g. changing from running a fixed mass merger tree to simulate the Milky Way, to a distribution of merger tree masses to simulate an entire population of galaxies).

Change files are simple XML files with a syntax described below. They can be included on the Galacticus command line after the primary parameter file, e.g.:

.. code-block:: none

   ./Galacticus.exe parameters.xml changes1.xml changes2.xml

which will cause Galacticus to first read the ``parameters.xml`` file, then apply changes from ``changes1.xml`` to it, and then to apply changes from ``changes2.xml``.

An example change file can be found `here <https://raw.githubusercontent.com/galacticusorg/galacticus/master/testSuite/parameters/changes.xml>`_, and looks as follows:

.. code-block:: none

   <?xml version="1.0" encoding="UTF-8"?>
   <!-- Change a parameter file -->
   <changes>

     <!-- Append a parameter to the top-level section -->
     <change type="append" path="">
       <galacticFilter value="haloIsolated"/>
     </change>

     <!-- Insert a parameter before a specific parameter -->
     <change type="insertBefore" path="nodeOperator/nodeOperator[@value='barInstability']">
       <nodeOperator value="indexShift"/>
       <nodeOperator value="null"      />
     </change>

     <!-- Insert a parameter after a specific parameter -->
     <change type="insertAfter" path="nodeOperator/nodeOperator[@value='blackHolesSeed']">
       <nodeOperator value="indexBranchTip"/>
     </change>

     <!-- Insert a parameter after a specific parameter (the last parameter in this case) -->
     <change type="insertAfter" path="nodeOperator/nodeOperator[@value='blackHolesCGMHeating']">
       <nodeOperator value="indexLastHost"/>
     </change>

     <!-- Remove a parameter -->
     <change type="remove" path="nodeOperator/nodeOperator[@value='blackHolesWinds']">
     </change>

     <!-- Replace a parameter -->
     <change type="replace" path="nodeOperator/nodeOperator[@value='stellarFeedbackSpheroids']/stellarFeedbackOutflows/stellarFeedbackOutflows">
       <stellarFeedbackOutflows value="vlctyMxSclng">
         <fraction value="0.015"/>
         <exponentVelocity value="3.5"/>
       </stellarFeedbackOutflows>
     </change>

     <!-- Replace a parameter with a different parameter-->
     <change type="replace" path="nodeOperator/nodeOperator[@value='stellarFeedbackSpheroids']/stellarFeedbackOutflows/stellarFeedbackOutflows" target="nodeOperator/nodeOperator[@value='stellarFeedbackDisks']/stellarFeedbackOutflows/stellarFeedbackOutflows"/>

     <!-- Replace or append a parameter -->
     <change type="replaceOrAppend" path="starFormationRateSpheroids/starFormationTimescale/efficiency">
       <efficiency value="0.06"/>
     </change>

     <!-- Encapsulate a parameter -->
     <change type="encapsulate" path="cosmologicalMassVariance">
       <cosmologicalMassVariance value="scaled">
         <scale value="0.689"/>
       </cosmologicalMassVariance>
     </change>

     <!-- Update a parameter value -->
     <change type="update" path="nodeOperator/nodeOperator[@value='barInstability']/galacticDynamicsBarInstability/stabilityThresholdGaseous" value="0.75">
     </change>

   </changes>

A parameter change file should contain one or more change elements. Each change element must have a ``path`` attribute which gives an `XPath <http://www.w3schools.com/Xml/xpath_syntax.asp>`_ expression that identifies a parameter in the parameter file that will be acted upon by the change. (Note that currently only a limited set of XPath functionality is supported - matching of element names and single attributes only.) An empty path indicates the top-level section of the parameter file.

Each change element must also have a ``type`` attribute which specifies the type of change to be made. Supported types are:

``append``:
   This will append all parameters enclosed within the change element to the children of the parameter specified by the ``path`` attribute.

``insertBefore``:
   This will insert all parameters enclosed within the change element before the parameter specified by the ``path`` attribute.

``insertAfter``:
   This will insert all parameters enclosed within the change element after the parameter specified by the ``path`` attribute.

``replace``:
   This will replace the parameter specified by the ``path`` attribute with all parameters enclosed within the change element.

``replaceWith``:
   This will replace the parameter specified by the ``path`` attribute with the parameter (and all children) from the parameter file at the path specified by the ``target`` attribute (note that the parameter identified by the ``target`` attribute is not removed)

``replaceOrAppend``:
   This will replace the parameter specified by the ``path`` attribute, if it exists, with all parameters enclosed within the change element. If the parameter does not exist, instead all parameters enclosed within the change element will be appended to the children of the parent element.

``encapsulate``:
   This will take the parameter specified by the ``path`` attribute, replace it by all parameters enclosed within the change element, and then reinsert the parameter specified by the ``path`` attribute as a child of the first element enclosed within the change element. This allows an existing parameter to be encapsulated inside a new parameter.

``remove``:
   This will remove the parameter specified by the ``path`` attribute.

``update``:
   This will update the value of the parameter specified by the ``path`` attribute with that given by the ``value`` attribute of the change element.

.. _manual-sec-outputFile:

General Structure of Output File
--------------------------------

Figure :numref:`{number} <fig-glcOutputFileStructure>` shows the structure of a typical Galacticus output file. Galacticus can perform many different tasks, and output a wide variety of different data. The following describes the structure of an output file resulting from the "``evolveForests``" :galacticus-class:`taskClass` (i.e. the usual mode of running Galacticus where it is asked to form and evolve a population of galaxies within a set of merger trees) using the "``standard``" :galacticus-class:`mergerTreeOutputterClass`.

The various groups and subgroups are described below.

.. code-block:: none
   :name: fig-glcOutputFileStructure
   :caption: Structure of a Galacticus HDF5 output file. ``<treeCount>`` is the total number of merger trees present in a given output, and ``<nodeCount`` is the total number of nodes (in all trees) present in an output.

   galacticus.hdf5
    |
    +-> UUID                                     Attribute {1}
    |
    +-> Build                                    Group
    |    |
    |    +-> FoX_library_version                 Attribute {1}
    |    +-> GSL_library_version                 Attribute {1}
    |    +-> HDF5_library_version                Attribute {1}
    |    +-> make_CCOMPILER                      Attribute {1}
    |    +-> make_CCOMPILER_VERSION              Attribute {1}
    |    +-> make_CFLAGS                         Attribute {1}
    |    +-> make_CPPCOMPILER                    Attribute {1}
    |    +-> make_CPPCOMPILER_VERSION            Attribute {1}
    |    +-> make_CPPFLAGS                       Attribute {1}
    |    +-> make_FCCOMPILER                     Attribute {1}
    |    +-> make_FCCOMPILER_VERSION             Attribute {1}
    |    +-> make_FCFLAGS                        Attribute {1}
    |    +-> make_FCFLAGS_NOOPT                  Attribute {1}
    |    +-> make_MODULETYPE                     Attribute {1}
    |    +-> make_PREPROCESSOR                   Attribute {1}
    |    +-> sourceChangeSetBundle               Dataset   {1}
    |    +-> sourceChangeSetMerge                Dataset   {1}
    |
    +-> Outputs                                  Group
    |    |
    |    +-> Output1                             Group
    |    |    |
    |    |    +-> nodeData                       Group
    |    |    |     |
    |    |    |     +-> nodeProperty1            Dataset {<nodeCount>}
    |    |    |     +-> ...                      Dataset {<nodeCount>}
    |    |    |     +-> ...                      Dataset {<nodeCount>}
    |    |    |     +-> ...                      Dataset {<nodeCount>}
    |    |    |     +-> nodePropertyN            Dataset {<nodeCount>}
    |    |    |
    |    |    +-> mergerTreeCount                Dataset {<treeCount>}
    |    |    |
    |    |    +-> mergerTreeIndex                Dataset {<treeCount>}
    |    |    |
    |    |    +-> mergerTreeSeed                 Dataset {<treeCount>}
    |    |    |
    |    |    +-> mergerTreeStartIndex           Dataset {<treeCount>}
    |    |    |
    |    |    +-> mergerTreeWeight               Dataset {<treeCount>}
    |    |    |
    |    |    +-> mergerTree1                    Group              [optional]
    |    |    |     |
    |    |    |     +-> nodeProperty1            Reference
    |    |    |     +-> ...                      Reference
    |    |    |     +-> ...                      Reference
    |    |    |     +-> ...                      Reference
    |    |    |     +-> nodePropertyN            Reference
    |    |    |
    |    |    x-> ...                            Group              [optional]
    |    |    x-> ...                            Group              [optional]
    |    |    x-> ...                            Group              [optional]
    |    |    x-> mergerTree<treeCount>          Group              [optional]
    |    |    |
    |    |    +-> outputType                     Attribute {1}
    |    |    +-> outputExpansionFactor          Attribute {1}
    |    |    +-> outputTime                     Attribute {1}
    |    |
    |    x-> Output2                             Group
    |
    +-> Filters                                  Group
    |    |
    |    +-> name                                Dataset   {<filterCount>}
    |    +-> wavelengthEffective                 Dataset   {<filterCount>}
    |
    +-> Parameters                               Group
    |    |
    |    +-> inputParameter1                     Attribute {}
    |    +-> ...                                 Attribute {}
    |    +-> ...                                 Attribute {}
    |    +-> ...                                 Attribute {}
    |    +-> inputParameterN                     Attribute {}
    |    +-> inputParameter1                     Group
    |         |
    |         +-> subInputParameter1             Attribute {}
    |         +-> ...                            Attribute {}
    |         +-> subInputParameterN             Attribute {}
    |    x-> ...                                 Attribute {}
    |    x-> ...                                 Attribute {}
    |    x-> ...                                 Attribute {}
    |    x-> inputParameterN                     Group
    |
    +-> Version                                  Group
         |
         +-> buildTime                           Attribute {1}
         +-> memoryUsageMaximum                  Attribute {1}
         +-> runStartTime                        Attribute {1}
         +-> runEndTime                          Attribute {1}
         +-> runDuration                         Attribute {1}
         +-> gitBranch                           Attribute {1}
         +-> gitHash                             Attribute {1}
         +-> gitHashDatasets                     Attribute {1}
         +-> runByName                           Attribute {1}
         +-> runByEmail                          Attribute {1}

.. _manual-sec-UUID:

UUID
~~~~

The UUID (`Universally Unique Identifier <https://secure.wikimedia.org/wikipedia/en/wiki/Universally_unique_identifier>`_) is a unique identifier assigned to each Galacticus model that is run. It allows identification of a given model and can be referenced from, for example, an external database.

.. _manual-sec-BuildInformation:

Build Information
~~~~~~~~~~~~~~~~~

Galacticus automatically stores various information about how it was built in the ``Build`` group attributes. Currently, included attributes consist of:

``FoX_library_version``
   The version number of the FoX library;

``GSL_library_version``
   The version number of the GSL library;

``HDF5_library_version``
   The version number of the HDF5 library;

``make_CCOMPILER``
   The C compiler command used;

``make_CCOMPILER_VERSION``
   The C compiler version information;

``make_CFLAGS``
   The flags passed to the C compiler;

``make_CPPCOMPILER``
   The C++ compiler command used;

``make_CPPCOMPILER_VERSION``
   The C++ compiler version information;

``make_CPPFLAGS``
   The flags passed to the C++ compiler;

``make_FCCOMPILER``
   The Fortran compiler command used;

``make_FCCOMPILER_VERSION``
   The Fortran compiler version information;

``make_FCFLAGS``
   The flags passed to the Fortran compiler;

``make_FCFLAGS_NOOPT``
   The flags passed to the Fortran compiler for unoptimized compiles;

``make_MODULETYPE``
   The Fortran module type identifier string;

``make_PREPROCESSOR``
   The preprocessor command used.

Additionally, two datasets are included which store details of the Galacticus source changeset. ``sourceChangeSetBundle`` contains the output of "``git bundle create HEAD ôrigin``", that is, it contains a Git archive that incorporates any changes made to the current branch relative to the main Galacticus branch. ``sourceChangeSetDiff`` contains the output of "``git diff``", that is, all differences between the source code in the working directory and that which has been committed to Git. Used together, these two datasets allow the precise source code used to run the model to be recovered from the main branch Galacticus source.

Filters
~~~~~~~

For each broadband filter used in the Galacticus model run an entry is added to the datasets in this group. Currently, two datasets are generated:

``name``
   The name of each filter used.

``wavelengthEffective``
   The effective wavelength, :math:`\lambda_\mathrm{eff}` (defined as :math:`\lambda_\mathrm{eff}=\left. \int_0^\infty \lambda R(\lambda) \mathrm{d}\lambda \right/ \int_0^\infty R(\lambda) \mathrm{d}\lambda`, where :math:`R(\lambda)` is the filter response) of the filter in Å.

.. _manual-sec-outputFile-parametersGroup:

Parameters
~~~~~~~~~~

The ``Parameters`` group contains a record of all parameter values (either input or default) that were used for this Galacticus run. The group contains a long list of attributes, each attribute named for the corresponding parameter and with a single entry giving the value of that parameter. If a parameter has subparameters, a group is created having the same name as the parameter, which will contain attributes corresponding to each subparameter. In cases where a parameter appears more than once in a given node of the parameter tree,it will be output with "``[N]``" appended to its name, where "``N``" is an integer indicating the instance of the parameter.

Version
~~~~~~~

The ``Version`` group contains a record of the Galacticus version used for this model, storing the  Git commit branch and hash (if the code is being maintained using  Git, otherwise a value of "``unknown``" is entered) in the attributes ``gitBranch`` and ``gitHash`` respectively, along with the time at which the executable was built as ``buildTime``. If the ``datasets`` path is a  Git repo then the hash of the checked-out commit is stored as ``gitHashDatasets`` (if ``datasets`` is not a  Git repo then a value of "``unknown``" is entered instead). Additionally, the times at which the model run started and ended are stored as ``runStartTime`` and ``runEndTime``, with the duration of the run (in seconds) stored as ``runDuration``, and the maximum memory used by the model (in bytes) is stored as ``memoryUsageMaximum``. If the ``galacticusConfig.xml`` file (see Section :galacticus-ref:`ConfigFile`) is present and contains contact details, the name and e-mail address of the person who ran the model are stored as ``runByName`` and ``runByEmail`` respectively.

Outputs
~~~~~~~

The ``Outputs`` group contains one or more sub-groups corresponding to the output times requested from Galacticus. Each sub-group contains the following information:

``outputType`` *(attribute)*
   The type of the output---common types are ``snapshot`` (for outputs containing halos/galaxies at a single snapshot in time) and ``lightcone`` (for outputs containing halos/galaxies on an observer's past lightcone);

``outputTime`` *(attribute)*
   The cosmic time (in Gyr) at this output;

``outputExpansionFactor`` *(attribute)*
   The expansion factor at this output;

``nodeData``
   A group of node properties as described below.

``mergerTree`` subgroups *(optional)*
   A set of ``mergerTree`` groups as described below.

Output is controlled by parameters given within the ``mergerTreeOutput`` section of the parameter file. Current options are:

``outputMergerTrees``
   If ``true`` then each merger tree is output to the relevant sub-group at each output time (see Section :galacticus-ref:`nodeDataGroup`). Otherwise merger trees are not output. [Default: ``true``.]

``outputReferences``
   If ``true`` then an HDF5 reference dataset is written for each merger tree subgroup (see Section :galacticus-ref:`mergerTreeSubgroups`). [Default: ``false``.]

``galacticFilter``
   A :galacticus-class:`galacticFilterClass` which is applied to each node in the tree to determine whether or not it should be output. By combining multiple filters it is possible to construct arbitrarily complex criteria for output. [Default: ``always``.]

.. _manual-sec-nodeDataGroup:

nodeData group
^^^^^^^^^^^^^^

The ``nodeData`` group contains all data from nodes in all merger trees. The group consists of a collection of datasets each of which lists a property of all nodes in the trees which exist at the output time. Where relevant, each dataset contains an attribute, ``units``, which is a composite HDF5 data-type containing the units of the dataset in the SI system as ``unitsInSI``, a human-readable description of the units as ``description``, an ``astropy.units``-parseable description of the units as ``quantity``, and ``isComoving`` which is ``1`` if the units are in comoving coordinates, and ``0`` otherwise.

.. _manual-sec-mergerTreeDatasets:

mergerTree datasets
^^^^^^^^^^^^^^^^^^^

To allow locating of nodes belonging to a given merger tree in the datasets in the ``nodeData`` group, the ``mergerTreeStartIndex`` and ``mergerTreeCount`` datasets list the starting index of each tree's nodes in the ``nodeData`` datasets, and the number of nodes belonging to each tree respectively. Additionally, the ``mergerTreeWeight`` dataset lists the ``volumeWeight`` property for each tree (see Section :galacticus-ref:`mergerTreeSubgroups`) which gives the weight (in Mpc\ :math:`^{-3}`) which should be assigned to this tree (and all nodes in it) to create a volume-averaged sample (see Section :galacticus-ref:`volumeLimitedSamples`). The ``mergerTreeIndex`` dataset gives the index of each tree stored in the ``nodeData`` datasets. Finally, the ``mergerTreeSeed`` dataset gives the seed used by this tree when generating random numbers.

.. _manual-sec-mergerTreeSubgroups:

mergerTree subgroups
^^^^^^^^^^^^^^^^^^^^

These subgroups will be present if the ``[mergerTreeOutputReferences]`` parameter is set to true. Each ``mergerTree`` subgroup contains HDF5 references to all data on a single merger tree. The group consists of a collection of scalar references each of which points to the appropriate region of the corresponding dataset in the ``nodeData`` group. Additionally, the ``volumeWeight`` attribute of this group gives the weight (in Mpc\ :math:`^{-3}`) which should be assigned to this tree (and all nodes in it) to create a volume-averaged sample. (A second attribute, ``units``, gives the units of ``volumeWeight`` in the SI system, along with human-readable, ``astropy.units``-parseable descriptions, and a boolean indicating that these are in comoving coordinates.)

.. _manual-sec-onTheFlyAnalyses:

On-the-fly Analyses
-------------------

In addition to simply outputting the properties of every galaxy formed, Galacticus can perform various analyses as it runs, and output the resulting quantities. For example, it can construct a galaxy stellar mass function matched in terms of binning, redshift range, survey geometry, and stellar mass uncertainties to an observed stellar mass function---outputting the resulting model expectation, along with the observed result, and the likelihood of the model given the data.

To perform such on-the-fly analysis, simply use the "``analyzer``" :galacticus-class:`mergerTreeOutputterClass`, by including the following in your parameter file:

.. code-block:: none

    <mergerTreeOutputter value="analyzer"/>

If you want to also keep the standard output of all galaxy properties, simply combine the "``standard``" and "``analyzer``" :galacticus-class:`mergerTreeOutputterClass`\ es:

.. code-block:: none

    <mergerTreeOutputter value="multi">
     <mergerTreeOutputter value="standard"/>
     <mergerTreeOutputter value="analyzer"/>
    </mergerTreeOutputter>

The analysis to be performed is determined by the :galacticus-class:`outputAnalysisClass`. For example, to compute the GAMA stellar mass function :cite:p:`baldry_galaxy_2012` you would include in your parameter file:

.. code-block:: none

    <outputAnalysis value="massFunctionStellarBaldry2012GAMA"/>

The resulting mass function is written to the output file in a group named "``analyses/massFunctionStellarBaldry2012GAMA``", which will contain the following datasets:

``massStellar``
   The stellar mass corresponding to each bin in the mass function;

``massFunction``
   The mass function expectation from Galacticus;

``massFunctionCovariance``
   The covariance matrix for the mass function expectation from Galacticus;

``massFunctionTarget``
   The observed mass function from :cite:t:`baldry_galaxy_2012`;

``massFunctionCovarianceTarget``
   The covariance matrix for the observed mass function from :cite:t:`baldry_galaxy_2012`.

The output group also contains a "``logLikelihood``" attribute which gives the log-likelihood of the model given this data. Additional attributes provide a description of the analysis that allows a simple plot comparing the Galacticus and observed mass functions to be made using the `Dendros <https://github.com/galacticusorg/dendros>`_ package. Dendros auto-discovers on-the-fly analyses in a Galacticus output file and produces comparison plots. Dendros is available on `PyPI <https://pypi.org/project/dendros/>`_ and can be installed via:

.. code-block:: none

   pip install dendros

Multiple on-the-fly analyses can be performed, by simply grouping them inside a "``multi``" :galacticus-class:`outputAnalysisClass`, e.g.:

.. code-block:: none

    <outputAnalysis value="multi">
     <outputAnalysis value="massFunctionStellarBaldry2012GAMA"/>
     <outputAnalysis value="massFunctionHIALFALFAMartin2010"  />
    </outputAnalysis>

which would then compute both the stellar mass function of :cite:t:`baldry_galaxy_2012` and the HI mass function of :cite:t:`martin_arecibo_2010`.

.. _manual-sec-volumeLimitedSamples:

Building Volume Limited Samples
-------------------------------

The ``mergerTreeWeight`` property (see Section :galacticus-ref:`mergerTreeDatasets`) property specifies the weight to be assigned to each merger tree in a model to construct a representative (i.e. volume limited) sample of galaxies. Galacticus does not typically generate every merger tree in a fixed volume of the Universe (as an N-body simulation might for example) as it's generally a waste of time to simulate millions of low mass halos and only a small number of high mass halos. The ``mergerTreeWeight`` factors correct for this sampling. If merger trees are being built, then the ``mergerTreeWeight``, :math:`w_i`, for each tree of mass :math:`M_i` (where the trees are ranked in order of increasing mass) is given by

.. math::

   w_i = \int_{M_\mathrm{min}}^{M_\mathrm{max}} n(M) \mathrm{d}M,

where :math:`n(M)` is the dark matter halo mass function and

.. math::

   M_\mathrm{min} &=& \sqrt{M_{i-1}M_i}, \\
   M_\mathrm{min} &=& \sqrt{M_i M_{i+1}}.

Suppose, for example, that we wish to construct a luminosity function of galaxies. In particular, we consider a luminosity bin :math:`k` which extends from :math:`L_k-\Delta k/2` to :math:`L_k+\Delta k/2`. If tree :math:`i` contains :math:`N_i` galaxies with luminosities :math:`l_{i,j}`, where :math:`j` runs from :math:`1` to :math:`N_i`, then the luminosity function in this bin is given by:

.. math::

   \phi_k = \sum_i \sum_{j=1}^{N_i} \left\{ \begin{array}{ll} w_i & \hbox{ if  } L_k-\Delta k/2 < l_{i,j} \le L_k+\Delta k/2 \\ 0 & \hbox{ otherwise.} \end{array} \right.

.. _manual-sec-RunningGrids:

Running Grids of Models
-----------------------

You can easily write your own scripts to generate parameter files and run Galacticus on these files. An example of such a script is ``scripts/aux/launch.py``. This script will loop over a sequence of parameter values, generate appropriate parameter files, run Galacticus using those parameters and analyze the results. This script currently supports running of Galacticus on a local machine, via a PBS queue (as multiple jobs or a single job), or on a `Condor <https://research.cs.wisc.edu/htcondor/>`_ cluster. To run the script simply enter:

.. code-block:: none

    ./scripts/aux/launch.py <runFile>

This will launch a single instance of the script. Multiple instances can be launched and will share the work load (i.e. they will not attempt to run a model which another instance is already running or has finished). If multiple instances are to be launched on multiple machines a command line option to ``launch.py`` can be used to ensure that they do not duplicate work. Adding ``-{}-instance 2:4`` for example will tell the script to run only the second model from each block of four models it finds. Launching for ``launch.py`` scripts on four different machines with ``-{}-instance 1:4``, ``-{}-instance 2:4``, ``-{}-instance 3:4`` and ``-{}-instance 4:4`` will then divide the models between those machines.

The ``runFile`` is an XML file with the following structure:

.. code-block:: none

   <parameterGrid>
    <modelRootDirectory>models.new</modelRootDirectory>
    <baseParameters>newBestParametersQuick.xml</baseParameters>
    <compressModels>no</compressModels>
    <splitModels>4</splitModels>

    <launchMethod>pbs</launchMethod>

     <local>
      <threadCount>3</threadCount>
      <ompThreads>4</ompThreads>
     </local>

    <condor>
     <galacticusDirectory>/home/condor/Galacticus/v0.9.3</galacticusDirectory>
     <universe>vanilla</universe>
     <environment>LD_LIBRARY_PATH=/usr/lib:/usr/lib64:/usr/local/lib</environment>
     <requirement>Memory &gt;= 1000 &amp;&amp; Memory &lt; 2000</requirement>
     <transferFile>{PWD}/myFile.data</transferFile>
     <wholeMachine>true</wholeMachine>
     <postSubmitSleepDuration>5</postSubmitSleepDuration>
     <jobWaitSleepDuration>10</jobWaitSleepDuration>
    </condor>

    <pbs>
     <scratchPath>/scratch/me</scratchPath>
     <wallTime>48:00:00</wallTime>
     <memory>3gb</memory>
     <ompThreads>8</ompThreads>
     <queue>standard</queue>
     <maxJobsInQueue>10</maxJobsInQueue>
     <mpiLaunch>yes</mpiLaunch>
     <mpiRun>/opt/openmpi/bin/mpirun</mpiRun>
     <environment>LD_LIBRARY_PATH=/home/me/software/Galacticus/Tools/lib64:$LD_LIBRARY_PATH</environment>
     <postSubmitSleepDuration>10</postSubmitSleepDuration>
     <jobWaitSleepDuration>60</jobWaitSleepDuration>
    </pbs>

    <monolithicPBS>
     <mpiLaunch>yes</mpiLaunch>
     <nodes>1</nodes>
     <threadsPerNode>12</threadsPerNode>
     <ompThreads>6</ompThreads>
     <jobWaitSleepDuration>60</jobWaitSleepDuration>
     <analyze>no</analyze>
     <environment>LD_LIBRARY_PATH=/home/me/software/Galacticus/Tools/lib64:$LD_LIBRARY_PATH</environment>
     <includePath>/my/include/path</includePath>
     <libraryPath>/opt/sgi/mpt/mpt-2.04/lib</libraryPath>
     <shell>csh</shell>
     <pbsCommand>source /usr/share/modules/init/csh</pbsCommand>
     <pbsCommand>module load mpi-sgi/2.04_64</pbsCommand>
    </monolithicPBS>

    <parameters>
     <label>modelLabel</label>
     <stabilityThresholdStellar value="1.1"/>
     <stabilityThresholdStellar value="0.9"/>
    </parameters>

    <parameters>
     <starFormationFeedbackDisks value="powerLaw">
      <exponent value="2.5"/>
      <exponent value="3.0"/>
     </starFormationFeedbackDisks>
     <starFormationFeedbackDisks value="creasey2012"/>
    </parameters>

    <parameters>
     <imfSelection value="fixed">
       <imfSelectionFixed value="Chabrier" parameterLevel="top"/>
       <imfSelectionFixed value="Salpeter" parameterLevel="top"/>
     </imfSelection>
     <imfSelection value="diskSpheroid">
       <imfSelectionDisk     value="Chabrier"  parameterLevel="top"/>
       <imfSelectionSpheroid value="Kennicutt" parameterLevel="top"/>
     </imfSelection>
    </parameters>

    <parameters>
      <coolingFunction value="summation">
        <coolingFunction value="atomicCIECloudy"             iterable="no"/>
        <coolingFunction value="CMBCompton"                  iterable="no"/>
        <coolingFunction value="molecularHydrogenGalliPalla" iterable="no"/>
      </coolingFunction>
    </parameters>

   </parameterGrid>

Each ``parameters`` block contains a list of parameters following the format used in standard Galacticus parameter files, with the difference that each parameter can appear multiple times, each time with a different ``value`` attribute, as is the case for ``stabilityThresholdStellar`` in the first ``parameters`` element in the above. A model will be run for all possible combinations of these values. For nested parameters with multiple values, all possible values of these parameters will be looped over when, and only when, the appropriate value of the containing parameter is being used. For example, in the second ``parameters`` element in the above example, models will be run with subparameter ``[exponent]``\ :math:`=`\ ``2.5`` and ``3.5`` for the ``starFormationFeedbackDisks`` element only when ``[starFormationFeedbackDisks]``\ :math:`=`\ ``powerLaw`` and not when ``[starFormationFeedbackDisks]``\ :math:`=`\ ``creasey2012``. It is also possible to specify that subparameter should be promoted to the top-level of the parameter file. In the third ``parameters`` element in the above example, ``imfSelectionFixed`` will take on values of ``Chabrier`` and ``Salpeter`` only when ``imfSelection``\ :math:`=`\ ``fixed``, and the ``imfSelectionFixed`` element will be promoted from a sub-parameter of ``imfSelection`` to the top-level of the parameter file due to the presence of the ``parameterLevel="top"`` attribute. Finally in some cases a parameter which appears multiple times is not to be iterated over. In the fourth ``parameters`` element in the above example, this is the case for the ``coolingFunction`` subparameters. The addition of an ``iterable="no"`` attribute specifies that these parameters are not to be iterated over, but simply left as they are.

Some variables, which are expanded at run time, are available. These include:

``%%galacticusOutputPath%%``
   This will be expanded to the output path of a model. Useful for specifying paths for any additional output.

By default, each model is output into a sequentially numbered directory within the ``./models`` directory. By default, these directories have the prefix ``galacticus``. This can be changed by including a ``label`` element inside a ``parameters`` block, in which case the content of the ``label`` element will be used as the prefix. This root directory can be modified by the optional ``modelRootDirectory`` element. Additionally, a set of base parameters can be read from a file specified by the ``baseParameters`` file---these will be read before each model is run and before any variations in parameters for the specific model are applied. As such, it defines the default model around which parameter variations occur. Additional options that may be present in the file (as elements within the ``parameterGrid`` element) are:

``doAnalysis``
   If set to "no" then no analysis scripts will be run on completed models, otherwise, they will be. Optionally, the analysis script to run can be specified via the ``analysisScript`` element;

``emailReport``
   If set to "yes" a report will be e-mailed to the address specified in ``galacticusConfig.xml`` when a model fails. Otherwise, the report will be written to standard output instead.

``compressModels``
   If "no" then models are not compressed after being run. Otherwise, the contents of the model output directory will be compressed using ``bzip2``.

``splitModels``
   If set to an integer larger than :math:`1`, each Galacticus model will be split into that number of jobs, and those jobs will be launched (using the selected method) independently. Once finished, the outputs from these split models will be merged back into a single model. This allows, for example, effectively distributing a single Galacticus model over multiple nodes of a PBS cluster.

The method by which to launch jobs must be specified in the ``launchMethod`` element. Currently available options are:

``local``
   The models will be run on the local machine. Two additional options can be specified within a ``local`` XML block:

   ``threadCount``
      The number of individual model threads to be launched.

   ``ompThreads``
      The number of OpenMP threads to be used by each model.

``pbs``
   Jobs will be submitted to a ``PBS`` batch queue system. The following options are available and can be specified within a ``pbs`` XML block:

   ``scratchPath``
      An optional path to which the model output will be written at run time. At the completion of each run, the data will be transferred to the usual output location. This is useful to avoid network I/O during run time;

   ``wallTime``
      A limit on the wall time allowed for each model (optional);

   ``memory``
      A limit on the memory allowed for each model (optional);

   ``ompThreads``
      The number of OpenMP threads to use for each model (optional). This is used to request an appropriate number of processors per node;

   ``queue``
      The name of the queue to submit the jobs to (optional);

   ``maxJobsInQueue``
      The maximum number of jobs to place in the queue. Additional jobs will be held and submitted once the number of jobs in the queue drops below this value (optional);

   ``mpiLaunch``
      If set to "``yes``" then the ``mpirun`` command will be used to launch a single copy of Galacticus (which may then spawn multiple OpenMP threads). If instead set to "``no``" then Galacticus is launch without the use of the ``mpirun`` command. Some systems will limit a code launched with ``mpirun`` to using just a single CPU (even if multiple OpenMP threads are spawned). In such cases, setting this option to "``no``" should permit multiple CPUs to be utilized.

   ``mpiRun``
      The path to the ``mpirun`` executable (optional---if not present, ``mpirun`` must be in ``PATH``);

   ``environment``
      Any settings here are set in each  PBS job in order to set appropriate environment variables on the machine where a job is executed;

   ``analyze``
      If set to "``yes``" then analysis (if any) will be performed as part of the PBS job. Otherwise, analysis is performed by the submitting machine.

   ``postSubmitSleepDuration``
      The time (in seconds) to wait after submitting each job. This prevents flooding the PBS queue manager with a large number of jobs in rapid succession.

   ``jobWaitSleepDuration``
      The time (in seconds) to sleep between successive checks of the PBS queue to see if any of the submitted jobs have finished.

``monolithicPBS``
   A single job will be submitted to a ``PBS`` batch queue system. This job will internally run multiple copies of Galacticus each with a different set of parameters. The following options are available and can be specified within a ``monolithicPBS`` XML block:

   ``nodes``
      The total number of nodes to use for the PBS job.

   ``threadsPerNode``
      The number of threads per node to use for the PBS job.

   ``ompThreads``
      The number of OpenMP threads to use for each model (optional). This is used to request an appropriate number of processors per node, and must be an factor of ``threadsPerNode``;

   ``scratchPath``
      An optional path to which the model output will be written at run time. At the completion of each run, the data will be transferred to the usual output location. This is useful to avoid network I/O during run time;

   ``wallTime``
      A limit on the wall time allowed for each model (optional);

   ``memory``
      A limit on the memory allowed for each model (optional);

   ``queue``
      The name of the queue to submit the jobs to (optional);

   ``mpiRun``
      The path to the ``mpirun`` executable (optional---if not present, ``mpirun`` must be in ``PATH``);

   ``environment``
      Any settings here are set in each  PBS job in order to set appropriate environment variables on the machine where a job is executed;

   ``analyze``
      If set to "``yes``" then analysis (if any) will be performed as part of the PBS job. Otherwise, analysis is performed by the submitting machine.

   ``postSubmitSleepDuration``
      The time (in seconds) to wait after submitting each job. This prevents flooding the PBS queue manager with a large number of jobs in rapid succession.

   ``jobWaitSleepDuration``
      The time (in seconds) to sleep between successive checks of the PBS queue to see if any of the submitted jobs have finished.

``condor``
   Jobs will be submitted to a ``Condor`` cluster. The following options are available and can be specified within a ``condor`` XML block:

   ``galacticusDirectory``
      When a Galacticus job is submitted to a  Condor cluster the Galacticus executable and the input parameter file are transferred to the machine where the job runs. Other files, such as data files, are not transferred. Therefore, they must be already present on any remote machine on which the job can run. This option specifies where a complete Galacticus installation can be found on the remote machine. If not present, it defaults to ``/home/condor/Galacticus/v0.9.0``;

   ``universe``
      Specifies to which  Condor universe jobs should be submitted. Allowed options are "vanilla" and "standard". If the standard universe is to be used then Galacticus must have been linked with ``condor_compile``---the ``Makefile`` allows this if the relevant lines are uncommented;

   ``environment``
      Any settings here are passed to  Condor's ``environment`` option in order to set appropriate environment variables on the machine where a job is executed;

   ``requirement``
      Any setting here is passed to  Condor's ``requirements`` option to specify requirements for each job. Multiple ``requirement`` entries will be combined (using logical and).

   ``transferFile``
      Any files listed here will be transferred the Condor worker (and so will be accessible from the path in which Galacticus is running). The macro ``{PWD}`` will be automatically expanded to the present working directory. Multiple ``transferFile`` entries can be given.

   ``wholeFile``
      Setting this option to ``true`` will add ``+RequiresWholeMachine = True`` to the Condor submit file. If Condor has been configured to allow jobs to take over a whole machine\ [#]_, this will cause jobs to do so. This is useful if you want to run OpenMP Galacticus on a Condor cluster.

   ``postSubmitSleepDuration``
      The time (in seconds) to wait after submitting each job. This prevents flooding the Condor queue manager with a large number of jobs in rapid succession.

   ``jobWaitSleepDuration``
      The time (in seconds) to sleep between successive checks of the Condor queue to see if any of the submitted jobs have finished.

In addition to the ``galacticus.hdf5`` output file, each model directory will contain a file ``parameters.xml`` which contains the parameters used to run the model and ``galacticus.log`` which contains any output from Galacticus during the run.

If present, the file ``galacticusConfig.xml``, described in Section :galacticus-ref:`ConfigFile`, is parsed for configuration options. If the ``contact`` element is present, the listed name and e-mail address will be used to determine who should receive error reports should a model crash. The error report will contain the host name of the computer running the model, the location of the model output and the log file (which may be incomplete if output is being buffered). Additionally, any core file produced will be stored in the model directory for later perusal, and the state files (see Section :galacticus-ref:`Restarting`) for the run can also be found in the model directory.

.. _manual-sec-ConfigFile:

Configuration File
------------------

The file ``galacticusConfig.xml``, if present in the current ``GALACTICUS_EXEC_PATH``, is used to configure Galacticus and provide useful information. If no such file is present in ``GALACTICUS_EXEC_PATH`` then the existence of ``$HOME/.galacticusConfig.xml`` is checked and that file used for configuration if it is present.

This configuration file should have the following structure:

.. code-block:: none

   <config>
     <contact>
       <name>My Name</name>
       <email>me@ivory.towers.edu</email>
     </contact>
     <email>
       <host>
         <name>myComputerHostName</name>
         <method>smtp</method>
         <host>smtp-server.ivory.towers.edu</host>
         <user>myUserName</user>
         <passwordFrom>kdewallet</passwordFrom>
       </host>
       <host>
         <name>default</name>
         <method>sendmail</method>
       </host>
     </email>
   </config>

The name and e-mail address in the ``contact`` section will be stored in any Galacticus models run---this helps track the provenance of the model. The ``email`` section determines how e-mail will be sent. Within this section, you can place one or more ``host`` elements, the ``name`` element of which specifies the host name of the computer to which these rules apply (the ``default`` host is used if no other match is found). For each host, the ``method`` element specifies how e-mail should be sent, either by ``sendmail`` or via ``smtp``. For SMTP transport (which currently supports SSL connections only), you must specify the ``host`` SMTP server, ``user`` name. The ``passwordFrom`` element specifies how the password for the SMTP log in should be obtained. If set to ``input`` then the user will be prompted for the password as needed.

Writing Data To a Temporary File
--------------------------------

When running Galacticus on a compute cluster it is often advantageous to have output written to a local scratch disk during run time and only moved to networked storage after the run is complete. (Otherwise, Galacticus will perform many small writes to networked storage which can result in extremely slow run times.) To do this, simply set the parameter ``[galacticusOutputScratchFileName]`` to the full path of a file to write to on local scratch space. During the run, data will be written to this file. After the run is finished, Galacticus will move this file to its permanent location as specified by the parameter ``[galacticusOutputFileName]``.

Error Handling
--------------

If something goes wrong and Galacticus fails, it will attempt to provide an error report, including a backtrace of where in the code the error happened, and will list any previously suppressed warning messages which may have occurred.

Backtraces
~~~~~~~~~~

You can also request a backtrace from the code at any time by sending signal ``SIGUSR1`` to a running Galacticus process. To do use, issue the command:

.. code-block:: none

   kill -SIGUSR1 <pid>

from the command line, where ``<pid>`` is the process ID of the running Galacticus process.

.. _manual-sec-Restarting:

Restarting A Crashed Run
------------------------

If Galacticus crashes, it can be useful to restart the calculation from just prior to the crash to speed the debugging process. Galacticus has functionality to store and retrieve the internal state of any modules and to recover this to permit such restarting. Currently, this is implemented with the ``build`` and ``read`` methods of merger tree construction, such that the internal state is stored prior to commencing the building or reading of each tree, thereby allowing a calculation to be restarted with the tree that crashed. More general store/retrieve behavior is planned for future releases.

To cause Galacticus to periodically store its internal state include the following input parameter:

.. code-block:: none

     <stateFileRoot value="galacticusState" />

This will cause the internal state to be stored to files ``galacticusState.state`` prior to commencing building each merger tree. Should a tree crash then replace this input parameter with:

.. code-block:: none

     <stateRetrieveFileRoot       value="galacticusState"/>
     <mergerTreeConstructor value="build"           >
      <treeBeginAt value="N"/>
     </mergerTreeConstructor>

where ``N`` is the number of the tree that crashed. This will cause calculations to begin with tree ``N`` and for the internal state to be recovered from the above mentioned files. The resulting tree and all galaxy formation calculations should therefore proceed just as in the original run (and so create the same crash condition).

OpenMP
~~~~~~

When running a model in parallel using OpenMP, a separate state file will be written for each thread, with the thread number appended to the end of each state file name. For debugging purposes, it is suggested that a crashed OpenMP run be restarted using just a single thread. To do this, change the appended thread number on the state files corresponding to the thread which crashed to 0 such that they will be used by the single thread when the run is restarted.

Processing Individual Merger Trees In Parallel
----------------------------------------------

By default, Galacticus utilizes the available parallel threads to process multiple merger trees simultaneously, with one tree processed by each thread. When the total number of trees to be processed is large, and there are not a small number of outlier trees with masses very much larger than the other trees, this approach generally results in good parallel efficiency.

However, in cases where a small number of trees are much more massive than any other (or are just slow to process for some other reason) it may be more efficient to have multiple parallel threads process each tree. To achieve this, set ``[treeEvolveSingleForest]``\ :math:`=`\ ``true``. In this case, trees are processed sequentially, with multiple threads assigned to each tree. To do this, a tree is broken up into a set of time slices, or "sections". The number of sections between each successive output (or between the earliest node in the tree and the first output) is specified by the ``[treeEvolveSingleForestSections]`` parameter. Individual branches of the tree within each section are assigned to parallel threads. *Note that this results in valid evolution only if the evolution of disjoint tree branches are independent of each other.*

.. [#] This approach allows a direct connection to be made between the structure of the input parameter XML file and the internal object hierarchy used by Galacticus, allowing very fine-grained control over the composition of Galacticus functionality. In particular it permits easy construction of objects which work by modifying results from other objects, such as the :galacticus-class:`darkMatterProfileConcentrationSchneider2015` model for dark matter halo concentrations.
.. [#] This functionality requires that ``libmatheval`` is installed.
.. [#] Note that this does not capture the full details of the correlations between parameters, since it uses just the covariance matrix. For a more accurate calculation the full Monte Carlo Markov Chains used in the WMAP-9 parameter fitting should be used instead.
.. [#] As described `here <https://www-auth.cs.wisc.edu/lists/htcondor-users/2009-January/msg00086.shtml>`_ for example.
