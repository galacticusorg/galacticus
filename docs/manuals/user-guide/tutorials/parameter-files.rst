Introduction to Galacticus Parameter Files
===========================================

File format
-----------

Galacticus parameter files are formatted as `XML <https://en.wikipedia.org/wiki/XML>`_ - a widely used and supported format. Most text editors will recognize XML (providing formatting and highlighting support), and most languages provide tools to work with and manipulate XML files (e.g. the `elementTree <https://docs.python.org/3/library/xml.etree.elementtree.html>`_ module in Python).

An extremely simple parameter file (taken from the `power spectrum tutorial <https://galacticus.readthedocs.io/en/latest/manuals/user-guide/tutorials/power-spectra.html>`_) might look like:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8"?>
   <!-- Parameters for tutorial on computing the power spectrum - https://github.com/galacticusorg/galacticus/wiki/Tutorial%3A-Power-spectra -->
   <parameters>
     <formatVersion>2</formatVersion>
     <lastModified revision="975fffd5aadf697707003f16bf9b4caee8ebe97d" time="2025-10-22T21:34:03" strict="true"/>

     <!-- Specify tasks to perform -->
     <task value="powerSpectra"/>

     <!-- Cosmological parameters -->
     <cosmologyFunctions value="matterLambda"/>
     <cosmologyParameters value="simple">
       <HubbleConstant value="70.20000"/>
       <OmegaMatter value=" 0.27250"/>
       <OmegaDarkEnergy value=" 0.72750"/>
       <OmegaBaryon value=" 0.04550"/>
       <temperatureCMB value=" 2.72548"/>
     </cosmologyParameters>

     <!-- Power spectrum options -->
     <transferFunction value="eisensteinHu1999">
       <neutrinoNumberEffective value="3.046"/>
       <neutrinoMassSummed value="0.000"/>
     </transferFunction>
     <powerSpectrumPrimordial value="powerLaw">
       <index value="0.961"/>
       <wavenumberReference value="1.000"/>
       <running value="0.000"/>
     </powerSpectrumPrimordial>
     <powerSpectrumPrimordialTransferred value="simple"/>
     <cosmologicalMassVariance value="filteredPower">
       <sigma_8 value="0.807"/>
       <tolerance value="1.0e-3"/>
     </cosmologicalMassVariance>

     <!-- Structure formation options -->
     <linearGrowth value="collisionlessMatter"/>

     <!-- Output options -->
     <outputFileName value="powerSpectrum.hdf5"/>
     <outputTimes value="list">
       <redshifts value="0.0 1.0"/>
     </outputTimes>

   </parameters>

The first few lines:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8"?>
   <!-- Parameters for tutorial on computing the power spectrum - https://github.com/galacticusorg/galacticus/wiki/Tutorial%3A-Power-spectra -->
   <parameters>
     <formatVersion>2</formatVersion>
     <lastModified revision="975fffd5aadf697707003f16bf9b4caee8ebe97d" time="2025-10-22T21:34:03" strict="true"/>

contain some header information.

The first line just indicates the XML version and encoding (in general, this never needs to be changed).

The line starting with ``<!--`` and ending with ``-->`` is a comment line - you can add comments anywhere in your file to help make it clear what the parameter means, why they were chosen, etc.

This is followed by a line:

.. code-block:: xml

   <parameters>

which marks the start of the actual parameters. This will be matched by a corresponding

.. code-block:: xml

   </parameters>

element at the end of the file.

The remaining two lines in this opening section:

.. code-block:: xml

   <formatVersion>2</formatVersion>
   <lastModified revision="975fffd5aadf697707003f16bf9b4caee8ebe97d" time="2025-10-22T21:34:03" strict="true"/>

give some useful metadata. They are optional but are usually included. The first, ``formatVersion``, identifies the version of the parameter file format used - it should currently be set to ``2`` (an older format, version ``1``, is no longer supported).

The second, ``lastModified``, can be used to indicate the version of Galacticus with which this parameter file is intended to be used. The ``revision`` element records the ``git`` commit hash corresponding to the last version of Galacticus for which this parameter file was updated to work with. When the parameter file is run, Galacticus will check if any more recent changes to the code may require changes to the parameter file (or if the parameter file corresponds to a later version of Galacticus than the one currently being used) - and will issue a warning if so. The ``time`` element is informational only, and indicates the time at which the parameter file was updated. Finally, the optional ``strict`` element, if set to ``true`` (if not present, it defaults to ``false``), causes any warnings about parameters (version mismatches as described above, or unknown parameter names, etc.) to be treated as fatal errors. This can be useful if you want to be sure that your parameter file has been updated before running it with a new version of Galacticus.

Parameters
----------

The actual parameters in the file all follow the same structure, for example:

.. code-block:: xml

   <cosmologyFunctions value="matterLambda"/>
   <countTrees         value="100"         />

defines a parameter named "``cosmologyFunctions``" to have a value of "``matterLambda``", and a parameter named "``countTrees``" to have a value of "``100``".

Checking your parameter file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can check that your parameter file has correct syntax using the ``validateParameters.py`` script. For example, to run this script on a parameter file named ``myParameterFile.xml`` you would do:

.. code-block:: bash

   ./scripts/aux/validateParameters.py myParameterFile.xml

Any problems with the parameters will be reported. As an example, try running the script on the (intentionally) incorrect parameter file ``testSuite/parameters/validation/duplicate-value-invalid.xml``:

.. code-block:: bash

   ./scripts/aux/validateParameters.py testSuite/parameters/validation/duplicate-value-invalid.xml

The output should be:

.. code-block:: text

   Parameter 'componentSpin' has multiple values

informing you that a parameter appears multiple times in the file.

Warnings
~~~~~~~~

When you run Galacticus, it checks that your parameter file is correctly formatted, and that it recognizes the parameters that you've included in there - reporting warnings if it finds any problems. These warnings are output in a section beginning:

.. code-block:: text

   -> WARNING: problems found with input parameters:

(with ``WARNING:`` colored in magenta if supported).

It is recommended to take warnings seriously, and to figure out how to fix your parameter file to remove the warnings. The presence of warnings is often an indicator that your parameter file won't do what you expect it to do.

Multiple copies
^^^^^^^^^^^^^^^

Warnings of the type:

.. code-block:: text

   multiple copies of parameter [componentDisk] present - only the first will be utilized

indicate that the named parameter ("``componentDisk``" in this case) appears more than once in a section of the parameter file. Most parameters are allowed to appear only once (see `here <#multiple-copies-1>`_ for an explanation of exceptions to this rule). The warning indicates that copies of the parameter after the first are simply being ignored. Often this means that you added a parameter to your file, but forgot to remove the old copy. If the old copy appears first in the file, your new copy will be ignored and you will get unexpected results. The solution is to remove additional copies of the parameter, retaining only the one that you want.

Empty value
^^^^^^^^^^^

Warnings of the type:

.. code-block:: text

   empty value for parameter [cosmologyFunctions]

indicate that the ``value`` element of the named parameter ("``cosmologyFunctions``" in this case) is empty, which is invalid. Either enter a value for the parameter or, if you want the parameter to take its default value (see `here <#cardinality>`_) simply remove the parameter from the file entirely.

Ambiguous value
^^^^^^^^^^^^^^^

Warnings of the type:

.. code-block:: text

   ambiguous value for parameter [countTrees]

indicate that you have mixed the new and old-styles of parameters, for example:

.. code-block:: xml

   <countTrees value="10"> <!-- New-style is to define the value via the value attribute here. -->
     <value>20</value>    <!-- Old-style is to define the value in a separate value element.  -->
   </countTrees>

which leads to an ambiguous choice of value for this parameter. The old-style, while still supported, is deprecated and should not be used, so replace the above with:

.. code-block:: xml

   <countTrees value="10"/>

Unrecognized parameter
^^^^^^^^^^^^^^^^^^^^^^^

Warnings of the type:

.. code-block:: text

   unrecognized parameter [cosmologyFnuctions] (did you mean [cosmologyFunctions]?)

indicate that the named parameter ("``cosmologyFnuctions``" in this case) was not recognized as an allowed parameter name. Sometimes this can be due to a typo - Galacticus will try to guess what you actually meant and offer a suggestion if it can. In the above, it guessed that you actually meant "``cosmologyFunctions``" - correcting this will remove the warning.

Where to find information about available parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The best resource to get information about available parameters, their allowed values, etc. is the `Galacticus Physics <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ document. Specifically, the section on `Functions <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ section details the available parameters, their types, cardinality, and default values for all allowed classes.

Types
~~~~~

The value must have the correct "type" for each parameter. For example, you can not set a value of "``one hundred``" for a parameter that expects an integer value. Allowed parameter types are:

* Strings, e.g. "``matterLambda``", "``starFormation``";
* Integers, e.g. "``123``", "``-45``";
* Reals, e.g. "``5.432e-8``";
* Booleans, e.g. "``true``", "``false``".

In the `Functions <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ section of the documentation a parameter description defines the cardinality, for example, the `haloMassFunctionShethTormen <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ class details the following parameters:

   *Parameters*

   ``[a]`` (real; 0,1) {``0.707e0``} The parameter :math:`a` in the Sheth et al. [2001] halo mass function fit.

   ``[p]`` (real; 0,1) {``0.3e0``} The parameter :math:`p` in the Sheth et al. [2001] halo mass function fit.

   ``[normalization]`` (real; 0,1) {``0.3221836349e0``} The normalization parameter :math:`A` in the Sheth et al. [2001] halo mass function fit.

The word "real" (in the "``(real; 0,1)``") in each of the above parameters indicates that a real type value is expected.

Cardinality
~~~~~~~~~~~

The cardinality of a parameter (i.e. how many values you are allowed to put in the "``value``" attribute) is typically 1 - i.e. you must provide a single value for the parameter.

Many parameters also allow a cardinality of 0 - corresponding to the case where the parameter simply is not present in the parameter file. In such cases the parameter will be assigned a default value.

Some parameters allow multiple values. For example, in defining a set of output redshifts you might have a parameter:

.. code-block:: xml

   <redshifts value="0.0 1.0 2.0"/>

which has a cardinality of 3, and so defines three output redshifts.

In the `Functions <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ section of the documentation a parameter description defines the cardinality, for example, the `outputTimesList <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ class details the following parameters:

   *Parameters*

   ``[times]`` (real; 1..\*) A list of (space-separated) times at which Galacticus results should be output. Times need not be in any particular order.

   ``[redshifts]`` (real; 0..\*) {``[0.0e0]``} A list of (space-separated) redshifts at which Galacticus results should be output. Redshifts need not be in any particular order.

The cardinality information is given in the parentheses after the parameter name. For the ``times`` parameter, the cardinality is ``1..\*`` - indicating that at least one value is required (if the parameter is present at all), and that there is no upper limit to the number of values that can appear.

For the ``redshifts`` parameter the cardinality is ``0..*`` - indicating that zero, one, or any number of values can be given. If the ``redshifts`` parameter is not present (i.e. there are zero values provided for it) then the default value (given in curly braces in the above) of ``[0.0e0]`` will be used instead.

Correspondence with classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Internally, Galacticus (mostly!) utilizes an `object-oriented programming <https://en.wikipedia.org/wiki/Object-oriented_programming>`_ design. You don't need to understand object-oriented programming to use Galacticus, but it can be helpful to appreciate that the content of your parameter file translates directly into the internal objects that Galacticus uses.

Here's a simple example. One *class* of object used in Galacticus is the ``cosmologyFunctions`` class. Here, "*class*" just means a group of objects which all provide the same functions and interfaces, but which might implement different physics models. As the name implies, the ``cosmologyFunctions`` class provides functions related to cosmology (including, for example, the relation between time, :math:`t`, and expansion factor, :math:`a(t)`). The ``cosmologyFunctions`` class contains three member classes:

* ``cosmologyFunctionsMatterLambda``: Assumes a universe containing matter and a cosmological constant;
* ``cosmologyFunctionsMatterDarkEnergy``: Assumes a universe containing matter and dark energy with an arbitrary equation of state;
* ``cosmologyFunctionsStaticUniverse``: Assumes a static (non-expanding) universe.

You can find all available classes and members in the `Functions <https://galacticus.readthedocs.io/en/latest/physics/index.html>`_ section of the documentation.

You can control which member of a class is used by setting the appropriate parameter in your parameter file. For example:

.. code-block:: xml

   <cosmologyFunctions value="matterLambda"/>

would select the ``cosmologyFunctionsMatterLambda`` class. Galacticus will build a corresponding ``cosmologyFunctionsMatterLambda`` object internally and use it for cosmological calculations.

Nested parameters
~~~~~~~~~~~~~~~~~~

Many classes require other parameters to configure their behavior. These are specified by setting those parameters *inside* the parameter selecting the class. For example, the ``cosmologyFunctionsMatterDarkEnergy`` class assumes an equation of state :math:`w(a)=w_0+w_1 a (1-a)` the dark energy, where :math:`w_0` and :math:`w_1` are specified via the parameters ``darkEnergyEquationOfStateW0`` and ``darkEnergyEquationOfStateW1`` respectively. This would be specified as follows:

.. code-block:: xml

   <cosmologyFunctions value="matterDarkEnergy">
     <darkEnergyEquationOfStateW0 value="-1.0"/>
     <darkEnergyEquationOfStateW1 value="+0.0"/>
   </cosmologyFunctions>

Many classes also make use of other classes. For example, the ``cosmologyFunctionsMatterLambda`` class requires a helper object of the ``cosmologyParameters`` class to provide it with cosmological parameter values (:math:`\Omega_\mathrm{m}`, :math:`H_0`, etc.). We can create an object of the ``cosmologyParameters`` class and provide it to the ``cosmologyFunctions`` object like this:

.. code-block:: xml

   <cosmologyFunctions value="matterLambda">
     <cosmologyParameters value="simple">
       <HubbleConstant  value="67.36000"/>
       <OmegaMatter     value=" 0.31530"/>
       <OmegaDarkEnergy value=" 0.68470"/>
       <OmegaBaryon     value=" 0.04930"/>
       <temperatureCMB  value=" 2.72548"/>
     </cosmologyParameters>
   </cosmologyFunctions>

When some object is looking for a helper object in the parameter file, the helper is first looked for inside the object itself - in the above the helper object is ``cosmologyParameters``, and is located inside the ``cosmologyFunctions`` object that needs it. If the helper is not found, Galacticus proceeds to look for the helper object at the next higher level in the parameter file. For example in:

.. code-block:: xml

   <cosmologyParameters value="simple">
     <HubbleConstant  value="67.36000"/>
     <OmegaMatter     value=" 0.31530"/>
     <OmegaDarkEnergy value=" 0.68470"/>
     <OmegaBaryon     value=" 0.04930"/>
     <temperatureCMB  value=" 2.72548"/>
   </cosmologyParameters>

   <cosmologyFunctions value="matterLambda"/>

the ``cosmologyFunctions`` object, not finding a helper ``cosmologyParameters`` object inside itself, looks outside itself, find the ``cosmologyParameters`` object there, and uses it. This is helpful as, for example, many classes need the ``cosmologyParameters`` object as a helper - this approach allows us to define it only once and have it be used by many other classes.

If no helper object is found, a default choice is adopted.

Multiple copies
~~~~~~~~~~~~~~~

Most parameters allow only one copy of the parameter to appear in any given section of the parameter file (and `warnings <#warnings>`_ will be issued if multiple copies appear). Some parameters, however, allow multiple copies. Most often, these correspond to classes which sum over the results of other classes (e.g. a cooling function class that sums the cooling rates from multiple other cooling function classes), or applies a sequence of multiple operators.

A common example of this is the ``nodePropertyExtractor`` class, which is used to determine which properties to output to the output file. Frequently we want to output many different properties. The ``nodePropertyExtractorMulti`` class allows this, for example:

.. code-block:: xml

   <nodePropertyExtractor value="multi">
     <nodePropertyExtractor value="nodeIndices"          />
     <nodePropertyExtractor value="indicesTree"          />
     <nodePropertyExtractor value="virialProperties"     />
     <nodePropertyExtractor value="radiusVelocityMaximum"/>
     <nodePropertyExtractor value="velocityMaximum"      />
   </nodePropertyExtractor>

will extract the properties provided by the ``nodePropertyExtractorNodeIndices``, ``nodePropertyExtractorIndicesTree``, ``nodePropertyExtractorVirialProperties``, ``nodePropertyExtractorRadiusVelocityMaximum``, and ``nodePropertyExtractorVelocityMaximum`` classes, and concatenate them for output.

Parameters in the output file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All parameters used while running a Galacticus model are written to the output file. This includes parameters set to their default values (and, therefore, did not appear in the original parameter file). Values are stored as attributes in the ``Parameters`` group of the output ``hdf5`` file. Nested parameters are stored in a group within the ``Parameters`` group, named for their parent parameter.

So, for example, the following parameter file structure:

.. code-block:: xml

   <cosmologyParameters value="simple">
     <HubbleConstant  value="70.2    "/>
     <OmegaMatter     value=" 0.2725 "/>
     <OmegaDarkEnergy value=" 0.7275 "/>
     <OmegaBaryon     value=" 0.0455 "/>
     <temperatureCMB  value=" 2.72548"/>
   </cosmologyParameters>

   <cosmologyFunctions value="matterLambda"/>

would appear in the output file as:

.. code-block:: text

      ATTRIBUTE "cosmologyFunctions" {
         DATATYPE  H5T_STRING {
            STRSIZE 12;
            STRPAD H5T_STR_SPACEPAD;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         }
         DATASPACE  SCALAR
         DATA {
         (0): "matterLambda"
         }
      }
      ATTRIBUTE "cosmologyParameters" {
         DATATYPE  H5T_STRING {
            STRSIZE 6;
            STRPAD H5T_STR_SPACEPAD;
            CSET H5T_CSET_ASCII;
            CTYPE H5T_C_S1;
         }
         DATASPACE  SCALAR
         DATA {
         (0): "simple"
         }
      }
      GROUP "cosmologyFunctions" {
      }
      GROUP "cosmologyParameters" {
         ATTRIBUTE "HubbleConstant" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 70.2
            }
         }
         ATTRIBUTE "OmegaBaryon" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 0.0455
            }
         }
         ATTRIBUTE "OmegaDarkEnergy" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 0.7275
            }
         }
         ATTRIBUTE "OmegaMatter" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 0.2725
            }
         }
         ATTRIBUTE "temperatureCMB" {
            DATATYPE  H5T_IEEE_F64LE
            DATASPACE  SCALAR
            DATA {
            (0): 2.72548
            }
         }
      }

Advanced topics
---------------

Math evaluation
~~~~~~~~~~~~~~~

Parameter values of type real and cardinality 1 can be written as equations using basic mathematical operations and referencing other (`type <#types>`_ real, `cardinality <#cardinality>`_ 1) parameters as variables.

This functionality is provided via `libmatheval <https://github.com/galacticusorg/libmatheval>`_. This library automatically compiled into the pre-compiled versions of Galacticus. If you compile your own Galacticus then you must have ``libmatheval`` available, otherwise math expressions in the parameter file will cause a run-time error.

Math expressions in a parameter are indicated by starting the ``value`` element with an ``=``, for example:

.. code-block:: xml

   <massHalo value="=10.0^13.5"/>

which would cause the ``massHalo`` parameter to evaluate to :math:`10^{13.5} \approx 3.16 \times 10^{13}`.

Basic mathematical functions are supported - full details are given in the ``libmatheval`` `documentation <https://galacticusorg.github.io/libmatheval/doc/evaluator_005fcreate.html#Description>`_.

It is also possible to reference the values of other parameters to use as variable, by enclosing their names in square brackets. For example:

.. code-block:: xml

   <parameters>

     <cosmologyParameters value="simple">
       <HubbleConstant  value="67.36000"/>
       <OmegaMatter     value=" 0.31530"/>
       <OmegaDarkEnergy value="=1.0-[./OmegaMatter]"/>
       <OmegaBaryon     value=" 0.04930"/>
       <temperatureCMB  value=" 2.72548"/>
     </cosmologyParameters>

   </parameters>

Here, the value of ``OmegaDarkEnergy`` is set to equal ``1.0`` minus the value of ``OmegaMatter`` - thereby enforcing a flat universe, i.e. :math:`\Omega_\mathrm{m}+\Omega_\Lambda=1`.

Multiple levels of the parameter structure can be navigated by separating parameter names by ``/``.

Paths to parameters can be relative to the location of the current parameter by prefixing the parameter name with ``./`` as in the above, and you can navigate to higher levels in the parameter hierarchy using ``../`` (i.e. ``.`` and ``..`` work similarly to their counterparts in the Unix ``cd`` command for navigating directories). If neither ``./`` or ``../`` is found at the start of the path, the path is assumed to start from the top of the parameter file.

So, in the above example each of:

#. ``OmegaDarkEnergy value="=1.0-[./OmegaMatter]"/>``
#. ``OmegaDarkEnergy value="=1.0-[cosmologyParameters/OmegaMatter]"/>``
#. ``OmegaDarkEnergy value="=1.0-[./../cosmologyParameters/OmegaMatter]"/>``

would achieve the same result (although the last one is a contrived example!).

Referencing other parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it can be useful to "reuse" an entire section of parameters elsewhere in the parameter file. For example, if you need to create the exact same object more than once, but want to avoid having to update both copies if you decide to change some parameter value. This can be achieved by using parameter IDs and ID references. For example:

.. code-block:: xml

   <cosmologicalMassVariance value="peakBackgroundSplit">
     <cosmologicalMassVariance value="filteredPower" id="myCMV">
       <sigma_8 value="0.812"/>
     </cosmologicalMassVariance>
   </cosmologicalMassVariance>

   <mergerTreeBuildController value="constrained">

     <mergerTreeBranchingProbabilityUnconstrained value="PCH">
       <cosmologicalMassVariance idRef="myCMV"/>
     <mergerTreeBranchingProbabilityUnconstrained/>

     <mergerTreeBranchingProbabilityConstrained value="gnrlzdPrssSchchtr">
       <cosmologicalMassVariance idRef="myCMV"/>
     <mergerTreeBranchingProbabilityConstrained/>

   </mergerTreeBuildController>

In this example, we want to use a peak-background split model for :math:`\sigma(M)` in the main part of our parameter file, but need the regular :math:`\sigma(M)` for calculations related to merger tree construction. To avoid having to copy the definition of ``cosmologicalMassVariance`` and its parameter ``sigma_8`` three times, we add an attribute ``id="myCMV"`` to the first instance of the ``cosmologicalMassVariance`` parameter. Then, in subsequent cases we can reference it using this ID:

.. code-block:: xml

   <cosmologicalMassVariance idRef="myCMV"/>

This means that, should we decide to change the value of :math:`\sigma_8` we need change it in only one location.

Note that, internally in Galacticus, these three instances of the ``cosmologicalMassVariance`` object are *not* copies of each other - instead they are *the exact same object*, reused in these three separate instances.

Ignoring warnings
~~~~~~~~~~~~~~~~~~

As described in the `Warnings <#warnings>`_ section, Galacticus will emit warnings about any potential problems that it finds in a parameter file. In general, the recommendation is to fix these potential problems!

However, if you are 100% sure that the parameter file is correct and doing what you intended, you can silence the warning for any parameter by adding an attribute ``ignoreWarnings="true"``. For example:

.. code-block:: xml

   <myIntentionallyEmptyParameter value="" ignoreWarnings="true">
     <someParameter value="123"/>
   </myIntentionallyEmptyParameter>

This would normally emit a warning:

.. code-block:: text

   MM: -> WARNING: problems found with input parameters:
   MM:     empty value for parameter [myIntentionallyEmptyParameter]
   MM: <-

but with the inclusion of the ``ignoreWarnings="true"`` attribute this warning is silenced.

Changing Parameters
~~~~~~~~~~~~~~~~~~~~

Galacticus allows parameter files to be modified before they are run, using a set of changes specified in one or more "change files". This can be useful to allow, for example, one base parameter file to be defined, and then to make small changes to allow different models/calculations to be performed (e.g. changing from running a fixed mass merger tree to simulate the Milky Way, to a distribution of merger tree masses to simulate an entire population of galaxies).

Change files are simple XML files with a syntax described below. They can be included on the Galacticus command line after the primary parameter file, e.g.:

.. code-block:: bash

   ./Galacticus.exe parameters.xml changes1.xml changes2.xml

which will cause Galacticus to first read the ``parameters.xml`` file, then apply changes from ``changes1.xml`` to it, and then to apply changes from ``changes2.xml``.

An example change file can be found `here <https://raw.githubusercontent.com/galacticusorg/galacticus/master/testSuite/parameters/changes.xml>`_, and looks as follows:

.. code-block:: xml

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
     <change type="replaceWith" path="nodeOperator/nodeOperator[@value='stellarFeedbackSpheroids']/stellarFeedbackOutflows/stellarFeedbackOutflows" target="nodeOperator/nodeOperator[@value='stellarFeedbackDisks']/stellarFeedbackOutflows/stellarFeedbackOutflows"/>

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
     <change type="update" path="nodeOperator/nodeOperator[@value='barInstability']/galacticDynamicsBarInstability/stabilityThresholdGaseous" value="0.75" append="false">
     </change>

   </changes>

A parameter change file should contain one or more ``change`` elements. Each ``change`` element must have a ``path`` attribute which gives an `XPath <http://www.w3schools.com/Xml/xpath_syntax.asp>`_ expression that identifies a parameter in the parameter file that will be acted upon by the change. (Note that currently only a limited set of XPath functionality is supported - matching of element names and single attributes only.) An empty ``path`` indicates the top-level section of the parameter file.

Each ``change`` element must also have a ``type`` attribute which specifies the type of change to be made. Supported types are:

* ``append``: This will append all parameters enclosed within the ``change`` element to the children of the parameter specified by the ``path`` attribute.
* ``insertBefore``: This will insert all parameters enclosed within the ``change`` element before the parameter specified by the ``path`` attribute.
* ``insertAfter``: This will insert all parameters enclosed within the ``change`` element after the parameter specified by the ``path`` attribute.
* ``replace``: This will replace the parameter specified by the ``path`` attribute with all parameters enclosed within the ``change`` element.
* ``replaceWith``: This will replace the parameter specified by the ``path`` attribute with the parameter (and all children) from the parameter file at the path specified by the ``target`` attribute (note that the parameter identified by the ``target`` element is *not* removed)
* ``replaceOrAppend``: This will replace the parameter specified by the path attribute, if it exists, with all parameters enclosed within the change element. If the parameter does not exist, instead, all parameters enclosed within the change element will be appended to the children of the parent element.
* ``encapsulate``: This will take the parameter specified by the ``path`` attribute, replace it with all parameters enclosed within the change element, and then reinsert the parameter specified by the ``path`` attribute as a child of the first element enclosed within the change element. This allows an existing parameter to be encapsulated inside a new parameter.
* ``remove``: This will remove the parameter specified by the ``path`` attribute.
* ``update``: This will update the value of the parameter specified by the ``path`` attribute with that given by the ``value`` attribute of the ``change`` element - if the optional ``append`` attribute is set to ``true`` then the value is appended to the current value (instead of replacing it).

Outputting changed parameter files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to have Galacticus output a new parameter file after processing any ``XInclude`` elements (to include other parameter files) and applying any change files as described above. The new parameter file will therefore be a standalone file that can be run without reference to any included files or changes files, and will maintain the structure and comments from the original files.

To output a processed parameter file use the ``--output-processed-parameters`` option. For example:

.. code-block:: bash

   ./Galacticus.exe myParameters.xml myChanges.xml --output-processed-parameters processedParameters.xml

This will parse the ``myParameters.xml`` file, apply changes from the ``myChanges.xml`` file and then output the modified parameters to ``processedParameters.xml``.
