.. _manual-sec-coding:

Coding Galacticus
=================

Node Class Hierarchy Builder
----------------------------

The hierarchy of classes describing a merger tree node and its galaxy (i.e. ``treeNode``, ``nodeComponent``, ``nodeComponentBasic``, ``nodeComponentBasicStandard``, etc.) are built automatically to realize the set of component implementations and properties specified in source files.

.. _manual-sec-nodeBuilderVariableDefinitions:

Variable Definitions
~~~~~~~~~~~~~~~~~~~~

Throughout the node objects builder code, Fortran variables are defined using a common specification. This is a Python ``dict`` containing the following entries:

``intrinsic``
   a string specifying the intrinsic Fortran type (``integer``, ``logical``, ``real``, ``double precision``, ``complex``, ``double complex``, ``type``, ``class``, ``procedure``);

``type``
   *[optional, except for type and class intrinsics]* a string specifying the type of the variables;

``attributes``
   a list of strings specifying all attributes of the variables;

``variables``
   a list of strings giving the names of the variables.

For example, rank-1, long integer arguments which will not be modified in their function would be specified as:

.. code-block:: none

   {
       'intrinsic':  'integer',
       'type':       'kind_int8',
       'attributes': ['dimension(:)', 'intent(in)'],
       'variables':  ['argument1', 'argument2'],
   }

The ``build`` Data Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``build`` data structure (a Python ``dict``) is used to accumulate all objects, variables, interfaces, and functions needed to build the class hierarchy used to represent nodes and components in Galacticus. The build pipeline lives under ``python/Galacticus/Build/Components/`` and is driven by ``scripts/build/buildCode.py``; sub-modules register hooks that mutate ``build`` during a sequence of named phases (``preValidate``, ``default``, ``gather``, ``scatter``, ``postValidate``, ``content``, ``types``, ``interfaces``, ``functions``). At the end of the node objects build process, this dict is processed to generate the required Fortran code, which is appended to ``build['content']``. The following subsections describe how to add information to it.

Component Classes and Implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``dict`` of structures defining all component classes, keyed by class name, is provided as ``build['componentClasses']``. Each component class structure has the form:

.. code-block:: none

   {
       'name':    <name of the class>,
       'members': <dict of structures defining all member implementations of the class>,
   }

Each member implementation structure has the form:

.. code-block:: none

   {
       'name':       <name of the implementation>,
       'properties': <dict of structures defining all properties of this implementation>,
   }

.. _manual-sec-buildHierarchyTypes:

Types
^^^^^

Definitions of derived types are accumulated in the ``dict`` ``build['types']``, with the key of each entry being the derived type's name. Each derived type structure has the form:

.. code-block:: none

   {
       'name':           <name of the derived type>,
       'comment':        <a description of the type>,
       'isPublic':       <True if the class should have public visibility, False otherwise>,
       'dataContent':    <list of variable definitions containing all data for this class>,
       'boundFunctions': [
           {
               'type':        <'procedure'|'generic'>,
               'descriptor':  <a function descriptor dict> [optional],
               'name':        <name of the bound method>,
               'function':    <name of function to bind to> [optional; not required if descriptor is provided],
               'description': <LaTeX-syntax description of the method> [optional; not required if descriptor is provided],
               'returnType':  <LaTeX-syntax return type of the method> [optional; not required if descriptor is provided],
               'arguments':   <LaTeX-syntax definition of all arguments to the method> [optional; not required if descriptor is provided],
           },
       ],
   }

The preferred approach is to provide a function descriptor (see Section :galacticus-ref:`buildHierarchyFunctions`), and to omit the ``function``, ``description``, ``returnType``, and ``arguments`` entries (which will be inferred from the function descriptor). The ``function``, ``description``, ``returnType``, and ``arguments`` entries will eventually be deprecated in favor of the ``descriptor`` entry.

.. _manual-sec-buildHierarchyInterfaces:

Interfaces
^^^^^^^^^^

Definitions of interfaces are accumulated in the ``dict`` ``build['interfaces']``, with the key of each entry being the interface's name. Each interface structure has the form:

.. code-block:: none

   {
       'name':      <name of the interface>,
       'comment':   <a description of the interface>,
       'intrinsic': <the intrinsic type of the function (or 'void' for a subroutine)>,
       'data':      <list of variable definitions containing all arguments for this interface>,
   }

.. _manual-sec-buildHierarchyFunctions:

Functions
^^^^^^^^^

Definitions of functions are accumulated either in the list ``build['functions']`` or are included as the ``descriptor`` in a derived type structure (this is the preferred method for functions that *are* bound to a derived type; see Section :galacticus-ref:`buildHierarchyTypes`). Each function structure has the form:

.. code-block:: none

   {
       'type':        <the type of function>,
       'name':        <function name>,
       'description': <LaTeX-syntax description of the function>,
       'modules':     <list of names of modules required by the function> [optional],
       'variables':   <list of variable definitions required by the function> [optional],
       'content':     <the code of the function (excluding opener, closer, and variable definitions)>,
   }

Module-scope Variables
^^^^^^^^^^^^^^^^^^^^^^

Any module-scope variables can be appended to the list ``build['variables']`` (e.g.\ ``build['variables'].append({...})``), using the usual definition format described in Section :galacticus-ref:`nodeBuilderVariableDefinitions`.

.. _manual-sec-sourceTreePreprocessor:

Galacticus Preprocessor Directives
----------------------------------

Galacticus has its own preprocessor for Fortran source files. This preprocesses parses each source file into an internal tree representation, performs various manipulations on that tree, and then outputs the preprocessed file for compilation. The preprocessor is used to automate and standardize many common tasks, through the inclusion of directives into the source code. Directives are specified in comment lines beginning ``!#``, and are written in XML. The remainder of this section describes the various preprocessor functionalities, and gives examples of their usage.

Source Code Introspection
~~~~~~~~~~~~~~~~~~~~~~~~~

The source code introspection functionality allows automated generation of information about the source code. Specifically:

* The directive ``{introspection:location}`` in source code will be replaced with a character string giving a backtrace of the current location in the code, including any function, module, and file (including line number).

Function Attributes
~~~~~~~~~~~~~~~~~~~

Galacticus allows the specification of attributes for functions which alter the way the compiler treats them. Function attributes are specified as:

.. code-block:: none

   !$GLC function attributes {attributes} :: {functionNames}

where ``{attributes}`` is a space-separated list of attributes, and ``{functionNames}`` is a space-separated list of functions to apply these functions to. Function attribute directives may appear anywhere in a file containing the named function, but it is good practice to locate them immediately before the function.

Currently, the supported attributes are:

``unused``
   The function is marked as being *possibly* unused. Compiler warnings about unused functions will be suppressed for this function.

Source Digest
~~~~~~~~~~~~~

The ``sourceDigest`` directive will generate an :term:`MD5 hash` hash of the source code of the file in which the directive is placed, along with the source code of any files upon which it depends. This can be useful in generating unique labels (e.g. to use as suffixes in file names) which automatically update if the source code is modified. To generate a source digest simply use:

.. code-block:: none

     !# <sourceDigest name="mySourceDigest"/>

A source digest will be generated and stored as a ``character(len=22)`` variable called ``mySourceDigest``.

Object Builder
~~~~~~~~~~~~~~

When constructing instances of a class from a provided parameter set, a common pattern is to need to construct other objects based on those parameters, which will be used by the instance. For example, a transfer function class might require a cosmological parameters object for its operation. In such cases, we often want to use the default instance of the required class unless a different instance is explicitly specified. The ``objectBuilder`` directive automates this process.

As an example, the following constructor requires an instance of the ``cosmologyParameters`` class, which is passes to an internal constructor:

.. code-block:: none

     function myConstructorParameters(parameters)
       use Input_Parameters2
       implicit none
       type(myClass                  )                :: myConstructorParameters
       type(inputParameters          ), intent(in   ) :: parameters
       class(cosmologyParametersClass), pointer       :: cosmologyParameters_

       !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
       myConstructorParameters=myConstructorInternal(cosmologyParameters_)
       return
     end function myConstructorParameters

The ``objectBuilder`` directive will assign a member of the class specified by the ``class`` attribute to the variable specified by the ``name`` attribute. If the parameter set specified by the ``parameters`` attribute contains an explicit definition of the relevant class, that definition will be used to construct the instance. Otherwise, the parent parameter set will be checked for a definition of the relevant class, and so on. If no definition is found when the root parameter set is found the default instance will be used. Note that objects created by an ``objectBuilder`` directive are directly associated with the element in the input parameters XML document from which they were created. Therefore, if a later ``objectBuilder`` requires the object from that same XML element, it will reuse the one previously created.

The ``objectBuilder`` directive by default searches for a parameter named ``{class}`` where ``{class}`` is the value of the ``class`` attribute in the directive. An alternative name may be specified via the addition of a ``parameterName`` attribute.

In cases where an explicit ``parameterName`` attribute is given, it is possible to specify an explicit default for the object if no such named parameter is found. (Where no explicit ``parameterName`` attribute is given the global default of the class is used.) A default is specified by adding a ``default`` element to the directive, which should contain default specification for the object (and any subobjects) in the usual format. For example:

.. code-block:: none

    !# <objectBuilder class="massDistribution" parameterName="diskMassDistribution" name="diskMassDistribution" source="globalParameters">
    !#  <default>
    !#   <diskMassDistribution value="exponentialDisk">
    !#    <dimensionless value="true"/>
    !#   </diskMassDistribution>
    !#  </default>
    !# </objectBuilder>

When a class accepts a list of objects of the same type (e.g.\ to implement a summation over multiple instances), the ``copy`` attribute may be used to build multiple copies in a loop. The value of ``copy`` is a loop variable (optionally with explicit bounds as a Fortran ``do`` statement). For example:

.. code-block:: none

    !# <objectBuilder class="mergerTreeEvolveTimestep" name="mergerTreeEvolveTimestep_%mergerTreeEvolveTimestep_" source="parameters" copy="i"/>

or with explicit bounds:

.. code-block:: none

    !# <objectBuilder class="stellarPopulationSpectraPostprocessor" name="postprocessors(i)%stellarPopulationSpectraPostprocessor_" source="parameters" copy="i=1,countPostprocessors"/>

In the first case, the loop index ``i`` is used as the ``copyInstance`` argument when retrieving the object, allowing the parameter node to hold multiple instances. In the second case, the loop runs from 1 to ``countPostprocessors`` explicitly.

Object Destructor
~~~~~~~~~~~~~~~~~

This directive can (and should) be used to destroy objects built by the ``objectBuilder`` directive. The ``objectDestructor`` directive automates the process of deciding if these objects should be destroyed (or merely have pointers to them nullified). As an example, the following destructor destroys two associated objects:

.. code-block:: none

     subroutine simpleDestructor(self)
       implicit none
       type(powerSpectrumPrimordialTransferredSimple), intent(inout) :: self

       !# <objectDestructor name="self%transferFunction_"       />
       !# <objectDestructor name="self%powerSpectrumPrimordial_"/>
      return
     end subroutine simpleDestructor

Constructor Assignments
~~~~~~~~~~~~~~~~~~~~~~~

A common requirement in object constructors is to assign the values of arguments to the constructor to corresponding entries in the object. The ``constructorAssign`` directive performs this assignment for a list of comma-separated variables. For example:

.. code-block:: none

    function stellarMassConstructorInternal(massThreshold)
       implicit none
       type            (galacticFilterStellarMass)                :: stellarMassConstructorInternal
       double precision                           , intent(in   ) :: massThreshold
       !# <constructorAssign variables="massThreshold"/>
       return
     end function stellarMassConstructorInternal

will cause the value of the ``massThreshold`` argument to ``stellarMassConstructorInternal%massThreshold``. If an argument name is prefixed with ``*`` in the variables list, pointer assignment is used instead of standard assignment.

State Storing
~~~~~~~~~~~~~

Galacticus supports storing its internal state to file to allow `restarts <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Usage.pdf\#sec.Restarting>`_. Code to store and restore the internal state of objects of a given class can be generated automatically through use of the ``stateStorable`` directive. An example is:

.. code-block:: none

    !# <stateStorable class="table">
    !#  <table1DGeneric>
    !#   <restoreTo variables="reset" state=".true."/>
    !#   <exclude variables="staticData"/>
    !#  </table1DGeneric>
    !#  <table2DLinLinLin>
    !#   <restoreTo variables="resetX, resetY" state=".true."/>
    !#  </table2DLinLinLin>
    !# </stateStorable>

This specifies that the ``table`` class (and all child classes) can and should be stored to file as part of the representation of the internal state. Code to store and restore all data associated with any object of this class (as well as restoring polymorphic objects to the correct type) will be generated. If certain variables of the class or subclass should be restored to specific values this can be specified through a ``restoreTo`` element placed within an element with the name of the class or subclass (e.g. the ``table1DGeneric`` in the above example). The ``restoreTo`` element should specify a comma-separated list of one or more variables to set in its ``variables`` attribute, and the state to which they should be restored in its ``state`` attribute. Any variables which should be excluded from state store/restore (e.g. if their values are known to be determined statically at construction) can be specified via a ``exclude`` element---a list of variables to exclude should be given as a comma-separated list in its ``variables`` attribute.

.. _manual-sec-eventHooks:

Event Hooks
~~~~~~~~~~~

Galacticus provides an event hook infrastructure which allows ``functionClass`` objects to hook into events triggered elsewhere in the code. There are two components to the event hook system: the definition of an event hook point in the code, and the attachment of functions to that event hook from within a ``functionClass`` implementation.

Defining an Event Hook Point
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An event hook point is defined using the ``eventHook`` directive. This directive is placed at the location in the code where the event should be triggered, and specifies the interface (i.e.\ the arguments) that any hooked functions must implement. For example:

.. code-block:: none

    !![
    <eventHook name="branchJumpPostProcess">
      <import>
        <module name="Galacticus_Nodes" symbols="treeNode"/>
      </import>
      <interface>
        type(treeNode), intent(inout), pointer :: node
      </interface>
      <callWith>node</callWith>
    </eventHook>
    !!]

The ``name`` attribute gives the name of the event hook. The ``interface`` element defines the Fortran declarations of the arguments to be passed to any hooked functions. The ``callWith`` element gives the list of actual arguments to pass when the event is triggered. The optional ``import`` element can be used to specify any types or symbols that should be imported into the abstract interface of hooked functions.

When the preprocessor encounters this directive, it generates code to call all functions that have been attached to this event hook. The event hook infrastructure (defined in ``Events_Hooks``) handles the list of attached functions and ensures each is called in the correct order (respecting any dependencies specified at attachment time), and with the correct OpenMP thread binding.

A simpler variant, ``eventHookStatic``, is available for cases where the set of hooked functions is fixed at compile time (i.e.\ the functions are known statically and do not change at runtime). In this case, the directive takes only a ``name`` attribute, and the hooked functions are identified by a matching directive of the same name, which specifies the function to call via a ``function`` attribute.

.. _manual-sec-autoHook:

Attaching Functions to an Event Hook
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ``functionClass`` implementation can attach one or more of its functions to event hooks by implementing an ``autoHook`` method. This method is automatically called whenever a new instance is created via deep copy (see Section :galacticus-ref:`functionClassAll`), and also on newly-constructed instances. The method should call the ``attach`` method on the relevant ``eventHookUnspecified`` object (or a typed variant) from the ``Events_Hooks`` module.

For example, the following ``autoHook`` implementation attaches the ``betaProfileCalculationReset`` function to the ``calculationResetEvent`` event hook:

.. code-block:: none

     subroutine betaProfileAutoHook(self)
       use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
       implicit none
       class(coolingRadiusBetaProfile), intent(inout) :: self

       call calculationResetEvent%attach(self,betaProfileCalculationReset, &
            & openMPThreadBindingAllLevels,label='coolingRadiusBetaProfile')
       return
     end subroutine betaProfileAutoHook

The ``attach`` method takes the following arguments:

``object_``
   The object to be passed as the first argument to the hooked function (typically ``self``).

``function_``
   The function to attach to the event hook.

``openMPThreadBinding``
   *(optional)* Specifies how the hooked function is bound to OpenMP threads. Possible values are ``openMPThreadBindingNone`` (default; the function is always called regardless of which thread triggers the event), ``openMPThreadBindingAtLevel`` (the function is only called if the current OpenMP level and thread number match), and ``openMPThreadBindingAllLevels`` (the function is called if the thread number matches at all OpenMP levels at or above the level at which the function was hooked).

``label``
   *(optional)* A label used for diagnostics and dependency resolution.

``dependencies``
   *(optional)* An array of ``dependency`` objects (of type ``dependencyExact`` or ``dependencyRegEx``) specifying ordering constraints relative to other hooks with matching labels.

Each ``functionClass`` implementation that provides an ``autoHook`` method must also detach its functions from event hooks in its destructor, using the ``detach`` method:

.. code-block:: none

     subroutine betaProfileDestructor(self)
       use :: Events_Hooks, only : calculationResetEvent
       implicit none
       type(coolingRadiusBetaProfile), intent(inout) :: self

       if (calculationResetEvent%isAttached(self,betaProfileCalculationReset)) &
          & call calculationResetEvent%detach(self,betaProfileCalculationReset)
       return
     end subroutine betaProfileDestructor

The ``isAttached`` method should always be called before ``detach`` to guard against attempting to detach a hook that was never attached. The hooked function is identified by the same function pointer used in the ``attach`` call.

Conditional Call
~~~~~~~~~~~~~~~~

In some instances it is useful to be able to call a function with different combinations of optional arguments depending on certain conditions. (For example, if some combinations of optional arguments are mutually exclusive.) This can be achieved using the ``conditionalCall`` directive. An example is:

.. code-block:: none

     !# <conditionalCall>
     !#  <call>self=massDistributionBetaProfile(beta{conditions})</call>
     !#  <argument name="densityNormalization" value="densityNormalization" parameterPresent="parameters"/>
     !#  <argument name="mass"                 value="mass"                 parameterPresent="parameters"/>
     !#  <argument name="outerRadius"          value="outerRadius"          parameterPresent="parameters"/>
     !#  <argument name="coreRadius"           value="coreRadius"           parameterPresent="parameters"/>
     !#  <argument name="dimensionless"        value="dimensionless"        parameterPresent="parameters"/>
     !# </conditionalCall>
     !# <inputParametersValidate source="parameters"/>

The ``call`` element specifies the function call, and contains the special sequence ``{conditions}`` which will be replaced with the conditionally-present arguments. One or more ``argument`` elements should specify the various arguments which should be included in the call. For each such element the ``name`` attribute specifies the name of the dummy argument in the called function, the ``value`` attribute specifies the value (or variables) to pass this this dummy argument. A condition for inclusion of the argument must also be specified. In the above, the special ``parameterPresent`` condition is used. The argument will be included in the call if a parameter with a name matching the ``name`` attribute exists in the parameter set named in the ``parameterPresent`` attribute. Alternatively a ``condition`` attribute can be given. An argument is included in the call if the expression given in the ``condition`` attribute evaluates to true.

Code will be generated to call the function with all possible combinations of arguments.

Optional Arguments
~~~~~~~~~~~~~~~~~~

Fortran supports optional arguments to functions, but does not provide for a default value if those arguments are not present. The ``optionalArgument`` directive allows a default value to be specified. In the following example, a default value is defined for the ``units`` argument:

.. code-block:: none

     double precision function simpleHubbleConstant(self,units)
       implicit none
       class  (cosmologyParametersSimple), intent(inout)           :: self
       integer                           , intent(in   ), optional :: units
       !# <optionalArgument name="units" defaultsTo="hubbleUnitsStandard" />

       select case (units_)
       case (hubbleUnitsStandard)
          ! Return the value using the default units.
       case ....
          ! Return the value using some other units.
       end select
       return
     end function simpleHubbleConstant

The ``optionalArgument`` directive should appear after variable declarations and before any attempt to use the optional argument, and should have two attributes, ``name`` and ``defaultsTo`` which give the name of the argument variable and its default value respectively. The preprocessor will add a new variable with the same name plus an underscore suffix, and will ensure that it is initialized to the default value if the optional variable is not present, otherwise setting it to the value of the optional variable.

Note that for optional arguments that are ``intent(out)`` or ``intent(inout)`` the preprocessor currently *does not* ensure that the value of the new variable is copied back to the argument prior to exit from the function.

Enumerations
~~~~~~~~~~~~

The ``enumeration`` directive allows specification of an enumeration (a set of labels), and (optionally) functions to decode such a label from user input. The following example illustrates this usage:

.. code-block:: none

   module Cosmology_Parameters

     !# <enumeration>
     !#  <name>hubbleUnits</name>
     !#  <description>Specifies the units for the Hubble constant.</description>
     !#  <visibility>public</visibility>
     !#  <validator>yes</validator>
     !#  <encodeFunction>yes</encodeFunction>
     !#  <entry label="standard" />
     !#  <entry label="time"     />
     !#  <entry label="littleH"  />
     !# </enumeration>

   contains

   subroutine Test_Enumeration()
       use ISO_Varying_String
       implicit none
       class(cosmologyParametersClass), pointer :: cosmologyParameters_
       cosmologyParameters_ => cosmologyParameters()
       write (0,*) "Enumeration contains ",hubbleUnitsCount," entries from ",hubbleUnitsMin," to ",hubbleUnitsMax
       write (0,*) "Hubble constant in little-h units is: ",cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
       if (enumerationHubbleUnitsEncode('hubbleUnitsStandard') == hubbleUnitsStandard) then
         write (0,*) "Enumeration decoding succeeded"
       else
         write (0,*) "Enumeration decoding failed"
       end if
       if (enumerationHubbleUnitsEncode('standard',includesPrefix=.false.) == hubbleUnitsStandard) then
         write (0,*) "Enumeration decoding succeeded"
       else
         write (0,*) "Enumeration decoding failed"
         end if
         write (0,*) "Name from time value is '",char(enumerationHubbleUnitsDecode(hubbleUnitsTime,includePrefix=.true.)),"'"
       return
     end subroutine Test_Enumeration

   end module Cosmology_Parameters

Enumerations must be defined in the declaration section of a ``module``. The encoding and decoding functions will only be generated if the ``encode`` element is present and has content ``yes``. The enumeration variables are given ``public`` visibility by default---this can be overridden using the ``visibility`` element. If the ``validator`` element is present and set to ``yes`` then a function is created which will return true if the given value is a valid one for the enumeration. Additionally, if the ``validator`` element is present and set to ``yes`` variables are created giving a count of the number of entries in the enumeration, along with the minimum and maximum values in the enumeration. Note that the ``description`` element is used to generate an entry for the enumeration in the document, and so should be written in LaTeX\ syntax.

Input Parameters
~~~~~~~~~~~~~~~~

The ``inputParameter`` directive reads an input parameter and assigns the appropriate value to the given variable. ``inputParameter`` directives must occur within the main body of a function, subroutine, or program. A default value can be specified if desired. The following example illustrates this usage:

.. code-block:: none

     subroutine simpleParametersRead()
       implicit none
       double precision :: hubbleConstant
       !# <inputParameter>
       !#   <name>HubbleConstant</name>
       !#   <source>myParameters</source>
       !#   <variable>hubbleConstant</variable>
       !#   <defaultValue>69.7d0</defaultValue>
       !#   <defaultSource>(\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultSource>
       !#   <description>The present day value of the Hubble parameter in units of km/s/Mpc.</description>
       !#   <type>real</type>
       !#   <cardinality>0..1</cardinality>
       !# </inputParameter>
       if (hubbleConstant < 0.0d0) write (0,*) "The universe is collapsing!"
       return
     end subroutine simpleParametersRead

In this case, the value of the ``HubbleConstant`` parameter is assigned to the ``hubbleConstant`` variable, with a default of :math:`69.7` if no value was specified in the input parameter file. If the ``source`` element is present, the parameter will be read from the named ``inputParameters`` set, otherwise the parameter will be read from the top-level of the parameters file. The ``defaultSource``, ``<description>The``, ``type``, and ``cardinality`` elements are used only for adding an entry for the input parameter to the documentation, and so should be written in LaTeX\ syntax.

It is also possible to specify a set of parameter which iterate over names defined by other directives. The following example would read one parameter named "``fileNameForXXXXXXIMF``" where "``XXXXXX``" equals the ``name`` element each ``imfRegisterName`` directive:

.. code-block:: none

       !# <inputParameter>
       !#   <iterator>fileNameFor(#imfRegisterName->name)IMF</iterator>
       !#   <source>parameters</source>
       !#   <variable>fileNames(IMF_Index("$1"))</variable>
       !#   <defaultValue>Galacticus_Input_Path()//"data/SSP_Spectra_imf$1.hdf5"</defaultValue>
       !#   <description>The name of the file of stellar populations to use for the named \gls{imf}.</description>
       !#   <type>string</type>
       !#   <cardinality>0..1</cardinality>
     !# </inputParameter>

Input Parameter Lists
~~~~~~~~~~~~~~~~~~~~~

The ``inputParameterList`` directive will construct a ``varying_string`` array containing the names of all input parameters which are defined in the unit in which the directive appears. The name of the array is specified by the ``label`` attribute of the ``inputParameterList`` directive. If a ``source`` attribute is specified in the directive then only parameters being read from the named variable will be included in the list, otherwise any parameters read will be included. Such a list can be used to validate the names of parameters passed to a function for example.

Function Classes
~~~~~~~~~~~~~~~~

Most\ [#]_ of the internal functionality within Galacticus is provided by "function classes". These are classes (in the object oriented sense) which model some particular physical entity or concept (e.g. the underlying cosmological model) and provide one or more functions associated with that entity or concept. A default implementation of each function class can be selected at run-time, allowing for simple user-defined control of model behavior. A function class is specified by a ``functionClass`` directive, together with one or more implementations of the function class specified by their own directives.

An example of a function class directive, which defines a class for cosmological parameters (which in this case, for simplicity, consists of just the Hubble constant) is given below:

.. code-block:: none

     !# <functionClass>
     !#  <name>cosmologyParameters</name>
     !#  <descriptiveName>Cosmological Parameters</descriptiveName>
     !#  <description>Object providing various cosmological parameters.</description>
     !#  <default>simple</default>
     !#  <defaultThreadPrivate>no</defaultThreadPrivate>
     !#  <stateful>no</stateful>
     !#  <calculationReset>no</calculationReset>
     !#  <method name="HubbleConstant" >
     !#   <description>Return the Hubble constant at the present day. The optional \mono{units} argument specifies if the return value should be in units of km/s/Mpc (hubbleUnitsStandard), Gyr$^{-1}$ (hubbleUnitsTime), or 100 km/s/Mpc (hubbleUnitsLittleH).</description>
     !#   <type>double precision</type>
     !#   <pass>yes</pass>
     !#   <argument>integer, intent(in   ), optional :: units</argument>
     !#  </method>
     !# </functionClass>

The directive should contain the following elements:

``name``
   The name of this function class.

``descriptiveName``
   A descriptive name for the function class, suitable for inclusion in the documentation.

``description``
   A description of the purpose of this function class. This description will be included into the documentation so should be written in LaTeX\ syntax.

``default``
   The default implementation to use for this function class if no choice is made in the input parameter file.

``defaultThreadPrivate``
   *(optional)* If present and set to ``yes`` then the default implementation of this function class will be made OpenMP thread private. Otherwise, the default implementation is shared between threads. A thread private default implementation can be useful if the function class may need to generate look-up tables unique to each thread on the fly for example.

``calculationReset``
   *(optional)* If present and set to ``yes`` then the default implementation of the function class is assumed to possibly want to reset its calculations when the active :term:`node` changes (see Section :galacticus-ref:`autoHook`). In this case, an additional method is generated for the function class: ``calculationReset`` with interface:

   .. code-block:: none

          subroutine calculationReset(self,thisNode)
            class  (functionClassBaseName), intent(inout)          :: self
            type   (treeNode             ), intent(inout), pointer :: thisNode
          end subroutine calculationReset

   This method has a null implementation for the base class of the function class, but can be overridden to reset calculations of any given implementation.

``method``
   Each ``method`` element defines a method which the function class will support, the name of which is given by a ``name`` attribute. The method definition must contain the following elements:

   ``description``
      A description of this method (in LaTeX\ syntax) suitable for inclusion into the documentation.

   ``type``
      The type of the function (e.g. ``double precision``; use ``void`` for a subroutine).

   ``pass``
      If ``yes`` pass the object that the method was called on as the first argument.

   ``argument``
      *(optional)* Zero or more declarations (in standard Fortran syntax) for the method arguments (there is no need to specify the declaration for the object upon which the method was called in the case where this object is passed to the method function).

Implementations of this function class must be declared with a directive having the same name as the function class (i.e. ``cosmologyParameters`` in the example above). The directive must give the name of the implementation as an attribute, and a description of the implementation suitable for inclusion into the documentation. Each implementation should be placed in a separate file---the preprocessor will find these files and merge the implementations and function class definition into a single file for compilation. Each implementation should declare a class which extends the basic function class, constructor interfaces, and any module-scope data required by the class. The implementation should also define all necessary functions required by the class (separated from the declarations by a ``contains`` keyword). An example is given below:

.. code-block:: none

     !# <cosmologyParameters name="cosmologyParametersSimple">
     !#  <description>Provides the Hubble constant: $H_0$.</description>
     !# </cosmologyParameters>
     type, extends(cosmologyParametersClass) :: cosmologyParametersSimple
        private
        double precision :: HubbleConstantValue
      contains
        final     ::                    simpleDestructor
        procedure :: HubbleConstant  => simpleHubbleConstant
     end type cosmologyParametersSimple

     interface cosmologyParametersSimple
        module procedure simpleDefaultConstructor
        module procedure simpleConstructor
     end interface cosmologyParametersSimple

   contains

     function simpleDefaultConstructor()
       implicit none
       type(cosmologyParametersSimple) :: simpleDefaultConstructor

       ! Construct an instance of this class using a value of the Hubble constant read from the input parameter file.
       !# <inputParameter>
       !#   <name>H_0</name>
       !#   <variable>simpleDefaultConstructor%HubbleConstantValue</variable>
       !#   <defaultValue>69.7d0</defaultValue>
       !#   <defaultSource>(\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultSource>
       !#   <description>The present day value of the Hubble parameter in units of km/s/Mpc.</description>
       !#   <type>real</type>
       !#   <cardinality>0..1</cardinality>
       !# </inputParameter>
       return
     end function simpleDefaultConstructor

     function simpleConstructor(HubbleConstant)
       implicit none
       type            (cosmologyParametersSimple)                :: simpleConstructor
       double precision                           , intent(in   ) :: HubbleConstant

       ! Construct an instance of this class using a value of the Hubble constant provided directly.
       simpleConstructor%HubbleConstantValue=HubbleConstant
       return
     end function simpleConstructor

     elemental subroutine simpleDestructor(self)
       implicit none
       type(cosmologyParametersSimple), intent(inout) :: self

       ! Do any clean-up required by this class when an instance goes out-of-scope.
       return
     end subroutine simpleDestructor

     double precision function simpleHubbleConstant(self,units)
       implicit none
       class  (cosmologyParametersSimple), intent(inout)           :: self
       integer                           , intent(in   ), optional :: units

       ! Do whatever is necessary to return the Hubble constant in the appropriate units.
       return
     end function simpleHubbleConstant

In the above example, we define a "simple" implementation of the cosmologyParameters class. Key points are:

Name:
   The name should always be prefixed with the function class name. In this case, we have a ``simple`` implementation of the ``cosmologyParameters`` function class, and so our name is ``cosmologyParametersSimple``.

Extends:
   The base class for the function class is always the function class name suffixed with ``Class``, in this case ``cosmologyParametersClass``. Implementations must always be extensions of either this base class, or of another implementation.

Procedures:
   The implementation must define procedures for all methods of the function class, *except* for where a method specified a ``code`` element in the function class directive (in which case a procedure for the method may still be optionally defined). Specification of a ``final`` function is encouraged.

Constructors:
   The implementation must specify at least one constructor, which takes no arguments (usually known as the default constructor). This constructor must create an instance of the implementation, setting any parameters from the input parameter file as necessary. Additional constructors may be defined as required.

Procedures:
   All required procedures (including constructors and destructors) should be given after a line containing the ``contains`` keyword. Galacticus coding policy is that all procedures associated with an implementation should be prefixed with the implementation name, ``simple`` in this case.

Functionality to store and restore the state (see `here <https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Usage.pdf\#sec.Restarting>`_) of classes built via a ``functionClass`` directive are automatically built. If variables of a given implementation should be restored to a specific state, this can be specified by adding a ``restoreTo`` element to the directive declaring the implementation. The ``restoreTo`` element should specify a comma-separated list of one or more variables to set in its ``variables`` attribute, and the state to which they should be restored in its ``state`` attribute. Any variables which should be excluded from state store/restore (e.g. if their values are known to be determined statically at construction) can be specified via a ``exclude`` element---a list of variables to exclude should be given as a comma-separated list in its ``variables`` attribute.

.. _manual-sec-deepCopy:

Deep Copy Infrastructure
^^^^^^^^^^^^^^^^^^^^^^^^

Deep copy functionality is generated automatically for all classes built via the ``functionClass`` directive. The generated ``deepCopy`` method copies all data members of an object into a destination object of the same type. In most cases this happens automatically, with the preprocessor detecting and appropriately handling pointer members, allocatable arrays, and nested ``functionClass`` objects. After a deep copy is made, ``autoHook`` (see Section :galacticus-ref:`autoHook`) is called on the destination object to re-establish any event hook attachments.

In some cases, the default deep copy behavior must be customized. This is done by adding a ``deepCopy`` element to the directive declaring the ``functionClass`` implementation. For example:

.. code-block:: none

    !![
    <coolingRadius name="coolingRadiusBetaProfile">
     ...
     <deepCopy>
      <functionClass variables="radiation"/>
     </deepCopy>
    </coolingRadius>
    !!]

The ``deepCopy`` element may contain the following sub-elements:

``functionClass``
   A ``functionClass`` element with a ``variables`` attribute specifying a comma-separated list of member variable names identifies member variables that are ``functionClass`` objects but are *not* declared as ``class(…), pointer`` (e.g.\ they are declared as ``type(…), pointer`` or are not pointers). The deep copy infrastructure will perform a full deep copy of these objects (recursively calling ``deepCopy`` on them) and will call ``autoHook`` on the copies.

``ignore``
   An ``ignore`` element with a ``variables`` attribute specifying a comma-separated list of variable names instructs the preprocessor to skip those variables during deep copy. This is useful for members whose values are managed externally (e.g.\ elements of a linked list that are populated separately) and should not be replicated during a deep copy.

``increment``
   An ``increment`` element with a ``variables`` attribute specifying a comma-separated list of variable names causes those variables to be incremented in the destination object after the copy. The optional ``atomic`` attribute (set to ``yes``) may be used to make the increment an OpenMP atomic operation.

In addition to the deep copy functionality generated within ``functionClass`` objects, the preprocessor provides three standalone directives for use in other contexts:

``deepCopyReset``
   Resets the internal ``copiedSelf`` pointer on one or more ``functionClass`` objects before beginning a new sequence of deep copies. Usage:

   .. code-block:: none

       !![
       <deepCopyReset variables="self%myObject_"/>
       !!]

``deepCopy``
   Performs a deep copy of a source ``functionClass`` object into a destination object. Usage:

   .. code-block:: none

       !![
       <deepCopy source="self%myObject_" destination="newObject_"/>
       !!]

   The ``autoHook`` method is called on the destination object after the copy.

``deepCopyFinalize``
   Resets the ``copiedSelf`` pointer after a sequence of deep copies is complete. Usage:

   .. code-block:: none

       !![
       <deepCopyFinalize variables="self%myObject_"/>
       !!]

These three directives are typically used in sequence when manually deep-copying a ``functionClass`` member: first ``deepCopyReset`` to prepare the source, then ``deepCopy`` to perform the copy, then ``deepCopyFinalize`` to clean up.

Generic Programming
~~~~~~~~~~~~~~~~~~~

The preprocessor supports generic programming by allowing generic types to be defined, which are automatically expanded to a set of specific types. Significant flexibility is provided to allow control over how each specific type is handled. A generic type is specific via a ``generic`` directive such as:

.. code-block:: none

     !# <generic identifier="Type">
     !#  <instance label="Logical"        intrinsic="logical"                         outputConverter="regEx@\textbrokenbar@(.*)@\textbrokenbar@char($1)@\textbrokenbar@"/>
     !#  <instance label="Integer"        intrinsic="integer"                         outputConverter="regEx@\textbrokenbar@(.*)@\textbrokenbar@$1@\textbrokenbar@"      />
     !#  <instance label="Double"         intrinsic="double precision"                outputConverter="regEx@\textbrokenbar@(.*)@\textbrokenbar@$1@\textbrokenbar@"      />
     !#  <instance label="LogicalRank1"   intrinsic="logical          , dimension(:)" outputConverter="regEx@\textbrokenbar@(.*)@\textbrokenbar@char($1)@\textbrokenbar@"/>
     !#  <instance label="IntegerRank1"   intrinsic="integer          , dimension(:)" outputConverter="regEx@\textbrokenbar@(.*)@\textbrokenbar@$1@\textbrokenbar@"      />
     !#  <instance label="DoubleRank1"    intrinsic="double precision , dimension(:)" outputConverter="regEx@\textbrokenbar@(.*)@\textbrokenbar@$1@\textbrokenbar@"      />
     !# </generic>

The above defines a generic type, which will be identified using the label "``Type``". The directive contains several ``instance`` elements, each of which specifies a specific type which should be implemented for the generic type. Each instance can contain an arbitrary number of attributes which specify strings or regular expressions which will be used to construct the specific implementation.

A generic directive applies to the entire unit within which it is scoped. The preprocessor will examine every element within that unit. If a generic tag (see below) is found in the opening of any subunit, that entire subunit is copied once for each instance, and any generic tags replaced with the appropriate content from the ``instance`` element. Where a generic tag is found in a non-opening line (and that line is not contained within a subunit whose opener *does* contain a generic tag), the line itself is replicated in the same way.

An example of the usage of generic tags using the above generic directive is:

.. code-block:: none

     type :: exampleType
        private
        .
        .
        .
      contains
        final     ::        exampleTypeDestroy
        procedure ::        exampleTypeSet{Type@\textbrokenbar@label}
        generic   :: set => exampleTypeSet{Type@\textbrokenbar@label}
     end type exampleType

   contains

     subroutine exampleTypeSet{Type@\textbrokenbar@label}(self,setValue)
       implicit none
       class           (exampleType), intent(in   ) :: self
       {Type@\textbrokenbar@intrinsic}             , intent(in   ) :: setValue

       {Type@\textbrokenbar@match@\textbrokenbar@^Logical@\textbrokenbar@! Do something to set a logical value.@\textbrokenbar@! Do something different to set a numerical value.@\textbrokenbar@}
       write (0,*) "Value is: ",{Type@\textbrokenbar@outputConverter@\textbrokenbar@setValue}
     end subroutine exampleTypeSet{Type@\textbrokenbar@label}

In this example, the ``exampleType`` class is defined to have a generic ``set`` method. The presence of the ``{Type@¦@label}`` generic tag will cause those lines to be replicated with the tag replaced by the content of the ``label`` attribute of each instance of the generic type. In the contained subroutine, a generic tag appears in the opener. As such, the entire subroutine will be replicated once for each instance of the generic type, and the generic tags replaced as appropriate.

When a generic instance attribute begins with ``regEx``, matching generic tags are handled differently. In particular, a match-and-replace regular expression is applied to the third element of the generic tag (elements are separated by broken vertical bar characters). The match and replace components of the regular expression are defined in the instance attribute, once again separated by broken vertical bars.

Finally, the special generic tag ``match`` acts as a ternary operator. If the regular expression specified in the third element matches the ``label`` attribute of a specific instance, the generic tag is replaced with its fourth element, otherwise it is replaced with its fifth element.

Numerical Tools
---------------

Galacticus provides a variety of tools to solve basic numerical problems. These can be found in files ``source/numerical.*``. Galacticus makes use of the `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_ for many of these tools, but typically provides a higher-level wrapper around those functions, providing a cleaner interface and, in some cases, additional functionality.

Finding Roots of Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~

Tools for solving equations of the form :math:`f(x)=0` are provided by the ``rootFinder`` object (available via the ``Root_Finder`` module). Typical use of this object is as follows:

.. code-block:: none

   ! Import the module.
   use Root_Finder
   ...
   ! Create a rootFinder object - make it OpenMP threadprivate so it can be used
   ! simultaneously by all threads.
   type(rootFinder), save :: finder
   !$omp threadprivate(finder)
   ...
   ! Check if our root finder has been initialized.
   if (.not.finder%isInitialized()) then
     ! Specify the function that evaluates f(x).
     call finder%rootFunction   (myRootFunction                     )
     ! Specify the type of root-finding algorithm - this is optional (Brent's
     ! method will be used by default).
     call finder%type           (GSL_Root_fSolver_Brent             )
     ! Specify the tolerances to use in finding the root. Both arguments are
     !optional - values of 1.0d-10 will be used for both absolute and relative
     ! tolerance by default.
     call finder%tolerance      (toleranceAbsolute,toleranceRelative)
     ! Specify how the initially provided range can be expanded to bracket the
     ! root. This is optional - if not provided no range expansion will be attempted.
     call finder%rangeExpand                                               &
          &  (                                                             &
          &   rangeExpandDownward          =0.5d0                        , &
          &   rangeExpandUpward            =2.0d0                        , &
          &   rangeExpandType              =rangeExpandMultiplicative    , &
          &   rangeDownwardLimit           =1.0d-3                       , &
          &   rangeUpwardLimit             =1.0d+3                       , &
          &   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
          &   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative
          &  )
   end if
   x=finder%find(rootGuess=1.0d0)
   .
   .
   .
   double precision function myRootFunction(x)
     implicit none
     double precision, intent(in   ) :: x
     ...
     return
   end function myRootFunction

The above example begins by importing the ``Root_Finder`` module and then creating a ``rootFinder`` object called ``finder``. This is made OpenMP ``threadprivate`` so that it may be used simultaneously by all threads. The first step is to initialize ``finder``---the ``isInitialized`` method tells us if this has already happened. The most important step is to specify the function that will evaluate :math:`f(x)`. This is done via the ``rootFunction`` method---once done, the ``rootFinder`` object is marked as initialized (and the ``isInitialized`` method will return ``true``). All other initialization steps are optional. In this example, we use the ``type`` method to specify that the ``Brent`` algorithm should be used for root finding. Any valid `GSL-supported root finding algorithm <http://www.gnu.org/software/gsl/manual/html_node/Root-Bracketing-Algorithms.html>`_ can be used. We then use the ``tolerance`` method to specify both the absolute and relative tolerances in the :math:`x` variable that must be attained to declare the root to be found. Both arguments are optional---default values of :math:`10^{-10}` will be used if either tolerance is not specified.

The final step of initialization is to call the ``rangeExpand`` method. This specifies how the initial guessed value or range for :math:`x` should be expanded to bracket the root. If you plan to always specify an initial range, and know that it will always bracket the root, you do not need to specify how the range should be expanded. In this case we've specified that range expansion is multiplicative---that is, the lower and upper values of :math:`x` defining the range will be multiplied by fixed factors until the root is bracketed---via the ``rangeExpandType=rangeExpandMultiplicative`` option. Alternatively, additive expansion is possible using ``rangeExpandType=rangeExpandAdditive``. The factors by which to multiply the lower and upper bounds of the range (or the factor to add in the case of additive expansion) are specified by the ``rangeExpandDownward`` and ``rangeExpandUpward`` options. It is possible to specify absolute lower/upper limits to the range via the ``rangeDownwardLimit`` and ``rangeUpwardLimit`` options. The range will not be expanded beyond these limits---if the root cannot be bracketed without exceeding these limits an error condition will occur. Finally, it is possible to indicate the expected sign of :math:`f(x)` at the lower and/or upper limits via the ``rangeExpandDownwardSignExpect`` and ``rangeExpandUpwardSignExpect`` options. Valid settings are ``rangeExpandSignExpectNegative``, ``rangeExpandSignExpectPositive``, and ``rangeExpandSignExpectNone`` (the default---implying that there is no expectation for the sign). If the sign of :math:`f(x)` is specified, then range expansion will stop once the expected sign is found. This can often improve efficiency, by allowing the range expander to expand the range in only one direction, resulting in a narrower range in which to search for the root.

Finally, we use the ``find`` method to return the value of the root. The first argument to ``find`` is the name of the function that evaluates :math:`f(x)`. Additionally, we must supply either ``rootGuess`` (a scalar value guess to use as the initial value for both the lower and upper values of the range---note that range expansion must be allowed in this case), or ``rootRange`` (a two-element array to use as the initial lower and upper values of the range bracketing the root).

The function evaluating :math:`f(x)` must have a form compatible with that shown for ``myRootFunction`` in the above example.

.. _manual-sec-codeUniqueLabels:

Computation Dependencies and Data Files
---------------------------------------

In many situations, some module in Galacticus might want to perform a calculation and then store the results to a file so that they can be reused later. A good example is the :galacticus-class:`transferFunctionCAMB` transfer function model, which computes a transfer function using  CAMB and stores this function in a file so that it can be re-read next time, avoiding the need to recompute the transfer function. A problem arises in such cases as the calculation may depend on the values of parameters (in our example, the transfer function will depend on cosmological parameters for example). We would like to record which parameter values this calculation refers to, perhaps encoding these into the file name, so that we can reuse these data in a future run only if the parameter values are unchanged. Given the modular nature of Galacticus it is impossible to know in advance which parameters will be relevant (e.g. does the cosmological parameter implementation have a parameter that describes a time varying equation of state for dark energy?).

To address this problem, Galacticus provides a mechanism to generate a unique descriptor for a given object. This descriptor encodes the parameter used to construct the object, and recursively includes the parameters used to construct any other object which is composited. A long-form (human readable) descriptor is returned by the ``descriptor`` method associated with all ``functionClass`` objects. Additionally, the ``hashedDescriptor`` method will return an MD5 hash of the descriptor, which will be unique (up to collisions) and can be used to identify the object both internally and, for example, when used as a suffix to file names. If the optional ``includeSourceDigest`` argument is set to true in the ``hashedDescriptor`` method then the hashed descriptor will include a hash of the source code of the object (and all composited objects) such that the descriptor will change should the source code be changed.

.. _manual-sec-Optimization:

Optimization
------------

In designing Galacticus, we opted for simplicity and clarity over speed. However, there are numerous parts of the code where optimization has been performed without a significant loss of clarity. In this section we discuss some of the techniques used.

Unique IDs and Stored Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frequently, a given property of a node may be required in many different aspects of the calculation. For example, the dark matter halo virial radius is used extensively in several distinct calculations within Galacticus. Frequently such calculations are performed for the same node, with the same properties several times\ [#]_. Obviously this is inefficient. It can be advantageous in such cases to store the result of a calculation and, if the function is called again with the same unchanged node to simply return the stored value. Galacticus facilitates this by two features.

The first feature is the "unique ID"---an integer number assigned to each node in Galacticus and which uniquely identifies a node (i.e. no two nodes processed in a Galacticus run will have the same unique ID). This number, which can be retrieved using the ``uniqueID`` property of a tree node, can be recorded each time a function is called. If called again for a node with the same unique ID as the previous call, the function can simply return the same answer as on the previous call.

The second feature accounts for the fact that the properties of a node will change, so even if a function is called on a node with the same unique ID it may occasionally need to recompute its result. Galacticus provides a calculation reset task (see Section :galacticus-ref:`autoHook`). All such tasks are performed just prior to the computation of derivatives for a node being evolved. A function can register a calculation reset task and use it to flag that it must update its calculations even if called again with the same node.

Global Functions
----------------

In very exceptional circumstances it is necessary to subvert the module hierarchy used by Galacticus to permit one module to call a function in a higher level module\ [#]_. Examples of where this approach is necessary usually involve initial bootstrapping (i.e. to establish halo density contrasts, which requires knowledge of the halo density profile, which in turn requires knowledge of the halo density contrast…).

Galacticus provides for global functions which facilitate this---specifically, it is possible to generate function pointers to higher level functions which are accessible via a very low-level module. To create a globally callable copy of a function use add the follow prior to the function definition:

.. code-block:: none

     !# <functionGlobal>
     !#  <unitName>myFunction</unitName>
     !#  <type>double precision</type>
     !#  <arguments>double precision , intent(in   ) :: mass, time</arguments>
     !# </functionGlobal>

Here, ``myFunction`` is the name of the function to make global, while the ``type`` and ``arguments`` (of which there may be more than one) elements are used to generate a suitable interface for the function. At run-time, a pointer to this function is then available from the ``Functions_Global`` module, named ``myFunction_``. Note that ``Functions_Global_Set`` provided by the ``Functions_Global_Utilities`` module must be called once to initialize these global function pointers prior to their use.

Galacticus Metadata
-------------------

Galacticus can collect metadata on its own activity. This is useful for profiling the code for example.

OpenMP Critical Section Wait Times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Galacticus makes use of :term:`OpenMP` for parallel operation. :term:`OpenMP` ``critical`` sections are used throughout to limit access to parts of the code which must be executed in serial. Threads will block at each ``critical`` section if another thread is currently within it. Galacticus can profile the amount of time spent waiting at each named :term:`OpenMP` ``critical`` section across all threads. To enable this profiling Galacticus should be compiled with ``-DOMPPROFILE`` added to the compile options (in ``GALACTICUS_FCFLAGS``). This will cause profiling instructions to be added to the source code prior to compilation. When this instrumented executable is run an ``metaData/openMP/`` group will be added to the output file. This group contains two datasets, ``criticalSectionNames`` and ``criticalSectionWaitTimes``. The first lists the names of all named :term:`OpenMP` ``critical`` sections, while the second lists the total number of seconds (across all threads) spent waiting at each section. A script is provided to analyze this metadata:

.. code-block:: none

   ./scripts/aux/openMPCriticalWaitProfile.py <modelFileName>

This script will analyze the :term:`OpenMP` ``critical`` section wait time metadata in the named file, reporting the total time spent waiting at critical sections followed by a rank ordered list of the top ten sections by wait time. This can be useful for assessing whether optimization might help to reduce :term:`OpenMP` ``critical`` section wait times.

Enumerations
------------

Enumerations are used to communicate options to many functions in Galacticus. All available enumerations, along with their members, are described below.

Defined Constants
-----------------

Galacticus defines numerous constants, including mathematical constants (e.g. :math:`\pi`), physical constants (e.g. the speed of light), unit conversions (e.g. Angstroms to meters), and prefixes (e.g. "kilo", "mega", etc.). These should be used whenever a constant is needed in the code---it is bad practice to use the numerical value of a constant directly in the code\ [#]_.

All defined constants are described below, along with references to their source, and the name of the module in which the constant is defined. To import a constant into a function, you would add a ``use`` statement. For example, for the constant :math:`\pi`, you would use:

.. code-block:: none

      use :: Numerical_Constants_Math, only : Pi

after which you can use this constant, e.g.:

.. code-block:: none

      volumeSphere=4.0d0*Pi/3.0d0*radiusSphere**3

Indicating Units of Defined Constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The unit system in which a constant is defined should be indicated by a suffix starting with an underscore. The default, which requires no suffix, is the SI system. So, for example, the defined constant ``massSolar`` (which has no suffix) will be in SI units of kilograms. The defined constant ``lymanSeriesLimitWavelengthHydrogen_atomic`` is in the "atomic" unit system, so will be in units of Angstroms.

Currently defined unit systems are:

no suffix:
   With no suffix, the defined constant is in the `SI <https://en.wikipedia.org/wiki/International_System_of_Units>`_ unit system.

``_internal``:
   The "internal" suffix indicates a defined constant is in Galacticus's internal unit system of :math:`\mathrm{M}_\odot`, Mpc, Gyr, and km/s.

``_atomic``:
   The "atomic" suffix indicates atomic units (Å).

``_cgs``:
   The "cgs" suffix indicates a defined constant is in the `CGS <https://en.wikipedia.org/wiki/Centimetre%E2%80%93gram%E2%80%93second_system_of_units>`_ unit system.

.. _manual-sec-definedConstants:

Available Defined Constants
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Classes
-------

The return type of each method, and the interfaces (i.e. the types and names of its arguments) are specified for each method of each object. A "``void``" return type indicates a subroutine.

.. _manual-sec-functionClassAll:

Methods available to all ``functionClass``\ es
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Although not methods directly bound to the ``functionClass`` class, the following methods will be automatically created for all child classes of ``functionClass`` when they are built by Galacticus.

   ``autoHook`` Insert any event hooks required by this object.

   * Return type: ``void``
   * Interface: ``()``

   ``deepCopy``  Perform a deep copy of the object. Here, ``<functionClassBase >`` refers to the base class of the created ``functionClass`` child.

   * Return type: ``void``
   * Interface: ``(destination)``

     ``class(<functionClassBase>), intent(inout) :: destination``

   ``deepCopyReset`` Reset deep copy pointers in this object and any objects that it uses.

   * Return type: ``void``
   * Interface: ``()``

   ``descriptor`` Return an input parameter list descriptor which could be used to recreate this object.

   * Return type: ``void``
   * Interface: ``(descriptor, includeClass)``

     ``type(inputParameters), intent(inout) :: descriptor``

     ``logical, intent(in   ), optional :: includeClass``

   ``hashedDescriptor`` Return a hash of the descriptor for this object, optionally include the source code digest in the hash.

   * Return type: ``type(varying_string)``
   * Interface: ``(descriptor, includeClass)``

     ``logical, intent(in   ), optional :: includeSourceDigest``

   ``objectType`` Return the type of the object.

   * Return type: ``type(varying_string)``
   * Interface: ``()``

   ``stateRestore`` Restore the state of this object from file.

   * Return type: ``void``
   * Interface: ``(stateFile, gslFile, stateOperationID)``

     ``integer, intent(in   ) :: stateFile``

     ``type(c_ptr), intent(in   ) :: gslStateFile``

     ``integer(c_size_t), intent(in   ) :: stateOperationID``

   ``stateStore`` Store the state of this object to file.

   * Return type: ``void``
   * Interface: ``(stateFile, gslFile, stateOperationID)``

     ``integer, intent(in   ) :: stateFile``

     ``type(c_ptr), intent(in   ) :: gslStateFile``

     ``integer(c_size_t), intent(in   ) :: stateOperationID``

.. [#] At this time the Galacticus code base is being transitioned to use this approach.
.. [#] For example, Galacticus's ODE solver will fix the properties of a node and then request that derivatives of all properties be computed. Some functions will then be called multiple times for the same node with unchanged properties.
.. [#] This usually arises because circular dependencies would arise if the called function were placed in a lower level module.
.. [#] Both because it is prone to mistakes (the more times a numerical value is used directly, the more chances there are for typos and other errors), and because using a named constant makes it much easier to understand *what* the code is doing and *why*.
