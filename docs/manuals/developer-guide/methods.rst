Adding New Classes
==================

.. _manual-sec-CodeDirectives:

Code Directives
---------------

Galacticus is designed to be flexible and extensible, allowing you to add new classes and functionality without having to hack the code extensively. To achieve this it makes much use of embedded code directives which, for example, explain to the build system how a particular subroutine or function connects into the Galacticus code. Such code directives are enclosed in blocks delimited by ``!![`` and ``!!]``, and take the form of short blocks of XML. For example, a typical code directive might look like:

.. code-block:: none

    !![
    <accretionDisks name="accretionDisksADAF">
      <description>An accretion disk class for ADAF (advection-dominated accretion flow) disks.</description>
    </accretionDisks>
    !!]

This directive appears just prior to the definition of a type which implements an advection-dominated accretion flow (ADAF) accretion disk. The ``accretionDisks`` tag explains to the Galacticus build system that this file contains an implementation of the black hole accretion disks class. The build system will then merge this implementation into the ``accretionDisks`` function class (see Section :galacticus-ref:`functionClass`), making it selectable at run time. A complete reference of the general-purpose preprocessor directives is given in Section :galacticus-ref:`sourceTreePreprocessor`; this chapter describes the directives used to define node components and to register task functions.

.. _manual-sec-ComponentMassTypes:

Identifying Components and Mass Types
-------------------------------------

Many functions can be applied to different components or groups of components and to different types of mass within a node. In general, these functions make use of a set of label defined in the :ref:`Galactic_Structure_Options <module-galactic_structure_options>` module. Components are identified by a ``componentType`` label which can take on the following values:

``componentTypeAll``
   All components are matched;

``componentTypeDisk``
   Only disk components are matched;

``componentTypeSpheroid``
   Only spheroid components are matched.

``componentTypeBlackHole``
   Only black hole components are matched.

``componentTypeHotHalo``
   Only hot halo components are matched.

``componentTypeDarkHalo``
   Only dark matter halo components are matched.

Types of mass are identified by a ``massType`` which can take one of the following values:

``massTypeAll``
   All mass is included;

``massTypeDark``
   Only dark matter is included;

``massTypeBaryonic``
   Only baryonic mass is included;

``massTypeGalactic``
   Only galactic mass is included.

``massTypeGaseous``
   Only gaseous mass is included.

``massTypeStellar``
   Only stellar mass is included.

``massTypeBlackHole``
   Only black hole mass is included.

Components
----------

This section describes the internal structure of node components, and how a component is implemented.

Component Structure
~~~~~~~~~~~~~~~~~~~

Each node in the merger tree consists of an arbitrary number of "components", each of which can actually be an array, allowing multiple components of a given class. Each component represents a specific class of object, which could be a dark matter halo, a galactic disk or a black hole etc. A component of each class may be of one or more different implementations of that component class. Component classes are extensions of the ``nodeComponent`` base class, while each implementation is an extension of its component class (or, sometimes, of another implementation of that same class). Each component implementation type consists of a set of data\ [#]_, representing the properties (mass, size etc.) of the component, along with the rates of change (and ODE solver tolerances) for any properties which are evolvable. Additionally, each component contains a large number of methods (functions) which can be used to access its properties, query its interfaces and which are used internally to perform ODE evolution, output etc. The ``nodeComponent`` base class and all classes derived from it are built automatically at compile time by the ``componentBuilder`` source-tree process hook (which drives the generator package ``python/Galacticus/Build/Components``); the generated code is grafted into ``objects/nodes/_class.F90`` during preprocessing, so take a look in ``work/build/objects/nodes/_class.p.F90`` if you want to see it.

Extending Components
~~~~~~~~~~~~~~~~~~~~

It is possible to create a component which extends an existing component (see the discussion of the ``extends`` element in Section :galacticus-ref:`ComponentDefinition`). This capability is intended to allow new properties to be added to a component without having to create a whole new copy of the component. It is *not* intended to allow changes in the way in which the component is evolved through the halo hierarchy. (With the exception that rules to describe how the newly added properties will evolve through the halo hierarchy can be added of course.)

A simple example of this extension capability can be found in the :ref:`scaleShape <manual-sec-darkmatterprofilescaleshape>` dark matter profile component, which extends the :ref:`scale <manual-sec-darkmatterprofilescale>` dark matter profile component. In this case, the ``scaleShape`` component adds a new property, ``shape``, and specifies how it is to be initialized, evolved, output, and change by node promotion events. It *does not* affect how the ``scale`` property, inherited from the ``scale`` dark matter profile component, is evolved.

.. _manual-sec-ComponentImplement:

Implementing a New Component
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementing a new component involves writing some modules and functions which contain a definition of the component and, if necessary, handle initialization, creation, evolution, and responses to any events. Frequently, the easiest way to make a new component is to copy a previously existing one and modify it as needed. Details of the various functions that component modules must perform are given below. By convention, a component implementation lives in its own directory, ``objects/nodes/components/<component>/<implementation>/`` (with ``<component>`` and ``<implementation>`` acting as placeholders for the component class and implementation names, e.g. ``disk/standard``), and is split into three or four files, although some components might not need all of these files. These files are named as follows:

``objects/nodes/components/<component>/<implementation>/_class.F90``
   The primary file which describes the component and its properties, and which contains functions that manipulate the component as it evolves through a merger tree (ODE rates, behavior during mergers, etc.);

``objects/nodes/components/<component>/<implementation>/bound_functions.Inc``
   Contains functions which will be bound to the component object (i.e. the ``nodeComponent<Class><Implementation>`` class), and so will be available as type bound procedures. Generally, these functions will include any which get or set values of properties in the component, those which return information about its internal state (such as a :galacticus-class:`massDistributionClass` object describing the structure of the component), and any other functions which we may want to be overridden by extensions to the component.

``objects/nodes/components/<component>/<implementation>/data.F90``
   Contains any data which may need to be shared between the above two files. This might contain parameters which control some property of the component that is the same for all instances (e.g. if spheroids are modeled as Sérsic profiles all with the same value of the Sérsic index, that value might be placed into this file).

``objects/nodes/components/<component>/<implementation>/structure.F90``
   Contains any functions which implement the structure (e.g. density, rotation curve) of the component and which cannot be placed in ``objects/nodes/components/<component>/<implementation>/bound_functions.Inc`` due to dependencies on modules which in turn depend on the ``Galacticus_Nodes`` module.

In general, the ``_class.F90`` file is the place for the component definition and functions which process the component during tree evolution (including output), while ``bound_functions.Inc`` is intended for functions which record or report the internal state of the component.

.. _manual-sec-ComponentDefinition:

Component Definition
^^^^^^^^^^^^^^^^^^^^

Component definition itself takes the form of an embedded XML document. The following example illustrates such a document:

.. code-block:: none

     !![
     <component>
      <class>disk</class>
      <name>exponential</name>
      <isDefault>yes</isDefault>
      <properties>
       <property>
         <name>isInitialized</name>
         <type>logical</type>
         <rank>0</rank>
         <attributes isSettable="true" isGettable="true" isEvolvable="false" />
       </property>
       <property>
         <name>massStellar</name>
         <type>real</type>
         <rank>0</rank>
         <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
         <output unitsInSI="massSolar" unitsDescription="Solar masses" unitsQuantity="solMass" comment="Mass of stars in the exponential disk."/>
       </property>
       <property>
         <name>abundancesStellar</name>
         <type>abundances</type>
         <rank>0</rank>
         <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
         <output unitsInSI="massSolar" unitsDescription="Solar masses" unitsQuantity="solMass" comment="Mass of metals in the stellar phase of the exponential disk."/>
       </property>
       <property>
         <name>massGas</name>
         <type>real</type>
         <rank>0</rank>
         <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" isNonNegative="true" />
         <output unitsInSI="massSolar" unitsDescription="Solar masses" unitsQuantity="solMass" comment="Mass of gas in the exponential disk."/>
       </property>
       <property>
         <name>coolingMass</name>
         <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
         <type>real</type>
         <rank>0</rank>
       </property>
       <property>
         <name>halfMassRadius</name>
         <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
         <type>real</type>
         <rank>0</rank>
         <getFunction>Node_Component_Disk_Exponential_Half_Mass_Radius</getFunction>
       </property>
       <property>
         <name>luminositiesStellar</name>
         <type>real</type>
         <rank>1</rank>
         <attributes isSettable="true" isGettable="true" isEvolvable="true" isNonNegative="true" />
         <classDefault modules="Stellar_Population_Properties_Luminosities" count="Stellar_Population_Luminosities_Count()">0.0d0</classDefault>
         <output labels="':'//Stellar_Population_Luminosities_Name({i})" count="Stellar_Population_Luminosities_Count()" condition="Stellar_Population_Luminosities_Output({i},time)" modules="Stellar_Population_Properties_Luminosities" unitsInSI="luminosityZeroPointAB" unitsDescription="AB-magnitude zero point" unitsQuantity="4.465920e17 W/Hz" comment="Luminosity of disk stars."/>
       </property>
      </properties>
      <bindings>
       <binding method="attachPipes" function="Node_Component_Disk_Exponential_Attach_Pipes" type="void" />
      </bindings>
      <functions>objects/nodes/components/disk/exponential/custom_methods.inc</functions>
     </component>
     !!]

The elements of this document have the following meaning:

``class``
   *[Required]* Specifies the component class of which this is an implementation.

``name``
   *[Required]* Specifies the name of this specific implementation.

``extends``
   *[Optional]* If present, this element must contain ``class`` and ``name`` elements which specify the type of component which should be extended. The component then automatically inherits all properties and type-bound functions of the extended type.

``isDefault``
   *[Required]* Specifies whether or not this should be the default implementation of this class. Note that only one implementation of each class can be declared to be the default. If no implementation of a given class is declared to be the default then the (automatically generated) ``null`` implementation will be made the default.

``properties``
   *[Optional]* Contains an array of ``property`` elements which specify the properties of this implementation. Each member ``property`` has the following structure:

   ``name``
      *[Required]* The name of the property.

   ``type``
      *[Required]* The type (one of ``real``, ``integer``, ``logical``, ``history``, ``abundances``, ``chemicals``, or ``keplerOrbit`` at present) of the property.

   ``rank``
      *[Required]* The rank of this property (currently ``0`` for a scalar or ``1`` for a 1-D array).

   ``attributes``
      *[Required]* Attributes of this property:

      ``isSettable``
         If ``true`` then the value of this property can be set directly.

      ``isGettable``
         If ``true`` then the value of this property can be got directly.

      ``isEvolvable``
         If ``true`` this property evolves as part of the Galacticus :term:`ODE` system.

      ``createIfNeeded``
         If ``true`` then any attempt to get, set, or adjust the rate of this property will cause the component to be created if it does not already exist. This is useful if the component should be created in response to mass transfer from some other component for example.

      ``isDeferred``
         Contains a "``:``" separated list which can contain ``get``, ``set``, and ``rate``. The methods present in this list will not have functions bound to them at compile time. Instead a function will be created which allows a function to be bound to these methods at run time. For example:

         .. code-block:: none

              call myComponent%massFunction    (My_Component_Mass_Get_Function)
              call myComponent%massSetFunction (My_Component_Mass_Set_Function)
              call myComponent%massRateFunction(My_Component_Mass_Rate_Function)

         Additionally, a method is created which returns true or false depending on whether the method has been attached to a function yet, e.g.

         .. code-block:: none

             myComponent%massIsAttached    ()
             myComponent%massSetIsAttached ()
             myComponent%massRateIsAttached()

   ``output``
      *[Optional]* If present, the property will be included in the Galacticus output file. The following attributes control the details of that output:

      ``unitsInSI``
         The units of the output quantity in the SI system.

      ``unitsDescription``
         A human-readable description of the units of the output quantity.

      ``unitsQuantity``
         An ``astropy.units``-parseable description of the units of the output quantity.

      ``comment``
         A comment to be included with the HDF5 dataset for this property.

      ``condition``
         A statement which must evaluate to ``true`` or ``false`` and which will be used to determine if the property will be output. The present output time for is available as ``time``. In the case of an array property the construct "``{i}``" can be used to pass the index of the element for which the condition should be evaluated.

      ``modules``
         A comma-separated list of any modules required to perform the output (e.g. modules which contain functions or values that are used).

      Additional attributes are required for array properties:

      ``labels``
         This can be an array, declared as "``[L_1,…,L_N]``", specifying the suffix to be added to the property name for each component of the array in the output, or a function which returns the suffix. In the case of a function the construct "``{i}``" can be used to pass the index of the element for which the suffix is required.

      ``count``
         A statement which evaluates the the number of elements to be output (i.e. the length of the array).

   ``isVirtual``
      *[Optional]* If present and set to "``true``", this property is a virtual property. A virtual property has no data associated with it and must supply its own functions for getting, setting and adjusting its rate of change (if allowed by the property's attributes). Virtual properties are used for quantities which are derived from actual properties of the component implementation (for example, a star formation rate could be a virtual property if it is derived from an actual gas mass property) or for adjusting the rates of several actual properties simultaneously.

   ``isNonNegative``
      *[Optional]* If present and set to "``true``", this property should always be non-negative (i.e. zero or positive). This is typically the case for quantities such as masses for example. This attribute is informative only---it may or may not be taking into account by the class responsible for evolving the component properties. For example, the :galacticus-class:`mergerTreeNodeEvolverStandard` class will evolve properties marked as non-negative in such a way as to ensure they remain non-negative, but only if its parameter ``[enforceNonNegativity]=true``.

   ``getFunction``
      *[Optional]* Specifies the function to be used for getting the value of the property, overriding the default get function. The function must be included in the :ref:`Galacticus_Nodes <module-galacticus_nodes>` module by use of the ``functions`` element described below. Note that this function, by virtue of its privileged access to the internal structure of node components, can access the value of the data associated with the property using:

      .. code-block:: none

         myComponent%<property>Data%value

   ``setFunction``
      *[Optional]* The same as ``getFunction`` but defines a function to set the value of the property.

   ``classDefault``
      *[Optional]* Specifies the default value for this property if the component class has not been created (i.e. has no specific implementation yet). The content of this element gives the default value (which can be a scalar, an array, a function, etc.). Additional, optional attributes control the use of this element:

      ``modules``
         Specifies a comma-separated list of modules which are required to set the default values (e.g. modules which contain the value or function to be used).

      ``count``
         For array properties whose size is not known at compile-time, it is possible to specify a function which will return the appropriate size of the array at run-time. The scalar default value given in the ``classDefault`` element will then be replicated the appropriate number of times.

``bindings``
   *[Optional]* Contains an array of ``binding`` elements which specify functions to bind to this implementation. Each member ``binding`` has the following structure:

   ``method``
      The name of the bound method, such that the function can be accessed using

      .. code-block:: none

          myComponent%<method>(...)

   ``function``
      The function to which the method should be bound. (This function must be included in the :ref:`Galacticus_Nodes <module-galacticus_nodes>` module by use of the ``functions`` element described below.

   ``type``
      The type of function.

``functions``
   *[Optional]* Contains the name of a file which will be included into the :ref:`Galacticus_Nodes <module-galacticus_nodes>` module. This file can contain functions which will be bound to this implementation. By virtue of being included in the :ref:`Galacticus_Nodes <module-galacticus_nodes>` module these functions have privileged access to the internal structure of all node component objects.

Component Initialization
^^^^^^^^^^^^^^^^^^^^^^^^

Initialization of a component module (if necessary, for example, to read parameters or allocate workspace) can occur at a number of different points in the execution of Galacticus. Providing initialization occurs in advance of any calculations then any point is acceptable. One possibility is simply to call an initialization function at the head of all functions defined in the component module. This initialization function should return immediately if it has already been called (to avoid duplicate initialization). Another option is to use a :galacticus-class:`mergerTreeOperator` to perform initialization just before merger trees are constructed (the initialization function must again return immediately if it has been previously called).

Alternatively, a component may register initialization (and uninitialization) functions using the following directives. Each is placed immediately before the subroutine that it registers, and names that subroutine in a ``function`` attribute (these directives, like all of the task directives described in this chapter, attach their function to a static event hook point---see Section :galacticus-ref:`eventHooks`). For example:

.. code-block:: none

    !![
    <nodeComponentThreadInitializationTask function="Node_Component_Black_Hole_Standard_Thread_Initialize"/>
    !!]

``nodeComponentInitializationTask``
   The function is called once when node components are first initialized. It receives a single argument, ``parameters``, of type ``type(inputParameters), intent(inout)``, from which any required parameters can be read.

``nodeComponentThreadInitializationTask``
   The function is called by each thread prior to merger tree evolution (with the same interface as above), and can therefore be used to perform any "per thread" initialization---for example, constructing threadprivate objects, or attaching functions to event hooks.

``nodeComponentThreadUninitializationTask``
   The function (which takes no arguments) is called by each thread once merger tree evolution is complete, and should undo any per-thread initialization (e.g. detaching functions from event hooks, and releasing objects).

Component Access, Creation and Destruction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a node is created, it initially contains no components. A component must therefore create itself on the fly as needed. Typically, a component is first created when an attempt is made to set a property value, or to adjust the rate of change of a property value or in response to some event (e.g. a satellite component may be created in response to a node merging with a larger node). Requests for property values frequently *do not* require that the component exist, as a zero value can often be returned instead\ [#]_.

To access a component from a node, use:

.. code-block:: none

    myComponent => thisNode%<class>([instance=<N>,autoCreate=<create>])

where ``class`` is the component class required, the optional ``instance`` argument requests a specific instance of the component (relevant if the node contains more than one of a particular component, e.g. if it contains two supermassive black holes for example; if no ``instance`` is specified the first instance will be returned), and the ``autoCreate`` option specifies whether or not the component should be automatically created (assuming it does not already exist). ``autoCreate``\ :math:`=`\ ``true`` should be used to create components initially.

A component of a node can be destroyed using:

.. code-block:: none

   call thisNode%<class>Destroy()

.. _manual-sec-ComponentMethods:

Component Implementations
^^^^^^^^^^^^^^^^^^^^^^^^^

Component implementations optionally provide functions to get and set their properties (and to set the rate of change of evolvable properties) so that other components and functions within Galacticus to can interact with them in a way that is independent of the specific component implementation chosen. To permit this, Galacticus creates functions for each property to access it in all permitted ways. For example, the ``exponential`` implementation of the ``disk`` component class has a "``massStellar``" property defined by:

.. code-block:: none

    <method>
      <name>massStellar</name>
      <type>real</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </method>

This causes Galacticus to define several functions bound to the ``nodeComponentDisk`` class:

``massStellarIsSettable``
   Returns ``true`` if this property is settable;

``massStellarIsGettable``
   Returns ``true`` if this property is gettable;

``massStellarSet``
   Sets the value of this property to the supplied argument;

``massStellarGet``
   Gets the value of this property;

``massStellarRate``
   Cumulates its argument to the rate of change of this property;

``massStellarScale``
   Sets the absolute scale for this property used in ODE error control;

along with several others used internally for output, serialization etc.

.. _manual-sec-ComponentEvolution:

Component Evolution
^^^^^^^^^^^^^^^^^^^

All component properties which have an ``isEvolvable`` attribute set to ``true`` are included in Galacticus's ODE solver as the node is evolved forward in time. As described in Section :galacticus-ref:`ComponentMethods`, Galacticus will create two functions that permit the rate of change of a property adjusted and for the absolute scale used in ODE error control to be set.

When evolving ODEs the ODE solver aims to keep the error on property :math:`i` below

.. math::

   D_i = \epsilon_\mathrm{abs} s_i + \epsilon_\mathrm{rel} |y_i|,

where :math:`epsilon_\mathrm{abs}=`\ ``[odeToleranceAbsolute]``, :math:`epsilon_\mathrm{rel}=`\ ``[odeToleranceRelative]``, :math:`y_i` is the value of property :math:`i` and :math:`s_i` is a scaling factor which controls the absolute tolerance for this property. By default, :math:`s_i=1`, but this can be changed for a component utilizing the ``scaleSetTask`` directive. This allows a function to be called in which the component sets suitable scale factors for each of its properties prior to any ODE evolution being carried out. This can be very useful, for example, in cases where two components are coupled. Consider a case where a disk is transferring material to a spheroid via a bar instability. If the disk is orders of magnitude more massive that the spheroid then the rate of mass transfer can be very high (i.e. :math:`\dot{y}/y` for the spheroid will be large). With just a relative tolerance (i.e. the :math:`\epsilon_\mathrm{rel} |y_i|` term) this would require very short timesteps for the spheroid. However, in such cases we don't care about such tiny tolerances for the spheroid (since it will grow to be substantially more massive). Therefore, it may be appropriate to set :math:`s_i` to be equal to the sum of the disk and spheroid properties for example. The scale set directive and associated subroutine should follow this template:

.. code-block:: none

     !![
     <scaleSetTask function="Node_Component_Disk_Exponential_Scale_Set"/>
     !!]
     subroutine Node_Component_Disk_Exponential_Scale_Set(node)
       implicit none
       type (treeNode         ), pointer, intent(inout) :: node
       class(nodeComponentDisk), pointer                :: disk

       ! Get the disk component.
       disk => node%disk()
       ! Check if an exponential disk component exists.
       select type (disk)
       class is (nodeComponentDiskExponential)
         ...
         call disk%massStellarScale(massScale)
         ...
       end select
       return
     end subroutine Node_Component_Disk_Exponential_Scale_Set

Sensible choices for the :math:`s_i` factors can significantly speed-up execution of Galacticus.

Other Component Task Directives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several further task directives allow a component to hook into specific points of the evolution cycle. Each takes the same form as ``scaleSetTask`` above---the directive is placed immediately before the subroutine to be registered, and names it in a ``function`` attribute:

``preEvolveTask``
   The function is called for each node before differential evolution of that node begins. The interface is ``(node)`` with ``node`` of type ``type(treeNode), intent(inout), pointer``. This is useful, for example, to ensure a component has been initialized before evolution starts.

``inactiveSetTask``
   The function is called (when an ODE solver making use of inactive-property optimizations is in use) before evolution of each node, and should mark as inactive any properties of the component which are known to not affect the evolution of other (active) properties, by calling their ``JacobianZero`` methods. The interface is the same as for ``preEvolveTask``.

``stateStoreTask`` and ``stateRetrieveTask``
   The functions are called when the internal state is written to (or restored from) file to allow :ref:`restarts <manual-sec-restarting>`, and should store (or restore) any state held by the component (typically making use of the ``stateStore``/``stateRestore`` directives---see Section :galacticus-ref:`stateStorable`). The interface is ``(stateFile,gslStateFile,stateOperationID)`` where ``stateFile`` is an ``integer``, ``gslStateFile`` is a ``type(c_ptr)``, and ``stateOperationID`` is an ``integer(c_size_t)``, all ``intent(in)``.

Evolution Interrupts
^^^^^^^^^^^^^^^^^^^^

It is often necessary to interrupt the smooth ODE evolution of a node in Galacticus. This can happen if, for example, a galaxy mergers with another galaxy (in which case the merger must be processed prior to further evolution) or if a component must be created before evolution can continue. The rate adjust and rate compute subroutines allow for interrupts to be flagged via their ``interrupt`` and ``interruptProcedure`` arguments. If an interrupt is required then ``interrupt`` should be set to true, while ``interruptProcedure`` should be set to point to a procedure which will handle the interrupt. Then, providing no other interrupt occurred earlier, the evolution will be stopped and the interrupt procedure called before evolution is continued.

An interrupt procedure should have the form:

.. code-block:: none

     subroutine My_Interrupt_Procedure(thisNode)
       implicit none
       type(treeNode), pointer, intent(inout) :: thisNode

       ! Do whatever needs to be done to handle the interrupt.

       return
     end subroutine My_Interrupt_Procedure

Existing Classes
----------------

Function Classes
~~~~~~~~~~~~~~~~

Functions implement basic calculations (e.g. computing the power spectrum). Additional implementations of a ``functionClass`` are added using the a directive with the name of that ``functionClass``. The implementation should be placed in a file containing the directive. For example, for the ``mergerTreeTimestep`` ``functionClass`` the file must containing a directive of the form

.. code-block:: none

   !![
   <mergerTreeTimestep name="mergerTreeTimestepMyImplementation">
     <description>A short description of the implementation.</description>
   </mergerTreeTimestep>
   !!]

where ``MyImplementation`` is an appropriate name for the implementation. This file should be treated as a regular Fortran submodule, but without the initial ``submodule`` and final ``end submodule`` lines. That is, it may contain ``use`` statements and variable declarations prior to the ``contains`` line, and should contain all functions required by the implementation after that line. It is standard, but not required, that function names should be prefixed with the name of the implementation (e.g. "``MyImplementation``" in the above example). The file *must* define a type that extends the named ``functionClass`` class (or extends another type which is itself an extension of the named ``functionClass`` class), containing any data needed by the implementation along with type-bound functions required by the implementation.

Events
~~~~~~

Events are triggered during merger tree evolution. Examples are when a node needs to be promoted to its parent node, or when a minor node merges with its parent. Code responds to such events by attaching functions to the corresponding event hooks (see Section :galacticus-ref:`eventHooks`)---for example, the ``nodePromotion`` event (triggered when a primary progenitor reaches the time of its parent halo), the ``haloFormation`` event (triggered when a halo is deemed to have formed, or reformed), the ``satelliteHostChange`` event (triggered when a satellite node moves to a new host), and the ``mergerTreeExtraOutput`` event (triggered for each node at each output time, allowing extra, non-standard output to be written). Attachment is done at run time via the event's ``attach`` method (typically from a ``functionClass`` object's ``autoHook`` method, or from a node component's thread initialization task).

.. note::
   In older versions of Galacticus these events were handled by dedicated directives (``nodePromotionTask``, ``haloFormationTask``, ``satelliteHostChangeTask``, ``mergerTreeExtraOutputTask``, ``hdfPreCloseTask``, and the ``mergerTreeOutputPropertyCount``/``mergerTreeOutputNames``/``mergerTreeOutputTask`` family). These directives no longer exist---use the event hooks described above instead (for output tasks, use a :galacticus-class:`nodePropertyExtractorClass` implementation, or attach to the ``mergerTreeExtraOutput`` or ``outputFileClose`` events).

Tasks
~~~~~

Tasks are any processing which must be performed at some specific point in the execution of Galacticus (e.g. before or after evolving a node, or when the output file is opened or closed). A function is registered as a task by a directive placed immediately before it, naming it in a ``function`` attribute. Each such directive corresponds to a static event hook point (an ``eventHookStatic`` directive of the same name---see Section :galacticus-ref:`eventHooks`) at which all registered functions are called. An optional ``after`` attribute (naming another registered function) may be used to constrain the order in which the registered functions are called.

The component-related task directives (``nodeComponentInitializationTask``, ``nodeComponentThreadInitializationTask``, ``nodeComponentThreadUninitializationTask``, ``scaleSetTask``, ``preEvolveTask``, ``inactiveSetTask``, ``postStepTask``, ``stateStoreTask``, and ``stateRetrieveTask``) are described in Section :galacticus-ref:`ComponentImplement`. In addition, the following global task directives are available:

``outputFileOpen``
   The function is called immediately after the Galacticus output HDF5 file is opened (useful for writing one-time datasets, e.g. version or build information). The function takes no arguments.

``outputFileClose``
   The function is called immediately prior to closing the Galacticus output HDF5 file (typically to write accumulated data to that file). The function takes no arguments.

``universePostEvolveTask``
   The function is called once evolution of all merger trees is complete. The function takes no arguments.

Merger Tree Initialization Tasks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Additional tasks to be performed during merger tree initialization can be added using the ``mergerTreeInitializeTask`` directive. For example, the ``standard`` basic component uses this directive as follows:

.. code-block:: none

     !![
     <mergerTreeInitializeTask function="Halo_Mass_Accretion_Rate"/>
     !!]

Here, ``Halo_Mass_Accretion_Rate`` is the name of a subroutine which will be called to perform whatever initialization is required. The subroutine must have the following form:

.. code-block:: none

      subroutine Merger_Tree_Initialize_Task(node)
       implicit none
       type(treeNode), pointer, intent(inout) :: node
       .
       .
       .
       return
     end subroutine Merger_Tree_Initialize_Task

where ``node`` is the node to be initialized. The subroutine will be called once for each node in the tree.

Post-step Tasks
^^^^^^^^^^^^^^^

Additional methods for post-step tasks (i.e. things that should be done after each ODE solver step when evolving a node differentially) can be added using the ``postStepTask`` directive. For example, the standard hot halo component adds a task as follows:

.. code-block:: none

     !![
     <postStepTask function="Node_Component_Hot_Halo_Standard_Post_Step"/>
     !!]

Here, ``Node_Component_Hot_Halo_Standard_Post_Step`` is the name of a subroutine which will be called to perform whatever tasks are required. The subroutine must have the following form:

.. code-block:: none

      subroutine Post_Step_Task(node,status)
       implicit none
       type   (treeNode), intent(inout), pointer :: node
       integer          , intent(inout)          :: status
       .
       .
       .
       return
     end subroutine Post_Step_Task

where ``node`` is the node for which tasks should be performed. If any change is made to the state of the node then ``status`` should be set equal to ``GSL_Failure``. Tasks typically involve cleaning up after differential evolution.

.. _manual-sec-radius-solver:

Radius Solver Tasks
^^^^^^^^^^^^^^^^^^^

Galactic radii solver functions (see :galacticus-class:`galacticStructureSolver`) need to be able to interact with the components of a tree node to

#. Determine whether the node is physically plausible (and hence whether radii should be solved for at all);
#. Determine which components want a radius to be solved for, and get and set the properties of those components.

The ``radiusSolverPlausibility`` and ``radiusSolverTask`` directives facilitate this. A component which has a radius to be solved for should include directives of the form:

.. code-block:: none

    !![
    <radiusSolverPlausibility function="Component_Radius_Solver_Plausibility"/>
    !!]

and

.. code-block:: none

    !![
    <radiusSolverTask function="Component_Radius_Solver"/>
    !!]

where ``Component_Radius_Solver_Plausibility`` is the name of a subroutine which will specify whether or not the component is physically plausible for radius solving (e.g. has non-negative mass) and should have the following form:

.. code-block:: none

    subroutine Component_Radius_Solver_Plausibility(node)
       implicit none
       type(treeNode), intent(inout) :: node
       .
       .
       .
       return
    end subroutine Component_Radius_Solver_Plausibility

which should set ``node%isPhysicallyPlausible`` (and/or ``node%isSolvable``) to false if the component is not physically plausible, but should otherwise leave them unchanged. Additionally, ``Component_Radius_Solver`` is the name of a subroutine which will supply the necessary information about the node, and which should have the following form:

.. code-block:: none

    subroutine Component_Radius_Solver(node,componentActive,component,specificAngularMomentumRequired,specificAngularMomentum, &
         &                             Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
       implicit none
       type            (treeNode                    ), intent(inout)          :: node
       logical                                       , intent(  out)          :: componentActive
       type            (enumerationComponentTypeType), intent(  out)          :: component
       logical                                       , intent(in   )          :: specificAngularMomentumRequired
       double precision                              , intent(  out)          :: specificAngularMomentum
       procedure       (...                         ), intent(  out), pointer :: Radius_Get,Velocity_Get
       procedure       (...                         ), intent(  out), pointer :: Radius_Set,Velocity_Set
       .
       .
       .
       return
    end subroutine Component_Radius_Solver

When called, the subroutine should set ``componentActive`` to indicate whether or not this node contains an active component of the type, and ``component`` to the corresponding component type (e.g. ``componentTypeDisk``; see Section :galacticus-ref:`ComponentMassTypes`). If the component is active, it should also set ``specificAngularMomentum`` to reflect the specific angular momentum (in km s\ :math:`^{-1}` Mpc) of the component (at whatever point in its profile the radius is required; this is required only if ``specificAngularMomentumRequired`` is true) and should point the four procedure pointers to routines which get and set the radius and circular velocity properties of the component (which should have the standard form for component get and set methods). It is acceptable for the set procedures to point to dummy routines.

The galactic structure radii solver routines will use this information to determine (and set) the radius and circular velocity of the component. An advantage of this approach is that different radii solver methods can all use this same system, ensuring that just a single interface is needed in each component.

Analytic Solvers
----------------

Galacticus typically solves the system of :term:`ODE` which describe the evolution of a galaxy using a numerical solver. However, in some situations analytic solutions may be available. Galacticus allows for the use of analytic solvers in such cases. Currently available analytic solvers are described below.

Simple Disk Satellites
~~~~~~~~~~~~~~~~~~~~~~

This solver, which can be activated by setting ``[diskVerySimpleUseAnalyticSolver]``\ :math:`=`\ ``true`` is applicable to satellite systems in the very simple disk component. In particular, the following conditions must be met for it to be valid:

* The hot halo component must be of type ``verySimple`` or ``verySimpleDelayed``;
* The satellite component must be of type ``verySimple``;
* The disk component must be of type ``verySimple``;
* Properties of all other components must not change for satellite galaxies (with the exception of the ``time`` property of the basic component), or else the component must be ``null``;
* Disk star formation and outflow rates must scale as :math:`M_\mathrm{gas}/\tau` where :math:`M_\mathrm{gas}` is the instantaneous mass of gas in the galaxy disk, and :math:`\tau` is a fixed timescale for any given satellite galaxy.

Under these conditions, the flow of mass between gas, stellar, and outflowed phases is analytically solvable.

Subsystems
----------

This section describes some of the subsystems within Galacticus that support various physical entities or processes.

.. _manual-sec-KeplerOrbits:

Kepler Orbits
~~~~~~~~~~~~~

The ``keplerOrbit`` object (provided by the ``Kepler_Orbits_Structure`` module) stores the parameters of a single Keplerian orbit. It internally handles computation of additional/alternate orbital parameters once an orbit has been fully defined. Currently, the orientation of orbits (i.e. the unit vector normal to the orbital plane and the argument of periapsis) is not tracked. As such, orbits are fully defined by three parameters (in addition to the masses of the orbiting bodies). The following limitations presently apply to the ``keplerOrbit`` object:

* If an orbit is overdefined (i.e. if more than three parameters are set manually) no checking is performed to ensure that the parameters are consistent with a Keplerian orbit;
* Not all interconversions between parameters are supported\ [#]_. If a conversion cannot be performed, an error message will be given.

A ``keplerOrbit`` object can be reset by calling the ``reset()`` method, and its defined/undefined status can be tested with the ``isDefined()`` method or asserted with the ``assertIsDefined()`` method. The following orbital parameters are supported, each method returning the value of the parameter and a corresponding method suffixed with ``Set`` can be used to set the parameter: ``radius``, ``velocityRadial``, ``velocityTangential``, ``energy``, ``angularMomentum``, ``eccentricity``, ``semiMajorAxis``, ``radiusPericenter``, ``radiusApocenter``. Additionally, the masses of the orbiting bodies are provided by the ``hostMass()`` and ``reducedMassSpecific()``\ :math:`=M_\mathrm{host}/(M_\mathrm{host}+M_\mathrm{satellite})` methods. Finally, the ``velocityScale()`` method returns :math:`\mathrm{G}M_\mathrm{host}/r` where :math:`r` is the radius of the orbit.

.. _manual-sec-ChemicalSubsystem:

Chemicals
~~~~~~~~~

The chemicals subsystem provides both a interface to a database of known chemicals (allowing their physical properties to be queried) and a structure to store abundances/masses/etc. of the set of chemicals being tracked in Galacticus. The name "chemicals" is used to denote any chemical species that might be involved in reactions, including molecules, atoms, atomic and molecular ions and electrons.

Chemical Database
^^^^^^^^^^^^^^^^^

The file ``data/Chemical_Database.cml`` contains a database of chemicals that can currently be used by Galacticus. It uses a simplified version of the `Chemical Markup Language <http://www.xml-cml.org>`_ to describe chemicals. An excerpt from the database is shown below:

.. code-block:: none

    <list>
     <chemical>
      <id>MolecularHydrogenAnion</id>
      <formalCharge>-1</formalCharge>
      <atomArray>
        <atom>
         <id>1</id>
         <elementType>H</elementType>
        </atom>
        <atom>
         <id>2</id>
         <elementType>H</elementType>
        </atom>
      </atomArray>
      <bondArray>
        <bond>
         <atomRefs2>1 2</atomRefs2>
         <order>1</order>
        </bond>
      </bondArray>
     </chemical>
     .
     .
     .
    </list>

The database contains a ``list`` of chemicals, each contained within a ``chemical`` element. The ``id`` element provides a label for the chemical (usually a descriptive name with no white space). The ``formalCharge`` element gives the charge of the chemical in units of the elementary charge. The chemical is then describe by a list of atoms and bonds inside ``atomArray`` and ``bondArray`` elements respectively. The ``atomArray`` can contain any number of ``atom`` elements, which should describe each atom in the chemical giving it a unique ``id`` number and an ``elementType``, which is the short one or two letter label for the element (e.g. H, Ni, etc.). The ``bondArray`` should contain a ``bond`` entry for each atomic bond, which itself contains a ``atomRefs2`` element giving the IDs of the two atoms participating in the bond and an ``order`` element which gives the order of the bond (e.g. "1" for a single bond).

Chemical Structure
^^^^^^^^^^^^^^^^^^

Within Galacticus a chemical is represented using the ``chemicalStructure`` type which is provided by the ``Chemical_Structures`` module. A ``chemicalStructure`` object can be assigned a particular chemical by retrieving that chemical from the database using:

.. code-block:: none

    call myChemical%retrieve("chemicalID")

where ``chemicalID`` is the ID of the chemical in the database. Any chemical can be exported to a CML file using

.. code-block:: none

    call myChemical%export(fileName)

where ``fileName`` gives the name of the file to which to export.

Once assigned a chemical, basic properties such as mass and charge (in atomic units) can be accessed using ``myChemical%mass`` and ``myChemical%charge`` respectively. The mass is computed from the known atomic masses of the constituent atoms of the chemical.

Chemical Abundances
^^^^^^^^^^^^^^^^^^^

Within Galacticus a set of abundances (or masses, or densities…) for all chemicals being tracked, as specified by the ``[chemicalsToTrack]`` input parameter, is stored within a ``chemicalAbundancesStructure`` type, as provided by the ``Chemical_Abundances_Structure`` module. The structure provides interfaces for setting and retrieving the abundance of a given chemical species, to pack/unpack all chemicals to/from an array, to convert from mass-weighted to number-weighted quantities and to multiply and divide the chemicals abundances by a given amount. Additionally, the ``Chemical_Abundances_Structure`` module provides functions which provide a count of the number of chemicals tracked, to look up the index of a chemical array representation from its name, and to retrieve the name of a given chemical.

.. _manual-sec-Coordinates:

Coordinates
~~~~~~~~~~~

The :galacticus-class:`coordinate` class describes a position in three-dimensional space. Each extension of this class (currently, :galacticus-class:`coordinateCartesian`, :galacticus-class:`coordinateCylindrical`, and :galacticus-class:`coordinateSpherical`) supply methods to convert to and from Cartesian coordinates. The assignment operator (``=``) is overloaded such that coordinate objects of any class can be assigned to any other class and conversion to the appropriate coordinate system will happen automatically. A function accepting a :galacticus-class:`coordinate` object can therefore convert it to, for example, spherical coordinates simply using

.. code-block:: none

    class(coordinate         ), intent(in) :: coordinates
    type (coordinateSpherical)             :: coordinatesSpherical
    coordinatesSpherical=coordinates

and thereby allow a position to be passed to it in any coordinate system.

.. [#] Data objects in components can be real, integer, boolean or of derived type. For derived types, currently ``history``, ``abundances``, ``chemicals``, and ``keplerOrbit`` are supported. Adding additional derived types is possible, providing that the type supports the required methods for output, serialization, etc. Data objects can currently be scalar or rank-1 arrays.
.. [#] Or some other value if a ``classDefault`` has been specified (see Section :galacticus-ref:`ComponentImplement`).
.. [#] The ``keplerOrbit`` object works by trying to convert to the combination radius, radial and tangential velocities. Once these are defined, all other parameters can be computed. However, for orbits defined in terms of other parameters, the ``keplerOrbit`` object does not know how to convert from every such combination of parameters.
