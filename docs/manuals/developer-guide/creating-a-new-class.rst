.. _manual-sec-creating-a-new-class:

Creating a New Class
====================

Overview of Classes
-------------------

Galacticus uses an `object-oriented programming <https://en.wikipedia.org/wiki/Object-oriented_programming>`_ model. This means that most of the functionality in Galacticus is provided by "`classes <https://en.wikipedia.org/wiki/Class-based_programming>`_".

The easiest way to think about these as used by Galacticus is that each of these represents a "class" of calculation or physical process (e.g. a class might represent the rate of star formation in a galaxy).

The class itself defines *only* the interfaces of the functions (also known as "*methods*") needed for that calculation - the arguments to those functions, and the types of result they return. For example, in the case of a star formation rate class, it might define a method "``rate``" which takes as an argument a ``component`` object (representing, for example, the disk or spheroid of a galaxy) and returns a single number, the star formation rate in that component.

Importantly, the class itself *does not* actually do the calculation - all it does is define the interface to the method that *will* do the calculation. This allows us to separate interface and implementation - any other part of the Galacticus code can now make use of the star formation rate class without having to know anything about how the star formation rate is computed - it just passes the appropriate ``component`` to the star formation rate class and gets back an answer.

The actual calculation (of star formation rate in this case) is carried out by an "*implementation*" of the class - an object that provides the actual functions implementing the methods that the class defined. There can be many different implementations of each class, which can be selected between via Galacticus' parameter files. This is where the flexibility and modularity of Galacticus arises - you can switch between different implementations of any calculation simply by changing the parameter file.

Creating a new class
--------------------

As you might expect, creating a new class involves first creating the class itself, and then creating an implementation (a class without an implementation is not very useful). These go in two separate files, which we'll see examples of below. For the purposes of this tutorial, we're going to use as an example the star formation timescale class.

We'll begin by looking at the file that defines the class itself.

The class file
~~~~~~~~~~~~~~

Below is the content of the file `source/star_formation/timescales/_class.F90 <https://github.com/galacticusorg/galacticus/blob/master/source/star_formation/timescales/_class.F90>`_. Source files are organized into a directory hierarchy under the ``source/`` folder, in which each concept (here "star formation" and "timescales" are the two concepts relevant to this class) becomes a directory level, and words within a concept are separated by ``_``. The abstract base class for a set of implementations lives in a file named ``_class.F90`` within the class's directory (here ``source/star_formation/timescales/``).

.. code-block:: fortran

   !!{
   Contains a module which provides a class that implements timescales for star formation.
   !!}

   module Star_Formation_Timescales
     !!{
     Provides a class that implements calculations of timescales for star formation.
     !!}
     use :: Galacticus_Nodes, only : nodeComponent
     private

     !![
     <functionClass>
      <name>starFormationTimescale</name>
      <descriptiveName>Timescales for star formation</descriptiveName>
      <description>
       Class providing models of timescales for star formation.
      </description>
      <default>dynamicalTime</default>
      <method name="timescale" >
       <description>Returns the timescale (in Gyr) for star formation in the provided ``component``.</description>
       <type>double precision</type>
       <pass>yes</pass>
       <argument>class(nodeComponent), intent(inout) :: component</argument>
      </method>
     </functionClass>
     !!]

   end module Star_Formation_Timescales

We will now go through each component of this file one at a time.

Opening comment
^^^^^^^^^^^^^^^

This section:

.. code-block:: fortran

   !!{
   Provides a class that implements timescales for star formation.
   !!}

is just a comment that tells someone reading this file what to expect to find in this file. Content inside ``!!{`` ``!!}`` blocks in Galacticus should be in reStructuredText format - it will be included into our ReadTheDocs documentation.

Module opener
^^^^^^^^^^^^^

Every class is placed inside its own `module <https://fortran-lang.org/learn/best_practices/modules_programs/>`_ - this allows it to be imported into other parts of Galacticus wherever it's needed. The module is begun by this line:

.. code-block:: fortran

   module Star_Formation_Timescales

We give the module a descriptive name - usually this is built from the concepts in the file's directory path, with words capitalized, separated by ``_``, and pluralized. So, ``source/star_formation/timescales/_class.F90`` becomes ``Star_Formation_Timescales``

The module is closed at the end of the file by the line:

.. code-block:: fortran

   end module Star_Formation_Timescales

Module comment
^^^^^^^^^^^^^^

Next, we have:

.. code-block:: fortran

     !!{
     Provides a class that implements calculations of timescales for star formation.
     !!}

This is another comment, this time associated with the ``module`` instead of the entire file. Once again, it just provides a brief description of what this module does.

Imports
^^^^^^^

If the class is going to make use of any other objects or classes, we need to import them here. In our example case we want to be able to pass a ``nodeComponent`` object as an argument to one of the class methods. So, we import that with a `use <https://fortran-lang.org/learn/best_practices/modules_programs/>`_ statement.

.. code-block:: fortran

     use :: Galacticus_Nodes, only : nodeComponent

Note that we import ``only`` the object we want. The ``Galacticus_Nodes`` module contains *many* objects and variables. It's good practice to only import the one we want to avoid any possible conflicts with objects/variables that we may import from other modules.

Exports
^^^^^^^

The next line:

.. code-block:: fortran

     private

declares that all content of this module is private by default - it cannot be exported to other parts of the code unless explicitly made public. Setting this default is also good practice - it prevents the rest of the code from accessing module content that it shouldn't. Of course, our module needs to export *something* in order to be useful! Galacticus knows that it needs to exploit that class that we are about to create in this module, and will do so automatically.

The class!
^^^^^^^^^^

Finally we get to define our class. This is enclosed inside a ``!![`` ``!!]`` section. Code inside such sections is interpreted as `XML <https://en.wikipedia.org/wiki/XML>`_, and is used when Galacticus is built to automatically generate lots of extra code.

This section begins with:

.. code-block:: xml

     <functionClass>

which simply tells Galacticus that we want to create a new class (referred to as a ``functionClass`` internally).

Next, we give a name for our class:

.. code-block:: xml

      <name>starFormationTimescale</name>
      <descriptiveName>Timescales for star formation</descriptiveName>

The ``name`` element here should be a simple, descriptive name. It should contain the same "concepts" used in naming the file and module, but uses ``camelCaseConvention`` where the first word is lower case, and the rest have their first letter uppercased - there are no separators between words. The ``descriptiveName`` element gives a descriptive name for this class - this will appear in the documentation so should be something that a reader will readily understand.

Next we provide a description for the class - this will be incorporated into the ReadTheDocs documentation.

.. code-block:: xml

      <description>
       Class providing models of timescales for star formation.
      </description>

Every class has a default implementation. If a class is needed by Galacticus, but the parameter file does not specify which implementation to use, the default implementation is chosen. In this case, we make the ``dynamicalTime`` implementation the default option.

.. code-block:: xml

      <default>dynamicalTime</default>

We now get to define any methods that we want our class to have. Our example class has just one, but we can include as many as we want (just add additional ``method`` sections below the first). Our ``method`` section looks like this:

.. code-block:: xml

      <method name="timescale" >
       <description>Returns the timescale (in Gyr) for star formation in the provided ``component``.</description>
       <type>double precision</type>
       <pass>yes</pass>
       <argument>class(nodeComponent), intent(inout) :: component</argument>
      </method>

In the first line, we give a ``name`` for our method - in this case that name is ``timescale``. This name is how other parts of the code will call this method. So, it should be descriptive of what the method computes. In this case, our method computes the timescale for star formation - so ``timescale`` is an appropriate name.

Next we include a ``description`` - this should be in reStructuredText format and will be incorporated into the ReadTheDocs documentation - it should explain clearly what the method does (and it's good to include details of what units are expected for arguments and the return value).

The ``type`` element specifies the type of quantity that the method should return. In our example case, we want to return a timescale, so a single ``double precision`` value. So, we just set the ``type`` to ``double precision``. Any valid Fortran type can be included here. For methods that don't return any value, set this to ``void`` - these will be implemented as ``subroutines`` in the generated code.

Next, the ``pass`` arguments specifies whether the object of this class (i.e. the specific object created from whichever implementation of the class was chosen at run time) should be passed as the first argument to the method. Almost always this should be ``yes``.

Finally, we need to specify what arguments our method expects. Arguments are given by ``argument`` elements - you can include multiple of these if needed. In our case we have:

.. code-block:: xml

       <argument>class(nodeComponent), intent(inout) :: component</argument>

This specifies that we want an argument of the ``nodeComponent`` ``class`` (which is used to represent components of galaxies in Galacticus). The ``intent(inout)`` specifies that the argument must have a defined value when passed to the method, and may be modified (and returned) by the method. Last, we give the name of the argument, ``component`` in this case. Any valid Fortran argument definition is allowed here.

Last, we close the ``functionClass`` element.

.. code-block:: xml

     </functionClass>

That's all we need to define our class. When you compile Galacticus, the ``functionClass`` definition is used to generate the class itself in Fortran code. Additionally, a lot of useful infrastructure associated with your class is automatically built - including functions that allow duplication of class objects across parallel threads, serialization/deserialization of the class to/from files, descriptor strings used to ensure tabulated results are unique, and functionality to allow the class to be built from the contents of the parameter file provided at run time. You should never need to think about this infrastructure - it's created automatically and should "just work".

Next, we'll look at creating an implementation of this class.

The implementation file
~~~~~~~~~~~~~~~~~~~~~~~

Below is the content of the file `source/star_formation/timescales/halo_scaling.F90 <https://github.com/galacticusorg/galacticus/blob/master/source/star_formation/timescales/halo_scaling.F90>`_. The file is placed in the same directory as the corresponding class file (``source/star_formation/timescales/``), and is named for the concept specific to this implementation - ``halo_scaling`` in this example.

.. code-block:: fortran

     !!{
     Implementation of a timescale for star formation which scales with the circular velocity of the host halo.
     !!}

     use :: Cosmology_Functions    , only : cosmologyFunctionsClass
     use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

     !![
     <starFormationTimescale name="starFormationTimescaleHaloScaling">
      <description>
       A star formation timescale class in which the timescale scales with halo properties. Specifically,
       \begin{equation}
        \tau_\star = \tau_\mathrm{\star,0} \left( {V_\mathrm{vir} \over 200\hbox{km/s}} \right)^{\alpha_\star} (1+z)^{\beta_\star},
       \end{equation}
       where $\tau_\mathrm{\star,0}=${\normalfont \ttfamily [timescale]}, $\alpha_\star=${\normalfont \ttfamily
       [exponentVelocityVirial]}, and $\beta_\star=${\normalfont \ttfamily [exponentRedshift]}.
      </description>
     </starFormationTimescale>
     !!]
     type, extends(starFormationTimescaleClass) :: starFormationTimescaleHaloScaling
        !!{
        Implementation of a haloScaling timescale for star formation.
        !!}
        private
        class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_           => null()
        class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_          => null()
        double precision                                    :: expansionFactorFactorPrevious          , exponentVelocityVirial , &
             &                                                 exponentRedshift                       , timescaleNormalization , &
             &                                                 timescaleStored                        , velocityPrevious       , &
             &                                                 velocityFactorPrevious                 , expansionFactorPrevious, &
             &                                                 timeScale_
        logical                                             :: timescaleComputed
        integer         (kind_int8                        ) :: lastUniqueID
      contains
        !![
        <methods>
          <method description="Reset memoized calculations." method="calculationReset" />
        </methods>
        !!]
        final     ::                     haloScalingDestructor
        procedure :: autoHook         => haloScalingAutoHook
        procedure :: timescale        => haloScalingTimescale
        procedure :: calculationReset => haloScalingCalculationReset
     end type starFormationTimescaleHaloScaling

     interface starFormationTimescaleHaloScaling
        !!{
        Constructors for the {\normalfont \ttfamily haloScaling} timescale for star formation class.
        !!}
        module procedure haloScalingConstructorParameters
        module procedure haloScalingConstructorInternal
     end interface starFormationTimescaleHaloScaling

     double precision, parameter :: velocityVirialNormalization=200.0d0

   contains

     function haloScalingConstructorParameters(parameters) result(self)
       !!{
       Constructor for the {\normalfont \ttfamily haloScaling} timescale for star formation feedback class which takes a
       parameter set as input.
       !!}
       use :: Input_Parameters, only : inputParameter, inputParameters
       implicit none
       type            (starFormationTimescaleHaloScaling)                :: self
       type            (inputParameters                  ), intent(inout) :: parameters
       class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
       class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
       double precision                                                   :: timescale           , exponentVelocityVirial, &
            &                                                                exponentRedshift

       !![
       <inputParameter>
         <name>timescale</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The timescale for star formation in the halo scaling timescale model.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>exponentVelocityVirial</name>
         <defaultValue>0.0d0</defaultValue>
         <description>The exponent of virial velocity in the timescale for star formation in the halo scaling timescale model.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>exponentRedshift</name>
         <defaultValue>0.0d0</defaultValue>
         <description>The exponent of redshift in the timescale for star formation in the halo scaling timescale model.</description>
         <source>parameters</source>
       </inputParameter>
       <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
       <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
       !!]
       self=starFormationTimescaleHaloScaling(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_)
       !![
       <inputParametersValidate source="parameters"/>
       <objectDestructor name="cosmologyFunctions_" />
       <objectDestructor name="darkMatterHaloScale_"/>
       !!]
       return
     end function haloScalingConstructorParameters

     function haloScalingConstructorInternal(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_) result(self)
       !!{
       Internal constructor for the {\normalfont \ttfamily haloScaling} timescale for star formation class.
       !!}
       implicit none
       type            (starFormationTimescaleHaloScaling)                        :: self
       double precision                                   , intent(in   )         :: timescale           , exponentVelocityVirial, &
            &                                                                        exponentRedshift
       class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
       class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
       !![
       <constructorAssign variables="exponentVelocityVirial, exponentRedshift, *cosmologyFunctions_, *darkMatterHaloScale_"/>
       !!]

       self%lastUniqueID                 =-1_kind_int8
       self%timescaleComputed            =.false.
       self%velocityPrevious             =-1.0d0
       self%velocityFactorPrevious       =-1.0d0
       self%expansionFactorPrevious      =-1.0d0
       self%expansionFactorFactorPrevious=-1.0d0
       ! Compute the normalization of the timescale.
       self%timeScale_            =+timescale
       self%timeScaleNormalization=+timescale                                                &
            &                      /velocityVirialNormalization**self%exponentVelocityVirial
       return
     end function haloScalingConstructorInternal

     subroutine haloScalingAutoHook(self)
       !!{
       Attach to the calculation reset event.
       !!}
       use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
       implicit none
       class(starFormationTimescaleHaloScaling), intent(inout) :: self

       call calculationResetEvent%attach(self,haloScalingCalculationReset,openMPThreadBindingAllLevels,label='starFormationTimescaleHaloScaling')
       return
     end subroutine haloScalingAutoHook

     subroutine haloScalingDestructor(self)
       !!{
       Destructor for the {\normalfont \ttfamily haloScaling} timescale for star formation class.
       !!}
       use :: Events_Hooks, only : calculationResetEvent
       implicit none
       type(starFormationTimescaleHaloScaling), intent(inout) :: self

       !![
       <objectDestructor name="self%cosmologyFunctions_" />
       <objectDestructor name="self%darkMatterHaloScale_"/>
       !!]
       if (calculationResetEvent%isAttached(self,haloScalingCalculationReset)) call calculationResetEvent%detach(self,haloScalingCalculationReset)
       return
     end subroutine haloScalingDestructor

     subroutine haloScalingCalculationReset(self,node,uniqueID)
       !!{
       Reset the halo scaling star formation timescale calculation.
       !!}
       use :: Galacticus_Nodes, only : treeNode
       use :: Kind_Numbers    , only : kind_int8
       implicit none
       class  (starFormationTimescaleHaloScaling), intent(inout) :: self
       type   (treeNode                         ), intent(inout) :: node
       integer(kind_int8                        ), intent(in   ) :: uniqueID
       !$GLC attributes unused :: node

       self%timescaleComputed=.false.
       self%lastUniqueID     =uniqueID
       return
     end subroutine haloScalingCalculationReset

     double precision function haloScalingTimescale(self,component) result(timescale)
       !!{
       Returns the timescale (in Gyr) for star formation in the given {\normalfont \ttfamily component} in the halo scaling
       timescale model.
       !!}
       use :: Galacticus_Nodes, only : nodeComponentBasic
       implicit none
       class           (starFormationTimescaleHaloScaling), intent(inout) :: self
       class           (nodeComponent                    ), intent(inout) :: component
       class           (nodeComponentBasic               ), pointer       :: basic
       double precision                                                   :: expansionFactor, velocityVirial

       ! Check if node differs from previous one for which we performed calculations.
       if (component%hostNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(component%hostNode,component%hostNode%uniqueID())
       ! Compute the timescale if necessary.
       if (.not.self%timescaleComputed) then
          ! Get virial velocity and expansion factor.
          basic           => component%hostNode%basic                               (                    )
          velocityVirial  =  self              %darkMatterHaloScale_%velocityVirial (component%hostNode  )
          expansionFactor =  self              %cosmologyFunctions_ %expansionFactor(basic    %time    ())
          ! Compute the velocity factor.
          if (velocityVirial /= self%velocityPrevious) then
              self%velocityPrevious      =velocityVirial
              self%velocityFactorPrevious=velocityVirial**self%exponentVelocityVirial
          end if
          ! Compute the expansion-factor factor.
          if (expansionFactor /= self%expansionFactorPrevious) then
             self%expansionFactorPrevious      =      expansionFactor
             self%expansionFactorFactorPrevious=1.0d0/expansionFactor**self%exponentRedshift
          end if
          ! Computed the timescale.
          self%timescaleStored=+self%timeScaleNormalization        &
               &               *self%velocityFactorPrevious        &
               &               *self%expansionFactorFactorPrevious
          ! Record that the timescale is now computed.
          self%timescaleComputed=.true.
       end if
       ! Return the stored timescale.
       timescale=self%timescaleStored
       return
     end function haloScalingTimescale

We'll now go through each part of this file step-by-step.

Opening comment
^^^^^^^^^^^^^^^

Once again, we start with a comment that describes what this file contains:

.. code-block:: fortran

     !!{
     Implementation of a timescale for star formation which scales with the circular velocity of the host halo.

     !!}

This is in reStructuredText format, and gets incorporated into our ReadTheDocs documentation.

Implementation definition
^^^^^^^^^^^^^^^^^^^^^^^^^

We next include a section that tells Galacticus that we're going to include an implementation of a class in this file:

.. code-block:: xml

     <starFormationTimescale name="starFormationTimescaleHaloScaling">
      <description>
       A star formation timescale class in which the timescale scales with halo properties. Specifically,

       .. math::

          \tau_\star = \tau_\mathrm{\star,0} \left( {V_\mathrm{vir} \over 200\hbox{km/s}} \right)^{\alpha_\star} (1+z)^{\beta_\star},

       where :math:`\tau_\mathrm{\star,0}=` ``[timescale]``, :math:`\alpha_\star=` ``[exponentVelocityVirial]``, and
       :math:`\beta_\star=` ``[exponentRedshift]``.
      </description>
     </starFormationTimescale>

This is in XML format (because it's inside a ``!![`` ``!!]`` section). The opening element of the XML *must* match the name of the class, as you defined it in the class file (``starFormationTimescale`` in this case), and then we give the name for this implementation, which should always start with the class time, followed by something descriptive of the concepts of this implementation. In our example, we take the class name ``starFormationTimescale`` and add ``HaloScaling`` to describe this implementation, giving us ``starFormationTimescaleHaloScaling``.

The remainder of this section is a ``description`` element, which should give a detailed description of this implementation - the physics it implements, any relevant equations and references. This should be in reStructuredText format, as it will be incorporated into the ReadTheDocs documentation.

The implementation ``type``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we must define a ``type`` for our implementation. This is the data type that will represent our implementation internally. It will contain any variables that our implementation might need (e.g. to store values of parameters, or tables of results that it needs to refer to, helper objects (e.g. root finders, interpolators), and pointers to other classes that it wants to make use of.

.. code-block:: fortran

     type, extends(starFormationTimescaleClass) :: starFormationTimescaleHaloScaling
        !!{
        Implementation of a haloScaling timescale for star formation.
        !!}
        private
        class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_           => null()
        class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_          => null()
        double precision                                    :: expansionFactorFactorPrevious          , exponentVelocityVirial , &
             &                                                 exponentRedshift                       , timescaleNormalization , &
             &                                                 timescaleStored                        , velocityPrevious       , &
             &                                                 velocityFactorPrevious                 , expansionFactorPrevious, &
             &                                                 timeScale_
        logical                                             :: timescaleComputed
        integer         (kind_int8                        ) :: lastUniqueID
      contains
        !![
        <methods>
          <method description="Reset memoized calculations." method="calculationReset" />
        </methods>
        !!]
        final     ::                     haloScalingDestructor
        procedure :: autoHook         => haloScalingAutoHook
        procedure :: timescale        => haloScalingTimescale
        procedure :: calculationReset => haloScalingCalculationReset
     end type starFormationTimescaleHaloScaling

We begin by opening the type, and giving it a name:

.. code-block:: fortran

     type, extends(starFormationTimescaleClass) :: starFormationTimescaleHaloScaling

Note that the ``extends(starFormationTimescaleClass)`` states that this implementation is an extension of the ``starFormationTimescale`` class that we created in the first file. The internal object representing this class is always named in this way - the name that we gave for our class, ``starFormationTimescale``, plus the suffix ``Class``. Then, the ``starFormationTimescaleHaloScaling`` is the name for our implementation. This *must* begin with the name of the class, ``starFormationTimescale`` in this case, followed by words describing the concepts of this implementation, all following ``camelCaseFormat``.

.. tip::

   *Advanced:* It's possible to have one implementation be an extension of another - reusing most of the functionality of that first implementation, but adding/overriding some behavior. In such cases instead of ``extends(starFormationTimescaleClass)`` we would extend the first implementation. So, for example, if we wanted to create a second implementation that extended our current one, we would use ``extends(starFormationTimescaleHaloScaling)``.

Next we have a description of the implementation - this can be brief - it is included into the documentation on this implementation object:

.. code-block:: fortran

        !!{
        Implementation of a haloScaling timescale for star formation.
        !!}

After this we declare any variables which we want to belong to our implementation. Every time an object is created from our implementation it will contains these variables - this is the "state" or "content" of the implementation.

.. code-block:: fortran

        private
        class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_           => null()
        class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_          => null()
        double precision                                    :: expansionFactorFactorPrevious          , exponentVelocityVirial , &
             &                                                 exponentRedshift                       , timescaleNormalization , &
             &                                                 timescaleStored                        , velocityPrevious       , &
             &                                                 velocityFactorPrevious                 , expansionFactorPrevious, &
             &                                                 timeScale_
        logical                                             :: timescaleComputed
        integer         (kind_int8                        ) :: lastUniqueID

The first line, ``private``, states that all content should be private to the object - the rest of the code can't directly see or modify this content. This is good practice - it allows us to keep the internal state of our implementation fully under its own control.

The remaining lines declare variables and pointers needed by our implementation, following standard Fortran syntax. We won't go into the specifics of the variables for this particular implementation in this tutorial, but we'll note a few things:

* Entries such as ``class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()`` are pointers to other classes. Often, an implementation will need to make use of some other class to perform part of its calculation. Here, for example, the implementation needs access to functions related to cosmology (the relation between time and redshift specifically). This is done by getting a pointer to such an object (we'll see how this works below). Note the ``null()`` - this ensures that the pointer is set to be null by default - this avoids potential errors that can occur if a pointer was otherwise left pointing to some random piece of memory. These pointers to other classes are typically given the same name as that class, but with a trailing underscore.
* This implementation stores some variables, e.g. ``exponentVelocityVirial``, ``exponentRedshift``, that store the parameters of the model that it implements. We'll see how these are determined from the parameter file below.
* Other variables, such as ``timescaleNormalization`` are quantities that the implementation computes once, and stores for continued re-use. This is often done to speed up calculations (by avoiding computing the same quantity many times).
* This implementation also has variables such as ``velocityPrevious`` which are used for `memoization <https://en.wikipedia.org/wiki/Memoization>`_ in which we store the result of the previous calculation of some function in case we are asked for that same value again - allowing us to re-use the stored result (again, to speed things up).
* Also worth noting is the ``timeScale_`` variable. In this example this stores the normalization of the star formation timescale. Since it has the same name as the ``timeScale`` method, we add a trailing underscore to distinguish it.

After defining the content of the implementation, we move on to defining the methods that it implements. These are separated from the content by the line

.. code-block:: fortran

      contains

The methods section then begins with descriptions of any methods unique to this implementation (i.e. that were not already defined in the parent class):

.. code-block:: fortran

        !![
        <methods>
          <method description="Reset memoized calculations." method="calculationReset" />
        </methods>
        !!]

Most implementations won't have any unique methods, in which case this section can be omitted. The descriptions are placed with a ``methods`` XML section, in which each method is defined by a ``method`` element, with two attributes, ``method`` which gives the name of the method, and ``description`` that describes what that method does. These descriptions are incorporated into the ReadTheDocs documentation.

After this we define the methods - in general this just connect the name of the method to the function that provides the functionality:

.. code-block:: fortran

        final     ::                     haloScalingDestructor
        procedure :: autoHook         => haloScalingAutoHook
        procedure :: timescale        => haloScalingTimescale
        procedure :: calculationReset => haloScalingCalculationReset

In this example we have:

* ``final``: The ``final`` line is a special case which gives the function that should be called to destroy objects of this type (we'll see more about this below).
* ``autoHook``: Another special case - the ``autoHook`` function is automatically called by Galacticus after an object is created or copied. We'll see more on this below.
* ``timescale``: Finally(!) this is the method we defined in the parent class - that will actually do the work of computing a star formation timescale. In this case, we specify that this method is implemented by the function ``haloScalingTimescale``
* ``calculationReset``: Last we have a new method, unique to this implementation. We'll learn more about this new method below.

Last, we close the type:

.. code-block:: fortran

     end type starFormationTimescaleHaloScaling

List of constructors
^^^^^^^^^^^^^^^^^^^^

After defining our implementation type, we need to give a list of functions that can be used to build objects of this type - these are known as `constructors <https://en.wikipedia.org/wiki/Constructor_(object-oriented_programming)>`_. We'll see more about how constructors work below, for now, we just list the names of the constructor function:

.. code-block:: fortran

     interface starFormationTimescaleHaloScaling
        !!{
        Constructors for the {\normalfont \ttfamily haloScaling} timescale for star formation class.
        !!}
        module procedure haloScalingConstructorParameters
        module procedure haloScalingConstructorInternal
     end interface starFormationTimescaleHaloScaling

This list is placed inside an ``interface`` section with a name matched to the name of our implementation. The ``interface`` section contains a short description, followed by the list of constructors, each on a line starting with ``module procedure``.

Implementation-scope parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It can often be useful to define parameters that are shared between multiple functions in our class. This is the place to include them. In this example we have:

.. code-block:: fortran

     double precision, parameter :: velocityVirialNormalization=200.0d0

which defines a velocity scale at which the star formation timescale model is normalized.

End of the definition, beginning of the functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We're now done with the preamble (defining the type, listing constructors, etc.). The line:

.. code-block:: fortran

   contains

indicates that we've reached the end of this, and are ready to start the actual functions of our implementation.

Constructors
^^^^^^^^^^^^

The first functions we define are constructors. `Above <#list-of-constructors>`_, in our list of constructors, we named two constructor functions. Each of these must return as its result an object of our implementation type, ``starFormationTimescaleHaloScaling``, fully initialized so that it is ready to use. Constructor functions are named by starting with the non-class part of our implementation name (so, take ``starFormationTimescaleHaloScaling``, remove the ``starFormationTimescale``, leaving just ``haloScaling``).

Two constructors is typical, but only the first is actually required. In our example, the two constructor functions are:

#. ``haloScalingConstructorParameters`` - This one is required for every implementation - it will build an object of our implementation from the parameter file given to Galacticus (e.g. reading values of parameters from it).
#. ``haloScalingConstructorInternal`` - This one is not required, but is often useful to have. It can be used internally to build objects of our type as needed. Typically, all relevant parameters needed to set up the object are passed to this constructor as arguments.

We can now look at the first constructor function:

.. code-block:: fortran

     function haloScalingConstructorParameters(parameters) result(self)
       !!{
       Constructor for the {\normalfont \ttfamily haloScaling} timescale for star formation feedback class which takes a
       parameter set as input.
       !!}
       use :: Input_Parameters, only : inputParameter, inputParameters
       implicit none
       type            (starFormationTimescaleHaloScaling)                :: self
       type            (inputParameters                  ), intent(inout) :: parameters
       class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
       class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
       double precision                                                   :: timescale           , exponentVelocityVirial, &
            &                                                                exponentRedshift

       !![
       <inputParameter>
         <name>timescale</name>
         <defaultValue>1.0d0</defaultValue>
         <description>The timescale for star formation in the halo scaling timescale model.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>exponentVelocityVirial</name>
         <defaultValue>0.0d0</defaultValue>
         <description>The exponent of virial velocity in the timescale for star formation in the halo scaling timescale model.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>exponentRedshift</name>
         <defaultValue>0.0d0</defaultValue>
         <description>The exponent of redshift in the timescale for star formation in the halo scaling timescale model.</description>
         <source>parameters</source>
       </inputParameter>
       <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
       <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
       !!]
       self=starFormationTimescaleHaloScaling(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_)
       !![
       <inputParametersValidate source="parameters"/>
       <objectDestructor name="cosmologyFunctions_" />
       <objectDestructor name="darkMatterHaloScale_"/>
       !!]
       return
     end function haloScalingConstructorParameters

We begin by defining the function name and interface:

.. code-block:: fortran

     function haloScalingConstructorParameters(parameters) result(self)

Here we define the arguments to this function - for this required constructor there must be precisely one argument, called ``parameters``. We also define the name of the variable that will hold the result of the constructor (the object of our implementation that it builds) - this should always be called ``self``, so we have ``result(self)``.

Next we have the usual description for inclusion into the ReadTheDocs document. This is followed by any module imports that we might need:

.. code-block:: fortran

       use :: Input_Parameters, only : inputParameters

Here we just import the ``inputParameters`` class from the ``Input_Parameters`` module. This is the type of our ``parameters`` argument, so we must include this.

We then define all variables needed by our constructor function:

.. code-block:: fortran

       implicit none
       type            (starFormationTimescaleHaloScaling)                :: self
       type            (inputParameters                  ), intent(inout) :: parameters
       class           (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
       class           (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
       double precision                                                   :: timescale           , exponentVelocityVirial, &
            &                                                                exponentRedshift

The ``implicit none`` just enforces that we must explicitly declare all variables (this is good practice). Then we have:

.. code-block:: fortran

       type            (starFormationTimescaleHaloScaling)                :: self
       type            (inputParameters                  ), intent(inout) :: parameters

The first of these lines defines our ``self`` object - the thing we will construct. It must be defined as a ``type()`` of our implementation type, ``starFormationTimescaleHaloScaling``. The second line defines the argument - for this constructor that is always the ``parameters`` argument and *must* have the form given above.

After these, we typically have any variables that we need to get from the parameter file so that we can build our object. These will include any pointers to other classes (e.g. ``cosmologyFunctions_`` here), and parameters to be read from the parameter file. Note that the names of these should match the name of the corresponding entry in the `implementation type <#the-implementation-type>`_.

.. note::

   The reason for this is that Galacticus will automatically match names in this constructor to those in the type, and automatically generate functions that allow the parameters of the object to be serialized and deserialized in various ways.

.. tip::

   Note that here we define a variable ``timescale``, but in the implementation type we called this variable ``timescale_`` (with a trailing underscore to avoid conflicting with the name of a method). Galacticus knows to look for these trailing underscores and will correctly match up these two variables.

We then begin getting parameter values and other class objects as needed from the parameter file. To get a parameter we include code like this (enclosed inside a ``!![`` ``!!]`` section):

.. code-block:: xml

       <inputParameter>
         <name>exponentRedshift</name>
         <defaultValue>0.0d0</defaultValue>
         <description>The exponent of redshift in the timescale for star formation in the halo scaling timescale model.</description>
         <source>parameters</source>
       </inputParameter>

This ``inputParameter`` directive tells Galacticus to get a value from the parameter file and assign it to a variable. In this case we ask for the parameter ``exponentRedshift`` which will be assigned to the variable of the same name, and we specify (via the ``source`` element) that this parameter should be obtained from the ``parameters`` object that our function was given as an argument. We also give a description of the parameter (which, as usual, goes into the ReadTheDocs documentation). We can also choose to specify a default value - here we chose ``0.0d0``. If there is an obvious or sensible choice for a default, set one - if there isn't, leave this element out. If no default is given and the parameter file doesn't provide a value, Galacticus will write an error message and stop.

To get pointers to objects of other classes we can use the ``objectBuilder`` directive, like this:

.. code-block:: xml

       <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>

We simply tell it what class of object we want, ``cosmologyFunctions``, what variable to assign this to, ``cosmologyFunctions_`` (note the trailing underscore), and from where to get this object, the ``parameters`` object.

Now that we have the parameters and other objects that we need, we need to actually construct our own object. A typical way to do this is to just call the "internal" constructor function, and make it do the work. The advantage to this is that we then don't have to duplicate the logic of setting up our object - we just do it all in the internal constructor function. So, we have:

.. code-block:: fortran

       self=starFormationTimescaleHaloScaling(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_)

where we call the constructor, ``starFormationTimescaleHaloScaling`` (note that when calling a constructor, we just use the name of the implementation itself - the compiler will figure out which constructor function we actually want here), giving it all of the arguments it requires (see below), and assign the result to our ``self`` object - which is now fully initialized.

.. tip::

   If you get the arguments to the internal constructor wrong (out of order, wrong type, wrong rank, etc.) you'll get a not so helpful error message from the compiler. If that happens, try temporarily calling the internal constructor function directly by replacing ``starFormationTimescaleHaloScaling`` in the above with ``haloScalingConstructorInternal`` - you'll probably get a more informative error message.

We now just have to clean up, and return our resulting object. First, the clean up:

.. code-block:: xml

       <inputParametersValidate source="parameters"/>
       <objectDestructor name="cosmologyFunctions_" />
       <objectDestructor name="darkMatterHaloScale_"/>

The first of these lines adds some automatic checking of the parameters provided to this function. If the parameter file included some parameter that was not recognized by this function, Galacticus will emit a warning.

The next two lines destroy the two class objects that we needed to build our own object. It may seem strange to destroy these objects - they're still needed by our own object of course. Galacticus use a `reference counting <https://en.wikipedia.org/wiki/Reference_counting>`_ approach, so it knows that, even though we've told it to destroy these objects, they're in use elsewhere, so it doesn't actually destroy them. It will keep them around until the last reference to them is removed, only then will it destroy the object.

Now we just return our result, and end the function.

.. code-block:: fortran

       return
     end function haloScalingConstructorParameters

Our second constructor is the *internal* constructor - this refers to a constructor that we can use internally to build an object of our implementation - we already did this in our *parameters* constructor above. The internal constructor looks like this:

.. code-block:: fortran

     function haloScalingConstructorInternal(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_) result(self)
       !!{
       Internal constructor for the {\normalfont \ttfamily haloScaling} timescale for star formation class.
       !!}
       implicit none
       type            (starFormationTimescaleHaloScaling)                        :: self
       double precision                                   , intent(in   )         :: timescale           , exponentVelocityVirial, &
            &                                                                        exponentRedshift
       class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
       class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
       !![
       <constructorAssign variables="exponentVelocityVirial, exponentRedshift, *cosmologyFunctions_, *darkMatterHaloScale_"/>
       !!]

       self%lastUniqueID                 =-1_kind_int8
       self%timescaleComputed            =.false.
       self%velocityPrevious             =-1.0d0
       self%velocityFactorPrevious       =-1.0d0
       self%expansionFactorPrevious      =-1.0d0
       self%expansionFactorFactorPrevious=-1.0d0
       ! Compute the normalization of the timescale.
       self%timeScale_            =+timescale
       self%timeScaleNormalization=+timescale                                                &
            &                      /velocityVirialNormalization**self%exponentVelocityVirial
       return
     end function haloScalingConstructorInternal

It begins by defining the function arguments and return value:

.. code-block:: fortran

   function haloScalingConstructorInternal(timescale,exponentVelocityVirial,exponentRedshift,cosmologyFunctions_,darkMatterHaloScale_) result(self)

For this constructor we pass, as arguments, all of the quantities and objects that our object will need. We also specify that the object will be constructed and return as ``self``.

There is then the usual brief description text, after which we define the types of all the arguments:

.. code-block:: fortran

       implicit none
       type            (starFormationTimescaleHaloScaling)                        :: self
       double precision                                   , intent(in   )         :: timescale           , exponentVelocityVirial, &
            &                                                                        exponentRedshift
       class           (cosmologyFunctionsClass          ), intent(in   ), target :: cosmologyFunctions_
       class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_

As usual, we begin with ``implicit none`` to assert that all variables must be declared. The ``self`` object must be of ``type(starFormationTimescaleHaloScaling)`` in any constructor, so we declare that here. Then we give the types of all of the arguments.

A common thing that we want to do in a constructor is to take the arguments given to us, and assign them to the variables inside of our ``self`` object. Galacticus provides a directive, ``constructorAssign``, to do this for us. It takes less code than writing these assignments out one at a time, and it also does some reference counting behind the scenes to ensure that we keep track of any other class objects being used.

.. code-block:: xml

       <constructorAssign variables="exponentVelocityVirial, exponentRedshift, *cosmologyFunctions_, *darkMatterHaloScale_"/>

This tells Galacticus to assign each of the named ``variables`` to the corresponding object in ``self``. For the pointers to other objects, we add a ``*`` prefix - this tells Galacticus to do pointer assignment (instead of copying the object) and handles reference counting.

We can now do any other initialization that may be needed for our object. In our example we have:

.. code-block:: fortran

       self%lastUniqueID                 =-1_kind_int8
       self%timescaleComputed            =.false.
       self%velocityPrevious             =-1.0d0
       self%velocityFactorPrevious       =-1.0d0
       self%expansionFactorPrevious      =-1.0d0
       self%expansionFactorFactorPrevious=-1.0d0

which just initializes the variables used for memoization to some unphysical values (good practice since if we attempt to use them before updating them we will then get an error message [which is better than not getting an error message and instead getting wrong results]), and:

.. code-block:: fortran

       ! Compute the normalization of the timescale.
       self%timeScale_            =+timescale
       self%timeScaleNormalization=+timescale                                                &
            &                      /velocityVirialNormalization**self%exponentVelocityVirial

which precomputes (and stores in the ``self`` object) a quantity that we will need to use every time our object is used to compute a star formation timescale. Since this value will be the same every time, it makes things slightly faster to compute it once, store the result, and re-use that stored result when needed.

We then just return our constructed object and finish the function.

.. code-block:: fortran

       return
     end function haloScalingConstructorInternal

The ``autoHook`` method
^^^^^^^^^^^^^^^^^^^^^^^

The ``autoHook`` method is a special that gets called automatically whenever an object is created or copied. This can be useful to perform actions on a new object. In our example case, the function linked to the ``autoHook`` method looks like this:

.. code-block:: fortran

     subroutine haloScalingAutoHook(self)
       !!{
       Attach to the calculation reset event.
       !!}
       use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
       implicit none
       class(starFormationTimescaleHaloScaling), intent(inout) :: self

       call calculationResetEvent%attach(self,haloScalingCalculationReset,openMPThreadBindingAllLevels,label='starFormationTimescaleHaloScaling')
       return
     end subroutine haloScalingAutoHook

It takes a single argument, ``self``, which must be declared with:

.. code-block:: fortran

       class(starFormationTimescaleHaloScaling), intent(inout) :: self

This is the object that was created/copied. In this example, we're using the ``autoHook`` method to attach our object to the ``calculationReset`` event. We won't go into detail about this event (or events in general) here, but, briefly, the ``calculationReset`` event is called every time Galacticus updates a node in a merger tree. Since our implementation uses memoization to store and re-use results of its calculations it needs to know when those stored results are no longer valid. The ``calculationReset`` event lets it know that. We tell the ``calculationReset`` event to call our function `haloScalingCalculationReset <#other-functions>`_ as needed.

Destructor
^^^^^^^^^^

We created our object with the constructors described above. When our object is no longer needed, we need to destroy it. That is handled by the `destructor <https://en.wikipedia.org/wiki/Destructor_(computer_programming)>`_. This is a function that cleans up anything in the object before the memory associated with it is destroyed. In our example, the destructor looks like this:

.. code-block:: fortran

     subroutine haloScalingDestructor(self)
       !!{
       Destructor for the {\normalfont \ttfamily haloScaling} timescale for star formation class.
       !!}
       use :: Events_Hooks, only : calculationResetEvent
       implicit none
       type(starFormationTimescaleHaloScaling), intent(inout) :: self

       !![
       <objectDestructor name="self%cosmologyFunctions_" />
       <objectDestructor name="self%darkMatterHaloScale_"/>
       !!]
       if (calculationResetEvent%isAttached(self,haloScalingCalculationReset)) call calculationResetEvent%detach(self,haloScalingCalculationReset)
       return
     end subroutine haloScalingDestructor

Destructors have a single ``self`` argument, declared as:

.. code-block:: fortran

       type(starFormationTimescaleHaloScaling), intent(inout) :: self

In our example, there are two things we need to do in our destructor:

#. Detach from the ``calculationReset`` event (that we attached to in via the `autoHook <#the-autohook-method>`_ method);
#. Destroy any objects that we kept pointers to in our object. This is done using the ``objectDestructor`` directive. For example:

.. code-block:: xml

       <objectDestructor name="self%cosmologyFunctions_" />

tells Galacticus that we are done with the object pointed to by ``self%cosmologyFunctions_``. It will update the reference count to this object, and, if it decides the object is no longer in use anywhere, it will destroy it.

Other functions
^^^^^^^^^^^^^^^

In our example implementation we also have a function ``haloScalingCalculationReset`` that will be called by the ``calculationReset`` event (see `above <#the-autohook-method>`_ for details). We won't go into details about this function, but note that it simply resets the state of our implementation's memoization variables, forcing them to be recomputed as needed.

The physics!
^^^^^^^^^^^^

Finally, we get to implementing the actual physics of this class. (Object-oriented programming can be very verbose - we did a lot of work above just to get to the point of writing this simple function that actually does the physics. The advantage of this however is that we then have a lot of flexibility in how we piece together different parts of the model.)

In our example, our function looks like this:

.. code-block:: fortran

     double precision function haloScalingTimescale(self,component) result(timescale)
       !!{
       Returns the timescale (in Gyr) for star formation in the given {\normalfont \ttfamily component} in the halo scaling
       timescale model.
       !!}
       use :: Galacticus_Nodes, only : nodeComponentBasic
       implicit none
       class           (starFormationTimescaleHaloScaling), intent(inout) :: self
       class           (nodeComponent                    ), intent(inout) :: component
       class           (nodeComponentBasic               ), pointer       :: basic
       double precision                                                   :: expansionFactor, velocityVirial

       ! Check if node differs from previous one for which we performed calculations.
       if (component%hostNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(component%hostNode,component%hostNode%uniqueID())
       ! Compute the timescale if necessary.
       if (.not.self%timescaleComputed) then
          ! Get virial velocity and expansion factor.
          basic           => component%hostNode%basic                               (                    )
          velocityVirial  =  self              %darkMatterHaloScale_%velocityVirial (component%hostNode  )
          expansionFactor =  self              %cosmologyFunctions_ %expansionFactor(basic    %time    ())
          ! Compute the velocity factor.
          if (velocityVirial /= self%velocityPrevious) then
              self%velocityPrevious      =velocityVirial
              self%velocityFactorPrevious=velocityVirial**self%exponentVelocityVirial
          end if
          ! Compute the expansion-factor factor.
          if (expansionFactor /= self%expansionFactorPrevious) then
             self%expansionFactorPrevious      =      expansionFactor
             self%expansionFactorFactorPrevious=1.0d0/expansionFactor**self%exponentRedshift
          end if
          ! Computed the timescale.
          self%timescaleStored=+self%timeScaleNormalization        &
               &               *self%velocityFactorPrevious        &
               &               *self%expansionFactorFactorPrevious
          ! Record that the timescale is now computed.
          self%timescaleComputed=.true.
       end if
       ! Return the stored timescale.
       haloScalingTimescale=self%timescaleStored
       return
     end function haloScalingTimescale

We start by opening the function:

.. code-block:: fortran

     double precision function haloScalingTimescale(self,component) result(timescale)

The type of the function (``double precision``) and the names of the arguments *must* match precisely with those defined in our `class <#the-class>`_. The name of the ``result`` can be whatever you want, but typically it makes sense to set this to match the method name, i.e. ``result(timescale)`` in this example. The declarations of the ``self`` argument must be:

.. code-block:: fortran

       class           (starFormationTimescaleHaloScaling), intent(inout) :: self

and other arguments must match their definitions in our `class <#the-class>`_:

.. code-block:: fortran

       class           (nodeComponent                    ), intent(inout) :: component

The remainder of the function does whatever it needs to do to evaluate this particular model of the star formation timescale. We won't discuss this in detail in this tutorial, but a few things are worth noticing:

#. To make use of memoization we do:

   .. code-block:: fortran

      ! Check if node differs from previous one for which we performed calculations.
      if (component%hostNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(component%hostNode,component%hostNode%uniqueID())
      ! Compute the timescale if necessary.
      if (.not.self%timescaleComputed) then

   which first checks if the ``component`` we were passed belongs to the same node for which we've already computed and memoized results. It does this by comparing the ``uniqueID()`` of the node (which, as the name suggests, is guaranteed to be unique to each node in any Galacticus run). It then checks if the timescale has already been computed and, if not, computes it, recording at the end that the memoized results are now initialized:

   .. code-block:: fortran

      ! Record that the timescale is now computed.
      self%timescaleComputed=.true.

#. We make use of the other class objects that we kept pointers to, for example:

   .. code-block:: fortran

      velocityVirial  =  self              %darkMatterHaloScale_%velocityVirial (component%hostNode  )
      expansionFactor =  self              %cosmologyFunctions_ %expansionFactor(basic    %time    ())

   in which we compute the virial velocity of the node to which the ``component`` belongs, and compute the expansion factor of the universe at the time at which this node exists.
#. The final result is returned by assigning it from the memoized value:

   .. code-block:: fortran

      ! Return the stored timescale.
      timescale=self%timescaleStored

Using your new class
~~~~~~~~~~~~~~~~~~~~

You can now use your new class anywhere in Galacticus. For example, some other class might want to know the timescale for star formation. To do so, it would import the class:

.. code-block:: fortran

   use :: Star_Formation_Timescales, only : starFormationTimescaleClass

add a pointer to it inside its own implementation type:

.. code-block:: fortran

   class(starFormationTimescaleClass), pointer :: starFormationTimescale_ => null()

and build an instance of our new class using:

.. code-block:: xml

   <objectBuilder class="starFormationTimescale" name="starFormationTimescale_" source="parameters"/>

in just the same way as our example implementation above made use of the ``cosmologyFunctions`` class.

Then, when you compile Galacticus it will automatically discover your new class, build it, and link it in to the resulting ``Galacticus.exe`` executable. It will then recognize your new class, and its parameters, in the parameter file, e.g.:

.. code-block:: xml

   <starFormationTimescale value="haloScaling">
     <exponentRedshift       value="0.0"/>
     <exponentVelocityVirial value="0.0"/>
     <timescale              value="1.0"/>
   </starFormationTimescale>

Recursive classes
~~~~~~~~~~~~~~~~~~

Occasionally a class composites a member of its *own* class---either directly, or via another, mutually-compositing class. If no such member is provided explicitly in the parameter file, the build searches up the parameter tree, re-discovers the object it is *currently* constructing, and tries to build it again---an unbounded recursion. Galacticus detects this and aborts with an informative error (see `issue #397 <https://github.com/galacticusorg/galacticus/issues/397>`_), listing the build stack and telling you to provide the required object explicitly.

For a small number of classes, however, such a construction cycle is *legitimate and bounded*: the physics guarantees that re-entry terminates. (Examples include :galacticus-class:`virialDensityContrastPercolation`, whose percolation objects composite a :galacticus-class:`darkMatterHaloScaleVirialDensityContrastDefinition` that resolves back to the percolation object.) To allow such a class to participate in a bounded cycle, add ``recursive="yes"`` to its ``functionClass`` directive:

.. code-block:: fortran

   !![
   <starFormationTimescale name="starFormationTimescaleMyRecursive" recursive="yes">
    <description>A star formation timescale class that composites a bounded member of its own class.</description>
   </starFormationTimescale>
   !!]

That single attribute is all that is required. When the build re-enters the node currently under construction, the factory returns a generated *shim*---a lightweight ``<name>Recursive`` object that holds only a weak pointer back to the real object under construction and forwards every method call to it. All of the supporting machinery is generated automatically: the shim type and its method forwarders, its ``deepCopy``/``stateStore``/``descriptor`` behaviour, and the weak (uncounted) back-reference that keeps reference counting sound so no memory leak results. You do **not** need to write any ``recursiveSelf`` pointers, per-method guards, or custom ``deepCopy`` code---this boilerplate was removed under `issue #695 <https://github.com/galacticusorg/galacticus/issues/695>`_.

.. warning::

   Only mark a class ``recursive="yes"`` if re-entry into it is genuinely bounded. The opt-in is deliberately *not* automatic: for a class whose cycle does *not* terminate (e.g. a filter that wraps itself), automatic shimming would convert the clean construction-time error above into an infinite recursion at *run* time, with no diagnostic. Leaving such a class un-opted-in preserves the informative construction-time abort.
