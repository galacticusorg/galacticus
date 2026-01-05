!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

  use :: Cosmological_Density_Field, only : cosmologicalMassVariance   , cosmologicalMassVarianceClass
  use :: Cosmology_Parameters      , only : cosmologyParameters        , cosmologyParametersClass
  use :: Power_Spectra_Primordial  , only : powerSpectrumPrimordial    , powerSpectrumPrimordialClass
  use :: Stateful_Types            , only : statefulDouble             , statefulInteger              , statefulLogical
  use :: Transfer_Functions        , only : transferFunction           , transferFunctionClass
  use :: Merger_Tree_Data_Structure, only : enumerationPropertyTypeType, enumerationMetaDataTypeType  , enumerationMergerTreeFormatType

  type :: mergerTreeMetaData
     !!{
     Type used for metadata in merger tree files.
     !!}
     type(enumerationMetaDataTypeType) :: type
     type(varying_string             ) :: name, content
  end type mergerTreeMetaData

  type :: propertyColumn
     !!{
     Type used to specify which property to read from which column,
     !!}
     type   (enumerationPropertyTypeType) :: property
     integer                              :: column
     type   (statefulDouble             ) :: conversionFactor
  end type propertyColumn

  !![
  <task name="taskMergerTreeFileBuilder">
   <description>
    This task will build a merger tree file in the format described
    \href{https://github.com/galacticusorg/galacticus/wiki/Merger-Tree-File-Format}{here} from merger tree descriptions in other
    formats, such as ASCII output from an SQL database. An example of how the builder can be used can be found in this
    \href{https://github.com/galacticusorg/galacticus/wiki/Tutorial:-Building-merger-trees-from-ASCII-files}{tutorial}.
    
    The builder is flexible, and therefore requires many parameters to control how it processes input files, all of which are
    described below.
    
    The builder reads from an ASCII file containing one halo per line. The {\normalfont \ttfamily [property]} sub-parameters allow
    specification of which properties are present in the file, and in which column. The builder can optionally also read a file of
    associated particle data---this can be used to assign positions to orphan halos if desired.
    
    The merger tree file builder can currently export in one of two formats:
    \begin{description}
    \item [{\normalfont \ttfamily galacticus}] merger trees are exported in \glc's native format described in detail
      \href{https://github.com/galacticusorg/galacticus/wiki/Merger-Tree-File-Format}{here};
    \item [{\normalfont \ttfamily irate}] merger trees are exported in the \href{https://irate-format.readthedocs.io/en/stable/formatspec.html}{\normalfont
        \ttfamily IRATE} format.
    \end{description}
    
    Properties to read from the file are specified through multiple {\normalfont \ttfamily property} sub-parameter sections, which take the form:
    \begin{verbatim}
     &lt;property>
      &lt;name             value="propertyName"/>
      &lt;column           value="columnNumber"/>
      &lt;conversionFactor value="1.0e0"       />
     &lt;/property>
    \end{verbatim}
    where {\normalfont \ttfamily [name]} is the property name (see below), {\normalfont \ttfamily [column]} is the column number
    (starting from 1) from which to read the property, and the optional {\normalfont \ttfamily [conversionFactor]} specifies an
    additional factor by which the property should be multiplied to place it into the correct internal units for \glc\footnote{The
      units for masses, lengths, and velocities in the input file are specified in their own parameter sub-sections. Conversion from
      these units to \glc's internal units is performed automatically. However, sometimes the input data may have inconsistent units
      between columns (e.g. positions in units of Mpc, but scale radii in units of kpc). In such cases this additional conversion
      factor can be applied to bring all quantities into a consistent unit system.}.
    
    Recognized property names are
    \begin{description}
    \item [{\normalfont \ttfamily treeIndex}] A unique ID number for the tree to which this node belongs;
    \item [{\normalfont \ttfamily nodeIndex}] An ID (unique within the tree) for this node;
    \item [{\normalfont \ttfamily descendantIndex}] The ID of the node's descendant node;
    \item [{\normalfont \ttfamily hostIndex}] The ID of the larger halo in which this node is hosted (equal to the node's own ID if
    the node is self-hosting);
    \item [{\normalfont \ttfamily redshift}] The redshift of the node;
    \item [{\normalfont \ttfamily nodeMass}] The mass of the node;
    \item [{\normalfont \ttfamily particleCount}] The number of particles in the node;
    \item [{\normalfont \ttfamily positionX}] The $x$-position of the node (if present, both $y$ and $z$ components must also be
    present);
    \item [{\normalfont \ttfamily positionY}] The $y$-position of the node (if present, both $x$ and $z$ components must also be
    present);
    \item [{\normalfont \ttfamily positionZ}] The $z$-position of the node (if present, both $x$ and $y$ components must also be
    present);
    \item [{\normalfont \ttfamily velocityX}] The $x$-velocity of the node (if present, both $y$ and $z$ components must also be
    present);
    \item [{\normalfont \ttfamily velocityY}] The $y$-velocity of the node (if present, both $x$ and $z$ components must also be
    present);
    \item [{\normalfont \ttfamily velocityZ}] The $z$-velocity of the node (if present, both $x$ and $y$ components must also be
    present);
    \item [{\normalfont \ttfamily spinX}] The $x$ component of the node's spin parameter (if present, both $y$ and $z$ components must
    also be present; cannot be present if spin magnitude is given);
     \item [{\normalfont \ttfamily spinY}] The $y$ component of the node's spin parameter (if present, both $x$ and $z$ components
    must also be present; cannot be present if spin magnitude is given);
     \item [{\normalfont \ttfamily spinZ}] The $z$ component of the node's spin parameter (if present, both $x$ and $y$ components
    must also be present; cannot be present if spin magnitude is given);
     \item [{\normalfont \ttfamily spin}] The magnitude of the node's spin parameter (cannot be present if spin vector components are
    given);
     \item [{\normalfont \ttfamily angularMomentumX}] The $x$-component of the node's angular momentum (if present, both $y$ and $z$
    components must also be present; cannot be present if angular momentum magnitude is given);
     \item [{\normalfont \ttfamily angularMomentumY}] The $y$-component of the node's angular momentum (if present, both $x$ and $z$
    components must also be present; cannot be present if angular momentum magnitude is given);
     \item [{\normalfont \ttfamily angularMomentumZ}] The $z$-component of the node's angular momentum (if present, both $x$ and $y$
    components must also be present; cannot be present if angular momentum magnitude is given);
     \item [{\normalfont \ttfamily angularMomentum}] The magnitude of the node's angular momentum (cannot be present if angular
    momentum vector components are given);
     \item [{\normalfont \ttfamily halfMassRadius}] The half-mass radius of the node;
     \item [{\normalfont \ttfamily mostBoundParticleIndex}] The index of the most bound particle in this node.
    \end{description}
    Not all properties must be specified---any required properties that are not specified will result in an error. Likewise, some
    properties, if present, require that other properties also be present. For example, if any of the position properties is given
    then all three positions are required.
    
    Properties of particles to read from the (optional) particle data file are specified through multiple {\normalfont \ttfamily
      particleProperty} sub-parameter sections, which take the form:
    \begin{verbatim}
     &lt;particleProperty>
      &lt;name   value="propertyName"/>
      &lt;column value="columnNumber"/>
     &lt;/particleProperty>
    \end{verbatim}
    where {\normalfont \ttfamily [name]} is the particle property name (see below), {\normalfont \ttfamily [column]} is the column
    number (starting from 1) from which to read the particle property.
    
    Recognized particle property names are
    \begin{description}
     \item [{\normalfont \ttfamily particleIndex}] A unique ID for the particle;
     \item [{\normalfont \ttfamily redshift}] The redshift of the particle;
     \item [{\normalfont \ttfamily nodeMass}] The mass of the particle;
     \item [{\normalfont \ttfamily particleCount}] The number of particles in the particle;
     \item [{\normalfont \ttfamily positionX}] The $x$-position of the particle (if present, both $y$ and $z$ components must also be present);
     \item [{\normalfont \ttfamily positionY}] The $y$-position of the particle (if present, both $x$ and $z$ components must also be present);
     \item [{\normalfont \ttfamily positionZ}] The $z$-position of the particle (if present, both $x$ and $y$ components must also be present);
     \item [{\normalfont \ttfamily velocityX}] The $x$-velocity of the particle (if present, both $y$ and $z$ components must also be present);
     \item [{\normalfont \ttfamily velocityY}] The $y$-velocity of the particle (if present, both $x$ and $z$ components must also be present);
     \item [{\normalfont \ttfamily velocityZ}] The $z$-velocity of the particle (if present, both $x$ and $y$ components must also be present).
    \end{description}
    
    The units used in the files are specified via the {\normalfont \ttfamily unitsMass}, {\normalfont \ttfamily unitsLength}, and
    {\normalfont \ttfamily unitsVelocity} sub-parameter sections. These have the following form:
    \begin{verbatim}
     &lt;unitsMass>
      &lt;name                value="Solar masses/h"/>
      &lt;unitsInSI           value="1.99e30"       />
      &lt;hubbleExponent      value="-1"            />
      &lt;scaleFactorExponent value=" 0"            />
     &lt;/unitsMass>
    \end{verbatim}
    where {\normalfont \ttfamily [name]} is a human-readable name for the units, {\normalfont \ttfamily [unitsInSI]} gives the units
    in the SI system, {\normalfont \ttfamily [hubbleExponent]} specifies the power to which $h$ appears in the units and {\normalfont
      \ttfamily [scaleFactorExponent]} specifies the number of powers of the expansion factor by which the quantity should be
    multiplied to place it into physical units.
    
    Finally, arbitrary metadata can be added to the file (which can be useful to record, for example, the origin of the data, or
    details of the simulation, or halo finder used). Metadata is specified via {\normalfont \ttfamily metaData} sub-parameter
    sections. These have the following form:
    \begin{verbatim}
     &lt;metaData>
      &lt;name    value="simulationCode"/>
      &lt;content value="Gadget2"       />
      &lt;type    value="simulation"    />
     &lt;/metaData>
    \end{verbatim}
    where {\normalfont \ttfamily [name]} is a name for this metadatum, {\normalfont \ttfamily [content]} is the value for the
    metadatum (integer, floating point, and text content is allowed), and {\normalfont \ttfamily [type]} specifies the type of
    metadatum, and must be one of:
    \begin{description}
     \item [{\normalfont \ttfamily generic}] Add to the generic {\normalfont \ttfamily metaData} group;
     \item [{\normalfont \ttfamily cosmology}] Add to the {\normalfont \ttfamily cosmology} group;
     \item [{\normalfont \ttfamily simulation}] Add to the {\normalfont \ttfamily simulation} group;
     \item [{\normalfont \ttfamily groupFinder}] Add to the {\normalfont \ttfamily groupFinder} group;
     \item [{\normalfont \ttfamily treeBuilder}] Add to the {\normalfont \ttfamily treeBuilder} group;
     \item [{\normalfont \ttfamily provenance}] Add to the {\normalfont \ttfamily provenance} group.
    \end{description}
   </description>
    <descriptorSpecial>descriptorSpecial</descriptorSpecial>
  </task>
  !!]
  type, extends(taskClass) :: taskMergerTreeFileBuilder
     !!{
     Implementation of a task which computes and outputs the halo mass function and related quantities.
     !!}
     private
     class           (cosmologyParametersClass       ), pointer                   :: cosmologyParameters_       => null()
     class           (cosmologicalMassVarianceClass  ), pointer                   :: cosmologicalMassVariance_  => null()
     class           (powerSpectrumPrimordialClass   ), pointer                   :: powerSpectrumPrimordial_   => null()
     class           (transferFunctionClass          ), pointer                   :: transferFunction_          => null()
     type            (propertyColumn                 ), allocatable, dimension(:) :: properties                          , particleProperties
     type            (mergerTreeMetaData             ), allocatable, dimension(:) :: metaData
     type            (varying_string                 )                            :: outputFileName                      , inputFileName                   , &
          &                                                                          particlesFileName
     type            (enumerationMergerTreeFormatType)                            :: outputFormat
     double precision                                                             :: unitsMassInSI                       , unitsLengthInSI                 , &
          &                                                                          unitsVelocityInSI
     integer                                                                      :: unitsMassHubbleExponent             , unitsMassScaleFactorExponent    , &
          &                                                                          unitsLengthHubbleExponent           , unitsLengthScaleFactorExponent  , &
          &                                                                          unitsVelocityHubbleExponent         , unitsVelocityScaleFactorExponent
     type            (varying_string                 )                            :: unitsMassName                       , unitsLengthName                 , &
          &                                                                          unitsVelocityName
     type            (statefulInteger                )                            :: dummyHostId
     type            (statefulDouble                 )                            :: massParticle
     type            (statefulLogical                )                            :: haloMassesIncludeSubhalos           , includesHubbleFlow              , &
          &                                                                          positionsArePeriodic
     logical                                                                      :: traceParticles                      , columnHeaders
     character       (len=1                          )                            :: columnSeparator
   contains
     !![
     <methods>
       <method method="descriptorSpecial" description="Handle adding special parameters to the descriptor."/>
     </methods>
     !!]
     final     ::                       mergerTreeFileBuilderDestructor
     procedure :: perform            => mergerTreeFileBuilderPerform
     procedure :: requiresOutputFile => mergerTreeFileBuilderRequiresOutputFile
     procedure :: descriptorSpecial  => mergerTreeFileBuilderDescriptorSpecial
  end type taskMergerTreeFileBuilder

  interface taskMergerTreeFileBuilder
     !!{
     Constructors for the \refClass{taskMergerTreeFileBuilder} task.
     !!}
     module procedure mergerTreeFileBuilderConstructorParameters
     module procedure mergerTreeFileBuilderConstructorInternal
  end interface taskMergerTreeFileBuilder

contains

  function mergerTreeFileBuilderConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskMergerTreeFileBuilder} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters                , only : inputParameter                   , inputParameters
    use :: Merger_Tree_Data_Structure      , only : enumerationMergerTreeFormatEncode, enumerationMetaDataTypeEncode, enumerationPropertyTypeEncode, unitsLength, &
          &                                         unitsMass                        , unitsVelocity
    use :: Numerical_Constants_Astronomical, only : massSolar                        , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    type            (taskMergerTreeFileBuilder    )                              :: self
    type            (inputParameters              ), intent(inout)               :: parameters
    class           (cosmologyParametersClass     ), pointer                     :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass), pointer                     :: cosmologicalMassVariance_
    class           (powerSpectrumPrimordialClass ), pointer                     :: powerSpectrumPrimordial_
    class           (transferFunctionClass        ), pointer                     :: transferFunction_
    type            (propertyColumn               ), allocatable  , dimension(:) :: properties                 , particleProperties
    type            (mergerTreeMetaData           ), allocatable  , dimension(:) :: metaData
    type            (inputParameters              ), allocatable                 :: subParameters
    type            (varying_string               )                              :: inputFileName              , outputFileName                  , &
         &                                                                          outputFormat               , name                            , &
         &                                                                          particlesFileName
    integer                                                                      :: i                          , column                          , &
         &                                                                          countProperties            , countMetaData
    type            (statefulInteger              )                              :: dummyHostId
    type            (statefulDouble               )                              :: massParticle
    type            (statefulLogical              )                              :: haloMassesIncludeSubhalos  , includesHubbleFlow              , &
         &                                                                          positionsArePeriodic
    double precision                                                             :: unitsMassInSI              , unitsLengthInSI                 , &
         &                                                                          unitsVelocityInSI
    integer                                                                      :: unitsMassHubbleExponent    , unitsMassScaleFactorExponent    , &
         &                                                                          unitsLengthHubbleExponent  , unitsLengthScaleFactorExponent  , &
         &                                                                          unitsVelocityHubbleExponent, unitsVelocityScaleFactorExponent
    type            (varying_string               )                              :: unitsMassName              , unitsLengthName                 , &
         &                                                                          unitsVelocityName          , metaDataType                    , &
         &                                                                          metaDataContent
    logical                                                                      :: columnHeaders
    type            (varying_string               )                              :: columnSeparator

    !![
    <inputParameter>
      <name>inputFileName</name>
      <description>The name of the file from which to read merger tree data.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('particlesFileName')) then
       !![
       <inputParameter>
         <name>particlesFileName</name>
         <description>The name of the file from which to read particle data.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    else
       particlesFileName=''
    end if
    !![
    <inputParameter>
      <name>outputFileName</name>
      <description>The name of the file to which to write merger tree data.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputFormat</name>
      <description>The format to use for merger tree output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>columnHeaders</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, the file is assumed to contain a single line of column headers, which will be skipped.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>columnSeparator</name>
      <defaultValue>var_str(',')</defaultValue>
      <description>The separator for columns.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (columnSeparator == "space") columnSeparator=" "
    massParticle%isSet=parameters%isPresent('massParticle')
    if (massParticle%isSet) then
       !![
       <inputParameter>
         <name>massParticle</name>
         <description>The mass of the simulation particle.</description>
         <source>parameters</source>
         <variable>massParticle%value</variable>
       </inputParameter>
       !!]
    end if
    dummyHostId%isSet=parameters%isPresent('dummyHostId')
    if (dummyHostId%isSet) then
       !![
       <inputParameter>
         <name>dummyHostId</name>
         <description>If present, specifies the dummy host ID for self-hosting halos. Otherwise, self-hosting halos have {\normalfont \ttfamily hostIndex == nodeIndex}.</description>
         <source>parameters</source>
         <variable>dummyHostId%value</variable>
       </inputParameter>
       !!]
    end if
    haloMassesIncludeSubhalos%isSet=parameters%isPresent('haloMassesIncludeSubhalos')
    if (haloMassesIncludeSubhalos%isSet) then
       !![
       <inputParameter>
         <name>haloMassesIncludeSubhalos</name>
         <description>Specifies whether or not halo masses include the masses of their subhalos.</description>
         <source>parameters</source>
         <variable>haloMassesIncludeSubhalos%value</variable>
       </inputParameter>
       !!]
    end if
    includesHubbleFlow%isSet=parameters%isPresent('includesHubbleFlow')
    if (includesHubbleFlow%isSet) then
       !![
       <inputParameter>
         <name>includesHubbleFlow</name>
         <description>Specifies whether or not Hubble flow is included in velocities.</description>
         <source>parameters</source>
         <variable>includesHubbleFlow%value</variable>
       </inputParameter>
       !!]
    end if
    positionsArePeriodic%isSet=parameters%isPresent('positionsArePeriodic')
    if (positionsArePeriodic%isSet) then
       !![
       <inputParameter>
         <name>positionsArePeriodic</name>
         <description>Specifies whether or not positions are periodic.</description>
         <source>parameters</source>
         <variable>positionsArePeriodic%value</variable>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="powerSpectrumPrimordial"  name="powerSpectrumPrimordial_"  source="parameters"/>
    <objectBuilder class="transferFunction"         name="transferFunction_"         source="parameters"/>
    !!]
    countProperties=parameters%copiesCount('property',requireValue=.false.)
    allocate(properties(countProperties))
    do i=1,countProperties
       allocate(subParameters)
       subParameters=parameters%subParameters('property',requireValue=.false.,copyInstance=i)
       !![
       <inputParameter>
         <name>name</name>
         <description>The name of the property to read.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>column</name>
         <description>The column from which to read the property.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       properties(i)%property              =enumerationPropertyTypeEncode(char(name),includesPrefix=.false.)
       properties(i)%column                =column
       properties(i)%conversionFactor%isSet=subParameters%isPresent('conversionFactor')
       if (properties(i)%conversionFactor%isSet) then
          !![
          <inputParameter>
            <name>conversionFactor</name>
            <variable>properties(i)%conversionFactor%value</variable>
            <description>An additional conversion factor to apply to the property to get it into the correct units.</description>
	    <source>subParameters</source>
	    <type>real</type>
	    <cardinality>1</cardinality>
          </inputParameter>
          !!]
       end if
       deallocate(subParameters)
    end do
    countProperties=parameters%copiesCount('particleProperty',requireValue=.false.,zeroIfNotPresent=.true.)
    allocate(particleProperties(countProperties))
    do i=1,countProperties
       allocate(subParameters)
       subParameters=parameters%subParameters('particleProperty',requireValue=.false.,copyInstance=i)
       !![
       <inputParameter>
         <name>name</name>
         <description>The name of the particle property to read.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>column</name>
         <description>The column from which to read the particle property.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       deallocate(subParameters)
       particleProperties(i)%property=enumerationPropertyTypeEncode(char(name),includesPrefix=.false.)
       particleProperties(i)%column  =column
    end do
    countMetaData=parameters%copiesCount('metaData',requireValue=.false.,zeroIfNotPresent=.true.)
    allocate(metaData(countMetaData))
    do i=1,countMetaData
       allocate(subParameters)
       subParameters=parameters%subParameters('metaData',requireValue=.false.,copyInstance=i)
       !![
       <inputParameter>
         <name>name</name>
         <description>The name of the metadata.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>content</name>
         <variable>metaDataContent</variable>
         <description>The value of the metadata.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>type</name>
         <variable>metaDataType</variable>
         <description>The metadata type.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       deallocate(subParameters)
       metaData(i)%name   =name
       metaData(i)%content=metaDataContent
       metaData(i)%type   =enumerationMetaDataTypeEncode(char(metaDataType),includesPrefix=.false.)
    end do
    unitsMassInSI               =massSolar
    unitsMassHubbleExponent     =0
    unitsMassScaleFactorExponent=0
    unitsMassName               =''
    if (parameters%isPresent('unitsMass',requireValue=.false.)) then
       allocate(subParameters)
       subParameters=parameters%subParameters('unitsMass',requireValue=.false.)
       !![
       <inputParameter>
         <name>unitsInSI</name>
         <variable>unitsMassInSI</variable>
         <description>The mass unit in the SI system.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>hubbleExponent</name>
         <variable>unitsMassHubbleExponent</variable>
         <description>The exponent of the ``little-$h$'' Hubble parameter needed to convert the masses to little-$h$-free units.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>scaleFactorExponent</name>
         <variable>unitsMassScaleFactorExponent</variable>
         <description>The exponent of the cosmological scale factor needed to convert the masses to physical units.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>name</name>
         <variable>unitsMassName</variable>
         <description>A human-readable name for the units of mass.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       deallocate(subParameters)
    end if
    unitsLengthInSI               =megaParsec
    unitsLengthHubbleExponent     =0
    unitsLengthScaleFactorExponent=0
    unitsLengthName               =''
    if (parameters%isPresent('unitsLength',requireValue=.false.)) then
       allocate(subParameters)
       subParameters=parameters%subParameters('unitsLength',requireValue=.false.)
       !![
       <inputParameter>
         <name>unitsInSI</name>
         <variable>unitsLengthInSI</variable>
         <description>The length unit in the SI system.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>hubbleExponent</name>
         <variable>unitsLengthHubbleExponent</variable>
         <description>The exponent of the ``little-$h$'' Hubble parameter needed to convert the lengths to little-$h$-free units.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>scaleFactorExponent</name>
         <variable>unitsLengthScaleFactorExponent</variable>
         <description>The exponent of the cosmological scale factor needed to convert the lengths to physical units.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>name</name>
         <variable>unitsLengthName</variable>
         <description>A human-readable name for the units of length.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       deallocate(subParameters)
    end if
    unitsVelocityInSI               =kilo
    unitsVelocityHubbleExponent     =0
    unitsVelocityScaleFactorExponent=0
    unitsVelocityName               =''
    if (parameters%isPresent('unitsVelocity',requireValue=.false.)) then
       allocate(subParameters)
       subParameters=parameters%subParameters('unitsVelocity',requireValue=.false.)
       !![
       <inputParameter>
         <name>unitsInSI</name>
         <variable>unitsVelocityInSI</variable>
         <description>The velocity unit in the SI system.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>hubbleExponent</name>
         <variable>unitsVelocityHubbleExponent</variable>
         <description>The exponent of the ``little-$h$'' Hubble parameter needed to convert the velocities to little-$h$-free units.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>scaleFactorExponent</name>
         <variable>unitsVelocityScaleFactorExponent</variable>
         <description>The exponent of the cosmological scale factor needed to convert the velocities to physical units.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>name</name>
         <variable>unitsVelocityName</variable>
         <description>A human-readable name for the units of velocity.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
       deallocate(subParameters)
    end if
    self=taskMergerTreeFileBuilder(inputFileName,particlesFileName,outputFileName,enumerationMergerTreeFormatEncode(char(outputFormat),includesPrefix=.false.),columnHeaders,char(columnSeparator),properties,particleProperties,metaData,massParticle,dummyHostId,haloMassesIncludeSubhalos,includesHubbleFlow,positionsArePeriodic,unitsMassInSI,unitsMassHubbleExponent,unitsMassScaleFactorExponent,unitsMassName,unitsLengthInSI,unitsLengthHubbleExponent,unitsLengthScaleFactorExponent,unitsLengthName,unitsVelocityInSI,unitsVelocityHubbleExponent,unitsVelocityScaleFactorExponent,unitsVelocityName,cosmologyParameters_,cosmologicalMassVariance_,powerSpectrumPrimordial_,transferFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="powerSpectrumPrimordial_" />
    <objectDestructor name="transferFunction_"        />
    !!]
    return
  end function mergerTreeFileBuilderConstructorParameters

  function mergerTreeFileBuilderConstructorInternal(inputFileName,particlesFileName,outputFileName,outputFormat,columnHeaders,columnSeparator,properties,particleProperties,metaData,massParticle,dummyHostId,haloMassesIncludeSubhalos,includesHubbleFlow,positionsArePeriodic,unitsMassInSI,unitsMassHubbleExponent,unitsMassScaleFactorExponent,unitsMassName,unitsLengthInSI,unitsLengthHubbleExponent,unitsLengthScaleFactorExponent,unitsLengthName,unitsVelocityInSI,unitsVelocityHubbleExponent,unitsVelocityScaleFactorExponent,unitsVelocityName,cosmologyParameters_,cosmologicalMassVariance_,powerSpectrumPrimordial_,transferFunction_) result(self)
    !!{
    Constructor for the \refClass{taskMergerTreeFileBuilder} task class which takes a parameter set as input.
    !!}
    implicit none
    type            (taskMergerTreeFileBuilder      )                              :: self
    class           (cosmologyParametersClass       ), intent(in   ), target       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass  ), intent(in   ), target       :: cosmologicalMassVariance_
    class           (powerSpectrumPrimordialClass   ), intent(in   ), target       :: powerSpectrumPrimordial_
    class           (transferFunctionClass          ), intent(in   ), target       :: transferFunction_
    type            (varying_string                 ), intent(in   )               :: inputFileName              , outputFileName                  , &
         &                                                                            particlesFileName
    type            (statefulInteger                ), intent(in   )               :: dummyHostId
    type            (statefulDouble                 ), intent(in   )               :: massParticle
    type            (statefulLogical                ), intent(in   )               :: haloMassesIncludeSubhalos  , includesHubbleFlow              , &
         &                                                                            positionsArePeriodic
    type            (enumerationMergerTreeFormatType), intent(in   )               :: outputFormat
    type            (propertyColumn                 ), intent(in   ), dimension(:) :: properties                 , particleProperties
    double precision                                 , intent(in   )               :: unitsMassInSI              , unitsLengthInSI                 , &
         &                                                                            unitsVelocityInSI
    type            (mergerTreeMetaData             ), intent(in   ), dimension(:) :: metaData
    integer                                          , intent(in   )               :: unitsMassHubbleExponent    , unitsMassScaleFactorExponent    , &
         &                                                                            unitsLengthHubbleExponent  , unitsLengthScaleFactorExponent  , &
         &                                                                            unitsVelocityHubbleExponent, unitsVelocityScaleFactorExponent
    type            (varying_string                 ), intent(in   )               :: unitsMassName              , unitsLengthName                 , &
         &                                                                            unitsVelocityName
    logical                                          , intent(in   )               :: columnHeaders
    character       (len=1                          ), intent(in   )               :: columnSeparator
    !![
    <constructorAssign variables="inputFileName, particlesFileName, outputFileName, outputFormat, columnHeaders, columnSeparator, properties, particleProperties, metaData, massParticle, dummyHostId, haloMassesIncludeSubhalos, includesHubbleFlow, positionsArePeriodic, unitsMassInSI, unitsMassHubbleExponent, unitsMassScaleFactorExponent, unitsMassName, unitsLengthInSI, unitsLengthHubbleExponent, unitsLengthScaleFactorExponent, unitsLengthName, unitsVelocityInSI, unitsVelocityHubbleExponent, unitsVelocityScaleFactorExponent, unitsVelocityName, *cosmologyParameters_, *cosmologicalMassVariance_, *powerSpectrumPrimordial_, *transferFunction_"/>
    !!]

    self%traceParticles=particlesFileName /= ''
    return
  end function mergerTreeFileBuilderConstructorInternal

  subroutine mergerTreeFileBuilderDestructor(self)
    !!{
    Destructor for the \refClass{taskMergerTreeFileBuilder} task class.
    !!}
    implicit none
    type(taskMergerTreeFileBuilder), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%powerSpectrumPrimordial_" />
    <objectDestructor name="self%transferFunction_"        />
    !!]
    return
  end subroutine mergerTreeFileBuilderDestructor

  subroutine mergerTreeFileBuilderPerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use :: Cosmology_Parameters      , only : hubbleUnitsLittleH
    use :: Dates_and_Times           , only : Formatted_Date_and_Time
    use :: Display                   , only : displayIndent          , displayUnindent
    use :: Error          , only : errorStatusSuccess
    use :: HDF5                      , only : hsize_t
    use :: Input_Parameters          , only : inputParameters
    use :: Merger_Tree_Data_Structure, only : mergerTreeData         , metaDataTypeCosmology, metaDataTypeProvenance, metaDataTypeSimulation, &
          &                                   unitsLength            , unitsMass            , unitsVelocity
    implicit none
    class           (taskMergerTreeFileBuilder), intent(inout), target   :: self
    integer                                    , intent(  out), optional :: status
    integer                                                              :: hdfChunkSize           =1024, hdfCompressionLevel       =9
    type            (mergerTreeData           )                          :: mergerTrees
    type            (inputParameters          )                          :: descriptorPowerSpectrum     , descriptorTransferFunction
    integer                                                              :: i                           , ioStatus                    , &
         &                                                                  metaDataValueInteger
    double precision                                                     :: metaDataValueFloat
    character       (len=1024                 )                          :: metaDataText

    call displayIndent('Begin task: merger tree file builder')
    ! Initialize the data structure.
    call mergerTrees%reset()
    ! Set columns to read.
    do i=1,size(self%properties)
       call mergerTrees%setPropertyColumn(self%properties(i)%property,self%properties(i)%column)
       ! Specify properties that need conversion due to inconsistent units.
       if (self%properties(i)%conversionFactor%isSet) call mergerTrees%setConversionFactor(self%properties(i)%property,self%properties(i)%conversionFactor%value)
    end do
    ! Dummy host ID.
    if (self%dummyHostId              %isSet) &
         & call mergerTrees%setDummyHostId          (self%dummyHostId              %value)
    ! Set subhalo mass definition.
    if (self%haloMassesIncludeSubhalos%isSet) &
         & call mergerTrees%setIncludesSubhaloMasses(self%haloMassesIncludeSubhalos%value)
    ! Set status of Hubble flow inclusion in velocities.
    if (self%includesHubbleFlow       %isSet) &
         & call mergerTrees%setIncludesHubbleFlow   (self%includesHubbleFlow       %value)
    ! Set position periodicity.
    if (self%positionsArePeriodic     %isSet) &
         & call mergerTrees%setPositionsArePeriodic (self%positionsArePeriodic     %value)
    ! Particle mass.
    if (self%massParticle             %isSet) &
         & call mergerTrees%setParticleMass         (self%massParticle             %value)
    ! Read in the data.
    call mergerTrees%readASCII(char(self%inputFileName),columnHeaders=self%columnHeaders,commentCharacter="#",separator=self%columnSeparator)
    ! Specify that we do not want to create individual merger tree reference datasets.
    call mergerTrees%makeReferences          (.false.)
    ! Specify that trees are self-contained (i.e. nodes never move from one tree to another).
    call mergerTrees%setSelfContained        (.true. )
    call mergerTrees%setUnits(unitsMass    ,unitsInSI=self%unitsMassInSI    ,hubbleExponent=self%unitsMassHubbleExponent    ,scaleFactorExponent=self%unitsMassScaleFactorExponent    ,name=char(self%unitsMassName    ))
    call mergerTrees%setUnits(unitsLength  ,unitsInSI=self%unitsLengthInSI  ,hubbleExponent=self%unitsLengthHubbleExponent  ,scaleFactorExponent=self%unitsLengthScaleFactorExponent  ,name=char(self%unitsLengthName  ))
    call mergerTrees%setUnits(unitsVelocity,unitsInSI=self%unitsVelocityInSI,hubbleExponent=self%unitsVelocityHubbleExponent,scaleFactorExponent=self%unitsVelocityScaleFactorExponent,name=char(self%unitsVelocityName))
    ! Get descriptors.
    descriptorPowerSpectrum   =inputParameters()
    descriptorTransferFunction=inputParameters()
    call self%powerSpectrumPrimordial_%descriptor(descriptorPowerSpectrum   ,includeClass=.true.)
    call self%transferFunction_       %descriptor(descriptorTransferFunction,includeClass=.true.)
    ! Set metadata.
    call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaMatter'       ,     self%cosmologyParameters_      %OmegaMatter      (                  ) )
    call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaBaryon'       ,     self%cosmologyParameters_      %OmegaBaryon      (                  ) )
    call mergerTrees%addMetadata(metaDataTypeCosmology ,'OmegaLambda'       ,     self%cosmologyParameters_      %OmegaDarkEnergy  (                  ) )
    call mergerTrees%addMetadata(metaDataTypeCosmology ,'HubbleParam'       ,     self%cosmologyParameters_      %HubbleConstant   (hubbleUnitsLittleH) )
    call mergerTrees%addMetadata(metaDataTypeCosmology ,'sigma_8'           ,     self%cosmologicalMassVariance_ %sigma8           (                  ) )
    call mergerTrees%addMetadata(metaDataTypeCosmology ,'powerSpectrumIndex',char(     descriptorPowerSpectrum   %serializeToString(                  )))
    call mergerTrees%addMetadata(metaDataTypeCosmology ,'transferFunction'  ,char(     descriptorTransferFunction%serializeToString(                  )))
    call mergerTrees%addMetadata(metaDataTypeProvenance,'fileBuiltBy'       ,     'Galacticus'                                                          )
    call mergerTrees%addMetadata(metaDataTypeProvenance,'fileTimestamp'     ,char(Formatted_Date_and_Time                          (                  )))
    do i=1,size(self%metaData)
       ! Determine the metadata type.
       metaDataText=self%metaData(i)%content
       !! Attempt reading as an integer.
       read (metaDataText,*,iostat=ioStatus) metaDataValueInteger
       if (ioStatus == 0) then
          call mergerTrees%addMetadata(self%metaData(i)%type,char(self%metaData(i)%name),     metaDataValueInteger )
          cycle
       end if
       !! Attempt reading as a float.
       read (metaDataText,*,iostat=ioStatus) metaDataValueFloat
       if (ioStatus == 0) then
          call mergerTrees%addMetadata(self%metaData(i)%type,char(self%metaData(i)%name),     metaDataValueFloat   )
          cycle
       end if
       !! Assume a string.
       call    mergerTrees%addMetadata(self%metaData(i)%type,char(self%metaData(i)%name),trim(metaDataText        ))
    end do
    ! Process particles if necessary.
    if (self%traceParticles) then
       do i=1,size(self%particleProperties)
          call mergerTrees%setParticlePropertyColumn(self%particleProperties(i)%property,self%particleProperties(i)%column)
       end do
       call mergerTrees%readParticlesASCII(char(self%particlesFileName),columnHeaders=self%columnHeaders,commentCharacter="#",separator=self%columnSeparator)
    end if
    ! Export the trees.
    call mergerTrees%export(char(self%outputFileName),self%outputFormat,int(hdfChunkSize,kind=hsize_t),hdfCompressionLevel)
    ! Done.
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: merger tree file builder' )
    return
  end subroutine mergerTreeFileBuilderPerform

  logical function mergerTreeFileBuilderRequiresOutputFile(self)
    !!{
    Specifies that this task does not require the main output file.
    !!}
    implicit none
    class(taskMergerTreeFileBuilder), intent(inout) :: self
    !$GLC attributes unused :: self

    mergerTreeFileBuilderRequiresOutputFile=.false.
    return
  end function mergerTreeFileBuilderRequiresOutputFile

  subroutine mergerTreeFileBuilderDescriptorSpecial(self,descriptor)
    !!{
    Add special parameters to the descriptor.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: ISO_Varying_String        , only : char
    use :: String_Handling           , only : String_Join
    use :: Merger_Tree_Data_Structure, only : enumerationMetaDataTypeDecode, enumerationPropertyTypeDecode
    implicit none
    class    (taskMergerTreeFileBuilder), intent(inout) :: self
    type     (inputParameters          ), intent(inout) :: descriptor
    type     (inputParameters          )                :: subParameters
    integer                                             :: i
    character(len=12                   )                :: column       , conversionFactor
    
    if (allocated(self%properties)) then
       do i=1,size(self%properties)
          call descriptor%addParameter('property','')
          subParameters=descriptor%subParameters('property',copyInstance=i)
          write (column,'(i12)') self%properties(i)%column
          call    subParameters%addParameter('name'            ,char(enumerationPropertyTypeDecode(self%properties(i)%property,includePrefix=.false.)))
          call    subParameters%addParameter('column'          ,trim(column                                                                          ))
          if (self%properties(i)%conversionFactor%isSet) then
             write (conversionFactor,'(e12.6)') self%properties(i)%conversionFactor%value
             call subParameters%addParameter('conversionFactor',trim(conversionFactor                                                                ))
          end if
       end do
    end if
    if (allocated(self%particleProperties)) then
       do i=1,size(self%particleProperties)
          call descriptor%addParameter('particleProperty','')
          subParameters=descriptor%subParameters('particleProperty',copyInstance=i)
           write (column,'(i12)') self%particleProperties(i)%column
           call subParameters%addParameter('name'  ,char(enumerationPropertyTypeDecode(self%particleProperties(i)%property,includePrefix=.false.)))
           call subParameters%addParameter('column',trim(column                                                                                  ))
        end do
    end if
    if (allocated(self%metaData)) then
       do i=1,size(self%metaData)
          call descriptor%addParameter('metaData','')
          subParameters=descriptor%subParameters('metaData',copyInstance=i)
          call subParameters%addParameter('name'   ,char(                              self%metaData(i)%name                          ))
          call subParameters%addParameter('content',char(                              self%metaData(i)%content                       ))
          call subParameters%addParameter('type'   ,char(enumerationMetaDataTypeDecode(self%metaData(i)%type   ,includePrefix=.false.)))
        end do
    end if
    return
  end subroutine mergerTreeFileBuilderDescriptorSpecial
