!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  use :: Cosmological_Density_Field, only : cosmologicalMassVariance, cosmologicalMassVarianceClass
  use :: Cosmology_Parameters      , only : cosmologyParameters     , cosmologyParametersClass
  use :: Power_Spectra_Primordial  , only : powerSpectrumPrimordial , powerSpectrumPrimordialClass
  use :: Stateful_Types            , only : statefulDouble          , statefulInteger              , statefulLogical
  use :: Transfer_Functions        , only : transferFunction        , transferFunctionClass

  type :: mergerTreeMetaData
     !% Type used for metadata in merger tree files.
     integer                 :: type
     type   (varying_string) :: name, content
  end type mergerTreeMetaData

  type :: propertyColumn
     !% Type used to specify which property to read from which column,
     integer                 :: property        , column
     type   (statefulDouble) :: conversionFactor
  end type propertyColumn

  !# <task name="taskMergerTreeFileBuilder">
  !#  <description>A task which computes and outputs the halo mass function and related quantities.</description>
  !# </task>
  type, extends(taskClass) :: taskMergerTreeFileBuilder
     !% Implementation of a task which computes and outputs the halo mass function and related quantities.
     private
     class           (cosmologyParametersClass     ), pointer                   :: cosmologyParameters_
     class           (cosmologicalMassVarianceClass), pointer                   :: cosmologicalMassVariance_
     class           (powerSpectrumPrimordialClass ), pointer                   :: powerSpectrumPrimordial_
     class           (transferFunctionClass        ), pointer                   :: transferFunction_
     type            (propertyColumn               ), allocatable, dimension(:) :: properties                 , particleProperties
     type            (mergerTreeMetaData           ), allocatable, dimension(:) :: metaData
     type            (varying_string               )                            :: outputFileName             , inputFileName                   , &
          &                                                                        particlesFileName
     integer                                                                    :: outputFormat
     double precision                                                           :: unitsMassInSI              , unitsLengthInSI                 , &
          &                                                                        unitsVelocityInSI
     integer                                                                    :: unitsMassHubbleExponent    , unitsMassScaleFactorExponent    , &
          &                                                                        unitsLengthHubbleExponent  , unitsLengthScaleFactorExponent  , &
          &                                                                        unitsVelocityHubbleExponent, unitsVelocityScaleFactorExponent
     type            (varying_string               )                            :: unitsMassName              , unitsLengthName                 , &
          &                                                                        unitsVelocityName
     type            (statefulInteger              )                            :: dummyHostId
     type            (statefulDouble               )                            :: massParticle
     type            (statefulLogical              )                            :: haloMassesIncludeSubhalos  , includesHubbleFlow              , &
          &                                                                        positionsArePeriodic
     logical                                                                    :: traceParticles             , columnHeaders
     character       (len=1                        )                            :: columnSeparator
   contains
     final     ::                       mergerTreeFileBuilderDestructor
     procedure :: perform            => mergerTreeFileBuilderPerform
     procedure :: requiresOutputFile => mergerTreeFileBuilderRequiresOutputFile
  end type taskMergerTreeFileBuilder

  interface taskMergerTreeFileBuilder
     !% Constructors for the {\normalfont \ttfamily mergerTreeFileBuilder} task.
     module procedure mergerTreeFileBuilderConstructorParameters
     module procedure mergerTreeFileBuilderConstructorInternal
  end interface taskMergerTreeFileBuilder

contains

  function mergerTreeFileBuilderConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily mergerTreeFileBuilder} task class which takes a parameter set as input.
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

    !# <inputParameter>
    !#   <name>inputFileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the file from which to read merger tree data.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    if (parameters%isPresent('particlesFileName')) then
       !# <inputParameter>
       !#   <name>particlesFileName</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The name of the file from which to read particle data.</description>
       !#   <source>parameters</source>
       !#   <type>string</type>
       !# </inputParameter>
    else
       particlesFileName=''
    end if
    !# <inputParameter>
    !#   <name>outputFileName</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The name of the file to which to write merger tree data.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>outputFormat</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The format to use for merger tree output.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>columnHeaders</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>If true, the file is assumed to contain a single line of column headers, which will be skipped.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>columnSeparator</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str(',')</defaultValue>
    !#   <description>The separator for columns.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    if (columnSeparator == "space") columnSeparator=" "
    massParticle%isSet=parameters%isPresent('massParticle')
    if (massParticle%isSet) then
       !# <inputParameter>
       !#   <name>massParticle</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The mass of the simulation particle.</description>
       !#   <source>parameters</source>
       !#   <variable>massParticle%value</variable>
       !#   <type>real</type>
       !# </inputParameter>
    end if
    dummyHostId%isSet=parameters%isPresent('dummyHostId')
    if (dummyHostId%isSet) then
       !# <inputParameter>
       !#   <name>dummyHostId</name>
       !#   <cardinality>1</cardinality>
       !#   <description>If present, specifies the dummy host ID for self-hosting halos. Otherwise, self-hosting halos have {\normalfont \ttfamily hostIndex == nodeIndex}.</description>
       !#   <source>parameters</source>
       !#   <variable>dummyHostId%value</variable>
       !#   <type>real</type>
       !# </inputParameter>
    end if
    haloMassesIncludeSubhalos%isSet=parameters%isPresent('haloMassesIncludeSubhalos')
    if (haloMassesIncludeSubhalos%isSet) then
       !# <inputParameter>
       !#   <name>haloMassesIncludeSubhalos</name>
       !#   <cardinality>1</cardinality>
       !#   <description>Specifies whether or not halo masses include the masses of their subhalos.</description>
       !#   <source>parameters</source>
       !#   <variable>haloMassesIncludeSubhalos%value</variable>
       !#   <type>boolean</type>
       !# </inputParameter>
    end if
    if (includesHubbleFlow%isSet) then
       !# <inputParameter>
       !#   <name>includesHubbleFlow</name>
       !#   <cardinality>1</cardinality>
       !#   <description>Specifies whether or not Hubble flow is included in velocities.</description>
       !#   <source>parameters</source>
       !#   <variable>includesHubbleFlow%value</variable>
       !#   <type>boolean</type>
       !# </inputParameter>
    end if
    if (positionsArePeriodic%isSet) then
       !# <inputParameter>
       !#   <name>positionsArePeriodic</name>
       !#   <cardinality>1</cardinality>
       !#   <description>Specifies whether or not positions are periodic.</description>
       !#   <source>parameters</source>
       !#   <variable>positionsArePeriodic%value</variable>
       !#   <type>boolean</type>
       !# </inputParameter>
    end if
    !# <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !# <objectBuilder class="powerSpectrumPrimordial"  name="powerSpectrumPrimordial_"  source="parameters"/>
    !# <objectBuilder class="transferFunction"         name="transferFunction_"         source="parameters"/>
    countProperties=parameters%copiesCount('property',requireValue=.false.)
    allocate(properties(countProperties))
    do i=1,countProperties
       allocate(subParameters)
       subParameters=parameters%subParameters('property',requireValue=.false.,copyInstance=i)
       !# <inputParameter>
       !#   <name>name</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The name of the property to read.</description>
       !#   <source>subParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>column</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The column from which to read the property.</description>
       !#   <source>subParameters</source>
       !#   <type>integer</type>
       !# </inputParameter>
       properties(i)%property              =enumerationPropertyTypeEncode(char(name),includesPrefix=.false.)
       properties(i)%column                =column
       properties(i)%conversionFactor%isSet=subParameters%isPresent('conversionFactor')
       if (properties(i)%conversionFactor%isSet) then
          !# <inputParameter>
          !#   <name>conversionFactor</name>
          !#   <variable>properties(i)%conversionFactor%value</variable>
          !#   <cardinality>1</cardinality>
          !#   <description>An additional conversion factor to apply to the property to get it into the correct units.</description>
          !#   <source>subParameters</source>
          !#   <type>float</type>
          !# </inputParameter>
       end if
       deallocate(subParameters)
    end do
    countProperties=parameters%copiesCount('particleProperty',requireValue=.false.,zeroIfNotPresent=.true.)
    allocate(particleProperties(countProperties))
    do i=1,countProperties
       allocate(subParameters)
       subParameters=parameters%subParameters('particleProperty',requireValue=.false.,copyInstance=i)
       !# <inputParameter>
       !#   <name>name</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The name of the particle property to read.</description>
       !#   <source>subParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>column</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The column from which to read the particle property.</description>
       !#   <source>subParameters</source>
       !#   <type>integer</type>
       !# </inputParameter>
       deallocate(subParameters)
       particleProperties(i)%property=enumerationPropertyTypeEncode(char(name),includesPrefix=.false.)
       particleProperties(i)%column  =column
    end do
    countMetaData=parameters%copiesCount('metaData',requireValue=.false.,zeroIfNotPresent=.true.)
    allocate(metaData(countMetaData))
    do i=1,countMetaData
       allocate(subParameters)
       subParameters=parameters%subParameters('metaData',requireValue=.false.,copyInstance=i)
       !# <inputParameter>
       !#   <name>name</name>
       !#   <cardinality>1</cardinality>
       !#   <description>The name of the metadata.</description>
       !#   <source>subParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>content</name>
       !#   <variable>metaDataContent</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The value of the metadata.</description>
       !#   <source>subParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>type</name>
       !#   <variable>metaDataType</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The metadata type.</description>
       !#   <source>subParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
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
       !# <inputParameter>
       !#   <name>unitsInSI</name>
       !#   <variable>unitsMassInSI</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The mass unit in the SI system.</description>
       !#   <source>subParameters</source>
       !#   <type>float</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>hubbleExponent</name>
       !#   <variable>unitsMassHubbleExponent</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The exponent of the ``little-$h$'' Hubble parameter needed to convert the masses to little-$h$-free units.</description>
       !#   <source>subParameters</source>
       !#   <type>integer</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>scaleFactorExponent</name>
       !#   <variable>unitsMassScaleFactorExponent</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The exponent of the cosmological scale factor needed to convert the masses to physical units.</description>
       !#   <source>subParameters</source>
       !#   <type>integer</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>name</name>
       !#   <variable>unitsMassName</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>A human-readable name for the units of mass.</description>
       !#   <source>subParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
       deallocate(subParameters)
    end if
    unitsLengthInSI               =megaParsec
    unitsLengthHubbleExponent     =0
    unitsLengthScaleFactorExponent=0
    unitsLengthName               =''
    if (parameters%isPresent('unitsLength',requireValue=.false.)) then
       allocate(subParameters)
       subParameters=parameters%subParameters('unitsLength',requireValue=.false.)
       !# <inputParameter>
       !#   <name>unitsInSI</name>
       !#   <variable>unitsLengthInSI</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The length unit in the SI system.</description>
       !#   <source>subParameters</source>
       !#   <type>float</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>hubbleExponent</name>
       !#   <variable>unitsLengthHubbleExponent</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The exponent of the ``little-$h$'' Hubble parameter needed to convert the lengthes to little-$h$-free units.</description>
       !#   <source>subParameters</source>
       !#   <type>integer</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>scaleFactorExponent</name>
       !#   <variable>unitsLengthScaleFactorExponent</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The exponent of the cosmological scale factor needed to convert the lengthes to physical units.</description>
       !#   <source>subParameters</source>
       !#   <type>integer</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>name</name>
       !#   <variable>unitsLengthName</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>A human-readable name for the units of length.</description>
       !#   <source>subParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
       deallocate(subParameters)
    end if
    unitsVelocityInSI               =kilo
    unitsVelocityHubbleExponent     =0
    unitsVelocityScaleFactorExponent=0
    unitsVelocityName               =''
    if (parameters%isPresent('unitsVelocity',requireValue=.false.)) then
       allocate(subParameters)
       subParameters=parameters%subParameters('unitsVelocity',requireValue=.false.)
       !# <inputParameter>
       !#   <name>unitsInSI</name>
       !#   <variable>unitsVelocityInSI</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The velocity unit in the SI system.</description>
       !#   <source>subParameters</source>
       !#   <type>float</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>hubbleExponent</name>
       !#   <variable>unitsVelocityHubbleExponent</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The exponent of the ``little-$h$'' Hubble parameter needed to convert the velocities to little-$h$-free units.</description>
       !#   <source>subParameters</source>
       !#   <type>integer</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>scaleFactorExponent</name>
       !#   <variable>unitsVelocityScaleFactorExponent</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>The exponent of the cosmological scale factor needed to convert the velocities to physical units.</description>
       !#   <source>subParameters</source>
       !#   <type>integer</type>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>name</name>
       !#   <variable>unitsVelocityName</variable>
       !#   <cardinality>1</cardinality>
       !#   <description>A human-readable name for the units of velocity.</description>
       !#   <source>subParameters</source>
       !#   <type>string</type>
       !# </inputParameter>
       deallocate(subParameters)
    end if
    self=taskMergerTreeFileBuilder(inputFileName,particlesFileName,outputFileName,enumerationMergerTreeFormatEncode(char(outputFormat),includesPrefix=.false.),columnHeaders,char(columnSeparator),properties,particleProperties,metaData,massParticle,dummyHostId,haloMassesIncludeSubhalos,includesHubbleFlow,positionsArePeriodic,unitsMassInSI,unitsMassHubbleExponent,unitsMassScaleFactorExponent,unitsMassName,unitsLengthInSI,unitsLengthHubbleExponent,unitsLengthScaleFactorExponent,unitsLengthName,unitsVelocityInSI,unitsVelocityHubbleExponent,unitsVelocityScaleFactorExponent,unitsVelocityName,cosmologyParameters_,cosmologicalMassVariance_,powerSpectrumPrimordial_,transferFunction_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"     />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    !# <objectDestructor name="powerSpectrumPrimordial_" />
    !# <objectDestructor name="transferFunction_"        />
    return
  end function mergerTreeFileBuilderConstructorParameters

  function mergerTreeFileBuilderConstructorInternal(inputFileName,particlesFileName,outputFileName,outputFormat,columnHeaders,columnSeparator,properties,particleProperties,metaData,massParticle,dummyHostId,haloMassesIncludeSubhalos,includesHubbleFlow,positionsArePeriodic,unitsMassInSI,unitsMassHubbleExponent,unitsMassScaleFactorExponent,unitsMassName,unitsLengthInSI,unitsLengthHubbleExponent,unitsLengthScaleFactorExponent,unitsLengthName,unitsVelocityInSI,unitsVelocityHubbleExponent,unitsVelocityScaleFactorExponent,unitsVelocityName,cosmologyParameters_,cosmologicalMassVariance_,powerSpectrumPrimordial_,transferFunction_) result(self)
    !% Constructor for the {\normalfont \ttfamily mergerTreeFileBuilder} task class which takes a parameter set as input.
    implicit none
    type            (taskMergerTreeFileBuilder    )                              :: self
    class           (cosmologyParametersClass     ), intent(in   ), target       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass), intent(in   ), target       :: cosmologicalMassVariance_
    class           (powerSpectrumPrimordialClass ), intent(in   ), target       :: powerSpectrumPrimordial_
    class           (transferFunctionClass        ), intent(in   ), target       :: transferFunction_
    type            (varying_string               ), intent(in   )               :: inputFileName              , outputFileName                  , &
         &                                                                          particlesFileName
    type            (statefulInteger              ), intent(in   )               :: dummyHostId
    type            (statefulDouble               ), intent(in   )               :: massParticle
    type            (statefulLogical              ), intent(in   )               :: haloMassesIncludeSubhalos  , includesHubbleFlow              , &
         &                                                                          positionsArePeriodic
    integer                                        , intent(in   )               :: outputFormat
    type            (propertyColumn               ), intent(in   ), dimension(:) :: properties                 , particleProperties
    double precision                               , intent(in   )               :: unitsMassInSI              , unitsLengthInSI                 , &
         &                                                                          unitsVelocityInSI
    type            (mergerTreeMetaData           ), intent(in   ), dimension(:) :: metaData
    integer                                        , intent(in   )               :: unitsMassHubbleExponent    , unitsMassScaleFactorExponent    , &
         &                                                                          unitsLengthHubbleExponent  , unitsLengthScaleFactorExponent  , &
         &                                                                          unitsVelocityHubbleExponent, unitsVelocityScaleFactorExponent
    type            (varying_string               ), intent(in   )               :: unitsMassName              , unitsLengthName                 , &
         &                                                                          unitsVelocityName
    logical                                        , intent(in   )               :: columnHeaders
    character       (len=1                        ), intent(in   )               :: columnSeparator
    !# <constructorAssign variables="inputFileName, particlesFileName, outputFileName, outputFormat, columnHeaders, columnSeparator, properties, particleProperties, metaData, massParticle, dummyHostId, haloMassesIncludeSubhalos, includesHubbleFlow, positionsArePeriodic, unitsMassInSI, unitsMassHubbleExponent, unitsMassScaleFactorExponent, unitsMassName, unitsLengthInSI, unitsLengthHubbleExponent, unitsLengthScaleFactorExponent, unitsLengthName, unitsVelocityInSI, unitsVelocityHubbleExponent, unitsVelocityScaleFactorExponent, unitsVelocityName, *cosmologyParameters_, *cosmologicalMassVariance_, *powerSpectrumPrimordial_, *transferFunction_"/>

    self%traceParticles=particlesFileName /= ''
    return
  end function mergerTreeFileBuilderConstructorInternal

  subroutine mergerTreeFileBuilderDestructor(self)
    !% Destructor for the {\normalfont \ttfamily mergerTreeFileBuilder} task class.
    implicit none
    type(taskMergerTreeFileBuilder), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"     />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    !# <objectDestructor name="self%powerSpectrumPrimordial_" />
    !# <objectDestructor name="self%transferFunction_"        />
    return
  end subroutine mergerTreeFileBuilderDestructor

  subroutine mergerTreeFileBuilderPerform(self,status)
    !% Compute and output the halo mass function.
    use :: Cosmology_Parameters      , only : hubbleUnitsLittleH
    use :: Dates_and_Times           , only : Formatted_Date_and_Time
    use :: Galacticus_Display        , only : Galacticus_Display_Indent, Galacticus_Display_Unindent
    use :: Galacticus_Error          , only : errorStatusSuccess
    use :: HDF5                      , only : hsize_t
    use :: Input_Parameters          , only : inputParameters
    use :: Merger_Tree_Data_Structure, only : mergerTreeData           , metaDataTypeCosmology      , metaDataTypeProvenance, metaDataTypeSimulation, &
          &                                   unitsLength              , unitsMass                  , unitsVelocity
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

    call Galacticus_Display_Indent('Begin task: merger tree file builder')
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
    ! Read in the data.
    call mergerTrees%readASCII(char(self%inputFileName),columnHeaders=self%columnHeaders,commentCharacter="#",separator=self%columnSeparator)
    ! Specify that we do not want to create individual merger tree reference datasets.
    call mergerTrees%makeReferences          (.false.)
    ! Specify that trees are self-contained (i.e. nodes never move from one tree to another).
    call mergerTrees%setSelfContained        (.true. )
    call mergerTrees%setUnits(unitsMass    ,unitsInSI=self%unitsMassInSI    ,hubbleExponent=self%unitsMassHubbleExponent    ,scaleFactorExponent=self%unitsMassScaleFactorExponent    ,name=char(self%unitsMassName    ))
    call mergerTrees%setUnits(unitsLength  ,unitsInSI=self%unitsLengthInSI  ,hubbleExponent=self%unitsLengthHubbleExponent  ,scaleFactorExponent=self%unitsLengthScaleFactorExponent  ,name=char(self%unitsLengthName  ))
    call mergerTrees%setUnits(unitsVelocity,unitsInSI=self%unitsVelocityInSI,hubbleExponent=self%unitsVelocityHubbleExponent,scaleFactorExponent=self%unitsVelocityScaleFactorExponent,name=char(self%unitsVelocityName))
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
    ! Get descriptors.
    descriptorPowerSpectrum   =inputParameters()
    descriptorTransferFunction=inputParameters()
    call self%powerSpectrumPrimordial_%descriptor(descriptorPowerSpectrum   ,includeMethod=.true.)
    call self%transferFunction_       %descriptor(descriptorTransferFunction,includeMethod=.true.)
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
    call Galacticus_Display_Unindent('Done task: merger tree file builder' )
    return
  end subroutine mergerTreeFileBuilderPerform

  logical function mergerTreeFileBuilderRequiresOutputFile(self)
    !% Specifies that this task does not require the main output file.
    implicit none
    class(taskMergerTreeFileBuilder), intent(inout) :: self
    !$GLC attributes unused :: self

    mergerTreeFileBuilderRequiresOutputFile=.false.
    return
  end function mergerTreeFileBuilderRequiresOutputFile
