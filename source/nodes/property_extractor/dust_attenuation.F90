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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a property extractor which applies dust attenuation to the luminosities of other extractors.
  !!}

  use :: Dust_Attenuations, only : dustAttenuationClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorDustAttenuation" docformat="rst">
   <description>
   A property extractor which applies dust attenuation to the luminosities produced by other property extractors.

   Each child extractor is asked to decompose its luminosity into parcels of emission---one wavelength, one component,
   one range of stellar population ages---at a resolution which the ``[dustAttenuation]`` object itself specifies. The
   attenuator returns a transmission factor for each parcel, and the child recombines the attenuated parcels into its
   own output. A child which can not decompose its output, or whose component the attenuator refuses, is rejected.

   Attenuated properties are named for the child with ``:dustAttenuated:`` and the name of the attenuator appended, so
   they may coexist with their unattenuated counterparts. Setting ``[outputUnattenuated]`` emits those counterparts
   alongside, which is the usual way to measure the effect of the dust model. Setting ``[outputSum]`` additionally
   emits the sum of the attenuated luminosities of all children: dust must be applied to each component separately,
   since components are attenuated differently, so recovering a galaxy-wide total requires summing afterwards, and
   this does that in one place rather than in post-processing.

   Only children producing scalar or tuple properties are supported at present. Extractors returning arrays---spectral
   energy distributions in particular---will be supported once they implement decomposition.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorMulti) :: nodePropertyExtractorDustAttenuation
     !!{RST
     A property extractor which applies dust attenuation to the luminosities of other extractors.
     !!}
     private
     class  (dustAttenuationClass), pointer :: dustAttenuation_   => null()
     logical                                :: outputUnattenuated          , outputSum
   contains
     !![
     <methods docformat="rst">
       <method method="countChildElements" description="Return the number of properties produced by a child extractor."/>
       <method method="attenuate"          description="Return the attenuated properties of a child extractor."        />
     </methods>
     !!]
     final     ::                       dustAttenuationDestructor
     procedure :: countChildElements => dustAttenuationCountChildElements
     procedure :: attenuate          => dustAttenuationAttenuate
     procedure :: elementCount       => dustAttenuationElementCount
     procedure :: extractDouble      => dustAttenuationExtractDouble
     procedure :: names              => dustAttenuationNames
     procedure :: descriptions       => dustAttenuationDescriptions
     procedure :: unitsInSI          => dustAttenuationUnitsInSI
     procedure :: units              => dustAttenuationUnits
     procedure :: ranks              => dustAttenuationRanks
  end type nodePropertyExtractorDustAttenuation

  interface nodePropertyExtractorDustAttenuation
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorDustAttenuation` property extractor class.
     !!}
     module procedure dustAttenuationConstructorParameters
     module procedure dustAttenuationConstructorInternal
  end interface nodePropertyExtractorDustAttenuation

contains

  function dustAttenuationConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorDustAttenuation` property extractor class which takes
    a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorDustAttenuation)                :: self
    type   (inputParameters                     ), intent(inout) :: parameters
    class  (dustAttenuationClass                ), pointer       :: dustAttenuation_
    logical                                                      :: outputUnattenuated, outputSum

    self%nodePropertyExtractorMulti=nodePropertyExtractorMulti(parameters)
    !![
    <inputParameter docformat="rst">
      <name>outputUnattenuated</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, the unattenuated properties of each child extractor are emitted alongside the attenuated ones.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputSum</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, the sum of the attenuated properties over all child extractors is emitted as an additional property.
      Useful for recovering a galaxy-wide total from luminosities which had to be attenuated component by component.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="dustAttenuation" name="dustAttenuation_" source="parameters"/>
    <inputParametersValidate source="parameters" multiParameters="nodePropertyExtractor"/>
    !!]
    self%dustAttenuation_   => dustAttenuation_
    !![
    <referenceCountIncrement owner="self" object="dustAttenuation_"/>
    !!]
    self%outputUnattenuated =  outputUnattenuated
    self%outputSum          =  outputSum
    call dustAttenuationValidate(self)
    !![
    <objectDestructor name="dustAttenuation_"/>
    !!]
    return
  end function dustAttenuationConstructorParameters

  function dustAttenuationConstructorInternal(dustAttenuation_,outputUnattenuated,outputSum,extractors) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorDustAttenuation` property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorDustAttenuation)                        :: self
    class  (dustAttenuationClass                ), intent(in   ), target :: dustAttenuation_
    logical                                      , intent(in   )         :: outputUnattenuated, outputSum
    type   (multiExtractorList                  ), intent(in   )         :: extractors
    !![
    <constructorAssign variables="outputUnattenuated, outputSum, *dustAttenuation_"/>
    !!]

    self%nodePropertyExtractorMulti=nodePropertyExtractorMulti(extractors)
    call dustAttenuationValidate(self)
    return
  end function dustAttenuationConstructorInternal

  subroutine dustAttenuationValidate(self)
    !!{RST
    Check that every child extractor can be attenuated, reporting any that can not at construction rather than part
    way through a run.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type (nodePropertyExtractorDustAttenuation), intent(inout) :: self
    type (multiExtractorList                  ), pointer       :: extractor_

    extractor_ => self%extractors
    do while (associated(extractor_))
       if (.not.extractor_%extractor_%supportsAttenuation())                                                   &
            & call Error_Report(                                                                               &
            &                   'property extractor "'                                                      // &
            &                   extractor_%extractor_%objectType()                                          // &
            &                   '" does not support dust attenuation - it can not decompose its output into'// &
            &                   ' parcels of emission'                                                      // &
            &                   {introspection:location}                                                       &
            &                  )
       select type (child => extractor_%extractor_)
       class is (nodePropertyExtractorScalar)
          ! Supported.
       class is (nodePropertyExtractorTuple )
          ! Supported.
       class default
          call Error_Report(                                                                                       &
               &            'property extractor "'                                                              // &
               &            extractor_%extractor_%objectType()                                                  // &
               &            '" produces properties of a rank which dust attenuation does not yet support - only'// &
               &            ' scalar and tuple extractors may be attenuated'                                    // &
               &            {introspection:location}                                                               &
               &           )
       end select
       extractor_ => extractor_%next
    end do
    return
  end subroutine dustAttenuationValidate

  subroutine dustAttenuationDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorDustAttenuation` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorDustAttenuation), intent(inout) :: self

    !![
    <objectDestructor name="self%dustAttenuation_"/>
    !!]
    return
  end subroutine dustAttenuationDestructor

  integer function dustAttenuationCountChildElements(self,extractor_,time) result(countElements)
    !!{RST
    Return the number of properties produced by a child extractor.
    !!}
    implicit none
    class           (nodePropertyExtractorDustAttenuation), intent(inout) :: self
    class           (nodePropertyExtractorClass          ), intent(inout) :: extractor_
    double precision                                      , intent(in   ) :: time
    !$GLC attributes unused :: self

    countElements=0
    select type (extractor_)
    class is (nodePropertyExtractorScalar)
       countElements=1
    class is (nodePropertyExtractorTuple )
       countElements=extractor_%elementCount(time)
    end select
    return
  end function dustAttenuationCountChildElements

  integer function dustAttenuationElementCount(self,elementType,time) result(elementCount)
    !!{RST
    Return the number of properties emitted, counting attenuated properties, optionally their unattenuated
    counterparts, and optionally their sum.
    !!}
    implicit none
    class           (nodePropertyExtractorDustAttenuation), intent(inout) :: self
    type            (enumerationElementTypeType          ), intent(in   ) :: elementType
    double precision                                      , intent(in   ) :: time
    type            (multiExtractorList                  ), pointer       :: extractor_
    integer                                                               :: countChild

    elementCount=0
    ! Only double-precision properties are produced.
    if (elementType /= elementTypeDouble) return
    extractor_ => self%extractors
    do while (associated(extractor_))
       countChild  =self%countChildElements(extractor_%extractor_,time)
       elementCount=elementCount+countChild
       if (self%outputUnattenuated) elementCount=elementCount+countChild
       extractor_  => extractor_%next
    end do
    if (self%outputSum) elementCount=elementCount+1
    return
  end function dustAttenuationElementCount

  function dustAttenuationAttenuate(self,extractor_,node,time,instance) result(values)
    !!{RST
    Return the attenuated properties of a child extractor: decompose its luminosity into parcels, attenuate each, and
    let the child recombine them.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                                      , allocatable  , dimension(:) :: values
    class           (nodePropertyExtractorDustAttenuation), intent(inout)               :: self
    class           (nodePropertyExtractorClass          ), intent(inout), target       :: extractor_
    type            (treeNode                            ), intent(inout), target       :: node
    double precision                                      , intent(in   )               :: time
    type            (multiCounter                        ), intent(inout), optional     :: instance
    type            (luminosityDecomposition             )                              :: decomposition
    double precision                                      , allocatable  , dimension(:) :: transmission
    integer                                                                             :: i
    !$GLC attributes unused :: instance

    decomposition=extractor_%decompose(node,time,self%dustAttenuation_%request())
    ! The attenuator may refuse the component a parcel came from -- most refuse a luminosity summed over components,
    ! which can not be attenuated meaningfully. The child's component is not visible until it decomposes, so this is
    ! checked here rather than at construction.
    do i=1,decomposition%countTerms()
       if (.not.self%dustAttenuation_%supportsComponent(decomposition%descriptors(i)%componentType))                &
            & call Error_Report(                                                                                    &
            &                   'the dust attenuation model refuses the component of a parcel of emission from "'// &
            &                   extractor_%objectType()                                                          // &
            &                   '" - a luminosity summed over components can not be attenuated'                  // &
            &                   {introspection:location}                                                            &
            &                  )
    end do
    allocate(transmission(decomposition%countTerms()))
    if (decomposition%countTerms() > 0) transmission=self%dustAttenuation_%transmission(node,decomposition%descriptors)
    call extractor_%recompose(decomposition,transmission,values)
    return
  end function dustAttenuationAttenuate

  function dustAttenuationExtractDouble(self,node,time,instance,ranks) result(properties)
    !!{RST
    Extract the attenuated properties.
    !!}
    use :: Error     , only : Error_Report
    use :: Poly_Ranks, only : polyRankDouble
    implicit none
    type            (polyRankDouble                      )                         , allocatable, dimension(:) :: properties
    class           (nodePropertyExtractorDustAttenuation), intent(inout)                                      :: self
    type            (treeNode                            ), intent(inout)                                      :: node
    double precision                                      , intent(in   )                                      :: time
    type            (multiCounter                        ), intent(inout), optional                            :: instance
    integer                                               , intent(  out), optional, allocatable, dimension(:) :: ranks
    type            (multiExtractorList                  ), pointer                                            :: extractor_
    double precision                                                               , allocatable, dimension(:) :: valuesAttenuated, valuesRaw
    integer                                                                                                    :: offset          , countChild, &
         &                                                                                                        i
    double precision                                                                                           :: sumAttenuated

    allocate(properties(self%elementCount(elementTypeDouble,time)))
    if (present(ranks)) then
       allocate(ranks(self%elementCount(elementTypeDouble,time)))
       ranks=0
    end if
    offset        =  0
    sumAttenuated =  0.0d0
    extractor_    => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       if (self%outputUnattenuated) then
          select type (child => extractor_%extractor_)
          class is (nodePropertyExtractorScalar)
             allocate(valuesRaw(1))
             valuesRaw(1)=child%extract(node,     instance)
          class is (nodePropertyExtractorTuple )
             valuesRaw   =child%extract(node,time,instance)
          end select
          do i=1,countChild
             properties(offset+i)=polyRankDouble(valuesRaw(i))
          end do
          offset=offset+countChild
          deallocate(valuesRaw)
       end if
       valuesAttenuated=self%attenuate(extractor_%extractor_,node,time,instance)
       if (size(valuesAttenuated) /= countChild)                                       &
            & call Error_Report(                                                       &
            &                   'the number of attenuated properties returned by "' // &
            &                   extractor_%extractor_%objectType()                  // &
            &                   '" does not match the number of properties it emits'// &
            &                   {introspection:location}                               &
            &                  )
       do i=1,countChild
          properties(offset+i)=polyRankDouble(valuesAttenuated(i))
          sumAttenuated       =sumAttenuated+valuesAttenuated(i)
       end do
       offset=offset+countChild
       deallocate(valuesAttenuated)
       extractor_ => extractor_%next
    end do
    if (self%outputSum) then
       offset            =offset+1
       properties(offset)=polyRankDouble(sumAttenuated)
    end if
    return
  end function dustAttenuationExtractDouble

  subroutine dustAttenuationNames(self,elementType,time,names)
    !!{RST
    Return the names of the properties emitted. Attenuated properties take the name of the property they attenuate,
    with the name of the attenuator appended, so that the two may coexist in the same output.
    !!}
    implicit none
    class           (nodePropertyExtractorDustAttenuation), intent(inout)                            :: self
    type            (enumerationElementTypeType          ), intent(in   )                            :: elementType
    double precision                                      , intent(in   )                            :: time
    type            (varying_string                      ), intent(inout), dimension(:), allocatable :: names
    type            (varying_string                      )               , dimension(:), allocatable :: namesChild
    type            (multiExtractorList                  ), pointer                                  :: extractor_
    type            (varying_string                      )                                           :: suffix
    integer                                                                                          :: offset     , countChild, &
         &                                                                                               i

    allocate(names(self%elementCount(elementType,time)))
    if (elementType /= elementTypeDouble) return
    suffix     =  ":dustAttenuated:"//self%dustAttenuation_%objectType(short=.true.)
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       call dustAttenuationChildNames(self,extractor_%extractor_,time,namesChild)
       if (self%outputUnattenuated) then
          do i=1,countChild
             names(offset+i)=namesChild(i)
          end do
          offset=offset+countChild
       end if
       do i=1,countChild
          names(offset+i)=namesChild(i)//suffix
       end do
       offset=offset+countChild
       deallocate(namesChild)
       extractor_ => extractor_%next
    end do
    if (self%outputSum) then
       offset       =offset+1
       names(offset)="luminosityTotal"//suffix
    end if
    return
  end subroutine dustAttenuationNames

  subroutine dustAttenuationChildNames(self,extractor_,time,names)
    !!{RST
    Return the names of the properties of a child extractor.
    !!}
    implicit none
    type            (nodePropertyExtractorDustAttenuation), intent(inout)                            :: self
    class           (nodePropertyExtractorClass          ), intent(inout)                            :: extractor_
    double precision                                      , intent(in   )                            :: time
    type            (varying_string                      ), intent(inout), dimension(:), allocatable :: names
    !$GLC attributes unused :: self

    select type (extractor_)
    class is (nodePropertyExtractorScalar)
       allocate(names(1))
       names(1)=extractor_%name()
    class is (nodePropertyExtractorTuple )
       call extractor_%names(time,names)
    end select
    return
  end subroutine dustAttenuationChildNames

  subroutine dustAttenuationDescriptions(self,elementType,time,descriptions)
    !!{RST
    Return descriptions of the properties emitted.
    !!}
    implicit none
    class           (nodePropertyExtractorDustAttenuation), intent(inout)                            :: self
    type            (enumerationElementTypeType          ), intent(in   )                            :: elementType
    double precision                                      , intent(in   )                            :: time
    type            (varying_string                      ), intent(inout), dimension(:), allocatable :: descriptions
    type            (varying_string                      )               , dimension(:), allocatable :: descriptionsChild
    type            (multiExtractorList                  ), pointer                                  :: extractor_
    type            (varying_string                      )                                           :: suffix
    integer                                                                                          :: offset           , countChild, &
         &                                                                                              i

    allocate(descriptions(self%elementCount(elementType,time)))
    if (elementType /= elementTypeDouble) return
    suffix     =  ", attenuated by dust using the '"//self%dustAttenuation_%objectType(short=.true.)//"' model"
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       select type (child => extractor_%extractor_)
       class is (nodePropertyExtractorScalar)
          allocate(descriptionsChild(1))
          descriptionsChild(1)=child%description()
       class is (nodePropertyExtractorTuple )
          call child%descriptions(time,descriptionsChild)
       end select
       if (self%outputUnattenuated) then
          do i=1,countChild
             descriptions(offset+i)=descriptionsChild(i)
          end do
          offset=offset+countChild
       end if
       do i=1,countChild
          descriptions(offset+i)=descriptionsChild(i)//suffix
       end do
       offset=offset+countChild
       deallocate(descriptionsChild)
       extractor_ => extractor_%next
    end do
    if (self%outputSum) then
       offset            =offset+1
       descriptions(offset)="Sum over all child extractors of the luminosity"//suffix
    end if
    return
  end subroutine dustAttenuationDescriptions

  function dustAttenuationUnitsInSI(self,elementType,time) result(unitsInSI)
    !!{RST
    Return the units, in the SI system, of the properties emitted. Attenuation is dimensionless, so an attenuated
    property carries the units of the property it attenuates.
    !!}
    implicit none
    double precision                                      , allocatable  , dimension(:) :: unitsInSI
    class           (nodePropertyExtractorDustAttenuation), intent(inout)               :: self
    type            (enumerationElementTypeType          ), intent(in   )               :: elementType
    double precision                                      , intent(in   )               :: time
    double precision                                      , allocatable  , dimension(:) :: unitsChild
    type            (multiExtractorList                  ), pointer                     :: extractor_
    integer                                                                             :: offset     , countChild

    allocate(unitsInSI(self%elementCount(elementType,time)))
    if (elementType /= elementTypeDouble) return
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       select type (child => extractor_%extractor_)
       class is (nodePropertyExtractorScalar)
          allocate(unitsChild(1))
          unitsChild(1)=child%unitsInSI()
       class is (nodePropertyExtractorTuple )
          unitsChild   =child%unitsInSI(time)
       end select
       if (self%outputUnattenuated) then
          unitsInSI(offset+1:offset+countChild)=unitsChild
          offset                               =offset+countChild
       end if
       unitsInSI(offset+1:offset+countChild)=unitsChild
       offset                               =offset+countChild
       deallocate(unitsChild)
       extractor_ => extractor_%next
    end do
    ! The summed property takes the units of the first child, all children being luminosities in the same units.
    if (self%outputSum) then
       offset           =offset+1
       if (offset > 1) then
          unitsInSI(offset)=unitsInSI(offset-1)
       else
          unitsInSI(offset)=1.0d0
       end if
    end if
    return
  end function dustAttenuationUnitsInSI

  function dustAttenuationUnits(self,elementType,time) result(units)
    !!{RST
    Return metadata describing the units of the properties emitted.
    !!}
    implicit none
    type            (unitType                            ), allocatable  , dimension(:) :: units
    class           (nodePropertyExtractorDustAttenuation), intent(inout)               :: self
    type            (enumerationElementTypeType          ), intent(in   )               :: elementType
    double precision                                      , intent(in   )               :: time
    type            (unitType                            ), allocatable  , dimension(:) :: unitsChild
    type            (multiExtractorList                  ), pointer                     :: extractor_
    integer                                                                             :: offset     , countChild, &
         &                                                                                 i

    allocate(units(self%elementCount(elementType,time)))
    if (elementType /= elementTypeDouble) return
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       ! `units` is defined on the rank-specific extractor classes rather than on the base class, and with differing
       ! signatures, so the child's class must be resolved before it can be called.
       select type (child => extractor_%extractor_)
       class is (nodePropertyExtractorScalar)
          allocate(unitsChild(1))
          unitsChild(1)=child%units(    )
       class is (nodePropertyExtractorTuple )
          unitsChild   =child%units(time)
       end select
       if (self%outputUnattenuated) then
          do i=1,countChild
             units(offset+i)=unitsChild(i)
          end do
          offset=offset+countChild
       end if
       do i=1,countChild
          units(offset+i)=unitsChild(i)
       end do
       offset=offset+countChild
       deallocate(unitsChild)
       extractor_ => extractor_%next
    end do
    if (self%outputSum) then
       offset      =offset+1
       if (offset > 1) units(offset)=units(offset-1)
    end if
    return
  end function dustAttenuationUnits

  function dustAttenuationRanks(self,elementType,time) result(ranks)
    !!{RST
    Return the ranks of the properties emitted. Only scalar and tuple children are supported, so every property is of
    rank zero.
    !!}
    implicit none
    integer                                               , allocatable  , dimension(:) :: ranks
    class           (nodePropertyExtractorDustAttenuation), intent(inout)               :: self
    type            (enumerationElementTypeType          ), intent(in   )               :: elementType
    double precision                                      , intent(in   )               :: time

    allocate(ranks(self%elementCount(elementType,time)))
    ranks=0
    return
  end function dustAttenuationRanks
