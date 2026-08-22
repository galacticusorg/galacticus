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
   emits the sum of the attenuated luminosities of all children, and ``[outputSumOnly]`` emits *only* that sum. Dust
   must be applied to each component separately, since components are attenuated differently, so recovering a
   galaxy-wide total requires summing afterwards; these options do that here rather than in post-processing.
   ``[outputSumOnly]`` in particular is what lets a consumer which expects a single value per galaxy, such as an
   output analysis, use this extractor at all.

   Children producing scalar, tuple, or array properties may all be attenuated. The sum is formed elementwise, so
   spectra are added wavelength by wavelength and tuples element by element, which requires that the children agree on
   how many elements they emit.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorMulti) :: nodePropertyExtractorDustAttenuation
     !!{RST
     A property extractor which applies dust attenuation to the luminosities of other extractors.
     !!}
     private
     class  (dustAttenuationClass), pointer :: dustAttenuation_   => null()
     logical                                :: outputUnattenuated          , outputSum , &
          &                                       outputSumOnly
     type   (varying_string      )          :: sumName
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
     procedure :: metaData           => dustAttenuationMetaData
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
    logical                                                      :: outputUnattenuated, outputSum, &
         &                                                                outputSumOnly
    type   (varying_string                      )                :: sumName

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

      The sum is formed *elementwise*, so summing spectra adds them wavelength by wavelength, and summing tuples adds
      corresponding elements. Children must therefore emit the same number of elements as one another, and it is an
      error if they do not.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputSumOnly</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, *only* the summed property is emitted, and the per-child properties from which it is formed are
      suppressed. This is what a consumer expecting a single property needs---an output analysis, for example, which
      requires one value per galaxy---since dust has to be applied to each component separately and the total
      recovered afterwards. Setting this implies ``[outputSum]``, and is incompatible with ``[outputUnattenuated]``.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>sumName</name>
      <defaultValue>var_str('luminosityTotal')</defaultValue>
      <description>
      The name given to the summed property emitted when ``[outputSum]`` is set. It must be set explicitly where more
      than one dust attenuation extractor using the same attenuator emits a sum, since the two would otherwise emit
      properties of the same name.
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
    self%outputSumOnly      =  outputSumOnly
    self%sumName            =  sumName
    call dustAttenuationValidate(self)
    !![
    <objectDestructor name="dustAttenuation_"/>
    !!]
    return
  end function dustAttenuationConstructorParameters

  function dustAttenuationConstructorInternal(dustAttenuation_,outputUnattenuated,outputSum,outputSumOnly,sumName,extractors) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorDustAttenuation` property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorDustAttenuation)                        :: self
    class  (dustAttenuationClass                ), intent(in   ), target :: dustAttenuation_
    logical                                      , intent(in   )         :: outputUnattenuated, outputSum, &
         &                                                                  outputSumOnly
    type   (varying_string                      ), intent(in   )         :: sumName
    type   (multiExtractorList                  ), intent(in   )         :: extractors
    !![
    <constructorAssign variables="outputUnattenuated, outputSum, outputSumOnly, sumName, *dustAttenuation_"/>
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

    if (self%outputSumOnly) then
       if (self%outputUnattenuated)                                                                                 &
            & call Error_Report(                                                                                    &
            &                   '[outputUnattenuated] and [outputSumOnly] are contradictory: the first asks for the'// &
            &                   ' per-child properties to be emitted, the second for them to be suppressed'         // &
            &                   {introspection:location}                                                             &
            &                  )
       ! Emitting only the sum requires that the sum be formed.
       self%outputSum=.true.
    end if
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
       class is (nodePropertyExtractorScalar  )
          ! Supported.
       class is (nodePropertyExtractorTuple   )
          ! Supported.
          class is (nodePropertyExtractorArray)
             ! Supported.
       class default
          call Error_Report(                                                                                   &
               &            'property extractor "'                                                          // &
               &            extractor_%extractor_%objectType()                                              // &
               &            '" produces properties of a rank which dust attenuation does not support - only'// &
               &            ' scalar, tuple, and array extractors may be attenuated'                        // &
               &            {introspection:location}                                                           &
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
    class is (nodePropertyExtractorArray )
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
    ! Where only the sum is emitted, that is the entirety of the output.
    if (self%outputSumOnly) then
       elementCount=1
       return
    end if
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
    type            (polyRankDouble                      )                         , allocatable, dimension(:  ) :: properties
    class           (nodePropertyExtractorDustAttenuation), intent(inout)                                        :: self
    type            (treeNode                            ), intent(inout)                                        :: node
    double precision                                      , intent(in   )                                        :: time
    type            (multiCounter                        ), intent(inout), optional                              :: instance
    integer                                               , intent(  out), optional, allocatable, dimension(:  ) :: ranks
    type            (multiExtractorList                  ), pointer                                              :: extractor_
    double precision                                                               , allocatable, dimension(:  ) :: valuesAttenuated, valuesRaw , &
         &                                                                                                          sumAttenuated
    double precision                                                               , allocatable, dimension(:,:) :: valuesRank1
    integer                                                                                                      :: offset          , countChild, &
         &                                                                                                          i               , countRows , &
         &                                                                                                          lengthElement

    allocate(properties(self%elementCount(elementTypeDouble,time)))
    if (present(ranks)) then
       ranks=self%ranks(elementTypeDouble,time)
    end if
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       if (self%outputUnattenuated) then
          select type (child => extractor_%extractor_)
          class is (nodePropertyExtractorScalar)
             allocate(valuesRaw(1))
             valuesRaw(1)=child%extract(node,     instance)
             do i=1,countChild
                properties(offset+i)=polyRankDouble(valuesRaw(i))
             end do
             deallocate(valuesRaw)
          class is (nodePropertyExtractorTuple )
             valuesRaw   =child%extract(node,time,instance)
             do i=1,countChild
                properties(offset+i)=polyRankDouble(valuesRaw(i))
             end do
             deallocate(valuesRaw)
          class is (nodePropertyExtractorArray )
             valuesRank1 =child%extract(node,time,instance)
             do i=1,countChild
                properties(offset+i)=polyRankDouble(valuesRank1(:,i))
             end do
             deallocate(valuesRank1)
          end select
          if (.not.self%outputSumOnly) offset=offset+countChild
       end if
       valuesAttenuated=self%attenuate(extractor_%extractor_,node,time,instance)
       ! The child returns one value per output *element*. For a scalar or tuple child that is one per property; for
       ! an array child, whose properties are themselves arrays, it is `size` values per property.
       countRows=1
       select type (child => extractor_%extractor_)
       class is (nodePropertyExtractorArray)
          countRows=int(child%size(time))
       end select
       if (size(valuesAttenuated) /= countChild*countRows)                                     &
            & call Error_Report(                                                               &
            &                   'the number of attenuated values returned by "'             // &
            &                   extractor_%extractor_%objectType()                          // &
            &                   '" does not match the number of output elements it declares'// &
            &                   {introspection:location}                                       &
            &                  )
       ! `recompose` returns the child's output elements flattened in column-major order, so for an array child the
       ! i-th property occupies the i-th contiguous block of `countRows` entries.
       if (.not.self%outputSumOnly) then
          do i=1,countChild
             if (countRows == 1) then
                properties(offset+i)=polyRankDouble(valuesAttenuated( i                           ))
             else
                properties(offset+i)=polyRankDouble(valuesAttenuated((i-1)*countRows+1:i*countRows))
             end if
          end do
       end if
       ! Accumulate the sum over children, if requested. Children must agree on the length of their properties for the
       ! sum to mean anything.
       if (self%outputSum) then
          lengthElement=size(valuesAttenuated)
          if (.not.allocated(sumAttenuated)) then
             allocate(sumAttenuated(lengthElement))
             sumAttenuated=0.0d0
          end if
          if (size(sumAttenuated) /= lengthElement)                                                             &
               & call Error_Report(                                                                             &
               &                   'children of this extractor emit properties of differing length, so their'// &
               &                   ' sum is not defined - unset [outputSum], or wrap children which match'   // &
               &                   {introspection:location}                                                     &
               &                  )
          sumAttenuated=sumAttenuated+valuesAttenuated
       end if
       if (.not.self%outputSumOnly) offset=offset+countChild
       deallocate(valuesAttenuated)
       extractor_ => extractor_%next
    end do
    if (self%outputSum) then
       if (.not.allocated(sumAttenuated)) allocate(sumAttenuated(0))
       offset=offset+1
       if (size(sumAttenuated) == 1) then
          properties(offset)=polyRankDouble(sumAttenuated(1))
       else
          properties(offset)=polyRankDouble(sumAttenuated   )
       end if
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
          if (.not.self%outputSumOnly) offset=offset+countChild
       end if
       if (.not.self%outputSumOnly) then
          do i=1,countChild
             names(offset+i)=namesChild(i)//suffix
          end do
       end if
       if (.not.self%outputSumOnly) offset=offset+countChild
       deallocate(namesChild)
       extractor_ => extractor_%next
    end do
    if (self%outputSum) then
       offset       =offset+1
       names(offset)=self%sumName//suffix
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
    class is (nodePropertyExtractorArray )
       ! Note the argument order: the array class takes the output array first, unlike the tuple class.
       call extractor_%names(names,time)
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
       class is (nodePropertyExtractorArray )
          call child%descriptions(descriptionsChild,time)
       end select
       if (self%outputUnattenuated) then
          do i=1,countChild
             descriptions(offset+i)=descriptionsChild(i)
          end do
          if (.not.self%outputSumOnly) offset=offset+countChild
       end if
       if (.not.self%outputSumOnly) then
          do i=1,countChild
             descriptions(offset+i)=descriptionsChild(i)//suffix
          end do
       end if
       if (.not.self%outputSumOnly) offset=offset+countChild
       deallocate(descriptionsChild)
       extractor_ => extractor_%next
    end do
    if (self%outputSum) then
       offset            =offset+1
       descriptions(offset)="Sum over all child extractors of "//self%sumName//suffix
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
    double precision                                                                    :: unitsSum
    logical                                                                             :: unitsSumSet

    allocate(unitsInSI(self%elementCount(elementType,time)))
    if (elementType /= elementTypeDouble) return
    offset      =  0
    unitsSum    =  1.0d0
    unitsSumSet =  .false.
    extractor_  => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       select type (child => extractor_%extractor_)
       class is (nodePropertyExtractorScalar)
          allocate(unitsChild(1))
          unitsChild(1)=child%unitsInSI()
       class is (nodePropertyExtractorTuple )
          unitsChild   =child%unitsInSI(time)
       class is (nodePropertyExtractorArray )
          unitsChild   =child%unitsInSI(time)
       end select
       if (.not.self%outputSumOnly) then
          if (self%outputUnattenuated) then
             unitsInSI(offset+1:offset+countChild)=unitsChild
             offset                               =offset+countChild
          end if
          unitsInSI(offset+1:offset+countChild)=unitsChild
          offset                               =offset+countChild
       end if
       ! Remember the units of the first child. The summed property carries them, and can not read them back from an
       ! emitted property: where only the sum is emitted, there is no earlier property to read.
       if (.not.unitsSumSet .and. size(unitsChild) > 0) then
          unitsSum   =unitsChild(1)
          unitsSumSet=.true.
       end if
       deallocate(unitsChild)
       extractor_ => extractor_%next
    end do
    ! The summed property takes the units of its children, all of which are luminosities in the same units.
    if (self%outputSum) then
       offset           =offset+1
       unitsInSI(offset)=unitsSum
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

    type            (unitType                            )                              :: unitsSum
    logical                                                                             :: unitsSumSet
    allocate(units(self%elementCount(elementType,time)))
    if (elementType /= elementTypeDouble) return
    offset      =  0
    unitsSumSet =  .false.
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
       class is (nodePropertyExtractorArray )
          unitsChild   =child%units(time)
       end select
       if (self%outputUnattenuated) then
          do i=1,countChild
             units(offset+i)=unitsChild(i)
          end do
          if (.not.self%outputSumOnly) offset=offset+countChild
       end if
       if (.not.self%outputSumOnly) then
          do i=1,countChild
             units(offset+i)=unitsChild(i)
          end do
       end if
       if (.not.self%outputSumOnly) offset=offset+countChild
       ! Remember the units of the first child. The summed property carries them, and can not read them back from an
       ! emitted property: where only the sum is emitted, there is no earlier property to read.
       if (.not.unitsSumSet .and. size(unitsChild) > 0) then
          unitsSum   =unitsChild(1)
          unitsSumSet=.true.
       end if
       deallocate(unitsChild)
       extractor_ => extractor_%next
    end do
    if (self%outputSum) then
       offset      =offset+1
       units(offset)=unitsSum
    end if
    return

  end function dustAttenuationUnits

  function dustAttenuationRanks(self,elementType,time) result(ranks)
    !!{RST
    Return the ranks of the properties emitted: rank one for the properties of an array child, and rank zero
    otherwise. Attenuation does not change the shape of a property, so an attenuated property has the rank of the
    property it attenuates.
    !!}
    implicit none
    integer                                               , allocatable  , dimension(:) :: ranks
    class           (nodePropertyExtractorDustAttenuation), intent(inout)               :: self
    type            (enumerationElementTypeType          ), intent(in   )               :: elementType
    double precision                                      , intent(in   )               :: time
    type            (multiExtractorList                  ), pointer                     :: extractor_
    integer                                                                             :: offset    , countChild, &
         &                                                                                 rankChild , countRows , &
         &                                                                                 lengthSum

    allocate(ranks(self%elementCount(elementType,time)))
    ranks=0
    if (elementType /= elementTypeDouble) return
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       rankChild =0
       select type (child => extractor_%extractor_)
       class is (nodePropertyExtractorArray)
          rankChild=1
       end select
       if (.not.self%outputSumOnly) then
          if (self%outputUnattenuated) then
             ranks(offset+1:offset+countChild)=rankChild
             offset                           =offset+countChild
          end if
          ranks(offset+1:offset+countChild)=rankChild
          offset                           =offset+countChild
       end if
       extractor_ => extractor_%next
    end do
    ! The summed property is formed elementwise, so it has one value per output element of a child. Its rank must
    ! therefore be derived exactly as `extractDouble` derives its shape -- from that element count -- and not from
    ! whether a child happens to be of the array class: a tuple child emitting several elements yields a rank one sum,
    ! while an array child of unit length yields a rank zero one.
    if (self%outputSum) then
       lengthSum  =  0
       extractor_ => self%extractors
       if (associated(extractor_)) then
          countChild=self%countChildElements(extractor_%extractor_,time)
          countRows =1
          select type (child => extractor_%extractor_)
          class is (nodePropertyExtractorArray)
             countRows=int(child%size(time))
          end select
          lengthSum=countChild*countRows
       end if
       offset=offset+1
       if (lengthSum > 1) then
          ranks(offset)=1
       else
          ranks(offset)=0
       end if
    end if
    return
  end function dustAttenuationRanks

  subroutine dustAttenuationMetaData(self,node,elementType,time,iProperty,metaDataRank0,metaDataRank1)
    !!{RST
    Populate metadata for the properties emitted.

    The inherited implementation maps a property index onto a child by walking the child list, which is wrong here:
    this extractor emits its children's properties twice when ``[outputUnattenuated]`` is set, not at all when
    ``[outputSumOnly]`` is set, and appends a summed property belonging to no child. Attaching a child's metadata to
    that sum would be actively misleading---the sum of two emission lines is not emitted at either line's
    wavelength---so the summed property is given none.
    !!}
    implicit none
    class           (nodePropertyExtractorDustAttenuation), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    type            (enumerationElementTypeType          ), intent(in   ) :: elementType
    double precision                                      , intent(in   ) :: time
    integer                                               , intent(in   ) :: iProperty
    type            (doubleDictionary                    ), intent(inout) :: metaDataRank0
    type            (rank1DoubleDictionary               ), intent(inout) :: metaDataRank1
    type            (multiExtractorList                  ), pointer       :: extractor_
    integer                                                               :: offset    , countChild

    if (elementType /= elementTypeDouble) return
    ! Where only the sum is emitted, there is no child property to describe.
    if (self%outputSumOnly                ) return
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       countChild=self%countChildElements(extractor_%extractor_,time)
       ! The unattenuated copy of a property, where emitted, carries the same metadata as the attenuated one.
       if (self%outputUnattenuated) then
          if (iProperty > offset .and. iProperty <= offset+countChild) then
             call dustAttenuationChildMetaData(extractor_%extractor_,node,iProperty-offset,metaDataRank0,metaDataRank1)
             return
          end if
          offset=offset+countChild
       end if
       if    (iProperty > offset .and. iProperty <= offset+countChild) then
          call dustAttenuationChildMetaData(extractor_%extractor_,node,iProperty-offset,metaDataRank0,metaDataRank1)
          return
       end if
       offset     =  offset+countChild
       extractor_ => extractor_%next
    end do
    ! Anything beyond the children is the summed property, which describes no single child and so carries no metadata.
    return
  end subroutine dustAttenuationMetaData

  subroutine dustAttenuationChildMetaData(extractor_,node,indexProperty,metaDataRank0,metaDataRank1)
    !!{RST
    Populate metadata from a child extractor for one of its properties.
    !!}
    implicit none
    class  (nodePropertyExtractorClass ), intent(inout) :: extractor_
    type   (treeNode                   ), intent(inout) :: node
    integer                             , intent(in   ) :: indexProperty
    type   (doubleDictionary           ), intent(inout) :: metaDataRank0
    type   (rank1DoubleDictionary      ), intent(inout) :: metaDataRank1

    select type (extractor_)
    class is (nodePropertyExtractorScalar)
       call extractor_%metaData(node,              metaDataRank0,metaDataRank1)
    class is (nodePropertyExtractorTuple )
       call extractor_%metaData(node,indexProperty,metaDataRank0,metaDataRank1)
    class is (nodePropertyExtractorArray )
       call extractor_%metaData(node,indexProperty,metaDataRank0,metaDataRank1)
    end select
    return
  end subroutine dustAttenuationChildMetaData
