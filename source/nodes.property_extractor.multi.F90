!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  Implements a multi node property extractor class.
  !!}

  use :: Hashes, only : doubleHash, rank1DoubleHash

  type, public :: multiExtractorList
     class(nodePropertyExtractorClass), pointer :: extractor_ => null()
     type (multiExtractorList        ), pointer :: next       => null()
  end type multiExtractorList

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMulti">
   <description>A multi output extractor property extractor class.</description>
   <linkedList type="multiExtractorList" variable="extractors" next="next" object="extractor_" objectType="nodePropertyExtractorClass"/>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorClass) :: nodePropertyExtractorMulti
     !!{
     A multi property extractor output extractor class, which concatenates properties from any number of other property
     extractors.
     !!}
     private
     type(multiExtractorList), pointer :: extractors => null()
   contains
     !![
     <methods>
       <method description="Return a description of the columns."                                        method="columnDescriptions"/>
       <method description="Return the number of properties in the tuple."                               method="elementCount"      />
       <method description="Extract the double properties from the given {\normalfont \ttfamily node}."  method="extractDouble"     />
       <method description="Extract the integer properties from the given {\normalfont \ttfamily node}." method="extractInteger"    />
       <method description="Return the names of the properties extracted."                               method="names"             />
       <method description="Return descriptions of the properties extracted."                            method="descriptions"      />
       <method description="Return the units of the properties extracted in the SI system."              method="unitsInSI"         />
       <method description="Return the ranks of the properties extracted."                               method="ranks"             />
       <method description="Populate a hash with meta-data for the property."                            method="metaData"          />
     </methods>
     !!]
     final     ::                       multiDestructor
     procedure :: columnDescriptions => multiColumnDescriptions
     procedure :: elementCount       => multiElementCount
     procedure :: extractDouble      => multiExtractDouble
     procedure :: extractInteger     => multiExtractInteger
     procedure :: names              => multiNames
     procedure :: descriptions       => multiDescriptions
     procedure :: unitsInSI          => multiUnitsInSI
     procedure :: ranks              => multiRanks
     procedure :: addInstances       => multiAddInstances
     procedure :: metaData           => multiMetaData
  end type nodePropertyExtractorMulti

  interface nodePropertyExtractorMulti
     !!{
     Constructors for the ``multi'' output extractor class.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface nodePropertyExtractorMulti

  !![
  <enumeration>
   <name>elementType</name>
   <description>Enumeration of extracted property element types.</description>
   <visibility>public</visibility>
   <entry label="integer"/>
   <entry label="double" />
  </enumeration>
  !!]

contains

  function multiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``multi'' output extractor property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorMulti)                :: self
    type   (inputParameters           ), intent(inout) :: parameters
    type   (multiExtractorList        ), pointer       :: extractor_
    integer                                            :: i

    self      %extractors => null()
    extractor_            => null()
    do i=1,parameters%copiesCount('nodePropertyExtractor',zeroIfNotPresent=.true.)
       if (associated(extractor_)) then
          allocate(extractor_%next)
          extractor_ => extractor_%next
       else
          allocate(self%extractors)
          extractor_ => self%extractors
       end if
       !![
       <objectBuilder class="nodePropertyExtractor" name="extractor_%extractor_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="nodePropertyExtractor"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(extractors) result(self)
    !!{
    Internal constructor for the ``multi'' output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMulti)                         :: self
    type(multiExtractorList        ), target , intent(in   ) :: extractors
    type(multiExtractorList        ), pointer                :: extractor_

    self      %extractors => extractors
    extractor_            => extractors
    do while (associated(extractor_))
       !![
       <referenceCountIncrement owner="extractor_" object="extractor_"/>
       !!]
       extractor_ => extractor_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily multi} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMulti), intent(inout) :: self
    type(multiExtractorList        ), pointer       :: extractor_, extractorNext

    if (associated(self%extractors)) then
       extractor_ => self%extractors
       do while (associated(extractor_))
          extractorNext => extractor_%next
          !![
          <objectDestructor name="extractor_%extractor_"/>
          !!]
          deallocate(extractor_)
          extractor_ => extractorNext
       end do
    end if
    return
  end subroutine multiDestructor

  integer function multiElementCount(self,elementType,time)
    !!{
    Return the number of elements in the multiple property extractors.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorMulti), intent(inout) :: self
    type            (enumerationElementTypeType), intent(in   ) :: elementType
    double precision                            , intent(in   ) :: time
    type            (multiExtractorList        ), pointer       :: extractor_

    multiElementCount =  0
    extractor_        => self%extractors
    do while (associated(extractor_))
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) multiElementCount=multiElementCount+1
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) multiElementCount=multiElementCount+extractor_%elementCount(            time)
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) multiElementCount=multiElementCount+1
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) multiElementCount=multiElementCount+extractor_%elementCount(            time)
       class is (nodePropertyExtractorArray        )
          if (elementType == elementTypeDouble ) multiElementCount=multiElementCount+extractor_%elementCount(            time)
       class is (nodePropertyExtractorList         )
          if (elementType == elementTypeDouble ) multiElementCount=multiElementCount+extractor_%elementCount(                )
       class is (nodePropertyExtractorIntegerList  )
          if (elementType == elementTypeInteger) multiElementCount=multiElementCount+extractor_%elementCount(                )
       class is (nodePropertyExtractorList2D       )
          if (elementType == elementTypeDouble ) multiElementCount=multiElementCount+extractor_%elementCount(                )
       class is (nodePropertyExtractorMulti        )
          multiElementCount                                       =multiElementCount+extractor_%elementCount(elementType,time)
       class default
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       extractor_ => extractor_%next
    end do
    return
  end function multiElementCount

  function multiExtractDouble(self,node,time,instance,ranks)
    !!{
    Implement a multi output extractor.
    !!}
    use :: Error     , only : Error_Report
    use :: Poly_Ranks, only : polyRankDouble
    implicit none
    type            (polyRankDouble            )                         , allocatable, dimension(:    ) :: multiExtractDouble
    class           (nodePropertyExtractorMulti), intent(inout)                                          :: self
    type            (treeNode                  ), intent(inout)                                          :: node
    double precision                            , intent(in   )                                          :: time
    type            (multiCounter              ), intent(inout), optional                                :: instance
    integer                                     , intent(  out), optional, allocatable, dimension(:    ) :: ranks
    type            (multiExtractorList        ), pointer                                                :: extractor_
    double precision                                                     , allocatable, dimension(:    ) :: rank0
    double precision                                                     , allocatable, dimension(:,:  ) :: rank1
    double precision                                                     , allocatable, dimension(:,:,:) :: rank2
    integer                                                              , allocatable, dimension(:    ) :: ranks_
    integer                                                                                              :: offset            , elementCount, &
         &                                                                                                  i

    elementCount=self%elementCount(elementTypeDouble,time)
    allocate(multiExtractDouble(elementCount))
    if (present(ranks)) allocate(ranks(self%elementCount(elementTypeDouble,time)))
    if (elementCount == 0) return
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          elementCount=1
          multiExtractDouble(offset+1)=polyRankDouble(extractor_%extract(node,instance))
          if (present(ranks)) ranks(offset+1)=0
       class is (nodePropertyExtractorTuple        )
          elementCount=extractor_%elementCount(time)
          rank0=extractor_%extract(node,time,instance)
          if (elementCount > 0) then
             do i=1,elementCount
                multiExtractDouble(offset+i)=polyRankDouble(rank0(i))
                if (present(ranks)) ranks(offset+i)=0
             end do
             deallocate(rank0)
          end if
       class is (nodePropertyExtractorArray        )
          elementCount=extractor_%elementCount(time)
          if (elementCount > 0) then
             rank1=extractor_%extract(node,time,instance)
             do i=1,elementCount
                multiExtractDouble(offset+i)=polyRankDouble(rank1(:,i))
                if (present(ranks)) ranks(offset+i)=1
             end do
             deallocate(rank1)
          end if
       class is (nodePropertyExtractorList         )
          elementCount=extractor_%elementCount()
          if (elementCount > 0) then
             rank1=extractor_%extract(node     ,instance)
             do i=1,elementCount
                multiExtractDouble(offset+i)=polyRankDouble(rank1(:,i))
                if (present(ranks)) ranks(offset+i)=-1
             end do
             deallocate(rank1)
          end if
       class is (nodePropertyExtractorList2D       )
          elementCount=extractor_%elementCount()
          if (elementCount > 0) then
             rank2=extractor_%extract(node     ,instance)
             do i=1,elementCount
                multiExtractDouble(offset+i)=polyRankDouble(rank2(:,:,i))
                if (present(ranks)) ranks(offset+i)=-2
             end do
             deallocate(rank2)
          end if
       class is (nodePropertyExtractorMulti        )
          elementCount=extractor_%elementCount(elementTypeDouble,time)
          if (elementCount > 0) then
             !![
	     <conditionalCall>
	       <call>multiExtractDouble(offset+1:offset+elementCount)=extractor_%extractDouble(node,time,instance{conditions})</call>
	       <argument name="ranks" value="ranks_" condition="present(ranks)"/>
	     </conditionalCall>
             !!]
             if (present(ranks)) then
                ranks(offset+1:offset+elementCount)=ranks_
                deallocate(ranks_)
             end if
          end if
       class is (nodePropertyExtractorIntegerScalar)
          elementCount=0
       class is (nodePropertyExtractorIntegerTuple )
          elementCount=0
       class is (nodePropertyExtractorIntegerList  )
          elementCount=0
       class default
          elementCount=0
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiExtractDouble

  function multiExtractInteger(self,node,time,instance,ranks)
    !!{
    Implement a multi output extractor.
    !!}
    use :: Error     , only : Error_Report
    use :: Poly_Ranks, only : polyRankInteger
    implicit none
    type            (polyRankInteger           )                         , allocatable, dimension(:  ) :: multiExtractInteger
    class           (nodePropertyExtractorMulti), intent(inout)                                        :: self
    type            (treeNode                  ), intent(inout)                                        :: node
    double precision                            , intent(in   )                                        :: time
    type            (multiCounter              ), intent(inout), optional                              :: instance
    integer                                     , intent(  out), optional, allocatable, dimension(:  ) :: ranks
    type            (multiExtractorList        ), pointer                                              :: extractor_
    integer         (kind_int8                 )                         , allocatable, dimension(:  ) :: rank0
    integer         (kind_int8                 )                         , allocatable, dimension(:,:) :: rank1
    integer                                                              , allocatable, dimension(:  ) :: ranks_
    integer                                                                                            :: offset             , elementCount, &
         &                                                                                                i

    elementCount=self%elementCount(elementTypeInteger,time)
    allocate(multiExtractInteger(elementCount))
    if (present(ranks)) allocate(ranks(self%elementCount(elementTypeInteger,time)))
    if (elementCount == 0) return
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          elementCount=0
       class is (nodePropertyExtractorTuple        )
          elementCount=0
       class is (nodePropertyExtractorArray        )
          elementCount=0
       class is (nodePropertyExtractorList         )
          elementCount=0
       class is (nodePropertyExtractorList2D       )
          elementCount=0
       class is (nodePropertyExtractorIntegerScalar)
          elementCount=1
          multiExtractInteger(offset+1)=polyRankInteger(extractor_%extract(node,time,instance))
          if (present(ranks)) ranks(offset+1)=0
       class is (nodePropertyExtractorIntegerTuple )
          elementCount=extractor_%elementCount(time)
          if (elementCount > 0) then
             rank0=extractor_%extract(node,time,instance)
             do i=1,elementCount
                multiExtractInteger(offset+i)=polyRankInteger(rank0(i))
                if (present(ranks)) ranks(offset+i)=0
             end do
             deallocate(rank0)
          end if
       class is (nodePropertyExtractorIntegerList  )
          elementCount=extractor_%elementCount()
          if (elementCount > 0) then
             rank1=extractor_%extract(node     ,instance)
             do i=1,elementCount
                multiExtractInteger(offset+i)=polyRankInteger(rank1(:,i))
                if (present(ranks)) ranks(offset+i)=-1
             end do
             deallocate(rank1)
          end if
       class is (nodePropertyExtractorMulti        )
          elementCount=extractor_%elementCount(elementTypeInteger,time)
          if (elementCount > 0) then
             !![
	     <conditionalCall>
	       <call>multiExtractInteger(offset+1:offset+elementCount)=extractor_%extractInteger(node,time,instance{conditions})</call>
	       <argument name="ranks" value="ranks_" condition="present(ranks)"/>
	     </conditionalCall>
             !!]
             if (present(ranks)) then
                ranks(offset+1:offset+elementCount)=ranks_
                deallocate(ranks_)
             end if
          end if
       class default
          elementCount=0
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiExtractInteger

  subroutine multiAddInstances(self,node,instance)
    !!{
    Implement adding of instances to a multi output extractor.
    !!}
    implicit none
    class(nodePropertyExtractorMulti), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    type (multiCounter              ), intent(inout) :: instance
    type (multiExtractorList        ), pointer       :: extractor_

    extractor_ => self%extractors
    do while (associated(extractor_))
       call extractor_%extractor_%addInstances(node,instance)
       extractor_ => extractor_%next
    end do
    return
  end subroutine multiAddInstances

  subroutine multiNames(self,elementType,time,names)
    !!{
    Return the names of the multiple properties.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorMulti), intent(inout)                             :: self
    type            (enumerationElementTypeType), intent(in   )                             :: elementType
    double precision                            , intent(in   )                             :: time
    type            (varying_string            ), intent(inout), dimension(:) , allocatable :: names
    type            (varying_string            )               , dimension(:) , allocatable :: namesTmp
    type            (multiExtractorList        ), pointer                                   :: extractor_
    integer                                                                                 :: offset    , elementCount
    
    allocate(names(self%elementCount(elementType,time)))
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       elementCount=0
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) then
             elementCount                              =1
             names       (offset+1:offset+elementCount)=extractor_%name (    )
          end if
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             if (elementCount > 0) then
                call extractor_%names(time,namesTmp)
                names(offset+1:offset+elementCount)=namesTmp
                deallocate(namesTmp)
             end if
          end if
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) then
             elementCount                              =1
             names       (offset+1:offset+elementCount)=extractor_%name (    )
          end if
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount(time)
             if (elementCount > 0) then
                call extractor_%names(time,namesTmp)
                names(offset+1:offset+elementCount)=namesTmp
                deallocate(namesTmp)
             end if
          end if
       class is (nodePropertyExtractorArray        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time) 
             if (elementCount > 0) then
                call extractor_%names(namesTmp,time)
                names(offset+1:offset+elementCount)=namesTmp
                deallocate(namesTmp)
             end if
          end if
       class is (nodePropertyExtractorList         )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             if (elementCount > 0) then
                call extractor_%names(namesTmp     )
                names(offset+1:offset+elementCount)=namesTmp
                deallocate(namesTmp)
             end if
          end if
       class is (nodePropertyExtractorIntegerList  )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount()
             if (elementCount > 0) then
                call extractor_%names(namesTmp     )
                names(offset+1:offset+elementCount)=namesTmp
                deallocate(namesTmp)
             end if
          end if
       class is (nodePropertyExtractorList2D       )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             if (elementCount > 0) then
                call extractor_%names(namesTmp     )
                names(offset+1:offset+elementCount)=namesTmp
                deallocate(namesTmp)
             end if
          end if
       class is (nodePropertyExtractorMulti        )
          elementCount=extractor_%elementCount(elementType,time)
          if (elementCount > 0) then
             call extractor_%names(elementType,time,namesTmp)
             names(offset+1:offset+elementCount)=namesTmp
             deallocate(namesTmp)
          end if
       class default
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end subroutine multiNames

  subroutine multiColumnDescriptions(self,elementType,i,time,descriptions,values,valuesDescription,valuesUnitsInSI)
    !!{
    Return column descriptions of the multiple properties.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorMulti), intent(inout)                             :: self
    type            (enumerationElementTypeType), intent(in   )                             :: elementType
    integer                                     , intent(in   )                             :: i
    double precision                            , intent(in   )                             :: time
    type            (varying_string            ), intent(  out), dimension(:) , allocatable :: descriptions
    double precision                            , intent(  out), dimension(:) , allocatable :: values
    type            (varying_string            ), intent(  out)                             :: valuesDescription
    double precision                            , intent(  out)                             :: valuesUnitsInSI
    type            (multiExtractorList        ), pointer                                   :: extractor_
    integer                                                                                 :: elementCount     , offset

    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       elementCount=0
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) then
             elementCount=1
             if (offset+elementCount >= i) then
                allocate(descriptions(0))
                return
             end if
          end if
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             if (offset+elementCount >= i) then
                allocate(descriptions(0))
                return
             end if
          end if
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) then
             elementCount=1
             if (offset+elementCount >= i) then
                allocate(descriptions(0))
                return
             end if
          end if
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount(time)
             if (offset+elementCount >= i) then
                allocate(descriptions(0))
                return
             end if
          end if
       class is (nodePropertyExtractorArray        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             if (offset+elementCount >= i) then
                call extractor_%columnDescriptions(descriptions,values,valuesDescription,valuesUnitsInSI,time)
                return
             end if
          end if
       class is (nodePropertyExtractorList         )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             if (offset+elementCount >= i) then
                allocate(descriptions(0))
                return
             end if
          end if
       class is (nodePropertyExtractorIntegerList  )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount()
             if (offset+elementCount >= i) then
                allocate(descriptions(0))
                return
             end if
          end if
       class is (nodePropertyExtractorList2D       )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             if (offset+elementCount >= i) then
                allocate(descriptions(0))
                return
             end if
          end if
       class is (nodePropertyExtractorMulti        )
          elementCount=extractor_%elementCount(elementType,time)
          if (offset+elementCount >= i) then
             call extractor_%columnDescriptions(elementType,i-offset,time,descriptions,values,valuesDescription,valuesUnitsInSI)
             return
          end if
       class default
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end subroutine multiColumnDescriptions

  subroutine multiDescriptions(self,elementType,time,descriptions)
    !!{
    Return the descriptions of the multiple properties.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorMulti), intent(inout)                             :: self
    type            (enumerationElementTypeType), intent(in   )                             :: elementType
    double precision                            , intent(in   )                             :: time
    type            (varying_string            ), intent(inout), dimension(:) , allocatable :: descriptions
    type            (varying_string            )               , dimension(:) , allocatable :: descriptionsTmp
    type            (multiExtractorList        ), pointer                                   :: extractor_
    integer                                                                                 :: offset         , elementCount

    allocate(descriptions(self%elementCount(elementType,time)))
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       elementCount=0
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) then
             elementCount=1
             descriptions       (offset+1:offset+elementCount)=extractor_%description (    )
          end if
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             if (elementCount > 0) then
                call extractor_%descriptions(time,descriptionsTmp)
                descriptions(offset+1:offset+elementCount)=descriptionsTmp
                deallocate(descriptionsTmp)
             end if
          end if
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) then
             elementCount=1
             descriptions       (offset+1:offset+elementCount)=extractor_%description (    )
          end if
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount(time)
             if (elementCount > 0) then
                call extractor_%descriptions(time,descriptionsTmp)
                descriptions(offset+1:offset+elementCount)=descriptionsTmp
                deallocate(descriptionsTmp)
             end if
          end if
       class is (nodePropertyExtractorArray        )
          if (elementType == elementTypeDouble) then
             elementCount=extractor_%elementCount(time)
             if (elementCount > 0) then
                call extractor_%descriptions(descriptionsTmp,time)
                descriptions(offset+1:offset+elementCount)=descriptionsTmp
                deallocate(descriptionsTmp)
             end if
          end if
       class is (nodePropertyExtractorList         )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             if (elementCount > 0) then
                call extractor_%descriptions(descriptionsTmp)
                descriptions(offset+1:offset+elementCount)=descriptionsTmp
                deallocate(descriptionsTmp)
             end if
          end if
        class is (nodePropertyExtractorIntegerList )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount()
             if (elementCount > 0) then
                call extractor_%descriptions(descriptionsTmp)
                descriptions(offset+1:offset+elementCount)=descriptionsTmp
                deallocate(descriptionsTmp)
             end if
          end if
       class is (nodePropertyExtractorList2D       )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             if (elementCount > 0) then
                call extractor_%descriptions(descriptionsTmp)
                descriptions(offset+1:offset+elementCount)=descriptionsTmp
                deallocate(descriptionsTmp)
             end if
          end if
       class is (nodePropertyExtractorMulti        )
          elementCount=extractor_%elementCount(elementType,time)
          if (elementCount > 0) then
             call extractor_%descriptions(elementType,time,descriptionsTmp)
             descriptions(offset+1:offset+elementCount)=descriptionsTmp
             deallocate(descriptionsTmp)
          end if
       class default
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end subroutine multiDescriptions

  function multiUnitsInSI(self,elementType,time)
    !!{
    Return the units of the multiple properties in the SI system.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                            , dimension(:) , allocatable :: multiUnitsInSI
    class           (nodePropertyExtractorMulti), intent(inout)              :: self
    type            (enumerationElementTypeType), intent(in   )              :: elementType
    double precision                            , intent(in   )              :: time
    type            (multiExtractorList        ), pointer                    :: extractor_
    integer                                                                  :: offset        , elementCount

    allocate(multiUnitsInSI(self%elementCount(elementType,time)))
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       elementCount=0
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) then
             elementCount=1
             multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(    )
          end if
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(time)
          end if
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) then
             elementCount=1
             multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(    )
          end if
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount(time)
             multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(time)
          end if
       class is (nodePropertyExtractorArray        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(time)
          end if
       class is (nodePropertyExtractorList         )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(    )
          end if
       class is (nodePropertyExtractorIntegerList  )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount()
             multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(    )
          end if
       class is (nodePropertyExtractorList2D       )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(    )
          end if
       class is (nodePropertyExtractorMulti        )
          elementCount=extractor_%elementCount(elementType,time)
          if (elementCount > 0)                                                                      &
               & multiUnitsInSI(offset+1:offset+elementCount)=extractor_%unitsInSI(elementType,time)
       class default
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiUnitsInSI

  function multiRanks(self,elementType,time)
    !!{
    Return the ranks of the multiple properties. Negative values indicate variable length ranks.
    !!}
    use :: Error, only : Error_Report
    implicit none
    integer                                     , dimension(:) , allocatable :: multiRanks
    class           (nodePropertyExtractorMulti), intent(inout)              :: self
    type            (enumerationElementTypeType), intent(in   )              :: elementType
    double precision                            , intent(in   )              :: time
    type            (multiExtractorList        ), pointer                    :: extractor_
    integer                                                                  :: offset     , elementCount

    allocate(multiRanks(self%elementCount(elementType,time)))
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       elementCount=0
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) then
             elementCount=1
             multiRanks(offset+1:offset+elementCount)=0
          end if
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             multiRanks(offset+1:offset+elementCount)=0
          end if
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) then
             elementCount=1
             multiRanks(offset+1:offset+elementCount)=0
          end if
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount(time)
             multiRanks(offset+1:offset+elementCount)=0
          end if
       class is (nodePropertyExtractorArray        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             multiRanks(offset+1:offset+elementCount)=1
          end if
       class is (nodePropertyExtractorList        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             multiRanks(offset+1:offset+elementCount)=-1
          end if
       class is (nodePropertyExtractorIntegerList )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount()
             multiRanks(offset+1:offset+elementCount)=-1
          end if
       class is (nodePropertyExtractorList2D      )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             multiRanks(offset+1:offset+elementCount)=-2
          end if
       class is (nodePropertyExtractorMulti       )
          elementCount=extractor_%elementCount(elementType,time)
          if (elementCount > 0) multiRanks(offset+1:offset+elementCount)=extractor_%ranks(elementType,time)
       class default
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiRanks

  subroutine multiMetaData(self,node,elementType,time,iProperty,metaDataRank0,metaDataRank1)
    !!{
    Populate multiple property meta-data.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (nodePropertyExtractorMulti), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    type            (enumerationElementTypeType), intent(in   ) :: elementType
    double precision                            , intent(in   ) :: time
    integer                                     , intent(in   ) :: iProperty
    type            (doubleHash                ), intent(inout) :: metaDataRank0
    type            (rank1DoubleHash           ), intent(inout) :: metaDataRank1
    type            (multiExtractorList        ), pointer       :: extractor_
    integer                                                     :: offset       , elementCount

    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       elementCount=0
              select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) then
             elementCount=1
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node                                  ,metaDataRank0,metaDataRank1)
          end if
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node                 ,iProperty-offset,metaDataRank0,metaDataRank1)
          end if
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) then
             elementCount=1
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node                                  ,metaDataRank0,metaDataRank1)
          end if
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount(time)
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node                 ,iProperty-offset,metaDataRank0,metaDataRank1)
          end if
       class is (nodePropertyExtractorArray        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node                 ,iProperty-offset,metaDataRank0,metaDataRank1)
          end if
       class is (nodePropertyExtractorList         )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node                                  ,metaDataRank0,metaDataRank1)
          end if
       class is (nodePropertyExtractorIntegerList  )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount()
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node                                  ,metaDataRank0,metaDataRank1)
          end if
       class is (nodePropertyExtractorList2D       )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount()
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node            ,time                 ,metaDataRank0,metaDataRank1)
          end if
       class is (nodePropertyExtractorMulti        )
          elementCount=extractor_%elementCount(elementType,time)
          if (elementCount > 0) then
             if (offset+1 <= iProperty .and. offset+elementCount >= iProperty) call extractor_%metaData(node,elementType,time,iProperty-offset,metaDataRank0,metaDataRank1)
          end if
       class default
          call Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end subroutine multiMetaData
