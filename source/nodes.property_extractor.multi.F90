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

  !% Contains a module which implements a multi node property extractor class.

  type, public :: multiExtractorList
     class(nodePropertyExtractorClass), pointer :: extractor_
     type (multiExtractorList        ), pointer :: next       => null()
  end type multiExtractorList

  !# <nodePropertyExtractor name="nodePropertyExtractorMulti">
  !#  <description>A multi output extractor property extractor class.</description>
  !#  <deepCopy>
  !#   <linkedList type="multiExtractorList" variable="extractors" next="next" object="extractor_" objectType="nodePropertyExtractorClass"/>
  !#  </deepCopy>
  !# </nodePropertyExtractor>
  type, extends(nodePropertyExtractorClass) :: nodePropertyExtractorMulti
     !% A multi property extractor output extractor class, which concatenates properties from any number of other property
     !% extractors.
     private
     type(multiExtractorList), pointer :: extractors => null()
   contains
     !# <methods>
     !#   <method description="Return the number of properties in the tuple." method="elementCount" pass="yes" />
     !#   <method description="Extract the double properties from the given {\normalfont \ttfamily node}." method="extractDouble" pass="yes" />
     !#   <method description="Extract the integer properties from the given {\normalfont \ttfamily node}." method="extractInteger" pass="yes" />
     !#   <method description="Return the names of the properties extracted." method="names" pass="yes" />
     !#   <method description="Return descriptions of the properties extracted." method="descriptions" pass="yes" />
     !#   <method description="Return the units of the properties extracted in the SI system." method="unitsInSI" pass="yes" />
     !# </methods>
     final     ::                   multiDestructor
     procedure :: elementCount   => multiElementCount
     procedure :: extractDouble  => multiExtractDouble
     procedure :: extractInteger => multiExtractInteger
     procedure :: names          => multiNames
     procedure :: descriptions   => multiDescriptions
     procedure :: unitsInSI      => multiUnitsInSI
     procedure :: addInstances   => multiAddInstances
     procedure :: type           => multiType
  end type nodePropertyExtractorMulti

  interface nodePropertyExtractorMulti
     !% Constructors for the ``multi'' output extractor class.
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface nodePropertyExtractorMulti

  !# <enumeration>
  !#  <name>elementType</name>
  !#  <description>Enumeration of extracted property element types.</description>
  !#  <visibility>public</visibility>
  !#  <entry label="integer"/>
  !#  <entry label="double" />
  !# </enumeration>

contains

  function multiConstructorParameters(parameters) result(self)
    !% Constructor for the ``multi'' output extractor property extractor class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorMulti)                :: self
    type   (inputParameters           ), intent(inout) :: parameters
    type   (multiExtractorList        ), pointer       :: extractor_
    integer                                            :: i

    self      %extractors => null()
    extractor_            => null()
    do i=1,parameters%copiesCount('nodePropertyExtractorMethod',zeroIfNotPresent=.true.)
       if (associated(extractor_)) then
          allocate(extractor_%next)
          extractor_ => extractor_%next
       else
          allocate(self%extractors)
          extractor_ => self%extractors
       end if
       !# <objectBuilder class="nodePropertyExtractor" name="extractor_%extractor_" source="parameters" copy="i" />
    end do
    return
  end function multiConstructorParameters

  function multiConstructorInternal(extractors) result(self)
    !% Internal constructor for the ``multi'' output extractor property extractor class.
    implicit none
    type(nodePropertyExtractorMulti)                         :: self
    type(multiExtractorList        ), target , intent(in   ) :: extractors
    type(multiExtractorList        ), pointer                :: extractor_

    self      %extractors => extractors
    extractor_            => extractors
    do while (associated(extractor_))
       !# <referenceCountIncrement owner="extractor_" object="extractor_"/>
       extractor_ => extractor_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !% Destructor for the {\normalfont \ttfamily multi} output extractor property extractor class.
    implicit none
    type(nodePropertyExtractorMulti), intent(inout) :: self
    type(multiExtractorList        ), pointer       :: extractor_, extractorNext

    if (associated(self%extractors)) then
       extractor_ => self%extractors
       do while (associated(extractor_))
          extractorNext => extractor_%next
          !# <objectDestructor name="extractor_%extractor_"/>
          deallocate(extractor_)
          extractor_ => extractorNext
       end do
    end if
    return
  end subroutine multiDestructor

  integer function multiElementCount(self,elementType,time)
    !% Return the number of elements in the multiple property extractors.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (nodePropertyExtractorMulti), intent(inout) :: self
    integer                                     , intent(in   ) :: elementType
    double precision                            , intent(in   ) :: time
    type            (multiExtractorList        ), pointer       :: extractor_

    multiElementCount =  0
    extractor_        => self%extractors
    do while (associated(extractor_))
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) multiElementCount=multiElementCount+1
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) multiElementCount=multiElementCount+extractor_%elementCount(time)
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) multiElementCount=multiElementCount+1
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) multiElementCount=multiElementCount+extractor_%elementCount(time)
       class default
          call Galacticus_Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       extractor_ => extractor_%next
    end do
    return
  end function multiElementCount

  function multiExtractDouble(self,node,time,instance)
    !% Implement a multi output extractor.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                            , dimension(:) , allocatable :: multiExtractDouble
    class           (nodePropertyExtractorMulti), intent(inout)              :: self
    type            (treeNode                  ), intent(inout)              :: node
    double precision                            , intent(in   )              :: time
    type            (multiCounter              ), intent(inout), optional    :: instance
    type            (multiExtractorList        ), pointer                    :: extractor_
    integer                                                                  :: offset            , elementCount

    allocate(multiExtractDouble(self%elementCount(elementTypeDouble,time)))
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          elementCount=1
          multiExtractDouble(offset+1:offset+elementCount)=extractor_%extract(node     ,instance)
       class is (nodePropertyExtractorTuple        )
          elementCount=extractor_%elementCount(time)
          multiExtractDouble(offset+1:offset+elementCount)=extractor_%extract(node,time,instance)
       class is (nodePropertyExtractorIntegerScalar)
          elementCount=0
       class is (nodePropertyExtractorIntegerTuple )
          elementCount=0
       class default
          elementCount=0
          call Galacticus_Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiExtractDouble

  function multiExtractInteger(self,node,time,instance)
    !% Implement a multi output extractor.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    integer         (kind_int8                 ), dimension(:) , allocatable :: multiExtractInteger
    class           (nodePropertyExtractorMulti), intent(inout)              :: self
    type            (treeNode                  ), intent(inout)              :: node
    double precision                            , intent(in   )              :: time
    type            (multiCounter              ), intent(inout), optional    :: instance
    type            (multiExtractorList        ), pointer                    :: extractor_
    integer                                                                  :: offset             , elementCount

    allocate(multiExtractInteger(self%elementCount(elementTypeInteger,time)))
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          elementCount=0
       class is (nodePropertyExtractorTuple        )
          elementCount=0
       class is (nodePropertyExtractorIntegerScalar)
          elementCount=1
          multiExtractInteger(offset+1:offset+elementCount)=extractor_%extract(node,time,instance)
       class is (nodePropertyExtractorIntegerTuple )
          elementCount=extractor_%elementCount(time)
          multiExtractInteger(offset+1:offset+elementCount)=extractor_%extract(node,time,instance)
       class default
          elementCount=0
          call Galacticus_Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiExtractInteger

  subroutine multiAddInstances(self,node,instance)
    !% Implement adding of instances to a multi output extractor.
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

  function multiNames(self,elementType,time)
    !% Return the names of the multiple properties.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (varying_string            ), dimension(:) , allocatable :: multiNames
    class           (nodePropertyExtractorMulti), intent(inout)              :: self
    integer                                     , intent(in   )              :: elementType
    double precision                            , intent(in   )              :: time
    type            (multiExtractorList        ), pointer                    :: extractor_
    integer                                                                  :: offset    , elementCount

    allocate(multiNames(self%elementCount(elementType,time)))
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       elementCount=0
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) then
             elementCount=1
             multiNames(offset+1:offset+elementCount)=extractor_%name (    )
          end if
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             multiNames(offset+1:offset+elementCount)=extractor_%names(time)
          end if
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) then
             elementCount=1
             multiNames(offset+1:offset+elementCount)=extractor_%name (    )
          end if
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount(time)
             multiNames(offset+1:offset+elementCount)=extractor_%names(time)
          end if
       class default
          call Galacticus_Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiNames

  function multiDescriptions(self,elementType,time)
    !% Return the descriptions of the multiple properties.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (varying_string            ), dimension(:) , allocatable :: multiDescriptions
    class           (nodePropertyExtractorMulti), intent(inout)              :: self
    integer                                     , intent(in   )              :: elementType
    double precision                            , intent(in   )              :: time
    type            (multiExtractorList        ), pointer                    :: extractor_
    integer                                                                  :: offset           , elementCount

    allocate(multiDescriptions(self%elementCount(elementType,time)))
    offset     =  0
    extractor_ => self%extractors
    do while (associated(extractor_))
       elementCount=0
       select type (extractor_ => extractor_%extractor_)
       class is (nodePropertyExtractorScalar       )
          if (elementType == elementTypeDouble ) then
             elementCount=1
             multiDescriptions(offset+1:offset+elementCount)=extractor_%description (    )
          end if
       class is (nodePropertyExtractorTuple        )
          if (elementType == elementTypeDouble ) then
             elementCount=extractor_%elementCount(time)
             multiDescriptions(offset+1:offset+elementCount)=extractor_%descriptions(time)
          end if
       class is (nodePropertyExtractorIntegerScalar)
          if (elementType == elementTypeInteger) then
             elementCount=1
             multiDescriptions(offset+1:offset+elementCount)=extractor_%description (    )
          end if
       class is (nodePropertyExtractorIntegerTuple )
          if (elementType == elementTypeInteger) then
             elementCount=extractor_%elementCount(time)
             multiDescriptions(offset+1:offset+elementCount)=extractor_%descriptions(time)
          end if
       class default
          call Galacticus_Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiDescriptions

  function multiUnitsInSI(self,elementType,time)
    !% Return the units of the multiple properties in the SI system.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                            , dimension(:) , allocatable :: multiUnitsInSI
    class           (nodePropertyExtractorMulti), intent(inout)              :: self
    integer                                     , intent(in   )              :: elementType
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
       class default
          call Galacticus_Error_Report('unsupported property extractor type'//{introspection:location})
       end select
       offset     =  offset         +elementCount
       extractor_ => extractor_%next
    end do
    return
  end function multiUnitsInSI

  integer function multiType(self)
    !% Return the type of the multi property.
    use :: Output_Analyses_Options, only : outputAnalysisPropertyTypeLinear
    implicit none
    class(nodePropertyExtractorMulti), intent(inout) :: self
    !$GLC attributes unused :: self

    multiType=outputAnalysisPropertyTypeLinear
    return
  end function multiType
