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

!!{
Implements a node operator class that shifts preset named properties at node promotion.
!!}

  !![
  <nodeOperator name="nodeOperatorPresetNamedShift">
   <description>A node operator class that shifts preset named properties at node promotion.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorPresetNamedShift
     !!{
     A node operator class that shifts preset named properties at node promotion.
     !!}
     private
     logical                            :: presetsIdentified
     integer, allocatable, dimension(:) :: presetIndicesFloats, presetIndicesIntegers
   contains
     procedure :: nodePromote => presetNamedShiftNodePromote
  end type nodeOperatorPresetNamedShift

  interface nodeOperatorPresetNamedShift
     !!{
     Constructors for the \refClass{nodeOperatorPresetNamedShift} node operator class.
     !!}
     module procedure presetNamedShiftConstructorParameters
     module procedure presetNamedShiftConstructorInternal
  end interface nodeOperatorPresetNamedShift

contains

  function presetNamedShiftConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorPresetNamedShift} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorPresetNamedShift)                :: self
    type(inputParameters             ), intent(inout) :: parameters
    
    self=nodeOperatorPresetNamedShift()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function presetNamedShiftConstructorParameters

  function presetNamedShiftConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorPresetNamedShift} node operator class.
    !!}
    implicit none
    type(nodeOperatorPresetNamedShift) :: self

    self%presetsIdentified=.false.
    return
  end function presetNamedShiftConstructorInternal

  subroutine presetNamedShiftNodePromote(self,node)
    !!{
    Act on node promotion.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: ISO_Varying_String, only : extract           , operator(==)
    implicit none
    class  (nodeOperatorPresetNamedShift), intent(inout)             :: self
    type   (treeNode                    ), intent(inout)             :: node
    class  (nodeComponentBasic          ), pointer                   :: basic              , basicParent
    integer                              , dimension(:), allocatable :: indexMetaProperty
    type   (varying_string              )                            :: nameMetaProperty
    integer                                                          :: countMetaProperties, i          , &
         &                                                              ii

    ! Get the basic components of the two nodes.
    basic       => node       %basic()
    basicParent => node%parent%basic()
    ! Check if we need to first identify preset meta-properties for promotion.
    if (.not.self%presetsIdentified) then
       ! Identify preset meta-properties.
       self%presetsIdentified=.true.
       !! Float preset meta-properties.
       countMetaProperties=basic%countFloatRank0MetaProperties      ()
       allocate(indexMetaProperty(countMetaProperties))
       if (countMetaProperties > 0) then
          ii=0
          do i=1,countMetaProperties
             nameMetaProperty=basic%nameFloatRank0MetaProperty      (i)
             if (extract(nameMetaProperty,1,13) == "basic:preset:") then
                ii=ii+1
                indexMetaProperty(ii)=i
             end if
          end do
          allocate(self%presetIndicesFloats  (ii))
          self%presetIndicesFloats  =indexMetaProperty(1:ii)
          deallocate(indexMetaProperty)
       else
          allocate(self%presetIndicesFloats  (0 ))
       end if
       !! Integer preset meta-properties.
       countMetaProperties=basic%countLongIntegerRank0MetaProperties()
       allocate(indexMetaProperty(countMetaProperties))
       if (countMetaProperties > 0) then
          ii=0
          do i=1,countMetaProperties
             nameMetaProperty=basic%nameLongIntegerRank0MetaProperty(i)
             if (extract(nameMetaProperty,1,13) == "basic:preset:") then
                ii=ii+1
                indexMetaProperty(ii)=i
             end if
          end do
          allocate(self%presetIndicesIntegers(ii))
          self%presetIndicesIntegers=indexMetaProperty(1:ii)
          deallocate(indexMetaProperty)
       else
          allocate(self%presetIndicesIntegers(0 ))
       end if
    end if
    ! Promote float meta-properties.
    if (size(self%presetIndicesFloats) > 0) then
       do i=1,size(self%presetIndicesFloats)
          call basic%      floatRank0MetaPropertySet(self%presetIndicesFloats  (i),basicParent%      floatRank0MetaPropertyGet(self%presetIndicesFloats  (i)))
       end do
    end if
    ! Promote integer meta-properties.
    if (size(self%presetIndicesIntegers) > 0) then
       do i=1,size(self%presetIndicesIntegers)
          call basic%longIntegerRank0MetaPropertySet(self%presetIndicesIntegers(i),basicParent%longIntegerRank0MetaPropertyGet(self%presetIndicesIntegers(i)))
       end do
    end if
    return
  end subroutine presetNamedShiftNodePromote
