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
  Implementation of a merger tree masses class which constructs a subset of another class.
  !!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesSubset">
   <description>A merger tree masses class which constructs a subset of another class.</description>
  </mergerTreeBuildMasses>
  !!]
  type, extends(mergerTreeBuildMassesClass) :: mergerTreeBuildMassesSubset
     !!{
     Implementation of a merger tree masses class which constructs a subset of another class.
     !!}
     private
     class  (mergerTreeBuildMassesClass), pointer :: mergerTreeBuildMasses_ => null()
     integer(c_size_t                  )          :: subsetBegin                     , subsetEnd
   contains
     final     ::              subsetDestructor
     procedure :: construct => subsetConstruct
  end type mergerTreeBuildMassesSubset

  interface mergerTreeBuildMassesSubset
     module procedure subsetConstructorParameters
     module procedure subsetConstructorInternal
  end interface mergerTreeBuildMassesSubset

contains

  function subsetConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassesSubset} merger tree masses class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeBuildMassesSubset)                :: self
    type   (inputParameters            ), intent(inout) :: parameters
    class  (mergerTreeBuildMassesClass ), pointer       :: mergerTreeBuildMasses_
    integer(c_size_t                   )                :: subsetBegin           , subsetEnd

    !![
    <inputParameter>
      <name>subsetBegin</name>
      <description>The first entry in the subset (indexed from 1).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>subsetEnd</name>
      <description>The last entry in the subset (indexed from 1).</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerTreeBuildMasses" name="mergerTreeBuildMasses_" source="parameters"/>
    !!]
    self=mergerTreeBuildMassesSubset(subsetBegin,subsetEnd,mergerTreeBuildMasses_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function subsetConstructorParameters

  function subsetConstructorInternal(subsetBegin,subsetEnd,mergerTreeBuildMasses_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassesSubset} merger tree masses class.
    !!}
    implicit none
    type   (mergerTreeBuildMassesSubset)                        :: self
    class  (mergerTreeBuildMassesClass ), target, intent(in   ) :: mergerTreeBuildMasses_
    integer(c_size_t                   )        , intent(in   ) :: subsetBegin           , subsetEnd
    !![
    <constructorAssign variables="subsetBegin, subsetEnd, *mergerTreeBuildMasses_"/>
    !!]

    return
  end function subsetConstructorInternal

  subroutine subsetDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeBuildMassesSubset} merger tree masses class.
    !!}
    implicit none
    type(mergerTreeBuildMassesSubset), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerTreeBuildMasses_"/>
    !!]
    return
  end subroutine subsetDestructor

  subroutine subsetConstruct(self,time,mass,massMinimum,massMaximum,weight)
    !!{
    Construct a set of merger tree masses by taking a subset from another set.
    !!}
    implicit none
    class           (mergerTreeBuildMassesSubset), intent(inout)                            :: self
    double precision                             , intent(in   )                            :: time
    double precision                             , intent(  out), allocatable, dimension(:) :: mass          , weight        , &
         &                                                                                     massMinimum   , massMaximum
    double precision                                            , allocatable, dimension(:) :: massSet       , weightSet     , &
         &                                                                                     massMinimumSet, massMaximumSet
    integer         (c_size_t                    )                                          :: subsetSize

    subsetSize=self%subsetEnd-self%subsetBegin+1
    call self%mergerTreeBuildMasses_%construct(time,massSet,massMinimumSet,massMaximumSet,weightSet)
    if (allocated(massSet       )) then
       allocate(mass       (subsetSize))
       mass       =massSet       (self%subsetBegin:self%subsetEnd)
    end if
    if (allocated(massMinimumSet)) then
       allocate(massMinimum(subsetSize))
       massMinimum=massMinimumSet(self%subsetBegin:self%subsetEnd)
    end if
    if (allocated(massMaximumSet)) then
       allocate(massMaximum(subsetSize))
       massMaximum=massMaximumSet(self%subsetBegin:self%subsetEnd)
    end if
    if (allocated(weightSet     )) then
       allocate(weight     (subsetSize))
       weight     =weightSet      (self%subsetBegin:self%subsetEnd)
    end if
    return
  end subroutine subsetConstruct
