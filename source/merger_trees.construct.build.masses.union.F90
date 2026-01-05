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
  Implementation of a merger tree masses class which constructs the union of other classes.
  !!}

  type, public :: mergerTreeBuildMassesList
     class(mergerTreeBuildMassesClass), pointer :: mergerTreeBuildMasses_ => null()
     type (mergerTreeBuildMassesList ), pointer :: next                   => null()
  end type mergerTreeBuildMassesList

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesUnion">
   <description>A merger tree masses class which constructs the union of other classes.</description>
   <linkedList type="mergerTreeBuildMassesList" variable="mergerTreeBuildMasses_" next="next" object="mergerTreeBuildMasses_" objectType="mergerTreeBuildMassesClass"/>
  </mergerTreeBuildMasses>
  !!]
  type, extends(mergerTreeBuildMassesClass) :: mergerTreeBuildMassesUnion
     !!{
     Implementation of a merger tree masses class which constructs the union of other classes.
     !!}
     private
     type(mergerTreeBuildMassesList), pointer :: mergerTreeBuildMasses_ => null()
   contains
     final     ::              unionDestructor
     procedure :: construct => unionConstruct
  end type mergerTreeBuildMassesUnion

  interface mergerTreeBuildMassesUnion
     module procedure unionConstructorParameters
     module procedure unionConstructorInternal
  end interface mergerTreeBuildMassesUnion

contains

  function unionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassesUnion} merger tree masses class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (mergerTreeBuildMassesUnion)                :: self
    type   (inputParameters           ), intent(inout) :: parameters
    type   (mergerTreeBuildMassesList ), pointer       :: mergerTreeBuildMasses_
    integer                                            :: i

    self%mergerTreeBuildMasses_ => null()
    mergerTreeBuildMasses_      => null()
    do i=1,parameters%copiesCount('mergerTreeBuildMasses',zeroIfNotPresent=.true.)
       if (associated(mergerTreeBuildMasses_)) then
          allocate(mergerTreeBuildMasses_%next)
          mergerTreeBuildMasses_ => mergerTreeBuildMasses_%next
       else
          allocate(self%mergerTreeBuildMasses_)
          mergerTreeBuildMasses_ => self                  %mergerTreeBuildMasses_
       end if
       mergerTreeBuildMasses_%mergerTreeBuildMasses_ => mergerTreeBuildMasses(parameters,i)
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="mergerTreeBuildMasses"/>
    !!]
    return
  end function unionConstructorParameters

  function unionConstructorInternal(mergerTreeBuildMasses_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassesUnion} merger tree masses class.
    !!}
    implicit none
    type(mergerTreeBuildMassesUnion)                        :: self
    type(mergerTreeBuildMassesList ), target, intent(in   ) :: mergerTreeBuildMasses_
    !![
    <constructorAssign variables="*mergerTreeBuildMasses_"/>
    !!]

    return
  end function unionConstructorInternal

  elemental subroutine unionDestructor(self)
    !!{
    Destructor for the merger tree mergerTreeBuildMasses function class.
    !!}
    implicit none
    type(mergerTreeBuildMassesUnion), intent(inout) :: self
    type(mergerTreeBuildMassesList ), pointer       :: mergerTreeBuildMasses_, mergerTreeBuildMassesNext

    if (associated(self%mergerTreeBuildMasses_)) then
       mergerTreeBuildMasses_ => self%mergerTreeBuildMasses_
       do while (associated(mergerTreeBuildMasses_))
          mergerTreeBuildMassesNext => mergerTreeBuildMasses_   %next
          deallocate(mergerTreeBuildMasses_%mergerTreeBuildMasses_)
          deallocate(mergerTreeBuildMasses_                       )
          mergerTreeBuildMasses_    => mergerTreeBuildMassesNext
       end do
    end if
    return
  end subroutine unionDestructor

  subroutine unionConstruct(self,time,mass,massMinimum,massMaximum,weight)
    !!{
    Construct a set of merger tree masses by sampling from a distribution.
    !!}
    use            :: Error            , only : Error_Report
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    use            :: Sorting          , only : sortIndex
    implicit none
    class           (mergerTreeBuildMassesUnion), intent(inout)                            :: self
    double precision                            , intent(in   )                            :: time
    double precision                            , intent(  out), allocatable, dimension(:) :: mass                  , weight        , &
         &                                                                                    massMinimum           , massMaximum
    type            (mergerTreeBuildMassesList  ), pointer                                 :: mergerTreeBuildMasses_
    double precision                                           , allocatable, dimension(:) :: massMember            , massTmp       , &
         &                                                                                    massMinimumMember     , massMinimumTmp, &
         &                                                                                    massMaximumMember     , massMaximumTmp, &
         &                                                                                    weightMember          , weightTmp
    integer         (c_size_t                   )              , allocatable, dimension(:) :: rankMass
    logical                                                                                :: useWeight
    integer         (c_size_t                   )                                          :: treeCount             , i

    useWeight              =  .false.
    treeCount              =  0_c_size_t
    mergerTreeBuildMasses_ => self%mergerTreeBuildMasses_
    do while (associated(mergerTreeBuildMasses_))
       call mergerTreeBuildMasses_%mergerTreeBuildMasses_%construct(time,massMember,massMinimumMember,massMaximumMember,weightMember)
       if (allocated(mass)) then
          if (useWeight.and..not.allocated(weightMember)) &
               & call Error_Report('all members must provide weights, or mass intervals - a mixture of the two is not permitted'//{introspection:location})
          call move_alloc(mass,     massTmp                  )
          allocate       (mass(size(massTmp)+size(massMember)))
          mass(1:size(massTmp))=massTmp
          deallocate(massTmp)
          if (useWeight) then
             call move_alloc(weight,     weightTmp                    )
             allocate       (weight(size(weightTmp)+size(weightMember)))
             weight(1:size(weightTmp))=weightTmp
             deallocate(weightTmp)
          else
             call move_alloc(massMinimum,     massMinimumTmp                         )
             call move_alloc(massMaximum,     massMaximumTmp                         )
             allocate       (massMinimum(size(massMinimumTmp)+size(massMinimumMember)))
             allocate       (massMaximum(size(massMaximumTmp)+size(massMaximumMember)))
             massMinimum(1:size(massMinimumTmp))=massMinimumTmp
             massMaximum(1:size(massMaximumTmp))=massMaximumTmp
             deallocate(massMinimumTmp)
             deallocate(massMaximumTmp)
          end if
       else
          useWeight=allocated(weightMember)
          allocate(mass,mold=massMaximumMember)
          if (useWeight) then
             allocate(weight     ,mold=     weightMember)
          else
             allocate(massMinimum,mold=massMinimumMember)
             allocate(massMaximum,mold=massMaximumMember)
          end if
       end if
       mass          (treeCount+1_c_size_t:treeCount+size(       massMember))=       massMember
       if (useWeight) then
          weight     (treeCount+1_c_size_t:treeCount+size(     weightMember))=     weightMember
       else
          massMinimum(treeCount+1_c_size_t:treeCount+size(massMinimumMember))=massMinimumMember
          massMaximum(treeCount+1_c_size_t:treeCount+size(massMaximumMember))=massMaximumMember
       end if
       treeCount=treeCount+size(massMember)
       mergerTreeBuildMasses_ => mergerTreeBuildMasses_%next
    end do
    ! Avoid overlapping mass intervals.
    if (.not.useWeight) then
       allocate(rankMass(size(mass)))
       rankMass=sortIndex(mass)
       do i=1,size(rankMass)
          if (i > 1 .and. massMinimum(rankMass(i)) < massMaximum(rankMass(i-1))) then
             massMinimum(rankMass(i  ))=sqrt(                     &
                  &                          +mass(rankMass(i-1)) &
                  &                          *mass(rankMass(i  )) &
                  &                         )
             massMaximum(rankMass(i-1))=sqrt(                     &
                  &                          +mass(rankMass(i-1)) &
                  &                          *mass(rankMass(i  )) &
                  &                         )
          end if
       end do

    end if
    return
  end subroutine unionConstruct
