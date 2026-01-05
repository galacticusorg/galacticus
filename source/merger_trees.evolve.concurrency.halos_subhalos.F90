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
  Implements a merger tree evolution concurrency class that separates halo and subhalo evolution.
  !!}

  !![
  <mergerTreeEvolveConcurrency name="mergerTreeEvolveConcurrencyHalosSubhalos">
    <description>
    A merger tree evolution concurrency class which separates evolution of halos and their subhalos.
    </description>
  </mergerTreeEvolveConcurrency>
  !!]
  type, extends(mergerTreeEvolveConcurrencyClass) :: mergerTreeEvolveConcurrencyHalosSubhalos
     !!{
     Implementation of a merger tree evolution concurrency class which separates evolution of halos and their subhalos.
     !!}
     private
     integer :: depthHierarchyMaximum
   contains
     procedure :: initializeTree     => halosSubhalosInitializeTree
     procedure :: countPhases        => halosSubhalosCountPhases
     procedure :: includeInEvolution => halosSubhalosIncludeInEvolution
  end type mergerTreeEvolveConcurrencyHalosSubhalos

  interface mergerTreeEvolveConcurrencyHalosSubhalos
     !!{
     Constructors for the \refClass{mergerTreeEvolveConcurrencyHalosSubhalos} merger tree evolution concurrency model.
     !!}
     module procedure threadedConstructorParameters
     module procedure threadedConstructorInternal
  end interface mergerTreeEvolveConcurrencyHalosSubhalos

contains

  function threadedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveConcurrencyHalosSubhalos} merger tree evolution concurrency model class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(mergerTreeEvolveConcurrencyHalosSubhalos)                :: self
    type(inputParameters                         ), intent(inout) :: parameters

    self=mergerTreeEvolveConcurrencyHalosSubhalos()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function threadedConstructorParameters

  function threadedConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeEvolveConcurrencyHalosSubhalos} merger tree evolution concurrency model class.
    !!}
    implicit none
    type(mergerTreeEvolveConcurrencyHalosSubhalos) :: self

    self%depthHierarchyMaximum=-1
    return
  end function threadedConstructorInternal

  subroutine halosSubhalosInitializeTree(self)
    !!{
    Initialize state for evolution of a new tree.
    !!}
    implicit none
    class(mergerTreeEvolveConcurrencyHalosSubhalos), intent(inout) :: self

    self%depthHierarchyMaximum=-1
    return
  end subroutine halosSubhalosInitializeTree

  integer function halosSubhalosCountPhases(self) result(countPhases)
    !!{
    Return the number of evolution phases needed for this tree.
    !!}
    implicit none
    class(mergerTreeEvolveConcurrencyHalosSubhalos), intent(inout) :: self

    countPhases=max(1,self%depthHierarchyMaximum+1)
    return
  end function halosSubhalosCountPhases

  logical function halosSubhalosIncludeInEvolution(self,evolutionPhase,node) result(includeInEvolution)
    !!{
    Return true if the {\normalfont \ttfamily node} should be included in this phase of evolution.
    !!}
    implicit none
    class  (mergerTreeEvolveConcurrencyHalosSubhalos), intent(inout)          :: self
    integer                                          , intent(in   )          :: evolutionPhase
    type   (treeNode                                ), intent(inout), target  :: node
    integer                                                                   :: depthHierarchy, depthHierarchyTarget
    type   (treeNode                                )               , pointer :: nodeWork

    includeInEvolution   =  .false.
    depthHierarchy       =                 0
    depthHierarchyTarget =  evolutionPhase-1
    nodeWork             => node
    do while (nodeWork%isSatellite())
       ! Exit early if we have already exceeded the target hierarchy depth (except on the first phase of evolution in which we need to record the maximum hierarchy depth reached).
       if     (                                        &
            &   depthHierarchy >= depthHierarchyTarget &
            &  .and.                                   &
            &   evolutionPhase /= 1                    &
            & ) return
       depthHierarchy =  depthHierarchy       +1
       nodeWork       => nodeWork      %parent
    end do
    if (evolutionPhase == 1) self%depthHierarchyMaximum=max(self%depthHierarchyMaximum,depthHierarchy)
    includeInEvolution=(depthHierarchy == depthHierarchyTarget)
    return
  end function halosSubhalosIncludeInEvolution
