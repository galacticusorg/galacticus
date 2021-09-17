!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a null star formation histories class.
!!}

  !![
  <starFormationHistory name="starFormationHistoryNull">
   <description>A null star formation histories class.</description>
  </starFormationHistory>
  !!]
  type, extends(starFormationHistoryClass) :: starFormationHistoryNull
     !!{
     A null star formation histories class.
     !!}
     private
   contains
     procedure :: create                => nullCreate
     procedure :: rate                  => nullRate
     procedure :: output                => nullOutput
     procedure :: scales                => nullScales
     procedure :: metallicityBoundaries => nullMetallicityBoundaries
  end type starFormationHistoryNull

  interface starFormationHistoryNull
     !!{
     Constructors for the ``null'' star formation history class.
     !!}
     module procedure nullConstructorParameters
  end interface starFormationHistoryNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``null'' star formation history class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(starFormationHistoryNull)                :: self
    type(inputParameters         ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=starFormationHistoryNull()
    return
  end function nullConstructorParameters

  subroutine nullCreate(self,node,historyStarFormation,timeBegin)
    !!{
    Create the history required for storing star formation history.
    !!}
    implicit none
    class           (starFormationHistoryNull), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    type            (history                 ), intent(inout) :: historyStarFormation
    double precision                          , intent(in   ) :: timeBegin
    !$GLC attributes unused :: self, node, historyStarFormation, timeBegin

    ! Do nothing.
    return
  end subroutine nullCreate

  subroutine nullRate(self,node,historyStarFormation,abundancesFuel,rateStarFormation)
    !!{
    Set the rate the star formation history for {\normalfont \ttfamily node}.
    !!}
    implicit none
    class           (starFormationHistoryNull), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    type            (history                 ), intent(inout) :: historyStarFormation
    type            (abundances              ), intent(in   ) :: abundancesFuel
    double precision                          , intent(in   ) :: rateStarFormation
    !$GLC attributes unused :: self, node, historyStarFormation, abundancesFuel, rateStarFormation

    ! Ensure the history does not exist.
    call historyStarFormation%destroy()
    return
  end subroutine nullRate

  subroutine nullOutput(self,node,nodePassesFilter,historyStarFormation,indexOutput,indexTree,componentType,treeLock)
    !!{
    Output the star formation history for {\normalfont \ttfamily node}.
    !!}
    implicit none
    class  (starFormationHistoryNull), intent(inout)         :: self
    type   (treeNode                ), intent(inout), target :: node
    logical                          , intent(in   )         :: nodePassesFilter
    type   (history                 ), intent(inout)         :: historyStarFormation
    integer(c_size_t                ), intent(in   )         :: indexOutput
    integer(kind=kind_int8          ), intent(in   )         :: indexTree
    integer                            intent(in   )         :: componentType
    type   (ompLock                 ), intent(inout)         :: treeLock
    !$GLC attributes unused :: self, node, nodePassesFilter, historyStarFormation, indexOutput, indexTree, componentType, treeLock

    ! Do nothing.
    return
  end subroutine nullOutput

  subroutine nullScales(self,historyStarFormation,massStellar,abundancesStellar)
    !!{
    Set the scalings for error control on the absolute values of star formation histories.
    !!}
    implicit none
    class           (starFormationHistoryNull), intent(inout) :: self
    double precision                          , intent(in   ) :: massStellar
    type            (abundances              ), intent(in   ) :: abundancesStellar
    type            (history                 ), intent(inout) :: historyStarFormation
    !$GLC attributes unused :: self, historyStarFormation, massStellar, abundancesStellar

    ! Do nothing.
    return
  end subroutine nullScales

  function nullMetallicityBoundaries(self)
    !!{
    Return the boundaries of the metallicities used in this tabulation.
    !!}
    implicit none
    double precision                          , allocatable  , dimension(:) :: nullMetallicityBoundaries
    class           (starFormationHistoryNull), intent(inout)               :: self

    allocate(nullMetallicityBoundaries(0:1))
    nullMetallicityBoundaries(0:1)=[0.0d0,huge(0.0d0)]
    return
  end function nullMetallicityBoundaries
