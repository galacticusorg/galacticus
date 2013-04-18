!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a very simple model of mass movements during satellite mergers.

module Satellite_Merging_Mass_Movements_Very_Simple
  !% Implements a very simple model of mass movements during satellite mergers.
  use Satellite_Merging_Mass_Movements_Descriptors
  implicit none
  private
  public :: Satellite_Merging_Mass_Movements_Very_Simple_Initialize

contains

  !# <satelliteMergingMassMovementsMethod>
  !#  <unitName>Satellite_Merging_Mass_Movements_Very_Simple_Initialize</unitName>
  !# </satelliteMergingMassMovementsMethod>
  subroutine Satellite_Merging_Mass_Movements_Very_Simple_Initialize(satelliteMergingMassMovementsMethod,Satellite_Merging_Mass_Movement_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingMassMovementsMethod
    procedure(Satellite_Merging_Mass_Movement_Very_Simple),          pointer, intent(inout) :: Satellite_Merging_Mass_Movement_Get

    if (satelliteMergingMassMovementsMethod == 'verySimple') Satellite_Merging_Mass_Movement_Get => Satellite_Merging_Mass_Movement_Very_Simple
    return
  end subroutine Satellite_Merging_Mass_Movements_Very_Simple_Initialize

  subroutine Satellite_Merging_Mass_Movement_Very_Simple(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo,mergerIsMajor)
    !% Determine where stars and gas move as the result of a merger event using a simple algorithm.
    use Galacticus_Nodes
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer       , intent(  out)           :: gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo
    logical       , intent(  out)           :: mergerIsMajor

    mergerIsMajor  =.false.
    gasMovesTo     =movesToDisk
    starsMoveTo    =movesToDisk
    hostGasMovesTo =doesNotMove
    hostStarsMoveTo=doesNotMove

    return
  end subroutine Satellite_Merging_Mass_Movement_Very_Simple

end module Satellite_Merging_Mass_Movements_Very_Simple
