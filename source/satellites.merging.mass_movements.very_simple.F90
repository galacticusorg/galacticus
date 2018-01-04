!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a very simple model of mass movements during satellite mergers.

module Satellite_Merging_Mass_Movements_Very_Simple
  !% Implements a very simple model of mass movements during satellite mergers.
  use Satellite_Merging_Mass_Movements_Descriptors
  implicit none
  private
  public :: Satellite_Merging_Mass_Movements_Very_Simple_Initialize

  ! Mass ratio above which a merger is considered to be "major".
  double precision :: majorMergerMassRatio

contains

  !# <satelliteMergingMassMovementsMethod>
  !#  <unitName>Satellite_Merging_Mass_Movements_Very_Simple_Initialize</unitName>
  !# </satelliteMergingMassMovementsMethod>
  subroutine Satellite_Merging_Mass_Movements_Very_Simple_Initialize(satelliteMergingMassMovementsMethod,Satellite_Merging_Mass_Movement_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                             ), intent(in   )          :: satelliteMergingMassMovementsMethod
    procedure(Satellite_Merging_Mass_Movement_Very_Simple), intent(inout), pointer :: Satellite_Merging_Mass_Movement_Get

    if (satelliteMergingMassMovementsMethod == 'verySimple') then
       Satellite_Merging_Mass_Movement_Get => Satellite_Merging_Mass_Movement_Very_Simple
       !# <inputParameter>
       !#   <name>majorMergerMassRatio</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>0.0d0</defaultValue>
       !#   <description>The mass ratio above which mergers are considered to be ``major'' in the very simple merger mass movements method.</description>
       !#   <source>globalParameters</source>
       !#   <type>real</type>
       !# </inputParameter>
    end if
    return
  end subroutine Satellite_Merging_Mass_Movements_Very_Simple_Initialize

  subroutine Satellite_Merging_Mass_Movement_Very_Simple(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo,mergerIsMajor)
    !% Determine where stars and gas move as the result of a merger event using a simple algorithm.
    use Galacticus_Nodes
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    type            (treeNode), intent(inout) :: thisNode
    integer                   , intent(  out) :: gasMovesTo   , hostGasMovesTo, hostStarsMoveTo, starsMoveTo
    logical                   , intent(  out) :: mergerIsMajor
    type            (treeNode), pointer       :: hostNode
    double precision                          :: hostMass     , satelliteMass
    
    ! Determine if this merger is considered major.
    if      (majorMergerMassRatio <= 0.0d0) then
       mergerIsMajor=.true.
    else if (majorMergerMassRatio >  1.0d0) then
       mergerIsMajor=.false.
    else
       ! Find the node to merge with.
       hostNode => thisNode%mergesWith()
       ! Find the baryonic masses of the two galaxies.
       satelliteMass=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic)
       hostMass     =Galactic_Structure_Enclosed_Mass(hostNode,massType=massTypeGalactic)    
       mergerIsMajor=satelliteMass >= majorMergerMassRatio*hostMass
    end if
    gasMovesTo     =movesToDisk
    starsMoveTo    =movesToDisk
    hostGasMovesTo =doesNotMove
    hostStarsMoveTo=doesNotMove
    return
  end subroutine Satellite_Merging_Mass_Movement_Very_Simple

end module Satellite_Merging_Mass_Movements_Very_Simple
