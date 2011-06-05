!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements a simple model of mass movements during satellite mergers.

module Satellite_Merging_Mass_Movements_Simple
  !% Implements a simple model of mass movements during satellite mergers.
  use Satellite_Merging_Mass_Movements_Descriptors
  private
  public :: Satellite_Merging_Mass_Movements_Simple_Initialize

  ! Mass ratio above which a merger is considered to be "major".
  double precision :: majorMergerMassRatio

  ! Location to which gas from satellite galaxy in minor merger is moved.
  integer          :: minorMergerGasMovesToValue

contains

  !# <satelliteMergingMassMovementsMethod>
  !#  <unitName>Satellite_Merging_Mass_Movements_Simple_Initialize</unitName>
  !# </satelliteMergingMassMovementsMethod>
  subroutine Satellite_Merging_Mass_Movements_Simple_Initialize(satelliteMergingMassMovementsMethod,Satellite_Merging_Mass_Movement_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingMassMovementsMethod
    procedure(),          pointer, intent(inout) :: Satellite_Merging_Mass_Movement_Get
    character(len=10)                            :: minorMergerGasMovesTo

    if (satelliteMergingMassMovementsMethod == 'simple') then
       Satellite_Merging_Mass_Movement_Get => Satellite_Merging_Mass_Movement_Simple
       !@ <inputParameter>
       !@   <name>majorMergerMassRatio</name>
       !@   <defaultValue>0.3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass ratio above which mergers are considered to be ``major'' in the simple merger mass movements method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("majorMergerMassRatio",majorMergerMassRatio,defaultValue=0.3d0)
       !@ <inputParameter>
       !@   <name>minorMergerGasMovesTo</name>
       !@   <defaultValue>disk</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The component to which satellite galaxy gas moves to as a result of a minor merger.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("minorMergerGasMovesTo",minorMergerGasMovesTo,defaultValue="disk")
       select case (trim(minorMergerGasMovesTo))
       case ("disk")
          minorMergerGasMovesToValue=movesToDisk
       case ("spheroid")
          minorMergerGasMovesToValue=movesToSpheroid
       case default
          call Galacticus_Error_Report('Satellite_Merging_Mass_Movements_Simple_Initialize','unrecognized location for minor merger satellite gas')
       end select
    end if
    return
  end subroutine Satellite_Merging_Mass_Movements_Simple_Initialize

  subroutine Satellite_Merging_Mass_Movement_Simple(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo)
    !% Return orbital velocities of a satellite selected at random from the fitting function found by \cite{benson_orbital_2005}.
    use Tree_Nodes
    use Tree_Node_Methods
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(out)             :: gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo
    type(treeNode),                pointer  :: hostNode
    double precision                        :: satelliteMass,hostMass

    ! Get the host node.
    hostNode => thisNode%parentNode

    ! Find the baryonic masses of the two galaxies.
    satelliteMass=Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic)
    hostMass     =Galactic_Structure_Enclosed_Mass(hostNode,massType=massTypeGalactic)

    ! Decide if the mass ratio is large enough to trigger a major merger.
    if (satelliteMass >= majorMergerMassRatio*hostMass) then
       gasMovesTo     =movesToSpheroid
       starsMoveTo    =movesToSpheroid
       hostGasMovesTo =movesToSpheroid
       hostStarsMoveTo=movesToSpheroid
    else
       gasMovesTo     =minorMergerGasMovesToValue
       starsMoveTo    =movesToSpheroid
       hostGasMovesTo =doesNotMove
       hostStarsMoveTo=doesNotMove
    end if

    return
  end subroutine Satellite_Merging_Mass_Movement_Simple

end module Satellite_Merging_Mass_Movements_Simple
