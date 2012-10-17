!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which determines how mass is moved around as a consequence of a satellite merging event.

module Satellite_Merging_Mass_Movements
  !% Determines how mass is moved around as a consequence of a satellite merging event.
  use ISO_Varying_String
  implicit none
  private
  public :: Satellite_Merging_Mass_Movement_Store, Satellite_Merging_Mass_Movement

  ! Flag to indicate if this module has been initialized.  
  logical                                             :: satelliteMergingMassMovementsInitialized=.false.

  ! Name of mass movement method used.
  type(varying_string)                                :: satelliteMergingMassMovementsMethod

  ! Pointer to the subroutine that returns descriptors for mass movement.
  procedure(Satellite_Merging_Mass_Movement), pointer :: Satellite_Merging_Mass_Movement_Get => null()
 
contains

  !# <satelliteMergerTask>
  !#  <unitName>Satellite_Merging_Mass_Movement_Store</unitName>
  !# </satelliteMergerTask>
  subroutine Satellite_Merging_Mass_Movement_Store(thisNode)
    !% Compute and store the mass movement descriptors for this satellite merger.
    use Tree_Nodes
    use Satellite_Merging_Mass_Movements_Descriptors
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode

    call Satellite_Merging_Mass_Movement(thisNode,thisMergerGasMovesTo,thisMergerStarsMoveTo,thisHostGasMovesTo,thisHostStarsMoveTo,thisMergerIsMajor)
    return
  end subroutine Satellite_Merging_Mass_Movement_Store

  subroutine Satellite_Merging_Mass_Movement(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo,mergerIsMajor)
    !% Returns descriptors of how gas and stars move as the result of a satellite merger.
    use Tree_Nodes
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteMergingMassMovementsMethod" type="moduleUse">
    include 'satellites.merging.mass_movements.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer  :: thisNode
    integer,        intent(out)             :: gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo
    logical,        intent(out)             :: mergerIsMajor

    if (.not.satelliteMergingMassMovementsInitialized) then
       !$omp critical(satelliteMergingMassMovementsInitialize)
       if (.not.satelliteMergingMassMovementsInitialized) then
          ! Get the satellite merging mass movement method parameter.
          !@ <inputParameter>
          !@   <name>satelliteMergingMassMovementsMethod</name>
          !@   <defaultValue>simple</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@    Selects the method to be used for deciding mass movements during satellite mergers.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteMergingMassMovementsMethod',satelliteMergingMassMovementsMethod,defaultValue='simple')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="satelliteMergingMassMovementsMethod" type="code" action="subroutine">
          !#  <subroutineArgs>satelliteMergingMassMovementsMethod,Satellite_Merging_Mass_Movement_Get</subroutineArgs>
          include 'satellites.merging.mass_movements.inc'
          !# </include>
          if (.not.associated(Satellite_Merging_Mass_Movement_Get)) call Galacticus_Error_Report('Satellite_Merging_Mass_Movement','method ' &
               &//char(satelliteMergingMassMovementsMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          satelliteMergingMassMovementsInitialized=.true.
       end if
       !$omp end critical(satelliteMergingMassMovementsInitialize)
    end if

    ! Call the routine to get the descriptors.
    call Satellite_Merging_Mass_Movement_Get(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo,mergerIsMajor)

    return
  end subroutine Satellite_Merging_Mass_Movement
  
end module Satellite_Merging_Mass_Movements
