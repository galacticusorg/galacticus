!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which determines how mass is moved around as a consequence of a satellite merging event.

module Satellite_Merging_Mass_Movements
  !% Determines how mass is moved around as a consequence of a satellite merging event.
  use ISO_Varying_String
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

    call Satellite_Merging_Mass_Movement(thisNode,thisMergerGasMovesTo,thisMergerStarsMoveTo,thisHostGasMovesTo,thisHostStarsMoveTo)
    return
  end subroutine Satellite_Merging_Mass_Movement_Store

  subroutine Satellite_Merging_Mass_Movement(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo)
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

    ! Call the routine to get the descriptors.
    call Satellite_Merging_Mass_Movement_Get(thisNode,gasMovesTo,starsMoveTo,hostGasMovesTo,hostStarsMoveTo)

    return
  end subroutine Satellite_Merging_Mass_Movement
  
end module Satellite_Merging_Mass_Movements
