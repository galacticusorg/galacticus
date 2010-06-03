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






!% Contains a module that implements calculations of the cooling radius.

module Cooling_Radii
  !% Implements calculations of the cooling radius.
  use ISO_Varying_String
  use Tree_Nodes
  !# <include directive="coolingRadiusMethod" type="moduleUse">
  include 'cooling.cooling_radius.modules.inc'
  !# </include>
  private
  public :: Cooling_Radius, Cooling_Radius_Growth_Rate

  ! Flag to indicate if this module has been initialized.  
  logical              :: coolingRadiusInitialized=.false.

  ! Name of cooling radius available method used.
  type(varying_string) :: coolingRadiusMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Cooling_Radius_Get_Template), pointer :: Cooling_Radius_Get => null()
  procedure(Cooling_Radius_Get_Template), pointer :: Cooling_Radius_Growth_Rate_Get => null()
  abstract interface
     double precision function Cooling_Radius_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Cooling_Radius_Get_Template
  end interface
  
contains

  subroutine Cooling_Radius_Initialize
    !% Initialize the cooling radius module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    !$omp critical(Cooling_Radius_Initialization) 
    ! Initialize if necessary.
    if (.not.coolingRadiusInitialized) then
       ! Get the cooling radius method parameter.
       !@ <inputParameter>
       !@   <name>coolingRadiusMethod</name>
       !@   <defaultValue>simple</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of the cooling radius.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingRadiusMethod',coolingRadiusMethod,defaultValue='simple')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="coolingRadiusMethod" type="code" action="subroutine">
       !#  <subroutineArgs>coolingRadiusMethod,Cooling_Radius_Get,Cooling_Radius_Growth_Rate_Get</subroutineArgs>
       include 'cooling.cooling_radius.inc'
       !# </include>
       if (.not.(associated(Cooling_Radius_Get).and.associated(Cooling_Radius_Growth_Rate_Get))) call&
            & Galacticus_Error_Report('Cooling_Radius','method ' //char(coolingRadiusMethod)//' is unrecognized')
       coolingRadiusInitialized=.true.
    end if
    !$omp end critical(Cooling_Radius_Initialization) 
    return
  end subroutine Cooling_Radius_Initialize

  double precision function Cooling_Radius(thisNode)
    !% Return the cooling radius for {\tt thisNode} (in units of Mpc).
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Cooling_Radius_Initialize

    ! Get the cooling radius using the selected method.
    Cooling_Radius=Cooling_Radius_Get(thisNode)

    return
  end function Cooling_Radius

  double precision function Cooling_Radius_Growth_Rate(thisNode)
    !% Return the rate at which the cooling radius grows for {\tt thisNode} (in units of Mpc/Gyr).
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Cooling_Radius_Initialize

    ! Get the cooling radius using the selected method.
    Cooling_Radius_Growth_Rate=Cooling_Radius_Growth_Rate_Get(thisNode)

    return
  end function Cooling_Radius_Growth_Rate

end module Cooling_Radii
