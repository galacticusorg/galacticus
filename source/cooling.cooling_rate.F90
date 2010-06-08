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






!% Contains a module that implements calculations of the cooling rate.

module Cooling_Rates
  !% Implements calculations of the cooling rate.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Cooling_Rate

  ! Flag to indicate if this module has been initialized.  
  logical              :: coolingRateInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: coolingRateMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Cooling_Rate_Get_Template), pointer :: Cooling_Rate_Get => null()
  abstract interface
     double precision function Cooling_Rate_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Cooling_Rate_Get_Template
  end interface
  
contains

  double precision function Cooling_Rate(thisNode)
    !% Return the cooling rate for {\tt thisNode} (in units of $M_\odot$ Gyr$^{-1}$).
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="coolingRateMethod" type="moduleUse">
    include 'cooling.cooling_rate.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    !$omp critical(Cooling_Rate_Initialization) 
    ! Initialize if necessary.
    if (.not.coolingRateInitialized) then
       ! Get the cooling rate method parameter.
       !@ <inputParameter>
       !@   <name>coolingRateMethod</name>
       !@   <defaultValue>White + Frenk</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used when computing the cooling rate.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingRateMethod',coolingRateMethod,defaultValue='White + Frenk')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="coolingRateMethod" type="code" action="subroutine">
       !#  <subroutineArgs>coolingRateMethod,Cooling_Rate_Get</subroutineArgs>
       include 'cooling.cooling_rate.inc'
       !# </include>
       if (.not.associated(Cooling_Rate_Get)) call Galacticus_Error_Report('Cooling_Rate','method ' &
            &//char(coolingRateMethod)//' is unrecognized')
       coolingRateInitialized=.true.
    end if
    !$omp end critical(Cooling_Rate_Initialization) 

    ! Get the cooling rate using the selected method.
    Cooling_Rate=Cooling_Rate_Get(thisNode)

    return
  end function Cooling_Rate

end module Cooling_Rates
