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

!% Contains a module that implements calculations of the time available for cooling.

module Cooling_Times_Available
  use ISO_Varying_String
  use Tree_Nodes
  !# <include directive="coolingTimeAvailableMethod" type="moduleUse">
  include 'cooling.time_available.modules.inc'
  !# </include>
  implicit none
  private
  public :: Cooling_Time_Available, Cooling_Time_Available_Increase_Rate

  ! Flag to indicate if this module has been initialized.  
  logical              :: coolingTimeAvailableInitialized=.false.

  ! Name of cooling time available method used.
  type(varying_string) :: coolingTimeAvailableMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Cooling_Time_Available_Get_Template), pointer :: Cooling_Time_Available_Get => null()
  procedure(Cooling_Time_Available_Get_Template), pointer :: Cooling_Time_Available_Increase_Rate_Get => null()
  interface Cooling_Time_Available_Get_Template
     double precision function Cooling_Time_Available_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Cooling_Time_Available_Get_Template
  end interface
  
contains

  subroutine Cooling_Time_Available_Initialize
    !% Initializes the cooling time available module.
    use Galacticus_Error
    use Input_Parameters
    implicit none
    
    ! Initialize if necessary.
    if (.not.coolingTimeAvailableInitialized) then
       !$omp critical(Cooling_Time_Available_Initialization) 
       if (.not.coolingTimeAvailableInitialized) then
          ! Get the cooling time available method parameter.
          !@ <inputParameter>
          !@   <name>coolingTimeAvailableMethod</name>
          !@   <defaultValue>White-Frenk1991</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used when computing the time available for cooling.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('coolingTimeAvailableMethod',coolingTimeAvailableMethod,defaultValue='White-Frenk1991')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="coolingTimeAvailableMethod" type="code" action="subroutine">
          !#  <subroutineArgs>coolingTimeAvailableMethod,Cooling_Time_Available_Get,Cooling_Time_Available_Increase_Rate_Get</subroutineArgs>
          include 'cooling.time_available.inc'
          !# </include>
          if (.not.(associated(Cooling_Time_Available_Get).and.associated(Cooling_Time_Available_Increase_Rate_Get))) call&
               & Galacticus_Error_Report('Cooling_Time_Available','method ' //char(coolingTimeAvailableMethod)//' is unrecognized')
          coolingTimeAvailableInitialized=.true.
       end if
       !$omp end critical(Cooling_Time_Available_Initialization) 
    end if
    return
  end subroutine Cooling_Time_Available_Initialize

  double precision function Cooling_Time_Available(thisNode)
    !% Return the time available for cooling in {\tt thisNode}.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize if necessary.
    call Cooling_Time_Available_Initialize

    ! Get the cooling time using the selected method.
    Cooling_Time_Available=Cooling_Time_Available_Get(thisNode)

    return
  end function Cooling_Time_Available

  double precision function Cooling_Time_Available_Increase_Rate(thisNode)
    !% Return the rate at which the time available for cooling increases in {\tt thisNode}.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize if necessary.
    call Cooling_Time_Available_Initialize

    ! Get the cooling time using the selected method.
    Cooling_Time_Available_Increase_Rate=Cooling_Time_Available_Increase_Rate_Get(thisNode)

    return
  end function Cooling_Time_Available_Increase_Rate

end module Cooling_Times_Available
