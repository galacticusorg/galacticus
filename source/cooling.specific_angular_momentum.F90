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

!% Contains a module that implements calculations of the specific angular momentum of cooling gas.

module Cooling_Specific_Angular_Momenta
  !% Implements calculations of the specific angular momentum of cooling gas.
  use ISO_Varying_String
  use Galacticus_Nodes
  !# <include directive="coolingSpecificAngularMomentumMethod" type="moduleUse">
  include 'cooling.specific_angular_momentum.modules.inc'
  !# </include>
  implicit none
  private
  public :: Cooling_Specific_Angular_Momentum

  ! Flag to indicate if this module has been initialized.  
  logical                                               :: coolingAngularMomentumInitialized    =.false.  
  
  ! Name of cooling radius available method used.                                                                                                     
  type     (varying_string                   )          :: coolingSpecificAngularMomentumMethod           
  
  ! Pointer to the function that actually does the calculation.                                                                                                     
  procedure(Cooling_Specific_Angular_Momentum), pointer :: Cooling_Specific_Angular_Momentum_Get=>null()  
                                                                                                       
contains

  subroutine Cooling_Specific_Angular_Momentum_Initialize
    !% Initialize the specific angular momentum of cooling gas module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.coolingAngularMomentumInitialized) then
       !$omp critical(Cooling_Specific_Angular_Momentum_Initialization) 
       if (.not.coolingAngularMomentumInitialized) then
          ! Get the cooling radius method parameter.
          !@ <inputParameter>
          !@   <name>coolingSpecificAngularMomentumMethod</name>
          !@   <defaultValue>constantRotation</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of the specific angular momentum of cooling gas.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('coolingSpecificAngularMomentumMethod',coolingSpecificAngularMomentumMethod,defaultValue='constantRotation')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="coolingSpecificAngularMomentumMethod" type="functionCall" functionType="void">
          !#  <functionArgs>coolingSpecificAngularMomentumMethod,Cooling_Specific_Angular_Momentum_Get</functionArgs>
          include 'cooling.specific_angular_momentum.inc'
          !# </include>
          if (.not.associated(Cooling_Specific_Angular_Momentum_Get)) call&
               & Galacticus_Error_Report('Cooling_Specific_Angular_Momentum','method ' //char(coolingSpecificAngularMomentumMethod)//' is unrecognized')
          coolingAngularMomentumInitialized=.true.
       end if
       !$omp end critical(Cooling_Specific_Angular_Momentum_Initialization) 
    end if
    return
  end subroutine Cooling_Specific_Angular_Momentum_Initialize

  double precision function Cooling_Specific_Angular_Momentum(thisNode,radius)
    !% Return the specific angular momentum (in units of km/s Mpc) of cooling gas in {\tt thisNode}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: radius    
    
    ! Initialize the module.                                                             
    call Cooling_Specific_Angular_Momentum_Initialize

    ! Get the cooling radius using the selected method.
    Cooling_Specific_Angular_Momentum=Cooling_Specific_Angular_Momentum_Get(thisNode,radius)

    return
  end function Cooling_Specific_Angular_Momentum

end module Cooling_Specific_Angular_Momenta
