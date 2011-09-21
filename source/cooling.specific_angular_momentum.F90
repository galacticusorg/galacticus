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


!% Contains a module that implements calculations of the specific angular momentum of cooling gas.

module Cooling_Specific_Angular_Momenta
  !% Implements calculations of the specific angular momentum of cooling gas.
  use ISO_Varying_String
  use Tree_Nodes
  !# <include directive="coolingSpecificAngularMomentumMethod" type="moduleUse">
  include 'cooling.specific_angular_momentum.modules.inc'
  !# </include>
  implicit none
  private
  public :: Cooling_Specific_Angular_Momentum

  ! Flag to indicate if this module has been initialized.  
  logical              :: coolingAngularMomentumInitialized=.false.

  ! Name of cooling radius available method used.
  type(varying_string) :: coolingSpecificAngularMomentumMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Cooling_Specific_Angular_Momentum_Get_Template), pointer :: Cooling_Specific_Angular_Momentum_Get => null()
  abstract interface
     double precision function Cooling_Specific_Angular_Momentum_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Cooling_Specific_Angular_Momentum_Get_Template
  end interface
  
contains

  subroutine Cooling_Specific_Angular_Momentum_Initialize
    !% Initialize the specific angular momentum of cooling gas module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    !$omp critical(Cooling_Specific_Angular_Momentum_Initialization) 
    ! Initialize if necessary.
    if (.not.coolingAngularMomentumInitialized) then
       ! Get the cooling radius method parameter.
       !@ <inputParameter>
       !@   <name>coolingSpecificAngularMomentumMethod</name>
       !@   <defaultValue>constantRotation</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of the specific angular momentum of cooling gas.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingSpecificAngularMomentumMethod',coolingSpecificAngularMomentumMethod,defaultValue='constantRotation')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="coolingSpecificAngularMomentumMethod" type="code" action="subroutine">
       !#  <subroutineArgs>coolingSpecificAngularMomentumMethod,Cooling_Specific_Angular_Momentum_Get</subroutineArgs>
       include 'cooling.specific_angular_momentum.inc'
       !# </include>
       if (.not.associated(Cooling_Specific_Angular_Momentum_Get)) call&
            & Galacticus_Error_Report('Cooling_Specific_Angular_Momentum','method ' //char(coolingSpecificAngularMomentumMethod)//' is unrecognized')
       coolingAngularMomentumInitialized=.true.
    end if
    !$omp end critical(Cooling_Specific_Angular_Momentum_Initialization) 
    return
  end subroutine Cooling_Specific_Angular_Momentum_Initialize

  double precision function Cooling_Specific_Angular_Momentum(thisNode)
    !% Return the specific angular momentum (in units of km/s Mpc) of cooling gas in {\tt thisNode}.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Cooling_Specific_Angular_Momentum_Initialize

    ! Get the cooling radius using the selected method.
    Cooling_Specific_Angular_Momentum=Cooling_Specific_Angular_Momentum_Get(thisNode)

    return
  end function Cooling_Specific_Angular_Momentum

end module Cooling_Specific_Angular_Momenta
