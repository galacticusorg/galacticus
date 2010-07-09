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
