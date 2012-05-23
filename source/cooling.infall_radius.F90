!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module that implements calculations of the infall radius for cooling calculations.

module Cooling_Infall_Radii
  !% Implements calculations of the infall radius for cooling calculations.
  use ISO_Varying_String
  use Tree_Nodes
  !# <include directive="infallRadiusMethod" type="moduleUse">
  include 'cooling.infall_radius.modules.inc'
  !# </include>
  implicit none
  private
  public :: Infall_Radius, Infall_Radius_Growth_Rate

  ! Flag to indicate if this module has been initialized.  
  logical              :: infallRadiusInitialized=.false.

  ! Name of cooling radius available method used.
  type(varying_string) :: infallRadiusMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Infall_Radius_Get_Template), pointer :: Infall_Radius_Get             => null()
  procedure(Infall_Radius_Get_Template), pointer :: Infall_Radius_Growth_Rate_Get => null()
  abstract interface
     double precision function Infall_Radius_Get_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Infall_Radius_Get_Template
  end interface
  
contains

  subroutine Infall_Radius_Initialize
    !% Initialize the infall radius module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.infallRadiusInitialized) then
       !$omp critical(Infall_Radius_Initialization) 
       if (.not.infallRadiusInitialized) then
          ! Get the infall radius method parameter.
          !@ <inputParameter>
          !@   <name>infallRadiusMethod</name>
          !@   <defaultValue>coolingRadius</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of the infall radius for cooling calculations.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('infallRadiusMethod',infallRadiusMethod,defaultValue='coolingRadius')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="infallRadiusMethod" type="code" action="subroutine">
          !#  <subroutineArgs>infallRadiusMethod,Infall_Radius_Get,Infall_Radius_Growth_Rate_Get</subroutineArgs>
          include 'cooling.infall_radius.inc'
          !# </include>
          if (.not.(associated(Infall_Radius_Get).and.associated(Infall_Radius_Growth_Rate_Get))) call&
               & Galacticus_Error_Report('Infall_Radius','method ' //char(infallRadiusMethod)//' is unrecognized')
          infallRadiusInitialized=.true.
       end if
       !$omp end critical(Infall_Radius_Initialization) 
    end if
    return
  end subroutine Infall_Radius_Initialize

  double precision function Infall_Radius(thisNode)
    !% Return the infall radius for {\tt thisNode} (in units of Mpc).
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Infall_Radius_Initialize

    ! Get the cooling radius using the selected method.
    Infall_Radius=Infall_Radius_Get(thisNode)

    return
  end function Infall_Radius

  double precision function Infall_Radius_Growth_Rate(thisNode)
    !% Return the rate at which the infall radius grows for {\tt thisNode} (in units of Mpc/Gyr).
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Infall_Radius_Initialize

    ! Get the cooling radius using the selected method.
    Infall_Radius_Growth_Rate=Infall_Radius_Growth_Rate_Get(thisNode)

    return
  end function Infall_Radius_Growth_Rate

end module Cooling_Infall_Radii
