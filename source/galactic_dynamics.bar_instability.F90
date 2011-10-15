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


!% Contains a module which implements calculations of bar instability in galactic disks.

module Galactic_Dynamics_Bar_Instabilities
  !% Implements calculations of bar instability in galactic disks.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Bar_Instability_Timescale

  ! Flag to indicate if this module has been initialized.  
  logical              :: barInstabilitiesInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: barInstabilityMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Bar_Instability_Template), pointer :: Bar_Instability_Timescale_Get => null()
  abstract interface
     double precision function Bar_Instability_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Bar_Instability_Template
  end interface

contains

  subroutine Galactic_Dynamics_Bar_Instability_Initialize
    !% Initialize the bar instability module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="barInstabilityMethod" type="moduleUse">
    include 'galactic_dynamics.bar_instability.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Galactic_Dynamics_Bar_Instability_Initialize) 
    ! Initialize if necessary.
    if (.not.barInstabilitiesInitialized) then
       ! Get the halo spin distribution method parameter.
       !@ <inputParameter>
       !@   <name>barInstabilityMethod</name>
       !@   <defaultValue>ELN</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for bar instability calculations.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('barInstabilityMethod',barInstabilityMethod,defaultValue='ELN')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="barInstabilityMethod" type="code" action="subroutine">
       !#  <subroutineArgs>barInstabilityMethod,Bar_Instability_Timescale_Get</subroutineArgs>
       include 'galactic_dynamics.bar_instability.inc'
       !# </include>
       if (.not.associated(Bar_Instability_Timescale_Get)) &
            & call Galacticus_Error_Report('Galactic_Dynamics_Bar_Instability_Initialize','method ' //char(barInstabilityMethod)//' is unrecognized')
       barInstabilitiesInitialized=.true.
    end if
    !$omp end critical(Galactic_Dynamics_Bar_Instability_Initialize) 

    return
  end subroutine Galactic_Dynamics_Bar_Instability_Initialize

  double precision function Bar_Instability_Timescale(thisNode)
    !% Returns a timescale on which the bar instability depletes material from a disk into a pseudo-bulge. A negative value
    !% indicates no instability.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Galactic_Dynamics_Bar_Instability_Initialize

    ! Get the timescale using the selected method.
    Bar_Instability_Timescale=Bar_Instability_Timescale_Get(thisNode)

    return
  end function Bar_Instability_Timescale

end module Galactic_Dynamics_Bar_Instabilities
