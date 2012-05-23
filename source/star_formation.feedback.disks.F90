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


!% Contains a module which implements calculations of feedback from star formation in disks.

module Star_Formation_Feedback_Disks
  !% Implements calculations of feedback from star formation in disks.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Star_Formation_Feedback_Disk_Outflow_Rate
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationFeedbackDisksInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationFeedbackDisksMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Feedback_Disk_Outflow_Rate_Template), pointer :: Star_Formation_Feedback_Disk_Outflow_Rate_Get => null()
  interface Star_Formation_Feedback_Disk_Outflow_Rate_Template
     double precision function Star_Formation_Feedback_Disk_Outflow_Rate_Template(thisNode,starFormationRate,energyInputRate)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: starFormationRate,energyInputRate
     end function Star_Formation_Feedback_Disk_Outflow_Rate_Template
  end interface

contains

  subroutine Star_Formation_Feedback_Disks_Initialize
    !% Initialize the disk star formation feedback module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationFeedbackDisksMethod" type="moduleUse">
    include 'star_formation.feedbacks.disks.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.starFormationFeedbackDisksInitialized) then
       !$omp critical(Star_Formation_Feedback_Disks_Initialization) 
       if (.not.starFormationFeedbackDisksInitialized) then
          ! Get the disk star formation feedback method parameter.
          !@ <inputParameter>
          !@   <name>starFormationFeedbackDisksMethod</name>
          !@   <defaultValue>powerLaw</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of \gls{sne} feedback in disks.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('starFormationFeedbackDisksMethod',starFormationFeedbackDisksMethod,defaultValue='powerLaw')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="starFormationFeedbackDisksMethod" type="code" action="subroutine">
          !#  <subroutineArgs>starFormationFeedbackDisksMethod,Star_Formation_Feedback_Disk_Outflow_Rate_Get</subroutineArgs>
          include 'star_formation.feedbacks.disks.inc'
          !# </include>
          if (.not.associated(Star_Formation_Feedback_Disk_Outflow_Rate_Get)) call Galacticus_Error_Report('Star_Formation_Feedback_Disks'&
               &,'method ' //char(starFormationFeedbackDisksMethod)//' is unrecognized')
          starFormationFeedbackDisksInitialized=.true.
       end if
       !$omp end critical(Star_Formation_Feedback_Disks_Initialization) 
    end if
    return
  end subroutine Star_Formation_Feedback_Disks_Initialize

  double precision function Star_Formation_Feedback_Disk_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate due to star formation in the disk component of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: starFormationRate,energyInputRate

    ! Initialize the module.
    call Star_Formation_Feedback_Disks_Initialize

    ! Get the energy using the selected method.
    Star_Formation_Feedback_Disk_Outflow_Rate=Star_Formation_Feedback_Disk_Outflow_Rate_Get(thisNode,starFormationRate,energyInputRate)

    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate
  
end module Star_Formation_Feedback_Disks
