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


!% Contains a module which implements calculations of expulsive feedback from star formation in spheroids.

module Star_Formation_Feedback_Expulsion_Spheroids
  !% Implements calculations of expulsive feedback from star formation in spheroids.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationExpulsiveFeedbackSpheroidsInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationExpulsiveFeedbackSpheroidsMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Expulsive_Feedback_Spheroid_Rate_Template), pointer :: Star_Formation_Expulsive_Feedback_Spheroid_Rate_Get => null()
  abstract interface
     double precision function Star_Formation_Expulsive_Feedback_Spheroid_Rate_Template(thisNode,starFormationRate,energyInputRate)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: starFormationRate,energyInputRate
     end function Star_Formation_Expulsive_Feedback_Spheroid_Rate_Template
  end interface

contains

  subroutine Star_Formation_Expulsive_Feedback_Spheroids_Initialize
    !% Initialize the spheroid star formation expulsive feedback module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationExpulsiveFeedbackSpheroidsMethod" type="moduleUse">
    include 'star_formation.feedbacks_expulsive.spheroids.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Star_Formation_Expulsive_Feedback_Spheroids_Initialization) 
    ! Initialize if necessary.
    if (.not.starFormationExpulsiveFeedbackSpheroidsInitialized) then
       ! Get the spheroid star formation expulsive feedback method parameter.
       !@ <inputParameter>
       !@   <name>starFormationExpulsiveFeedbackSpheroidsMethod</name>
       !@   <defaultValue>null</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of expulsive \SNe\ feedback in spheroids.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationExpulsiveFeedbackSpheroidsMethod',starFormationExpulsiveFeedbackSpheroidsMethod,defaultValue='null')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="starFormationExpulsiveFeedbackSpheroidsMethod" type="code" action="subroutine">
       !#  <subroutineArgs>starFormationExpulsiveFeedbackSpheroidsMethod,Star_Formation_Expulsive_Feedback_Spheroid_Rate_Get</subroutineArgs>
       include 'star_formation.feedback_expulsive.spheroids.inc'
       !# </include>
       if (.not.associated(Star_Formation_Expulsive_Feedback_Spheroid_Rate_Get)) call Galacticus_Error_Report('Star_Formation_Expulsive_Feedback_Spheroids'&
            &,'method ' //char(starFormationExpulsiveFeedbackSpheroidsMethod)//' is unrecognized')
       starFormationExpulsiveFeedbackSpheroidsInitialized=.true.
    end if
    !$omp end critical(Star_Formation_Expulsive_Feedback_Spheroids_Initialization) 

    return
  end subroutine Star_Formation_Expulsive_Feedback_Spheroids_Initialize

  double precision function Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
    !% Returns the expulsive outflow rate due to star formation in the spheroid component of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: starFormationRate,energyInputRate

    ! Initialize the module.
    call Star_Formation_Expulsive_Feedback_Spheroids_Initialize

    ! Get the energy using the selected method.
    Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate=Star_Formation_Expulsive_Feedback_Spheroid_Rate_Get(thisNode,starFormationRate,energyInputRate)

    return
  end function Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate
  
end module Star_Formation_Feedback_Expulsion_Spheroids
