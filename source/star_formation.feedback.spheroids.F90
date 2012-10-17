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

!% Contains a module which implements calculations of feedback from star formation in spheroids.

module Star_Formation_Feedback_Spheroids
  !% Implements calculations of feedback from star formation in spheroids.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Star_Formation_Feedback_Spheroid_Outflow_Rate
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationFeedbackSpheroidsInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationFeedbackSpheroidsMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Feedback_Spheroid_Outflow_Rate_Template), pointer :: Star_Formation_Feedback_Spheroid_Outflow_Rate_Get => null()
  interface Star_Formation_Feedback_Spheroid_Outflow_Rate_Template
     double precision function Star_Formation_Feedback_Spheroid_Outflow_Rate_Template(thisNode,starFormationRate,energyInputRate)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: starFormationRate,energyInputRate
     end function Star_Formation_Feedback_Spheroid_Outflow_Rate_Template
  end interface

contains

  subroutine Star_Formation_Feedback_Spheroids_Initialize
    !% Initialize the spheroid star formation feedback module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationFeedbackSpheroidsMethod" type="moduleUse">
    include 'star_formation.feedbacks.spheroids.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.starFormationFeedbackSpheroidsInitialized) then
       !$omp critical(Star_Formation_Feedback_Spheroids_Initialization) 
       if (.not.starFormationFeedbackSpheroidsInitialized) then
          ! Get the spheroid star formation feedback method parameter.
          !@ <inputParameter>
          !@   <name>starFormationFeedbackSpheroidsMethod</name>
          !@   <defaultValue>powerLaw</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of \gls{sne} feedback in spheroids.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('starFormationFeedbackSpheroidsMethod',starFormationFeedbackSpheroidsMethod,defaultValue='powerLaw')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="starFormationFeedbackSpheroidsMethod" type="code" action="subroutine">
          !#  <subroutineArgs>starFormationFeedbackSpheroidsMethod,Star_Formation_Feedback_Spheroid_Outflow_Rate_Get</subroutineArgs>
          include 'star_formation.feedbacks.spheroids.inc'
          !# </include>
          if (.not.associated(Star_Formation_Feedback_Spheroid_Outflow_Rate_Get)) call Galacticus_Error_Report('Star_Formation_Feedback_Spheroids'&
               &,'method ' //char(starFormationFeedbackSpheroidsMethod)//' is unrecognized')
          starFormationFeedbackSpheroidsInitialized=.true.
       end if
       !$omp end critical(Star_Formation_Feedback_Spheroids_Initialization) 
    end if
    return
  end subroutine Star_Formation_Feedback_Spheroids_Initialize

  double precision function Star_Formation_Feedback_Spheroid_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate due to star formation in the spheroid component of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: starFormationRate,energyInputRate

    ! Initialize the module.
    call Star_Formation_Feedback_Spheroids_Initialize

    ! Get the energy using the selected method.
    Star_Formation_Feedback_Spheroid_Outflow_Rate=Star_Formation_Feedback_Spheroid_Outflow_Rate_Get(thisNode,starFormationRate,energyInputRate)

    return
  end function Star_Formation_Feedback_Spheroid_Outflow_Rate
  
end module Star_Formation_Feedback_Spheroids
