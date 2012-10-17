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

!% Contains a module which implements calculations of expulsive feedback from star formation in disks.

module Star_Formation_Feedback_Expulsion_Disks
  !% Implements calculations of expulsive feedback from star formation in disks.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationExpulsiveFeedbackDisksInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationExpulsiveFeedbackDisksMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Expulsive_Feedback_Disk_Rate_Template), pointer :: Star_Formation_Expulsive_Feedback_Disk_Rate_Get => null()
  abstract interface
     double precision function Star_Formation_Expulsive_Feedback_Disk_Rate_Template(thisNode,starFormationRate,energyInputRate)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: starFormationRate,energyInputRate
     end function Star_Formation_Expulsive_Feedback_Disk_Rate_Template
  end interface

contains

  subroutine Star_Formation_Expulsive_Feedback_Disks_Initialize
    !% Initialize the disk star formation expulsive feedback module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationExpulsiveFeedbackDisksMethod" type="moduleUse">
    include 'star_formation.feedbacks_expulsive.disks.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.starFormationExpulsiveFeedbackDisksInitialized) then
       !$omp critical(Star_Formation_Expulsive_Feedback_Disks_Initialization) 
       if (.not.starFormationExpulsiveFeedbackDisksInitialized) then
          ! Get the disk star formation expulsive feedback method parameter.
          !@ <inputParameter>
          !@   <name>starFormationExpulsiveFeedbackDisksMethod</name>
          !@   <defaultValue>null</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for calculations of expulsive \gls{sne} feedback in disks.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('starFormationExpulsiveFeedbackDisksMethod',starFormationExpulsiveFeedbackDisksMethod,defaultValue='null')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="starFormationExpulsiveFeedbackDisksMethod" type="code" action="subroutine">
          !#  <subroutineArgs>starFormationExpulsiveFeedbackDisksMethod,Star_Formation_Expulsive_Feedback_Disk_Rate_Get</subroutineArgs>
          include 'star_formation.feedback_expulsive.disks.inc'
          !# </include>
          if (.not.associated(Star_Formation_Expulsive_Feedback_Disk_Rate_Get)) call Galacticus_Error_Report('Star_Formation_Expulsive_Feedback_Disks'&
               &,'method ' //char(starFormationExpulsiveFeedbackDisksMethod)//' is unrecognized')
          starFormationExpulsiveFeedbackDisksInitialized=.true.
       end if
       !$omp end critical(Star_Formation_Expulsive_Feedback_Disks_Initialization) 
    end if

    return
  end subroutine Star_Formation_Expulsive_Feedback_Disks_Initialize

  double precision function Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate(thisNode,starFormationRate,energyInputRate)
    !% Returns the expulsive outflow rate due to star formation in the disk component of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: starFormationRate,energyInputRate

    ! Initialize the module.
    call Star_Formation_Expulsive_Feedback_Disks_Initialize

    ! Get the energy using the selected method.
    Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate=Star_Formation_Expulsive_Feedback_Disk_Rate_Get(thisNode,starFormationRate,energyInputRate)

    return
  end function Star_Formation_Expulsive_Feedback_Disk_Outflow_Rate
  
end module Star_Formation_Feedback_Expulsion_Disks
