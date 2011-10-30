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


!% Contains a module which implements a power-law outflow rate due to star formation feedback in galactic disks.

module Star_Formation_Feedback_Disks_Power_Law
  !% Implements a power-law outflow rate due to star formation feedback in galactic disks.
  use Tree_Nodes
  implicit none
  private
  public :: Star_Formation_Feedback_Disks_Power_Law_Initialize

  ! Parameters of the feedback model.
  double precision :: diskOutflowVelocity,diskOutflowExponent
  
contains

  !# <starFormationFeedbackDisksMethod>
  !#  <unitName>Star_Formation_Feedback_Disks_Power_Law_Initialize</unitName>
  !# </starFormationFeedbackDisksMethod>
  subroutine Star_Formation_Feedback_Disks_Power_Law_Initialize(starFormationFeedbackDisksMethod,Star_Formation_Feedback_Disk_Outflow_Rate_Get)
    !% Initializes the ``power law'' disk star formation feedback module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: starFormationFeedbackDisksMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Feedback_Disk_Outflow_Rate_Get
    
    if (starFormationFeedbackDisksMethod == 'powerLaw') then
       Star_Formation_Feedback_Disk_Outflow_Rate_Get => Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law
       ! Get parameters of for the feedback calculation.
       !@ <inputParameter>
       !@   <name>diskOutflowVelocity</name>
       !@   <defaultValue>250</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity scale at which the \SNe-driven outflow rate equals the star formation rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowVelocity',diskOutflowVelocity,defaultValue=250.0d0)
       !@ <inputParameter>
       !@   <name>diskOutflowExponent</name>
       !@   <defaultValue>3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity scaling of the \SNe-driven outflow rate in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('diskOutflowExponent',diskOutflowExponent,defaultValue=  3.0d0)
    end if
    return
  end subroutine Star_Formation_Feedback_Disks_Power_Law_Initialize

  double precision function Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law(thisNode,starFormationRate,energyInputRate)
    !% Returns the outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic disk of {\tt thisNode}. The outflow
    !% rate is given by
    !% \begin{equation}
    !% \dot{M}_{\rm outflow} = \left({V_{\rm disk,outflow} \over V_{\rm disk}}\right)^{\alpha_{\rm disk,outflow}},
    !% \end{equation}
    !% where $V_{\rm disk,outflow}$(={\tt diskOutflowVelocity}) is the velocity scale at which outflow rate equals star formation
    !% rate and $\alpha_{\rm disk,outflow}$(={\tt diskOutflowExponent}) controls the scaling with velocity. Note that the velocity
    !% $V_{\rm disk}$ is whatever characteristic value returned by the disk method. This scaling is functionally similar to that
    !% adopted by \cite{cole_hierarchical_2000}, but that they specifically used the circular velocity at half-mass radius.
    use Tree_Nodes
    use Numerical_Constants_Units
    use Stellar_Feedback
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: starFormationRate,energyInputRate
    double precision                         :: diskVelocity,outflowRateToStarFormationRate

    ! Get disk circular velocity.
    diskVelocity=Tree_Node_Disk_Velocity(thisNode)

    ! Check for zero velocity disk.
    if (diskVelocity <= 0.0d0) then
       Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law=0.0d0 ! No well defined answer in this case.
    else
       outflowRateToStarFormationRate=(diskOutflowVelocity/diskVelocity)**diskOutflowExponent
       Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law=outflowRateToStarFormationRate*energyInputRate&
            &/feedbackEnergyInputAtInfinityCanonical
    end if
    return
  end function Star_Formation_Feedback_Disk_Outflow_Rate_Power_Law
  
end module Star_Formation_Feedback_Disks_Power_Law
