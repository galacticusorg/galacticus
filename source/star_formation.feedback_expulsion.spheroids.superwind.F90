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


!% Contains a module which implements a ``superwind'' expulsive outflow rate (as in \citep{baugh_can_2005}) due to star
!% formation feedback in galactic spheroids.

module Star_Formation_Expulsive_Feedback_Spheroids_Superwind
  !% Implements a``superwind'' outflow rate (as in \citep{baugh_can_2005}) due to star formation feedback in
  !% galactic spheroids.
  use Tree_Nodes
  implicit none
  private
  public :: Star_Formation_Expulsive_Feedback_Spheroids_SW_Initialize

  ! Parameters of the feedback model.
  double precision :: spheroidSuperwindMassLoading,spheroidSuperwindVelocity
  
contains

  !# <starFormationExpulsiveFeedbackSpheroidsMethod>
  !#  <unitName>Star_Formation_Expulsive_Feedback_Spheroids_SW_Initialize</unitName>
  !# </starFormationExpulsiveFeedbackSpheroidsMethod>
  subroutine Star_Formation_Expulsive_Feedback_Spheroids_SW_Initialize(starFormationExpulsiveFeedbackSpheroidsMethod&
       &,Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Get)
    !% Initializes the ``superwind'' spheroid star formation expulsive feedback module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),                 intent(in)    :: starFormationExpulsiveFeedbackSpheroidsMethod
    procedure(double precision), pointer, intent(inout) :: Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Get
    
    if (starFormationExpulsiveFeedbackSpheroidsMethod == 'superwind') then
       Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_Get => Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_SW
       ! Get parameters of for the feedback calculation.
       !@ <inputParameter>
       !@   <name>spheroidSuperwindMassLoading</name>
       !@   <defaultValue>2</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass loading of the spheroid superwind.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidSuperwindMassLoading',spheroidSuperwindMassLoading,defaultValue=2.0d0)
       !@ <inputParameter>
       !@   <name>spheroidSuperwindVelocity</name>
       !@   <defaultValue>200 km/s</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity scale of the spheroid superwind.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidSuperwindVelocity',spheroidSuperwindVelocity,defaultValue=200.0d0)
    end if
    return
  end subroutine Star_Formation_Expulsive_Feedback_Spheroids_SW_Initialize

  double precision function Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_SW(thisNode,starFormationRate,energyInputRate)
    !% Returns the expulsive outflow rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic spheroid of {\tt thisNode}. The outflow
    !% rate is given by
    !% \begin{equation}
    !% \dot{M}_{\rm outflow} = f_{\rm SW,0} \left\{ \begin{array}{ll} 1 & \hbox{ if } V_{\rm spheroid} < V_{\rm spheroid,SW} \\ (V_{\rm spheroid,SW}/V_{\rm spheroid})^2 &  \hbox{ if } V_{\rm spheroid} \ge V_{\rm spheroid,SW} \end{array} \right. ,
    !% \end{equation}
    !%  where $V_{\rm spheroid,SW}=${\tt
    !% [spheroidSuperwindVelocity]} and $f_{\rm SW,0}=${\tt [spheroidSuperwindMassLoading]}. Note that the velocity $V_{\rm
    !% spheroid}$ is whatever characteristic value returned by the spheroid method. This scaling is functionally similar to
    !% that adopted by \cite{cole_hierarchical_2000} and \cite{baugh_can_2005}, except that they specifically used the
    !% circular velocity at half-mass radius.
    use Tree_Nodes
    use Numerical_Constants_Units
    use Stellar_Feedback
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: starFormationRate,energyInputRate
    double precision                         :: spheroidVelocity,outflowRateToStarFormationRate

    ! Get spheroid circular velocity.
    spheroidVelocity=Tree_Node_Spheroid_Velocity(thisNode)

    ! Check for zero velocity spheroid.
    if (spheroidVelocity <= 0.0d0) then
       Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_SW=0.0d0 ! No well defined answer in this case.
    else
       outflowRateToStarFormationRate=spheroidSuperwindMassLoading*(spheroidSuperwindVelocity/max(spheroidVelocity,spheroidSuperwindVelocity))**2
       Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_SW=outflowRateToStarFormationRate*energyInputRate &
            &/feedbackEnergyInputAtInfinityCanonical
    end if
    return
  end function Star_Formation_Expulsive_Feedback_Spheroid_Outflow_Rate_SW
  
end module Star_Formation_Expulsive_Feedback_Spheroids_Superwind
