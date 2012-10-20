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

!% Contains a module which implements a ``superwind'' expulsive outflow rate (as in \citep{baugh_can_2005}) due to star
!% formation feedback in galactic spheroids.

module Star_Formation_Expulsive_Feedback_Spheroids_Superwind
  !% Implements a``superwind'' outflow rate (as in \citep{baugh_can_2005}) due to star formation feedback in
  !% galactic spheroids.
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
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('spheroidSuperwindMassLoading',spheroidSuperwindMassLoading,defaultValue=2.0d0)
       !@ <inputParameter>
       !@   <name>spheroidSuperwindVelocity</name>
       !@   <defaultValue>200 km/s</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity scale of the spheroid superwind.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
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
