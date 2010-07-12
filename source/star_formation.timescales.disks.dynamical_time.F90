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


!% Contains a module which implements a a dynamical time-based star formation timescale for galactic disks.

module Star_Formation_Timescale_Disks_Dynamical_Time
  !% Implements a a dynamical time-based star formation timescale for galactic disks.
  use Tree_Nodes
  private
  public :: Star_Formation_Timescale_Disks_Dynamical_Time_Initialize

  ! Parameters of the timescale model.
  double precision :: starFormationDiskEfficiency,starFormationDiskVelocityExponent,starFormationDiskMinimumTimescale
  
contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Dynamical_Time_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Dynamical_Time_Initialize(starFormationTimescaleDisksMethod,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``dynamical time'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: starFormationTimescaleDisksMethod
    procedure(),          pointer, intent(inout) :: Star_Formation_Timescale_Disk_Get
    
    if (starFormationTimescaleDisksMethod == 'dynamical time') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Dynamical_Time
       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationDiskEfficiency</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The efficiency of star formation in disks for the dynamical time method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationDiskEfficiency'      ,starFormationDiskEfficiency      ,defaultValue= 0.01d0)
       !@ <inputParameter>
       !@   <name>starFormationDiskVelocityExponent</name>
       !@   <defaultValue>$-1.5$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity exponent for star formation in disks for the dynamical time method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationDiskVelocityExponent',starFormationDiskVelocityExponent,defaultValue=-1.500d0)
       !@ <inputParameter>
       !@   <name>starFormationDiskMinimumTimescale</name>
       !@   <defaultValue>$10^{-3}$ Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum timescale for star formation in disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationDiskMinimumTimescale',starFormationDiskMinimumTimescale,defaultValue=1.0d-3)
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Dynamical_time_Initialize

  double precision function Star_Formation_Timescale_Disk_Dynamical_Time(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt thisNode}. The timescale is given by
    !% \begin{equation}
    !% \tau_\star = \epsilon_\star^{-1} \tau_{\rm dynamical, disk} \left( {V_{\rm disk} \over 200\hbox{km/s}} \right)^{\alpha_\star},
    !% \end{equation}
    !% where $\epsilon_\star$(={\tt starFormationDiskEfficiency}) is a star formation efficiency and $\alpha_\star$(={\tt
    !% starFormationDiskVelocityExponent}) controls the scaling with velocity. Note that $\tau_{\rm dynamical,disk}=R_{\rm
    !% disk}/V_{\rm disk}$ where the radius and velocity are whatever characteristic values returned by the disk method. This
    !% scaling is functionally similar to that adopted by \cite{cole_hierarchical_2000}, but that they specifically used the
    !% half-mass radius and circular velocity at that radius.
    use Tree_Nodes
    use Numerical_Constants_Units
    use Tree_Node_Methods
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, parameter              :: velocityZeroPoint=200.0d0 ! (km/s)
    double precision                         :: diskVelocity,dynamicalTime

    ! Get disk circular velocity.
    diskVelocity=Tree_Node_Disk_Velocity(thisNode)

    ! Check for zero velocity disk.
    if (diskVelocity <= 0.0d0) then
       Star_Formation_Timescale_Disk_Dynamical_Time=0.0d0 ! No well defined answer in this case.
    else
       ! Get the dynamical time in Gyr.
       dynamicalTime=Mpc_per_km_per_s_To_Gyr*Tree_Node_Disk_Radius(thisNode)/diskVelocity
       
       ! Compute the star formation timescale using a simple scaling factor.
       Star_Formation_Timescale_Disk_Dynamical_Time=max(dynamicalTime*(diskVelocity/velocityZeroPoint)&
            &**starFormationDiskVelocityExponent/starFormationDiskEfficiency,starFormationDiskMinimumTimescale)
    end if
    return
  end function Star_Formation_Timescale_Disk_Dynamical_Time
  
end module Star_Formation_Timescale_Disks_Dynamical_Time
