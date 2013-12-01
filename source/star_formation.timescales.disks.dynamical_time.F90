!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a a dynamical time-based star formation timescale for galactic disks.

module Star_Formation_Timescale_Disks_Dynamical_Time
  !% Implements a a dynamical time-based star formation timescale for galactic disks.
  implicit none
  private
  public :: Star_Formation_Timescale_Disks_Dynamical_Time_Initialize

  ! Parameters of the timescale model.
  double precision :: starFormationDiskEfficiency      , starFormationDiskMinimumTimescale, &
       &              starFormationDiskVelocityExponent

contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Dynamical_Time_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Dynamical_Time_Initialize(starFormationTimescaleDisksMethod,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``dynamical time'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Nodes
    use Array_Utilities
    implicit none
    type     (varying_string                               ), intent(in   )          :: starFormationTimescaleDisksMethod
    procedure(Star_Formation_Timescale_Disk_Dynamical_Time ), intent(inout), pointer :: Star_Formation_Timescale_Disk_Get

    if (starFormationTimescaleDisksMethod == 'dynamicalTime') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Dynamical_Time
       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationDiskEfficiency</name>
       !@   <defaultValue>0.01</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The efficiency of star formation in disks for the dynamical time method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationDiskEfficiency'      ,starFormationDiskEfficiency      ,defaultValue= 0.01d0)
       !@ <inputParameter>
       !@   <name>starFormationDiskVelocityExponent</name>
       !@   <defaultValue>$-1.5$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The velocity exponent for star formation in disks for the dynamical time method.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationDiskVelocityExponent',starFormationDiskVelocityExponent,defaultValue=-1.500d0)
       !@ <inputParameter>
       !@   <name>starFormationDiskMinimumTimescale</name>
       !@   <defaultValue>$10^{-3}$ Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum timescale for star formation in disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationDiskMinimumTimescale',starFormationDiskMinimumTimescale,defaultValue=1.0d-3)
       ! Check that required properties are gettable.
       if     (                                                                                                       &
            &  .not.(                                                                                                 &
            &         defaultDiskComponent%velocityIsGettable()                                                       &
            &        .and.                                                                                            &
            &         defaultDiskComponent%  radiusIsGettable()                                                       &
            &       )                                                                                                 &
            & ) call Galacticus_Error_Report                                                                          &
            &        (                                                                                                &
            &         'Star_Formation_Timescale_Disks_Dynamical_Time_Initialize'                                    , &
            &         'disk component must have gettable radius and velocity properties.'//                           &
            &         Galacticus_Component_List(                                                                      &
            &                                   'disk'                                                              , &
            &                                    defaultDiskComponent%velocityAttributeMatch(requireGettable=.true.)  &
            &                                   .intersection.                                                        &
            &                                    defaultDiskComponent%  radiusAttributeMatch(requireGettable=.true.)  &
            &                                  )                                                                      &
            &        )
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Dynamical_Time_Initialize

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
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode         ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk)               , pointer :: thisDiskComponent
    double precision                   , parameter              :: velocityZeroPoint=200.0d0                !   (km/s)
    double precision                                            :: diskVelocity             , dynamicalTime

    ! Get the disk.
    thisDiskComponent => thisNode%disk()

    ! Get disk circular velocity.
    diskVelocity=thisDiskComponent%velocity()

    ! Check for zero velocity disk.
    if (diskVelocity <= 0.0d0) then
       Star_Formation_Timescale_Disk_Dynamical_Time=0.0d0 ! No well defined answer in this case.
    else if (starFormationDiskEfficiency == 0.0d0) then
       ! No star formation occurs if the efficiency is zero.
       Star_Formation_Timescale_Disk_Dynamical_Time=0.0d0
    else
       ! Get the dynamical time in Gyr.
       dynamicalTime=Mpc_per_km_per_s_To_Gyr*thisDiskComponent%radius()/diskVelocity

       ! Compute the star formation timescale using a simple scaling factor.
       Star_Formation_Timescale_Disk_Dynamical_Time=max(dynamicalTime*(diskVelocity/velocityZeroPoint)&
            &**starFormationDiskVelocityExponent/starFormationDiskEfficiency,starFormationDiskMinimumTimescale)
    end if
    return
  end function Star_Formation_Timescale_Disk_Dynamical_Time

end module Star_Formation_Timescale_Disks_Dynamical_Time
