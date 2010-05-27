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
       !@   <defaultValue>0.005</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The efficiency of star formation in disks for the dynamical time method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationDiskEfficiency'      ,starFormationDiskEfficiency      ,defaultValue= 0.005d0)
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
