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

!% Contains a module which implements a star formation timescale for galactic disks which scales with halo virial velocity and
!% redshift.

module Star_Formation_Timescale_Disks_Halo_Scaling
  !% Implements a star formation timescale for galactic disks which scales with halo virial velocity and
  !% redshift.
  use Galacticus_Nodes
  use Kind_Numbers
  implicit none
  private
  public :: Star_Formation_Timescale_Disks_Halo_Scaling_Initialize, Star_Formation_Timescale_Disks_Halo_Scaling_Reset

  ! Parameters of the timescale model.
  double precision                 :: starFormationTimescaleDisksHaloScalingRedshiftExponent              , starFormationTimescaleDisksHaloScalingTimescale, &
       &                              starFormationTimescaleDisksHaloScalingVirialVelocityExponent

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8) :: lastUniqueID                                                =-1
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not timescale has already been computed for this node.
  logical                          :: timescaleComputed                                           =.false.
  !$omp threadprivate(timescaleComputed)
  ! Stored values of the timescale.
  double precision                 :: timeScaleStored
  !$omp threadprivate(timescaleStored)
contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Halo_Scaling_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Halo_Scaling_Initialize(starFormationTimescaleDisksMethod,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``halo scaling'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type     (varying_string                            ), intent(in   )          :: starFormationTimescaleDisksMethod
    procedure(Star_Formation_Timescale_Disk_Halo_Scaling), intent(inout), pointer :: Star_Formation_Timescale_Disk_Get

    if (starFormationTimescaleDisksMethod == 'haloScaling') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Halo_Scaling
       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksHaloScalingTimescale</name>
       !@   <defaultValue>1 Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The timescale for star formation in the halo scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksHaloScalingTimescale',starFormationTimescaleDisksHaloScalingTimescale,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksHaloScalingVirialVelocityExponent</name>
       !@   <defaultValue>$0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of virial velocity in the timescale for star formation in the halo scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksHaloScalingVirialVelocityExponent',starFormationTimescaleDisksHaloScalingVirialVelocityExponent,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksHaloScalingRedshiftExponent</name>
       !@   <defaultValue>$0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of redshift in the timescale for star formation in the halo scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksHaloScalingRedshiftExponent',starFormationTimescaleDisksHaloScalingRedshiftExponent,defaultValue=0.0d0)
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Halo_Scaling_Initialize

  double precision function Star_Formation_Timescale_Disk_Halo_Scaling(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt thisNode} in the halo scaling timescale model.
    use Cosmology_Functions
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    double precision                    , parameter              :: virialVelocityNormalization=200.0d0
    double precision                                             :: expansionFactor                    , virialVelocity

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= lastUniqueID) call Star_Formation_Timescale_Disks_Halo_Scaling_Reset(thisNode)

    ! Compute the timescale if necessary.
    if (.not.timescaleComputed) then

       ! Get virial velocity and expansion factor.
       virialVelocity =Dark_Matter_Halo_Virial_Velocity(thisNode                 )
       expansionFactor=Expansion_Factor                (thisBasicComponent%time())

       ! Return the timescale.
       timescaleStored=                                                                                                   &
            &  starFormationTimescaleDisksHaloScalingTimescale                                                            &
            & *(virialVelocity/virialVelocityNormalization)**starFormationTimescaleDisksHaloScalingVirialVelocityExponent &
            & /expansionFactor**starFormationTimescaleDisksHaloScalingRedshiftExponent

       ! Record that the timescale is now computed.
       timescaleComputed=.true.
    end if

    ! Return the stored timescale.
    Star_Formation_Timescale_Disk_Halo_Scaling=timescaleStored
    return
  end function Star_Formation_Timescale_Disk_Halo_Scaling

  !# <calculationResetTask>
  !# <unitName>Star_Formation_Timescale_Disks_Halo_Scaling_Reset</unitName>
  !# </calculationResetTask>
  subroutine Star_Formation_Timescale_Disks_Halo_Scaling_Reset(thisNode)
    !% Reset the halo scaling disk star formation timescale calculation.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    timescaleComputed=.false.
    lastUniqueID     =thisNode%uniqueID()
    return
  end subroutine Star_Formation_Timescale_Disks_Halo_Scaling_Reset

end module Star_Formation_Timescale_Disks_Halo_Scaling
