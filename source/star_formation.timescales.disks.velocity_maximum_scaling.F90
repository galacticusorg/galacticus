! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a star formation timescale for galactic disks which scales with halo maximum velocity and
!% redshift.

module Star_Formation_Timescale_Disks_VlctyMxSclng
  !% Implements a star formation timescale for galactic disks which scales with halo maximum velocity and
  !% redshift.
  use Galacticus_Nodes
  use Kind_Numbers
  implicit none
  private
  public :: Star_Formation_Timescale_Disks_VlctyMxSclng_Initialize, Star_Formation_Timescale_Disks_VlctyMxSclng_Reset

  ! Parameters of the timescale model.
  double precision                             :: starFormationTimescaleDisksVlctyMxSclngRedshiftExponent              , starFormationTimescaleDisksVlctyMxSclngTimescale, &
       &                                          starFormationTimescaleDisksVlctyMxSclngVelocityExponent

  ! Record of unique ID of node which we last computed results for.
  integer         (kind=kind_int8)             :: lastUniqueID                                                =-1
  !$omp threadprivate(lastUniqueID)
  ! Record of whether or not timescale has already been computed for this node.
  logical                                      :: timescaleComputed                                           =.false.
  !$omp threadprivate(timescaleComputed)
  ! Stored values of the timescale.
  double precision                             :: timeScaleStored
  !$omp threadprivate(timescaleStored)
  ! Normalization of the timescale.
  double precision                            :: timeScaleNormalization
  ! Normalization for velocities. 
  double precision                , parameter :: velocityNormalization=200.0d0

contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_VlctyMxSclng_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_VlctyMxSclng_Initialize(starFormationTimescaleDisksMethod,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``halo scaling'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type            (varying_string                                        ), intent(in   )          :: starFormationTimescaleDisksMethod 
    procedure       (Star_Formation_Timescale_Disk_VlctyMxSclng), intent(inout), pointer :: Star_Formation_Timescale_Disk_Get 

    if (starFormationTimescaleDisksMethod == 'velocityMaximumScaling') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_VlctyMxSclng
       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksVelocityMaximumScalingTimescale</name>
       !@   <defaultValue>1 Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The timescale for star formation in the velocity maximum scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksVelocityMaximumScalingTimescale',starFormationTimescaleDisksVlctyMxSclngTimescale,defaultValue=1.0d0)
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksVelocityMaximumScalingVelocityExponent</name>
       !@   <defaultValue>$0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of virial velocity in the timescale for star formation in the velocity maximum scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksVelocityMaximumScalingVelocityExponent',starFormationTimescaleDisksVlctyMxSclngVelocityExponent,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksVelocityMaximumScalingRedshiftExponent</name>
       !@   <defaultValue>$0$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The exponent of redshift in the timescale for star formation in the velocity maximum scaling timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksVelocityMaximumScalingRedshiftExponent',starFormationTimescaleDisksVlctyMxSclngRedshiftExponent,defaultValue=0.0d0)
       ! Compute the normalization of the timescale.
       timeScaleNormalization= starFormationTimescaleDisksVlctyMxSclngTimescale                               &
       &                      /velocityNormalization**starFormationTimescaleDisksVlctyMxSclngVelocityExponent
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_VlctyMxSclng_Initialize

  double precision function Star_Formation_Timescale_Disk_VlctyMxSclng(node)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt node} in the velocity maximum scaling timescale model.
    use Cosmology_Functions
    use Dark_Matter_Profiles
    implicit none
    type            (treeNode               ), intent(inout), pointer :: node
    class           (nodeComponentBasic     )               , pointer :: basic
    class           (cosmologyFunctionsClass)               , pointer :: cosmologyFunctions_
    class           (darkMatterProfileClass )               , pointer :: darkMatterProfile_
    double precision                         , save                   :: velocityMaximumPrevious=-1.0d0, velocityFactorPrevious       =-1.0d0
    !$omp threadprivate(velocityMaximumPrevious,velocityFactorPrevious)
    double precision                         , save                   :: expansionFactorPrevious=-1.0d0, expansionFactorFactorPrevious=-1.0d0
    !$omp threadprivate(expansionFactorPrevious,expansionFactorFactorPrevious)
    double precision                                                  :: expansionFactor               , velocityMaximum

    ! Get the basic component.
    basic => node%basic()
    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= lastUniqueID) call Star_Formation_Timescale_Disks_VlctyMxSclng_Reset(node)
    ! Compute the timescale if necessary.
    if (.not.timescaleComputed) then
       ! Get the default cosmology functions object.
       cosmologyFunctions_ => cosmologyFunctions()
       darkMatterProfile_  => darkMatterProfile ()
       ! Get virial velocity and expansion factor.
       velocityMaximum=darkMatterProfile_ %circularVelocityMaximum(node        )
       expansionFactor=cosmologyFunctions_%expansionFactor        (basic%time())
       ! Compute the velocity factor.
       if (velocityMaximum /= velocityMaximumPrevious) then
          velocityMaximumPrevious      =      velocityMaximum
          velocityFactorPrevious       =      velocityMaximum**starFormationTimescaleDisksVlctyMxSclngVelocityExponent
       end if
       ! Compute the expansion-factor factor.
       if (expansionFactor /= expansionFactorPrevious) then
          expansionFactorPrevious      =      expansionFactor
          expansionFactorFactorPrevious=1.0d0/expansionFactor**starFormationTimescaleDisksVlctyMxSclngRedshiftExponent
       end if
       ! Return the timescale.
       timescaleStored=                      &
            &  timeScaleNormalization        &
            & *velocityFactorPrevious        &
            & *expansionFactorFactorPrevious
       ! Record that the timescale is now computed.
       timescaleComputed=.true.
    end if
    ! Return the stored timescale.
    Star_Formation_Timescale_Disk_VlctyMxSclng=timescaleStored
    return
  end function Star_Formation_Timescale_Disk_VlctyMxSclng

  !# <calculationResetTask>
  !# <unitName>Star_Formation_Timescale_Disks_VlctyMxSclng_Reset</unitName>
  !# </calculationResetTask>
  subroutine Star_Formation_Timescale_Disks_VlctyMxSclng_Reset(node)
    !% Reset the velocity maximum scaling disk star formation timescale calculation.
    implicit none
    type(treeNode), intent(inout), pointer :: node

    timescaleComputed=.false.
    lastUniqueID     =node%uniqueID()
    return
  end subroutine Star_Formation_Timescale_Disks_VlctyMxSclng_Reset

end module Star_Formation_Timescale_Disks_VlctyMxSclng
