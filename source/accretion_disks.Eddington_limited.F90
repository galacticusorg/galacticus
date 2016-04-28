!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Contains a module which implements calculations of properties of accretion disks which ignore the details of accretion physics
!% and assume that black hole jets have a power that is a fixed fraction of the Eddington luminosity.

module Accretion_Disks_Eddington
  !% Implements calculations of properties of accretion disks which ignore the details of accretion physics
  !% and assume that black hole jets have a power that is a fixed fraction of the Eddington luminosity.
  implicit none
  private
  public :: Accretion_Disks_Eddington_Initialize

  ! Efficiency parameters for the accretion disk model.
  double precision :: accretionDiskJetPowerEddington, accretionDiskRadiativeEfficiencyEddington

contains

  !# <accretionDisksMethod>
  !#  <unitName>Accretion_Disks_Eddington_Initialize</unitName>
  !# </accretionDisksMethod>
  subroutine Accretion_Disks_Eddington_Initialize(accretionDisksMethod,Accretion_Disk_Radiative_Efficiency_Get&
       &,Black_Hole_Spin_Up_Rate_Get,Accretion_Disk_Jet_Power_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                               ), intent(in   )          :: accretionDisksMethod
    procedure(Accretion_Disk_Radiative_Efficiency_Eddington), intent(inout), pointer :: Accretion_Disk_Radiative_Efficiency_Get
    procedure(Black_Hole_Spin_Up_Rate_Eddington            ), intent(inout), pointer :: Black_Hole_Spin_Up_Rate_Get
    procedure(Accretion_Disk_Jet_Power_Eddington           ), intent(inout), pointer :: Accretion_Disk_Jet_Power_Get

    if (accretionDisksMethod == 'eddingtonLimited') then
       Accretion_Disk_Radiative_Efficiency_Get => Accretion_Disk_Radiative_Efficiency_Eddington
       Black_Hole_Spin_Up_Rate_Get             => Black_Hole_Spin_Up_Rate_Eddington
       Accretion_Disk_Jet_Power_Get            => Accretion_Disk_Jet_Power_Eddington
       !@ <inputParameter>
       !@   <name>accretionDiskJetPowerEddington</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The jet power produced by an Eddington limited accretion disk in units of the Eddington luminosity.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionDiskJetPowerEddington",accretionDiskJetPowerEddington,defaultValue=0.1d0)
       !@ <inputParameter>
       !@   <name>accretionDiskRadiativeEfficiencyEddington</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    The radiative efficiency of an Eddington limited accretion disk.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter("accretionDiskRadiativeEfficiencyEddington",accretionDiskRadiativeEfficiencyEddington,defaultValue=0.1d0)
    end if
    return
  end subroutine Accretion_Disks_Eddington_Initialize

  double precision function Accretion_Disk_Radiative_Efficiency_Eddington(thisBlackHole,massAccretionRate)
    !% Computes the radiative efficiency for an Eddington-limited accretion disk.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    !GCC$ attributes unused :: thisBlackHole, massAccretionRate
    
    Accretion_Disk_Radiative_Efficiency_Eddington=accretionDiskRadiativeEfficiencyEddington
    return
  end function Accretion_Disk_Radiative_Efficiency_Eddington

  double precision function Accretion_Disk_Jet_Power_Eddington(blackHole,massAccretionRate)
    !% Computes the jet power for an Eddington-limited accretion disk. The jet power is held at
    !% a fixed fraction of the Eddington accretion rate, independent of the mass accretion rate
    !% onto the black hole, until the mass accretion rate falls below some small fraction of the
    !% Eddington rate, at which point the jet power is smoothly reduced to zero as mass
    !% accretion rate approaches zero. This avoids the unphysical behavior of having finite jet
    !% power at zero accretion rate.
    use Galacticus_Nodes
    use Black_Hole_Fundamentals
    use Numerical_Constants_Physical
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                        , intent(in   ) :: massAccretionRate
    double precision                        , parameter     :: fractionalRateCutOff                    =0.01d0
    double precision                                        :: massAccretionRateEddingtonLimited              , &
         &                                                     massAccretionRateEddingtonLimitedReduced       , &
         &                                                     fractionalRate

    massAccretionRateEddingtonLimited=+accretionDiskJetPowerEddington                 &
         &                            *Black_Hole_Eddington_Accretion_Rate(blackHole)
    if (massAccretionRate > fractionalRateCutOff*massAccretionRateEddingtonLimited) then
       massAccretionRateEddingtonLimitedReduced=massAccretionRateEddingtonLimited
    else
       fractionalRate                          =+massAccretionRate                 &
            &                                   /fractionalRateCutOff              &
            &                                   /massAccretionRateEddingtonLimited
       massAccretionRateEddingtonLimitedReduced=+massAccretionRateEddingtonLimited &
            &                                   *(                                 &
            &                                     +3.0d0*fractionalRate**2         &
            &                                     -2.0d0*fractionalRate**3         &
            &                                   )
    end if
    Accretion_Disk_Jet_Power_Eddington=+massAccretionRateEddingtonLimitedReduced &
         &                             *(                                        &
         &                              +speedLight                              &
         &                              /kilo                                    &
         &                             )**2
    return
  end function Accretion_Disk_Jet_Power_Eddington

  double precision function Black_Hole_Spin_Up_Rate_Eddington(thisBlackHole,massAccretionRate)
    !% Computes the spin up rate of the black hole in {\normalfont \ttfamily thisNode} due to accretion from an Eddington-limited accretion disk.
    !% disk. This is always zero, as no physical model is specified for this accretion disk method.
    use Galacticus_Nodes
    implicit none
    class           (nodeComponentBlackHole), intent(inout) :: thisBlackHole
    double precision                        , intent(in   ) :: massAccretionRate
    !GCC$ attributes unused :: thisBlackHole, massAccretionRate
    
    Black_Hole_Spin_Up_Rate_Eddington=0.0d0
    return
  end function Black_Hole_Spin_Up_Rate_Eddington

end module Accretion_Disks_Eddington
