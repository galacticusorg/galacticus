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

!% Contains a module which implements simple mass loss rates from disks due to tidal stripping.

module Tidal_Stripping_Mass_Loss_Rate_Disks_Simple
  !% Implements simple mass loss rates from disks due to tidal stripping.
  use Galacticus_Nodes
  implicit none
  private
  public :: Tidal_Stripping_Mass_Loss_Rate_Disks_Simple_Init

  ! Parameter controlling the maximum mass loss fraction per dynamical time.
  double precision :: tidalStrippingMassLossRateDiskSimpleFractionalRateMaximum

contains

  !# <tidalStrippingMassLossRateDisksMethod>
  !#  <unitName>Tidal_Stripping_Mass_Loss_Rate_Disks_Simple_Init</unitName>
  !# </tidalStrippingMassLossRateDisksMethod>
  subroutine Tidal_Stripping_Mass_Loss_Rate_Disks_Simple_Init(tidalStrippingMassLossRateDisksMethod,Tidal_Stripping_Mass_Loss_Rate_Disk_Get)
    !% Initializes the ``simple'' tidal stripping mass loss rate from disks module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                            ), intent(in   )          :: tidalStrippingMassLossRateDisksMethod
    procedure(Tidal_Stripping_Mass_Loss_Rate_Disk_Simple), intent(inout), pointer :: Tidal_Stripping_Mass_Loss_Rate_Disk_Get

    if (tidalStrippingMassLossRateDisksMethod == 'simple') then
       Tidal_Stripping_Mass_Loss_Rate_Disk_Get => Tidal_Stripping_Mass_Loss_Rate_Disk_Simple
       !@ <inputParameter>
       !@   <name>tidalStrippingMassLossRateDiskSimpleFractionalRateMaximum</name>
       !@   <defaultValue>$10$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum fractional mass loss rate per dynamical time in the simple model of mass loss from disks due to tidal stripping.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('tidalStrippingMassLossRateDiskSimpleFractionalRateMaximum',tidalStrippingMassLossRateDiskSimpleFractionalRateMaximum,defaultValue=10.0d0)
    end if
    return
  end subroutine Tidal_Stripping_Mass_Loss_Rate_Disks_Simple_Init

  double precision function Tidal_Stripping_Mass_Loss_Rate_Disk_Simple(thisNode)
    !% Computes the mass loss rate from disks due to tidal stripping assuming a simple model. Specifically, the mass loss
    !% rate is
    !% \begin{equation}
    !% \dot{M} = -\alpha M/\tau_{\rm disk},
    !% \end{equation}
    !% where
    !% \begin{equation}
    !% \alpha = F_{\rm tidal}/F_{\rm gravity},
    !% \end{equation}
    !% $F_{\rm tidal}=\mathcal{F}_{\rm tidal} r_{1/2}$, $\mathcal{F}_{\rm tidal}$ is the tidal field from the host halo (see \S\ref{sec:SatelliteTidalFields}), and
    !% \begin{equation}
    !% F_{\rm gravity} = V_{1/2}^2(r_{1/2})/r_{1/2}
    !% \end{equation}
    !% is the gravitational restoring force in the disk at the half-mass radius, $r_{\rm 1/2}$.
    use Galacticus_Nodes
    use Satellites_Tidal_Fields
    use Galactic_Structure_Options
    use Galactic_Structure_Rotation_Curves
    use Numerical_Constants_Math
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode         ), intent(inout), pointer :: thisNode
    class           (nodeComponentDisk)               , pointer :: thisDisk
    double precision                                            :: forceGravitational, forceTidal, massLossRateFractional, &
         &                                                         radiusHalfMass    , tidalField, timeDynamical         , &
         &                                                         velocityRotation

    ! Assume no mass loss rate due to tidal by default.
    Tidal_Stripping_Mass_Loss_Rate_Disk_Simple=0.0d0
    ! Return immediately if this is not a satellite.
    if (.not.thisNode%isSatellite()) return
    ! Get the disk component.
    thisDisk           => thisNode%disk()
    ! Get the tidal field due to the host halo.
    tidalField         =  Satellite_Tidal_Field(thisNode)
    ! Get the disk half-mass radius.
    radiusHalfMass     =  thisDisk%halfMassRadius()
    ! Get the tidal force exerted at the half-mass radius.
    forceTidal         =  tidalField*radiusHalfMass
    ! Return if the tidal field is compressive.
    if (forceTidal <= 0.0d0) return
    ! Compute the disk rotation curve.
    velocityRotation=Galactic_Structure_Rotation_Curve(                &
         &                                             thisNode      , &
         &                                             radiusHalfMass  &
         &                                            )
    ! Compute the gravitational restoring force in the disk midplane.
    forceGravitational=velocityRotation**2/radiusHalfMass
    ! Return zero rate if the gravitational force is zero.
    if (forceGravitational <= 0.0d0) return
    ! Compute the mass loss fraction per dynamical time.
    if (forceTidal < tidalStrippingMassLossRateDiskSimpleFractionalRateMaximum*forceGravitational) then
       massLossRateFractional=forceTidal/forceGravitational
    else
       massLossRateFractional=tidalStrippingMassLossRateDiskSimpleFractionalRateMaximum
    end if
    ! Compute the dynamical time.
    timeDynamical         =(megaParsec/kilo/gigaYear)*thisDisk%radius()/thisDisk%velocity()
    ! Compute the mass loss rate.
    Tidal_Stripping_Mass_Loss_Rate_Disk_Simple=massLossRateFractional*(thisDisk%massGas()+thisDisk%massStellar())/timeDynamical
    return
  end function Tidal_Stripping_Mass_Loss_Rate_Disk_Simple

end module Tidal_Stripping_Mass_Loss_Rate_Disks_Simple
