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

!% Contains a module which implements a calculation of hot halo ram pressure timescales based on an estimate of the acceleration
!% due to ram pressure.

module Hot_Halo_Ram_Pressure_Timescales_Ram_Pressure_Accel
  !% Implements a calculation of hot halo ram pressure timescales based on an estimate of the acceleration due to ram pressure.
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Timescales_Ram_Pressure_Accel_Initialize

contains

  !# <hotHaloRamPressureStrippingTimescaleMethod>
  !#  <unitName>Hot_Halo_Ram_Pressure_Timescales_Ram_Pressure_Accel_Initialize</unitName>
  !# </hotHaloRamPressureStrippingTimescaleMethod>
  subroutine Hot_Halo_Ram_Pressure_Timescales_Ram_Pressure_Accel_Initialize(hotHaloRamPressureStrippingTimescaleMethod,Hot_Halo_Ram_Pressure_Timescale_Get)
    !% Initializes the ``ram pressure acceleration'' hot halo ram pressure stripping timescale module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                    ), intent(in   )          :: hotHaloRamPressureStrippingTimescaleMethod
    procedure(Hot_Halo_Ram_Pressure_Timescale_Ram_Pressure_Accel), intent(inout), pointer :: Hot_Halo_Ram_Pressure_Timescale_Get

    if (hotHaloRamPressureStrippingTimescaleMethod == 'ramPressureAcceleration')                      &
         & Hot_Halo_Ram_Pressure_Timescale_Get => Hot_Halo_Ram_Pressure_Timescale_Ram_Pressure_Accel
    return
  end subroutine Hot_Halo_Ram_Pressure_Timescales_Ram_Pressure_Accel_Initialize

  double precision function Hot_Halo_Ram_Pressure_Timescale_Ram_Pressure_Accel(thisNode)
    !% Computes the hot halo ram pressure stripping timescale, based on the acceleration due to ram pressure forces. This
    !% timescale is approximated as \citep{roediger_ram_2007} $\tau \approx \sqrt{2 r_{\rm outer} \Sigma_{\rm outer} / P_{\rm
    !% ram}}$, where $r_{\rm outer}$ is the current outer radius of the hot halo, $\Sigma_{\rm outer}$ is the surface density at
    !% that radius, and $P_{\rm ram}$ is the ram pressure force (per unit area). The surface density is approximated as
    !% $\Sigma_{\rm outer} \approx r_{\rm outer} \rho_{\rm outer}$, where $\rho_{\rm outer}$ is the density at the outer radius.
    use Galacticus_Nodes
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    use Hot_Halo_Density_Profile
    use Hot_Halo_Ram_Pressure_Forces
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    class           (nodeComponentHotHalo)               , pointer :: thisHotHalo
    type            (treeNode            )               , pointer :: hostNode
    double precision                      , parameter              :: timescaleInfinite       =1.0d30
    double precision                      , parameter              :: velocityStrippingMaximum=1.0d01
    double precision                                               :: outerRadius                   , densityAtOuterRadius       , &
         &                                                            forceRamPressure              , surfaceDensityAtOuterRadius

    ! Evaluate surface density and ram pressure force.
    thisHotHalo                 => thisNode   %hotHalo        (                    )
    outerRadius                 =  thisHotHalo%outerRadius    (                    )
    densityAtOuterRadius        =  Hot_Halo_Density           (thisNode,outerRadius)
    forceRamPressure            =  Hot_Halo_Ram_Pressure_Force(thisNode            )
    surfaceDensityAtOuterRadius =  outerRadius*densityAtOuterRadius
    ! Exit with infinite timescale for zero ram pressure force.
    if (forceRamPressure <= 0.0d0) then
       Hot_Halo_Ram_Pressure_Timescale_Ram_Pressure_Accel=timescaleInfinite
       return
    end if
    ! For zero density or radius, return the halo dynamical time (the timescale should be irrelevant in such cases anyway).
    if (outerRadius <= 0.0d0 .or. surfaceDensityAtOuterRadius <= 0.0d0) then
       Hot_Halo_Ram_Pressure_Timescale_Ram_Pressure_Accel=Dark_Matter_Halo_Dynamical_Timescale(thisNode)
       return
    end if
    ! Find the hosting node.
    hostNode => thisNode%parent
    ! Evaluate the timescale.
    Hot_Halo_Ram_Pressure_Timescale_Ram_Pressure_Accel         &
         & =(                                                  &
         &    megaParsec                                       &
         &   /kilo                                             &
         &   /gigaYear                                         &
         &  )                                                  &
         &  *max(                                              &
         &       sqrt(                                         &
         &            2.0d0                                    &
         &            *outerRadius                             &
         &            *surfaceDensityAtOuterRadius             &
         &            /forceRamPressure                        &
         &           )                                       , &
         &         outerRadius                                 &
         &        /velocityStrippingMaximum                    &
         &        /Dark_Matter_Halo_Virial_Velocity(hostNode)  &
         &      )
    return
  end function Hot_Halo_Ram_Pressure_Timescale_Ram_Pressure_Accel

end module Hot_Halo_Ram_Pressure_Timescales_Ram_Pressure_Accel
