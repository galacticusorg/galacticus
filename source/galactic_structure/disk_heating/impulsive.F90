!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implementation of a galactic disk heating class which heats disks by tidal shocks in the impulse approximation.
  !!}

  !![
  <galacticDiskHeating name="galacticDiskHeatingImpulsive" docformat="rst">
    <description>
    A galactic disk heating class which heats disks by tidal shocks in the impulse approximation, following
    :cite:t:`carr_migration_2022`. A tidal shock, described by the time integral of the tidal tensor, :math:`G_{ij}`, imparts a
    velocity :math:`\Delta \mathbf{v} = -\mathbf{G} \cdot \mathbf{x}` to a star at position :math:`\mathbf{x}`. For a star at
    radius :math:`R` and azimuth :math:`\phi` in the plane of the disk (taken to be the :math:`x`-:math:`y` plane) this changes
    its specific angular momentum by

    .. math::

       \Delta j = -R^2 \left[ G_{xy} \cos 2\phi + \frac{G_{yy}-G_{xx}}{2} \sin 2\phi \right],

    which averages to zero around every ring of the disk. To first order the disk therefore retains its angular momentum, while
    the guiding radii of its stars are spread out. Averaging over each ring, and over random orientations of the disk relative
    to the tidal field, gives :math:`\langle \Delta j^2 \rangle = (3/5) Q R^4` for a traceless tidal tensor (the coefficient
    becomes approximately :math:`0.56` for the trace typical of the tidal field of a host halo), where :math:`Q` is the
    ``tidalHeatingNormalized`` property of the satellite (such that the mean specific energy injected at radius :math:`r` is
    :math:`Q r^2`; see :galacticus-class:`nodeOperatorSatelliteTidalHeating`). Assuming the self-gravity of the disk to be
    negligible, and its rotation curve to be flat, :math:`\Delta j = V \Delta R`, so a ring at radius :math:`R` is spread in
    radius with a dispersion :math:`\sigma_R = \sqrt{3Q/5} R^2/V`. The heating parameter is therefore

    .. math::

       s = \sqrt{\frac{3Q}{5}} \frac{R_\mathrm{d}}{V_\mathrm{d}},

    where :math:`R_\mathrm{d}` is the scale length of the disk, and :math:`V_\mathrm{d}` is its characteristic circular
    velocity (the ``velocity`` property of the disk component, which for the standard disk is the circular velocity at the radius
    used by the galactic structure solver, :math:`R_\mathrm{d}` for the default settings and an untruncated disk). Only this
    spreading of the disk in radius is modeled. The vertical heating (thickening) of the disk is not.

    This requires that the tidal heating of satellites is evolved (by :galacticus-class:`nodeOperatorSatelliteTidalHeating`),
    as otherwise :math:`Q` remains zero. An error is reported if it is not.
    </description>
  </galacticDiskHeating>
  !!]
  type, extends(galacticDiskHeatingClass) :: galacticDiskHeatingImpulsive
     !!{RST
     A galactic disk heating class which heats disks by tidal shocks in the impulse approximation.
     !!}
     private
   contains
     procedure :: isActive => impulsiveIsActive
     procedure :: heating  => impulsiveHeating
  end type galacticDiskHeatingImpulsive

  interface galacticDiskHeatingImpulsive
     !!{RST
     Constructors for the :galacticus-class:`galacticDiskHeatingImpulsive` galactic disk heating class.
     !!}
     module procedure impulsiveConstructorParameters
  end interface galacticDiskHeatingImpulsive

contains

  function impulsiveConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticDiskHeatingImpulsive` galactic disk heating class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticDiskHeatingImpulsive)                :: self
    type(inputParameters             ), intent(inout) :: parameters

    self=galacticDiskHeatingImpulsive()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function impulsiveConstructorParameters

  logical function impulsiveIsActive(self) result(isActive)
    !!{RST
    Return true, as this class heats disks.
    !!}
    implicit none
    class(galacticDiskHeatingImpulsive), intent(inout) :: self
    !$GLC attributes unused :: self

    isActive=.true.
    return
  end function impulsiveIsActive

  double precision function impulsiveHeating(self,tidalHeatingNormalized,radius,velocity) result(heating)
    !!{RST
    Return the heating parameter, :math:`s=\sqrt{3Q/5} R_\mathrm{d}/V_\mathrm{d}`, of the disk.
    !!}
    use :: Display                       , only : displayGreen                  , displayReset
    use :: Error                         , only : Error_Report
    use :: Satellite_Tidal_Heating_Status, only : satelliteTidalHeatingIsEvolved
    implicit none
    class           (galacticDiskHeatingImpulsive), intent(inout) :: self
    double precision                              , intent(in   ) :: tidalHeatingNormalized, radius, &
         &                                                           velocity
    !$GLC attributes unused :: self

    if (.not.satelliteTidalHeatingIsEvolved)                                                                                               &
         & call Error_Report(                                                                                                              &
         &                   'heating of disks requires that the tidal heating of satellites is evolved'//char(10)//                       &
         &                   displayGreen()//'  HELP:'//displayReset()//' include a `satelliteTidalHeating` node operator in your model'// &
         &                   {introspection:location}                                                                                      &
         &                  )
    if (tidalHeatingNormalized > 0.0d0 .and. radius > 0.0d0 .and. velocity > 0.0d0) then
       heating=+sqrt(                       &
            &        +3.0d0                 &
            &        /5.0d0                 &
            &        *tidalHeatingNormalized&
            &       )                       &
            &  *radius                      &
            &  /velocity
    else
       heating=+0.0d0
    end if
    return
  end function impulsiveHeating
