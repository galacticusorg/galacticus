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
  Implements the dust extinction curve of :cite:t:`calzetti_dust_2000`.
  !!}

  !![
  <dustExtinctionCurve name="dustExtinctionCurveCalzetti2000" docformat="rst">
   <description>
   The dust attenuation curve of :cite:t:`calzetti_dust_2000`, derived empirically from the ultraviolet-to-optical
   spectra of local starburst galaxies. Equation (4) of that work gives :math:`k(\lambda)`, which is divided here by
   :math:`R_\mathrm{V}=4.05` (their equation 5) to normalize to the :math:`V`-band.

   The fit is defined only over :math:`0.12\,\mu\mathrm{m} &lt; \lambda \le 2.20\,\mu\mathrm{m}`. Outside that range this
   implementation returns zero, i.e. *no* attenuation---the behavior of the original ``stellarSpectraDustAttenuation``
   implementation, retained here for consistency. Note that this is not the same as clamping the curve to its boundary
   value, and can matter when the curve is evaluated across a broad spectrum rather than through an optical filter;
   use ``wavelengthRange`` to detect the boundary.
   </description>
  </dustExtinctionCurve>
  !!]
  type, extends(dustExtinctionCurveClass) :: dustExtinctionCurveCalzetti2000
     !!{RST
     The dust extinction curve of :cite:t:`calzetti_dust_2000`.
     !!}
     private
   contains
     !![
     <methods docformat="rst">
       <method method="RvValue" description="Return the total-to-selective extinction ratio assumed by this curve."/>
     </methods>
     !!]
     procedure :: attenuationRelative => calzetti2000AttenuationRelative
     procedure :: wavelengthRange     => calzetti2000WavelengthRange
     procedure :: RvValue             => calzetti2000RvValue
  end type dustExtinctionCurveCalzetti2000

  interface dustExtinctionCurveCalzetti2000
     !!{RST
     Constructors for the :galacticus-class:`dustExtinctionCurveCalzetti2000` dust extinction curve class.
     !!}
     module procedure calzetti2000ConstructorParameters
  end interface dustExtinctionCurveCalzetti2000

  ! Range of validity of the fit, in microns (Calzetti et al. 2000; eqn. 4).
  double precision, parameter :: calzetti2000WavelengthMinimumMicrons=0.12d0
  double precision, parameter :: calzetti2000WavelengthMaximumMicrons=2.20d0
  ! Total-to-selective extinction ratio (Calzetti et al. 2000; eqn. 5).
  double precision, parameter :: calzetti2000Rv                      =4.05d0

contains

  function calzetti2000ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustExtinctionCurveCalzetti2000` dust extinction curve class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(dustExtinctionCurveCalzetti2000)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=dustExtinctionCurveCalzetti2000()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function calzetti2000ConstructorParameters

  double precision function calzetti2000AttenuationRelative(self,wavelength) result(attenuationRelative)
    !!{RST
    Return the relative dust opacity for the extinction curve of :cite:t:`calzetti_dust_2000`.
    !!}
    use :: Numerical_Constants_Units, only : micronsToAngstroms
    implicit none
    class           (dustExtinctionCurveCalzetti2000), intent(inout) :: self
    double precision                                 , intent(in   ) :: wavelength
    double precision                                                 :: wavelengthMicrons, kappa
    !$GLC attributes unused :: self

    ! Eqn. (4) of Calzetti et al. (2000).
    wavelengthMicrons=wavelength/micronsToAngstroms
    if      (wavelengthMicrons > calzetti2000WavelengthMinimumMicrons .and. wavelengthMicrons <= 0.63d0                             ) then
       kappa=+2.659d0                                &
            & *(                                     &
            &   -2.156d0                             &
            &   +1.509d0/wavelengthMicrons           &
            &   -0.198d0/wavelengthMicrons**2        &
            &   +0.011d0/wavelengthMicrons**3        &
            &  )                                     &
            & +calzetti2000Rv
    else if (wavelengthMicrons > 0.63d0                              .and. wavelengthMicrons <= calzetti2000WavelengthMaximumMicrons) then
       kappa=+2.659d0                                &
            & *(                                     &
            &   -1.857d0                             &
            &   +1.040d0/wavelengthMicrons           &
            &  )                                     &
            & +calzetti2000Rv
    else
       ! Outside the range of the fit - see the class description.
       kappa=+0.0d0
    end if
    attenuationRelative=+kappa            &
         &              /calzetti2000Rv
    return
  end function calzetti2000AttenuationRelative

  subroutine calzetti2000WavelengthRange(self,wavelengthMinimum,wavelengthMaximum)
    !!{RST
    Return the range of wavelengths over which the :cite:t:`calzetti_dust_2000` extinction curve is defined.
    !!}
    use :: Numerical_Constants_Units, only : micronsToAngstroms
    implicit none
    class           (dustExtinctionCurveCalzetti2000), intent(inout) :: self
    double precision                                 , intent(  out) :: wavelengthMinimum, wavelengthMaximum
    !$GLC attributes unused :: self

    wavelengthMinimum=calzetti2000WavelengthMinimumMicrons*micronsToAngstroms
    wavelengthMaximum=calzetti2000WavelengthMaximumMicrons*micronsToAngstroms
    return
  end subroutine calzetti2000WavelengthRange

  double precision function calzetti2000RvValue(self) result(Rv)
    !!{RST
    Return the total-to-selective extinction ratio, :math:`R_\mathrm{V}`, assumed by the
    :cite:t:`calzetti_dust_2000` curve. Provided as an accessor because each implementation of this class is generated
    into its own submodule, so a module-level constant declared here is not visible to sibling implementations---and
    :galacticus-class:`dustExtinctionCurveNoll2009`, which is defined as a modification of this curve, needs it.
    !!}
    implicit none
    class(dustExtinctionCurveCalzetti2000), intent(inout) :: self
    !$GLC attributes unused :: self

    Rv=calzetti2000Rv
    return
  end function calzetti2000RvValue
