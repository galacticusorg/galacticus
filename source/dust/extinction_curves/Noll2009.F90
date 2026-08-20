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
  Implements the modified starburst dust attenuation curve of :cite:t:`noll_analysis_2009`.
  !!}

  !![
  <dustExtinctionCurve name="dustExtinctionCurveNoll2009" docformat="rst">
   <description>
   The modified starburst attenuation curve of :cite:t:`noll_analysis_2009`, which augments the
   :cite:t:`calzetti_dust_2000` curve with a variable-strength :math:`2175\,`Å absorption bump and an adjustable
   overall slope:

   .. math::

      \frac{k(\lambda)}{k_\mathrm{V}} = \left[ \frac{k_\mathrm{Cal}(\lambda)}{k_\mathrm{V}} + \frac{D(\lambda)}{R_\mathrm{V}} \right] \left( \frac{\lambda}{\lambda_\mathrm{V}} \right)^\delta,

   where :math:`D(\lambda)` is a Drude profile of amplitude ``bumpStrength``,

   .. math::

      D(\lambda) = \frac{E_\mathrm{b} (\lambda \Delta\lambda)^2}{(\lambda^2-\lambda_0^2)^2 + (\lambda \Delta\lambda)^2},

   with :math:`\lambda_0=2175\,`Å and :math:`\Delta\lambda=350\,`Å. Setting ``bumpStrength`` and ``slope`` to zero
   recovers :cite:t:`calzetti_dust_2000` exactly. This two-parameter family is the one most commonly varied when
   fitting attenuation curves to observed spectral energy distributions.

   Because the underlying :cite:t:`calzetti_dust_2000` curve returns zero outside :math:`0.12`--:math:`2.20\,\mu`m,
   so does this one.
   </description>
  </dustExtinctionCurve>
  !!]
  type, extends(dustExtinctionCurveClass) :: dustExtinctionCurveNoll2009
     !!{RST
     The modified starburst attenuation curve of :cite:t:`noll_analysis_2009`.
     !!}
     private
     type            (dustExtinctionCurveCalzetti2000) :: curveCalzetti2000
     double precision                                  :: bumpStrength     , slope
   contains
     procedure :: attenuationRelative => noll2009AttenuationRelative
     procedure :: wavelengthRange     => noll2009WavelengthRange
  end type dustExtinctionCurveNoll2009

  interface dustExtinctionCurveNoll2009
     !!{RST
     Constructors for the :galacticus-class:`dustExtinctionCurveNoll2009` dust extinction curve class.
     !!}
     module procedure noll2009ConstructorParameters
     module procedure noll2009ConstructorInternal
  end interface dustExtinctionCurveNoll2009

  ! Parameters of the Drude profile describing the 2175Å bump (Noll et al. 2009; eqn. 3), in Angstroms.
  double precision, parameter :: noll2009WavelengthBump=2175.0d0
  double precision, parameter :: noll2009WidthBump     = 350.0d0

contains

  function noll2009ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustExtinctionCurveNoll2009` dust extinction curve class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustExtinctionCurveNoll2009)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    double precision                                             :: bumpStrength, slope

    !![
    <inputParameter docformat="rst">
      <name>bumpStrength</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The amplitude :math:`E_\mathrm{b}` of the :math:`2175\,`Å bump in the :cite:t:`noll_analysis_2009` attenuation
      curve. Zero gives no bump, as in :cite:t:`calzetti_dust_2000`; the Milky Way value is approximately 3.5.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>slope</name>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The slope :math:`\delta` modifying the overall wavelength dependence of the :cite:t:`noll_analysis_2009`
      attenuation curve. Zero recovers the :cite:t:`calzetti_dust_2000` slope; negative values give a steeper curve.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustExtinctionCurveNoll2009(bumpStrength,slope)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function noll2009ConstructorParameters

  function noll2009ConstructorInternal(bumpStrength,slope) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustExtinctionCurveNoll2009` dust extinction curve class.
    !!}
    implicit none
    type            (dustExtinctionCurveNoll2009)                :: self
    double precision                             , intent(in   ) :: bumpStrength, slope
    !![
    <constructorAssign variables="bumpStrength, slope"/>
    !!]

    self%curveCalzetti2000=dustExtinctionCurveCalzetti2000()
    return
  end function noll2009ConstructorInternal

  double precision function noll2009AttenuationRelative(self,wavelength) result(attenuationRelative)
    !!{RST
    Return the relative dust opacity for the attenuation curve of :cite:t:`noll_analysis_2009`.
    !!}
    implicit none
    class           (dustExtinctionCurveNoll2009), intent(inout) :: self
    double precision                             , intent(in   ) :: wavelength
    double precision                                             :: drude

    ! Eqn. (3) of Noll et al. (2009). The Drude profile is expressed in units of k, so is divided by Rv to bring it
    ! into the same relative normalization as the Calzetti curve.
    drude              =+self%bumpStrength                              &
         &              *(wavelength*noll2009WidthBump)**2              &
         &              /(                                              &
         &                +(wavelength**2-noll2009WavelengthBump**2)**2 &
         &                +(wavelength   *noll2009WidthBump        )**2 &
         &               )
    ! Eqn. (5) of Noll et al. (2009).
    attenuationRelative=+(                                                        &
         &                +self%curveCalzetti2000%attenuationRelative(wavelength) &
         &                +drude                                                  &
         &                /self%curveCalzetti2000%RvValue()                       &
         &               )                                                        &
         &              *(                                                        &
         &                +wavelength                                             &
         &                /wavelengthVBand                                        &
         &               )**self%slope
    return
  end function noll2009AttenuationRelative

  subroutine noll2009WavelengthRange(self,wavelengthMinimum,wavelengthMaximum)
    !!{RST
    Return the range of wavelengths over which the :cite:t:`noll_analysis_2009` attenuation curve is defined, which is
    that of the underlying :cite:t:`calzetti_dust_2000` curve.
    !!}
    implicit none
    class           (dustExtinctionCurveNoll2009), intent(inout) :: self
    double precision                             , intent(  out) :: wavelengthMinimum, wavelengthMaximum

    call self%curveCalzetti2000%wavelengthRange(wavelengthMinimum,wavelengthMaximum)
    return
  end subroutine noll2009WavelengthRange
