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
  Implements an abstract class of spectra of thermal emission from dust given by a tabulation.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <dustEmissionSpectrum name="dustEmissionSpectrumTabulated" abstract="yes" docformat="rst">
   <description>
   An abstract spectrum of thermal emission from dust, given by a tabulation of its shape against wavelength. Between the
   tabulated wavelengths the spectrum is interpolated as a power law, and it is zero outside the tabulated range. It is
   normalized by integrating that interpolant analytically, so that it integrates to exactly the absorbed luminosity
   whatever wavelengths it is evaluated at. The mass of dust is not used: the shape is fixed by the tabulation.

   Concrete classes supply the tabulation, by calling ``initializeTabulation`` from their constructors.
   </description>
  </dustEmissionSpectrum>
  !!]
  type, abstract, extends(dustEmissionSpectrumClass) :: dustEmissionSpectrumTabulated
     !!{RST
     An abstract spectrum of thermal emission from dust given by a tabulation.
     !!}
     private
     ! The tabulated wavelengths (in Å), and the logarithm of the spectrum at each, normalized to unit luminosity (so in
     ! Hz⁻¹).
     double precision              , allocatable, dimension(:) :: wavelength            , spectrumLogarithmic
     type            (interpolator)                            :: interpolatorWavelength
   contains
     !![
     <methods docformat="rst">
       <method method="initializeTabulation" description="Set the tabulated spectrum, normalizing it to unit luminosity."/>
     </methods>
     !!]
     procedure :: initializeTabulation => tabulatedInitializeTabulation
     procedure :: luminosity           => tabulatedLuminosity
  end type dustEmissionSpectrumTabulated

contains

  subroutine tabulatedInitializeTabulation(self,wavelength,luminosityLogarithmic)
    !!{RST
    Set the tabulated spectrum from ``luminosityLogarithmic``, the logarithm of :math:`\nu L_\nu` (equivalently,
    :math:`\lambda L_\lambda`) in arbitrary units, at each of the ``wavelength`` (in Å). There must be at least two.

    The spectrum is converted to :math:`L_\nu`, which is proportional to :math:`\lambda \nu L_\nu`, and normalized to
    unit luminosity by integrating its power-law interpolant (see ``dustEmissionIntegralPowerLaw``).
    !!}
    use :: Error                       , only : Error_Report
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms
    use :: Table_Labels                , only : extrapolationTypeFix
    implicit none
    class           (dustEmissionSpectrumTabulated), intent(inout)               :: self
    double precision                               , intent(in   ), dimension(:) :: wavelength   , luminosityLogarithmic
    double precision                                                             :: normalization

    if (size(wavelength) < 2                           ) call Error_Report('at least two tabulated wavelengths are required'               //{introspection:location})
    if (size(wavelength) /= size(luminosityLogarithmic)) call Error_Report('`wavelength` and `luminosityLogarithmic` must be the same size'//{introspection:location})
    self%wavelength            =wavelength
    normalization              =dustEmissionIntegralPowerLaw(self%wavelength,luminosityLogarithmic)
    self%spectrumLogarithmic   =+luminosityLogarithmic                             &
         &                      +log(self%wavelength/speedLight/metersToAngstroms) &
         &                      -log(normalization)
    self%interpolatorWavelength=interpolator(log(self%wavelength),extrapolationType=extrapolationTypeFix)
    return
  end subroutine tabulatedInitializeTabulation

  function tabulatedLuminosity(self,wavelengths,luminosityAbsorbed,massDust,time) result(luminosity)
    !!{RST
    Return the luminosity per unit frequency of a tabulated dust emission spectrum, normalized to the luminosity
    absorbed.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (dustEmissionSpectrumTabulated), intent(inout)                               :: self
    double precision                               , intent(in   ), dimension(:                ) :: wavelengths
    double precision                               , intent(in   )                               :: luminosityAbsorbed, massDust, &
         &                                                                                          time
    double precision                                              , dimension(size(wavelengths)) :: luminosity
    double precision                                              , dimension(0:1              ) :: weights
    integer         (c_size_t                     )                                              :: j
    integer                                                                                      :: i
    !$GLC attributes unused :: massDust, time

    luminosity=0.0d0
    if (luminosityAbsorbed <= 0.0d0) return
    do i=1,size(wavelengths)
       if (wavelengths(i) < self%wavelength(1) .or. wavelengths(i) > self%wavelength(size(self%wavelength))) cycle
       call self%interpolatorWavelength%linearFactors(log(wavelengths(i)),j,weights)
       luminosity(i)=+luminosityAbsorbed                            &
            &        *exp(                                          &
            &             +weights(0)*self%spectrumLogarithmic(j  ) &
            &             +weights(1)*self%spectrumLogarithmic(j+1) &
            &            )
    end do
    return
  end function tabulatedLuminosity
