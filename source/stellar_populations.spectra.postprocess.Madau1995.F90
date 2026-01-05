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

  !!{
  An implementation of a spectrum postprocessor that applies the \cite{madau_radiative_1995} calculation of the attenuation of spectra by the intergalactic medium.
  !!}

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorMadau1995">
   <description>
    A stellar population postprocessor class that postprocesses spectra through absorption by the \gls{igm} using the results
    of \cite{madau_radiative_1995}.
   </description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorMadau1995
     !!{
     A spectrum postprocessor implementing the \cite{madau_radiative_1995} calculation of the attenuation of spectra by the intergalactic medium.
     !!}
     private
   contains
     procedure :: multiplier => madau1995Multiplier
  end type stellarPopulationSpectraPostprocessorMadau1995

  interface stellarPopulationSpectraPostprocessorMadau1995
     !!{
     Constructors for the \refClass{stellarPopulationSpectraPostprocessorMadau1995} stellar population spectra postprocessor class.
     !!}
     module procedure madau1995ConstructorParameters
  end interface stellarPopulationSpectraPostprocessorMadau1995

contains

  function madau1995ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSpectraPostprocessorMadau1995} stellar population spectra postprocessor class which takes a
    parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(stellarPopulationSpectraPostprocessorMadau1995)                :: self
    type(inputParameters                               ), intent(inout) :: parameters

    self=stellarPopulationSpectraPostprocessorMadau1995()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function madau1995ConstructorParameters

  double precision function madau1995Multiplier(self,wavelength,age,redshift)
    !!{
    Suppress the Lyman continuum in a spectrum.
    !!}
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    implicit none
    class           (stellarPopulationSpectraPostprocessorMadau1995), intent(inout)           :: self
    double precision                                                , intent(in   )           :: age                                                                                                                           , redshift                        , &
         &                                                                                       wavelength
    double precision                                                , dimension(9), parameter :: opticalDepthLymanLinesCoefficients=[0.00360d0,0.00170d0,0.00120d0,0.00093d0,0.00093d0,0.00093d0,0.00093d0,0.00093d0,0.00093d0]
    double precision                                                , dimension(9)            :: opticalDepthLymanLines
    integer                                                                                   :: iLine
    double precision                                                                          :: continuumFactor                                                                                                               , emissionFactor                  , &
         &                                                                                       opticalDepth                                                                                                                  , wavelengthObservedLymanContinuum
    !$GLC attributes unused :: self, age

    ! Check if this is a zero redshift case.
    if (redshift <= 0.0d0) then
       ! It is, so return no modification.
       madau1995Multiplier=1.0d0
       return
    else
       ! Compute the observed wavelength in units of the Lyman-continuum wavelength.
       wavelengthObservedLymanContinuum=wavelength*(1.0d0+redshift)/lymanSeriesLimitWavelengthHydrogen_atomic
       ! Compute contribution to optical depth from Lyman-series lines.
       forall(iLine=1:9)
          opticalDepthLymanLines(iLine)=opticalDepthLymanLinesCoefficients(iLine)*(wavelengthObservedLymanContinuum*(1.0d0-1.0d0&
               &/dble((iLine+1)**2)))**3.46d0
       end forall
       ! Cumulate optical depth from lines.
       opticalDepth=0.0d0
       do iLine=1,9
          if (wavelengthObservedLymanContinuum < (1.0d0+redshift)/(1.0d0-1.0d0/dble((iLine+1)**2))) opticalDepth=opticalDepth&
               &+opticalDepthLymanLines(iLine)
       end do
       ! Add in continuum optical depth.
       if (wavelengthObservedLymanContinuum < 1.0d0+redshift) then
          emissionFactor =1.0d0+redshift
          continuumFactor=wavelengthObservedLymanContinuum
          opticalDepth=opticalDepth+0.25d0*(continuumFactor**3)*(emissionFactor**0.46d0-continuumFactor**0.46d0)+9.4d0&
               &*(continuumFactor**1.5d0)*(emissionFactor**0.18d0-continuumFactor**0.18d0)-0.7d0*(continuumFactor**3)*(1.0d0&
               &/continuumFactor**1.32d0-1.0d0/emissionFactor**1.32d0)-0.023d0*(emissionFactor**1.68d0-continuumFactor**1.68d0)
       end if
       ! Compute attenuation from optical depth.
       madau1995Multiplier=exp(-opticalDepth)
    end if
    return
  end function madau1995Multiplier

