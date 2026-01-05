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
  An implementation of a spectrum postprocessor that suppresses the Lyman continuum.
  !!}

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorLycSuppress">
   <description>
    A stellar population spectrum postprocessor class that suppresses all emission in the Lyman continuum.
   </description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorLycSuppress
     !!{
     A stellar population spectrum postprocessor which completely suppresses the Lyman continuum.
     !!}
     private
   contains
     procedure :: multiplier          => lycSuppressMultiplier
     procedure :: isRedshiftDependent => lycSuppressIsRedshiftDependent
  end type stellarPopulationSpectraPostprocessorLycSuppress

  interface stellarPopulationSpectraPostprocessorLycSuppress
     !!{
     Constructors for the \refClass{stellarPopulationSpectraPostprocessorLycSuppress} stellar population spectra postprocessor class.
     !!}
     module procedure lycSuppressConstructorParameters
  end interface stellarPopulationSpectraPostprocessorLycSuppress

contains

  function lycSuppressConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSpectraPostprocessorLycSuppress} stellar population spectra postprocessor class which takes a
    parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(stellarPopulationSpectraPostprocessorLycSuppress)                :: self
    type(inputParameters                                 ), intent(inout) :: parameters

    self=stellarPopulationSpectraPostprocessorLycSuppress()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function lycSuppressConstructorParameters

  double precision function lycSuppressMultiplier(self,wavelength,age,redshift)
    !!{
    Suppress the Lyman continuum in a spectrum.
    !!}
    use :: Numerical_Constants_Atomic, only : lymanSeriesLimitWavelengthHydrogen_atomic
    implicit none
    class           (stellarPopulationSpectraPostprocessorLycSuppress), intent(inout) :: self
    double precision                                                  , intent(in   ) :: age , redshift, wavelength
    !$GLC attributes unused :: self, age, redshift

    if (wavelength < lymanSeriesLimitWavelengthHydrogen_atomic) then
       lycSuppressMultiplier=0.0d0
    else
       lycSuppressMultiplier=1.0d0
    end if
    return
  end function lycSuppressMultiplier

  logical function lycSuppressIsRedshiftDependent(self) result(isRedshiftDependent)
    !!{
    Return false indicating that the postprocessor is redshift independent.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorLycSuppress), intent(inout) :: self
    !$GLC attributes unused :: self

    isRedshiftDependent=.false.
    return
  end function lycSuppressIsRedshiftDependent
