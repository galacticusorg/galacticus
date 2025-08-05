!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of a spectrum postprocessor that does nothing.
  !!}

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorIdentity">
   <description>
    A stellar population postprocessing class which leaves the spectrum unchanged.
   </description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorIdentity
     !!{
     An identity spectrum postprocessor.
     !!}
     private
   contains
     procedure :: multiplier          => identityMultiplier
     procedure :: isRedshiftDependent => identityIsRedshiftDependent
  end type stellarPopulationSpectraPostprocessorIdentity

  interface stellarPopulationSpectraPostprocessorIdentity
     !!{
     Constructors for the \refClass{stellarPopulationSpectraPostprocessorIdentity} stellar population spectra postprocessor class.
     !!}
     module procedure identityConstructorParameters
  end interface stellarPopulationSpectraPostprocessorIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSpectraPostprocessorIdentity} stellar population spectra postprocessor class which takes a
    parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(stellarPopulationSpectraPostprocessorIdentity)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=stellarPopulationSpectraPostprocessorIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

  double precision function identityMultiplier(self,wavelength,age,redshift)
    !!{
    Perform an identity postprocessing on a spectrum.
    !!}
    implicit none
    class           (stellarPopulationSpectraPostprocessorIdentity), intent(inout) :: self
    double precision                                               , intent(in   ) :: age , redshift, wavelength
    !$GLC attributes unused :: self, age, redshift, wavelength

    identityMultiplier=1.0d0
    return
  end function identityMultiplier

  logical function identityIsRedshiftDependent(self) result(isRedshiftDependent)
    !!{
    Return false indicating that the postprocessor is redshift independent.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorIdentity), intent(inout) :: self
    !$GLC attributes unused :: self

    isRedshiftDependent=.false.
    return
  end function identityIsRedshiftDependent
