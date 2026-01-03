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
  An implementation of a spectrum postprocessor that keeps only unescaped populations.
  !!}

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorUnescaped">
   <description>Retains only unescaped stellar populations.</description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorUnescaped
     !!{
     An unescaped spectrum postprocessor.
     !!}
     private
     double precision :: timescale
   contains
     procedure :: multiplier          => unescapedMultiplier
     procedure :: isRedshiftDependent => unescapedIsRedshiftDependent
  end type stellarPopulationSpectraPostprocessorUnescaped

  interface stellarPopulationSpectraPostprocessorUnescaped
     !!{
     Constructors for the unescaped spectrum postprocessor class.
     !!}
     module procedure unescapedConstructorParameters
     module procedure unescapedConstructorInternal
  end interface stellarPopulationSpectraPostprocessorUnescaped

contains

  function unescapedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the unescaped spectrum postprocessor class which accepts a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarPopulationSpectraPostprocessorUnescaped)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    double precision                                                                :: timescale

    !![
    <inputParameter>
      <name>timescale</name>
      <defaultValue>1.0d-2</defaultValue>
      <description>The timescale for ``escape'' of stellar populations in the ``unescaped'' spectra postprocessing method.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarPopulationSpectraPostprocessorUnescaped(timescale)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function unescapedConstructorParameters

  function unescapedConstructorInternal(timescale) result(self)
    !!{
    Internal constructor for the unescaped spectrum postprocessor class.
    !!}
    implicit none
    type            (stellarPopulationSpectraPostprocessorUnescaped)                :: self
    double precision                                                , intent(in   ) :: timescale
    !![
    <constructorAssign variables="timescale"/>
    !!]

    return
  end function unescapedConstructorInternal

  double precision function unescapedMultiplier(self,wavelength,age,redshift)
    !!{
    Perform an unescaped postprocessing on a spectrum.
    !!}
    implicit none
    class           (stellarPopulationSpectraPostprocessorUnescaped), intent(inout) :: self
    double precision                                                , intent(in   ) :: age       , redshift, &
         &                                                                             wavelength
    !$GLC attributes unused :: redshift, wavelength

    unescapedMultiplier=exp(-age/self%timescale)
    return
  end function unescapedMultiplier

  logical function unescapedIsRedshiftDependent(self) result(isRedshiftDependent)
    !!{
    Return false indicating that the postprocessor is redshift independent.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorUnescaped), intent(inout) :: self
    !$GLC attributes unused :: self

    isRedshiftDependent=.false.
    return
  end function unescapedIsRedshiftDependent
