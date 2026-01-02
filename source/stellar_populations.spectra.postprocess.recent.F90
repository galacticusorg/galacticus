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
  An implementation of a spectrum postprocessor that keeps only recent populations.
  !!}

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorRecent">
   <description>
    A stellar population postprocessor class which suppresses all emission from populations older than {\normalfont \ttfamily
    [timeLimit]} (in Gyr).
   </description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorRecent
     !!{
     An recent spectrum postprocessor.
     !!}
     private
     double precision :: timeLimit
   contains
     procedure :: multiplier          => recentMultiplier
     procedure :: isRedshiftDependent => recentIsRedshiftDependent
  end type stellarPopulationSpectraPostprocessorRecent

  interface stellarPopulationSpectraPostprocessorRecent
     !!{
     Constructors for the recent spectrum postprocessor class.
     !!}
     module procedure recentConstructorParameters
     module procedure recentConstructorInternal
  end interface stellarPopulationSpectraPostprocessorRecent

contains

  function recentConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the recent spectrum postprocessor class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarPopulationSpectraPostprocessorRecent)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    double precision                                                             :: timeLimit

    !![
    <inputParameter>
      <name>timeLimit</name>
      <defaultValue>1.0d-2</defaultValue>
      <description>The maximum age of stellar populations to retain in the ``recent'' spectra postprocessing method.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarPopulationSpectraPostprocessorRecent(timeLimit)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function recentConstructorParameters

  function recentConstructorInternal(timeLimit) result(self)
    !!{
    Internal constructor for the recent spectrum postprocessor class.
    !!}
    implicit none
    type            (stellarPopulationSpectraPostprocessorRecent)                :: self
    double precision                                             , intent(in   ) :: timeLimit
    !![
    <constructorAssign variables="timeLimit"/>
    !!]

    return
  end function recentConstructorInternal

  double precision function recentMultiplier(self,wavelength,age,redshift)
    !!{
    Perform a recent postprocessing on a spectrum.
    !!}
    implicit none
    class           (stellarPopulationSpectraPostprocessorRecent), intent(inout) :: self
    double precision                                             , intent(in   ) :: age       , redshift, &
         &                                                                          wavelength
    !$GLC attributes unused :: redshift, wavelength

    if (age > self%timeLimit) then
       recentMultiplier=0.0d0
    else
       recentMultiplier=1.0d0
    end if
    return
  end function recentMultiplier

  logical function recentIsRedshiftDependent(self) result(isRedshiftDependent)
    !!{
    Return false indicating that the postprocessor is redshift independent.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorRecent), intent(inout) :: self
    !$GLC attributes unused :: self

    isRedshiftDependent=.false.
    return
  end function recentIsRedshiftDependent
