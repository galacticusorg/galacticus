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
  An implementation of a spectrum postprocessor that keeps only populations in a specified age window.
  !!}
  
  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorAgeWindow">
   <description>
    A stellar population postprocessor class which keeps only emission from populations with ages between {\normalfont \ttfamily [ageMinimum]} and {\normalfont \ttfamily [ageMaximum]}.
   </description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorAgeWindow
     !!{
     An age window spectrum postprocessor.
     !!}
     private
     double precision :: ageMinimum, ageMaximum
   contains
     procedure :: multiplier          => ageWindowMultiplier
     procedure :: isRedshiftDependent => ageWindowIsRedshiftDependent
  end type stellarPopulationSpectraPostprocessorAgeWindow

  interface stellarPopulationSpectraPostprocessorAgeWindow
     !!{
     Constructors for the ageWindow spectrum postprocessor class.
     !!}
     module procedure ageWindowConstructorParameters
     module procedure ageWindowConstructorInternal
  end interface stellarPopulationSpectraPostprocessorAgeWindow

contains

  function ageWindowConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the ageWindow spectrum postprocessor class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarPopulationSpectraPostprocessorAgeWindow)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    double precision                                                                :: ageMinimum, ageMaximum

    !![
    <inputParameter>
      <name>ageMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum age of stellar populations to retain.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ageMaximum</name>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>The maximum age of stellar populations to retain.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarPopulationSpectraPostprocessorAgeWindow(ageMinimum,ageMaximum)
    return
  end function ageWindowConstructorParameters

  function ageWindowConstructorInternal(ageMinimum,ageMaximum) result(self)
    !!{
    Internal constructor for the \refClass{stellarPopulationSpectraPostprocessorAgeWindow} spectrum postprocessor class.
    !!}
    implicit none
    type            (stellarPopulationSpectraPostprocessorAgeWindow)                :: self
    double precision                                                , intent(in   ) :: ageMinimum, ageMaximum
    !![
    <constructorAssign variables="ageMinimum, ageMaximum"/>
    !!]
    
    return
  end function ageWindowConstructorInternal

  double precision function ageWindowMultiplier(self,wavelength,age,redshift)
    !!{
    Perform a ageWindow postprocessing on a spectrum.
    !!}
    implicit none
    class           (stellarPopulationSpectraPostprocessorAgeWindow), intent(inout) :: self
    double precision                                                , intent(in   ) :: age       , redshift, &
         &                                                                             wavelength
    !$GLC attributes unused :: redshift, wavelength

    if     (                        &
         &   age >= self%ageMinimum &
         &  .and.                   &
         &   age <  self%ageMaximum &
         & ) then
       ageWindowMultiplier=1.0d0
    else
       ageWindowMultiplier=0.0d0
    end if
    return
  end function ageWindowMultiplier

  logical function ageWindowIsRedshiftDependent(self) result(isRedshiftDependent)
    !!{
    Return false indicating that the postprocessor is redshift independent.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorAgeWindow), intent(inout) :: self
    !$GLC attributes unused :: self

    isRedshiftDependent=.false.
    return
  end function ageWindowIsRedshiftDependent
