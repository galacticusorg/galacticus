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
  An implementation of a spectrum postprocessor that keeps only populations which remain in their birth clouds according to the
  model of \cite{lacey_unified_2016}.
  !!}

  !![
  <stellarPopulationSpectraPostprocessor name="stellarPopulationSpectraPostprocessorBirthCloudsLacey2016">
   <description>Retains only only populations which remain in their birth clouds according to the model of \cite{lacey_unified_2016}.</description>
  </stellarPopulationSpectraPostprocessor>
  !!]
  type, extends(stellarPopulationSpectraPostprocessorClass) :: stellarPopulationSpectraPostprocessorBirthCloudsLacey2016
     !!{
     An implementation of a spectrum postprocessor that keeps only populations which remain in their birth clouds according to the
     model of \cite{lacey_unified_2016}.
     !!}
     private
     double precision :: timescale
   contains
     procedure :: multiplier          => birthCloudsLacey2016Multiplier
     procedure :: isRedshiftDependent => birthCloudsLacey2016IsRedshiftDependent
  end type stellarPopulationSpectraPostprocessorBirthCloudsLacey2016

  interface stellarPopulationSpectraPostprocessorBirthCloudsLacey2016
     !!{
     Constructors for the birthCloudsLacey2016 spectrum postprocessor class.
     !!}
     module procedure birthCloudsLacey2016ConstructorParameters
     module procedure birthCloudsLacey2016ConstructorInternal
  end interface stellarPopulationSpectraPostprocessorBirthCloudsLacey2016

contains

  function birthCloudsLacey2016ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarPopulationSpectraPostprocessorBirthCloudsLacey2016} spectrum postprocessor class which accepts a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (stellarPopulationSpectraPostprocessorBirthCloudsLacey2016)                :: self
    type            (inputParameters                                          ), intent(inout) :: parameters
    double precision                                                                           :: timescale

    !![
    <inputParameter>
      <name>timescale</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The timescale for ``escape'' of stellar populations from their birth clouds.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=stellarPopulationSpectraPostprocessorBirthCloudsLacey2016(timescale)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function birthCloudsLacey2016ConstructorParameters

  function birthCloudsLacey2016ConstructorInternal(timescale) result(self)
    !!{
    Internal constructor for the \refClass{stellarPopulationSpectraPostprocessorBirthCloudsLacey2016} spectrum postprocessor class.
    !!}
    implicit none
    type            (stellarPopulationSpectraPostprocessorBirthCloudsLacey2016)                :: self
    double precision                                                           , intent(in   ) :: timescale
    !![
    <constructorAssign variables="timescale"/>
    !!]

    return
  end function birthCloudsLacey2016ConstructorInternal

  double precision function birthCloudsLacey2016Multiplier(self,wavelength,age,redshift)
    !!{
    Perform an birthCloudsLacey2016 postprocessing on a spectrum.
    !!}
    implicit none
    class           (stellarPopulationSpectraPostprocessorBirthCloudsLacey2016), intent(inout) :: self
    double precision                                                           , intent(in   ) :: age       , redshift, &
         &                                                                                        wavelength
    !$GLC attributes unused :: redshift, wavelength

    ! Implement equation (A5) of Lacey et al. (2016).
    if      (age <       self%timescale) then
       birthCloudsLacey2016Multiplier=+1.0d0
    else if (age < 2.0d0*self%timescale) then
       birthCloudsLacey2016Multiplier=+2.0d0-age/self%timescale
    else
       birthCloudsLacey2016Multiplier=+0.0d0
    end if
    return
  end function birthCloudsLacey2016Multiplier

  logical function birthCloudsLacey2016IsRedshiftDependent(self) result(isRedshiftDependent)
    !!{
    Return false indicating that the postprocessor is redshift independent.
    !!}
    implicit none
    class(stellarPopulationSpectraPostprocessorBirthCloudsLacey2016), intent(inout) :: self
    !$GLC attributes unused :: self

    isRedshiftDependent=.false.
    return
  end function birthCloudsLacey2016IsRedshiftDependent
