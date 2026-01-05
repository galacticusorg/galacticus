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
  Implements a class for the cosmic microwave background radiation field.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <radiationField name="radiationFieldCosmicMicrowaveBackground">
   <description>A radiation field class for the cosmic microwave background.</description>
  </radiationField>
  !!]
  type, extends(radiationFieldBlackBody) :: radiationFieldCosmicMicrowaveBackground
     !!{
     A radiation field class for the cosmic microwave background.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: time_
   contains
     final     ::                      cosmicMicrowaveBackgroundDestructor
     procedure :: time              => cosmicMicrowaveBackgroundTime
     procedure :: timeSet           => cosmicMicrowaveBackgroundTimeSet
     procedure :: timeDependentOnly => cosmicMicrowaveBackgroundTimeDependentOnly
  end type radiationFieldCosmicMicrowaveBackground

  interface radiationFieldCosmicMicrowaveBackground
     !!{
     Constructors for the \refClass{radiationFieldCosmicMicrowaveBackground} radiation field class.
     !!}
     module procedure cosmicMicrowaveBackgroundConstructorParameters
     module procedure cosmicMicrowaveBackgroundConstructorInternal
  end interface radiationFieldCosmicMicrowaveBackground

contains

  function cosmicMicrowaveBackgroundConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiationFieldCosmicMicrowaveBackground} radiation field class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (radiationFieldCosmicMicrowaveBackground)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function cosmicMicrowaveBackgroundConstructorParameters

  function cosmicMicrowaveBackgroundConstructorInternal(cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{radiationFieldCosmicMicrowaveBackground} radiation field class.
    !!}
    implicit none
    type (radiationFieldCosmicMicrowaveBackground)                        :: self
    class(cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    self%temperature_=-1.0d0 ! Initialize to an unphysical value.
    return
  end function cosmicMicrowaveBackgroundConstructorInternal

  subroutine cosmicMicrowaveBackgroundDestructor(self)
    !!{
    Destructor for the \refClass{radiationFieldCosmicMicrowaveBackground} radiation field class.
    !!}
    implicit none
    type(radiationFieldCosmicMicrowaveBackground), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine cosmicMicrowaveBackgroundDestructor

  double precision function cosmicMicrowaveBackgroundTime(self)
    !!{
    Return the time for which this radiation field is set.
    !!}
    implicit none
    class(radiationFieldCosmicMicrowaveBackground), intent(inout) :: self

    cosmicMicrowaveBackgroundTime=self%time_
    return
  end function cosmicMicrowaveBackgroundTime

  subroutine cosmicMicrowaveBackgroundTimeSet(self,time)
    !!{
    Set the time (and temperature) of the cosmic microwave background radiation field.
    !!}
    implicit none
    class           (radiationFieldCosmicMicrowaveBackground), intent(inout) :: self
    double precision                                         , intent(in   ) :: time

    self%time_       =                                                    time
    self%temperature_=self%cosmologyFunctions_%temperatureCMBEpochal(time=time)
    return
  end subroutine cosmicMicrowaveBackgroundTimeSet

  logical function cosmicMicrowaveBackgroundTimeDependentOnly(self)
    !!{
    Return true as this radiation field depends on time only.
    !!}
    implicit none
    class(radiationFieldCosmicMicrowaveBackground), intent(inout) :: self
    !$GLC attributes unused :: self

    cosmicMicrowaveBackgroundTimeDependentOnly=.true.
    return
  end function cosmicMicrowaveBackgroundTimeDependentOnly
