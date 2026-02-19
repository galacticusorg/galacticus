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
  Implements a class that switches on a radiation field at a specified time.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <radiationField name="radiationFieldSwitchOn">
   <description>
    A radiation field class that switches on another radiation field at a specified time.
   </description>
  </radiationField>
  !!]
  type, extends(radiationFieldIntergalacticBackground) :: radiationFieldSwitchOn
     !!{
     A radiation field class that switches on another radiation field at a specified time.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     class           (radiationFieldClass    ), pointer :: radiationField_     => null()
     double precision                                   :: redshiftSwitchOn             , timeSwitchOn, &
          &                                                time_
   contains
     final     ::                      switchOnDestructor
     procedure :: flux              => switchOnFlux
     procedure :: time              => switchOnTime
     procedure :: timeSet           => switchOnTimeSet
     procedure :: timeDependentOnly => switchOnTimeDependentOnly
  end type radiationFieldSwitchOn

  interface radiationFieldSwitchOn
     !!{
     Constructors for the \refClass{radiationFieldSwitchOn} radiation field class.
     !!}
     module procedure switchOnConstructorParameters
     module procedure switchOnConstructorInternal
  end interface radiationFieldSwitchOn

contains

  function switchOnConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiationFieldSwitchOn} radiation field class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (radiationFieldSwitchOn )                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass), pointer       :: cosmologyFunctions_
    class           (radiationFieldClass    ), pointer       :: radiationField_
    double precision                                         :: redshiftSwitchOn

    !![
    <inputParameter>
      <name>redshiftSwitchOn</name>
      <description>The redshift at which to switch on the radiation field.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="radiationField"     name="radiationField_"     source="parameters"/>
    !!]
    self=radiationFieldSwitchOn(cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftSwitchOn)),radiationField_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="radiationField_"    />
    !!]
    return
  end function switchOnConstructorParameters

  function switchOnConstructorInternal(timeSwitchOn,radiationField_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{radiationFieldSwitchOn} radiation field class.
    !!}
    implicit none
    type            (radiationFieldSwitchOn )                        :: self
    double precision                         , intent(in   )         :: timeSwitchOn
    class           (cosmologyFunctionsClass), intent(in   ), target :: cosmologyFunctions_
    class           (radiationFieldClass    ), intent(in   ), target :: radiationField_
    !![
    <constructorAssign variables="timeSwitchOn, *radiationField_, *cosmologyFunctions_"/>
    !!]

    self%redshiftSwitchOn=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeSwitchOn))
    return
  end function switchOnConstructorInternal

  subroutine switchOnDestructor(self)
    !!{
    Destructor for the \refClass{radiationFieldSwitchOn} radiation field class.
    !!}
    implicit none
    type(radiationFieldSwitchOn), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%radiationField_"    />
    !!]
    return
  end subroutine switchOnDestructor

  double precision function switchOnFlux(self,wavelength,node)
    !!{
    Return the flux in the radiation field.
    !!}
    implicit none
    class           (radiationFieldSwitchOn), intent(inout)  :: self
    double precision                        , intent(in   )  :: wavelength
    type            (treeNode              ), intent(inout)  :: node

    if (self%time_ < self%timeSwitchOn) then
       switchOnFlux=0.0d0
    else
       switchOnFlux=self%radiationField_%flux(wavelength,node)
    end if
    return
  end function switchOnFlux

  double precision function switchOnTime(self)
    !!{
    Set the time of the radiation field.
    !!}
    implicit none
    class(radiationFieldSwitchOn), intent(inout) :: self

    switchOnTime=self%time_
    return
  end function switchOnTime

  subroutine switchOnTimeSet(self,time)
    !!{
    Set the time of the radiation field.
    !!}
    implicit none
    class           (radiationFieldSwitchOn), intent(inout) :: self
    double precision                        , intent(in   ) :: time

    self%time_=time
    call self%radiationField_%timeSet(time)
    return
  end subroutine switchOnTimeSet

  logical function switchOnTimeDependentOnly(self)
    !!{
    Return whether the radiation field in this class depends only on time.
    !!}
    implicit none
    class(radiationFieldSwitchOn), intent(inout) :: self

    switchOnTimeDependentOnly=self%radiationField_%timeDependentOnly()
    return
  end function switchOnTimeDependentOnly
