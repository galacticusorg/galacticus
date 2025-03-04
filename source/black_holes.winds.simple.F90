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
  Implements a simple black hole winds class.
  !!}

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass

  !![
  <blackHoleWind name="blackHoleWindSimple">
   <description>
    A simple black hole winds model.
   </description>
  </blackHoleWind>
  !!]
  type, extends(blackHoleWindClass) :: blackHoleWindSimple
     !!{
     A simple black hole winds model.
     !!}
     private
     class           (blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
     double precision                                       :: efficiencyWind
   contains
     final     ::          simpleDestructor
     procedure :: power => simplePower
  end type blackHoleWindSimple
  
  interface blackHoleWindSimple
     !!{
     Constructors for the {\normalfont \ttfamily simple} black hole winds class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface blackHoleWindSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily simple} black hole winds class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleWindSimple        )                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (blackHoleAccretionRateClass), pointer       :: blackHoleAccretionRate_
    double precision                                             :: efficiencyWind
    
    !![
    <inputParameter>
      <name>efficiencyWind</name>
      <defaultValue>2.2157d-3</defaultValue>
      <description>The efficiency of the black hole accretion-driven wind.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    !!]
    self=blackHoleWindSimple(efficiencyWind,blackHoleAccretionRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(efficiencyWind,blackHoleAccretionRate_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily simple} node operator class.
    !!}
    implicit none
    type            (blackHoleWindSimple        )                        :: self
    class           (blackHoleAccretionRateClass), target, intent(in   ) :: blackHoleAccretionRate_
    double precision                                     , intent(in   ) :: efficiencyWind
    !![
    <constructorAssign variables="efficiencyWind, *blackHoleAccretionRate_"/>
    !!]

    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !!{
    Destructor for the simple black hole winds class.
    !!}
    implicit none
    type(blackHoleWindSimple), intent(inout) :: self
    
    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    !!]
    return
  end subroutine simpleDestructor
  
  double precision function simplePower(self,blackHole) result(power)
    !!{
    Compute the power of a black hole-driven wind that couples to the surrounding galaxy.
    !!}
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class           (blackHoleWindSimple   ), intent(inout) :: self
    class           (nodeComponentBlackHole), intent(inout) :: blackHole
    double precision                                        :: rateAccretionSpheroid, rateAccretionHotHalo, &
         &                                                     rateAccretion

    ! Compute the wind power.
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo)
    rateAccretion=+rateAccretionSpheroid &
         &        +rateAccretionHotHalo
    if (rateAccretion > 0.0d0) then
       power  =+self%efficiencyWind &
            &  *rateAccretion       &
            &  *speedLight**2       &
            &  /kilo      **2
    else
       power  =+0.0d0
    end if
    return
  end function simplePower
