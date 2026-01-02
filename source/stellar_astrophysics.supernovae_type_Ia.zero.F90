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
  Implements a supernovae type Ia class with no supernovae.
  !!}

  !![
  <supernovaeTypeIa name="supernovaeTypeIaZero">
   <description>
    A supernovae type Ia class which produces zero supernovae.
   </description>
  </supernovaeTypeIa>
  !!]
  type, extends(supernovaeTypeIaClass) :: supernovaeTypeIaZero
     !!{
     A supernovae type Ia class that produces zero supernovae.
     !!}
     private
   contains
     procedure :: massInitialRange => zeroMassInitialRange
     procedure :: number           => zeroNumber
     procedure :: yield            => zeroYield
  end type supernovaeTypeIaZero

  interface supernovaeTypeIaZero
     !!{
     Constructors for the \refClass{supernovaeTypeIaZero} supernovae type Ia class.
     !!}
     module procedure zeroConstructorParameters
  end interface supernovaeTypeIaZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{supernovaeTypeIaZero} supernovae type Ia class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(supernovaeTypeIaZero)                :: self
    type(inputParameters     ), intent(inout) :: parameters

   
    self=supernovaeTypeIaZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  subroutine zeroMassInitialRange(self,initialMassFunction_,age,metallicity,massInitialMinimum,massInitialMaximum)
    !!{
    Return the range of initial stellar masses contributing to the Type Ia population.
    !!}
    implicit none
     class           (supernovaeTypeIaZero    ), intent(inout) :: self
     class           (initialMassFunctionClass), intent(inout) :: initialMassFunction_
     double precision                          , intent(in   ) :: age                 , metallicity
     double precision                          , intent(  out) :: massInitialMinimum  , massInitialMaximum
     !$GLC attributes unused :: self, initialMassFunction_, age, metallicity

     ! The range is arbitrary as we have no Type Ias.
     massInitialMinimum=1.0d0
     massInitialMaximum=2.0d0
    return
  end subroutine zeroMassInitialRange
  
  double precision function zeroNumber(self,initialMassFunction_,initialMass,age,metallicity) result(number)
    !!{
    Compute the cumulative number of Type Ia supernovae which is always zero.
    !!}
    implicit none
    class           (supernovaeTypeIaZero    ), intent(inout), target :: self
    class           (initialMassFunctionClass), intent(inout), target :: initialMassFunction_
    double precision                          , intent(in   )         :: age                 , initialMass, &
         &                                                               metallicity
    !$GLC attributes unused :: self, initialMassFunction_, initialMass, age, metallicity
    
    number=0.0d0
    return
  end function zeroNumber

  double precision function zeroYield(self,initialMassFunction_,initialMass,age,metallicity,atomIndex) result(yield)
    !!{
    Compute the cumulative yield from Type Ia supernovae which is always zero.
    !!}
    implicit none
    class           (supernovaeTypeIaZero    ), intent(inout)           :: self
    class           (initialMassFunctionClass), intent(inout)           :: initialMassFunction_
    double precision                          , intent(in   )           :: age                 , initialMass, &
         &                                                                 metallicity
    integer                                   , intent(in   ), optional :: atomIndex
    !$GLC attributes unused :: self, initialMassFunction_, initialMass, age, metallicity, atomIndex
    
    yield=0.0d0    
    return
  end function zeroYield
