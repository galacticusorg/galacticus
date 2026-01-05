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
  Implementation of a cooling rate class which modifies another cooling rate by cutting off cooling in satellites.
  !!}


  !![
  <coolingRate name="coolingRateMultiplier">
   <description>A cooling rate class which modifies another cooling rate by multiplying the rate by a fixed value.</description>
  </coolingRate>
  !!]
  type, extends(coolingRateClass) :: coolingRateMultiplier
     !!{
     Implementation of cooling rate class which modifies another cooling rate by multiplying the rate by a fixed value.
     !!}
     private
     class           (coolingRateClass), pointer :: coolingRate_ => null()
     double precision                            :: multiplier
   contains
     final     ::         multiplierDestructor
     procedure :: rate => multiplierRate
  end type coolingRateMultiplier

  interface coolingRateMultiplier
     !!{
     Constructors for the cut off cooling rate class.
     !!}
     module procedure multiplierConstructorParameters
     module procedure multiplierConstructorInternal
  end interface coolingRateMultiplier

contains

  function multiplierConstructorParameters(parameters) result(self)
    !!{
    Constructor for the cut off cooling rate class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (coolingRateMultiplier)                :: self
    type            (inputParameters      ), intent(inout) :: parameters
    class           (coolingRateClass     ), pointer       :: coolingRate_
    double precision                                       :: multiplier

    !![
    <inputParameter>
     <name>multiplier</name>
     <source>parameters</source>
     <defaultValue>0.0d0</defaultValue>
     <description>The value by which cooling rates should be multiplied.</description>
    </inputParameter>
    <objectBuilder class="coolingRate" name="coolingRate_" source="parameters"/>
    !!]
    self=coolingRateMultiplier(multiplier,coolingRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="coolingRate_"/>
    !!]
    return
  end function multiplierConstructorParameters

  function multiplierConstructorInternal(multiplier,coolingRate_) result(self)
    !!{
    Internal constructor for the cut off cooling rate class.
    !!}
    type            (coolingRateMultiplier)                        :: self
    class           (coolingRateClass     ), intent(in   ), target :: coolingRate_
    double precision                       , intent(in   )         :: multiplier

    !![
    <constructorAssign variables="multiplier, *coolingRate_"/>
    !!]
    return
  end function multiplierConstructorInternal

  subroutine multiplierDestructor(self)
    !!{
    Destructor for the cut off cooling rate class.
    !!}
    implicit none
    type(coolingRateMultiplier), intent(inout) :: self

    !![
    <objectDestructor name="self%coolingRate_"/>
    !!]
    return
  end subroutine multiplierDestructor

  double precision function multiplierRate(self,node)
    !!{
    Returns the cooling rate (in $M_\odot$ Gyr$^{-1}$) in the hot atmosphere for a model in which this rate is multiplied by
    some fixed value.
    !!}
    implicit none
    class(coolingRateMultiplier), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node

    multiplierRate=+self%multiplier              &
         &         *self%coolingRate_%rate(node)
    return
  end function multiplierRate

