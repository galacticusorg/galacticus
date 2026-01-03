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
  Implements a \gls{cgm} heating class that sums over other heating rate classes.
  !!}

  !![
  <circumgalacticMediumHeating name="circumgalacticMediumHeatingSummation">
   <description>
    A \gls{cgm} heating class that sums over other heating rate classes.
   </description>
   <linkedList type="heaterList" variable="heaters" next="next" object="circumgalacticMediumHeating" objectType="circumgalacticMediumHeatingClass"/>
  </circumgalacticMediumHeating>
  !!]

  type, public :: heaterList
     class(circumgalacticMediumHeatingClass), pointer :: circumgalacticMediumHeating => null()
     type (heaterList                      ), pointer :: next                        => null()
  end type heaterList

  type, extends(circumgalacticMediumHeatingClass) :: circumgalacticMediumHeatingSummation
     !!{
     A \gls{cgm} heating class that sums over other heating rate classes.
     !!}
     private
     type(heaterList), pointer :: heaters => null()
   contains
     final     ::                summationDestructor
     procedure :: heatingRate => summationHeatingRate
  end type circumgalacticMediumHeatingSummation
  
  interface circumgalacticMediumHeatingSummation
     !!{
     Constructors for the \refClass{circumgalacticMediumHeatingSummation} class.
     !!}
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface circumgalacticMediumHeatingSummation

contains

  function summationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{circumgalacticMediumHeatingSummation} class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (circumgalacticMediumHeatingSummation)                :: self
    type   (inputParameters                     ), intent(inout) :: parameters
    type   (heaterList                          ), pointer       :: heater
    integer                                                      :: i

    heater => null()
    do i=1,parameters%copiesCount('circumgalacticMediumHeating',zeroIfNotPresent=.true.)
       if (associated(heater)) then
          allocate(heater%next)
          heater => heater%next
       else
          allocate(self%heaters)
          heater => self%heaters
       end if
       !![
       <objectBuilder class="circumgalacticMediumHeating" name="heater%circumgalacticMediumHeating" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="circumgalacticMediumHeating"/>
    !!]
    return
  end function summationConstructorParameters

  function summationConstructorInternal(heaters) result(self)
    !!{
    Internal constructor for the \refClass{circumgalacticMediumHeatingSummation} cooling function class.
    !!}
    implicit none
    type(circumgalacticMediumHeatingSummation)                        :: self
    type(heaterList                          ), target, intent(in   ) :: heaters
    type(heaterList                          ), pointer               :: heater_

    self   %heaters => heaters
    heater_         => heaters
    do while (associated(heater_))
       !![
       <referenceCountIncrement owner="heater_" object="circumgalacticMediumHeating"/>
       !!]
       heater_ => heater_%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !!{
    Destructor for the \refClass{circumgalacticMediumHeatingSummation} cooling function class.
    !!}
    implicit none
    type(circumgalacticMediumHeatingSummation), intent(inout) :: self
    type(heaterList                          ), pointer       :: heater_, heaterNext

    if (associated(self%heaters)) then
       heater_ => self%heaters
       do while (associated(heater_))
          heaterNext => heater_%next
          !![
          <objectDestructor name="heater_%circumgalacticMediumHeating"/>
          !!]
          deallocate(heater_)
          heater_ => heaterNext
       end do
    end if
    return
  end subroutine summationDestructor

  double precision function summationHeatingRate(self,node) result(rateHeating)
    !!{
    Compute the heating rate of the \gls{cgm}, assumed to be always summation.
    !!}
    implicit none
    class(circumgalacticMediumHeatingSummation), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    type(heaterList                           ), pointer       :: heater

    rateHeating =  0.0d0
    heater      => self%heaters
    do while (associated(heater))
       rateHeating=+rateHeating                                          &
            &      +heater%circumgalacticMediumHeating%heatingRate(node)
       heater => heater%next
    end do
    return
  end function summationHeatingRate
