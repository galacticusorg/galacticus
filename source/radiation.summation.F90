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
Implements a radiation field class which sums over other radiation fields.
!!}

  type, public :: radiationFieldList
     class(radiationFieldClass), pointer :: radiationField_ => null()
     type (radiationFieldList ), pointer :: next            => null()
  end type radiationFieldList

  !![
  <radiationField name="radiationFieldSummation">
   <description>A summation radiation field class.</description>
   <linkedList type="radiationFieldList" variable="radiationFields" next="next" object="radiationField_" objectType="radiationFieldClass"/>
  </radiationField>
  !!]
  type, extends(radiationFieldClass) :: radiationFieldSummation
     !!{
     A summation radiation field class.
     !!}
     private
     type(radiationFieldList), pointer :: radiationFields => null()
   contains
     !![
     <methods>
       <method description="Return a list of all sub-components." method="list" />
     </methods>
     !!]
     final     ::                      summationDestructor
     procedure :: flux              => summationFlux
     procedure :: list              => summationList
     procedure :: time              => summationTime
     procedure :: timeSet           => summationTimeSet
     procedure :: timeDependentOnly => summationTimeDependentOnly
  end type radiationFieldSummation

  interface radiationFieldSummation
     !!{
     Constructors for the {\normalfont \ttfamily summation} radiation field class.
     !!}
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface radiationFieldSummation

contains

  function summationConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily summation} radiation field class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (radiationFieldSummation)                :: self
    type   (inputParameters        ), intent(inout) :: parameters
    type   (radiationFieldList     ), pointer       :: radiationField_
    integer                                         :: i

    self     %radiationFields => null()
    radiationField_           => null()
    do i=1,parameters%copiesCount('radiationField',zeroIfNotPresent=.true.)
       if (associated(radiationField_)) then
          allocate(radiationField_%next)
          radiationField_ => radiationField_%next
       else
          allocate(self%radiationFields)
          radiationField_ => self%radiationFields
       end if
       !![
       <objectBuilder class="radiationField" name="radiationField_%radiationField_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="radiationField"/>
    !!]
    return
  end function summationConstructorParameters

  function summationConstructorInternal(radiationFields) result (self)
    !!{
    Internal constructor for the summation radiation field class.
    !!}
    implicit none
    type(radiationFieldSummation)                        :: self
    type(radiationFieldList     ), target, intent(in   ) :: radiationFields
    type(radiationFieldList     ), pointer               :: radiationField_

    self           %radiationFields => radiationFields
    radiationField_                 => radiationFields
    do while (associated(radiationField_))
       !![
       <referenceCountIncrement owner="radiationField_" object="radiationField_"/>
       !!]
       radiationField_ => radiationField_%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !!{
    Destructor for the summation radiation field class.
    !!}
    implicit none
    type(radiationFieldSummation), intent(inout) :: self
    type(radiationFieldList     ), pointer       :: radiationField_, radiationFieldNext

    if (associated(self%radiationFields)) then
       radiationField_ => self%radiationFields
       do while (associated(radiationField_))
          radiationFieldNext => radiationField_%next
          !![
          <objectDestructor name="radiationField_%radiationField_"/>
          !!]
          deallocate(radiationField_)
          radiationField_ => radiationFieldNext
       end do
    end if
    return
  end subroutine summationDestructor

  double precision function summationFlux(self,wavelength,node)
    !!{
    Implement a summation radiation field.
    !!}
    implicit none
    class           (radiationFieldSummation), intent(inout) :: self
    double precision                         , intent(in   ) :: wavelength
    type            (treeNode               ), intent(inout) :: node
    type            (radiationFieldList     ), pointer       :: radiationField_

    summationFlux   =  0.0d0
    radiationField_ => self%radiationFields
    do while (associated(radiationField_))
       summationFlux   =  +summationFlux                                         &
            &             +radiationField_%radiationField_%flux(wavelength,node)
       radiationField_ =>  radiationField_%next
    end do
    return
  end function summationFlux

  double precision function summationTime(self)
    !!{
    Return the time for which this radiation field is set.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (radiationFieldSummation), intent(inout) :: self
    type            (radiationFieldList     ), pointer       :: radiationField_
    double precision                                         :: time

    summationTime   =  -huge(0.0d0)  
    radiationField_ => self%radiationFields
    do while (associated(radiationField_))
       time=radiationField_%radiationField_%time()
       if (time /= summationTime .and. summationTime /= -huge(0.0d0)) &
            & call Error_Report('inconsistent times in radiation field components'//{introspection:location})
       summationTime   =  time
       radiationField_ => radiationField_%next
    end do
    return
  end function summationTime

  subroutine summationTimeSet(self,time)
    !!{
    Set the time (and temperature) of the cosmic microwave background radiation field.
    !!}
    implicit none
    class           (radiationFieldSummation), intent(inout) :: self
    double precision                         , intent(in   ) :: time
    type            (radiationFieldList     ), pointer       :: radiationField_

    radiationField_ => self%radiationFields
    do while (associated(radiationField_))
       call radiationField_%radiationField_%timeSet(time)
       radiationField_ =>  radiationField_%next
    end do
    return
  end subroutine summationTimeSet

  logical function summationTimeDependentOnly(self)
    !!{
    Return true as this radiation field depends on time only.
    !!}
    implicit none
    class(radiationFieldSummation), intent(inout) :: self
    type (radiationFieldList     ), pointer       :: radiationField_

    summationTimeDependentOnly =  .true.
    radiationField_            => self%radiationFields
    do while (associated(radiationField_))
       if (.not.radiationField_%radiationField_%timeDependentOnly()) then
          summationTimeDependentOnly=.false.
          return
       end if
       radiationField_ =>  radiationField_%next
    end do
    return
  end function summationTimeDependentOnly

  function summationList(self)
    !!{
    Return a list of all components for the {\normalfont \ttfamily summation} radiation field class.
    !!}
    implicit none
    class(radiationFieldSummation), intent(inout) :: self
    type (radiationFieldList     ), pointer       :: summationList

    summationList => self%radiationFields
    return
  end function summationList
