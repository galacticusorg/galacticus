!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Implements a radiation field class which sums over other radiation fields.

  type, public :: radiationFieldList
     class(radiationFieldClass), pointer :: radiationField_
     type (radiationFieldList ), pointer :: next            => null()
  end type radiationFieldList

  !# <radiationField name="radiationFieldSummation">
  !#  <description>A summation radiation field class.</description>
  !# </radiationField>
  type, extends(radiationFieldClass) :: radiationFieldSummation
     !% A summation radiation field class.
     private
     type(radiationFieldList), pointer :: radiationFields
   contains
     !@ <objectMethods>
     !@  <object>radiationFieldSummation</object>
     !@  <objectMethod>
     !@   <method>list</method>
     !@   <type>\void</type>
     !@   <arguments>\textcolor{red}{\textless *type(radiationFieldList)\textgreater}</arguments>
     !@   <description>Return a list of all sub-components.</description>
     !@  </objectMethod>
     !@ </objectMethods>
     final     ::             summationDestructor
     procedure :: flux     => summationFlux
     procedure :: deepCopy => summationDeepCopy
     procedure :: list     => summationList
  end type radiationFieldSummation

  interface radiationFieldSummation
     !% Constructors for the ``summation'' radiation field class.
     module procedure summationConstructorParameters
     module procedure summationConstructorInternal
  end interface radiationFieldSummation

contains

  function summationConstructorParameters(parameters) result (self)
    !% Constructor for the ``summation'' radiation field class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type   (radiationFieldSummation)                :: self
    type   (inputParameters        ), intent(inout) :: parameters
    type   (radiationFieldList     ), pointer       :: radiationField_
    integer                                         :: i

    self     %radiationFields => null()
    radiationField_           => null()
    do i=1,parameters%copiesCount('radiationFieldMethod',zeroIfNotPresent=.true.)
       if (associated(radiationField_)) then
          allocate(radiationField_%next)
          radiationField_ => radiationField_%next
       else
          allocate(self%radiationFields)
          radiationField_ => self%radiationFields
       end if
       !# <objectBuilder class="radiationField" name="radiationField_%radiationField_" source="parameters" copy="i" />
    end do
    return
  end function summationConstructorParameters

  function summationConstructorInternal(radiationFields) result (self)
    !% Internal constructor for the summation radiation field class.
    implicit none
    type(radiationFieldSummation)                        :: self
    type(radiationFieldList     ), target, intent(in   ) :: radiationFields
    type(radiationFieldList     ), pointer               :: radiationField_

    self           %radiationFields => radiationFields
    radiationField_                 => radiationFields
    do while (associated(radiationField_))
       !# <referenceCountIncrement owner="radiationField_" object="radiationField_"/>
       radiationField_ => radiationField_%next
    end do
    return
  end function summationConstructorInternal

  subroutine summationDestructor(self)
    !% Destructor for the summation radiation field class.
    implicit none
    type(radiationFieldSummation), intent(inout) :: self
    type(radiationFieldList     ), pointer       :: radiationField_, radiationFieldNext

    if (associated(self%radiationFields)) then
       radiationField_ => self%radiationFields
       do while (associated(radiationField_))
          radiationFieldNext => radiationField_%next
          !# <objectDestructor name="radiationField_%radiationField_"/>
          deallocate(radiationField_)
          radiationField_ => radiationFieldNext
       end do
    end if
    return
  end subroutine summationDestructor

  double precision function summationFlux(self,wavelength,node)
    !% Implement a summation radiation field.
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
  
  subroutine summationDeepCopy(self,destination)
    !% Perform a deep copy for the {\normalfont \ttfamily summation} radiation field class.
    use Galacticus_Error
    implicit none
    class(radiationFieldSummation), intent(inout) :: self
    class(radiationFieldClass    ), intent(  out) :: destination
    type (radiationFieldList     ), pointer       :: radiationField_   , radiationFieldDestination_, &
         &                                           radiationFieldNew_

    call self%radiationFieldClass%deepCopy(destination)
    select type (destination)
    type is (radiationFieldSummation)
       destination%radiationFields => null          ()
       radiationFieldDestination_  => null          ()
       radiationField_             => self%radiationFields
       do while (associated(radiationField_))
          allocate(radiationFieldNew_)
          if (associated(radiationFieldDestination_)) then
             radiationFieldDestination_%next       => radiationFieldNew_
             radiationFieldDestination_            => radiationFieldNew_             
          else
             destination          %radiationFields => radiationFieldNew_
             radiationFieldDestination_            => radiationFieldNew_
          end if
          allocate(radiationFieldNew_%radiationField_,mold=radiationField_%radiationField_)
          !# <deepCopy source="radiationField_%radiationField_" destination="radiationFieldNew_%radiationField_"/>
          radiationField_ => radiationField_%next
       end do       
    class default
       call Galacticus_Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine summationDeepCopy

  function summationList(self)
    !% Return a list of all components for the {\normalfont \ttfamily summation} radiation field class.
    implicit none
    class(radiationFieldSummation), intent(inout) :: self
    type (radiationFieldList     ), pointer       :: summationList

    summationList => self%radiationFields
    return
  end function summationList
