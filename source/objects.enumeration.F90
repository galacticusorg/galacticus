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
Contains a module which defines the base class for all {\normalfont \ttfamily enumeration} classes.
!!}

module Enumerations
  !!{
  Defines the base class for all {\normalfont \ttfamily enumeration} classes.
  !!}
  implicit none
  private
  public :: enumerationType

  type, abstract :: enumerationType
     !!{
     The base class for all {\normalfont \ttfamily enumeration} classes.
     !!}
     integer :: ID
   contains
     !![
     <methods>
       <method description="Test if two enumeration entries are not equal."               method="operator(/=)"  />
       <method description="Test if one enumeration is less than another."                method="operator(&lt;)" />
       <method description="Test if one enumeration is less or equal tothan another."     method="operator(&lt;=)"/>
       <method description="Test if one enumeration is greater than another."             method="operator(&gt;)" />
       <method description="Test if one enumeration is greater than or equal to another." method="operator(&gt;=)"/>
       <method description="Subtract an integer from an enumeration member."              method="subtract"       />
     </methods>
     !!]
     procedure ::                 enumerationIsNotEqual
     procedure ::                 enumerationLessThan
     procedure ::                 enumerationLessThanOrEqual
     procedure ::                 enumerationGreaterThan
     procedure ::                 enumerationGreaterThanOrEqual
     procedure ::                 enumerationSubtractionInteger
     generic   :: operator(/=) => enumerationIsNotEqual
     generic   :: operator(< ) => enumerationLessThan
     generic   :: operator(<=) => enumerationLessThanOrEqual
     generic   :: operator(> ) => enumerationGreaterThan
     generic   :: operator(>=) => enumerationGreaterThanOrEqual
     generic   :: subtract     => enumerationSubtractionInteger
  end type enumerationType

contains

  elemental logical function enumerationIsNotEqual(enumerationA,enumerationB)
    !!{
    Return true if {\normalfont \ttfamily enumerationA} is not equal to {\normalfont \ttfamily enumerationB}.
    !!}
    implicit none
    class(enumerationType), intent(in   ) :: enumerationA, enumerationB

    if (.not.same_type_as(enumerationA,enumerationB)) then
       enumerationIsNotEqual=.true.
    else
       enumerationIsNotEqual=enumerationA%ID /= enumerationB%ID
    end if
    return
  end function enumerationIsNotEqual

  elemental logical function enumerationLessThan(enumerationA,enumerationB)
    !!{
    Return true if {\normalfont \ttfamily enumerationA} is less than {\normalfont \ttfamily enumerationB}.
    !!}
    implicit none
    class(enumerationType), intent(in   ) :: enumerationA, enumerationB

    if (.not.same_type_as(enumerationA,enumerationB)) then
       enumerationLessThan=.false.
    else
       enumerationLessThan=enumerationA%ID < enumerationB%ID
    end if
    return
  end function enumerationLessThan

  elemental logical function enumerationLessThanOrEqual(enumerationA,enumerationB)
    !!{
    Return true if {\normalfont \ttfamily enumerationA} is less than or equal to {\normalfont \ttfamily enumerationB}.
    !!}
    implicit none
    class(enumerationType), intent(in   ) :: enumerationA, enumerationB

    if (.not.same_type_as(enumerationA,enumerationB)) then
       enumerationLessThanOrEqual=.false.
    else
       enumerationLessThanOrEqual=enumerationA%ID <= enumerationB%ID
    end if
    return
  end function enumerationLessThanOrEqual

  elemental logical function enumerationGreaterThan(enumerationA,enumerationB)
    !!{
    Return true if {\normalfont \ttfamily enumerationA} is greater than {\normalfont \ttfamily enumerationB}.
    !!}
    implicit none
    class(enumerationType), intent(in   ) :: enumerationA, enumerationB

    if (.not.same_type_as(enumerationA,enumerationB)) then
       enumerationGreaterThan=.false.
    else
       enumerationGreaterThan=enumerationA%ID > enumerationB%ID
    end if
    return
  end function enumerationGreaterThan

  elemental logical function enumerationGreaterThanOrEqual(enumerationA,enumerationB)
    !!{
    Return true if {\normalfont \ttfamily enumerationA} is greater than or equal to {\normalfont \ttfamily enumerationB}.
    !!}
    implicit none
    class(enumerationType), intent(in   ) :: enumerationA, enumerationB
    
    if (.not.same_type_as(enumerationA,enumerationB)) then
       enumerationGreaterThanOrEqual=.false.
    else
       enumerationGreaterThanOrEqual=enumerationA%ID >= enumerationB%ID
    end if
    return
  end function enumerationGreaterThanOrEqual

  elemental subroutine enumerationSubtractionInteger(enumerationA,enumerationBID)
    !!{
    Subtract an integer ID from {\normalfont \ttfamily enumerationA}.
    !!}
    implicit none
    class  (enumerationType), intent(inout) :: enumerationA
    integer                 , intent(in   ) :: enumerationBID
    
    enumerationA%ID=+enumerationA%ID &
         &          -enumerationBID
    return
  end subroutine enumerationSubtractionInteger
  
end module Enumerations
