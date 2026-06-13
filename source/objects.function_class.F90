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

!!{RST
Contains a module which defines the base class for all ``functionClass`` classes.
!!}

module Function_Classes
  !!{RST
  Defines the base class for all ``functionClass`` classes.
  !!}
  use :: ISO_Varying_String, only : assignment(=), varying_string
  implicit none
  private
  public :: functionClass

  type, abstract :: functionClass
     !!{RST
     The base class for all ``functionClass`` classes.
     !!}
     logical :: isDefaultOfClass=.false., reportOn_=.false.
     integer :: referenceCount  =0
   contains
     !![
     <methods docformat="rst">
       <method method="referenceCountReset"     description="Reset the reference count to this object to 0."                                  />
       <method method="referenceCountIncrement" description="Increment the reference count to this object."                                   />
       <method method="referenceCountDecrement" description="Decrement the reference count to this object and return the new reference count."/>
       <method method="isDefault"               description="Return true if this is the default object of this class."                        />
       <method method="reportOn"                description="Indicate that reference count changes to this object should be reported on."     />
     </methods>
     !!]
     procedure :: isDefault               => functionClassIsDefault
     procedure :: reportOn                => functionClassReportOn
     procedure :: referenceCountReset     => functionClassReferenceCountReset
     procedure :: referenceCountIncrement => functionClassReferenceCountIncrement
     procedure :: referenceCountDecrement => functionClassReferenceCountDecrement
  end type functionClass

contains

  subroutine functionClassReportOn(self)
    !!{RST
    Indicate that reference count changes to this object should be reported on.
    !!}
    use :: Display           , only : displayMessage
    use :: String_Handling   , only : operator(//)
    use :: ISO_Varying_String, only : operator(//)  , var_str
    implicit none
    class(functionClass), intent(inout) :: self

    self%reportOn_=.true.
    call displayMessage(var_str('reporting on functionClass object [loc:')//loc(self)//' ] references - count = '//self%referenceCount)
    call backtrace()
    return
  end subroutine functionClassReportOn

  logical function functionClassIsDefault(self)
    !!{RST
    Return true if this is the default object of this class.
    !!}
    implicit none
    class(functionClass), intent(in   ) :: self

    functionClassIsDefault=self%isDefaultOfClass
    return
  end function functionClassIsDefault

  subroutine functionClassReferenceCountReset(self)
    !!{RST
    Reset the reference count to this object to 0.
    !!}
    use :: Display           , only : displayMessage
    use :: String_Handling   , only : operator(//)
    use :: ISO_Varying_String, only : operator(//)  , var_str
    implicit none
    class(functionClass), intent(inout) :: self

    self%referenceCount=1
    if (self%reportOn_) then
       call displayMessage(var_str('reset functionClass object [loc:')//loc(self)//' ] references - count = '//self%referenceCount)
       call backtrace()
    end if
    return
  end subroutine functionClassReferenceCountReset

  subroutine functionClassReferenceCountIncrement(self)
    !!{RST
    Increment the reference count to this object.
    !!}
    use :: Display           , only : displayMessage
    use :: String_Handling   , only : operator(//)
    use :: ISO_Varying_String, only : operator(//)  , var_str
    implicit none
    class(functionClass), intent(inout) :: self

    self%referenceCount=self%referenceCount+1
    if (self%reportOn_) then
       call displayMessage(var_str('increment functionClass object [loc:')//loc(self)//' ] references - count = '//self%referenceCount)
       call backtrace()
    end if
    return
  end subroutine functionClassReferenceCountIncrement

  integer function functionClassReferenceCountDecrement(self)
    !!{RST
    Decrement the reference count to this object and return the new count.
    !!}
    use :: Display           , only : displayMessage
    use :: String_Handling   , only : operator(//)
    use :: ISO_Varying_String, only : operator(//)  , var_str
    implicit none
    class(functionClass), intent(inout) :: self

    self%referenceCount=self%referenceCount-1
    functionClassReferenceCountDecrement=self%referenceCount
    if (self%reportOn_) then
       call displayMessage(var_str('decrement functionClass object [loc:')//loc(self)//' ] references - count = '//self%referenceCount)
       call backtrace()
    end if
    return
  end function functionClassReferenceCountDecrement

end module Function_Classes
