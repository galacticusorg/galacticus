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

!% Contains a module which defines the base class for all {\normalfont \ttfamily functionClass} classes.

module Function_Classes
  !% Defines the base class for all {\normalfont \ttfamily functionClass} classes.
  implicit none
  private
  public :: functionClass

  type functionClass
     !% The base class for all {\normalfont \ttfamily functionClass} classes.
     logical :: isIndestructible=.false., isDefaultOfClass=.false.
   contains
     !@ <objectMethods>
     !@   <object>functionClass</object>
     !@   <objectMethod>
     !@     <method>isFinalizable</method>
     !@     <type>logical</type>
     !@     <arguments></arguments>
     !@     <description>Return true if this object can be finalized.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>makeIndestructible</method>
     !@     <type>void</type>
     !@     <arguments></arguments>
     !@     <description>Make this object non-finalizable.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isDefault</method>
     !@     <type>logical</type>
     !@     <arguments></arguments>
     !@     <description>Return true if this is the default object of this class.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: isFinalizable      => functionClassIsFinalizable
     procedure :: makeIndestructible => functionClassMakeIndestructible
     procedure :: isDefault          => functionClassIsDefault
  end type functionClass

contains

  logical function functionClassIsFinalizable(self)
    !% Returns true if the object can be finalized.
    class(functionClass), intent(in   ) :: self

    functionClassIsFinalizable=.not.self%isIndestructible
    return
  end function functionClassIsFinalizable
  
  subroutine functionClassMakeIndestructible(self)
    !% Make this object non-finalizable.
    class(functionClass), intent(inout) :: self

    self%isIndestructible=.true.
    return
  end subroutine functionClassMakeIndestructible
  
  logical function functionClassIsDefault(self)
    !% Return true if this is the default object of this class.
    class(functionClass), intent(in   ) :: self

    functionClassIsDefault=self%isDefaultOfClass
    return
  end function functionClassIsDefault
  
end module Function_Classes
