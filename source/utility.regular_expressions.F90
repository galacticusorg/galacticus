!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements regular expressions by wrappring the GNU C Library implementations.

! Specify an explicit dependence on the semaphores.o object file.
!: ./work/build/regular_expressions.o

module Regular_Expressions
  !% Implements regular expressions by wrappring the GNU C Library implementations.
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: regEx

  type :: regEx
     !% A regular expression object.
     type(c_ptr) :: r
   contains
     !@ <objectMethods>
     !@   <object>regEx</object>
     !@   <objectMethod>
     !@     <method>matches</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} string\argin</arguments>
     !@     <description>Return true if a regular expression matches the supplied {\tt string}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Destroy the regex.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::            Regular_Expression_Destructor
     procedure :: destroy => Regular_Expression_Destroy
     procedure :: matches => Regular_Expression_Match
  end type regEx

  interface regEx
     !% Constructor for regular expression object.
     module procedure Regular_Expression_Constructor
  end interface regEx

  interface
     function Regular_Expression_Construct_C(pattern) bind(c,name='Regular_Expression_Construct_C')
       !% Template for a C function that initializes a regular expression.
       import
       type     (c_ptr )        :: Regular_Expression_Construct_C
       character(c_char), value :: pattern
     end function Regular_Expression_Construct_C
  end interface

  interface
     subroutine Regular_Expression_Destruct_C(r) bind(c,name='Regular_Expression_Destruct_C')
       !% Template for a C function that destroys a regular expression.
       import
       type(c_ptr), value :: r
     end subroutine Regular_Expression_Destruct_C
  end interface

  interface
     integer function Regular_Expression_Match_C(r,string) bind(c,name='Regular_Expression_Match_C')
       !% Template for a C function that checks for a match with a regular expression.
       import
       type     (c_ptr ), value :: r
       character(c_char), value :: string
     end function Regular_Expression_Match_C
  end interface

contains

  function Regular_Expression_Constructor(regularExpression)
    !% Constructor for {\tt regEx} objects.
    implicit none
    type     (regEx)                :: Regular_Expression_Constructor
    character(len=*), intent(in   ) :: regularExpression
    
    Regular_Expression_Constructor%r=Regular_Expression_Construct_C(trim(regularExpression)//char(0))
    return
  end function Regular_Expression_Constructor

  subroutine Regular_Expression_Destructor(self)
    !% Destructor for {\tt regEx} objects.
    implicit none
    type(regEx), intent(inout) :: self
    
    if (C_Associated(self%r)) call Regular_Expression_Destruct_C(self%r)
    self%r=C_NULL_PTR
    return
  end subroutine Regular_Expression_Destructor

  subroutine Regular_Expression_Destroy(self)
    !% Destroy a {\tt regEx} object.
    implicit none
    class(regEx), intent(inout) :: self
    
    select type (self)
    type is (regEx)
       call Regular_Expression_Destructor(self)
    end select
    return
  end subroutine Regular_Expression_Destroy

  logical function Regular_Expression_Match(self,string)
    !% Returns true if a {\tt regEx} object matches the supplied {\tt string}.
    implicit none
    class    (regEx)                :: self
    character(len=*), intent(in   ) :: string
    integer                         :: status
 
    status=Regular_Expression_Match_C(self%r,trim(string)//char(0))
    Regular_Expression_Match=(status==0)
    return
  end function Regular_Expression_Match

end module Regular_Expressions
