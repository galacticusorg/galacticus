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

!% Contains a module which interfaces with low-level aspects of the GSL library.

! Specify an explicit dependence on the interface.GSL.C.o object file.
!: $(BUILDPATH)/interface.GSL.C.o


module Interface_GSL
  !% Interfaces with low-level aspects of the GSL library.
  use, intrinsic :: ISO_C_Binding, only : c_funptr, c_ptr
  private
  public :: gslFunction, gslFunctionDestroy, gslFunctionTemplate

  abstract interface
     !% Interface for {\normalfont \ttfamily gsl\_function} type. We ignore the {\normalfont \ttfamily params} argument as this is not used by \glc.
     double precision function gslFunctionTemplate(x)
       double precision, intent(in   ), value :: x
     end function gslFunctionTemplate
  end interface

   interface
      !% Interfaces to C functions.
      function gslFunctionConstructor(f) bind(c,name="gslFunctionConstructor")
        !% Interface to a C function which establishes a {\normalfont \ttfamily gsl\_function} type.
        import c_ptr, c_funptr
        type(c_ptr   )        :: gslFunctionConstructor
        type(c_funptr), value :: f
      end function gslFunctionConstructor

      subroutine gslFunctionDestructor(f) bind(c,name="gslFunctionDestructor")
        !% Interface to a C function which destroys a {\normalfont \ttfamily gsl\_function} type.
        import c_funptr
        type(c_funptr), value :: f
      end subroutine gslFunctionDestructor
   end interface

contains

  function gslFunction(f)
    !% Return a {\normalfont \ttfamily c\_ptr} object for the given function {\normalfont \ttfamily f}.
    use, intrinsic :: ISO_C_Binding, only : c_funloc
    implicit none
    type     (c_ptr              ) :: gslFunction
    procedure(gslFunctionTemplate) :: f

    gslFunction=gslFunctionConstructor(c_funloc(f))
    return
  end function gslFunction

  subroutine gslFunctionDestroy(f)
    !% Destroy a {\normalfont \ttfamily c\_ptr} to a {\normalfont \ttfamily gsl\_function} object.
    implicit none
    type(c_ptr), intent(in   ) :: f

    call gslFunctionDestructor(f)
    return
  end subroutine gslFunctionDestroy

end module Interface_GSL
