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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module which provides utilities for tasks which build external tools.
!!}

module Tasks_Build_Tool_Utilities
  !!{RST
  Provides utilities for tasks which build external tools.
  !!}
  implicit none
  private
  public :: Task_Build_Tool

  abstract interface
     subroutine buildToolInitializer(path,version,static)
       !!{RST
       Interface for procedures which build an external tool, returning the path to, and version of, the tool built.
       !!}
       use :: ISO_Varying_String, only : varying_string
       type   (varying_string), intent(  out)           :: path  , version
       logical                , intent(in   ), optional :: static
     end subroutine buildToolInitializer
  end interface

contains

  subroutine Task_Build_Tool(name,initializer,status)
    !!{RST
    Build the named external tool using the given procedure, reporting the version built and where it was installed.
    !!}
    use :: Display           , only : displayIndent     , displayMessage, displayUnindent
    use :: Error             , only : errorStatusSuccess
    use :: ISO_Varying_String, only : varying_string    , operator(//)
    implicit none
    character(len=*               ), intent(in   )           :: name
    procedure(buildToolInitializer)                          :: initializer
    integer                        , intent(  out), optional :: status
    type     (varying_string      )                          :: path       , version
#include "os.inc"

    call displayIndent  ('Begin task: '//name//' tool build')
    call initializer(                &
         &                  path   , &
         &                  version, &
#ifdef __APPLE__
         &           static=.false.  &
#else
         &           static=.true.   &
#endif
         &          )
    call displayMessage (name//' version '//version//' successfully built in: '//path)
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: '//name//' tool build')
    return
  end subroutine Task_Build_Tool

end module Tasks_Build_Tool_Utilities
