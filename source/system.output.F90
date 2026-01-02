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
Contains a module which interfaces with the system output.
!!}

! Specify an explicit dependence on the C interface files.
!: $(BUILDPATH)/isatty.o

module System_Output
  use, intrinsic :: ISO_C_Binding, only : c_int
  implicit none
  private
  public :: stdOutIsATTY

  interface
     function stdOutIsATTY_() bind(c,name='stdOutIsATTY_')
       !!{
       Template for a C function that determines if stdout is a TTY.
       !!}
       import
       integer(c_int) :: stdOutIsATTY_
     end function stdOutIsATTY_
  end interface

contains

  logical function stdOutIsATTY()
    !!{
    Return {\normalfont \ttfamily true} if stdout is a {\normalfont \ttfamily TTY}.
    !!}
    implicit none

    stdOutIsATTY=stdOutIsATTY_() == 1_c_int
    return
  end function stdOutIsATTY
  
end module System_Output
