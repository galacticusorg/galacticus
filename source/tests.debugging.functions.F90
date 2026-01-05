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
Contains functions used in debugging tests.
!!}

module Test_Debugging_Functions
  !!{
  Contains functions used in debugging tests.
  !!}
  private
  public :: dummyCaller
  
contains

  subroutine dummyCaller(caller)
    !!{
    Callable function to generate a backtrace.
    !!}
    use :: ISO_Varying_String, only : varying_string
    use :: Debugging         , only : getCaller
    implicit none
    type(varying_string), intent(  out) :: caller

    caller=getCaller()
    return
  end subroutine dummyCaller
  
end module Test_Debugging_Functions
