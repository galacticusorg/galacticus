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
Contains a module which provides basic arithmetic functions.
!!}

module Math_Arithmetic
  !!{
  Provides basic arithmetic functions, often handling floating point issues.
  !!}
  implicit none
  private
  public :: divideSafe

contains


  pure double precision function divideSafe(x,y)
    !!{
    Compute $x/y$ but avoiding floating point overflow. Where overflow would occur the result is limited to the largest
    representable value.    
    !!}
    implicit none
    double precision, intent(in   ) :: x, y

    if (exponent(x)-exponent(y) < maxExponent(1.0d0)) then
       divideSafe=+x             &
            &     /y
    else
       divideSafe=+huge(1.0d0  ) &
            &     *sign(1.0d0,x) &
            &     *sign(1.0d0,y)
    end if
    return
  end function divideSafe

end module Math_Arithmetic
