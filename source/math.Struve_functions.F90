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
Contains a module which implements Struve functions.
!!}

module Struve_Functions
  !!{
  Implements Struve functions.
  !!}
  implicit none
  private
  public :: Struve_Function_L1

contains

  double precision function Struve_Function_L1(x)
    !!{
    Evaluate and return the Struve $L_1$ function.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    call stvl1(x,Struve_Function_L1)
    return
  end function Struve_Function_L1

  ! Following are the original functions taken with permission from copyrighted work of:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special Functions,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  ! The source code was extracted from:
  !
  !    https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90

  subroutine stvl1 ( x, sl1 )

    !*****************************************************************************80
    !
    !! STVL1 computes the modified Struve function L1(x).
    !
    !  Licensing:
    !
    !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
    !    they give permission to incorporate this routine into a user program
    !    provided that the copyright is acknowledged.
    !
    !  Modified:
    !
    !    05 July 2012
    !
    !  Author:
    !
    !    Shanjie Zhang, Jianming Jin
    !
    !  Reference:
    !
    !    Shanjie Zhang, Jianming Jin,
    !    Computation of Special Functions,
    !    Wiley, 1996,
    !    ISBN: 0-471-11963-6,
    !    LC: QA351.C45.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument.
    !
    !    Output, real ( kind = 8 ) SL1, the function value.
    !
    implicit none

    real ( kind = 8 ) a1
    real ( kind = 8 ) bi1
    integer ( kind = 4 ) k
    integer ( kind = 4 ) km
    real ( kind = 8 ) pi
    real ( kind = 8 ) r
    real ( kind = 8 ) s
    real ( kind = 8 ) sl1
    real ( kind = 8 ) x

    pi = 3.141592653589793D+00
    r = 1.0D+00
    if ( x <= 20.0D+00 ) then
       s = 0.0D+00
       do k = 1, 60
          r = r * x * x / ( 4.0D+00 * k * k - 1.0D+00 )
          s = s + r
          if ( abs ( r / s ) < 1.0D-12 ) then
             exit
          end if
       end do

       sl1 = 2.0D+00 / pi * s

    else

       s = 1.0D+00
       km = int ( 0.50D+00 * x )
       km = min ( km, 25 )

       do k = 1, km
          r = r * ( 2.0D+00 * k + 3.0D+00 ) &
               * ( 2.0D+00 * k + 1.0D+00 ) / ( x * x )
          s = s + r
          if ( abs ( r / s ) < 1.0D-12 ) then
             exit
          end if
       end do

       sl1 = 2.0D+00 / pi * ( -1.0D+00 + 1.0D+00 &
            / ( x * x ) + 3.0D+00 * s / x**4 )
       a1 = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
       r = 1.0D+00
       bi1 = 1.0D+00
       do k = 1, 16
          r = -0.125D+00 * r &
               * ( 4.0D+00 - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
          bi1 = bi1 + r
          if ( abs ( r / bi1 ) < 1.0D-12 ) then
             exit
          end if
       end do

       sl1 = sl1 + a1 * bi1

    end if

    return
  end subroutine stvl1

end module Struve_Functions
