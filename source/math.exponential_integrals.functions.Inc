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

! Mathematical special function code taken from John Burkardt's collection:
!  http://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html
!  http://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.f90

  subroutine e1z ( z, ce1 )
  !*****************************************************************************80
  !
  !! E1Z computes the complex exponential integral E1(z).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
  !    they give permission to incorporate this routine into a user program
  !    provided that the copyright is acknowledged.
  !
  !  Modified:
  !
  !    16 July 2012
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
  !    Input, complex ( kind = 8 ) Z, the argument.
  !
  !    Output, complex ( kind = 8 ) CE1, the function value.
  !
  implicit none

  real ( kind = 8 ) a0
  complex ( kind = 8 ) ce1
  complex ( kind = 8 ) cr
  complex ( kind = 8 ) ct
  complex ( kind = 8 ) ct0
  real ( kind = 8 ) el
  integer ( kind = 4 ) k
  real ( kind = 8 ) pi
  real ( kind = 8 ) x
  complex ( kind = 8 ) z

  pi = 3.141592653589793D+00
  el = 0.5772156649015328D+00
  x = real ( z, kind = 8 )
  a0 = abs ( z )

  if ( a0 == 0.0D+00 ) then
     ce1 = cmplx ( 1.0D+300, 0.0D+00, kind = 8 )
  else if ( a0 <= 10.0D+00 .or. &
       ( x < 0.0D+00 .and. a0 < 20.0D+00 ) ) then
     ce1 = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     cr = cmplx ( 1.0D+00, 0.0D+00, kind = 8 )
     do k = 1, 150
        cr = - cr * k * z / ( k + 1.0D+00 )**2
        ce1 = ce1 + cr
        if ( abs ( cr ) <= abs ( ce1 ) * 1.0D-15 ) then
           exit
        end if
     end do

     ce1 = - el - log ( z ) + z * ce1

  else

     ct0 = cmplx ( 0.0D+00, 0.0D+00, kind = 8 )
     do k = 120, 1, -1
        ct0 = k / ( 1.0D+00 + k / ( z + ct0 ) )
     end do
     ct = 1.0D+00 / ( z + ct0 )

     ce1 = exp ( - z ) * ct
     if ( x <= 0.0D+00 .and. imag ( z ) == 0.0D+00 ) then
        ce1 = ce1 - pi * cmplx ( 0.0D+00, 1.0D+00, kind = 8 )
     end if

  end if

  return
end subroutine e1z
