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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which implements the Riemann zeta function.
!!}

! Add dependency on GSL library.
!; gsl

module Zeta_Functions
  !!{RST
  Implements the Riemann zeta function.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double
  implicit none
  private
  public :: Zeta_Function

  interface
     function gsl_sf_zeta(s) bind(c,name='gsl_sf_zeta')
       !!{RST
       Template for the GSL Riemann zeta function.
       !!}
       import
       real(c_double)        :: gsl_sf_zeta
       real(c_double), value :: s
     end function gsl_sf_zeta
  end interface

contains

  double precision function Zeta_Function(s)
    !!{RST
    Return the Riemann zeta function, :math:`\zeta(s)=\sum_{n=1}^\infty n^{-s}`, for real ``s`` not equal to one.
    !!}
    implicit none
    double precision, intent(in   ) :: s

    Zeta_Function=gsl_sf_zeta(s)
    return
  end function Zeta_Function

end module Zeta_Functions
