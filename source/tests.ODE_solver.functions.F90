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
Contains a module of ODEs for unit tests.
!!}

module Test_ODE_Solver_Functions
  !!{
  Contains ODEs for unit tests.
  !!}
  use :: Interface_GSL, only : GSL_Success
  implicit none
  private
  public :: ODE_Set_1     , ODE_Set_2       , &
       &    Jacobian_Set_1, Jacobian_Set_2  , &
       &                    Integrands_Set_2

contains

  integer function ODE_Set_1(x,y,dydx)
    !!{
    A set of ODEs for unit tests.
    !!}
    double precision              , intent(in   ) :: x
    double precision, dimension(:), intent(in   ) :: y
    double precision, dimension(:), intent(  out) :: dydx
    !$GLC attributes unused :: y

    dydx(1)=sin(x)
    ODE_Set_1=GSL_Success
    return
  end function ODE_Set_1

  integer function Jacobian_Set_1(x,y,dfdy,dfdx)
    !!{
    Jacobian for a set of ODEs for unit tests.
    !!}
    double precision              , intent(in   ) :: x
    double precision, dimension(:), intent(in   ) :: y
    double precision, dimension(:), intent(  out) :: dfdy, dfdx
    !$GLC attributes unused :: y

    dfdy(1)=0.0d0
    dfdx(1)=cos(x)
    Jacobian_Set_1=GSL_Success
    return
  end function Jacobian_Set_1

  integer function ODE_Set_2(x,y,dydx)
    !!{
    A set of ODEs for unit tests.
    !!}
    double precision              , intent(in   ) :: x
    double precision, dimension(:), intent(in   ) :: y
    double precision, dimension(:), intent(  out) :: dydx
    !$GLC attributes unused :: x

    dydx(1)=       y(2)
    dydx(2)=-1.0d0*y(1)
    ODE_Set_2=GSL_Success
    return
  end function ODE_Set_2

  integer function Jacobian_Set_2(x,y,dfdy,dfdx)
    !!{
    Jacobian for a set of ODEs for unit tests.
    !!}
    double precision              , intent(in   ) :: x
    double precision, dimension(:), intent(in   ) :: y
    double precision, dimension(:), intent(  out) :: dfdy, dfdx
    !$GLC attributes unused :: x, y

    dfdy(1  )=+0.0d0
    dfdy(2  )=+1.0d0
    dfdy(3  )=-1.0d0
    dfdy(4  )=+0.0d0
    dfdx(1:2)=+0.0d0
    Jacobian_Set_2=GSL_Success
    return
  end function Jacobian_Set_2

  subroutine Integrands_Set_2(x,y,dydx,z0,e,dzdx)
    !!{
    A set of integrands for unit tests.
    !!}
    double precision, intent(in   ), dimension(        : ) :: x
    double precision, intent(in   ), dimension(:,:) :: y   , dydx
    double precision, intent(in   ), dimension(:        ) :: z0
    logical         , intent(inout), dimension(        : ) :: e
    double precision, intent(  out), dimension(:,:) :: dzdx
    !$GLC attributes unused :: x, dydx, z0

    if (e(1)) dzdx(1,:)=1.0d0/sqrt(1.0d0+y(1,:)**2)
    if (e(2)) dzdx(2,:)=1.0d0/sqrt(1.0d0+y(2,:)**2)
    return
  end subroutine Integrands_Set_2

end module Test_ODE_Solver_Functions
