!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of error functions.

module Error_Functions
  !% Implements calculations of error functions.
  implicit none
  private
  public :: Error_Function, Error_Function_Complementary, Faddeeva

  interface Error_Function
     module procedure Error_Function_Real
     module procedure Error_Function_Complex
  end interface Error_Function

  interface Error_Function_Complementary
     module procedure Error_Function_Complementary_Real
     module procedure Error_Function_Complementary_Complex
  end interface Error_Function_Complementary

contains

  double precision function Error_Function_Real(argument)
    !% Computes the error function.
    use FGSL
    implicit none
    double precision, intent(in   ) :: argument

    Error_Function_Real=FGSL_SF_Erf(argument)
    return
  end function Error_Function_Real

  double precision function Error_Function_Complementary_Real(argument)
    !% Computes the complementary error function.
    use FGSL
    implicit none
    double precision, intent(in   ) :: argument

    Error_Function_Complementary_Real=FGSL_SF_ErfC(argument)
    return
  end function Error_Function_Complementary_Real
  
  elemental double complex function Error_Function_Complex(argument)
    !% Computes the complex complementary error function.
    implicit none
    double complex, intent(in   ) :: argument

    Error_Function_Complex=1.0d0-Error_Function_Complementary_Complex(argument)
    return
  end function Error_Function_Complex

  elemental double complex function Error_Function_Complementary_Complex(argument)
    !% Computes the complex complementary error function, using the algorithm of \cite{abrarov_rapid_2013}.
    use Numerical_Constants_Math
    implicit none
    double complex  , intent(in   ) :: argument
    integer         , parameter     :: N        =23
    double precision, parameter     :: sigma    = 2.0d0
    double precision, parameter     :: tauM     =12.0d0
    integer                         :: i
    double complex                  :: summation       , z
    double precision                :: A               , B

    if (real(argument) < 0.0d0) then
       z=-argument
    else
       z=+argument
    end if
    summation=dcmplx(0.0d0,0.0d0)
    do i=1,N
       A=                      &
            & +2.0d0           &
            & *tauM            &
            & *exp(            &
            &      +sigma  **2 &
            &      -dble(i)**2 &
            &      *Pi     **2 &
            &      /tauM   **2 &
            &     )            &
            & *cos(            &
            &      +2.0d0      &
            &      *dble(i)    &
            &      *Pi         &
            &      *sigma      &
            &      /tauM       &
            &     )
       B=                      &
            & +2.0d0           &
            & *dble(i)         &
            & *Pi              &
            & *exp(            &
            &      +sigma  **2 &
            &      -dble(i)**2 &
            &      *Pi     **2 &
            &      /tauM   **2 &
            &     )            &
            & *sin(            &
            &      +2.0d0      &
            &      *dble(i)    &
            &      *Pi         &
            &      *sigma      &
            &      /tauM       &
            &     )
       summation=                            &
            & +summation                     &
            & +(                             &
            &   +A                           &
            &   *(                           &
            &     +sigma                     &
            &     +z                         &
            &    )                           &
            &   +B                           &
            &  )                             &
            & /(                             &
            &   +dcmplx(dble(i)*Pi,0.0d0)**2 &
            &   +(                           &
            &     +tauM                      &
            &     *(                         &
            &       +sigma                   &
            &       +z                       &
            &      )                         &
            &    )                       **2 &
            &  )
    end do
    Error_Function_Complementary_Complex= &
         & +exp(-z**2)                    &
         & *(                             &
         &   +exp(sigma**2)               &
         &   /tauM                        &
         &   /(                           &
         &     +sigma                     &
         &     +z                         &
         &    )                           &
         &   +summation                   &
         &  )
    if (real(argument) < 0.0d0) Error_Function_Complementary_Complex=2.0d0-Error_Function_Complementary_Complex
    return
  end function Error_Function_Complementary_Complex

  elemental double complex function Faddeeva(argument)
    !% The \href{http://en.wikipedia.org/wiki/Faddeeva_function}{Fadeeva function}.
    implicit none
    double complex  , intent(in   ) :: argument
    
    Faddeeva=exp(-argument**2)*Error_Function_Complementary(-dcmplx(0.0d0,1.0d0)*argument)
    return
  end function Faddeeva
  
end module Error_Functions
