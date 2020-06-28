!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!: $(BUILDPATH)/pFq/pfq.new.o
!: $(BUILDPATH)/gslSpecFuncApprox/hyperg_2F1.o

! Add dependency on GSL library.
!; gsl

!% Contains a module which implements hypergeometric functions.

module Hypergeometric_Functions
  !% Implements hypergeometric functions.
  use, intrinsic :: ISO_C_Binding, only : c_double     , c_int
  use            :: Interface_GSL, only : gsl_sf_result, gsl_success
  implicit none
  private
  public :: Hypergeometric_1F1, Hypergeometric_2F1, Hypergeometric_pFq

  interface Hypergeometric_pFq
     module procedure :: Hypergeometric_pFq_Real
     module procedure :: Hypergeometric_pFq_Complex
  end interface Hypergeometric_pFq

  interface
     function gsl_sf_hyperg_2F1_approx_e(a,b,c,x,tol,result) bind(c,name='gsl_sf_hyperg_2F1_approx_e')
       !% Template for the GSL approximate hypergeometric 2F1 C function.
       import
       integer(c_int        )        :: gsl_sf_hyperg_2F1_approx_e
       real   (c_double     ), value :: a,b,c,x,tol
       type   (gsl_sf_result)        :: result
     end function gsl_sf_hyperg_2F1_approx_e

     function gsl_sf_hyperg_2F1_e(a,b,c,x,result) bind(c,name='gsl_sf_hyperg_2F1_e')
       !% Template for the GSL hypergeometric 2F1 C function.
       import
       integer(c_int        )        :: gsl_sf_hyperg_2F1_e
       real   (c_double     ), value :: a                  , b, &
            &                           c                  , x
       type   (gsl_sf_result)        :: result
     end function gsl_sf_hyperg_2F1_e

     function gsl_sf_hyperg_1F1(a,b,x) bind(c,name='gsl_sf_hyperg_1F1')
       !% Template for the GSL hypergeometric 1F1 C function.
       import
       real(c_double)        :: gsl_sf_hyperg_1F1
       real(c_double), value :: a                , b, &
            &                   x
     end function gsl_sf_hyperg_1F1
  end interface

  ! Error status.
  integer(c_int) :: statusActual
  !$omp threadprivate(statusActual)

contains

  double precision function Hypergeometric_1F1(a,b,x)
    !% Evaluate the $_1F_1(a_1;b_1;x)$ hypergeometric function.
    implicit none
    double precision, intent(in   ) :: a(1), b(1), x

    Hypergeometric_1F1=GSL_SF_Hyperg_1F1(a(1),b(1),x)
    return
  end function Hypergeometric_1F1

  double precision function Hypergeometric_2F1(a,b,x,status,error,toleranceRelative)
    !% Evaluate the $_2F_1(a_1,a_2;b_1;x)$ hypergeometric function.
    use :: Galacticus_Error, only : Galacticus_Error_Report, Galacticus_GSL_Error_Handler_Abort_Off, Galacticus_GSL_Error_Handler_Abort_On
    implicit none
    double precision               , intent(in   )           :: a(2)             , b(1), &
         &                                                      x
    integer         (c_int        ), intent(  out), optional :: status
    double precision               , intent(  out), optional :: error
    double precision               , intent(in   ), optional :: toleranceRelative
    type            (gsl_sf_result)                          :: gslResult

    ! Use our own error handler.
    if (present(status)) then
       call Galacticus_GSL_Error_Handler_Abort_Off()
       statusActual=GSL_Success
    end if
    ! GSL only evaluates this function for |x|<1.
    if (abs(x) <= 1.0d0) then
       ! |x|<1 so simply call the GSL function to compute the function.
       if (.not.present(toleranceRelative)) then
          statusActual=GSL_SF_Hyperg_2F1_E(a(1),a(2),b(1),x,gslResult)
          Hypergeometric_2F1=gslResult%val
          if (present(error)) error=gslResult%err
       else
          statusActual=GSL_SF_Hyperg_2F1_Approx_E(a(1),a(2),b(1),x,toleranceRelative,gslResult)
          Hypergeometric_2F1=gslResult%val
          if (present(error)) error=gslResult%err
       end if
    else if (x < -1.0d0) then
       ! x<-1 so use a Pfaff transformation to evaluate in terms of a hypergeometric function with |x|<1.
       if (.not.present(toleranceRelative)) then
          statusActual=GSL_SF_Hyperg_2F1_E(a(2),b(1)-a(1),b(1),x/(x-1.0d0),gslResult)
          Hypergeometric_2F1=gslResult%val/(1.0d0-x)**a(2)
          if (present(error)) error=gslResult%err/(1.0d0-x)**a(2)
       else
          statusActual=GSL_SF_Hyperg_2F1_Approx_E(a(2),b(1)-a(1),b(1),x/(x-1.0d0),toleranceRelative,gslResult)
          Hypergeometric_2F1=gslResult%val/(1.0d0-x)**a(2)
          if (present(error)) error=gslResult%err/(1.0d0-x)**a(2)
       end if
    else
       Hypergeometric_2F1=0.0d0
       call Galacticus_Error_Report('function cannot be evaluated for x>1'//{introspection:location})
    end if
    if (present(status)) then
       status=statusActual
       ! Reset error handler.
       call Galacticus_GSL_Error_Handler_Abort_On()
    else if (statusActual /= GSL_Success) then
       call Galacticus_Error_Report('GSL failed'//{introspection:location})
    end if
    return
  end function Hypergeometric_2F1

  double complex function Hypergeometric_pFq_Complex(a,b,x)
    !% Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x)$, using the algorithm of
    !% \cite{perger_numerical_1993}.
    implicit none
    double complex, intent(in   ), dimension(:) :: a    , b
    double complex, intent(in   )               :: x
    double complex                              :: PFQ
    integer                                     :: LNPFQ, IX, NSIGFIG

    LNPFQ  = 0
    IX     = 0
    NSIGFIG=10
    if (dreal(x) == 0.0d0) then
       Hypergeometric_pFq_Complex=1.0d0
    else
       Hypergeometric_pFq_Complex=PFQ(a,b,size(a),size(b),x,LNPFQ,IX,NSIGFIG)
    end if
    return
  end function Hypergeometric_pFq_Complex

  double precision function Hypergeometric_pFq_Real(a,b,x)
    !% Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x)$ for real arguments.
    implicit none
    double precision, intent(in   ), dimension(:) :: a  , b
    double precision, intent(in   )               :: x

    Hypergeometric_pFq_Real=real(Hypergeometric_pFq_Complex(dcmplx(a),dcmplx(b),dcmplx(x)))
    return
  end function Hypergeometric_pFq_Real

end module Hypergeometric_Functions
