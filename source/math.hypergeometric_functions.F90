!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements hypergeometric functions.

module Hypergeometric_Functions
  !% Implements hypergeometric functions.
  use            :: FGSL         , only : FGSL_SF_Hyperg_1F1, FGSL_SF_Hyperg_2F1_E, FGSL_Success, fgsl_int, &
          &                               fgsl_sf_result    , gsl_sf_result
  use, intrinsic :: ISO_C_Binding, only : c_double          , c_int
  implicit none
  private
  public :: Hypergeometric_1F1, Hypergeometric_2F1, Hypergeometric_pFq

  interface Hypergeometric_pFq
     module procedure :: Hypergeometric_pFq_Real
     module procedure :: Hypergeometric_pFq_Complex
  end interface Hypergeometric_pFq

  interface
     function gsl_sf_hyperg_2F1_approx_e(a,b,c,x,tol,result) bind(c,name='gsl_sf_hyperg_2F1_approx_e')
       !% Template for the GSL hypergeometric 2F1 C function.
       import
       integer(c_int        )        :: gsl_sf_hyperg_2F1_approx_e
       real   (c_double     ), value :: a,b,c,x,tol
       type   (gsl_sf_result)        :: result
     end function gsl_sf_hyperg_2F1_approx_e
  end interface

  interface
     function gsl_sf_hyperg_2F1_approx_series(a,b,c,x,tol,result) bind(c,name='gsl_sf_hyperg_2F1_approx_series')
       !% Template for the GSL series approximation of the hypergeometric 2F1 C function.
       import
       integer(c_int        )        :: gsl_sf_hyperg_2F1_approx_series
       real   (c_double     ), value :: a,b,c,x,tol
       type   (gsl_sf_result)        :: result
     end function gsl_sf_hyperg_2F1_approx_series
  end interface

  ! Error status.
  integer(fgsl_int) :: statusActual
  !$omp threadprivate(statusActual)

contains

  double precision function Hypergeometric_1F1(a,b,x)
    !% Evaluate the $_1F_1(a_1;b_1;x)$ hypergeometric function.
    implicit none
    double precision, intent(in   ) :: a(1), b(1), x

    Hypergeometric_1F1=FGSL_SF_Hyperg_1F1(a(1),b(1),x)
    return
  end function Hypergeometric_1F1

  double precision function Hypergeometric_2F1(a,b,x,status,error,toleranceRelative)
    !% Evaluate the $_2F_1(a_1,a_2;b_1;x)$ hypergeometric function.
    use :: Galacticus_Error, only : Galacticus_Error_Report, Galacticus_GSL_Error_Handler_Abort_Off, Galacticus_GSL_Error_Handler_Abort_On
    use :: Gamma_Functions , only : Gamma_Function
    implicit none
    double precision                , intent(in   )           :: a(2)              , b(1)             , &
         &                                                       x
    integer         (fgsl_int      ), intent(  out), optional :: status
    double precision                , intent(  out), optional :: error
    double precision                , intent(in   ), optional :: toleranceRelative
    double precision                , parameter               :: gsl_dbl_epsilon  =2.2204460492503131d-16
    double precision                                          :: toleranceActual
    double precision                                          :: prefactor1        , prefactor2
    logical                                                   :: aIsNegativeInt (2), bIsNegativeInt(1), &
         &                                                       baIsNegativeInt(2)
    type            (fgsl_sf_result)                          :: fgslResult
    type            ( gsl_sf_result)                          ::  gslResult

    ! Use our own error handler.
    if (present(status)) then
       call Galacticus_GSL_Error_Handler_Abort_Off()
       statusActual=FGSL_Success
    end if
    if (present(toleranceRelative)) then
       toleranceActual=toleranceRelative
    else
       toleranceActual=gsl_dbl_epsilon
    end if
    ! GSL only evaluates this function for |x|<1.
    if (abs(x) <= 1.0d0) then
       ! |x|<1 so simply call the GSL function to compute the function.
       if (.not.present(toleranceRelative) .and. abs(b(1)-a(1)-a(2)-dnint(b(1)-a(1)-a(2))) >= 1000.0d0*gsl_dbl_epsilon) then
          statusActual=FGSL_SF_Hyperg_2F1_E(a(1),a(2),b(1),x,fgslResult)
          Hypergeometric_2F1=fgslResult%val
          if (present(error)) error=fgslResult%err
       else
          statusActual=GSL_SF_Hyperg_2F1_Approx_E(a(1),a(2),b(1),x,toleranceActual,gslResult)
          Hypergeometric_2F1=gslResult%val
          if (present(error)) error=gslResult%err
       end if
    else if (x < -1.0d0) then
       ! x<-1 so transformation of variables is needed so that the function can be evaluated
       ! in terms of a hypergeometric function with |x|<1.
       ! First, check if any of the parameters is a nonpositive integer.
       aIsNegativeInt =(a      < 0.0d0 .and. a     -dnint(a     ) == 0.0d0)
       bIsNegativeInt =(b      < 0.0d0 .and. b     -dnint(b     ) == 0.0d0)
       baIsNegativeInt=(b(1)-a < 0.0d0 .and. b(1)-a-dnint(b(1)-a) == 0.0d0)
       if (any(aIsNegativeInt) .or. any(a == 0.0d0)) then
          ! Parameter a(1) or a(2) is a nonpositive integer, evaluate the finite series directly.
          statusActual=GSL_SF_Hyperg_2F1_Approx_Series(a(1),a(2),b(1),x,toleranceActual,gslResult)
          Hypergeometric_2F1=gslResult%val
          if (present(error)) error=gslResult%err
       else if (abs(b(1)-a(1)) == 0.0d0 .or. abs(b(1)-a(2)) == 0.0d0) then
          ! Algebraic formula is available.
          Hypergeometric_2F1=(1.0d0-x)**(b(1)-a(1)-a(2))
          statusActual=FGSL_Success
          if (present(error)) error=gsl_dbl_epsilon*abs(Hypergeometric_2F1)
       else if (any(baIsNegativeInt)) then
          ! Do an Euler transformation,
          ! 2F1(a,b,c,x)=(1-x)**(c-a-b)*2F1(c-a,c-b,c,x),
          ! so that the function can be evaluated from a finite series.
          ! See <http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/17/02/07/0001/>.
          statusActual=GSL_SF_Hyperg_2F1_Approx_Series(b(1)-a(1),b(1)-a(2),b(1),x,toleranceActual,gslResult)
          prefactor1=(1-x)**(b(1)-a(1)-a(2))
          Hypergeometric_2F1=gslResult%val*prefactor1
          if (present(error)) error=gslResult%err*abs(prefactor1)
       else if ((a(1)-a(2)-dnint(a(1)-a(2)) .ne. 0.0d0) .and. (.not.bIsNegativeInt(1)) .and. (b(1) .ne. 0.0d0)) then
          ! Use the transformation x -> 1/(1-x) to evaluate in terms of a hypergeometric function with |x|<1.
          ! See <http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/17/02/07/0007/>.
          prefactor1=+1.0d0                     &
               &     /(1.0d0-x)**a(1)           &
               &     *Gamma_Function(a(2)-a(1)) &
               &     *Gamma_Function(b(1)     ) &
               &     /Gamma_Function(a(2)     ) &
               &     /Gamma_Function(b(1)-a(1))
          prefactor2=+1.0d0                     &
               &     /(1.0d0-x)**a(2)           &
               &     *Gamma_Function(a(1)-a(2)) &
               &     *Gamma_Function(b(1)     ) &
               &     /Gamma_Function(a(1)     ) &
               &     /Gamma_Function(b(1)-a(2))
          ! First term.
          statusActual=GSL_SF_Hyperg_2F1_Approx_E(a(1),b(1)-a(2),a(1)-a(2)+1.0d0,1.0d0/(1.0d0-x),toleranceActual,gslResult)
          Hypergeometric_2F1=gslResult%val*prefactor1
          if (present(error)) error=gslResult%err*abs(prefactor1)
          ! Second term.
          if (statusActual == FGSL_Success) then
             statusActual=GSL_SF_Hyperg_2F1_Approx_E(a(2),b(1)-a(1),-a(1)+a(2)+1.0d0,1.0d0/(1.0d0-x),toleranceActual,gslResult)
             Hypergeometric_2F1=Hypergeometric_2F1+gslResult%val*prefactor2
             if (present(error)) error=error+gslResult%err*abs(prefactor2)
          end if
       else
          ! Use a Pfaff transformation to evaluate in terms of a hypergeometric function with |x|<1.
          if (.not.present(toleranceRelative) .and. abs(a(1)-a(2)-dnint(a(1)-a(2))) >= 1000.0d0*gsl_dbl_epsilon) then
             statusActual=FGSL_SF_Hyperg_2F1_E(a(2),b(1)-a(1),b(1),x/(x-1.0d0),fgslResult)
             prefactor1=1.0d0/(1.0d0-x)**a(2)
             Hypergeometric_2F1=fgslResult%val*prefactor1
             if (present(error)) error=fgslResult%err*abs(prefactor1)
          else
             statusActual=GSL_SF_Hyperg_2F1_Approx_E(a(2),b(1)-a(1),b(1),x/(x-1.0d0),toleranceActual,gslResult)
             prefactor1=1.0d0/(1.0d0-x)**a(2)
             Hypergeometric_2F1=gslResult%val*prefactor1
             if (present(error)) error=gslResult%err*abs(prefactor1)
          end if
       end if
    else
       Hypergeometric_2F1=0.0d0
       call Galacticus_Error_Report('function cannot be evaluated for x>1'//{introspection:location})
    end if
    if (present(status)) then
       status=statusActual
       ! Reset error handler.
       call Galacticus_GSL_Error_Handler_Abort_On()
    else if (statusActual /= FGSL_Success) then
       call Galacticus_Error_Report('GSL failed'//{introspection:location})
    end if
    return
  end function Hypergeometric_2F1

  double complex function Hypergeometric_pFq_Complex(a,b,x,toleranceRelative)
    !% Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x)$, using the algorithm of
    !% \cite{perger_numerical_1993}.
    implicit none
    double complex  , intent(in   ), dimension(:) :: a                , b
    double complex  , intent(in   )               :: x
    double precision, intent(in   ), optional     :: toleranceRelative
    double complex                                :: PFQ
    integer                                       :: LNPFQ            , IX, NSIGFIG

    LNPFQ  = 0
    IX     = 0
    if (present(toleranceRelative)) then
       NSIGFIG=ceiling(-(log10(toleranceRelative)))
    else
       NSIGFIG=10
    end if
    if (dreal(x) == 0.0d0) then
       Hypergeometric_pFq_Complex=1.0d0
    else
       Hypergeometric_pFq_Complex=PFQ(a,b,size(a),size(b),x,LNPFQ,IX,NSIGFIG)
    end if
    return
  end function Hypergeometric_pFq_Complex

  double precision function Hypergeometric_pFq_Real(a,b,x,toleranceRelative)
    !% Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x)$ for real arguments.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Gamma_Functions , only : Gamma_Function
    implicit none
    double precision, intent(in   ), dimension(:) :: a                , b
    double precision, intent(in   )               :: x
    double precision, intent(in   ), optional     :: toleranceRelative
    double precision                              :: toleranceActual
    double precision                              :: aP(size(a))      , bQ(size(b))
    double precision                              :: aMulti           , bMulti     , &
         &                                           aaMulti          , baMulti    , &
         &                                           term
    integer                                       :: p                , q          , &
         &                                           i                , j          , &
         &                                           k

    if (present(toleranceRelative)) then
       toleranceActual=toleranceRelative
    else
       toleranceActual=1.0d-10
    end if
    p=size(a)
    q=size(b)
    if (x >= -1.0d0) then
       Hypergeometric_pFq_Real=real(Hypergeometric_pFq_Complex(dcmplx(a),dcmplx(b),dcmplx(x),toleranceActual))
    else
       if (p == q+1) then
          ! Use the transformation x -> 1/x to evaluate in terms of a hypergeometric function with |x|<1.
          ! Note that the transformation is only valid when a(j), b(j), a(j)-a(k) (|j-k|>0) are not zero
          ! negative integers. For details, see
          ! <http://functions.wolfram.com/07.31.06.0019.01>.
          aMulti=1.0d0
          bMulti=1.0d0
          do i=1, p
             aMulti=aMulti*Gamma_Function(a(i))
          end do
          do i=1, q
             bMulti=bMulti*Gamma_Function(b(i))
          end do
          Hypergeometric_pFq_Real=0.0d0
          do k=1, p
             aaMulti=1.0d0
             baMulti=1.0d0
             do j=1, p
                if (j .ne. k) aaMulti=aaMulti*Gamma_Function(a(j)-a(k))
             end do
             do j=1, q
                ! Check whether b(j)-a(k) is a nonpositive integer.
                if (b(j)-a(k) <= 0.0d0 .and. (b(j)-a(k)-dnint(b(j)-a(k))) == 0.0d0) then
                   aaMulti=0.0d0
                   baMulti=1.0d0
                else
                   baMulti=baMulti*Gamma_Function(b(j)-a(k))
                end if
             end do
             aP(1 )=a(k)
             aP(2:)=a(k)-b+1.0d0
             do i=1, q
                if (i < k) then
                   bQ(i)=1.0d0-a(i  )+a(k)
                else
                   bQ(i)=1.0d0-a(i+1)+a(k)
                end if
             end do
             term=real(Hypergeometric_pFq_Complex(dcmplx(aP),dcmplx(bQ),dcmplx(1.0d0/x),toleranceActual))
             Hypergeometric_pFq_Real=+Hypergeometric_pFq_Real &
                  &                  +Gamma_Function(a(k))    &
                  &                  *aaMulti                 &
                  &                  /baMulti                 &
                  &                  /(-x)**a(k)              &
                  &                  *term
          end do
          Hypergeometric_pFq_Real=bMulti/aMulti*Hypergeometric_pFq_Real
       else
          Hypergeometric_pFq_Real=0.0d0
          call Galacticus_Error_Report('not implemented yet'//{introspection:location})
       end if
    end if
    return
  end function Hypergeometric_pFq_Real

end module Hypergeometric_Functions
