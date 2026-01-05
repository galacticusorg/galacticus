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

!: $(BUILDPATH)/pFq/pfq.new.o
!: $(BUILDPATH)/gslSpecFuncApprox/hyperg_2F1.o

! Add dependency on GSL library.
!; gsl

!!{
Contains a module which implements hypergeometric functions.
!!}

module Hypergeometric_Functions
  !!{
  Implements hypergeometric functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double     , c_int
  use            :: Interface_GSL, only : gsl_sf_result, gsl_success
  implicit none
  private
  public :: Hypergeometric_1F1            , Hypergeometric_2F1, Hypergeometric_pFq, Hypergeometric_pFq_Regularized, &
       &    Hypergeometric_2F1_Regularized

  interface Hypergeometric_pFq
     module procedure :: Hypergeometric_pFq_Real
     module procedure :: Hypergeometric_pFq_Complex
  end interface Hypergeometric_pFq

  interface
     function gsl_sf_hyperg_2F1_approx_e(a,b,c,x,tol,result) bind(c,name='gsl_sf_hyperg_2F1_approx_e')
       !!{
       Template for the GSL approximate hypergeometric 2F1 C function.
       !!}
       import
       integer(c_int        )        :: gsl_sf_hyperg_2F1_approx_e
       real   (c_double     ), value :: a,b,c,x,tol
       type   (gsl_sf_result)        :: result
     end function gsl_sf_hyperg_2F1_approx_e

     function gsl_sf_hyperg_2F1_e(a,b,c,x,result) bind(c,name='gsl_sf_hyperg_2F1_e')
       !!{
       Template for the GSL hypergeometric 2F1 C function.
       !!}
       import
       integer(c_int        )        :: gsl_sf_hyperg_2F1_e
       real   (c_double     ), value :: a                  , b, &
            &                           c                  , x
       type   (gsl_sf_result)        :: result
     end function gsl_sf_hyperg_2F1_e

     function gsl_sf_hyperg_1F1(a,b,x) bind(c,name='gsl_sf_hyperg_1F1')
       !!{
       Template for the GSL hypergeometric 1F1 C function.
       !!}
       import
       real(c_double)        :: gsl_sf_hyperg_1F1
       real(c_double), value :: a                , b, &
            &                   x
     end function gsl_sf_hyperg_1F1

     function gsl_sf_hyperg_2F1_approx_series(a,b,c,x,tol,result) bind(c,name='gsl_sf_hyperg_2F1_approx_series')
       !!{
       Template for the GSL series approximation of the hypergeometric 2F1 C function.
       !!}
       import
       integer(c_int        )        :: gsl_sf_hyperg_2F1_approx_series
       real   (c_double     ), value :: a,b,c,x,tol
       type   (gsl_sf_result)        :: result
     end function gsl_sf_hyperg_2F1_approx_series
  end interface

  ! Error status.
  integer(c_int) :: statusActual
  !$omp threadprivate(statusActual)

contains

  double precision function Hypergeometric_1F1(a,b,x)
    !!{
    Evaluate the $_1F_1(a_1;b_1;x)$ hypergeometric function.
    !!}
    implicit none
    double precision, intent(in   ) :: a(1), b(1), x

    Hypergeometric_1F1=GSL_SF_Hyperg_1F1(a(1),b(1),x)
    return
  end function Hypergeometric_1F1

  double precision function Hypergeometric_2F1(a,b,x,status,error,toleranceRelative)
    !!{
    Evaluate the $_2F_1(a_1,a_2;b_1;x)$ hypergeometric function.
    !!}
    use :: Error           , only : Error_Report  , GSL_Error_Handler_Abort_Off, GSL_Error_Handler_Abort_On
    use :: Gamma_Functions , only : Gamma_Function
    implicit none
    double precision               , intent(in   )           :: a(2)              , b(1)             , &
         &                                                      x
    integer         (c_int        ), intent(  out), optional :: status
    double precision               , intent(  out), optional :: error
    double precision               , intent(in   ), optional :: toleranceRelative
    double precision               , parameter               :: gsl_dbl_epsilon  =2.2204460492503131d-16
    double precision                                         :: toleranceActual
    double precision                                         :: prefactor1        , prefactor2
    logical                                                  :: aIsNegativeInt (2), bIsNegativeInt(1), &
         &                                                      baIsNegativeInt(2)
    type            (gsl_sf_result)                          :: gslResult

    ! Use our own error handler.
    if (present(status)) then
       call GSL_Error_Handler_Abort_Off()
       statusActual=GSL_Success
    end if
    if (present(toleranceRelative)) then
       toleranceActual=toleranceRelative
    else
       toleranceActual=gsl_dbl_epsilon
    end if
    ! GSL only evaluates this function for |x|<1.
    if (abs(x) <= 1.0d0) then
       ! |x|<1 so simply call the GSL function to compute the function (special cases are treated separately).
       if (.not.present(toleranceRelative) .and. abs(b(1)-a(1)-a(2)-dnint(b(1)-a(1)-a(2))) >= 1000.0d0*gsl_dbl_epsilon) then
          statusActual=GSL_SF_Hyperg_2F1_E(a(1),a(2),b(1),x,gslResult)
          Hypergeometric_2F1=gslResult%val
          if (present(error)) error=gslResult%err
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
          statusActual=GSL_Success
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
          if (statusActual == GSL_Success) then
             statusActual=GSL_SF_Hyperg_2F1_Approx_E(a(2),b(1)-a(1),-a(1)+a(2)+1.0d0,1.0d0/(1.0d0-x),toleranceActual,gslResult)
             Hypergeometric_2F1=Hypergeometric_2F1+gslResult%val*prefactor2
             if (present(error)) error=error+gslResult%err*abs(prefactor2)
          end if
       else
          ! Use a Pfaff transformation to evaluate in terms of a hypergeometric function with |x|<1 (special cases are treated separately).
          if (.not.present(toleranceRelative) .and. abs(a(1)-a(2)-dnint(a(1)-a(2))) >= 1000.0d0*gsl_dbl_epsilon) then
             statusActual=GSL_SF_Hyperg_2F1_E(a(2),b(1)-a(1),b(1),x/(x-1.0d0),gslResult)
             prefactor1=1.0d0/(1.0d0-x)**a(2)
             Hypergeometric_2F1=gslResult%val*prefactor1
             if (present(error)) error=gslResult%err*abs(prefactor1)
          else
             statusActual=GSL_SF_Hyperg_2F1_Approx_E(a(2),b(1)-a(1),b(1),x/(x-1.0d0),toleranceActual,gslResult)
             prefactor1=1.0d0/(1.0d0-x)**a(2)
             Hypergeometric_2F1=gslResult%val*prefactor1
             if (present(error)) error=gslResult%err*abs(prefactor1)
          end if
       end if
    else
       Hypergeometric_2F1=0.0d0
       call Error_Report('function cannot be evaluated for x>1'//{introspection:location})
    end if
    if (present(status)) then
       status=statusActual
       ! Reset error handler.
       call GSL_Error_Handler_Abort_On()
    else if (statusActual /= GSL_Success) then
       call Error_Report('GSL failed'//{introspection:location})
    end if
    return
  end function Hypergeometric_2F1

  double precision function Hypergeometric_2F1_Regularized(a,b,x)
    !!{
    Evaluate the regularized generalized hypergeometric function
    $_2F_1(a_1,a_2;b_1;x)/\Gamma(b_1)$ for real arguments.
    !!}
    implicit none
    double precision, intent(in   ), dimension(2) :: a
    double precision, intent(in   ), dimension(1) :: b
    double precision, intent(in   )               :: x

    Hypergeometric_2F1_Regularized=Hypergeometric_2F1(a,b,x)/Gamma(b(1))
    return
  end function Hypergeometric_2F1_Regularized

  double complex function Hypergeometric_pFq_Complex(a,b,x,toleranceRelative)
    !!{
    Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x)$, using the algorithm of
    \cite{perger_numerical_1993}.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    double complex  , intent(in   ), dimension(:) :: a                , b
    double complex  , intent(in   )               :: x
    double precision, intent(in   ), optional     :: toleranceRelative
    double complex                                :: PFQ
    integer                                       :: LNPFQ            , IX   , &
         &                                           NSIGFIG
    logical                                       :: a1is1            , a1is2, &
         &                                           a2is2            , a2is4, &
         &                                           b1is2            , b1Is3

    LNPFQ  = 0
    IX     = 0
    if (present(toleranceRelative)) then
       NSIGFIG=ceiling(-(log10(toleranceRelative)))
    else
       NSIGFIG=10
    end if
    if (dreal(x) == 0.0d0) then
       Hypergeometric_pFq_Complex=1.0d0
       return
    end if
    if (size(a) == 2 .and. size(b) == 1) then
       ! Special cases for ₂F₁.
       a1Is1=Values_Agree(real(a(1)),1.0d0,absTol=1.0d-6) .and. Values_Agree(imag(a(1)),0.0d0,absTol=1.0d-6)
       a1Is2=Values_Agree(real(a(1)),2.0d0,absTol=1.0d-6) .and. Values_Agree(imag(a(1)),0.0d0,absTol=1.0d-6)
       a2Is2=Values_Agree(real(a(2)),2.0d0,absTol=1.0d-6) .and. Values_Agree(imag(a(2)),0.0d0,absTol=1.0d-6)
       a2Is4=Values_Agree(real(a(2)),4.0d0,absTol=1.0d-6) .and. Values_Agree(imag(a(2)),0.0d0,absTol=1.0d-6)
       b1Is2=Values_Agree(real(b(1)),2.0d0,absTol=1.0d-6) .and. Values_Agree(imag(b(1)),0.0d0,absTol=1.0d-6)
       b1Is3=Values_Agree(real(b(1)),3.0d0,absTol=1.0d-6) .and. Values_Agree(imag(b(1)),0.0d0,absTol=1.0d-6)
       if (a1Is1 .and. a2Is2 .and. b1Is2) then
          ! ₂F₁([1,2],[2],x) = 1/(1-x).
          Hypergeometric_pFq_Complex=1.0d0/(1.0d0-x)
          return
       end if
       if (a1Is2 .and. a2Is2 .and. b1Is3) then
          ! ₂F₁([2,2],[3],x) = 2 (-x-log[1-x]+x log[1-x])/x²/(-1+x)
          Hypergeometric_pFq_Complex=2.0d0*(-x-log(1.0d0-x)+x*log(1.0d0-x))/x**2/(-1.0d0+x)
          return
       end if
       if (a1Is1 .and. a2Is4 .and. b1Is2) then
          ! ₂F₁([1,4],[2],x) = -(3 + x [x-3])/3/(x-1)³
          Hypergeometric_pFq_Complex=-(3.0d0+(x-3.0d0)*x)/3.0d0/(x-1.0d0)**3
          return
       end if
       if (a1Is2 .and. a2Is4 .and. b1Is3) then
          ! ₂F₁([2,4],[3],x) = (x-3)/3/(x-1)³
          Hypergeometric_pFq_Complex=(x-3.0d0)/3.0d0/(x-1.0d0)**3
          return
       end if
    end if
    Hypergeometric_pFq_Complex=PFQ(a,b,size(a),size(b),x,LNPFQ,IX,NSIGFIG)    
    return
  end function Hypergeometric_pFq_Complex

  double precision function Hypergeometric_pFq_Real(a,b,x,toleranceRelative,useAcceleration)
    !!{
    Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x)$ for real arguments.
    !!}
    use :: Error           , only : Error_Report
    use :: Gamma_Functions , only : Gamma_Function
    implicit none
    double precision, intent(in   ), dimension(:) :: a                    , b
    double precision, intent(in   )               :: x
    double precision, intent(in   ), optional     :: toleranceRelative
    logical         , intent(in   ), optional     :: useAcceleration
    double precision                              :: toleranceActual
    double precision                              :: aP(size(a))          , bQ(size(b))
    double precision                              :: aMulti               , bMulti     , &
         &                                           aaMulti              , baMulti    , &
         &                                           term
    integer                                       :: p                    , q          , &
         &                                           i                    , j          , &
         &                                           k
    logical                                       :: useAccelerationActual

    if (present(toleranceRelative)) then
       toleranceActual=toleranceRelative
    else
       toleranceActual=1.0d-10
    end if
    if (present(useAcceleration)) then
       useAccelerationActual= useAcceleration
    else
       useAccelerationActual=.false.
    end if
    p=size(a)
    q=size(b)
    if (x >= 1.0d0) then
       Hypergeometric_pFq_Real=real(Hypergeometric_pFq_Complex(dcmplx(a),dcmplx(b),dcmplx(x),toleranceActual))
    else if (x >= -1.0d0) then
       if (useAccelerationActual) then
          Hypergeometric_pFq_Real=Hypergeometric_pFq_approx_series(a,b,x,toleranceActual)
       else
          Hypergeometric_pFq_Real=real(Hypergeometric_pFq_Complex(dcmplx(a),dcmplx(b),dcmplx(x),toleranceActual))
       end if
    else
       if (p == q+1) then
          ! Use the transformation x -> 1/x to evaluate in terms of a hypergeometric function with |x|<1.
          ! Note that the transformation is only valid when a(j), b(j), a(j)-a(k) (|j-k|>0) are not zero
          ! negative integers. For details, see
          ! <http://functions.wolfram.com/07.31.06.0019.01>.
          ! First, check whether any of a(i) is a nonpositive integer.
          if (any(a <= 0.0d0 .and. a-dnint(a) == 0.0d0)) then
             Hypergeometric_pFq_Real=Hypergeometric_pFq_approx_series(a,b,x,toleranceActual)
             return
          end if
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
             if (useAccelerationActual) then
                term=Hypergeometric_pFq_approx_series(aP,bQ,1.0d0/x,toleranceActual)
             else
                term=real(Hypergeometric_pFq_Complex(dcmplx(aP),dcmplx(bQ),dcmplx(1.0d0/x),toleranceActual))
             end if
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
          call Error_Report('not implemented yet'//{introspection:location})
       end if
    end if
    return
  end function Hypergeometric_pFq_Real

  double precision function Hypergeometric_pFq_approx_series(a,b,x,toleranceRelative)
    !!{
    Evaluate the generalized hypergeometric function $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x)$ by direct summation.
    Shanks transformation (\cite{shanks_non_linear_1955}) is used to accelerate the calculations.
    !!}
    use :: Error             , only : Error_Report
    use :: Display           , only : displayMessage, verbosityLevelWarn
    use :: ISO_Varying_String, only : varying_string, assignment(=)     , operator(//)
    implicit none
    double precision                , intent(in   ), dimension(:) :: a                 , b
    double precision                , intent(in   )               :: x
    double precision                , intent(in   ), optional     :: toleranceRelative
    double precision                                              :: toleranceActual
    double precision                , parameter                   :: gsl_dbl_epsilon=2.2204460492503131d-16
    double precision                                              :: An                , Sn                , &
         &                                                           SnPrevious        , delta             , &
         &                                                           deltaPrevious     , deltaSn           , &
         &                                                           deltaSnPrevious   , denominator       , &
         &                                                           numerator         , correction        , &
         &                                                           correctionPrevious
    integer                         , parameter                   :: nMax=50000
    integer                                                       :: i
    type            (varying_string)                              :: message
    character       (len=13        )                              :: label

    if (present(toleranceRelative)) then
       toleranceActual=toleranceRelative
    else
       toleranceActual=gsl_dbl_epsilon
    end if
    An        =1.0d0
    Sn        =1.0d0
    delta     =1.0d0
    deltaSn   =0.0d0
    correction=0.0d0
    deltaPrevious     =1.0d0
    deltaSnPrevious   =0.0d0
    correctionPrevious=0.0d0
    do i=1,nMax
       SnPrevious=Sn
       if (i >= 2) deltaPrevious     =delta
       if (i >= 3) correctionPrevious=correction
       if (i >= 4) deltaSnPrevious   =deltaSn
       delta=delta*product(a+i-1.0d0)/product(b+i-1.0d0)*x/i
       An   =An+delta
       Sn   =An
       ! Perform a Shanks transformation on An.
       correction=0.0d0
       if (i >= 2) then
          denominator=delta-deltaPrevious
          if (denominator .ne. 0.0d0) then
             numerator = delta**2
             correction=-numerator/denominator
             Sn=An+correction
          end if
       end if
       ! If x<0, perform a second Shanks transformation on Sn. This is very useful when x<0, which will make the
       ! result converge much faster.
       if (x < 0.0d0) then
          if (i >= 3) then
             deltaSn=delta+correction-correctionPrevious
          end if
          if (i >= 4) then
             denominator=deltaSn-deltaSnPrevious
             if (denominator .ne. 0.0d0) then
                numerator=deltaSn**2
                Sn=Sn-numerator/denominator
             end if
             ! If the increment of the series in one step is larger than the asymptotic value of the series,
             ! roundoff error may prevent the result from achieving specified relative tolerance. This usually
             ! happens when x is close to -1. In such case, try to increase the tolerance.
             if (abs(delta/Sn) > 1.0d0) toleranceActual=max(abs(delta/Sn)*gsl_dbl_epsilon, toleranceActual)
          end if
       end if
       ! Check the convergence of the result.
       if (abs(Sn-SnPrevious) < toleranceActual*abs(Sn)) exit
    end do
    if (i <= nMax) then
       Hypergeometric_pFq_approx_series=Sn
       if (present(toleranceRelative)) then
          if (toleranceActual > toleranceRelative) then
             message='Warning: roundoff error prevents the result from achieving specified relative tolerance ['
             write (label,'(e12.6)') toleranceRelative
             message=message//trim(label)//']. A relative tolerance of ['
             write (label,'(e12.6)') toleranceActual
             message=message//trim(label)//'] is used instead.'//char(10)
             call displayMessage(message//{introspection:location},verbosityLevelWarn)
          end if
       end if
    else
       Hypergeometric_pFq_approx_series=0.0d0
       call Error_Report('maximum number of iterations is reached before the relative tolerence is satisfied. Try to increase the value of [nMax]'//{introspection:location})
    end if
    return
  end function Hypergeometric_pFq_approx_series

  double precision function Hypergeometric_pFq_Regularized(a,b,x)
    !!{
    Evaluate the regularized generalized hypergeometric function
    $_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;x)/[\Gamma(b_1)\ldots\Gamma(b_q)]$ for real arguments.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: a, b
    double precision, intent(in   )               :: x

    Hypergeometric_pFq_Regularized=real(Hypergeometric_pFq(dcmplx(a),dcmplx(b),dcmplx(x)))/product(Gamma(b))
    return
  end function Hypergeometric_pFq_Regularized

end module Hypergeometric_Functions
