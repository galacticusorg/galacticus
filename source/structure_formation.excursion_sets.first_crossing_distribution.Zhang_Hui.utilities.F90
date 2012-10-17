!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements various utility functions for excursion set first crossing distribution set calculations
!% that use variants of the methodology of \cite{zhang_random_2006}.

module Excursion_Sets_First_Crossing_Zhang_Hui_Utilities
  !% Implements various utility functions for excursion set first crossing distribution set calculations
  !% that use variants of the methodology of \cite{zhang_random_2006}.
  implicit none
  private
  public :: g_1,g_2,g_2_Integrated,Delta

  ! Module variables used in integrations.
  double precision :: varianceGlobal,timeGlobal,barrierGlobal,barrierGradientGlobal
  !$omp threadprivate(varianceGlobal,timeGlobal,barrierGlobal,barrierGradientGlobal)

contains
  
  double precision function g_1(variance,time)
    !% Returns the function $g_1(S)$ in the \cite{zhang_random_2006} algorithm for excursion set barrier crossing probabilities.
    use Excursion_Sets_Barriers
    use Math_Distributions_Gaussian
    implicit none
    double precision, intent(in) :: variance,time
    double precision             :: barrier

    barrier=Excursion_Sets_Barrier(variance,time)
    g_1=(barrier/variance-2.0d0*Excursion_Sets_Barrier_Gradient(variance,time))*Gaussian_Distribution(barrier,sqrt(variance))
    return
  end function g_1
  
  double precision function g_2(variance,variancePrimed,time)
    !% Returns the function $g_2(S,S^\prime)$ in the \cite{zhang_random_2006} algorithm for excursion set barrier crossing probabilities.
    use Excursion_Sets_Barriers
    use Math_Distributions_Gaussian
    implicit none  
    double precision, intent(in) :: variance,variancePrimed,time
    double precision             :: barrierPrimed
    double precision, save       :: variancePrevious=-1.0d0,timePrevious=-1.0d0,barrier,barrierGradient
    !$omp threadprivate(variancePrevious,barrier,barrierGradient)
    
    ! Compute the barriers.
    if (variance /= variancePrevious .or. time /= timePrevious) then
       variancePrevious=variance
       timePrevious    =time
       barrier         =Excursion_Sets_Barrier         (variance,time)
       barrierGradient =Excursion_Sets_Barrier_Gradient(variance,time)
    end if    
    barrierPrimed=Excursion_Sets_Barrier(variancePrimed,time)
    ! Compute the function.
    g_2=(2.0d0*barrierGradient-(barrier-barrierPrimed)/(variance-variancePrimed))*Gaussian_Distribution(barrier-barrierPrimed&
         &,sqrt(variance-variancePrimed))
    return
  end function g_2

  double precision function g_2_Integrated(variance,deltaVariance,time)
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Comparison
    use Numerical_Integration
    use Excursion_Sets_Barriers
    implicit none
    double precision,                 intent(in) :: variance,deltaVariance,time
    double precision,                 parameter  :: gradientChangeTolerance=1.0d-3
    double precision                             :: smallStep
    type(c_ptr)                                  :: parameterPointer
    type(fgsl_function)                          :: integrandFunction
    type(fgsl_integration_workspace)             :: integrationWorkspace
    
    varianceGlobal       =variance
    timeGlobal           =time
    barrierGlobal        =Excursion_Sets_Barrier         (varianceGlobal,timeGlobal)
    barrierGradientGlobal=Excursion_Sets_Barrier_Gradient(varianceGlobal,timeGlobal)
    
    ! Find a suitably small step in variance that allows us to compute the divergent part of the integral with an analytic
    ! approximation. The approximation used assumes that the barrier gradient, dB/dS, is constant, so find a step over which
    ! the gradient is constant to within a specified tolerance.
    smallStep=deltaVariance
    do while (Values_Differ(Excursion_Sets_Barrier_Gradient(variance-smallStep,timeGlobal),barrierGradientGlobal&
         &,relTol=gradientChangeTolerance))
       smallStep=0.5d0*smallStep
    end do
    ! Compute the non-divergent part of the integral numerically.
    g_2_Integrated=Integrate(variance-deltaVariance,variance-smallStep,g_2_Integrand_Zhang_Hui,parameterPointer&
         &,integrandFunction,integrationWorkspace ,toleranceAbsolute=1.0d-50,toleranceRelative=1.0d-6,hasSingularities=.true.&
         &,integrationRule=FGSL_Integ_Gauss15)
    ! Compute the divergent part of the integral with an analytic approximation.
    g_2_Integrated=g_2_Integrated+erf(barrierGradientGlobal*sqrt(0.5d0*smallStep))
    call Integrate_Done(integrandFunction,integrationWorkspace)
    return
  end function g_2_Integrated

  function g_2_Integrand_Zhang_Hui(variance,parameterPointer) bind(c)
    !% Integrand function used in computing $\Delta_{i,i}$ in the \cite{zhang_random_2006} algorithm for excursion set barrier
    !% crossing probabilities.
    use Excursion_Sets_Barriers
    use Math_Distributions_Gaussian
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: g_2_Integrand_Zhang_Hui
    real(c_double),   value :: variance
    type(c_ptr),      value :: parameterPointer
    real(c_double)          :: barrier

    if (variance >= varianceGlobal) then
       g_2_Integrand_Zhang_Hui=0.0d0
    else
       barrier=Excursion_Sets_Barrier(variance,timeGlobal)
       g_2_Integrand_Zhang_Hui=(2.0d0*barrierGradientGlobal-(barrierGlobal-barrier)/(varianceGlobal-variance))&
            &*Gaussian_Distribution(barrierGlobal-barrier,sqrt(varianceGlobal-variance))
    end if
    return
  end function g_2_Integrand_Zhang_Hui

  double precision function Delta(i,j,iVariance,jVariance,deltaVariance,time)
    !% Returns the factor $Delta_{i,j}$ in the \cite{zhang_random_2006} algorithm for excursion set barrier crossing probabilities.
    use Excursion_Sets_Barriers
    implicit none
    integer,          intent(in) :: i,j
    double precision, intent(in) :: iVariance,jVariance,deltaVariance,time
    double precision, save       :: timePrevious=-1.0d0,deltaPrevious
    integer,          save       :: iPrevious=-1,jPrevious=-1
    !$omp threadprivate(timePrevious,deltaPrevious,iPrevious,jPrevious)
 
    if (.not.(i == iPrevious .and. j == jPrevious .and. time == timePrevious)) then
       iPrevious   =i
       jPrevious   =j
       timePrevious=time
       if ( i == j ) then
          ! In this case integrate over the range to get an average value.
          deltaPrevious=0.5d0*g_2_Integrated(iVariance,deltaVariance,time)
       else
          ! Compute the appropriate factor.
          deltaPrevious=(deltaVariance/2.0d0)*g_2(iVariance,jVariance-deltaVariance/2.0d0,time)
       end if
    end if
    Delta=deltaPrevious
    return
  end function Delta

end module Excursion_Sets_First_Crossing_Zhang_Hui_Utilities
