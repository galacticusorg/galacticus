!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements a nonlinear power spectrum class in which the nonlinear power spectrum is computed using the
algorithm of \cite{peacock_non-linear_1996}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Power_Spectra      , only : powerSpectrumClass

  !![
  <powerSpectrumNonlinear name="powerSpectrumNonlinearPeacockDodds1996">
   <description>Provides a nonlinear power spectrum class in which the power spectrum is computed using the algorithm of \cite{peacock_non-linear_1996}.</description>
  </powerSpectrumNonlinear>
  !!]
  type, extends(powerSpectrumNonlinearClass) :: powerSpectrumNonlinearPeacockDodds1996
     !!{
     A nonlinear power spectrum class in which the power spectrum is computed using the algorithm of \cite{peacock_non-linear_1996}.
     !!}
     private
     double precision                         , dimension(2) :: waveNumberPrevious           , fNLPrevious
     double precision                                        :: timePrevious
     class           (cosmologyFunctionsClass), pointer      :: cosmologyFunctions_ => null()
     class           (powerSpectrumClass     ), pointer      :: powerSpectrum_      => null()
  contains
     final     ::               peacockDodds1996Destructor
     procedure :: value      => peacockDodds1996Value
  end type powerSpectrumNonlinearPeacockDodds1996

  interface powerSpectrumNonlinearPeacockDodds1996
     !!{
     Constructors for the \refClass{powerSpectrumNonlinearPeacockDodds1996} nonlinear power spectrum class.
     !!}
     module procedure peacockDodds1996ConstructorParameters
     module procedure peacockDodds1996ConstructorInternal
  end interface powerSpectrumNonlinearPeacockDodds1996

contains

  function peacockDodds1996ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the peacockDodds1996 nonlinear power spectrum class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (powerSpectrumNonlinearPeacockDodds1996)                        :: self
    type (inputParameters                       ), target, intent(inout) :: parameters
    class(cosmologyFunctionsClass               ), pointer               :: cosmologyFunctions_
    class(powerSpectrumClass                    ), pointer               :: powerSpectrum_

    ! Check and read parameters.
    ! Construct required objects.
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="powerSpectrum"      name="powerSpectrum_"      source="parameters"/>
    !!]
    ! Call the internal constructor.
    self=peacockDodds1996ConstructorInternal(cosmologyFunctions_,powerSpectrum_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="powerSpectrum_"     />
    !!]
    return
  end function peacockDodds1996ConstructorParameters

  function peacockDodds1996ConstructorInternal(cosmologyFunctions_,powerSpectrum_) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumNonlinearPeacockDodds1996} nonlinear power spectrum class.
    !!}
    implicit none
    type (powerSpectrumNonlinearPeacockDodds1996)                        :: self
    class(cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    class(powerSpectrumClass                    ), intent(in   ), target :: powerSpectrum_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *powerSpectrum_"/>
    !!]

    ! Initialize state.
    self%waveNumberPrevious=-1.0d0
    self%timePrevious      =-1.0d0
    return
  end function peacockDodds1996ConstructorInternal

  subroutine peacockDodds1996Destructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumNonlinearPeacockDodds1996} nonlinear power spectrum class.
    !!}
    implicit none
    type(powerSpectrumNonlinearPeacockDodds1996), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%powerSpectrum_"     />
    !!]
    return
  end subroutine peacockDodds1996Destructor

  double precision function peacockDodds1996Value(self,wavenumber,time)
    !!{
    Return a nonlinear power spectrum equal using the algorithm of \cite{peacock_non-linear_1996}.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(powerSpectrumNonlinearPeacockDodds1996), intent(inout) :: self
    double precision                             , intent(in   ) :: time                          , wavenumber
    integer                                      , parameter     :: iterationCountMaximum  =1000
    double precision                             , parameter     :: tolerance              =1.0d-3
    double precision                             , parameter     :: updateFraction         =0.5d+0
    logical                                                      :: converged
    integer                                                      :: iterationCount                , i
    double precision                                             :: A                             , B               , &
         &                                                          V                             , alpha           , &
         &                                                          beta                          , fNL             , &
         &                                                          fNLLastIteration              , g               , &
         &                                                          n                             , waveNumberLinear, &
         &                                                          x                             , logRatioBest    , &
         &                                                          omegaMatter                   , omegaDarkEnergy , &
         &                                                          omegaMatterFourSevenths

    ! Pre-compute quantities which depend only on time, and so will be constant throughout this calculation.
    omegaMatter              =self%cosmologyFunctions_%omegaMatterEpochal    (time)
    omegaDarkEnergy          =self%cosmologyFunctions_%omegaDarkEnergyEpochal(time)
    omegaMatterFourSevenths  =omegaMatter                                           **(4.0d0/7.0d0)
    ! Determine if we can use a previous estimate of fNL as our starting guess.
    fNL         =-1.0d0
    logRatioBest=log(2.0d0)
    if (time == self%timePrevious) then
       do i=1,2
          if (self%waveNumberPrevious(i) > 0.0d0 .and. abs(log(wavenumber/self%waveNumberPrevious(i))) < logRatioBest) then
             fNL         =self%fNLPrevious(i)
             logRatioBest=abs(log(wavenumber/self%waveNumberPrevious(i)))
          end if
       end do
    end if
    ! Make an initial guess that the nonlinear power spectrum equals the linear power spectrum.
    if (fNL < 0.0d0) fNL=self%powerSpectrum_%powerDimensionless(wavenumber,time)
    fNLLastIteration=fNL
    ! Iterate until a converged solution is found.
    converged     =.false.
    iterationCount=0
    do while (.not.converged .and. iterationCount < iterationCountMaximum)
       ! Find the corresponding linear wavenumber.
       waveNumberLinear=wavenumber/(1.0d0+fNL)**(1.0d0/3.0d0)
       ! Get the dimensionless linear power spectrum and its logarithmic slope.
       x=self%powerSpectrum_%powerDimensionless        (      waveNumberLinear,time)
       n=self%powerSpectrum_%powerLogarithmicDerivative(0.5d0*waveNumberLinear,time)
       ! Compute parameters of the Peacock & Dodds fitting function.
       A    = 0.482d0/(1.0d0+n/3.0d0)**0.947d0
       B    = 0.226d0/(1.0d0+n/3.0d0)**1.778d0
       alpha= 3.310d0/(1.0d0+n/3.0d0)**0.244d0
       beta = 0.862d0/(1.0d0+n/3.0d0)**0.287d0
       V    =11.550d0/(1.0d0+n/3.0d0)**0.423d0
       ! Compute growth factor using same fitting function as Peacock & Dodds (from Carroll, Press & Turner 1992).
       g    = +2.5d0                       &
            & *    omegaMatter             &
            & /(                           &
            &   +  omegaMatterFourSevenths &
            &   -  omegaDarkEnergy         &
            &   +(                         &
            &     +1.0d0                   &
            &     +0.5d0                   &
            &     *omegaMatter             &
            &    )                         &
            &   *(                         &
            &     +1.0d0                   &
            &     +omegaDarkEnergy         &
            &     /70.0d0                  &
            &    )                         &
            &  )
       ! Compute new estimate of non-linear power-spectrum.
       fNL=    x                                             &
            & *(                                             &
            &    (1.0d0+B*beta*x+(A*x)**(alpha*beta))        &
            &   /(1.0d0+((A*x)**alpha*g**3/V/sqrt(x))**beta) &
            &  )**(1.0d0/beta)
       ! Update our estimate using a mixture of new and old results (this avoids oscillating solutions by preventing large initial
       ! jumps in the estimate).
       fNL=updateFraction*fNL+(1.0d0-updateFraction)*fNLLastIteration
       ! Test for convergence.
       converged=(abs(fNL-fNLLastIteration) < tolerance*0.5d0*(abs(fNL)+abs(fNLLastIteration)))
       ! Move to next iteration.
       fNLLastIteration=fNL
       iterationCount  =iterationCount+1
    end do
    if (.not.converged) call Error_Report('nonlinear power spectrum calculation failed to converge'//{introspection:location})
    ! Store evaluations.
    self%timePrevious         =time
    self%waveNumberPrevious(1)=self%waveNumberPrevious(2)
    self%fNLPrevious       (1)=self%fNLPrevious       (2)
    self%waveNumberPrevious(2)=wavenumber
    self%fNLPrevious       (2)=fNL
    ! Convert to a dimensionful power spectrum.
    peacockDodds1996Value=(2.0d0*Pi)**3*fNL/4.0d0/Pi/wavenumber**3
    return
  end function peacockDodds1996Value
