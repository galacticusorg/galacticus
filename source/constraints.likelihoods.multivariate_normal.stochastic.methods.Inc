!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines methods for a multivariate normal likelihood function with stochastic sampling.

function likelihoodMultivariateNormalStochasticConstructor(means,covariance,realizationCount,realizationCountMinimum)
  !% Constructor for multivariate-normal likelihood class.
  implicit none
  type            (likelihoodMultivariateNormalStochastic)                                :: likelihoodMultivariateNormalStochasticConstructor
  double precision                                        , intent(in   ), dimension(:  ) :: means
  double precision                                        , intent(in   ), dimension(:,:) :: covariance
  integer                                                 , intent(in   )                 :: realizationCount, realizationCountMinimum

  ! Store means and covariance.
  likelihoodMultivariateNormalStochasticConstructor%means                  =means
  likelihoodMultivariateNormalStochasticConstructor%covariance             =covariance
  likelihoodMultivariateNormalStochasticConstructor%realizationCount       =realizationCount
  likelihoodMultivariateNormalStochasticConstructor%realizationCountMinimum=realizationCountMinimum
  ! Find the inverse of the covariance matrix.
  likelihoodMultivariateNormalStochasticConstructor%inverseCovariance=likelihoodMultivariateNormalStochasticConstructor%covariance%invert()
  return
end function likelihoodMultivariateNormalStochasticConstructor

double precision function likelihoodMultivariateNormalStochasticEvaluate(self,simulationState,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return the log-likelihood for a multivariate-normal likelihood function.
  use Constraints_State
  use Constraints_Convergence
  use Numerical_Constants_Math
  implicit none
  class           (likelihoodMultivariateNormalStochastic), intent(inout)               :: self
  class           (state                                 ), intent(in   )               :: simulationState
  class           (convergence                           ), intent(inout)               :: simulationConvergence
  double precision                                        , intent(in   )               :: temperature          , logLikelihoodCurrent, &
       &                                                                                   logPriorCurrent      , logPriorProposed
  type            (vector                                )                              :: stateVector          , difference             , &
       &                                                                                   stateStochasticVector
  type            (matrix                                )                              :: covarianceMatrix     , covarianceMatrixInverse
  double precision                                        , allocatable, dimension(:  ) :: stateStochastic      , stateTrue
  double precision                                        , allocatable, dimension(:,:) :: covarianceStochastic , covarianceFixed        , &
       &                                                                                   covarianceFull
  integer                                                                               :: i                    , j                      , &
       &                                                                                   realizationCount
  double precision                                                                      :: temperatureEffective , likelihoodEffective    , &
       &                                                                                   logDeterminant

  ! Get the current state vector.
  stateVector=simulationState%get()
  ! Make a stochastic realization of the state vector.
  allocate(     stateTrue      (simulationState%dimension()                            ))
  allocate(     stateStochastic(simulationState%dimension()                            ))
  allocate(covarianceStochastic(simulationState%dimension(),simulationState%dimension()))
  allocate(covarianceFixed     (simulationState%dimension(),simulationState%dimension()))
  stateTrue      =stateVector
  stateStochastic=0.0d0
  ! Determine number of realizations to make.
  realizationCount    =max(int(dble(self%realizationCount)/temperature     ),self%realizationCountMinimum)
  temperatureEffective=        dble(self%realizationCount)/dble(realizationCount)
  do i=1,realizationCount
     do j=1,simulationState%dimension()
        stateStochastic(j)=stateStochastic(j)+2.0d0*stateTrue(j)*self%pseudoRandomSequence%sample()/dble(realizationCount)
     end do
  end do
  stateStochasticVector=stateStochastic
  covarianceStochastic =0.0d0
  do j=1,simulationState%dimension()
     covarianceStochastic(j,j)=stateTrue(j)**2/3.0d0/dble(realizationCount)
  end do
  covarianceFixed        =self%covariance
  covarianceFull         =covarianceFixed*temperatureEffective+covarianceStochastic
  covarianceMatrix       =covarianceFull
  covarianceMatrixInverse=covarianceMatrix%invert                ()
  logDeterminant         =covarianceMatrix%logarithmicDeterminant()
  deallocate(     stateTrue      )
  deallocate(     stateStochastic)
  deallocate(covarianceStochastic)
  deallocate(covarianceFixed     )
  ! Construct the likelihood.
  difference         =stateStochasticVector-self%means
  likelihoodEffective=-0.5d0*(difference*(covarianceMatrixInverse*difference))
  ! Correct to unit temperature.
  likelihoodMultivariateNormalStochasticEvaluate= &
       &  temperatureEffective                    &
       & *likelihoodEffective                     &
       & +(                                       &
       &    temperatureEffective                  &
       &   -1.0d0                                 &
       &  )                                       &
       & *(                                       &
       &    0.5d0                                 &
       &   *simulationState%dimension()           &
       &   *log(                                  &
       &         2.0d0                            &
       &        *Pi                               &
       &       )                                  &
       &   +0.5d0                                 &
       &   *logDeterminant                        &
       &  )
  return
end function likelihoodMultivariateNormalStochasticEvaluate

subroutine likelihoodMultivariateNormalStochasticFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodMultivariateNormalStochastic), intent(inout) :: self
  
  return
end subroutine likelihoodMultivariateNormalStochasticFunctionChanged
