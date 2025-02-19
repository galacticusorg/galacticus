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
  Implementation of a posterior sampling likelihood class which implements a multivariate normal likelihood.
  !!}

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass
  
  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodMltiVrtNormalStochastic">
   <description>
    The likelihood is identical to that of the {\normalfont \ttfamily multivariateNormal} class, except that the likelihood function
    is evaluated stochastically. In addition to the parameter of the {\normalfont \ttfamily multivariateNormal} class, two additional
    parameters are required and are specified within the {\normalfont \ttfamily likelihood} element using:
    \begin{verbatim}
      &lt;realizationCount>4000&lt;/realizationCount>
      &lt;realizationCountMinimum>10&lt;/realizationCountMinimum>
    \end{verbatim}
    When evaluating the likelihood, the state vector is set equal to 
    \begin{equation}
     S^\prime_i = \sum_{j=1}^N {2 U(S_i) \over N},
    \end{equation}
    where $N=${\normalfont \ttfamily realizationCount} and $U(x)$ is a uniform random deviate in the range $0$ to $x$. This results in
    a variance in $S^\prime_i$ of $S_i^2/3N$. This variance is added to the covariance used in evaluating the likelihood. When
    evaluating the likelihood at a higher temperature the number of realizations is reduced (which increases the covariance, which has
    the same effect as increasing the temperature) to speed computation, and the likelihood corrected for this fact. The number of
    realizations is reduced to $N/T$, but never allowed to fall below {\normalfont \ttfamily realizationCountMinimum}.
   </description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodMultivariateNormal) :: posteriorSampleLikelihoodMltiVrtNormalStochastic
     !!{
     Implementation of a posterior sampling likelihood class which implements a multivariate likelihood.
     !!}
     private
     integer                                      :: realizationCount                , realizationCountMinimum
     class  (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
   contains
     final     ::                    multivariateNormalStochasticDestructor
     procedure :: evaluate        => multivariateNormalStochasticEvaluate
     procedure :: functionChanged => multivariateNormalStochasticFunctionChanged
  end type posteriorSampleLikelihoodMltiVrtNormalStochastic

  interface posteriorSampleLikelihoodMltiVrtNormalStochastic
     !!{
     Constructors for the {\normalfont \ttfamily multivariateNormalStochastic} posterior sampling convergence class.
     !!}
     module procedure multivariateNormalStochasticConstructorParameters
     module procedure multivariateNormalStochasticConstructorInternal
  end interface posteriorSampleLikelihoodMltiVrtNormalStochastic

contains

  function multivariateNormalStochasticConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily multivariateNormalStochastic} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodMltiVrtNormalStochastic)                              :: self
    type            (inputParameters                                 ), intent(inout)               :: parameters
    class           (randomNumberGeneratorClass                      ), pointer                     :: randomNumberGenerator_
    double precision                                                  , allocatable, dimension(:  ) :: means
    double precision                                                  , allocatable, dimension(:,:) :: covariance
    integer                                                                                         :: realizationCount, realizationCountMinimum

    allocate(means     (parameters%count('means')                          ))
    allocate(covariance(parameters%count('means'),parameters%count('means')))
    !![
    <inputParameter>
      <name>means</name>
      <description>The mean of the multivariate normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>covariance</name>
      <description>The covariance matrix for the of the multivariate normal distribution.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>realizationCount</name>
      <description>The number of realizations of the stochastic likelihood to compute at unit temperature.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>realizationCountMinimum</name>
      <description>The minimum number of realizations of the stochastic likelihood to compute at higher temperatures.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=posteriorSampleLikelihoodMltiVrtNormalStochastic(means,covariance,realizationCount,realizationCountMinimum,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function multivariateNormalStochasticConstructorParameters

  function multivariateNormalStochasticConstructorInternal(means,covariance,realizationCount,realizationCountMinimum,randomNumberGenerator_) result(self)
    !!{
    Constructor for ``multivariateNormalStochastic'' convergence class.
    !!}
    type            (posteriorSampleLikelihoodMltiVrtNormalStochastic)                                :: self
    double precision                                                  , intent(in   ), dimension(:  ) :: means
    double precision                                                  , intent(in   ), dimension(:,:) :: covariance
    integer                                                           , intent(in   )                 :: realizationCount, realizationCountMinimum
    class           (randomNumberGeneratorClass                      ), intent(in   ), target         :: randomNumberGenerator_
    !![
    <constructorAssign variables="realizationCount, realizationCountMinimum, *randomNumberGenerator_"/>
    !!]

    self%posteriorSampleLikelihoodMultivariateNormal=posteriorSampleLikelihoodMultivariateNormal(means,covariance)
    return
  end function multivariateNormalStochasticConstructorInternal

  subroutine multivariateNormalStochasticDestructor(self)
    !!{
    Destructor for the  {\normalfont \ttfamily multivariateNormalStochastic} model likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodMltiVrtNormalStochastic), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine multivariateNormalStochasticDestructor

  double precision function multivariateNormalStochasticEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for a multivariate-normal likelihood function.
    !!}
    use :: Linear_Algebra                , only : assignment(=)                  , operator(*)
    use :: Numerical_Constants_Math      , only : Pi
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodMltiVrtNormalStochastic), intent(inout), target         :: self
    class           (posteriorSampleStateClass                       ), intent(inout)                 :: simulationState
    type            (modelParameterList                              ), intent(inout), dimension(:)   :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass                 ), intent(inout)                 :: simulationConvergence
    double precision                                                  , intent(in   )                 :: temperature           , logLikelihoodCurrent    , &
         &                                                                                               logPriorCurrent       , logPriorProposed
    real                                                              , intent(inout)                 :: timeEvaluate
    double precision                                                  , intent(  out), optional       :: logLikelihoodVariance
    logical                                                           , intent(inout), optional       :: forceAcceptance
    type            (vector                                          )                                :: stateVector           , difference              , &
         &                                                                                               stateStochasticVector
    type            (matrix                                          )                                :: covarianceMatrix
    double precision                                                  , allocatable  , dimension(:  ) :: stateStochastic       , stateTrue
    double precision                                                  , allocatable  , dimension(:,:) :: covarianceStochastic  , covarianceFixed         , &
         &                                                                                               covarianceFull
    integer                                                                                           :: i                     , j                       , &
         &                                                                                               realizationCount
    double precision                                                                                  :: temperatureEffective  , likelihoodEffective     , &
         &                                                                                               logDeterminant
    !$GLC attributes unused :: self, logLikelihoodCurrent, logPriorCurrent, logPriorProposed, simulationConvergence, timeEvaluate, modelParametersInactive_, forceAcceptance

    ! Report zero variance, as this class is designed specifically for testing for unaccounted-for randomness in the likelihood
    ! estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Get the current state vector.
    stateVector=simulationState%get()
    ! Make a stochastic realization of the state vector.
    allocate(     stateTrue      (simulationState%dimension()                            ))
    allocate(     stateStochastic(simulationState%dimension()                            ))
    allocate(covarianceStochastic(simulationState%dimension(),simulationState%dimension()))
    allocate(covarianceFixed     (simulationState%dimension(),simulationState%dimension()))
    stateTrue      =stateVector
    stateStochastic=0.0d0
    ! Unmap the state vector.
    do i=1,size(stateTrue)
       stateTrue(i)=modelParametersActive_(i)%modelParameter_%unmap(stateTrue(i))
    end do
    ! Determine number of realizations to make.
    realizationCount    =max(int(dble(self%realizationCount)/     temperature     ),self%realizationCountMinimum)
    temperatureEffective=        dble(self%realizationCount)/dble(realizationCount)
    do i=1,realizationCount
       do j=1,simulationState%dimension()
          stateStochastic(j)=stateStochastic(j)+2.0d0*stateTrue(j)*self%randomNumberGenerator_%uniformSample()/dble(realizationCount)
       end do
    end do
    stateStochasticVector=stateStochastic
    covarianceStochastic =0.0d0
    do j=1,simulationState%dimension()
       covarianceStochastic(j,j)=stateTrue(j)**2/3.0d0/dble(realizationCount)
    end do
    covarianceFixed =self%covariance
    covarianceFull  =covarianceFixed*temperatureEffective+covarianceStochastic
    covarianceMatrix=covarianceFull
    logDeterminant  =covarianceMatrix%logarithmicDeterminant()
    deallocate(     stateTrue      )
    deallocate(     stateStochastic)
    deallocate(covarianceStochastic)
    deallocate(covarianceFixed     )
    ! Construct the likelihood.
    difference         =stateStochasticVector-self%means
    likelihoodEffective=-0.5d0*covarianceMatrix%covarianceProduct(difference)
    ! Correct to unit temperature.
    multivariateNormalStochasticEvaluate= &
         &  temperatureEffective          &
         & *likelihoodEffective           &
         & +(                             &
         &    temperatureEffective        &
         &   -1.0d0                       &
         &  )                             &
         & *(                             &
         &    0.5d0                       &
         &   *simulationState%dimension() &
         &   *log(                        &
         &         2.0d0                  &
         &        *Pi                     &
         &       )                        &
         &   +0.5d0                       &
         &   *logDeterminant              &
         &  )
    return
  end function multivariateNormalStochasticEvaluate

  subroutine multivariateNormalStochasticFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodMltiVrtNormalStochastic), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine multivariateNormalStochasticFunctionChanged
