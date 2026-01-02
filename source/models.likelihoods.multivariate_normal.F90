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
  Implementation of a posterior sampling likelihood class which implements a multivariate normal likelihood.
  !!}

  use :: Linear_Algebra, only : matrix, vector

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodMultivariateNormal">
   <description>
    The likelihood is a simple multivariate Gaussian, intended primarily for testing purposes. The distribution parameters are
    specified within the {\normalfont \ttfamily likelihood} element using:
    \begin{verbatim}
      &lt;mean>0.45 0.50&lt;/mean>
      &lt;covariance>
        &lt;row>1.0e-4 -0.9e-4&lt;/row>
        &lt;row>-0.9e-4 1.0e-4&lt;/row>
      &lt;/covariance>
    \end{verbatim}
    where the {\normalfont \ttfamily mean} element gives the mean vector of $N$ elements, and the {\normalfont \ttfamily covariance}
    element contains $N$ {\normalfont \ttfamily row} elements each containing a vector of $N$ elements giving a single row of the
    covariance matrix. The likelihood is then:
    \begin{equation}
    \log \mathcal{L} = - {1 \over 2} \Delta \mathcal{C}^{-1} \Delta^\mathrm{T},
    \end{equation}
    where $\Delta = \theta - \bar{\theta}$, $\theta$ is the state, $\bar{\theta}$ is the mean, and $\mathcal{C}$ is the covariance
    matrix.
   </description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodMultivariateNormal
     !!{
     Implementation of a posterior sampling likelihood class which implements a multivariate likelihood.
     !!}
     private
     type(vector) :: means
     type(matrix) :: covariance
   contains
     procedure :: evaluate        => multivariateNormalEvaluate
     procedure :: functionChanged => multivariateNormalFunctionChanged
     procedure :: descriptor      => multivariateNormalDescriptor
  end type posteriorSampleLikelihoodMultivariateNormal

  interface posteriorSampleLikelihoodMultivariateNormal
     !!{
     Constructors for the \refClass{posteriorSampleLikelihoodMultivariateNormal} posterior sampling convergence class.
     !!}
     module procedure multivariateNormalConstructorParameters
     module procedure multivariateNormalConstructorInternal
  end interface posteriorSampleLikelihoodMultivariateNormal

contains

  function multivariateNormalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleLikelihoodMultivariateNormal} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleLikelihoodMultivariateNormal)                              :: self
    type            (inputParameters                            ), intent(inout)               :: parameters
    double precision                                             , allocatable, dimension(:  ) :: means
    double precision                                             , allocatable, dimension(:,:) :: covariance

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
    !!]
    self=posteriorSampleLikelihoodMultivariateNormal(means,covariance)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function multivariateNormalConstructorParameters

  function multivariateNormalConstructorInternal(means,covariance) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleLikelihoodMultivariateNormal} convergence class.
    !!}
    use :: Linear_Algebra, only : assignment(=)
    implicit none
    type            (posteriorSampleLikelihoodMultivariateNormal)                                :: self
    double precision                                             , intent(in   ), dimension(:  ) :: means
    double precision                                             , intent(in   ), dimension(:,:) :: covariance
    !![
    <constructorAssign variables="means, covariance" allocate="no"/>
    !!]

    return
  end function multivariateNormalConstructorInternal

  double precision function multivariateNormalEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for a multivariate-normal likelihood function.
    !!}
    use :: Linear_Algebra                , only : assignment(=)                  , operator(*)
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodMultivariateNormal), intent(inout), target       :: self
    class           (posteriorSampleStateClass                  ), intent(inout)               :: simulationState
    type            (modelParameterList                         ), intent(inout), dimension(:) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass            ), intent(inout)               :: simulationConvergence
    double precision                                             , intent(in   )               :: temperature           , logLikelihoodCurrent    , &
         &                                                                                        logPriorCurrent       , logPriorProposed
    real                                                         , intent(inout)               :: timeEvaluate
    double precision                                             , intent(  out), optional     :: logLikelihoodVariance
    logical                                                      , intent(inout), optional     :: forceAcceptance
    double precision                                             , allocatable  , dimension(:) :: stateArray
    integer                                                                                    :: i
    type            (vector                                     )                              :: stateVector          , difference
    !$GLC attributes unused :: timeEvaluate, temperature, simulationConvergence, logPriorProposed, logPriorCurrent, logLikelihoodCurrent, modelParametersInactive_, forceAcceptance

    ! There is no variance in our likelihood estimate.
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Evaluate the likelihood.
    allocate(stateArray(simulationState%dimension()))
    stateArray =simulationState%get()
    do i=1,size(stateArray)
       stateArray(i)=modelParametersActive_(i)%modelParameter_%unmap(stateArray(i))
    end do
    stateVector               =stateArray
    difference                =stateVector-self%means
    multivariateNormalEvaluate=-0.5d0*self%covariance%covarianceProduct(difference)
    return
  end function multivariateNormalEvaluate

  subroutine multivariateNormalFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodMultivariateNormal), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine multivariateNormalFunctionChanged

  subroutine multivariateNormalDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    use :: Linear_Algebra  , only : assignment(=)
    implicit none
    class           (posteriorSampleLikelihoodMultivariateNormal), intent(inout)               :: self
    type            (inputParameters                            ), intent(inout)               :: descriptor
    logical                                                      , intent(in   ), optional     :: includeClass  , includeFileModificationTimes
    character       (len=18                                     )                              :: parameterLabel
    type            (inputParameters                            )                              :: parameters
    double precision                                             , allocatable, dimension(:  ) :: means
    double precision                                             , allocatable, dimension(:,:) :: covariance
    type            (varying_string                             )                              :: parameterCombined
    integer                                                                                    :: i                , j

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('posteriorSampleLikelihood','multivariateNormal')
    parameters       =descriptor%subparameters('posteriorSampleLikelihood')
    means            =self%means
    covariance       =self%covariance
    parameterCombined=""
    do i=1,size(means,dim=1)
       write (parameterLabel,'(e17.10)') means(i)
       if (i > 1) parameterCombined=parameterCombined//" "
       parameterCombined=parameterCombined//trim(adjustl(parameterLabel))
    end do
    call parameters%addParameter('means',trim(adjustl(parameterLabel)))
    do i=1,size(covariance,dim=1)
       do j=1,size(covariance,dim=2)
          write (parameterLabel,'(e17.10)') covariance(i,j)
          if (i > 1 .or. j > 1) parameterCombined=parameterCombined//" "
          parameterCombined=parameterCombined//trim(adjustl(parameterLabel))
       end do
    end do
    call parameters%addParameter('covariance',trim(adjustl(parameterLabel)))
    return
  end subroutine multivariateNormalDescriptor
