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
  Implementation of a posterior sampling state initializor class which initializes states using random draws from a Gaussian
  sphere centered on the mode of the prior distribution.  
  !!}

  !![
  <posteriorSampleStateInitialize name="posteriorSampleStateInitializeGaussianSphere">
    <description>
      A posterior sampling state initialization class which initializes states using random draws from a Gaussian sphere centered on
      the mode of the prior distribution.
   </description>
  </posteriorSampleStateInitialize>
  !!]
  type, extends(posteriorSampleStateInitializeClass) :: posteriorSampleStateInitializeGaussianSphere
     !!{
     Implementation of a posterior sampling state initialization class which initializes states using random draws from a Gaussian
     sphere centered on the mode of the prior distribution.
     !!}
     private
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     double precision                                      :: radiusSphere
     logical                                               :: radiusIsRelative
   contains
     final     ::                gaussianSphereDestructor
     procedure :: initialize  => gaussianSphereInitialize
  end type posteriorSampleStateInitializeGaussianSphere

  interface posteriorSampleStateInitializeGaussianSphere
     !!{
     Constructors for the \refClass{posteriorSampleStateInitializeGaussianSphere} posterior sampling state initialization class.
     !!}
     module procedure gaussianSphereConstructorParameters
     module procedure gaussianSphereConstructorInternal
  end interface posteriorSampleStateInitializeGaussianSphere

contains

  function gaussianSphereConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeGaussianSphere} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (posteriorSampleStateInitializeGaussianSphere)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass                  ), pointer       :: randomNumberGenerator_
    double precision                                                              :: radiusSphere
    logical                                                                       :: radiusIsRelative

    !![
    <inputParameter>
      <name>radiusSphere</name>
      <description>The radius of the Gaussian sphere.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusIsRelative</name>
      <description>If true, the radius of the sphere is assumed to be relative to the extent of the prior, otherwise it is assumed to be an absolute radius.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=posteriorSampleStateInitializeGaussianSphere(radiusSphere,radiusIsRelative,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function gaussianSphereConstructorParameters

  function gaussianSphereConstructorInternal(radiusSphere,radiusIsRelative,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{posteriorSampleStateInitializeGaussianSphere} posterior sampling state initialization class.
    !!}
    implicit none
    type            (posteriorSampleStateInitializeGaussianSphere)                        :: self
    double precision                                              , intent(in   )         :: radiusSphere
    logical                                                       , intent(in   )         :: radiusIsRelative
    class           (randomNumberGeneratorClass                  ), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="radiusSphere, radiusIsRelative, *randomNumberGenerator_"/>
    !!]
    
    return
  end function gaussianSphereConstructorInternal

  subroutine gaussianSphereDestructor(self)
    !!{
    Destructor for the \refClass{posteriorSampleStateInitializeGaussianSphere} posterior sampling state initialization class.
    !!}
    implicit none
    type(posteriorSampleStateInitializeGaussianSphere), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine gaussianSphereDestructor

  subroutine gaussianSphereInitialize(self,simulationState,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    !!{
    Initialize simulation state by drawing at random from the parameter priors.
    !!}
    use :: Models_Likelihoods_Constants, only : logImpossible
    implicit none
    class           (posteriorSampleStateInitializeGaussianSphere), intent(inout)                                    :: self
    class           (posteriorSampleStateClass                   ), intent(inout)                                    :: simulationState
    class           (posteriorSampleLikelihoodClass              ), intent(inout)                                    :: modelLikelihood
    type            (modelParameterList                          ), intent(inout), dimension(:                     ) :: modelParameters_
    double precision                                              , intent(  out)                                    :: timeEvaluatePrevious, logLikelihood       , &
         &                                                                                                              logPosterior
    double precision                                                             , dimension(size(modelParameters_)) :: state
    double precision                                                                                                 :: distributionMinimum  , distributionMaximum, &
         &                                                                                                              distributionMedian   , radius
    integer                                                                                                          :: j
    logical :: first
    !$GLC attributes unused ::  self, modelLikelihood

    ! No knowledge of evaluation time.
    timeEvaluatePrevious=-1.0d0
    ! We have no information about the likelihood of this state.
    logLikelihood=logImpossible
    logPosterior =logImpossible
    ! Initialize chain to some state vector.
    state=0.0d0
    do j=1,simulationState%dimension()
       ! Get the median, minimum, and maximum of the prior.
       distributionMinimum=modelParameters_(j)%modelParameter_%map(modelParameters_(j)%modelParameter_%priorMinimum(     ))
       distributionMaximum=modelParameters_(j)%modelParameter_%map(modelParameters_(j)%modelParameter_%priorMaximum(     ))
       distributionMedian =modelParameters_(j)%modelParameter_%map(modelParameters_(j)%modelParameter_%priorInvert (0.5d0))
       ! Determine the radius of the Gaussian.
       radius=self%radiusSphere
       if (self%radiusIsRelative) radius=radius*(distributionMaximum-distributionMinimum)
       ! Sample values until a value within the allowed prior range is found.
       first=.true.
       do while (                                &
            &     first                          &
            &    .or.                            &
            &     state(j) < distributionMinimum &
            &    .or.                            &
            &     state(j) > distributionMaximum &
            &   )
          first   =.false.
          state(j)=+self%randomNumberGenerator_%standardNormalSample() &
               &   *                            radius                 &
               &   +                            distributionMedian
       end do
    end do
    call simulationState%update(state,.false.,.false.)
  return
  end subroutine gaussianSphereInitialize
