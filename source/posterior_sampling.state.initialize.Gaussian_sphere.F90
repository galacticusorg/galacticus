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
     class           (randomNumberGeneratorClass), pointer                   :: randomNumberGenerator_ => null()
     double precision                                                        :: radiusSphere
     logical                                                                 :: radiusIsRelative                , usePriorMedian
     type            (varying_string            )                            :: position
     double precision                            , allocatable, dimension(:) :: stateInitial
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
    use :: String_Handling , only : String_Count_Words
    implicit none
    type            (posteriorSampleStateInitializeGaussianSphere)                            :: self
    type            (inputParameters                             ), intent(inout)             :: parameters
    class           (randomNumberGeneratorClass                  ), pointer                   :: randomNumberGenerator_
    double precision                                              , allocatable, dimension(:) :: stateInitial
    double precision                                                                          :: radiusSphere
    logical                                                                                   :: radiusIsRelative      , usePriorMedian
    type            (varying_string                              )                            :: position
    character(len=:), allocatable :: position_
    
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
    <inputParameter>
      <name>position</name>
      <description>The initial position for the sphere. If this is set to {\normalfont \ttfamily priorMedian}, then the sphere is placed at the median of the prior in each dimension. Otherwise, this must be a list of starting parameter values.</description>
      <source>parameters</source>
      <defaultValue>var_str('priorMedian')</defaultValue>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    usePriorMedian=position == 'priorMedian'
    if (.not.usePriorMedian) then
       allocate(stateInitial(String_Count_Words(position)))
       allocate(character(len=len(position)) :: position_)
       position_=position
       read (position_,*) stateInitial
    end if
    !![
    <conditionalCall>
      <call>self=posteriorSampleStateInitializeGaussianSphere(radiusSphere,radiusIsRelative,randomNumberGenerator_,usePriorMedian{conditions})</call>
      <argument name="stateInitial" value="stateInitial" condition=".not.usePriorMedian"/>
    </conditionalCall>
    !!]
    self%position=position
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function gaussianSphereConstructorParameters

  function gaussianSphereConstructorInternal(radiusSphere,radiusIsRelative,randomNumberGenerator_,usePriorMedian,stateInitial) result(self)
    !!{
    Internal constructor for the \refClass{posteriorSampleStateInitializeGaussianSphere} posterior sampling state initialization class.
    !!}
    implicit none
    type            (posteriorSampleStateInitializeGaussianSphere)                                        :: self
    double precision                                              , intent(in   )                         :: radiusSphere
    logical                                                       , intent(in   )                         :: radiusIsRelative
    class           (randomNumberGeneratorClass                  ), intent(in   ), target                 :: randomNumberGenerator_
    logical                                                       , intent(in   )                         :: usePriorMedian
    double precision                                              , intent(in   ), dimension(:), optional :: stateInitial
    !![
    <constructorAssign variables="radiusSphere, radiusIsRelative, *randomNumberGenerator_, usePriorMedian, stateInitial"/>
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
    logical                                                                                                          :: first
    !$GLC attributes unused :: modelLikelihood

    ! No knowledge of evaluation time.
    timeEvaluatePrevious=-1.0d0
    ! We have no information about the likelihood of this state.
    logLikelihood=logImpossible
    logPosterior =logImpossible
    ! If an initial state was provided, check it is of the correct size.
    if (.not.self%usePriorMedian .and. size(self%stateInitial) /= simulationState%dimension()) &
         & call Error_Report('initial state has the wrong size'//{introspection:location})
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
               &   *                            radius
          if (self%usePriorMedian) then
             state(j)=+     state        (j) &
                  &   +distributionMedian
          else
             state(j)=+     state        (j) &
                  &   +self%stateInitial (j)
          end if
       end do
    end do
    call simulationState%update(state,.false.,.false.)
  return
  end subroutine gaussianSphereInitialize
