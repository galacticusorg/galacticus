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
  Implementation of a posterior sampling state initializer class which initializes chains to the the maximum posterior point from
  a set of chain files, with a Gaussian sphere scatter around that point.
  !!}

  use :: ISO_Varying_String, only : varying_string

  !![
  <posteriorSampleStateInitialize name="posteriorSampleStateInitializePosteriorMaximumGaussianSphere">
    <description>
      Initializes chains to the the maximum posterior point from a set of chain files, with a Gaussian sphere scatter around that
      point.
    </description>
  </posteriorSampleStateInitialize>
  !!]
  type, extends(posteriorSampleStateInitializeClass) :: posteriorSampleStateInitializePosteriorMaximumGaussianSphere
     !!{
     Implementation of a posterior sampling state initialization class that initializes state to the maximum posterior point from
     a set of chain files, with a Gaussian sphere scatter around that point.
     !!}
     private
     type            (varying_string            )          :: logFileRoot
     class           (randomNumberGeneratorClass), pointer :: randomNumberGenerator_ => null()
     double precision                                      :: radiusSphere
     logical                                               :: radiusIsRelative
   contains
     final     ::                posteriorMaximumGaussianSphereDestructor
     procedure :: initialize  => posteriorMaximumGaussianSphereInitialize
  end type posteriorSampleStateInitializePosteriorMaximumGaussianSphere

  interface posteriorSampleStateInitializePosteriorMaximumGaussianSphere
     !!{
     Constructors for the {\normalfont \ttfamily posteriorMaximumGaussianSphere} posterior sampling state initialization class.
     !!}
     module procedure posteriorMaximumGaussianSphereConstructorParameters
     module procedure posteriorMaximumGaussianSphereConstructorInternal
  end interface posteriorSampleStateInitializePosteriorMaximumGaussianSphere

contains

  function posteriorMaximumGaussianSphereConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily posteriorMaximumGaussianSphere} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleStateInitializePosteriorMaximumGaussianSphere)                :: self
    type            (inputParameters                                             ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass                                  ), pointer       :: randomNumberGenerator_
    double precision                                                                              :: radiusSphere
    logical                                                                                       :: radiusIsRelative
    type            (varying_string                                              )                :: logFileRoot

    !![
    <inputParameter>
      <name>logFileRoot</name>
      <description>The root file name of the state files from which to resume.</description>
      <source>parameters</source>
    </inputParameter>
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
    self=posteriorSampleStateInitializePosteriorMaximumGaussianSphere(logFileRoot,radiusSphere,radiusIsRelative,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function posteriorMaximumGaussianSphereConstructorParameters

  function posteriorMaximumGaussianSphereConstructorInternal(logFileRoot,radiusSphere,radiusIsRelative,randomNumberGenerator_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily posteriorMaximumGaussianSphere} posterior sampling state initialization class.
    !!}
    implicit none
    type            (posteriorSampleStateInitializePosteriorMaximumGaussianSphere)                        :: self
    type            (varying_string                                              ), intent(in   )         :: logFileRoot
    double precision                                                              , intent(in   )         :: radiusSphere
    logical                                                                       , intent(in   )         :: radiusIsRelative
    class           (randomNumberGeneratorClass                                  ), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="logFileRoot, radiusSphere, radiusIsRelative, *randomNumberGenerator_"/>
    !!]

    return
  end function posteriorMaximumGaussianSphereConstructorInternal

  subroutine posteriorMaximumGaussianSphereDestructor(self)
    !!{
    Destructor for the  {\normalfont \ttfamily posteriorMaximumGaussianSphere} posterior sampling state initialization class.
    !!}
    implicit none
    type(posteriorSampleStateInitializePosteriorMaximumGaussianSphere), intent(inout) :: self

    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine posteriorMaximumGaussianSphereDestructor

  subroutine posteriorMaximumGaussianSphereInitialize(self,simulationState,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    !!{
    Initialize simulation state by reading parameter values from a parameter file.
    !!}
    use :: Models_Likelihoods_Constants, only : logImpossible
    use :: Posterior_Sampling_State    , only : posteriorSampleStateClass
    use :: File_Utilities              , only : File_Exists
    implicit none
    class           (posteriorSampleStateInitializePosteriorMaximumGaussianSphere), intent(inout)               :: self
    class           (posteriorSampleStateClass                                   ), intent(inout)               :: simulationState
    class           (posteriorSampleLikelihoodClass                              ), intent(inout)               :: modelLikelihood
    type            (modelParameterList                                          ), intent(inout), dimension(:) :: modelParameters_
    double precision                                                              , intent(  out)               :: timeEvaluatePrevious      , logLikelihood      , &
         &                                                                                                         logPosterior
    double precision                                                              , allocatable  , dimension(:) :: stateVector               , stateVectorMapped  , &
         &                                                                                                         stateVectorMappedPerturbed, stateVector_
    type            (varying_string                                              )                              :: logFileName
    integer                                                                                                     :: stateCount                , mpiRank            , &
         &                                                                                                         logFileUnit               , ioStatus           , &
         &                                                                                                         i                         , iChain
    double precision                                                                                            :: logPosterior_             , logLikelihood_     , &
         &                                                                                                         distributionMinimum       , distributionMaximum, &
         &                                                                                                         timeEvaluatePrevious_     , radius
    logical                                                                                                     :: converged
    character       (len=   4                                                    )                              :: label
    character       (len=4096                                                    )                              :: line
    logical                                                                                                     :: first
    !$GLC attributes unused :: modelLikelihood

    ! Assume we have no information about the likelihood of the state by default.
    logLikelihood        =logImpossible
    logPosterior         =logImpossible
    timeEvaluatePrevious_=0.0d0
    ! Allocate the state vector.
    allocate(stateVector_              (simulationState%dimension()))
    allocate(stateVector               (simulationState%dimension()))
    allocate(stateVectorMapped         (simulationState%dimension()))
    allocate(stateVectorMappedPerturbed(simulationState%dimension()))
    ! Read state from the log files.
    iChain=0
    do while (iChain >= 0)
       write (label,'(i4.4)') iChain
       logFileName=self%logFileRoot//'_'//trim(label)//'.log'
       if (.not.File_Exists(logFileName)) then
          iChain=-1
          exit
       end if
       open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
       ioStatus=0
       do while (ioStatus == 0)
          read (logFileUnit,'(a)',iostat=ioStatus) line
          if (ioStatus  /=  0 ) cycle
          if (line(1:1) == "#") cycle
          read (line,*) stateCount           , &
               &        mpiRank              , &
               &        timeEvaluatePrevious_, &
               &        converged            , &
               &        logPosterior_        , &
               &        logLikelihood_       , &
               &        stateVector_
          if (logPosterior_ <= logPosterior) cycle
          logPosterior        =logPosterior_
          logLikelihood       =logLikelihood_
          timeEvaluatePrevious=timeEvaluatePrevious_
          stateVector         =stateVector_
       end do
       close(logFileUnit)
       iChain=iChain+1
    end do
    ! Map the state.
    do i=1,size(stateVector)
       stateVectorMapped(i)=modelParameters_(i)%modelParameter_%map(stateVector(i))
    end do
    ! Apply Gaussian random perturbations around the state.
    stateVectorMappedPerturbed=stateVectorMapped
    do i=1,simulationState%dimension()
       ! Get the median, minimum, and maximum of the prior.
       distributionMinimum=modelParameters_(i)%modelParameter_%map(modelParameters_(i)%modelParameter_%priorMinimum())
       distributionMaximum=modelParameters_(i)%modelParameter_%map(modelParameters_(i)%modelParameter_%priorMaximum())
       ! Determine the radius of the Gaussian.
       radius=self%radiusSphere
       if (self%radiusIsRelative) radius=radius*(distributionMaximum-distributionMinimum)
       ! Sample values until a value within the allowed prior range is found.
       first=.true.
       do while (                                                     &
            &     first                                               &
            &    .or.                                                 &
            &     stateVectorMappedPerturbed(i) < distributionMinimum &
            &    .or.                                                 &
            &     stateVectorMappedPerturbed(i) > distributionMaximum &
            &   )
          first   =.false.
          stateVectorMappedPerturbed(i)=+self%randomNumberGenerator_%standardNormalSample( ) &
               &                        *                            radius                  &
               &                        +                            stateVectorMapped   (i)
       end do
    end do
    ! Set the simulation state.
    call simulationState%update(stateVectorMappedPerturbed,.false.,.false.)
    ! Reset likelihoods - we do not want to assume that these will be the same as in the prior model.
    logLikelihood       =logImpossible
    logPosterior        =logImpossible
    timeEvaluatePrevious=-1.0d0
    return
  end subroutine posteriorMaximumGaussianSphereInitialize
