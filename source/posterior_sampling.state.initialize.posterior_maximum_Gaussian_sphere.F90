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
     type            (varying_string                     )          :: logFileRoot
     class           (randomNumberGeneratorClass         ), pointer :: randomNumberGenerator_          => null()
     class           (posteriorSampleStateInitializeClass), pointer :: posteriorSampleStateInitialize_ => null()
     double precision                                               :: radiusSphere
     logical                                                        :: radiusIsRelative
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
    class           (posteriorSampleStateInitializeClass                         ), pointer       :: posteriorSampleStateInitialize_
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
    <objectBuilder class="randomNumberGenerator"          name="randomNumberGenerator_"          source="parameters"/>
    <objectBuilder class="posteriorSampleStateInitialize" name="posteriorSampleStateInitialize_" source="parameters"/>
    !!]
    self=posteriorSampleStateInitializePosteriorMaximumGaussianSphere(logFileRoot,radiusSphere,radiusIsRelative,posteriorSampleStateInitialize_,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function posteriorMaximumGaussianSphereConstructorParameters

  function posteriorMaximumGaussianSphereConstructorInternal(logFileRoot,radiusSphere,radiusIsRelative,posteriorSampleStateInitialize_,randomNumberGenerator_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily posteriorMaximumGaussianSphere} posterior sampling state initialization class.
    !!}
    implicit none
    type            (posteriorSampleStateInitializePosteriorMaximumGaussianSphere)                        :: self
    type            (varying_string                                              ), intent(in   )         :: logFileRoot
    double precision                                                              , intent(in   )         :: radiusSphere
    logical                                                                       , intent(in   )         :: radiusIsRelative
    class           (posteriorSampleStateInitializeClass                         ), intent(in   ), target :: posteriorSampleStateInitialize_
    class           (randomNumberGeneratorClass                                  ), intent(in   ), target :: randomNumberGenerator_
    !![
    <constructorAssign variables="logFileRoot, radiusSphere, radiusIsRelative, *posteriorSampleStateInitialize_, *randomNumberGenerator_"/>
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
    <objectDestructor name="self%randomNumberGenerator_"         />
    <objectDestructor name="self%posteriorSampleStateInitialize_"/>
    !!]
    return
  end subroutine posteriorMaximumGaussianSphereDestructor

  subroutine posteriorMaximumGaussianSphereInitialize(self,simulationState,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    !!{
    Initialize simulation state by reading parameter values from a parameter file.
    !!}
    use :: Models_Likelihoods_Constants, only : logImpossible
    use :: Posterior_Sampling_State    , only : posteriorSampleStateClass, posteriorSampleStateSimple
    use :: File_Utilities              , only : File_Exists
    implicit none
    class           (posteriorSampleStateInitializePosteriorMaximumGaussianSphere), intent(inout)               :: self
    class           (posteriorSampleStateClass                                   ), intent(inout)               :: simulationState
    class           (posteriorSampleLikelihoodClass                              ), intent(inout)               :: modelLikelihood
    type            (modelParameterList                                          ), intent(inout), dimension(:) :: modelParameters_
    double precision                                                              , intent(  out)               :: timeEvaluatePrevious , logLikelihood       , &
         &                                                                                                         logPosterior
    double precision                                                              , allocatable  , dimension(:) :: stateVector          , stateVectorPriorBest, &
         &                                                                                                         stateVectorPerturbed , stateVectorPrior
    integer                                                                       , allocatable  , dimension(:) :: stateMap
    type            (posteriorSampleStateSimple                                  ), allocatable                 :: simulationState__
    type            (varying_string                                              )                              :: logFileName
    integer                                                                                                     :: stateCount           , mpiRank             , &
         &                                                                                                         logFileUnit          , ioStatus            , &
         &                                                                                                         i                    , iChain              , &
         &                                                                                                         indexPrior           , iPass               , &
         &                                                                                                         countPrior
    double precision                                                                                            :: logPosterior_        , logLikelihood_      , &
         &                                                                                                         distributionMinimum  , distributionMaximum , &
         &                                                                                                         timeEvaluatePrevious_, radius
    logical                                                                                                     :: converged
    character       (len=   4                                                    )                              :: label
    character       (len=4096                                                    )                              :: line                 , namePrior
    logical                                                                                                     :: first
    !$GLC attributes unused :: modelLikelihood

    ! Get the state from our fallback initializor.
    allocate(simulationState__)
    simulationState__=posteriorSampleStateSimple(1)
    call simulationState__                                %parameterCountSet(simulationState  %dimension()                                                                                 )
    call self             %posteriorSampleStateInitialize_%initialize       (simulationState__            ,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    stateVector=simulationState__%get()
    deallocate(simulationState__)
    ! Assume we have no information about the likelihood of the state by default.
    logLikelihood        =logImpossible
    logPosterior         =logImpossible
    timeEvaluatePrevious_=0.0d0
    ! First read the chain file to establish the number and names of parameters in the prior simulation.
    logFileName=self%logFileRoot//'_0000.log'
    if (File_Exists(logFileName)) then
       countPrior=0
       do iPass=1,2
          open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted',ioStat=ioStatus)
          do while (ioStatus == 0)
             read (logFileUnit,'(a)',iostat=ioStatus) line
             if (      ioStatus           /=  0 ) cycle
             if (index(line         ,"=") ==  0 ) cycle
             if (      line    (1:1)      /= "#") exit
             if (      iPass              ==  1 ) countPrior=countPrior+1
             if (      iPass              ==  1 ) cycle
             ! On the second pass, determine which, if any, parameter in the current simulation this parameter corresponds to.
             !! Remove the leading "#".
             line=line(2:)
             !! Read the parameter index.
             read (line,*) indexPrior
             !! Read the parameter name.
             read (line(index(line,"`")+1:len_trim(line)-1),'(a)') namePrior
             indexPrior=indexPrior-6
             if (indexPrior > 0) then
                stateMap(indexPrior)=-1
                do i=1,size(stateVector)
                   if (trim(namePrior) == modelParameters_(i)%modelParameter_%name()) stateMap(indexPrior)=i
                end do
             end if
          end do
          close(logFileUnit)
          if (iPass == 1) then
             ! The first six lines of the chain file are for fixed output properties (not model parameters). Correct the count for
             ! that here.
             countPrior=countPrior-6
             allocate(stateMap(countPrior))
          end if
       end do
       ! Allocate the state vectors.
       allocate(stateVectorPrior    (                countPrior  ))
       allocate(stateVectorPriorBest(                countPrior  ))
       allocate(stateVectorPerturbed(simulationState%dimension ()))
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
                  &        stateVectorPrior
             if (logPosterior_ <= logPosterior) cycle
             logPosterior        =logPosterior_
             logLikelihood       =logLikelihood_
             timeEvaluatePrevious=timeEvaluatePrevious_
             stateVectorPriorBest=stateVectorPrior
          end do
          close(logFileUnit)
          iChain=iChain+1
       end do
       ! Copy state from the prior model into the new model, mapping as we go.
       do i=1,size(stateMap)
          if (stateMap(i) /= -1) stateVector(stateMap(i))=modelParameters_(stateMap(i))%modelParameter_%map(stateVectorPriorBest(i))
       end do
    end if
    ! Apply Gaussian random perturbations around the state.
    stateVectorPerturbed=stateVector
    do i=1,simulationState%dimension()
       ! Get the median, minimum, and maximum of the prior.
       distributionMinimum=modelParameters_(i)%modelParameter_%map(modelParameters_(i)%modelParameter_%priorMinimum())
       distributionMaximum=modelParameters_(i)%modelParameter_%map(modelParameters_(i)%modelParameter_%priorMaximum())
       ! Determine the radius of the Gaussian.
       radius=self%radiusSphere
       if (self%radiusIsRelative) radius=radius*(distributionMaximum-distributionMinimum)
       ! Sample values until a value within the allowed prior range is found.
       first=.true.
       do while (                                               &
            &     first                                         &
            &    .or.                                           &
            &     stateVectorPerturbed(i) < distributionMinimum &
            &    .or.                                           &
            &     stateVectorPerturbed(i) > distributionMaximum &
            &   )
          first   =.false.
          stateVectorPerturbed(i)=+self%randomNumberGenerator_%standardNormalSample( ) &
               &                  *                            radius                  &
               &                  +                            stateVector         (i)
       end do
    end do
    ! Set the simulation state.
    call simulationState%update(stateVectorPerturbed,.false.,.false.)
    ! Reset likelihoods - we do not want to assume that these will be the same as in the prior model.
    logLikelihood       =logImpossible
    logPosterior        =logImpossible
    timeEvaluatePrevious=-1.0d0
    return
  end subroutine posteriorMaximumGaussianSphereInitialize
