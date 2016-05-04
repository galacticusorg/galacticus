!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines methods for the ``differentialEvolution'' simulator class.

function simulatorDifferentialEvolutionConstructor(parameterPriors,randomDistributions,parameterMappings,modelLikelihood,simulationConvergence,simulationStoppingCriterion&
     &,simulationState,simulationStateInitializor,proposalSize,randomJump,stepsMaximum,acceptanceAverageCount,stateSwapCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot)
  !% Constructor for ``differentialEvolution'' simulator class.
  implicit none
  type     (simulatorDifferentialEvolution)                                      :: simulatorDifferentialEvolutionConstructor
  type     (prior                         ), intent(in   ), target, dimension(:) :: parameterPriors
  type     (distributionList              ), intent(in   ), target, dimension(:) :: randomDistributions
  type     (mappingList                   ), intent(in   ), target, dimension(:) :: parameterMappings
  class    (likelihood                    ), intent(in   ), target               :: modelLikelihood
  class    (convergence                   ), intent(in   ), target               :: simulationConvergence
  class    (stoppingCriterion             ), intent(in   ), target               :: simulationStoppingCriterion
  class    (state                         ), intent(in   ), target               :: simulationState
  class    (stateInitializor              ), intent(in   ), target               :: simulationStateInitializor
  class    (deProposalSize                ), intent(in   ), target               :: proposalSize
  class    (deRandomJump                  ), intent(in   ), target               :: randomJump
  integer                                  , intent(in   )                       :: stepsMaximum,acceptanceAverageCount,stateSwapCount,logFlushCount, reportCount
  character(len=*                         ), intent(in   )                       :: logFileRoot, interactionRoot
  logical                                  , intent(in   )                       :: sampleOutliers

  simulatorDifferentialEvolutionConstructor%parameterCount              =  size(parameterPriors)
  simulatorDifferentialEvolutionConstructor%parameterPriors             => parameterPriors
  simulatorDifferentialEvolutionConstructor%randomDistributions         => randomDistributions
  simulatorDifferentialEvolutionConstructor%parameterMappings           => parameterMappings
  simulatorDifferentialEvolutionConstructor%modelLikelihood             => modelLikelihood
  simulatorDifferentialEvolutionConstructor%simulationConvergence       => simulationConvergence
  simulatorDifferentialEvolutionConstructor%simulationStoppingCriterion => simulationStoppingCriterion
  simulatorDifferentialEvolutionConstructor%simulationState             => simulationState
  simulatorDifferentialEvolutionConstructor%simulationStateInitializor  => simulationStateInitializor
  simulatorDifferentialEvolutionConstructor%proposalSize                => proposalSize
  simulatorDifferentialEvolutionConstructor%randomJump                  => randomJump
  simulatorDifferentialEvolutionConstructor%stepsMaximum                =  stepsMaximum
  simulatorDifferentialEvolutionConstructor%acceptanceAverageCount      =  acceptanceAverageCount
  simulatorDifferentialEvolutionConstructor%stateSwapCount              =  stateSwapCount
  simulatorDifferentialEvolutionConstructor%logFileRoot                 =  logFileRoot
  simulatorDifferentialEvolutionConstructor%isInteractive               =  trim(interactionRoot) /= "none"
  simulatorDifferentialEvolutionConstructor%interactionRoot             =  trim(interactionRoot)
  simulatorDifferentialEvolutionConstructor%sampleOutliers              =  sampleOutliers
  simulatorDifferentialEvolutionConstructor%logFlushCount               =  logFlushCount
  simulatorDifferentialEvolutionConstructor%reportCount                 =  reportCount
  return
end function simulatorDifferentialEvolutionConstructor

subroutine simulatorDifferentialEvolutionSimulate(self)
  !% Perform a differential evolution simulation.
  use MPI_Utilities
  use Pseudo_Random
  use Galacticus_Error
  use Galacticus_Display
  use String_Handling
  use Constraints_Constants
  use Kind_Numbers
  use File_Utilities
  use System_Command
  implicit none
  class           (simulatorDifferentialEvolution), intent(inout)                    :: self
  integer                                         , dimension(                    2) :: chainPair
  double precision                                , dimension(self%parameterCount,2) :: statePair
  double precision                                , dimension(self%parameterCount  ) :: stateVector           , stateVectorProposed          , &
       &                                                                                stateVectorInteractive
  class           (state                         ), allocatable                      :: stateProposed
  real                                                                               :: timePreEvaluate       , timePostEvaluate             , &
       &                                                                                timeEvaluate          , timeEvaluatePrevious
  double precision                                                                   :: logLikelihoodProposed , logPosteriorProposed         , &
       &                                                                                logLikelihood         , logLikelihoodVariance        , &
       &                                                                                timeEvaluateInitial   , logLikelihoodVarianceProposed
  type            (pseudoRandom                  )                                   :: randomNumberGenerator
  type            (varying_string                )                                   :: logFileName           , message                      , &
       &                                                                                interactionFileName
  integer                                                                            :: logFileUnit           , convergedAtStep              , &
       &                                                                                convergenceFileUnit   , i                            , &
       &                                                                                ioStatus              , interactionFile
  logical                                                                               forceAcceptance
  character       (len=32                        )                                   :: label
  
  ! Check that we have sufficient chains for differential evolution.
  if (mpiSelf%count() < 3) call Galacticus_Error_Report('simulatorDifferentialEvolutionSimulate','at least 3 chains are required for differential evolution')
  ! Write start-up message.
  message="Process "//mpiSelf%rankLabel()//" [PID: "
  message=message//getPID()//"] is running on host '"//mpiSelf%hostAffinity()//"'"
  call Galacticus_Display_Message(message)
  ! Allocate a simple state object for the proposed state.
  allocate(stateSimple :: stateProposed)
  select type (stateProposed)
  type is (stateSimple)
     stateProposed=stateSimple(self%parameterCount,1)
  end select
  ! Initialize chain to some state vector.
  call self%simulationStateInitializor%initialize(self%simulationState,self%parameterPriors,self%parametermappings,timeEvaluateInitial)
  ! Evaluate the posterior in the initial state.
  timeEvaluate        =-1.0
  timeEvaluatePrevious=real(timeEvaluateInitial)
  call CPU_Time(timePreEvaluate )
  call self%posterior(self%simulationState,logImpossible,logImpossible,self%logPosterior,logLikelihood,logLikelihoodVariance,timeEvaluate,timeEvaluatePrevious)
  call CPU_Time(timePostEvaluate)
  if (timeEvaluate < 0.0) timeEvaluate=timePostEvaluate-timePreEvaluate
  timeEvaluatePrevious=timeEvaluate     
  ! Check for impossible state.
  if (self%logPosterior <= logImpossible) call Galacticus_Error_Report('simulatorDifferentialEvolutionSimulate','impossible initial state')
  ! Begin stepping.
  logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.log'
  open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
  self%isConverged=.false. 
  do while (                                                                                                                  &
       &          self%simulationState            %count(                                               ) < self%stepsMaximum &
       &    .and.                                                                                                             &
       &     .not.self%simulationStoppingCriterion%stop (self%simulationState,self%simulationConvergence)                     &
       &   )
     ! Pick two random processes to use for proposal.
     chainPair(1)=simulatorDifferentialChainSelect(self,randomNumberGenerator               )
     chainPair(2)=simulatorDifferentialChainSelect(self,randomNumberGenerator,[chainPair(1)])
     ! Receive states from selected chains.
     stateVector=self%simulationState%get()
     statePair  =mpiSelf%requestData(chainPair,stateVector)
     ! Generate proposal.
     stateVectorProposed= stateVector      &
          &              +self%stepSize()  &
          &              *(                &
          &                +statePair(:,1) &
          &                -statePair(:,2) &
          &               )
     ! Add random perturbations to the proposal.
     stateVectorProposed=stateVectorProposed+self%randomJump%sample(self%randomDistributions,self%simulationState)
     ! If simulation is interactive, check for any interaction file.
     forceAcceptance=.false.
     if (self%isInteractive) then
        ! Check if an interaction file exists.
        interactionFileName=self%interactionRoot//"_"//mpiSelf%rankLabel()
        if (File_Exists(interactionFileName)) then
           ! Read the file and validate.
           open(newUnit=interactionFile,file=char(interactionFileName),status='old',form='formatted',ioStat=ioStatus)
           if (ioStatus == 0) then
              read (interactionFile,*,ioStat=ioStatus) stateVectorInteractive
              if (ioStatus == 0) then
                 ! Copy the state to the proposed state vector.
                 stateVectorProposed=stateVectorInteractive
                 message="Chain "//mpiSelf%rankLabel()//" is being interactively moved to state:"
                 call Galacticus_Display_Indent(message)
                 ! Map parameters of interactively proposed state.
                 do i=1,size(stateVector)
                    write (label,*) stateVectorProposed(i)
                    message="State["
                    message=message//i//"] = "//trim(adjustl(label))
                    call Galacticus_Display_Message(message)
                    stateVectorProposed(i)=self%parameterMappings(i)%thisMapping%map(stateVectorProposed(i))
                 end do
                 call Galacticus_Display_Unindent('end')
                 ! Force acceptance of this state.
                 forceAcceptance=.true.
              else
                 message="WARNING: state proposed in interaction file '"//interactionFileName//"' cannot be read"
                 call Galacticus_Display_Message(message)
              end if
           else
              message="WARNING: unable to open interaction file '"//interactionFileName//"'"
              call Galacticus_Display_Message(message)
           end if
           close(interactionFile)
           ! Remove the interaction file.
           call System_Command_Do("rm -f "//interactionFileName)
        end if
     end if
     ! Store the proposed state vector.
     call stateProposed%update(stateVectorProposed,.false.,.false.)
     ! Evaluate likelihood.
     timeEvaluatePrevious=timeEvaluate
     timeEvaluate        =-1.0
     call CPU_Time(timePreEvaluate )
     call self%posterior(stateProposed,logLikelihood,self%logPosterior-logLikelihood,logPosteriorProposed,logLikelihoodProposed,logLikelihoodVarianceProposed,timeEvaluate,timeEvaluatePrevious)
     call CPU_Time(timePostEvaluate)
     if (timeEvaluate < 0.0) timeEvaluate=timePostEvaluate-timePreEvaluate
     ! Decide whether to take step.
     if     (                                                                                                                                       &
          &   self%acceptProposal(self%logPosterior,logPosteriorProposed,logLikelihoodVariance,logLikelihoodVarianceProposed,randomNumberGenerator) &
          &  .or.                                                                                                                                   &
          &   forceAcceptance                                                                                                                       &
          & ) then
        self%logPosterior    =logPosteriorProposed
        logLikelihood        =logLikelihoodProposed
        logLikelihoodVariance=logLikelihoodVarianceProposed
        stateVector          =stateVectorProposed
     else
        ! Step not accepted - retain the old estimate of work time.
        timeEvaluate         =timeEvaluatePrevious
     end if
     call self%update(stateVector)
     ! Unmap parameters and write to log file.
     do i=1,size(stateVector)
        stateVector(i)=self%parameterMappings(i)%thisMapping%unmap(stateVector(i))
     end do
     if (self%logging()) then
        write (logFileUnit,*) self   %simulationState%count(), &
             &                mpiSelf%rank                 (), &
             &                timeEvaluate                   , &
             &                self%isConverged               , &
             &                self%logPosterior              , &
             &                logLikelihood                  , &
             &                stateVector
        if (mod(self%simulationState%count(),self%logFlushCount) == 0) call flush(logFileUnit)
     end if
     ! Repeat.
     call mpiBarrier()
     ! Test for convergence.
     if (.not.self%isConverged) then
        self%isConverged=self%simulationConvergence%isConverged(self%simulationState,self%logPosterior)
        if (self%isConverged) then
           convergedAtStep=self%simulationState%count()
           if (mpiSelf%rank() == 0) then
              message='Converged after '
              message=message//convergedAtStep//' steps'
              call Galacticus_Display_Message(message)
              logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.convergence.log'
              open(newunit=convergenceFileUnit,file=char(logFileName),status='unknown',form='formatted',access='append')
              write (convergenceFileUnit,'(a,i8)') 'Converged at step: ',convergedAtStep
              call self%simulationConvergence%logReport(convergenceFileUnit)
              close(convergenceFileUnit)
           end if
        end if
     end if
  end do
  close(logFileUnit)
  return
end subroutine simulatorDifferentialEvolutionSimulate

subroutine simulatorDifferentialEvolutionUpdate(self,stateVector)
  !% Update the differential evolution simulator state.
  use MPI_Utilities
  implicit none
  class           (simulatorDifferentialEvolution), intent(inout)                                 :: self
  double precision                                , intent(in   ), dimension(self%parameterCount) :: stateVector
  logical                                         , allocatable  , dimension(:                  ) :: outlierMask
  integer                                                                                         :: i

  allocate(outlierMask(0:mpiSelf%count()-1))
  do i=0,mpiSelf%count()-1
     outlierMask(i)=self%simulationConvergence%chainIsOutlier(i)
  end do
  call self%simulationState%update(stateVector,self%logging(),self%simulationConvergence%isConverged(),outlierMask)
  return
end subroutine simulatorDifferentialEvolutionUpdate

integer function simulatorDifferentialChainSelect(self,randomNumberGenerator,blockedChains)
  !% Select a chain at random, optionally excluding blocked chains.
  use Pseudo_Random
  use MPI_Utilities
  implicit none
  class  (simulatorDifferentialEvolution), intent(inout)                         :: self
  type   (pseudoRandom                  ), intent(inout)                         :: randomNumberGenerator
  integer                                , intent(in   ), dimension(:), optional :: blockedChains
  logical                                                                        :: accept
  
  accept=.false.
  do while (.not.accept)
     simulatorDifferentialChainSelect=                                        &
          &      min(                                                         &
          &          int(                                                     &
          &               dble(mpiSelf%count())                               &
          &              *randomNumberGenerator%sample(mpiRankOffset=.true.)  &
          &             )                                                   , &
          &                    mpiSelf%count()                                &
          &          -1                                                       &
          &         )
     accept=.true.
     if (.not.self%sampleOutliers.and.self%isConverged)                                              &
          & accept=                                                                                  &
          &              self%simulationConvergence%chainIsOutlier(                  mpiSelf%rank()) &
          &        .or.                                                                              &
          &         .not.self%simulationConvergence%chainIsOutlier(simulatorDifferentialChainSelect)
     if (present(blockedChains)) accept=accept.and.all(simulatorDifferentialChainSelect /= blockedChains)
  end do
  return
end function simulatorDifferentialChainSelect

logical function simulatorDifferentialEvolutionLogging(self)
  !% Specifies whether or not the current state should be logged to file during differential evolution.
  implicit none
  class(simulatorDifferentialEvolution), intent(inout) :: self
  !GCC$ attributes unused :: self
  
  simulatorDifferentialEvolutionLogging=.true.
  return
end function simulatorDifferentialEvolutionLogging

subroutine simulatorDifferentialEvolutionPosterior(self,simulationState,logLikelihoodCurrent,logPriorCurrent,logPosterior,logLikelihood,logLikelihoodVariance,timeEvaluate,timeEvaluatePrevious)
  !% Return the log of the posterior for the current state.
  use MPI_Utilities
  use Sort
  use Galacticus_Error
  use Galacticus_Display
  implicit none
  class           (simulatorDifferentialEvolution), intent(inout)               :: self
  class           (state                         ), intent(inout)               :: simulationState
  double precision                                , intent(in   )               :: logLikelihoodCurrent  , logPriorCurrent
  double precision                                , intent(  out)               :: logLikelihood         , logPosterior            , &
       &                                                                           logLikelihoodVariance
  real                                            , intent(inout)               :: timeEvaluate
  real                                            , intent(in   )               :: timeEvaluatePrevious
  double precision                                , dimension(:  ), allocatable :: timesEvaluate         , nodeWork                , &
       &                                                                           stateVectorSelf       , timesEvaluateActual
  double precision                                , dimension(:,:), allocatable :: stateVectorWork
  integer                                         , dimension(:  ), allocatable :: processToProcess      , processFromProcess      , &
       &                                                                           timesEvaluateOrder    , nodeWorkOrder
  integer                                         , dimension(1  )              :: chainIndexSelf
  integer                                         , dimension(1,1)              :: chainIndexWork
  double precision                                , dimension(1,1)              :: logLikelihoodSelf     , logLikelihoodCurrentWork, &
       &                                                                           logPriorCurrentWork   , logPriorWork            , &
       &                                                                           timeEvaluateSelf
  double precision                                , dimension(1  )              :: logLikelihoodWork     , logLikelihoodCurrentSelf, &
       &                                                                           logPriorCurrentSelf   , logPriorSelf            , &
       &                                                                           timeEvaluateWork
  double precision                                                              :: logPrior              , timeEvaluateEffective
  integer                                                                       :: i                     , processTrial            , &
       &                                                                           nodeTrial
  type            (varying_string                )                              :: message
  character       (len=10                        )                              :: label
  
  ! Evaluate the proposed prior.
  logPrior=priorsEvaluateLog(self%parameterPriors,simulationState)
  ! Gather timing data from all chains.
  allocate(timesEvaluate      (0:mpiSelf%count()-1))
  allocate(timesEvaluateActual(0:mpiSelf%count()-1))
  allocate(processToProcess   (0:mpiSelf%count()-1))
  allocate(processFromProcess (0:mpiSelf%count()-1))
  allocate(timesEvaluateOrder (0:mpiSelf%count()-1))
  timeEvaluateEffective=timeEvaluatePrevious
  if     (                                                                     &
       &  .not.self%modelLikelihood%willEvaluate(                              &
       &                                         simulationState             , &
       &                                         self%parameterMappings      , &
       &                                         self%simulationConvergence  , &
       &                                         self%temperature          (), &
       &                                         logLikelihoodCurrent        , &
       &                                         logPriorCurrent             , &
       &                                         logPrior                      &
       &                                        )                              &
       & )                                                                     &
       & timeEvaluateEffective=0.0d0
  timesEvaluate=mpiSelf%gather(dble(timeEvaluateEffective))
  ! If previous time estimate is negative, don't do load balancing.
  if (mpiSelf%isMaster() .and. mod(self%simulationState%count(),self%reportCount) == 0) call Galacticus_Display_Indent('Load balancing report')
  if (any(timesEvaluate < 0.0d0)) then
     forall(i=0:mpiSelf%count()-1)
        processToProcess  (i)=i
        processFromProcess(i)=i
     end forall
     if (mpiSelf%isMaster() .and. mod(self%simulationState%count(),self%reportCount) == 0) call Galacticus_Display_Message('Not performing load balancing - missing work cost data')
  else
     ! Distribute tasks across nodes.
     timesEvaluateOrder=Sort_Index_Do(timesEvaluate)-1
     processToProcess=-1
     allocate(nodeWork     (mpiSelf%nodeCount()))
     allocate(nodeWorkOrder(mpiSelf%nodeCount()))
     nodeWork=0.0d0
     do i=mpiSelf%count()-1,0,-1
        nodeWorkOrder=Sort_Index_Do(nodeWork)
        do nodeTrial=1,mpiSelf%nodeCount()
           do processTrial=0,mpiSelf%count()-1
              if (mpiSelf%nodeAffinity(processTrial) == nodeWorkOrder(nodeTrial) .and. .not.any(processToProcess == processTrial)) then
                 processToProcess  (timesEvaluateOrder(i))=processTrial
                 processFromProcess(processTrial)=timesEvaluateOrder(i)
                 nodeWork(nodeWorkOrder(nodeTrial))=nodeWork(nodeWorkOrder(nodeTrial))+timesEvaluate(timesEvaluateOrder(i))
                 exit
              end if
           end do
           if (processToProcess(timesEvaluateOrder(i)) >= 0) exit
        end do
        if (processToProcess(timesEvaluateOrder(i)) < 0) call Galacticus_Error_Report('simulatorDifferentialEvolutionPosterior','failed to assign task to process')
     end do
     ! Report.
     if (mpiSelf%isMaster() .and. mod(self%simulationState%count(),self%reportCount) == 0) then
        call Galacticus_Display_Indent('Chain redistribution:')
        do i=0,mpiSelf%count()-1
           write (label,'(i4.4)') i
           message='Chain '//trim(label)//' -> process/node '
           write (label,'(i4.4)') processToProcess(i)
           message=message//trim(label)//'/'
           write (label,'(i4.4)') mpiSelf%nodeAffinity(processToProcess(i))
           message=message//trim(label)//' (work = '
           write (label,'(f9.2)') timesEvaluate(i)
           message=message//trim(label)//')'
           call Galacticus_Display_Message(message)
        end do
        call Galacticus_Display_Unindent('done')
        call Galacticus_Display_Indent('Node work loads:')
        do i=1,size(nodeWork)
           write (label,'(i4.4)') i
           message='Node '//trim(label)//': work = '
           write (label,'(f9.2)') nodeWork(i)
           message=message//trim(label)
           call Galacticus_Display_Message(message)
        end do
        call Galacticus_Display_Unindent('done')
     end if
  end if
  if (mpiSelf%isMaster() .and. mod(self%simulationState%count(),self%reportCount) == 0) call Galacticus_Display_Unindent('done')
  ! Get state vector, chain index, current likelihood, current prior and proposed prior.
  allocate(stateVectorSelf(self%parameterCount  ))
  allocate(stateVectorWork(self%parameterCount,1))
  stateVectorSelf         =simulationState%get       ()
  chainIndexSelf          =simulationState%chainIndex()
  logLikelihoodCurrentSelf=logLikelihoodCurrent
  logPriorCurrentSelf     =logPriorCurrent
  logPriorSelf            =logPrior
  stateVectorWork         =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),stateVectorSelf         )
  chainIndexWork          =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),chainIndexSelf          )
  logLikelihoodCurrentWork=mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),logLikelihoodCurrentSelf)
  logPriorCurrentWork     =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),logPriorCurrentSelf     )
  logPriorWork            =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),logPriorSelf            )
  ! Set state and chain index.
  call simulationState%update       (stateVectorWork(:,1),logState=.false.,isConverged=.false.)
  call simulationState%chainIndexSet( chainIndexWork(1,1)                                     )
  ! Evaluate the likelihood.
  logLikelihood=self%modelLikelihood%evaluate(simulationState,self%parameterMappings,self%simulationConvergence,self%temperature(),logLikelihoodCurrentWork(1,1),logPriorCurrentWork(1,1),logPriorWork(1,1),timeEvaluate,logLikelihoodVariance=logLikelihoodVariance)
  call mpiBarrier()
  ! Distribute likelihoods back to origins.
  logLikelihoodWork=                                                                         logLikelihood
  logLikelihoodSelf=mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     logLikelihoodWork      )
  logLikelihood    =                                                                         logLikelihoodSelf(1,1)
  ! Distribute likelihood variances back to origins.
  logLikelihoodWork    =                                                                         logLikelihoodVariance
  logLikelihoodSelf    =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     logLikelihoodWork      )
  logLikelihoodVariance=                                                                         logLikelihoodSelf(1,1)
  ! Distribute evaluation times back to origins.
  timeEvaluateWork =                                                                    dble(timeEvaluate          )
  timeEvaluateSelf =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     timeEvaluateWork       )
  timeEvaluate     =                                                                    real(timeEvaluateSelf (1,1))
  ! Restore state and chain index.
  call simulationState%update       (stateVectorSelf   ,logState=.false.,isConverged=.false.)
  call simulationState%chainIndexSet( chainIndexSelf(1)                                     )
  ! Compute the log posterior.
  logPosterior=logPrior+logLikelihood
  ! Gather actual evaluation times and report.
  timesEvaluateActual=mpiSelf%gather(dble(timeEvaluate))
  if (mpiSelf%isMaster() .and. mod(self%simulationState%count(),self%reportCount) == 0) then
     call Galacticus_Display_Indent('Node work done vs. expected:')
     do i=0,mpiSelf%count()-1
        write (label,'(i4.4)') i
        message='Node '//trim(label)//': work (actual/estimated) = '
        write (label,'(f9.2)') timesEvaluateActual(i)
        message=message//trim(label)//" / "
        write (label,'(f9.2)') timesEvaluate      (i)
        message=message//trim(label)
        call Galacticus_Display_Message(message)
     end do
     call Galacticus_Display_Unindent('done')
  end if
  return
end subroutine simulatorDifferentialEvolutionPosterior

double precision function simulatorDifferentialEvolutionStepSize(self)
  !% Return the step size parameter, $\gamma$, for a differential evolution step.
  implicit none
  class(simulatorDifferentialEvolution), intent(inout) :: self
  
  if (mod(self%simulationState%count(),self%stateSwapCount) == 0) then
     ! Every self%stateSwapCount steps, set gamma=1 to allow interchange of chains.
     simulatorDifferentialEvolutionStepSize=1.0d0
  else
     ! Otherwise, use the step-size algorithm.
     simulatorDifferentialEvolutionStepSize=self%                                               &
          &                                      proposalSize%                                  &
          &                                                   gamma(                            &
          &                                                         self%simulationState      , &
          &                                                         self%simulationConvergence  &
          &                                                        )
  end if
  return
end function simulatorDifferentialEvolutionStepSize

logical function simulatorDifferentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,logLikelihoodVariance,logLikelihoodVarianceProposed,randomNumberGenerator)
  !% Return whether or not to accept a proposal.
  use Pseudo_Random
  use MPI_Utilities
  implicit none
  class           (simulatorDifferentialEvolution), intent(inout) :: self
  double precision                                , intent(in   ) :: logPosterior         , logPosteriorProposed         , &
       &                                                             logLikelihoodVariance, logLikelihoodVarianceProposed
  type            (pseudoRandom                  ), intent(inout) :: randomNumberGenerator
  double precision                                                :: x
  !GCC$ attributes unused :: self, logLikelihoodVariance, logLikelihoodVarianceProposed
  
  ! Decide whether to take step.
  x=randomNumberGenerator%sample(mpiRankOffset=.true.)
  simulatorDifferentialEvolutionAcceptProposal=                           &
       &   logPosteriorProposed >      logPosterior                       &
       &  .or.                                                            &
       &   x                    < exp(-logPosterior+logPosteriorProposed)
  return
end function simulatorDifferentialEvolutionAcceptProposal

double precision function simulatorDifferentialEvolutionTemperature(self)
  !% Return the temperature.
  implicit none
  class(simulatorDifferentialEvolution), intent(inout) :: self
  !GCC$ attributes unused :: self

  simulatorDifferentialEvolutionTemperature=1.0d0
  return
end function simulatorDifferentialEvolutionTemperature
