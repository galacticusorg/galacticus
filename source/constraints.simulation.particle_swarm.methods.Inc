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

!% Defines methods for the ``particleSwarm'' simulator class.

function simulatorParticleSwarmConstructor(parameterPriors,parameterMappings,modelLikelihood,simulationConvergence,simulationStoppingCriterion&
     &,simulationState,simulationStateInitializor,stepsMaximum,logFileRoot,logFlushCount,reportCount,inertiaWeight,accelerationCoefficientPersonal,accelerationCoefficientGlobal,velocityCoefficient)
  !% Constructor for ``particleSwarm'' simulator class.
  implicit none
  type            (simulatorParticleSwarm)                                              :: simulatorParticleSwarmConstructor
  type            (prior                         ), intent(in   ), target, dimension(:) :: parameterPriors
  type            (mappingList                   ), intent(in   ), target, dimension(:) :: parameterMappings
  class           (likelihood                    ), intent(in   ), target               :: modelLikelihood
  class           (convergence                   ), intent(in   ), target               :: simulationConvergence
  class           (stoppingCriterion             ), intent(in   ), target               :: simulationStoppingCriterion
  class           (state                         ), intent(in   ), target               :: simulationState
  class           (stateInitializor              ), intent(in   ), target               :: simulationStateInitializor
  character       (len=*                         ), intent(in   )                       :: logFileRoot
  integer                                         , intent(in   )                       :: stepsMaximum                     , reportCount                    , &
       &                                                                                   logFlushCount
  double precision                                , intent(in   )                       :: inertiaWeight                    , accelerationCoefficientPersonal, &
       &                                                                                   accelerationCoefficientGlobal    , velocityCoefficient

  simulatorParticleSwarmConstructor%parameterCount                 =  size(parameterPriors)
  simulatorParticleSwarmConstructor%parameterPriors                => parameterPriors
  simulatorParticleSwarmConstructor%parameterMappings              => parameterMappings
  simulatorParticleSwarmConstructor%modelLikelihood                => modelLikelihood
  simulatorParticleSwarmConstructor%simulationConvergence          => simulationConvergence
  simulatorParticleSwarmConstructor%simulationStoppingCriterion    => simulationStoppingCriterion
  simulatorParticleSwarmConstructor%simulationState                => simulationState
  simulatorParticleSwarmConstructor%simulationStateInitializor     => simulationStateInitializor
  simulatorParticleSwarmConstructor%stepsMaximum                   =  stepsMaximum
  simulatorParticleSwarmConstructor%logFlushCount                  =  logFlushCount
  simulatorParticleSwarmConstructor%reportCount                    =  reportCount
  simulatorParticleSwarmConstructor%inertiaWeight                  =  inertiaWeight
  simulatorParticleSwarmConstructor%accelerationCoefficientPersonal=  accelerationCoefficientPersonal
  simulatorParticleSwarmConstructor%accelerationCoefficientGlobal  =  accelerationCoefficientGlobal
  simulatorParticleSwarmConstructor%velocityCoefficient            =  velocityCoefficient
  simulatorParticleSwarmConstructor%logFileRoot                    =  logFileRoot
  return
end function simulatorParticleSwarmConstructor

subroutine simulatorParticleSwarmSimulate(self)
  !% Perform a particle swarm simulation.
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
  class           (simulatorParticleSwarm), intent(inout)                               :: self
  double precision                        , dimension(self%parameterCount)              :: stateVector             , positionMinimum       , &
       &                                                                                   positionMaximum         , velocityParticle      , &
       &                                                                                   velocityMaximum         , stateBestPersonal     , &
       &                                                                                   stateBestGlobal
  double precision                        , dimension(:,:)                , allocatable :: stateVectors
  double precision                        , dimension(:  )                , allocatable :: logPosteriorsBestGlobal
  real                                                                                  :: timePreEvaluate         , timePostEvaluate      , &
       &                                                                                   timeEvaluate            , timeEvaluatePrevious
  double precision                                                                      :: timeEvaluateInitial     , logPosterior          , &
       &                                                                                   logPosteriorBestPersonal, logPosteriorBestGlobal
  type            (pseudoRandom          )                                              :: randomNumberGenerator
  type            (varying_string        )                                              :: logFileName             , message
  integer                                                                               :: logFileUnit             , convergedAtStep       , &
       &                                                                                   convergenceFileUnit     , i                     , &
       &                                                                                   ioStatus
  logical                                                                               :: isConverged
  
  ! Initialize the particle state vector.
  logPosterior=logImpossible
  do while (logPosterior <= logImpossible)
     ! Initialize particle to some state vector.
     call self%simulationStateInitializor%initialize(self%simulationState,self%parameterPriors,self%parametermappings,timeEvaluateInitial)
     ! Evaluate the posterior in the initial state.
     timeEvaluate        =-1.0
     timeEvaluatePrevious=real(timeEvaluateInitial)
     call CPU_Time(timePreEvaluate )
     call self%posterior(self%simulationState,logPosterior,timeEvaluate,timeEvaluatePrevious)
     call CPU_Time(timePostEvaluate)
     if (timeEvaluate < 0.0) timeEvaluate=timePostEvaluate-timePreEvaluate
     timeEvaluatePrevious=timeEvaluate
  end do
  ! Set the personal best state to the initial state.
  logPosteriorBestPersonal=logPosterior
  stateBestPersonal       =self%simulationState%get()
  ! Set global best state.
  allocate(stateVectors           (self%parameterCount,mpiSelf%count()))
  allocate(logPosteriorsBestGlobal(                    mpiSelf%count()))
  logPosteriorBestGlobal  =logPosterior
  stateBestGlobal         =self%simulationState%get()
  logPosteriorsBestGlobal =mpiSelf%gather(logPosteriorBestGlobal)
  stateVectors            =mpiSelf%gather(stateBestGlobal       )
  stateBestGlobal         =stateVectors(:,maxloc(logPosteriorsBestGlobal,dim=1))
  logPosteriorBestGlobal  =               maxval(logPosteriorsBestGlobal      )
  ! Set an initial velocity for the particle.
  do i=1,self%parameterCount     
     positionMinimum (i)=self%parameterPriors(i)%minimum()
     positionMaximum (i)=self%parameterPriors(i)%maximum()
     velocityMaximum (i)=self%velocityCoefficient*(positionMaximum(i)-positionMinimum(i))
     velocityParticle(i)=(2.0d0*randomNumberGenerator%sample()-1.0d0)*velocityMaximum(i)
  end do  
  ! Begin the simulation.
  logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.log'
  open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
  isConverged=.false. 
  do while (                                                                                                                  &
       &          self%simulationState            %count(                                               ) < self%stepsMaximum &
       &    .and.                                                                                                             &
       &     .not.self%simulationStoppingCriterion%stop (self%simulationState,self%simulationConvergence)                     &
       &   )
     ! Get the current particle state.
     stateVector=self%simulationState%get()
     ! Update the state vector
     stateVector=stateVector+velocityParticle
     do i=1,self%parameterCount
        if (stateVector(i) < positionMinimum(i)) then
           velocityParticle(i)=-velocityParticle(i)
           stateVector     (i)=+positionMinimum (i)
        end if
        if (stateVector(i) > positionMaximum(i)) then
           velocityParticle(i)=-velocityParticle(i)
           stateVector     (i)=+positionMaximum (i)
        end if
     end do
     call self%simulationState%update(stateVector,.true.,self%simulationConvergence%isConverged())
     ! Update the velocity.
     velocityParticle=+self%inertiaWeight                   &
          &           *velocityParticle                     &
          &           +self%accelerationCoefficientPersonal &
          &           *randomNumberGenerator%sample()       &
          &           *(                                    &
          &             +stateBestPersonal                  &
          &             -stateVector                        &
          &            )                                    &
          &           +self%accelerationCoefficientGlobal   &
          &           *randomNumberGenerator%sample()       &
          &           *(                                    &
          &             +stateBestGlobal                    &
          &             -stateVector                        &
          &            )
     where (velocityParticle > velocityMaximum)
        velocityParticle=+velocityMaximum
     end where
     where (velocityParticle < -velocityMaximum)
        velocityParticle=-velocityMaximum
     end where
     ! Evaluate posterior.
     timeEvaluatePrevious=timeEvaluate
     timeEvaluate        =-1.0
     call CPU_Time(timePreEvaluate )
     call self%posterior(self%simulationState,logPosterior,timeEvaluate,timeEvaluatePrevious)
     call CPU_Time(timePostEvaluate)
     if (timeEvaluate < 0.0) timeEvaluate=timePostEvaluate-timePreEvaluate
     ! Update personal best state.
     if (logPosterior > logPosteriorBestPersonal) then
        logPosteriorBestPersonal=logPosterior
        stateBestPersonal       =stateVector
     end if
     ! Update global best state.
     if (logPosterior > logPosteriorBestGlobal  ) then
        logPosteriorBestGlobal  =logPosterior
        stateBestGlobal         =stateVector
     end if
     logPosteriorsBestGlobal=mpiSelf%gather(logPosteriorBestGlobal)
     stateVectors           =mpiSelf%gather(stateBestGlobal         )
     stateBestGlobal        =stateVectors(:,maxloc(logPosteriorsBestGlobal,dim=1))
     logPosteriorBestGlobal =               maxval(logPosteriorsBestGlobal      )
     ! Unmap parameters and write to log file.
     do i=1,size(stateVector)
        stateVector(i)=self%parameterMappings(i)%thisMapping%unmap(stateVector(i))
     end do
     write (logFileUnit,*) self   %simulationState%count(), &
          &                mpiSelf%rank                 (), &
          &                timeEvaluate                   , &
          &                isConverged                    , &
          &                logPosterior                   , &
          &                stateVector
     if (mod(self%simulationState%count(),self%logFlushCount) == 0) call flush(logFileUnit)
    ! Repeat.
     call mpiBarrier()
     ! Test for convergence.
     if (.not.isConverged) then
        isConverged=self%simulationConvergence%isConverged(self%simulationState,logPosterior)
        if (isConverged) then
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
end subroutine simulatorParticleSwarmSimulate

subroutine simulatorParticleSwarmPosterior(self,simulationState,logPosterior,timeEvaluate,timeEvaluatePrevious)
  !% Return the log of the posterior for the current state.
  use Constraints_Constants
  use MPI_Utilities
  use Sort
  use Galacticus_Error
  use Galacticus_Display
  implicit none
  class           (simulatorParticleSwarm), intent(inout)               :: self
  class           (state                 ), intent(inout)               :: simulationState
  double precision                        , intent(  out)               :: logPosterior
  real                                    , intent(inout)               :: timeEvaluate
  real                                    , intent(in   )               :: timeEvaluatePrevious
  double precision                        , dimension(:  ), allocatable :: timesEvaluate             , nodeWork                , &
       &                                                                   stateVectorSelf           , timesEvaluateActual
  double precision                        , dimension(:,:), allocatable :: stateVectorWork
  integer                                 , dimension(:  ), allocatable :: processToProcess          , processFromProcess      , &
       &                                                                   timesEvaluateOrder        , nodeWorkOrder
  integer                                 , dimension(1  )              :: particleIndexSelf
  integer                                 , dimension(1,1)              :: particleIndexWork
  double precision                        , dimension(1,1)              :: logLikelihoodSelf         , logPriorWork            , &
       &                                                                   timeEvaluateSelf
  double precision                        , dimension(1  )              :: logLikelihoodWork         , logPriorSelf            , &
       &                                                                   timeEvaluateWork
  double precision                        , parameter                   :: temperature         =1.0d0
  double precision                                                      :: logPrior                  , timeEvaluateEffective   , &
       &                                                                   logLikelihood
  integer                                                               :: i                         , processTrial            , &
       &                                                                   nodeTrial
  type            (varying_string        )                              :: message
  character       (len=10                )                              :: label
  
  ! Evaluate the proposed prior.
  logPrior=priorsEvaluateLog(self%parameterPriors,simulationState)
  ! Gather timing data from all particles.
  allocate(timesEvaluate      (0:mpiSelf%count()-1))
  allocate(timesEvaluateActual(0:mpiSelf%count()-1))
  allocate(processToProcess   (0:mpiSelf%count()-1))
  allocate(processFromProcess (0:mpiSelf%count()-1))
  allocate(timesEvaluateOrder (0:mpiSelf%count()-1))
  timeEvaluateEffective=timeEvaluatePrevious
  if     (                                                                   &
       &  .not.self%modelLikelihood%willEvaluate(                            &
       &                                         simulationState           , &
       &                                         self%parameterMappings    , &
       &                                         self%simulationConvergence, &
       &                                         temperature               , &
       &                                         logImpossible             , &
       &                                         logImpossible             , &
       &                                         logPrior                    &
       &                                        )                            &
       & )                                                                   &
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
           if (processToProcess(timesevaluateorder(i)) >= 0) exit
        end do
        if (processToProcess(timesevaluateorder(i)) < 0) call Galacticus_Error_Report('simulatorParticleSwarmPosterior','failed to assign task to process')
     end do
     ! Report.
     if (mpiSelf%isMaster() .and. mod(self%simulationState%count(),self%reportCount) == 0) then
        call Galacticus_Display_Indent('Particle redistribution:')
        do i=0,mpiSelf%count()-1
           write (label,'(i4.4)') i
           message='Particle '//trim(label)//' -> process/node '
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
  ! Get state vector, particle index, and prior.
  allocate(stateVectorSelf(self%parameterCount  ))
  allocate(stateVectorWork(self%parameterCount,1))
  stateVectorSelf         =simulationState%get       ()
  particleIndexSelf       =simulationState%chainIndex()
  logPriorSelf            =logPrior
  stateVectorWork         =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),stateVectorSelf  )
  particleIndexWork       =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),particleIndexSelf)
  logPriorWork            =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),logPriorSelf     )
  ! Set state and particle index.
  call simulationState%update       (  stateVectorWork(:,1),logState=.false.,isConverged=.false.)
  call simulationState%chainIndexSet(particleIndexWork(1,1)                                     )
  ! Evaluate the likelihood.
  logLikelihood=self%modelLikelihood%evaluate(                                 &
       &                                      simulationState                , &
       &                                      self%parameterMappings         , &
       &                                      self%simulationConvergence     , &
       &                                      temperature                    , &
       &                                      logImpossible                  , &
       &                                      logImpossible                  , &
       &                                      logPriorWork              (1,1), &
       &                                      timeEvaluate                     &
       &                                     )
  call mpiBarrier()
  ! Distribute likelihoods back to origins.
  logLikelihoodWork    =                                                                         logLikelihood
  logLikelihoodSelf    =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     logLikelihoodWork         )
  logLikelihood        =                                                                         logLikelihoodSelf    (1,1)
  ! Distribute evaluation times back to origins.
  timeEvaluateWork    =                                                                     dble(timeEvaluate             )
  timeEvaluateSelf    =mpiSelf%requestData (processToProcess(mpiSelf%rank():mpiSelf%rank()),     timeEvaluateWork          )
  timeEvaluate        =                                                                     real(timeEvaluateSelf     (1,1))
  ! Restore state and particle index.
  call simulationState%update       (  stateVectorSelf   ,logState=.false.,isConverged=.false.)
  call simulationState%chainIndexSet(particleIndexSelf(1)                                     )
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
end subroutine simulatorParticleSwarmPosterior
