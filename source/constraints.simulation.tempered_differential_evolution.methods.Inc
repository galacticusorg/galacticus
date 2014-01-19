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

!% Defines methods for the ``differentialEvolution'' simulator class.

function simulatorTemperedDifferentialEvolutionConstructor(parameterPriors,randomDistributions,modelLikelihood,simulationConvergence &
          &,simulationState,proposalSize,proposalSizeTemperatureExponent,randomJump,stepsMaximum,stepsPostConvergence&
          &,acceptanceAverageCount,stateSwapCount,logFileRoot,temperatureMaximum ,untemperedStepCount,temperingLevelCount&
          &,stepsPerLevel)
  !% Constructor for ``temperedDifferentialEvolution'' simulator class.
  implicit none
  type            (simulatorTemperedDifferentialEvolution)                                      :: simulatorTemperedDifferentialEvolutionConstructor
  type            (prior                                 ), intent(in   ), target, dimension(:) :: parameterPriors
  type            (distributionList                      ), intent(in   ), target, dimension(:) :: randomDistributions
  class           (likelihood                            ), intent(in   ), target               :: modelLikelihood
  class           (convergence                           ), intent(in   ), target               :: simulationConvergence
  class           (state                                 ), intent(in   ), target               :: simulationState
  class           (deProposalSize                        ), intent(in   ), target               :: proposalSize
  class           (dePropSizeTempExp                     ), intent(in   ), target               :: proposalSizeTemperatureExponent
  class           (deRandomJump                          ), intent(in   ), target               :: randomJump
  double precision                                        , intent(in   )                       :: temperatureMaximum
  integer                                                 , intent(in   )                       :: stepsMaximum                   , stepsPostConvergence    , &
       &                                                                                           acceptanceAverageCount         , untemperedStepCount     , &
       &                                                                                           temperingLevelCount            , stepsPerLevel           , &
       &                                                                                           stateSwapCount
  character       (len=*                                 ), intent(in   )                       :: logFileRoot
  integer                                                                                       :: i

  simulatorTemperedDifferentialEvolutionConstructor%parameterCount                  =  size(parameterPriors)
  simulatorTemperedDifferentialEvolutionConstructor%parameterPriors                 => parameterPriors
  simulatorTemperedDifferentialEvolutionConstructor%randomDistributions             => randomDistributions
  simulatorTemperedDifferentialEvolutionConstructor%modelLikelihood                 => modelLikelihood
  simulatorTemperedDifferentialEvolutionConstructor%simulationConvergence           => simulationConvergence
  simulatorTemperedDifferentialEvolutionConstructor%simulationState                 => simulationState
  simulatorTemperedDifferentialEvolutionConstructor%proposalSize                    => proposalSize
  simulatorTemperedDifferentialEvolutionConstructor%proposalSizeTemperatureExponent => proposalSizeTemperatureExponent
  simulatorTemperedDifferentialEvolutionConstructor%randomJump                      => randomJump
  simulatorTemperedDifferentialEvolutionConstructor%stepsMaximum                    =  stepsMaximum
  simulatorTemperedDifferentialEvolutionConstructor%stepsPostConvergence            =  stepsPostConvergence
  simulatorTemperedDifferentialEvolutionConstructor%acceptanceAverageCount          =  acceptanceAverageCount
  simulatorTemperedDifferentialEvolutionConstructor%stateSwapCount                  =  stateSwapCount
  simulatorTemperedDifferentialEvolutionConstructor%logFileRoot                     =  logFileRoot
  simulatorTemperedDifferentialEvolutionConstructor%temperatureMaximum              =  temperatureMaximum
  simulatorTemperedDifferentialEvolutionConstructor%untemperedStepCount             =  untemperedStepCount
  simulatorTemperedDifferentialEvolutionConstructor%temperingLevelCount             =  temperingLevelCount
  simulatorTemperedDifferentialEvolutionConstructor%stepsPerLevel                   =  stepsPerLevel
  simulatorTemperedDifferentialEvolutionConstructor%temperingStep                   =  0
  simulatorTemperedDifferentialEvolutionConstructor%temperingLevelMonotonic         =  0
  allocate(simulatorTemperedDifferentialEvolutionConstructor%temperatures(temperingLevelCount))
  allocate(stateSimple :: simulatorTemperedDifferentialEvolutionConstructor%temperedStates(temperingLevelCount))
  do i=1,temperingLevelCount
     simulatorTemperedDifferentialEvolutionConstructor%temperatures(i) &
          &  =exp(                                                     &
          &        log (temperatureMaximum )                           &
          &       *dble(i                  )                           &
          &       /dble(temperingLevelCount)                           &
          &      )
     select type (s => simulatorTemperedDifferentialEvolutionConstructor%temperedStates(i))
     type is (stateSimple)
        s=stateSimple(size(parameterPriors),acceptanceAverageCount)
     end select
  end do
  return
end function simulatorTemperedDifferentialEvolutionConstructor

logical function simulatorTemperedDifferentialEvolutionLogging(self)
  !% Specifies whether or not the current state should be logged to file during differential evolution.
  implicit none
  class(simulatorTemperedDifferentialEvolution), intent(inout) :: self

  simulatorTemperedDifferentialEvolutionLogging=(self%temperingLevelMonotonic == 0)
  return
end function simulatorTemperedDifferentialEvolutionLogging

subroutine simulatorTemperedDifferentialEvolutionUpdate(self,stateVector)
  !% Update the differential evolution simulator state.
  use MPI_Utilities
  use Galacticus_Display
  use String_Handling
  use ISO_Varying_String
  implicit none
  class           (simulatorTemperedDifferentialEvolution), intent(inout)                                   :: self
  double precision                                        , intent(in   ), dimension(self%parameterCount  ) :: stateVector
  integer                                                                                                   :: i             , temperingLevelSaved
  logical                                                                                                   :: levelChanged
  double precision                                                                                          :: acceptanceRate, temperature        , &
       &                                                                                                       stepSize
  character       (len=30                                )                                                  :: label
  type            (varying_string                        )                                                  :: message

  ! Update the simulation state.
  call self%simulationState%update(stateVector,self%logging())
  ! Update tempering step count and level as necessary.
  levelChanged=.false.
  self%temperingStep=self%temperingStep+1
  if (self%temperingLevelMonotonic == 0) then
     ! Currently not tempering: check if we've run all required untempered steps.
     if (self%temperingStep > self%untemperedStepCount) then
        ! We have, switch to the first tempered level.
        self%temperingStep          =0
        self%temperingLevelMonotonic=1
        levelChanged                =.true.
     end if
  else
     ! Update the tempered level states.
     call self%temperedStates(self%level())%update(stateVector,.true.)
     ! Currently tempering: check if we've run all required steps at this level.
     if (self%temperingStep > self%stepsPerLevel) then
        ! We have, move to the next tempering level. Note that we run tempering levels from 1 to
        ! 2*temperingLevelCount-1 - levels above temperingLevelCount represent the cooling phase.
        self%temperingStep          =0
        self%temperingLevelMonotonic=self%temperingLevelMonotonic+1
        levelChanged                =.true.
        ! Check if we've finished tempering and switch back to untempered evolution is so.
        if (self%temperingLevelMonotonic > 2*self%temperingLevelCount-1) self%temperingLevelMonotonic=0 
     end if
  end if
  ! Check for change in level.
  if (levelChanged) then
     if (mpiSelf%isMaster().and.Galacticus_Verbosity_Level() >= verbosityInfo) then
        write (label,'(f8.1)') self%temperature()
        message='Tempering state: level='
        message=message//self%level()//'; temperature='//trim(label)
        call Galacticus_Display_Message(message)
     end if
     if (self%temperingLevelMonotonic == 0) then
        ! We've just returned to the untempered level. Report on acceptance rates in tempered levels.
        if (Galacticus_Verbosity_Level() >= verbosityInfo) then
           if (mpiSelf%isMaster()) then
              call Galacticus_Display_Indent('Acceptance rates in tempered levels')
              call Galacticus_Display_Message('Level Temperature  Gamma  Rate')
              call Galacticus_Display_Message('------------------------------')
           end if
           ! Store the current tempering level so that we can restore it below.
           temperingLevelSaved=self%temperingLevelMonotonic
           do i=1,self%temperingLevelCount
              self%temperingLevelMonotonic=i
              acceptanceRate=mpiSelf%average(self%temperedStates(i)%acceptanceRate())
              temperature   =                self                  %temperature   ()
              stepSize      =                self                  %stepSize      ()
              if (mpiSelf%isMaster())  then
                 write (label,'(2x,i3,4x,f8.1,1x,f6.3,1x,f5.3)') i,temperature,stepSize,acceptanceRate
                 call Galacticus_Display_Message(label)
              end if
           end do
           if (mpiSelf%isMaster()) call Galacticus_Display_Unindent('done')
           ! Restore the tempering level to the original.
           self%temperingLevelMonotonic=temperingLevelSaved
        end if
     else
        ! We're in a tempering level, update the state of the new level to the current state, without logging.
        call self%temperedStates(self%level())%update(stateVector,.false.)
     end if
  end if
  return
end subroutine simulatorTemperedDifferentialEvolutionUpdate

integer function simulatorTemperedDifferentialEvolutionLevel(self)
  !% Return the actual tempering level.
  implicit none
  class(simulatorTemperedDifferentialEvolution), intent(inout) :: self

  simulatorTemperedDifferentialEvolutionLevel=self%temperingLevelMonotonic
  if (self%temperingLevelMonotonic > self%temperingLevelCount)                     &
       & simulatorTemperedDifferentialEvolutionLevel= 2                            &
       &                                             *self%temperingLevelCount     &
       &                                             -self%temperingLevelMonotonic
  return
end function simulatorTemperedDifferentialEvolutionLevel

double precision function simulatorTemperedDifferentialEvolutionStepSize(self)
  !% Return the step size parameter, $\gamma$, for a differential evolution step.
  implicit none
  class           (simulatorTemperedDifferentialEvolution), intent(inout) :: self
  double precision                                                        :: gammaBoostFactor

  if (mod(self%simulationState%count(),self%stateSwapCount) == 0 .and. self%level() == 0) then
     ! Every self%stateSwapCount steps, set gamma=1 to allow interchange of chains.
     simulatorTemperedDifferentialEvolutionStepSize=1.0d0
  else
     gammaBoostFactor=self%temperature()**self%proposalSizeTemperatureExponent%exponent(self%temperedStates,self%temperatures,self%simulationState,self%simulationConvergence)
     simulatorTemperedDifferentialEvolutionStepSize=gammaBoostFactor*self%simulatorDifferentialEvolution%stepSize()
  end if
  return
end function simulatorTemperedDifferentialEvolutionStepSize

double precision function simulatorTemperedDifferentialEvolutionTemperature(self)
  !% Return the temperature.
  implicit none
  class(simulatorTemperedDifferentialEvolution), intent(inout) :: self

  if (self%level() == 0) then
     simulatorTemperedDifferentialEvolutionTemperature=1.0d0
  else
     simulatorTemperedDifferentialEvolutionTemperature=self%temperatures(self%level())
  end if
  return
end function simulatorTemperedDifferentialEvolutionTemperature

logical function simulatorTemperedDifferentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,randomNumberGenerator)
  !% Return whether or not to accept a proposal.
  use Pseudo_Random
  use MPI_Utilities
  implicit none
  class           (simulatorTemperedDifferentialEvolution), intent(inout) :: self
  double precision                                        , intent(in   ) :: logPosterior         , logPosteriorProposed
  type            (pseudoRandom                          ), intent(inout) :: randomNumberGenerator
  double precision                                                        :: x

  ! Decide whether to take step.
  x=randomNumberGenerator%sample(mpiRankOffset=.true.)
  simulatorTemperedDifferentialEvolutionAcceptProposal=      &
       &   logPosteriorProposed >      logPosterior          &
       &  .or.                                               &
       &   x                    < exp(                       &
       &                              (                      &
       &                               -logPosterior         &
       &                               +logPosteriorProposed &
       &                              )                      &
       &                              /self%temperature()    &
       &                             )
  return
end function simulatorTemperedDifferentialEvolutionAcceptProposal
