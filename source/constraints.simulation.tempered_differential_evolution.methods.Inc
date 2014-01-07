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
          &,simulationState,proposalSize,stepsMaximum,stepsPostConvergence,acceptanceAverageCount,logFileRoot,temperatureMaximum&
          &,untemperedStepCount,temperingLevelCount,stepsPerLevel,gammaTemperatureExponent)
  !% Constructor for ``differentialEvolution'' simulator class.
  implicit none
  type            (simulatorTemperedDifferentialEvolution)                                      :: simulatorTemperedDifferentialEvolutionConstructor
  type            (prior                                 ), intent(in   ), target, dimension(:) :: parameterPriors
  type            (distributionList                      ), intent(in   ), target, dimension(:) :: randomDistributions
  class           (likelihood                            ), intent(in   ), target               :: modelLikelihood
  class           (convergence                           ), intent(in   ), target               :: simulationConvergence
  class           (state                                 ), intent(in   ), target               :: simulationState
  class           (deProposalSize                        ), intent(in   ), target               :: proposalSize
  double precision                                        , intent(in   )                       :: temperatureMaximum    , gammaTemperatureExponent
  integer                                                 , intent(in   )                       :: stepsMaximum          , stepsPostConvergence    , &
       &                                                                                           acceptanceAverageCount, untemperedStepCount     , &
       &                                                                                           temperingLevelCount   , stepsPerLevel
  character       (len=*                                 ), intent(in   )                       :: logFileRoot

  simulatorTemperedDifferentialEvolutionConstructor%parameterCount           =  size(parameterPriors)
  simulatorTemperedDifferentialEvolutionConstructor%parameterPriors          => parameterPriors
  simulatorTemperedDifferentialEvolutionConstructor%randomDistributions      => randomDistributions
  simulatorTemperedDifferentialEvolutionConstructor%modelLikelihood          => modelLikelihood
  simulatorTemperedDifferentialEvolutionConstructor%simulationConvergence    => simulationConvergence
  simulatorTemperedDifferentialEvolutionConstructor%simulationState          => simulationState
  simulatorTemperedDifferentialEvolutionConstructor%proposalSize             => proposalSize
  simulatorTemperedDifferentialEvolutionConstructor%stepsMaximum             =  stepsMaximum
  simulatorTemperedDifferentialEvolutionConstructor%stepsPostConvergence     =  stepsPostConvergence
  simulatorTemperedDifferentialEvolutionConstructor%acceptanceAverageCount   =  acceptanceAverageCount
  simulatorTemperedDifferentialEvolutionConstructor%logFileRoot              =  logFileRoot
  simulatorTemperedDifferentialEvolutionConstructor%temperatureMaximum       =  temperatureMaximum
  simulatorTemperedDifferentialEvolutionConstructor%untemperedStepCount      =  untemperedStepCount
  simulatorTemperedDifferentialEvolutionConstructor%temperingLevelCount      =  temperingLevelCount
  simulatorTemperedDifferentialEvolutionConstructor%stepsPerLevel            =  stepsPerLevel
  simulatorTemperedDifferentialEvolutionConstructor%gammaTemperatureExponent =  gammaTemperatureExponent
  simulatorTemperedDifferentialEvolutionConstructor%temperature              =  1.0d0
  simulatorTemperedDifferentialEvolutionConstructor%temperingStep            =  0
  simulatorTemperedDifferentialEvolutionConstructor%temperingLevel           =  0
  return
end function simulatorTemperedDifferentialEvolutionConstructor

logical function simulatorTemperedDifferentialEvolutionLogging(self)
  !% Specifies whether or not the current state should be logged to file during differential evolution.
  implicit none
  class(simulatorTemperedDifferentialEvolution), intent(inout) :: self

  simulatorTemperedDifferentialEvolutionLogging=(self%temperingLevel == 0)
  return
end function simulatorTemperedDifferentialEvolutionLogging

double precision function simulatorTemperedDifferentialEvolutionPosterior(self,simulationState)
  !% Return the log of the posterior for the current state.
  implicit none
  class(simulatorTemperedDifferentialEvolution), intent(inout) :: self
  class(state                                 ), intent(in   ) :: simulationState

  simulatorTemperedDifferentialEvolutionPosterior=                         &
       &  self%simulatorDifferentialEvolution%posterior  (simulationState) &
       & /self                               %temperature
  return
end function simulatorTemperedDifferentialEvolutionPosterior

subroutine simulatorTemperedDifferentialEvolutionUpdate(self,stateVector)
  !% Update the differential evolution simulator state.
  use MPI_Utilities
  implicit none
  class           (simulatorTemperedDifferentialEvolution), intent(inout)                                   :: self
  double precision                                        , intent(in   ), dimension(self%parameterCount  ) :: stateVector
  integer                                                                                                   :: temperingLevelActual
  logical                                                                                                   :: report

  ! Call the parent class to store the state.
  call self%simulatorDifferentialEvolution%update(stateVector)
  ! Update tempering step count and level as necessary.
  report=.false.
  self%temperingStep=self%temperingStep+1
  if (self%temperingLevel == 0) then
     ! Currently not tempering: check if we've run all required untempered steps.
     if (self%temperingStep > self%untemperedStepCount) then
        ! We have, switch to the first tempered level.
        self%temperingStep =0
        self%temperingLevel=1
        report             =.true.
     end if
  else
     ! Currently tempering: check if we've run all required steps at this level.
     if (self%temperingStep > self%stepsPerLevel) then
        ! We have, move to the next tempering level. Note that we run tempering levels from 1 to
        ! 2*temperingLevelCount-1 - levels above temperingLevelCount represent the cooling phase.
        self%temperingStep =0
        self%temperingLevel=self%temperingLevel+1
        report             =.true.
       ! Check if we've finished tempering and switch back to untempered evolution is so.
        if (self%temperingLevel > 2*self%temperingLevelCount-1) self%temperingLevel=0 
     end if
  end if
  ! Compute the current temperature.
  temperingLevelActual=self%temperingLevel
  if (self%temperingLevel > self%temperingLevelCount)   &
       & temperingLevelActual= 2                        &
       &                      *self%temperingLevelCount &
       &                      -self%temperingLevel
  self%temperature=exp(                                &
       &                log (self%temperatureMaximum ) &
       &               *dble(temperingLevelActual    ) &
       &               /dble(self%temperingLevelCount) &
       &              )
  ! Report if necessary.
  if (mpiSelf%isMaster().and.report) &
       & write (0,*) 'Tempering state: level=',self%temperingLevel,'; temperature=',self%temperature
  return
end subroutine simulatorTemperedDifferentialEvolutionUpdate

double precision function simulatorTemperedDifferentialEvolutionStepSize(self)
  !% Return the step size parameter, $\gamma$, for a differential evolution step.
  implicit none
  class           (simulatorTemperedDifferentialEvolution), intent(inout) :: self
  double precision                                                        :: gammaBoostFactor

  gammaBoostFactor=self%temperature**self%gammaTemperatureExponent
  simulatorTemperedDifferentialEvolutionStepSize=gammaBoostFactor*self%simulatorDifferentialEvolution%stepSize()
  return
end function simulatorTemperedDifferentialEvolutionStepSize
