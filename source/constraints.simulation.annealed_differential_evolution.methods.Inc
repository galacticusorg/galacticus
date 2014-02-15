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

!% Defines methods for the ``annealedDifferentialEvolution'' simulator class.

function simulatorAnnealedDifferentialEvolutionConstructor(parameterPriors,randomDistributions,modelLikelihood,simulationConvergence &
          &,simulationState,simulationStateInitializor,proposalSize,randomJump,stepsMaximum,stepsPostConvergence&
          &,acceptanceAverageCount,stateSwapCount,logFileRoot,temperatureMaximum,temperatureLevelCount)
  !% Constructor for ``annealedDifferentialEvolution'' simulator class.
  implicit none
  type            (simulatorAnnealedDifferentialEvolution)                                      :: simulatorAnnealedDifferentialEvolutionConstructor
  type            (prior                                 ), intent(in   ), target, dimension(:) :: parameterPriors
  type            (distributionList                      ), intent(in   ), target, dimension(:) :: randomDistributions
  class           (likelihood                            ), intent(in   ), target               :: modelLikelihood
  class           (convergence                           ), intent(in   ), target               :: simulationConvergence
  class           (state                                 ), intent(in   ), target               :: simulationState
  class           (stateInitializor                      ), intent(in   ), target               :: simulationStateInitializor
  class           (deProposalSize                        ), intent(in   ), target               :: proposalSize
  class           (deRandomJump                          ), intent(in   ), target               :: randomJump
  double precision                                        , intent(in   )                       :: temperatureMaximum
  integer                                                 , intent(in   )                       :: stepsMaximum                   , stepsPostConvergence, &
       &                                                                                           acceptanceAverageCount         , stateSwapCount      , &
       &                                                                                           temperatureLevelCount
  character       (len=*                                 ), intent(in   )                       :: logFileRoot
  integer                                                                                       :: i

  simulatorAnnealedDifferentialEvolutionConstructor%parameterCount             =  size(parameterPriors)
  simulatorAnnealedDifferentialEvolutionConstructor%parameterPriors            => parameterPriors
  simulatorAnnealedDifferentialEvolutionConstructor%randomDistributions        => randomDistributions
  simulatorAnnealedDifferentialEvolutionConstructor%modelLikelihood            => modelLikelihood
  simulatorAnnealedDifferentialEvolutionConstructor%simulationConvergence      => simulationConvergence
  simulatorAnnealedDifferentialEvolutionConstructor%simulationState            => simulationState
  simulatorAnnealedDifferentialEvolutionConstructor%simulationStateInitializor => simulationStateInitializor
  simulatorAnnealedDifferentialEvolutionConstructor%proposalSize               => proposalSize
  simulatorAnnealedDifferentialEvolutionConstructor%randomJump                 => randomJump
  simulatorAnnealedDifferentialEvolutionConstructor%stepsMaximum               =  stepsMaximum
  simulatorAnnealedDifferentialEvolutionConstructor%stepsPostConvergence       =  stepsPostConvergence
  simulatorAnnealedDifferentialEvolutionConstructor%acceptanceAverageCount     =  acceptanceAverageCount
  simulatorAnnealedDifferentialEvolutionConstructor%stateSwapCount             =  stateSwapCount
  simulatorAnnealedDifferentialEvolutionConstructor%logFileRoot                =  logFileRoot
  simulatorAnnealedDifferentialEvolutionConstructor%temperatureMaximum         =  temperatureMaximum
  simulatorAnnealedDifferentialEvolutionConstructor%temperatureLevelCount      =  temperatureLevelCount
  simulatorAnnealedDifferentialEvolutionConstructor%temperatureLevelCurrent    =  temperatureLevelCount
  allocate(simulatorAnnealedDifferentialEvolutionConstructor%temperatures(temperatureLevelCount))
  if (temperatureLevelCount == 1) then
     simulatorAnnealedDifferentialEvolutionConstructor%temperatures(1)=1.0d0
  else
     do i=1,temperatureLevelCount
        simulatorAnnealedDifferentialEvolutionConstructor%temperatures(i) &
             &  =exp(                                                     &
             &        log (temperatureMaximum     )                       &
             &       *dble(i                    -1)                       &
             &       /dble(temperatureLevelCount-1)                       &
             &      )
     end do
  end if
  return
end function simulatorAnnealedDifferentialEvolutionConstructor

subroutine simulatorAnnealedDifferentialEvolutionUpdate(self,stateVector)
  !% Update the differential evolution simulator state.
  use MPI_Utilities
  use Galacticus_Display
  use ISO_Varying_String
  use String_Handling
  implicit none
  class           (simulatorAnnealedDifferentialEvolution), intent(inout)                                 :: self
  double precision                                        , intent(in   ), dimension(self%parameterCount) :: stateVector
  type            (varying_string                        )                                                :: message
  character       (len=7                                 )                                                :: label

  ! Check for convergence.
  if (self%simulationConvergence%isConverged(self%simulationState,self%logPosterior).and.self%temperatureLevelCurrent > 1) then
     ! Decrease the temperature level.
     self%temperatureLevelCurrent=self%temperatureLevelCurrent-1
     ! Reset state.
     call self%simulationState      %reset()
     call self%simulationConvergence%reset()
     self%isConverged=.false.
     ! Report.
     if (mpiSelf%isMaster()) then
        write (label,'(f7.2)') self%temperature()
        message="Annealing temperature is now "//label//" (level "
        message=message//self%temperatureLevelCurrent//" of "//self%temperatureLevelCount//")"
        call Galacticus_Display_Message(message)
     end if
  end if
  ! Update the simulation state.
  call self%simulationState%update(stateVector,self%logging())
  return
end subroutine simulatorAnnealedDifferentialEvolutionUpdate

double precision function simulatorAnnealedDifferentialEvolutionTemperature(self)
  !% Return the temperature.
  implicit none
  class(simulatorAnnealedDifferentialEvolution), intent(inout) :: self

  simulatorAnnealedDifferentialEvolutionTemperature=self%temperatures(self%temperatureLevelCurrent)
  return
end function simulatorAnnealedDifferentialEvolutionTemperature

logical function simulatorAnnealedDifferentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,randomNumberGenerator)
  !% Return whether or not to accept a proposal.
  use Pseudo_Random
  use MPI_Utilities
  implicit none
  class           (simulatorAnnealedDifferentialEvolution), intent(inout) :: self
  double precision                                        , intent(in   ) :: logPosterior         , logPosteriorProposed
  type            (pseudoRandom                          ), intent(inout) :: randomNumberGenerator
  double precision                                                        :: x

  ! Decide whether to take step.
  x=randomNumberGenerator%sample(mpiRankOffset=.true.)
  simulatorAnnealedDifferentialEvolutionAcceptProposal=      &
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
end function simulatorAnnealedDifferentialEvolutionAcceptProposal
