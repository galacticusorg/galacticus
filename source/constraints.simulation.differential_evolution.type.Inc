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


  type, public, extends(simulator) :: simulatorDifferentialEvolution
     !% Defines a simulator class implementing the differential evolution algorithm.
     private
     integer                                          :: parameterCount       , stepsMaximum          , &
          &                                              stepsPostConvergence , acceptanceAverageCount, &
          &                                              stateSwapCount
     type   (prior           ), pointer, dimension(:) :: parameterPriors
     type   (distributionList), pointer, dimension(:) :: randomDistributions
     class  (likelihood      ), pointer               :: modelLikelihood
     class  (convergence     ), pointer               :: simulationConvergence
     class  (state           ), pointer               :: simulationState
     class  (deProposalSize  ), pointer               :: proposalSize
     class  (deRandomJump    ), pointer               :: randomJump
     type   (varying_string  )                        :: logFileRoot
   contains
     !@ <objectMethods>
     !@   <object>simulatorDifferentialEvolution</object>
     !@   <objectMethod>
     !@     <method>logging</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if the simulator is currently logging state.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>posterior</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(state)\textgreater} simulationState/argin</arguments>
     !@     <description>Return the log of posterior probability for the given {\tt simulationState}.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>update</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ stateVector\argin</arguments>
     !@     <description>Update the simulator to the new {\tt stateVector} after a step.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: simulate       => simulatorDifferentialEvolutionSimulate
     procedure :: logging        => simulatorDifferentialEvolutionLogging
     procedure :: posterior      => simulatorDifferentialEvolutionPosterior
     procedure :: update         => simulatorDifferentialEvolutionUpdate
     procedure :: stepSize       => simulatorDifferentialEvolutionStepSize
     procedure :: acceptProposal => simulatorDifferentialEvolutionAcceptProposal
  end type simulatorDifferentialEvolution
  
  interface simulatorDifferentialEvolution
     !% Constructors for the ``differentialEvolution'' simulator class.
     module procedure simulatorDifferentialEvolutionConstructor
  end interface simulatorDifferentialEvolution
