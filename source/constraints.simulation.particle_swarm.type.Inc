!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  type, public, extends(simulator) :: simulatorParticleSwarm
     !% Defines a simulator class implementing the particle swarm algorithm.
     private
     integer                                                    :: parameterCount                 , stepsMaximum                 , &
          &                                                        reportCount                    , logFlushCount
     double precision                                           :: accelerationCoefficientPersonal, accelerationCoefficientGlobal, &
          &                                                        inertiaWeight                  , velocityCoefficient
     logical                                                    :: isInteractive                  , resume
     type            (varying_string   )                        :: logFileRoot                    , interactionRoot              , &
          &                                                        logFilePreviousRoot
     type            (prior            ), pointer, dimension(:) :: parameterPriors
     type            (mappingList      ), pointer, dimension(:) :: parameterMappings
     class           (likelihood       ), pointer               :: modelLikelihood
     class           (convergence      ), pointer               :: simulationConvergence
     class           (stoppingCriterion), pointer               :: simulationStoppingCriterion
     class           (state            ), pointer               :: simulationState
     class           (stateInitializor ), pointer               :: simulationStateInitializor
   contains
     !@ <objectMethods>
     !@   <object>simulatorParticleSwarm</object>
     !@   <objectMethod>
     !@     <method>posterior</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\textcolor{red}{\textless class(state)\textgreater} simulationState\argin</arguments>
     !@     <description>Return the log of posterior probability for the given {\normalfont \ttfamily simulationState}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: simulate  => simulatorParticleSwarmSimulate
     procedure :: posterior => simulatorParticleSwarmPosterior
  end type simulatorParticleSwarm
  
  interface simulatorParticleSwarm
     !% Constructors for the ``particleSwarm'' simulator class.
     module procedure simulatorParticleSwarmConstructor
  end interface simulatorParticleSwarm
