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


  type, public, extends(simulatorDifferentialEvolution) :: simulatorTemperedDifferentialEvolution
     !% Defines a simulator class implementing the tempered differential evolution algorithm.
     private
     integer                                                        :: untemperedStepCount            , temperingLevelCount     , stepsPerLevel
     integer                                                        :: temperingStep                  , temperingLevelMonotonic
     double precision                                               :: temperatureMaximum
     class           (dePropSizeTempExp), pointer                   :: proposalSizeTemperatureExponent
     double precision                   , allocatable, dimension(:) :: temperatures
     class           (state            ), allocatable, dimension(:) :: temperedStates
   contains
     !@ <objectMethods>
     !@   <object>simulatorTemperedDifferentialEvolution</object>
     !@   <objectMethod>
     !@     <method>level</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return the current tempering level.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>temperature</method>
     !@     <type>\doublezero</type>
     !@     <arguments></arguments>
     !@     <description>Return the current temperature.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: logging        => simulatorTemperedDifferentialEvolutionLogging
     procedure :: acceptProposal => simulatorTemperedDifferentialEvolutionAcceptProposal
     procedure :: update         => simulatorTemperedDifferentialEvolutionUpdate
     procedure :: stepSize       => simulatorTemperedDifferentialEvolutionStepSize
     procedure :: level          => simulatorTemperedDifferentialEvolutionLevel
     procedure :: temperature    => simulatorTemperedDifferentialEvolutionTemperature
  end type simulatorTemperedDifferentialEvolution
  
  interface simulatorTemperedDifferentialEvolution
     !% Constructors for the ``temperedDifferentialEvolution'' simulator class.
     module procedure simulatorTemperedDifferentialEvolutionConstructor
  end interface simulatorTemperedDifferentialEvolution
