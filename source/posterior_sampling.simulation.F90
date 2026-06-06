!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Contains a module which implements a class of posterior sampling simulators.
!!}

module Posterior_Sampling_Simulation
  !!{
  Implements a class of posterior sampling simulators.
  !!}
  private

  !![
  <functionClass>
   <name>posteriorSampleSimulation</name>
   <descriptiveName>Posterior Sampling Simulations</descriptiveName>
   <description>Class providing simulators for Bayesian posterior sampling---algorithms that draw samples from
    the joint posterior distribution of model parameters given observational constraints. Implementations
    orchestrate the Markov Chain Monte Carlo or other stochastic exploration strategy, managing the
    collection of walkers or chains, proposing new parameter vectors, evaluating the likelihood and prior,
    and deciding whether to accept or reject proposed moves. The default implementation uses the
    Differential Evolution MCMC algorithm for efficient high-dimensional sampling.</description>
   <default>differentialEvolution</default>
   <method name="simulate" >
    <description>Execute the full posterior sampling run, iterating the proposal, acceptance, and logging cycle until either the stopping criterion is met or the maximum step count is reached.</description>
    <type>void</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

end module Posterior_Sampling_Simulation
