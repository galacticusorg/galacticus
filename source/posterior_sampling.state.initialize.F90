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
Contains a module which implements a class for posterior sampling state initialization.
!!}

module Posterior_Sampling_State_Initialize
  !!{
  Implements a class for posterior sampling state initialization.
  !!}
  use :: Model_Parameters        , only : modelParameterList
  use :: Models_Likelihoods      , only : posteriorSampleLikelihoodClass
  use :: Posterior_Sampling_State, only : posteriorSampleStateClass
  private

  !![
  <functionClass>
   <name>posteriorSampleStateInitialize</name>
   <descriptiveName>Posterior Sampling State Initialization</descriptiveName>
   <description>Class providing state initialization for Bayesian posterior sampling simulations---the strategy
    for choosing the starting point of each Markov chain before the main sampling loop begins. Implementations
    draw initial parameter vectors from the prior (random prior sampling), from a multivariate Gaussian
    centered on a maximum-likelihood estimate, or from a user-supplied restart file. Proper initialization
    helps the chains reach the high-probability region quickly and avoids premature convergence or
    long burn-in periods.</description>
   <default>priorRandom</default>
   <method name="initialize" >
    <description>Set the initial parameter vector of the simulation state before the main sampling loop begins, also returning the initial log-likelihood, log-posterior, and elapsed evaluation time for the starting point.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>class           (posteriorSampleStateClass     ), intent(inout)               :: simulationState</argument>
    <argument>type            (modelParameterList            ), intent(inout), dimension(:) :: modelParameters_</argument>
    <argument>class           (posteriorSampleLikelihoodClass), intent(inout)               :: modelLikelihood</argument>
    <argument>double precision                                , intent(  out)               :: timeEvaluatePrevious, logLikelihood, logPosterior</argument>
   </method>
  </functionClass>
  !!]

end module Posterior_Sampling_State_Initialize
