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
Contains a module which implements a class that maintains during posterior sampling.
!!}

module Posterior_Sampling_State
  !!{
  Implements a class that maintains during posterior sampling.
  !!}
  private

  !![
  <functionClass>
   <name>posteriorSampleState</name>
   <descriptiveName>Posterior Sampling State</descriptiveName>
   <description>Class providing the state vector during Bayesian posterior sampling simulations---the current
    position of a Markov chain in parameter space, together with its chain index, step count, and parameter
    count. Implementations store and retrieve the parameter vector (via \mono{get}/\mono{set}), advance
    the step counter, and optionally maintain a history of chain positions for convergence diagnostics.
    The state is updated at each MCMC step and queried by the likelihood function to evaluate the
    model at the proposed parameter values.</description>
   <default>simple</default>
   <method name="parameterCountSet" >
    <description>Set the number of active parameters that this state vector will track, allocating internal storage for the parameter vector of the given dimension.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: parameterCount</argument>
   </method>
   <method name="chainIndex" >
    <description>Return the integer index (0-based) of the Markov chain that owns this state object, used to identify chains in multi-chain posterior sampling algorithms.</description>
    <type>integer</type>
    <pass>yes</pass>
    <code>
     posteriorSampleStateChainIndex=self%chainIndexValue
    </code>
   </method>
   <method name="chainIndexSet" >
    <description>Assign the integer chain index to this state object, identifying which Markov chain it belongs to in a multi-chain ensemble sampler.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: chainIndex</argument>
    <code>
     self%chainIndexValue=chainIndex
    </code>
   </method>
   <method name="count" >
    <description>Returns the total number of sampling steps that have been taken since the state was last reset, used for logging and convergence diagnostics.</description>
    <type>integer</type>
    <pass>yes</pass>
    <code>
     posteriorSampleStateCount=self%stepCount
    </code>
   </method>
   <method name="dimension" >
    <description>Returns the number of active parameters (dimension of the state vector) for this sampling state, equal to the value previously set by \mono{parameterCountSet}.</description>
    <type>integer</type>
    <pass>yes</pass>
    <code>
     posteriorSampleStateDimension=self%parameterCount
    </code>
   </method>
   <method name="reset" >
    <description>Reset the state object to its initial condition by zeroing the step counter and clearing any accumulated history, in preparation for a new sampling run.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
     self%stepCount=0
    </code>
   </method>
   <method name="get" >
    <description>Return the current parameter vector representing the position of this chain in the model parameter space at the most recently accepted sampling step.</description>
    <type>double precision, dimension(self%parameterCount)</type>
    <pass>yes</pass>
   </method>
   <method name="update" >
    <description>Advance the state to the new parameter vector \mono{stateNew}, optionally logging it to the chain history; \mono{isConverged} and \mono{outlierMask} are used to track whether post-convergence steps should be recorded.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), dimension(:)           :: stateNew</argument>
    <argument>logical         , intent(in   )                         :: logState</argument>
    <argument>logical         , intent(in   )                         :: isConverged</argument>
    <argument>logical         , intent(in   ), dimension(:), optional :: outlierMask</argument>
   </method>
   <method name="mean" >
    <description>Return the mean parameter vector computed over all stored steps in the chain history, providing a point estimate of the posterior mode or mean.</description>
    <type>double precision, dimension(self%parameterCount)</type>
    <pass>yes</pass>
   </method>
   <method name="variance" >
    <description>Return the per-parameter variance computed over all stored steps in the chain history, providing a measure of posterior width for each model parameter.</description>
    <type>double precision, dimension(self%parameterCount)</type>
    <pass>yes</pass>
   </method>
   <method name="acceptanceRate" >
    <description>Return the fraction of proposed moves that were accepted over the recent history of this chain, used to monitor and adaptively tune proposal distributions.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="restore" >
    <description>Replay a previously recorded state vector into the state history one step at a time, used when resuming a sampling run from a saved log file; \mono{first} signals the start of the restoration sequence.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), dimension(:) :: stateVector</argument>
    <argument>logical         , intent(in   )               :: first</argument>
   </method>
   <data>integer :: parameterCount, stepCount, chainIndexValue</data>
  </functionClass>
  !!]

end module Posterior_Sampling_State
