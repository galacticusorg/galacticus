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
Contains a module which implements a class that provides convergence criteria for posterior sampling simulations.
!!}

module Posterior_Sampling_Convergence
  !!{
  Implements a class that provides convergence criteria for posterior sampling simulations.
  !!}
  use :: Posterior_Sampling_State, only : posteriorSampleStateClass
  private

  !![
  <functionClass>
   <name>posteriorSampleConvergence</name>
   <descriptiveName>Posterior Sampling Convergence Criteria</descriptiveName>
   <description>Class providing convergence criteria for Bayesian posterior sampling simulations---diagnostics
    that assess whether the Markov chains have adequately explored the posterior distribution and can
    be declared converged. Methods return whether convergence has been reached (and at which step),
    reset the convergence calculation when chains are restarted, log convergence diagnostics to file,
    and identify outlier chains. The default implementation uses the Gelman-Rubin $\hat{R}$ statistic,
    which compares within-chain and between-chain variances across all active walkers.</description>
   <default>gelmanRubin</default>
   <method name="isConverged" >
    <description>Returns true if the posterior sampling is deemed to be converged.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>class           (posteriorSampleStateClass), intent(inout), optional :: simulationState</argument>
    <argument>double precision                           , intent(in   ), optional :: logLikelihood</argument>
   </method>
   <method name="convergedAtStep" >
    <description>Returns the simulation step index at which the posterior sampling was first declared converged, or a negative value if convergence has not yet been reached.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="reset" >
    <description>Reset the convergence diagnostic to its initial state, clearing any accumulated chain statistics, typically called when the sampler is restarted or the chain population is reinitialized.</description>
    <type>void</type>
    <pass>yes</pass>
   </method>
   <method name="logReport" >
    <description>Write a human-readable report of the current convergence state (e.g.\ the $\hat{R}$ statistic and per-parameter values) to the open file unit \mono{fileUnit}.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: fileUnit</argument>
   </method>
   <method name="stateIsOutlier" >
    <description>Return true if the specified state is deemed to be an outlier.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: stateIndex</argument>
   </method>
  </functionClass>
  !!]

end module Posterior_Sampling_Convergence
