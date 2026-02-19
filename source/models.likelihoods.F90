!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which implements a likelihood class for posterior sampling simulations.
!!}

module Models_Likelihoods
  !!{
  Implements a likelihood class for posterior sampling simulations.
  !!}
  use :: Model_Parameters              , only : modelParameterList
  use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
  use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
  private

  !![
  <functionClass>
   <name>posteriorSampleLikelihood</name>
   <descriptiveName>Posterior Sampling Likelihoods</descriptiveName>
   <description>Class providing likelihoods for posterior sampling simulations.</description>
   <method name="evaluate" >
    <description>Evaluate the likelihood.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>class           (posteriorSampleStateClass      ), intent(inout)               :: simulationState</argument>
    <argument>type            (modelParameterList             ), intent(inout), dimension(:) :: modelParametersActive_, modelParametersInactive_</argument>
    <argument>class           (posteriorSampleConvergenceClass), intent(inout)               :: simulationConvergence</argument>
    <argument>double precision                                 , intent(in   )               :: temperature, logLikelihoodCurrent, logPriorCurrent, logPriorProposed</argument>
    <argument>real                                             , intent(inout)               :: timeEvaluate</argument>
    <argument>double precision                                 , intent(  out), optional     :: logLikelihoodVariance</argument>
    <argument>logical                                          , intent(inout), optional     :: forceAcceptance</argument>
   </method>
   <method name="willEvaluate" >
    <description>Returns true if the likelihood will be evaluated for this state.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>class           (posteriorSampleStateClass      ), intent(inout)               :: simulationState</argument>
    <argument>type            (modelParameterList             ), intent(in   ), dimension(:) :: modelParameters_</argument>
    <argument>class           (posteriorSampleConvergenceClass), intent(inout)               :: simulationConvergence</argument>
    <argument>double precision                                 , intent(in   )               :: temperature, logLikelihoodCurrent, logPriorCurrent, logPriorProposed</argument>
    <code>
     !$GLC attributes unused :: self, simulationState, modelParameters_, simulationConvergence, temperature, logLikelihoodCurrent, logPriorCurrent, logPriorProposed
     posteriorSampleLikelihoodWillEvaluate=.true.
    </code>
   </method>
   <method name="functionChanged" >
    <description>Respond to possible changes in the likelihood function.</description>
    <type>void</type>
    <pass>yes</pass>
   </method>
   <method name="restore" >
    <description>Log a report on convergence state to the given file unit.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), dimension(:) :: simulationState</argument>
    <argument>double precision, intent(in   )               :: logLikelihood</argument>
    <code>
     !$GLC attributes unused :: self, simulationState, logLikelihood
    </code>
   </method>
  </functionClass>
  !!]

end module Models_Likelihoods
