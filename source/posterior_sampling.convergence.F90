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
   <description>Class providing convergence criteria for posterior sampling simulations.</description>
   <default>gelmanRubin</default>
   <method name="isConverged" >
    <description>Returns true if the posterior sampling is deemed to be converged.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>class           (posteriorSampleStateClass), intent(inout), optional :: simulationState</argument>
    <argument>double precision                           , intent(in   ), optional :: logLikelihood</argument>
   </method>
   <method name="convergedAtStep" >
    <description>Returns the step at which the sampling reached convergence.</description>
    <type>integer</type>
    <pass>yes</pass>
   </method>
   <method name="reset" >
    <description>Reset the convergence calculation.</description>
    <type>void</type>
    <pass>yes</pass>
   </method>
   <method name="logReport" >
    <description>Log a report on convergence state to the given file unit.</description>
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
