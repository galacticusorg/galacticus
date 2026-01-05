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
Contains a module that implements a class that returns lists of states at which to sample the posterior.
!!}

module Posterior_Sampling_State_Samples
  !!{
  Implements a class that returns lists of states at which to sample the posterior.
  !!}
  use :: Model_Parameters        , only : modelParameterList
  use :: Posterior_Sampling_State, only : posteriorSampleStateSimple
  private

  !![
  <functionClass>
   <name>posteriorSamples</name>
   <descriptiveName>Posterior Samples</descriptiveName>
   <description>Class providing lists of states at which to sample the posterior.</description>
   <default>priorGrid</default>
   <method name="samples" >
    <description>Return a list of states at which to sample the posterior.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type(posteriorSampleStateSimple), intent(inout), dimension(:), allocatable :: simulationStates</argument>
    <argument>type(modelParameterList        ), intent(inout), dimension(:)              :: modelParameters_</argument>
   </method>
  </functionClass>
  !!]

end module Posterior_Sampling_State_Samples
