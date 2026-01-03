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
Contains a module which implements class for the random jump component in differential evolution algorithms.
!!}

module Posterior_Sample_Differential_Random_Jump
  !!{
  Implements a class for the random jump component in differential evolution algorithms.
  !!}
  use :: Model_Parameters        , only : modelParameterList
  use :: Posterior_Sampling_State, only : posteriorSampleStateClass
  private

  !![
  <functionClass>
   <name>posteriorSampleDffrntlEvltnRandomJump</name>
   <descriptiveName>Posterior Sampling Differential Evolution Random Jumps</descriptiveName>
   <description>Class providing random jumps to be added to proposals in differential evolution posterior samplers.</description>
   <method name="sample" >
    <description>Sample from the jump distribution.</description>
    <type>double precision, dimension(size(modelParameters_))</type>
    <pass>yes</pass>
    <argument>type (modelParameterList       ), dimension(:), intent(in   ) :: modelParameters_</argument>
    <argument>class(posteriorSampleStateClass)              , intent(inout) :: simulationState</argument>
   </method>
  </functionClass>
  !!]

end module Posterior_Sample_Differential_Random_Jump
