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
Contains a module which implements a stopping criteria class for constraint simulations.
!!}

module Posterior_Sampling_Stopping_Criteria
  !!{
  Implements a stopping criteria class for constraint simulations.
  !!}
  use :: Posterior_Sampling_State, only : posteriorSampleStateClass
  private

  !![
  <functionClass>
   <name>posteriorSampleStoppingCriterion</name>
   <descriptiveName>Posterior Sampling Stopping Criteria</descriptiveName>
   <description>Class providing stopping criteria for posterior sampling simulations.</description>
   <default>never</default>
   <method name="stop" >
    <description>Returns true if the posterior sampling should stop.</description>
    <type>logical</type>
    <pass>yes</pass>
    <argument>class(posteriorSampleStateClass), intent(inout) :: simulationState</argument>
   </method>
  </functionClass>
  !!]

end module Posterior_Sampling_Stopping_Criteria
