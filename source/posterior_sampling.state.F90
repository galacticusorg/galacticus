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
   <description>Class providing state during posterior sampling simulations.</description>
   <default>simple</default>
   <data>integer :: parameterCount, stepCount, chainIndexValue</data>
   <method name="parameterCountSet" >
    <description>Set the number of parameters in this state.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: parameterCount</argument>
   </method>
   <method name="chainIndex" >
    <description>Return the index of the chain associated with this state.</description>
    <type>integer</type>
    <pass>yes</pass>
    <code>
     posteriorSampleStateChainIndex=self%chainIndexValue
    </code>
   </method>
   <method name="chainIndexSet" >
    <description>Set the index of the chain associated with this state.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>integer, intent(in   ) :: chainIndex</argument>
    <code>
     self%chainIndexValue=chainIndex
    </code>
   </method>
   <method name="count" >
    <description>Returns the number of steps in the current state.</description>
    <type>integer</type>
    <pass>yes</pass>
    <code>
     posteriorSampleStateCount=self%stepCount
    </code>
   </method>
   <method name="dimension" >
    <description>Returns the dimension of the state.</description>
    <type>integer</type>
    <pass>yes</pass>
    <code>
     posteriorSampleStateDimension=self%parameterCount
    </code>
   </method>
   <method name="reset" >
    <description>Reset the state object.</description>
    <type>void</type>
    <pass>yes</pass>
    <code>
     self%stepCount=0
    </code>
   </method>
   <method name="get" >
    <description>Get the current state vector.</description>
    <type>double precision, dimension(self%parameterCount)</type>
    <pass>yes</pass>
   </method>
   <method name="update" >
    <description>Update the state vector.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), dimension(:)           :: stateNew</argument>
    <argument>logical         , intent(in   )                         :: logState</argument>
    <argument>logical         , intent(in   )                         :: isConverged</argument>
    <argument>logical         , intent(in   ), dimension(:), optional :: outlierMask</argument>
   </method>
   <method name="mean" >
    <description>Return the mean state.</description>
    <type>double precision, dimension(self%parameterCount)</type>
    <pass>yes</pass>
   </method>
   <method name="variance" >
    <description>Return the variance in the state vector.</description>
    <type>double precision, dimension(self%parameterCount)</type>
    <pass>yes</pass>
   </method>
   <method name="acceptanceRate" >
    <description>Return the state acceptance rate.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="restore" >
    <description>Restore the state, one step at a time.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ), dimension(:) :: stateVector</argument>
    <argument>logical         , intent(in   )               :: first</argument>
   </method>
  </functionClass>
  !!]

end module Posterior_Sampling_State
