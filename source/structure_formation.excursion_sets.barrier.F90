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
Contains a module which provides a class that implements barriers for the excursion set problem.
!!}

module Excursion_Sets_Barriers
  !!{
  Provides a class that implements barriers for the excursion set problem.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>excursionSetBarrier</name>
   <descriptiveName>Excursion Set Barrier</descriptiveName>
   <description>
    Class providing barriers for the excursion set problem.
   </description>
   <default>criticalOverdensity</default>
   <method name="barrier" >
    <description>Return the barrier height at the given variance and time. The {\normalfont \ttfamily rateCompute} should be set to {\normalfont \ttfamily true} if the barrier is being used in a calculation of barrier crossing rates, and to {\normalfont \ttfamily false} otherwise.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: variance   , time</argument>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>logical                   , intent(in   ) :: rateCompute</argument>
   </method>
   <method name="barrierGradient" >
    <description>Return the gradient of the barrier with respect to variance at the given variance and time. The {\normalfont \ttfamily rateCompute} should be set to {\normalfont \ttfamily true} if the barrier is being used in a calculation of barrier crossing rates, and to {\normalfont \ttfamily false} otherwise.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: variance   , time</argument>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>logical                   , intent(in   ) :: rateCompute</argument>
   </method>
  </functionClass>
  !!]

  ! Enumeration for whether a remapping of the barrier should apply to cases where the excursion set is used for rate calculation, non-rate calculations, or both.
  !![
  <enumeration>
   <name>excursionSetRemap</name>
   <description>Specifies whether a remapping of the barrier should apply to cases where the excursion set is used for rate calculation, non-rate calculations, or both.</description>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <entry label="rates"   />
   <entry label="nonRates"/>
   <entry label="both"    />
  </enumeration>
  !!]

end module Excursion_Sets_Barriers
