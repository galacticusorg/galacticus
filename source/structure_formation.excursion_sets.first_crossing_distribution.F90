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
Contains a module which provides a class for first crossing distributions for excursion set calculations.
!!}


module Excursion_Sets_First_Crossings
  !!{
  Provides a class for first crossing distributions for excursion set calculations.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>excursionSetFirstCrossing</name>
   <descriptiveName>Excursion Set First Crossing Statistics</descriptiveName>
   <description>
    Class providing first crossing statistics for the excursion set problem.
   </description>
   <default>linearBarrier</default>
   <data>logical :: coordinatedMPI_=.false. </data>
   <method name="probability" >
    <description>Return the probability for a trajectory to make its first crossing of the barrier at the given {\normalfont \ttfamily variance} and {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: variance, time</argument>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="rate" >
    <description>Return the rate of first crossing for excursion sets beginning at the given {\normalfont \ttfamily variance} and {\normalfont \ttfamily time} to transition to a first crossing at the given {\normalfont \ttfamily varianceProgenitor}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: variance, varianceProgenitor, time</argument>
    <argument>type            (treeNode), intent(inout) :: node</argument>
   </method>
   <method name="rateNonCrossing" >
    <description>Return the rate of non-crossing for excursion sets beginning at the given {\normalfont \ttfamily variance} and {\normalfont \ttfamily time}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ) :: variance, massMinimum, time</argument>
    <argument>type            (treeNode), intent(inout) :: node                       </argument>
   </method>
   <method name="coordinatedMPI" >
    <description>Sets the state of coordination under MPI. If set to true then the object can assume that any calculations it performs are being performed identically by all other MPI processes. This permits, for example, coordinated tabulation of results across MPI processes.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>logical, intent(in   ) :: state</argument>
    <code>
      self%coordinatedMPI_=state
    </code>
   </method>
  </functionClass>
  !!]

end module Excursion_Sets_First_Crossings
