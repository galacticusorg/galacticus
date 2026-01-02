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
Contains a module which provides a class that implements profiling of merger tree evolution.
!!}

module Merger_Tree_Evolve_Profilers
  !!{
  Provides a class that implements profiling of merger tree evolution.
  !!}
  use, intrinsic :: ISO_C_Binding   , only : c_size_t
  use            :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>mergerTreeEvolveProfiler</name>
   <descriptiveName>Merger Tree Evolver Profiler</descriptiveName>
   <description>Class providing profilers for merger tree evolution.</description>
   <default>null</default>
   <method name="stepDescriptor" >
    <description>Provide a descriptor of the current step.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (varying_string), intent(in   )               :: descriptor                                                  </argument>
   </method>
   <method name="profile" >
    <description>Profile a differential evolution step.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>type            (treeNode      ), intent(in   )               :: node                                                        </argument>
    <argument>double precision                , intent(in   )               :: time            , timeStart   , timeEnd      , timeStep     </argument>
    <argument>integer         (c_size_t      ), intent(in   )               :: countEvaluations                                            </argument>
    <argument>logical                         , intent(in   )               :: interrupted                                                 </argument>
    <argument>integer         (c_size_t      ), intent(in   )               :: propertyIndex                                               </argument>
    <argument>type            (varying_string), intent(in   )               :: propertyName                                                </argument>
    <argument>double precision                , intent(in   ), dimension(:) :: propertyValue   , propertyRate, propertyScale, propertyError</argument>
    <argument>double precision                , intent(in   )               :: timeCPU                                                     </argument>
   </method>
  </functionClass>
  !!]
  
end module Merger_Tree_Evolve_Profilers
