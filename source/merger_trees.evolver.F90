!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which provides a class that implements evolution of merger trees.

module Merger_Trees_Evolve
  !% Provides a class that implements evolution of merger trees.
  use    Galacticus_Nodes, only : mergerTree
  !$ use OMP_Lib         , only : omp_lock_kind
  
  !# <functionClass>
  !#  <name>mergerTreeEvolver</name>
  !#  <descriptiveName>Merger Tree Evolvers</descriptiveName>
  !#  <description>Class providing evolvers for merger trees.</description>
  !#  <default>standard</default>
  !#  <stateful>yes</stateful>
  !#  <method name="evolve" >
  !#   <description>Evolve a merger tree.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>type            (mergerTree   ), target  , intent(inout) :: tree                           </argument>
  !#   <argument>double precision                         , intent(in   ) :: timeEnd                        </argument>
  !#   <argument>logical                                  , intent(  out) :: treeDidEvolve     , suspendTree</argument>
  !#   <argument>logical                                  , intent(in   ) :: deadlockReporting              </argument>
  !#   <argument>integer         (omp_lock_kind), optional, intent(inout) :: initializationLock             </argument>
  !#  </method>
  !# </functionClass>
  
end module Merger_Trees_Evolve
