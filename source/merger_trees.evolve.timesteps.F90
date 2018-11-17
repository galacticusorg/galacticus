!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a class for merger tree evolution timestepping.

module Merger_Tree_Timesteps
  !% Implements a class for merger tree evolution timestepping.
  use Galacticus_Nodes
  implicit none
  private
  public :: timestepTask

  !# <functionClass>
  !#  <name>mergerTreeEvolveTimestep</name>
  !#  <descriptiveName>Merger Tree Evolution Timesteps</descriptiveName>
  !#  <description>Class providing timestep control for merger tree evolution.</description>
  !#  <default>standard</default>
  !#  <method name="timeEvolveTo" >
  !#   <description>Return the time to which the node can be evolved.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>type     (treeNode      ), intent(inout), target  :: node    </argument>
  !#   <argument>procedure(timestepTask  ), intent(  out), pointer :: task    </argument>
  !#   <argument>class    (*             ), intent(  out), pointer :: taskSelf</argument>
  !#   <argument>logical                  , intent(in   )          :: report  </argument>
  !#   <argument>type     (treeNode      ), intent(  out), pointer :: lockNode</argument>
  !#   <argument>type     (varying_string), intent(  out)          :: lockType</argument>
  !#  </method>
  !# </functionClass>
  
  abstract interface
     subroutine timestepTask(self,tree,node,deadlockStatus)
       import mergerTree, treeNode
       class  (*         ), intent(inout)          :: self
       type   (mergerTree), intent(in   )          :: tree
       type   (treeNode  ), intent(inout), pointer :: node
       integer            , intent(inout)          :: deadlockStatus
     end subroutine timestepTask
  end interface

end module Merger_Tree_Timesteps
