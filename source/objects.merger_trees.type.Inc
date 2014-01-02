!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Defines an object class for merger trees.

  type, public :: mergerTree
     !% The merger tree object type.
     integer         (kind=kind_int8)          :: index
     type            (hdf5Object    )          :: hdf5Group
     double precision                          :: volumeWeight
     logical                                   :: initialized
     type            (treeNode      ), pointer :: baseNode    =>null()
     type            (mergerTree    ), pointer :: nextTree    =>null()
   contains
     ! Tree creation/destruction.
     !@ <objectMethods>
     !@   <object>mergerTree</object>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <description>Destroys the merger tree, including all nodes and their components.</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroyBranch</method>
     !@     <description>
     !@       Destroys a branch of a merger tree starting from the supplied node. All nodes and their components are destroyed.
     !@     </description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless *type(treeNode)\textgreater} thisNode\arginout</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>getNode</method>
     !@     <description>Returns a pointer to the node with given index in the merger tree, or a null pointer if no such node exists.</description>
     !@     <type>\textcolor{red}{\textless *type(treeNode)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless integer(kind\_int8)\textgreater} nodeIndex\argin</arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: destroy      =>Merger_Tree_Destroy
     procedure :: destroyBranch=>Merger_Tree_Destroy_Branch
     procedure :: getNode      =>Merger_Tree_Node_Get
  end type mergerTree

  ! Flag to indicate if the tree node create routines have been initialized.
  logical :: treeNodeCreateInitialized=.false.
