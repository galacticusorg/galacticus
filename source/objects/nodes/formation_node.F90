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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides access to the formation node of a node.
!!}

module Nodes_Formation_Node
  !!{RST
  Provides access to the formation node of a node.
  !!}
  private
  public :: nodeFormationRequired

contains

  function nodeFormationRequired(node,requiredBy) result(nodeFormation)
    !!{RST
    Return a pointer to the formation node of ``node``. Formation nodes exist only if a node operator creates them, so if there
    is none a fatal error is reported, naming the option (``requiredBy``) which required it and the node operator which would
    create it.
    !!}
    use :: Display         , only : displayGreen, displayReset
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type     (treeNode), pointer       :: nodeFormation
    type     (treeNode), intent(in   ) :: node
    character(len=*   ), intent(in   ) :: requiredBy

    nodeFormation => node%formationNode
    if (.not.associated(nodeFormation))                                                                             &
         & call Error_Report(                                                                                       &
         &                   'no formation node exists, but one is required by '//requiredBy//char(10)           // &
         &                   displayGreen()//'HELP:'//displayReset()                                             // &
         &                   ' formation nodes are created by, for example, the `nodeFormationTimeCole2000` node'// &
         &                   ' operator - add it to the list of node operators, or choose an option which does'  // &
         &                   ' not use the formation node'                                                       // &
         &                   {introspection:location}                                                               &
         &                  )
    return
  end function nodeFormationRequired

end module Nodes_Formation_Node
