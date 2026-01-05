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
Contains a module which handles resetting of calculations before a new or updated node is processed.
!!}

module Calculations_Resets
  !!{
  Handles resetting of calculations before a new or updated node is processed.
  !!}
  implicit none
  private
  public :: Calculations_Reset

contains

  !![
  <functionGlobal>
   <unitName>Calculations_Reset</unitName>
   <type>void</type>
   <module>Galacticus_Nodes, only : treeNode</module>
   <arguments>type(treeNode) , intent(inout) :: node</arguments>
  </functionGlobal>
  !!]
  subroutine Calculations_Reset(node)
    !!{
    Calls any routines required to reset all calculation for a new or updated node.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    !![
    <include directive="calculationResetTask" type="moduleUse">
    !!]
    include 'calculation_reset.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    type(treeNode), intent(inout) :: node

    !![
    <include directive="calculationResetTask" type="functionCall" functionType="void">
     <functionArgs>node,node%uniqueID()</functionArgs>
    !!]
    include 'calculation_reset.tasks.inc'
    !![
    </include>
    !!]

    !![
    <eventHook name="calculationReset">
     <callWith>node,node%uniqueID()</callWith>
    </eventHook>
    !!]
    return
  end subroutine Calculations_Reset

end module Calculations_Resets
