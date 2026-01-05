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
  Implements a node operator class that keeps the ``\gls{dmou}'' fixed between node promotion events.
  !!}

  !![
  <nodeOperator name="nodeOperatorDMOUninterpolated">
   <description>
    A node operator class that keeps the ``\gls{dmou}'' fixed between node promotion events.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDMOUninterpolated
     !!{
     A node operator class that keeps the ``\gls{dmou}'' fixed between node promotion events.
     !!}
     private
   contains
     procedure :: nodePromote => dmoUninterpolatedNodePromote
  end type nodeOperatorDMOUninterpolated
  
  interface nodeOperatorDMOUninterpolated
     !!{
     Constructors for the \refClass{nodeOperatorDMOUninterpolated} node operator class.
     !!}
     module procedure dmoUninterpolatedConstructorParameters
  end interface nodeOperatorDMOUninterpolated
  
contains
  
  function dmoUninterpolatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDMOUninterpolated} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDMOUninterpolated)                :: self
    type (inputParameters              ), intent(inout) :: parameters
 
    self=nodeOperatorDMOUninterpolated()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function dmoUninterpolatedConstructorParameters

  subroutine dmoUninterpolatedNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class    (nodeOperatorDMOUninterpolated), intent(inout) :: self
    type     (treeNode                     ), intent(inout) :: node
    type     (treeNode                     ), pointer       :: nodeParent
    class    (nodeComponentBasic           ), pointer       :: basicParent, basic
    !$GLC attributes unused :: self

    nodeParent  => node      %parent
    basic       => node      %basic ()
    basicParent => nodeParent%basic ()
    ! Ensure the two halos exist at the same time.
    if (basic%time() /= basicParent%time()) call Error_Report('node has not been evolved to its parent'//{introspection:location})
    ! Adjust the mass to that of the parent node.
    call basic%massSet(basicParent%mass())
    return
  end subroutine dmoUninterpolatedNodePromote
