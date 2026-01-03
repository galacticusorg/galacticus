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
Contains a module which implements a simple test of mapping a function over all components in a \gls{node}.
!!}

module Test_Nodes_Tasks
  !!{
  Implements a simple test of mapping a function over all components in a \gls{node}.
  !!}
  private
  public :: Test_Node_Task

  ! Records of which components have been seen in a test of function mapping.
  logical :: componentBasicStandardSeen    =.false.
  logical :: componentBlackHoleStandardSeen=.false.

contains

  subroutine Test_Node_Task(node)
    !!{
    Implements simple tests of mapping functions over all components in a \gls{node}.
    !!}
    use :: Display         , only : displayVerbositySet, verbosityLevelStandard
    use :: Galacticus_Nodes, only : nodeComponent      , nodeComponentBlackHole, treeNode
    use :: Unit_Tests      , only : Assert
    implicit none
    type            (treeNode     ), intent(inout) :: node
    procedure       (testVoidFunc ), pointer       :: myFuncVoid => testVoidFunc
    class           (nodeComponent), pointer       :: component

    ! Set verbosity level.
    call displayVerbositySet(verbosityLevelStandard)

    ! Create a black hole component.
    component => node%blackHole(autoCreate=.true.)

    ! Map a void function (subroutine) over all components.
    call node%mapVoid(myFuncVoid)
    call Assert('Map void function over all components',all([componentBasicStandardSeen,componentBlackHoleStandardSeen]),.true.)
    return
  end subroutine Test_Node_Task

  subroutine testVoidFunc(component)
    !!{
    A simple void function used in testing mapping over a function over all components.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponent
    use :: ISO_Varying_String, only : operator(==)
    implicit none
    class(nodeComponent), intent(inout) :: component

    if (component%type() == "nodeComponent:basic:standard"    ) componentBasicStandardSeen    =.true.
    if (component%type() == "nodeComponent:blackHole:standard") componentBlackHoleStandardSeen=.true.
    return
  end subroutine testVoidFunc

end module Test_Nodes_Tasks
