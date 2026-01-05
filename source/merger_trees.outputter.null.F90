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
  Implements a merger trees outputter class which does no output.
  !!}

  !![
  <mergerTreeOutputter name="mergerTreeOutputterNull">
   <description>A merger tree outputter which does no output.</description>
  </mergerTreeOutputter>
  !!]
  type, extends(mergerTreeOutputterClass) :: mergerTreeOutputterNull
     !!{
     Implementation of a merger tree outputter which does no output.
     !!}
     private
   contains
     procedure :: outputTree => nullOutputTree
     procedure :: outputNode => nullOutputNode
     procedure :: finalize   => nullFinalize
  end type mergerTreeOutputterNull

  interface mergerTreeOutputterNull
     !!{
     Constructors for the \refClass{mergerTreeOutputterNull} merger tree outputter.
     !!}
     module procedure nullConstructorParameters
  end interface mergerTreeOutputterNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOutputterNull} merger tree outputter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeOutputterNull)                :: self
    type(inputParameters        ), intent(inout) :: parameters

    self=mergerTreeOutputterNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  subroutine nullOutputTree(self,tree,indexOutput,time)
    !!{
    Perform no output.
    !!}
    implicit none
    class           (mergerTreeOutputterNull), intent(inout)         :: self
    type            (mergerTree             ), intent(inout), target :: tree
    integer         (c_size_t               ), intent(in   )         :: indexOutput
    double precision                         , intent(in   )         :: time
    !$GLC attributes unused :: self, tree, indexOutput, time

    return
  end subroutine nullOutputTree

  subroutine nullOutputNode(self,node,indexOutput)
    !!{
    Perform no output.
    !!}
    implicit none
    class           (mergerTreeOutputterNull), intent(inout) :: self
    type            (treeNode               ), intent(inout) :: node
    integer         (c_size_t               ), intent(in   ) :: indexOutput
    !$GLC attributes unused :: self, node, indexOutput

    return
  end subroutine nullOutputNode

  subroutine nullFinalize(self)
    !!{
    Finalize merger tree output.
    !!}
    implicit none
    class(mergerTreeOutputterNull), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine nullFinalize
