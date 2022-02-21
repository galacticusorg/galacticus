!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a merger tree build controller class which provides no control.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerUncontrolled">
   <description>A merger tree build controller class which provides no control.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerUncontrolled
     !!{
     A merger tree build controller class which provides no control.
     !!}
     private
  contains
     procedure :: control => controlUncontrolled
  end type mergerTreeBuildControllerUncontrolled

  interface mergerTreeBuildControllerUncontrolled
     !!{
     Constructors for the ``uncontrolled'' merger tree build controller class.
     !!}
     module procedure uncontrolledConstructorParameters
  end interface mergerTreeBuildControllerUncontrolled

contains

  function uncontrolledConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``uncontrolled'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(mergerTreeBuildControllerUncontrolled)                :: self
    type(inputParameters                      ), intent(inout) :: parameters
    
    self=mergerTreeBuildControllerUncontrolled()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function uncontrolledConstructorParameters

  logical function controlUncontrolled(self,node,treeWalker_)
    !!{
    Apply a set of filters to a {\normalfont \ttfamily tree} combined with ``all'' operations.
    !!}
    implicit none
    class(mergerTreeBuildControllerUncontrolled), intent(inout)          :: self
    type (treeNode                             ), intent(inout), pointer :: node
    class(mergerTreeWalkerClass                ), intent(inout)          :: treeWalker_
    !$GLC attributes unused :: self, node, treeWalker_

    ! Always return true as we never want to halt tree building.
    controlUncontrolled=.true.
    return
  end function controlUncontrolled
