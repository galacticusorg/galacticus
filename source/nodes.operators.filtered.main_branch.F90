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
  Implements a node operator class that applies only to main branch nodes during tree initialization only.
  !!}

  !![
  <nodeOperator name="nodeOperatorFilteredMainBranch">
    <description>
      A node operator class that applies only to main branch nodes during tree initialization only. This uses a fast algorithm to
      determine main branch status, so is more efficient that using the \refClass{nodeOperatorFiltered} class along with a
      \refClass{galacticFilterMainBranch} filter.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorMulti) :: nodeOperatorFilteredMainBranch
     !!{
     A node operator class that applies only to main branch nodes during tree initialization only.
     !!}
     private
     logical            :: isTreeInitialization=.false., invertFilter
     integer(kind_int8) :: uniqueIDMainBranch          , uniqueIDMainBranchRoot
   contains
     procedure :: nodeTreeInitialize => filteredMainBranchNodeTreeInitialize
     procedure :: isActive           => filteredMainBranchIsActive
  end type nodeOperatorFilteredMainBranch

  interface nodeOperatorFilteredMainBranch
     !!{
     Constructors for the \refClass{nodeOperatorFilteredMainBranch} node operator class.
     !!}
     module procedure filteredMainBranchConstructorParameters
     module procedure filteredMainBranchConstructorInternal
  end interface nodeOperatorFilteredMainBranch

contains

  function filteredMainBranchConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorFilteredMainBranch} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorFilteredMainBranch)                :: self
    type   (inputParameters               ), intent(inout) :: parameters

    self%nodeOperatorMulti=nodeOperatorMulti(parameters)
    !![
    <inputParameter>
      <name>invertFilter</name>
      <variable>self%invertFilter</variable>
      <defaultValue>.false.</defaultValue>
      <description>If true, the filter is inverted to pass only nodes \emph{not} on the main branch.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParametersValidate source="parameters" multiParameters="nodeOperator"/>
    !!]
    self%uniqueIDMainBranchRoot=-huge(1_kind_int8)
    self%uniqueIDMainBranch    =-huge(1_kind_int8)
    return
  end function filteredMainBranchConstructorParameters

  function filteredMainBranchConstructorInternal(processes,invertFilter) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorFilteredMainBranch} node operator class.
    !!}
    implicit none
    type   (nodeOperatorFilteredMainBranch)                        :: self
    type   (multiProcessList              ), intent(in   ), target :: processes
    logical                                , intent(in   )         :: invertFilter

    self%nodeOperatorMulti     =nodeOperatorMulti(processes)
    self%invertFilter          =invertFilter
    self%uniqueIDMainBranchRoot=-huge(1_kind_int8)
    self%uniqueIDMainBranch    =-huge(1_kind_int8)
    return
  end function filteredMainBranchConstructorInternal

  subroutine filteredMainBranchNodeTreeInitialize(self,node)
    !!{
    Perform node tree initialization.
    !!}
    implicit none
    class(nodeOperatorFilteredMainBranch), intent(inout), target  :: self
    type (treeNode                      ), intent(inout), target  :: node
    type (multiProcessList              )               , pointer :: process_

    self%isTreeInitialization=.true.
    if (self%isActive(node)) then
       process_ => self%processes
       do while (associated(process_))
          call process_%process_%nodeTreeInitialize(node)
          process_ => process_%next
       end do
    end if
    self%isTreeInitialization=.false.
    return
  end subroutine filteredMainBranchNodeTreeInitialize

  logical function filteredMainBranchIsActive(self,node) result(isActive)
    !!{
    Return true if the given {\normalfont \ttfamily node} is on the main branch of the tree. Here we assume that a depth-first
    walk of the tree is being performed. As such, we keep a record of the (unique ID of the) last main branch node seen (starting
    from the tip of the main branch). A subsequent node is then only on the main branch if its first child is that same last seen
    main branch node.
    !!}
    implicit none
    class(nodeOperatorFilteredMainBranch), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    type (treeNode                      ), pointer       :: nodeMainBranch

    isActive=.false.
    if (.not.self%isTreeInitialization) return
    if (node%hostTree%nodeBase%uniqueID() /= self%uniqueIDMainBranchRoot) then
       nodeMainBranch => node%hostTree%nodeBase
       do while (associated(nodeMainBranch%firstChild))
          nodeMainBranch => nodeMainBranch%firstChild
       end do
       self%uniqueIDMainBranchRoot=node          %hostTree%nodeBase%uniqueID()
       self%uniqueIDMainBranch    =nodeMainBranch                  %uniqueID()
    end if

    if (associated(node%firstChild)) then
       isActive=node%firstChild%uniqueID() == self%uniqueIDMainBranch
       if (isActive) self%uniqueIDMainBranch=node%uniqueID()
    else
       isActive=node           %uniqueID() == self%uniqueIDMainBranch
    end if
    if (self%invertFilter) isActive=.not.isActive
    if (.not.associated(node%parent)) self%uniqueIDMainBranchRoot=-huge(0_kind_int8)
    return
  end function filteredMainBranchIsActive
