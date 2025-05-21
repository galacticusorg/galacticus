!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a node operator class that simply initializes the dark matter profile.
  !!}
  
  !![
  <nodeOperator name="nodeOperatorDarkMatterProfileInitialize">
   <description>
    A node operator class that simply initializes the dark matter profile.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfileInitialize
     !!{
     A node operator class that simply initializes the dark matter profile in halos.
     !!}
     private
   contains
     procedure :: nodeTreeInitialize => darkMatterProfileInitializeNodeTreeInitialize
  end type nodeOperatorDarkMatterProfileInitialize
  
  interface nodeOperatorDarkMatterProfileInitialize
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfileInitialize} node operator class.
     !!}
     module procedure darkMatterProfileInitializeConstructorParameters
  end interface nodeOperatorDarkMatterProfileInitialize
  
contains
  
  function darkMatterProfileInitializeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileInitialize} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfileInitialize)                :: self
    type (inputParameters                        ), intent(inout) :: parameters

    self=nodeOperatorDarkMatterProfileInitialize()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkMatterProfileInitializeConstructorParameters

  subroutine darkMatterProfileInitializeNodeTreeInitialize(self,node)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileInitialize), intent(inout), target  :: self
    type (treeNode                               ), intent(inout), target  :: node
    class(nodeComponentDarkMatterProfile         )               , pointer :: darkMatterProfile

    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    return
  end subroutine darkMatterProfileInitializeNodeTreeInitialize
