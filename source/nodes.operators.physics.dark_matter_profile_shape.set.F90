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
  Implements a node operator class that sets dark matter profile shape parameter.
  !!}

  use :: Dark_Matter_Profiles_Shape, only : darkMatterProfileShapeClass
  
  !![
  <nodeOperator name="nodeOperatorDarkMatterProfileShapeSet">
   <description>
    A node operator class that sets dark matter profile shape parameter.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfileShapeSet
     !!{
     A node operator class that set the dark matter profile shape parameter.
     !!}
     private
     class(darkMatterProfileShapeClass), pointer :: darkMatterProfileShape_ => null()
   contains
     final     ::                       darkMatterProfileShapeSetConstructorDestructor
     procedure :: nodeTreeInitialize => darkMatterProfileShapeSetNodeTreeInitialize
     procedure :: nodePromote        => darkMatterProfileShapeSetNodePromote
  end type nodeOperatorDarkMatterProfileShapeSet
  
  interface nodeOperatorDarkMatterProfileShapeSet
     !!{
     Constructors for the \refClass{nodeOperatorDarkMatterProfileShapeSet} node operator class.
     !!}
     module procedure darkMatterProfileShapeSetConstructorParameters
     module procedure darkMatterProfileShapeSetConstructorInternal
  end interface nodeOperatorDarkMatterProfileShapeSet
  
contains
  
  function darkMatterProfileShapeSetConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileShapeSet} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfileShapeSet)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(darkMatterProfileShapeClass          ), pointer       :: darkMatterProfileShape_

    !![
    <objectBuilder class="darkMatterProfileShape" name="darkMatterProfileShape_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfileShapeSet(darkMatterProfileShape_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileShape_"/>
    !!]
    return
  end function darkMatterProfileShapeSetConstructorParameters

  function darkMatterProfileShapeSetConstructorInternal(darkMatterProfileShape_) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkMatterProfileShapeSet} node operator class which takes a parameter set as input.
    !!}
    implicit none
    type (nodeOperatorDarkMatterProfileShapeSet)                        :: self
    class(darkMatterProfileShapeClass          ), intent(in   ), target :: darkMatterProfileShape_
    !![
    <constructorAssign variables="*darkMatterProfileShape_"/>
    !!]

    return
  end function darkMatterProfileShapeSetConstructorInternal

  subroutine darkMatterProfileShapeSetConstructorDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorDarkMatterProfileShapeSet} node operator class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfileShapeSet), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileShape_"/>
    !!]
    return
  end subroutine darkMatterProfileShapeSetConstructorDestructor

  subroutine darkMatterProfileShapeSetNodeTreeInitialize(self,node)
    !!{
    Initialize dark matter profile shape parameters.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileShapeSet), intent(inout), target  :: self
    type (treeNode                             ), intent(inout), target  :: node
    class(nodeComponentDarkMatterProfile       )               , pointer :: darkMatterProfile
  
    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    call darkMatterProfile%shapeSet(self%darkMatterProfileShape_%shape(node))
    return
  end subroutine darkMatterProfileShapeSetNodeTreeInitialize
    
  subroutine darkMatterProfileShapeSetNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the shape
    parameter of {\normalfont \ttfamily node} to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileShapeSet), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile       ), pointer       :: darkMatterProfile, darkMatterProfileParent
    
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    call darkMatterProfile%shapeSet(darkMatterProfileParent%shape())
    return
  end subroutine darkMatterProfileShapeSetNodePromote
