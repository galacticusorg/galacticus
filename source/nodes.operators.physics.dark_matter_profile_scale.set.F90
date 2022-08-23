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
  Implements a node operator class that sets dark matter profile scale radius.
  !!}

  use :: Dark_Matter_Profile_Scales, only : darkMatterProfileScaleRadiusClass
  
  !![
  <nodeOperator name="nodeOperatorDarkMatterProfileScaleSet">
   <description>
    A node operator class that sets dark matter profile scale radius.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkMatterProfileScaleSet
     !!{
     A node operator class that sets the dark matter profile scale radius in halos.
     !!}
     private
     class(darkMatterProfileScaleRadiusClass), pointer :: darkMatterProfileScaleRadius_ => null()
   contains
     final     ::                       darkMatterProfileScaleSetConstructorDestructor
     procedure :: nodeTreeInitialize => darkMatterProfileScaleSetNodeTreeInitialize
     procedure :: nodePromote        => darkMatterProfileScaleSetNodePromote
  end type nodeOperatorDarkMatterProfileScaleSet
  
  interface nodeOperatorDarkMatterProfileScaleSet
     !!{
     Constructors for the {\normalfont \ttfamily darkMatterProfileScaleSet} node operator class.
     !!}
     module procedure darkMatterProfileScaleSetConstructorParameters
     module procedure darkMatterProfileScaleSetConstructorInternal
  end interface nodeOperatorDarkMatterProfileScaleSet
  
contains
  
  function darkMatterProfileScaleSetConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily darkMatterProfileScaleSet} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfileScaleSet)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(darkMatterProfileScaleRadiusClass    ), pointer       :: darkMatterProfileScaleRadius_

    !![
    <objectBuilder class="darkMatterProfileScaleRadius" name="darkMatterProfileScaleRadius_" source="parameters"/>
    !!]
    self=nodeOperatorDarkMatterProfileScaleSet(darkMatterProfileScaleRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileScaleRadius_"/>
    !!]
    return
  end function darkMatterProfileScaleSetConstructorParameters

  function darkMatterProfileScaleSetConstructorInternal(darkMatterProfileScaleRadius_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily darkMatterProfileScaleSet} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorDarkMatterProfileScaleSet)                        :: self
    class(darkMatterProfileScaleRadiusClass    ), intent(in   ), target :: darkMatterProfileScaleRadius_
    !![
    <constructorAssign variables="*darkMatterProfileScaleRadius_"/>
    !!]

    return
  end function darkMatterProfileScaleSetConstructorInternal

  subroutine darkMatterProfileScaleSetConstructorDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily darkMatterProfileScaleSet} dark matter halo profile scale radius class.
    !!}
    implicit none
    type(nodeOperatorDarkMatterProfileScaleSet), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileScaleRadius_"/>
    !!]
    return
  end subroutine darkMatterProfileScaleSetConstructorDestructor

  subroutine darkMatterProfileScaleSetNodeTreeInitialize(self,node)
    !!{
    Initialize dark matter profile scale radii.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileScaleSet), intent(inout), target  :: self
    type (treeNode                             ), intent(inout), target  :: node
    class(nodeComponentDarkMatterProfile       )               , pointer :: darkMatterProfile

    darkMatterProfile => node%darkMatterProfile(autoCreate=.true.)
    call darkMatterProfile%scaleSet(self%darkMatterProfileScaleRadius_%radius(node))
    return
  end subroutine darkMatterProfileScaleSetNodeTreeInitialize
    
  subroutine darkMatterProfileScaleSetNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the scale radius
    of {\normalfont \ttfamily node} to be that of its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorDarkMatterProfileScaleSet), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile       ), pointer       :: darkMatterProfile, darkMatterProfileParent
    
    darkMatterProfile       => node       %darkMatterProfile()
    darkMatterProfileParent => node%parent%darkMatterProfile()
    call darkMatterProfile%scaleSet(darkMatterProfileParent%scale())
    return
  end subroutine darkMatterProfileScaleSetNodePromote
