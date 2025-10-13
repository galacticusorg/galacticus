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
  Implements a node operator class that tracks the last-defined virial radius.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScale, darkMatterHaloScaleClass
  
  !![
  <nodeOperator name="nodeOperatorRadiusVirialLastDefined">
    <description>
      A node operator class that tracks the last-defined virial radius. Intended to be paired with the
      \refClass{nodePropertyExtractorRadiusVirialLastDefined} class to extract these times for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorRadiusVirialLastDefined
     !!{
     A node operator class that tracks the last-defined virial radius.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_      => null()
     integer                                    :: radiusVirialLastDefinedID
   contains
     final     ::                   radiusVirialLastDefinedDestructor
     procedure :: nodeInitialize => radiusVirialLastDefinedNodeInitialize
  end type nodeOperatorRadiusVirialLastDefined
  
  interface nodeOperatorRadiusVirialLastDefined
     !!{
     Constructors for the \refClass{nodeOperatorRadiusVirialLastDefined} node operator class.
     !!}
     module procedure radiusVirialLastDefinedConstructorParameters
     module procedure radiusVirialLastDefinedConstructorInternal
  end interface nodeOperatorRadiusVirialLastDefined
  
contains

  function radiusVirialLastDefinedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorRadiusVirialLastDefined} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorRadiusVirialLastDefined)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_

    !![
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=nodeOperatorRadiusVirialLastDefined(darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function radiusVirialLastDefinedConstructorParameters

  function radiusVirialLastDefinedConstructorInternal(darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorRadiusVirialLastDefined} node operator class.
    !!}
    implicit none
    type (nodeOperatorRadiusVirialLastDefined)                        :: self
    class(darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterHaloScale_"/>
    !!]

    !![
    <addMetaProperty component="basic" name="radiusVirialLastDefined" id="self%radiusVirialLastDefinedID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function radiusVirialLastDefinedConstructorInternal
  
  subroutine radiusVirialLastDefinedDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorRadiusVirialLastDefined} property extractor class.
    !!}
    implicit none
    type(nodeOperatorRadiusVirialLastDefined), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine radiusVirialLastDefinedDestructor

  subroutine radiusVirialLastDefinedNodeInitialize(self,node)
    !!{
    Initialize the last-defined virial radius for this node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorRadiusVirialLastDefined), intent(inout), target  :: self
    type (treeNode                           ), intent(inout), target  :: node
    type (treeNode                           )               , pointer :: nodeDescendant
    class(nodeComponentBasic                 )               , pointer :: basic

    ! Find the last-defined node.
    nodeDescendant => node
    do while (nodeDescendant%isPrimaryProgenitor())
       nodeDescendant => nodeDescendant%parent
    end do
    basic => node%basic()
    call basic%floatRank0MetaPropertySet(self%radiusVirialLastDefinedID,self%darkMatterHaloScale_%radiusVirial(nodeDescendant))
    return
  end subroutine radiusVirialLastDefinedNodeInitialize
