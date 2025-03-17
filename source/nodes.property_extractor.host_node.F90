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
Implements an output analysis property extractor class that extracts a property from the host node of the given node.
!!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorHostNode">
   <description>An output analysis property extractor class that extracts a property from the host node of the given node.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorHostNode
     !!{
     A property extractor output analysis class that extracts a property from the host node of the given node.
     !!}
     private
     class(nodePropertyExtractorScalar), pointer :: nodePropertyExtractor_ => null()
   contains
     final     ::                hostNodeDestructor
     procedure :: extract     => hostNodeExtract
     procedure :: name        => hostNodeName
     procedure :: description => hostNodeDescription
     procedure :: unitsInSI   => hostNodeUnitsInSI
  end type nodePropertyExtractorHostNode

  interface nodePropertyExtractorHostNode
     !!{
     Constructors for the ``hostNode'' node property extractor class.
     !!}
     module procedure hostNodeConstructorParameters
     module procedure hostNodeConstructorInternal
  end interface nodePropertyExtractorHostNode

contains

  function hostNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``hostNode'' node property extractor class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorHostNode)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(nodePropertyExtractorClass   ), pointer       :: nodePropertyExtractor_
    
    !![
    <objectBuilder class="nodePropertyExtractor" name="nodePropertyExtractor_" source="parameters"/>
    !!]
    select type (nodePropertyExtractor_)
    class is (nodePropertyExtractorScalar)
       self=nodePropertyExtractorHostNode(nodePropertyExtractor_)
    class default
       call Error_Report('extracted property must be a real scalar'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodePropertyExtractor_"/>
    !!]
    return
  end function hostNodeConstructorParameters

  function hostNodeConstructorInternal(nodePropertyExtractor_) result(self)
    !!{
    Internal constructor for the ``hostNode'' node property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorHostNode)                        :: self
    class(nodePropertyExtractorScalar  ), intent(in   ), target :: nodePropertyExtractor_
    !![
    <constructorAssign variables="*nodePropertyExtractor_"/>
    !!]

    return
  end function hostNodeConstructorInternal
  
  subroutine hostNodeDestructor(self)
    !!{
    Destructor for  the ``hostNode'' node property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorHostNode), intent(inout) :: self

    !![
    <objectDestructor name="self%nodePropertyExtractor_"/>
    !!]
    return
  end subroutine hostNodeDestructor
  
  double precision function hostNodeExtract(self,node,instance)
    !!{
    Implement a hostNode output analysis.
    !!}
    implicit none
    class(nodePropertyExtractorHostNode), intent(inout), target   :: self
    type (treeNode                     ), intent(inout), target   :: node
    type (multiCounter                 ), intent(inout), optional :: instance
    type (treeNode                     ), pointer                 :: nodeHost

    nodeHost => node
    do while (nodeHost%isSatellite())
       nodeHost => nodeHost%parent
    end do
    hostNodeExtract=self%nodePropertyExtractor_%extract(nodeHost,instance)
    return
  end function hostNodeExtract


  function hostNodeName(self)
    !!{
    Return the name of the hostNode property.
    !!}
    use :: String_Handling, only : String_Upper_Case_First
    implicit none
    type (varying_string               )                :: hostNodeName
    class(nodePropertyExtractorHostNode), intent(inout) :: self

    hostNodeName=var_str('host')//String_Upper_Case_First(char(self%nodePropertyExtractor_%name()))
    return
  end function hostNodeName

  function hostNodeDescription(self)
    !!{
    Return a description of the hostNode property.
    !!}
    implicit none
    type (varying_string               )                :: hostNodeDescription
    class(nodePropertyExtractorHostNode), intent(inout) :: self

    hostNodeDescription=self%nodePropertyExtractor_%description()//' (of the host node)'
    return
  end function hostNodeDescription

  double precision function hostNodeUnitsInSI(self)
    !!{
    Return the units of the hostNode property in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorHostNode), intent(inout) :: self

    hostNodeUnitsInSI=self%nodePropertyExtractor_%unitsInSI()
    return
  end function hostNodeUnitsInSI
