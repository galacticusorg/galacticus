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
Implements a last-defined virial radius property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusVirialLastDefined">
   <description>A last-defined virial radius property extractor class.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusVirialLastDefined
     !!{
     A last-defined virial radius infall property extractor class.
     !!}
     private
     integer :: radiusVirialLastDefinedID
   contains
     procedure :: extract     => radiusVirialLastDefinedExtract
     procedure :: name        => radiusVirialLastDefinedName
     procedure :: description => radiusVirialLastDefinedDescription
     procedure :: unitsInSI   => radiusVirialLastDefinedUnitsInSI
  end type nodePropertyExtractorRadiusVirialLastDefined

  interface nodePropertyExtractorRadiusVirialLastDefined
     !!{
     Constructors for the \refClass{nodePropertyExtractorRadiusVirialLastDefined} output analysis class.
     !!}
     module procedure radiusVirialLastDefinedConstructorParameters
     module procedure radiusVirialLastDefinedConstructorInternal
  end interface nodePropertyExtractorRadiusVirialLastDefined

contains

  function radiusVirialLastDefinedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorRadiusVirialLastDefined} output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorRadiusVirialLastDefined)                :: self
    type(inputParameters                             ), intent(inout) :: parameters
    
    self=nodePropertyExtractorRadiusVirialLastDefined()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function radiusVirialLastDefinedConstructorParameters

  function radiusVirialLastDefinedConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorRadiusVirialLastDefined} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorRadiusVirialLastDefined) :: self

    !![
    <addMetaProperty component="basic" name="radiusVirialLastDefined" id="self%radiusVirialLastDefinedID" isEvolvable="no" isCreator="no"/>
    !!]
    return
  end function radiusVirialLastDefinedConstructorInternal

  double precision function radiusVirialLastDefinedExtract(self,node,instance)
    !!{
    Implement a time of first infall property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(nodePropertyExtractorRadiusVirialLastDefined), intent(inout), target   :: self
    type (treeNode                                    ), intent(inout), target   :: node
    type (multiCounter                                ), intent(inout), optional :: instance
    class(nodeComponentBasic                          ), pointer                 :: basic
    !$GLC attributes unused :: self, instance

    basic                          => node %basic                    (                              )
    radiusVirialLastDefinedExtract =  basic%floatRank0MetaPropertyGet(self%radiusVirialLastDefinedID)
    return
  end function radiusVirialLastDefinedExtract

  function radiusVirialLastDefinedName(self)
    !!{
    Return the name of the time of first infall property.
    !!}
    implicit none
    type (varying_string                              )                :: radiusVirialLastDefinedName
    class(nodePropertyExtractorRadiusVirialLastDefined), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVirialLastDefinedName=var_str('darkMatterOnlyRadiusVirialLastDefined')
    return
  end function radiusVirialLastDefinedName

  function radiusVirialLastDefinedDescription(self)
    !!{
    Return a description of the time of first infall property.
    !!}
    implicit none
    type (varying_string                              )                :: radiusVirialLastDefinedDescription
    class(nodePropertyExtractorRadiusVirialLastDefined), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVirialLastDefinedDescription=var_str('Virial radius of the dark matter only node at the time it was last explictly defined [Mpc].')
    return
  end function radiusVirialLastDefinedDescription

  double precision function radiusVirialLastDefinedUnitsInSI(self)
    !!{
    Return the units of the time of first infall property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusVirialLastDefined), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusVirialLastDefinedUnitsInSI=megaParsec
    return
  end function radiusVirialLastDefinedUnitsInSI
