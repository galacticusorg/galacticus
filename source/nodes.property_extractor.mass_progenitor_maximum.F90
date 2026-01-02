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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassProgenitorMaximum">
   <description>
     A node property extractor which extracts the mass of the most massive progenitor of a node. Requires the
     \refClass{nodeOperatorMassProgenitorMaximum} node operator to be used to track the maximum progenitor mass.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassProgenitorMaximum
     !!{
     A property extractor which extracts the mass of the most massive progenitor of a node.
     !!}
     private
     integer :: massProgenitorMaximumID
   contains
     procedure :: extract     => massProgenitorMaximumExtract
     procedure :: name        => massProgenitorMaximumName
     procedure :: description => massProgenitorMaximumDescription
     procedure :: unitsInSI   => massProgenitorMaximumUnitsInSI
  end type nodePropertyExtractorMassProgenitorMaximum

  interface nodePropertyExtractorMassProgenitorMaximum
     !!{
     Constructors for the \refClass{nodePropertyExtractorMassProgenitorMaximum} output extractor class.
     !!}
     module procedure massProgenitorMaximumConstructorParameters
     module procedure massProgenitorMaximumConstructorInternal
  end interface nodePropertyExtractorMassProgenitorMaximum

contains

  function massProgenitorMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMassProgenitorMaximum} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorMassProgenitorMaximum)                :: self
    type(inputParameters                           ), intent(inout) :: parameters

    self=nodePropertyExtractorMassProgenitorMaximum()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massProgenitorMaximumConstructorParameters

  function massProgenitorMaximumConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMassProgenitorMaximum} output extractor property extractor class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type(nodePropertyExtractorMassProgenitorMaximum) :: self

    !![
    <addMetaProperty component="basic" name="massProgenitorMaximum" id="self%massProgenitorMaximumID" isEvolvable="no"/>
    !!]
    return
  end function massProgenitorMaximumConstructorInternal

  double precision function massProgenitorMaximumExtract(self,node,instance)
    !!{
    Implement a massProgenitorMaximum output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodePropertyExtractorMassProgenitorMaximum), intent(inout), target   :: self
    type (treeNode                                  ), intent(inout), target   :: node
    type (multiCounter                              ), intent(inout), optional :: instance
    class(nodeComponentBasic                        )               , pointer  :: basic
    !$GLC attributes unused :: instance

    basic                        => node %basic                    (                            )
    massProgenitorMaximumExtract =  basic%floatRank0MetaPropertyGet(self%massProgenitorMaximumID)
    return
  end function massProgenitorMaximumExtract

  function massProgenitorMaximumName(self)
    !!{
    Return the names of the {\normalfont \ttfamily massProgenitorMaximum} properties.
    !!}
    implicit none
    type (varying_string                            )                :: massProgenitorMaximumName
    class(nodePropertyExtractorMassProgenitorMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massProgenitorMaximumName=var_str('massProgenitorMaximum')
    return
  end function massProgenitorMaximumName

  function massProgenitorMaximumDescription(self)
    !!{
    Return the descriptions of the {\normalfont \ttfamily massProgenitorMaximum} properties.
    !!}
    implicit none
    type (varying_string                            )                :: massProgenitorMaximumDescription
    class(nodePropertyExtractorMassProgenitorMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massProgenitorMaximumDescription=var_str('The maximum mass of halo in which this node has ever been hosted. (Or $0$ if this node has never been a subhalo.)')
    return
  end function massProgenitorMaximumDescription

  double precision function massProgenitorMaximumUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily massProgenitorMaximum} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassProgenitorMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massProgenitorMaximumUnitsInSI=massSolar
    return
  end function massProgenitorMaximumUnitsInSI
