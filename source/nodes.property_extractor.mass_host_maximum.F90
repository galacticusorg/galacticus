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
  <nodePropertyExtractor name="nodePropertyExtractorMassHostMaximum">
   <description>
     A node property extractor which extracts the mass of the most massive node in which a node has been hosted. Requires the
     \refClass{nodeOperatorMassHostMaximum} node operator to be used to track the maximum host mass.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassHostMaximum
     !!{
     A property extractor which extracts the mass of the most massive node in which a node has been hosted.
     !!}
     private
     integer :: massHostMaximumID
   contains
     procedure :: extract     => massHostMaximumExtract
     procedure :: name        => massHostMaximumName
     procedure :: description => massHostMaximumDescription
     procedure :: unitsInSI   => massHostMaximumUnitsInSI
  end type nodePropertyExtractorMassHostMaximum

  interface nodePropertyExtractorMassHostMaximum
     !!{
     Constructors for the \refClass{nodePropertyExtractorMassHostMaximum} output extractor class.
     !!}
     module procedure massHostMaximumConstructorParameters
     module procedure massHostMaximumConstructorInternal
  end interface nodePropertyExtractorMassHostMaximum

contains

  function massHostMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMassHostMaximum} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorMassHostMaximum)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=nodePropertyExtractorMassHostMaximum()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massHostMaximumConstructorParameters

  function massHostMaximumConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMassHostMaximum} output extractor property extractor class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type(nodePropertyExtractorMassHostMaximum) :: self

    !![
    <addMetaProperty component="basic" name="massHostMaximum" id="self%massHostMaximumID" isEvolvable="no"/>
    !!]
    return
  end function massHostMaximumConstructorInternal

  double precision function massHostMaximumExtract(self,node,instance)
    !!{
    Implement a massHostMaximum output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodePropertyExtractorMassHostMaximum), intent(inout), target   :: self
    type (treeNode                            ), intent(inout), target   :: node
    type (multiCounter                        ), intent(inout), optional :: instance
    class(nodeComponentBasic                  )               , pointer  :: basic
    !$GLC attributes unused :: instance

    basic                  => node %basic                    (                      )
    massHostMaximumExtract =  basic%floatRank0MetaPropertyGet(self%massHostMaximumID)
    return
  end function massHostMaximumExtract

  function massHostMaximumName(self)
    !!{
    Return the names of the {\normalfont \ttfamily massHostMaximum} properties.
    !!}
    implicit none
    type (varying_string                      )                :: massHostMaximumName
    class(nodePropertyExtractorMassHostMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massHostMaximumName=var_str('massHostMaximum')
    return
  end function massHostMaximumName

  function massHostMaximumDescription(self)
    !!{
    Return the descriptions of the {\normalfont \ttfamily massHostMaximum} properties.
    !!}
    implicit none
    type (varying_string                      )                :: massHostMaximumDescription
    class(nodePropertyExtractorMassHostMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massHostMaximumDescription=var_str('The maximum mass of halo in which this node has ever been hosted. (Or $0$ if this node has never been a subhalo.)')
    return
  end function massHostMaximumDescription

  double precision function massHostMaximumUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily massHostMaximum} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassHostMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massHostMaximumUnitsInSI=massSolar
    return
  end function massHostMaximumUnitsInSI
