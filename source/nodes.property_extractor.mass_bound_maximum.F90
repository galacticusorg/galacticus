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

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassBoundMaximum">
   <description>
     A node property extractor which extracts the maximum bound mass of the node. Requires the
     \refClass{nodeOperatorMassBoundMaximum} node operator to be used to track the maximum bound mass.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorMassBoundMaximum
     !!{
     A property extractor which extracts the maximum bound mass of the node.
     !!}
     private
     integer :: massBoundMaximumID
   contains
     procedure :: extract     => massBoundMaximumExtract
     procedure :: name        => massBoundMaximumName
     procedure :: description => massBoundMaximumDescription
     procedure :: unitsInSI   => massBoundMaximumUnitsInSI
  end type nodePropertyExtractorMassBoundMaximum

  interface nodePropertyExtractorMassBoundMaximum
     !!{
     Constructors for the {\normalfont \ttfamily massBoundMaximum} output extractor class.
     !!}
     module procedure massBoundMaximumConstructorParameters
     module procedure massBoundMaximumConstructorInternal
  end interface nodePropertyExtractorMassBoundMaximum

contains

  function massBoundMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily massBoundMaximum} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorMassBoundMaximum)                :: self
    type(inputParameters                      ), intent(inout) :: parameters

    self=nodePropertyExtractorMassBoundMaximum()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massBoundMaximumConstructorParameters

  function massBoundMaximumConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily massBoundMaximum} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorMassBoundMaximum) :: self

    !![
    <addMetaProperty component="satellite" name="massBoundMaximum" id="self%massBoundMaximumID" isEvolvable="no"/>
    !!]
    return
  end function massBoundMaximumConstructorInternal

  double precision function massBoundMaximumExtract(self,node,instance)
    !!{
    Implement a massBoundMaximum output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite
    implicit none
    class(nodePropertyExtractorMassBoundMaximum), intent(inout), target   :: self
    type (treeNode                             ), intent(inout), target   :: node
    type (multiCounter                         ), intent(inout), optional :: instance
    class(nodeComponentBasic                   )               , pointer  :: basic
    class(nodeComponentSatellite               )               , pointer  :: satellite
    !$GLC attributes unused :: instance

    basic     => node%basic    ()
    satellite => node%satellite()
    select type (satellite)
    type is (nodeComponentSatellite)
       ! No satellite exists, so use the basic mass of the node.
       massBoundMaximumExtract=basic    %mass                     (                       )
    class default
       massBoundMaximumExtract=satellite%floatRank0MetaPropertyGet(self%massBoundMaximumID)
    end select
    return
  end function massBoundMaximumExtract
   
  function massBoundMaximumName(self)
    !!{
    Return the names of the {\normalfont \ttfamily massBoundMaximum} properties.
    !!}
    implicit none
    type (varying_string                       )                :: massBoundMaximumName
    class(nodePropertyExtractorMassBoundMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massBoundMaximumName=var_str('massBoundMaximum')
    return
  end function massBoundMaximumName

  function massBoundMaximumDescription(self)
    !!{
    Return the descriptions of the {\normalfont \ttfamily massBoundMaximum} properties.
    !!}
    implicit none
    type (varying_string                       )                :: massBoundMaximumDescription
    class(nodePropertyExtractorMassBoundMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massBoundMaximumDescription=var_str('The maximum bound mass ever attained by this halo.')
    return
  end function massBoundMaximumDescription

  double precision function massBoundMaximumUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily massBoundMaximum} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    class(nodePropertyExtractorMassBoundMaximum), intent(inout) :: self
    !$GLC attributes unused :: self

    massBoundMaximumUnitsInSI=massSolar
    return
  end function massBoundMaximumUnitsInSI
