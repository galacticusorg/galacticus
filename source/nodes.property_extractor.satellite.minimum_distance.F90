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
  <nodePropertyExtractor name="nodePropertyExtractorSatelliteMinimumDistance">
   <description>
     A node property extractor which extracts the minimum distance from the center that a satellite has ever reached in its
     current host halo.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorSatelliteMinimumDistance
     !!{     
     A node property extractor which extracts the minimum distance from the center that a satellite has ever reached in its
     current host halo.
     !!}
     private
     integer :: satelliteDistanceMinimumID
   contains
     procedure :: extract     => satelliteMinimumDistanceExtract
     procedure :: name        => satelliteMinimumDistanceName
     procedure :: description => satelliteMinimumDistanceDescription
     procedure :: unitsInSI   => satelliteMinimumDistanceUnitsInSI
  end type nodePropertyExtractorSatelliteMinimumDistance

  interface nodePropertyExtractorSatelliteMinimumDistance
     !!{
     Constructors for the \refClass{nodePropertyExtractorSatelliteMinimumDistance} output extractor class.
     !!}
     module procedure satelliteMinimumDistanceConstructorParameters
     module procedure satelliteMinimumDistanceConstructorInternal
  end interface nodePropertyExtractorSatelliteMinimumDistance

contains

  function satelliteMinimumDistanceConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorSatelliteMinimumDistance} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorSatelliteMinimumDistance)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=nodePropertyExtractorSatelliteMinimumDistance()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteMinimumDistanceConstructorParameters

  function satelliteMinimumDistanceConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorSatelliteMinimumDistance} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSatelliteMinimumDistance) :: self
    
    !![
    <addMetaProperty component="satellite" name="satelliteDistanceMinimum" id="self%satelliteDistanceMinimumID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function satelliteMinimumDistanceConstructorInternal

  double precision function satelliteMinimumDistanceExtract(self,node,instance)
    !!{
    Implement a satelliteMinimumDistance output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(nodePropertyExtractorSatelliteMinimumDistance), intent(inout), target   :: self
    type (treeNode                                     ), intent(inout), target   :: node
    type (multiCounter                                 ), intent(inout), optional :: instance
    class(nodeComponentSatellite                       )               , pointer  :: satellite
    !$GLC attributes unused :: instance

    if (node%isSatellite()) then
       satellite                       => node     %satellite                (                               )
       satelliteMinimumDistanceExtract =  satellite%floatRank0MetaPropertyGet(self%satelliteDistanceMinimumID)
    else
       satelliteMinimumDistanceExtract =  -1.0d0
    end if
    return
  end function satelliteMinimumDistanceExtract
  
  function satelliteMinimumDistanceName(self)
    !!{
    Return the names of the {\normalfont \ttfamily satelliteMinimumDistance} properties.
    !!}
    implicit none
    type (varying_string                               )                :: satelliteMinimumDistanceName
    class(nodePropertyExtractorSatelliteMinimumDistance), intent(inout) :: self
    !$GLC attributes unused :: self

    satelliteMinimumDistanceName=var_str('satelliteDistanceMinimum')
    return
  end function satelliteMinimumDistanceName

  function satelliteMinimumDistanceDescription(self)
    !!{
    Return the descriptions of the {\normalfont \ttfamily satelliteMinimumDistance} properties.
    !!}
    implicit none
    type (varying_string                               )                :: satelliteMinimumDistanceDescription
    class(nodePropertyExtractorSatelliteMinimumDistance), intent(inout) :: self
    !$GLC attributes unused :: self

    satelliteMinimumDistanceDescription=var_str('Minimum distance from the center of the current host ever reached.')
    return
  end function satelliteMinimumDistanceDescription

  double precision function satelliteMinimumDistanceUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily satelliteMinimumDistance} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorSatelliteMinimumDistance), intent(inout) :: self
    !$GLC attributes unused :: self

    satelliteMinimumDistanceUnitsInSI=megaParsec
    return
  end function satelliteMinimumDistanceUnitsInSI
