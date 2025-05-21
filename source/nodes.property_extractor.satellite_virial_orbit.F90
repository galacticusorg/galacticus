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
  <nodePropertyExtractor name="nodePropertyExtractorSatelliteVirialOrbit">
   <description>A property extractor class for {\normalfont \ttfamily keplerOrbit} objects.</description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorKeplerOrbit) :: nodePropertyExtractorSatelliteVirialOrbit
     !!{
     A property extractor for satellite node virial orbits.
     !!}
     private
   contains
     procedure :: extract => satelliteVirialOrbitExtract
  end type nodePropertyExtractorSatelliteVirialOrbit

  interface nodePropertyExtractorSatelliteVirialOrbit
     !!{
     Constructors for the \refClass{nodePropertyExtractorSatelliteVirialOrbit} extractor class.
     !!}
     module procedure satelliteVirialOrbitConstructorParameters
     module procedure satelliteVirialOrbitConstructorInternal
  end interface nodePropertyExtractorSatelliteVirialOrbit

contains

  function satelliteVirialOrbitConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorSatelliteVirialOrbit} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorSatelliteVirialOrbit)                              :: self
    type(inputParameters                          ), intent(inout)               :: parameters
    type(varying_string                           ), allocatable  , dimension(:) :: properties

    allocate(properties(parameters%count('properties')))
    !![
    <inputParameter>
      <name>properties</name>
      <source>parameters</source>
      <description>The set of properties of the orbit to output.</description>
    </inputParameter>
    !!]
    self=nodePropertyExtractorSatelliteVirialOrbit(properties)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteVirialOrbitConstructorParameters

  function satelliteVirialOrbitConstructorInternal(properties) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorSatelliteVirialOrbit} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSatelliteVirialOrbit)                              :: self
    type(varying_string                           ), intent(in   ), dimension(:) :: properties

    call self%initialize(properties,'satelliteVirialOrbit')
    return
  end function satelliteVirialOrbitConstructorInternal

  function satelliteVirialOrbitExtract(self,node,time,instance)
    !!{
    Implement a descendantNode output analysis.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    double precision                                           , dimension(:) , allocatable :: satelliteVirialOrbitExtract
    class           (nodePropertyExtractorSatelliteVirialOrbit), intent(inout), target      :: self
    type            (treeNode                                 ), intent(inout), target      :: node
    double precision                                           , intent(in   )              :: time
    type            (multiCounter                             ), intent(inout), optional    :: instance
    class           (nodeComponentSatellite                   )               , pointer     :: satellite
    type            (keplerOrbit                              )                             :: orbit
    !$GLC attributes unused :: time, instance

    satellite                   => node     %satellite       (     )
    orbit                       =  satellite%virialOrbit     (     )
    satelliteVirialOrbitExtract =  self     %extractFromOrbit(orbit)
    return
  end function satelliteVirialOrbitExtract
