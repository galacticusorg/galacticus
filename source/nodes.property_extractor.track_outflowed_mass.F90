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
  Implements a property extractor class for the mass and metal mass of gas outflowed to the \gls{cgm}.
  !!}
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorTrackOutflowedMass">
   <description>
    A property extractor class for the mass and metal mass of gas outflowed to the \gls{cgm}.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorTrackOutflowedMass
     !!{
     A property extractor class for the velocity dispersion at a set of radii.
     !!}
     private
     integer :: massOutflowedID, massMetalsOutflowedID
   contains
     procedure :: elementCount => trackOutflowedMassElementCount
     procedure :: extract      => trackOutflowedMassExtract
     procedure :: names        => trackOutflowedMassNames
     procedure :: descriptions => trackOutflowedMassDescriptions
     procedure :: unitsInSI    => trackOutflowedMassUnitsInSI
  end type nodePropertyExtractorTrackOutflowedMass

  interface nodePropertyExtractorTrackOutflowedMass
     !!{
     Constructors for the \refClass{nodePropertyExtractorTrackOutflowedMass} output analysis class.
     !!}
     module procedure trackOutflowedMassConstructorParameters
     module procedure trackOutflowedMassConstructorInternal
  end interface nodePropertyExtractorTrackOutflowedMass

contains

  function trackOutflowedMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorTrackOutflowedMass} property extractor class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodePropertyExtractorTrackOutflowedMass)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self=nodePropertyExtractorTrackOutflowedMass()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function trackOutflowedMassConstructorParameters

  function trackOutflowedMassConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorTrackOutflowedMass} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorTrackOutflowedMass) :: self

    !![
    <addMetaProperty component="hotHalo" name="massOutflowed"       id="self%massOutflowedID"       isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="hotHalo" name="massMetalsOutflowed" id="self%massMetalsOutflowedID" isEvolvable="yes" isCreator="no"/> 
    !!]
    return
  end function trackOutflowedMassConstructorInternal

  integer function trackOutflowedMassElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily trackOutflowedMass} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorTrackOutflowedMass), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    !$GLC attributes unused :: time

    trackOutflowedMassElementCount=2
    return
  end function trackOutflowedMassElementCount

  function trackOutflowedMassExtract(self,node,time,instance)
    !!{
    Implement a {\normalfont \ttfamily trackOutflowedMass} property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    double precision                                         , dimension(:) , allocatable :: trackOutflowedMassExtract
    class           (nodePropertyExtractorTrackOutflowedMass), intent(inout), target      :: self
    type            (treeNode                               ), intent(inout), target      :: node
    double precision                                         , intent(in   )              :: time
    type            (multiCounter                           ), intent(inout), optional    :: instance
    class           (nodeComponentHotHalo                   )               , pointer     :: hotHalo
    !$GLC attributes unused :: time, instance

    allocate(trackOutflowedMassExtract(2))
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Hot halo does not exist.
       trackOutflowedMassExtract=[                                                               & 
            &                     0.0d0                                                        , &
            &                     0.0d0                                                          &
            &                    ]
    class default
       trackOutflowedMassExtract=[                                                               &
            &                     hotHalo%floatRank0MetaPropertyGet(self%massOutflowedID      ), &
            &                     hotHalo%floatRank0MetaPropertyGet(self%massMetalsOutflowedID)  &            
            &                    ]
    end select
    return
  end function trackOutflowedMassExtract

  subroutine trackOutflowedMassNames(self,time,names)
    !!{
    Return the names of the {\normalfont \ttfamily trackOutflowedMass} properties.
    !!}
    implicit none
    class           (nodePropertyExtractorTrackOutflowedMass), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(2))
    names(1)=var_str('hotHaloMassOutflowTracked'      )
    names(2)=var_str('hotHaloMassMetalsOutflowTracked')
    return
  end subroutine trackOutflowedMassNames

  subroutine trackOutflowedMassDescriptions(self,time,descriptions)
    !!{
    Return descriptions of the {\normalfont \ttfamily trackOutflowedMass} property.
    !!}
    implicit none
    class           (nodePropertyExtractorTrackOutflowedMass), intent(inout)                             :: self
    double precision                                         , intent(in   )                             :: time
    type            (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(2))
    descriptions(1)=var_str('The mass of gas in the CGM which arrived via outflow from the central galaxy [M☉].'   )
    descriptions(2)=var_str('The mass of metals in the CGM which arrived via outflow from the central galaxy [M☉].')
    return
  end subroutine trackOutflowedMassDescriptions

  function trackOutflowedMassUnitsInSI(self,time)
    !!{
    Return the units of the {\normalfont \ttfamily trackOutflowedMass} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar
    implicit none
    double precision                                         , allocatable  , dimension(:) :: trackOutflowedMassUnitsInSI
    class           (nodePropertyExtractorTrackOutflowedMass), intent(inout)               :: self
    double precision                                         , intent(in   )               :: time
    !$GLC attributes unused :: time

    allocate(trackOutflowedMassUnitsInSI(2))
    trackOutflowedMassUnitsInSI=massSolar
    return
  end function trackOutflowedMassUnitsInSI
  
