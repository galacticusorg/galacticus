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
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMergersPhysical">
   <description>
     A node property extractor which extracts the physical properties of galaxy-galaxy mergers.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorGalaxyMergersPhysical
     !!{
     A property extractor which extracts the physical properties of galaxy-galaxy mergers.
     !!}
     private
     integer :: galaxyMergerHostMassGasID     , galaxyMergerHostMassStellarID     , &
          &     galaxyMergerSatelliteMassGasID, galaxyMergerSatelliteMassStellarID, &
          &     galaxyMergerTimeID 
   contains
     procedure :: elementCount => galaxyMergersPhysicalElementCount
     procedure :: extract      => galaxyMergersPhysicalExtract
     procedure :: names        => galaxyMergersPhysicalNames
     procedure :: descriptions => galaxyMergersPhysicalDescriptions
     procedure :: unitsInSI    => galaxyMergersPhysicalUnitsInSI
  end type nodePropertyExtractorGalaxyMergersPhysical

  interface nodePropertyExtractorGalaxyMergersPhysical
     !!{
     Constructors for the \refClass{nodePropertyExtractorGalaxyMergersPhysical} output extractor class.
     !!}
     module procedure galaxyMergersPhysicalConstructorParameters
     module procedure galaxyMergersPhysicalConstructorInternal
  end interface nodePropertyExtractorGalaxyMergersPhysical

contains

  function galaxyMergersPhysicalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorGalaxyMergersPhysical} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMergersPhysical)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMergersPhysical()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMergersPhysicalConstructorParameters

  function galaxyMergersPhysicalConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorGalaxyMergersPhysical} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGalaxyMergersPhysical) :: self
    
    !![
    <addMetaProperty component="basic" name="galaxyMergersTime"                id="self%galaxyMergerTimeID"                 rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="galaxyMergerHostMassGas"          id="self%galaxyMergerHostMassGasID"          rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="galaxyMergerSatelliteMassGas"     id="self%galaxyMergerSatelliteMassGasID"     rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="galaxyMergerHostMassStellar"      id="self%galaxyMergerHostMassStellarID"      rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="galaxyMergerSatelliteMassStellar" id="self%galaxyMergerSatelliteMassStellarID" rank="1" isCreator="no"/>
    !!]
    return
  end function galaxyMergersPhysicalConstructorInternal

  integer function galaxyMergersPhysicalElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergersPhysical), intent(inout) :: self

    galaxyMergersPhysicalElementCount=5
    return
  end function galaxyMergersPhysicalElementCount

  function galaxyMergersPhysicalExtract(self,node,instance) result(galaxyMergers)
    !!{
    Implement a galaxyMergersPhysical output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                                            , dimension(:,:), allocatable :: galaxyMergers
    class           (nodePropertyExtractorGalaxyMergersPhysical), intent(inout)               :: self
    type            (treeNode                                  ), intent(inout)               :: node
    type            (multiCounter                              ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic                        )                , pointer     :: basic
    double precision                                            , dimension(:  ), allocatable :: massesHostGas     , massesHostStellar     , &
         &                                                                                       massesSatelliteGas, massesSatelliteStellar, &
         &                                                                                       times
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: times, massesHostGas, massesSatelliteGas, massesHostStellar, massesSatelliteStellar
    
    basic                  => node %basic                    (                                       )
    times                  =  basic%floatRank1MetaPropertyGet(self%galaxyMergerTimeID                )
    massesHostGas          =  basic%floatRank1MetaPropertyGet(self%galaxyMergerHostMassGasID         )
    massesSatelliteGas     =  basic%floatRank1MetaPropertyGet(self%galaxyMergerSatelliteMassGasID    )
    massesHostStellar      =  basic%floatRank1MetaPropertyGet(self%galaxyMergerHostMassStellarID     )
    massesSatelliteStellar =  basic%floatRank1MetaPropertyGet(self%galaxyMergerSatelliteMassStellarID)
    allocate(galaxyMergers(size(times),5))
    galaxyMergers(:,1)=times
    galaxyMergers(:,2)=massesHostGas
    galaxyMergers(:,3)=massesSatelliteGas
    galaxyMergers(:,4)=massesHostStellar
    galaxyMergers(:,5)=massesSatelliteStellar
    return
  end function galaxyMergersPhysicalExtract
  
  subroutine galaxyMergersPhysicalNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily galaxyMergersPhysical} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergersPhysical), intent(inout)                             :: self
    type (varying_string                            ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(5))
    names(1)=var_str('galaxyMergersTime'                )
    names(2)=var_str('galaxyMergersMassHostGas'         )
    names(3)=var_str('galaxyMergersMassSatelliteGas'    )
    names(4)=var_str('galaxyMergersMassHostStellar'     )
    names(5)=var_str('galaxyMergersMassSatelliteStellar')
    return
  end subroutine galaxyMergersPhysicalNames

  subroutine galaxyMergersPhysicalDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily galaxyMergersPhysical} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyMergersPhysical), intent(inout)                             :: self
    type (varying_string                            ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(5))
    descriptions(1)=var_str('Times of galaxy-galaxy mergers.'                        )
    descriptions(2)=var_str('Host galaxy gas mass in galaxy-galaxy mergers.'         )
    descriptions(3)=var_str('Satellite galaxy gas mass in galaxy-galaxy mergers.'    )
    descriptions(4)=var_str('Host galaxy stellar mass in galaxy-galaxy mergers.'     )
    descriptions(5)=var_str('Satellite galaxy stellar mass in galaxy-galaxy mergers.')
    return
  end subroutine galaxyMergersPhysicalDescriptions

  function galaxyMergersPhysicalUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily galaxyMergersPhysical} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, massSolar
    implicit none
    double precision                                            , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorGalaxyMergersPhysical), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(5))
    unitsInSI(1)=gigaYear
    unitsInSI(2)=massSolar
    unitsInSI(3)=massSolar
    unitsInSI(4)=massSolar
    unitsInSI(5)=massSolar
    return
  end function galaxyMergersPhysicalUnitsInSI
