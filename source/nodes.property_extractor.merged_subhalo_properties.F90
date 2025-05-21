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

  use :: Kepler_Orbits, only : keplerOrbitCount
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorMergedSubhaloProperties">
   <description>
     A node property extractor which extracts properties of merged subhalo orbits.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorMergedSubhaloProperties
     !!{
     A property extractor which extracts properties of merged subhalo orbits.
     !!}
     private
     integer :: mergedSubhaloIDs(keplerOrbitCount)
   contains
     procedure :: elementCount => mergedSubhaloPropertiesElementCount
     procedure :: extract      => mergedSubhaloPropertiesExtract
     procedure :: names        => mergedSubhaloPropertiesNames
     procedure :: descriptions => mergedSubhaloPropertiesDescriptions
     procedure :: unitsInSI    => mergedSubhaloPropertiesUnitsInSI
  end type nodePropertyExtractorMergedSubhaloProperties

  interface nodePropertyExtractorMergedSubhaloProperties
     !!{
     Constructors for the \refClass{nodePropertyExtractorMergedSubhaloProperties} output extractor class.
     !!}
     module procedure mergedSubhaloPropertiesConstructorParameters
     module procedure mergedSubhaloPropertiesConstructorInternal
  end interface nodePropertyExtractorMergedSubhaloProperties

contains

  function mergedSubhaloPropertiesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorMergedSubhaloProperties} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorMergedSubhaloProperties)                :: self
    type(inputParameters                             ), intent(inout) :: parameters

    self=nodePropertyExtractorMergedSubhaloProperties()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function mergedSubhaloPropertiesConstructorParameters

  function mergedSubhaloPropertiesConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorMergedSubhaloProperties} output extractor property extractor class.
    !!}
    use :: Kepler_Orbits, only : keplerOrbitTimeInitial     , keplerOrbitMassSatellite, keplerOrbitMassHost, keplerOrbitRadius, &
         &                       keplerOrbitRadiusPericenter, keplerOrbitTimeCurrent
    implicit none
    type(nodePropertyExtractorMergedSubhaloProperties) :: self
    
    !![
    <addMetaProperty component="basic" name="mergedSubhaloTimeCurrent"      id="self%mergedSubhaloIDs(keplerOrbitTimeCurrent     %ID)" rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="mergedSubhaloTimeInitial"      id="self%mergedSubhaloIDs(keplerOrbitTimeInitial     %ID)" rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="mergedSubhaloMassSatellite"    id="self%mergedSubhaloIDs(keplerOrbitMassSatellite   %ID)" rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="mergedSubhaloMassHost"         id="self%mergedSubhaloIDs(keplerOrbitMassHost        %ID)" rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="mergedSubhaloRadius"           id="self%mergedSubhaloIDs(keplerOrbitRadius          %ID)" rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="mergedSubhaloRadiusPericenter" id="self%mergedSubhaloIDs(keplerOrbitRadiusPericenter%ID)" rank="1" isCreator="no"/>
    !!]
    return
  end function mergedSubhaloPropertiesConstructorInternal

  integer function mergedSubhaloPropertiesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorMergedSubhaloProperties), intent(inout) :: self

    mergedSubhaloPropertiesElementCount=6
    return
  end function mergedSubhaloPropertiesElementCount

  function mergedSubhaloPropertiesExtract(self,node,instance) result(propertiesOrbitalMergedSubhalos)
    !!{
    Implement a mergedSubhaloProperties output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    use :: Kepler_Orbits   , only : keplerOrbitTimeInitial, keplerOrbitMassSatellite, keplerOrbitMassHost, keplerOrbitRadiusPericenter, &
         &                          keplerOrbitRadius     , keplerOrbitTimeCurrent
    implicit none
    double precision                                              , dimension(:,:), allocatable :: propertiesOrbitalMergedSubhalos
    class           (nodePropertyExtractorMergedSubhaloProperties), intent(inout)               :: self
    type            (treeNode                                    ), intent(inout)               :: node
    type            (multiCounter                                ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic                          )                , pointer     :: basic
    double precision                                              , dimension(:  ), allocatable :: property
    integer                                                                                     :: i                              , ID
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: timesLastIsolated

    basic => node%basic()
    do i=1,6
       select case (i)
       case (1)
          ID=keplerOrbitTimeInitial     %ID
       case (2)
          ID=keplerOrbitTimeCurrent     %ID
       case (3)
          ID=keplerOrbitMassSatellite   %ID
       case (4)
          ID=keplerOrbitMassHost        %ID
       case (5)
          ID=keplerOrbitRadius          %ID
       case (6)
          ID=keplerOrbitRadiusPericenter%ID
       end select
       property=basic%floatRank1MetaPropertyGet(self%mergedSubhaloIDs(ID))
       if (.not.allocated(propertiesOrbitalMergedSubhalos)) allocate(propertiesOrbitalMergedSubhalos(size(property),6))
       propertiesOrbitalMergedSubhalos(:,i)=property
       deallocate(property)
    end do
    return
  end function mergedSubhaloPropertiesExtract
  
  subroutine mergedSubhaloPropertiesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily mergedSubhaloProperties} properties.
    !!}
    use :: Kepler_Orbits  , only : enumerationKeplerOrbitDecode, keplerOrbitTimeInitial, keplerOrbitMassSatellite, keplerOrbitMassHost, &
         &                         keplerOrbitRadiusPericenter ,  keplerOrbitRadius    , keplerOrbitTimeCurrent
    use :: String_Handling, only : String_Upper_Case_First
    implicit none
    class  (nodePropertyExtractorMergedSubhaloProperties), intent(inout)                             :: self
    type   (varying_string                              ), intent(inout), dimension(:) , allocatable :: names
    integer                                                                                          :: i    , ID
    !$GLC attributes unused :: self

    allocate(names(6))
    do i=1,6
        select case (i)
       case (1)
          ID=keplerOrbitTimeInitial     %ID
       case (2)
          ID=keplerOrbitTimeCurrent     %ID
       case (3)
          ID=keplerOrbitMassSatellite   %ID
       case (4)
          ID=keplerOrbitMassHost        %ID
       case (5)
          ID=keplerOrbitRadius          %ID
       case (6)
          ID=keplerOrbitRadiusPericenter%ID
       end select
       names(i)=var_str('mergedSubhalo')//String_Upper_Case_First(char(enumerationKeplerOrbitDecode(ID,includePrefix=.false.)))
    end do
    return
  end subroutine mergedSubhaloPropertiesNames

  subroutine mergedSubhaloPropertiesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily mergedSubhaloProperties} properties.
    !!}
    use :: Kepler_Orbits, only : enumerationKeplerOrbitDescription, keplerOrbitTimeInitial, keplerOrbitMassSatellite, keplerOrbitMassHost, &
         &                       keplerOrbitRadiusPericenter      , keplerOrbitRadius     , keplerOrbitTimeCurrent
    implicit none
    class  (nodePropertyExtractorMergedSubhaloProperties), intent(inout)                             :: self
    type   (varying_string                              ), intent(inout), dimension(:) , allocatable :: descriptions
    integer                                                                                          :: i           , ID
    !$GLC attributes unused :: self

    allocate(descriptions(6))
    do i=1,6
        select case (i)
       case (1)
          ID=keplerOrbitTimeInitial     %ID
       case (2)
          ID=keplerOrbitTimeCurrent     %ID
       case (3)
          ID=keplerOrbitMassSatellite   %ID
       case (4)
          ID=keplerOrbitMassHost        %ID
       case (5)
          ID=keplerOrbitRadius          %ID
       case (6)
          ID=keplerOrbitRadiusPericenter%ID
       end select
       descriptions(i)=var_str('Merged subhalos: ')//enumerationKeplerOrbitDescription(ID)
    end do
    return
  end subroutine mergedSubhaloPropertiesDescriptions

  function mergedSubhaloPropertiesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily mergedSubhaloProperties} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear, massSolar, megaParsec
    implicit none
    double precision                                              , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorMergedSubhaloProperties), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(6))
    unitsInSI=[gigaYear,massSolar,massSolar,megaParsec,megaParsec]
    return
  end function mergedSubhaloPropertiesUnitsInSI
