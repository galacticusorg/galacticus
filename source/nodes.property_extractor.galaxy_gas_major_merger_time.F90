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
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyGasMajorMergerTime">
   <description>
     A node property extractor which extracts the times of gas-mass-based major mergers for each galaxy.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorGalaxyGasMajorMergerTime
     !!{
     A property extractor which extracts the times of gas-mass-based major merger for each galaxy.
     !!}
     private
     integer :: galaxyGasMajorMergerTimeID
   contains
     procedure :: elementCount => galaxyGasMajorMergerTimeElementCount
     procedure :: extract      => galaxyGasMajorMergerTimeExtract
     procedure :: names        => galaxyGasMajorMergerTimeNames
     procedure :: descriptions => galaxyGasMajorMergerTimeDescriptions
     procedure :: unitsInSI    => galaxyGasMajorMergerTimeUnitsInSI
  end type nodePropertyExtractorGalaxyGasMajorMergerTime

  interface nodePropertyExtractorGalaxyGasMajorMergerTime
     !!{
     Constructors for the \refClass{nodePropertyExtractorGalaxyGasMajorMergerTime} output extractor class.
     !!}
     module procedure galaxyGasMajorMergerTimeConstructorParameters
     module procedure galaxyGasMajorMergerTimeConstructorInternal
  end interface nodePropertyExtractorGalaxyGasMajorMergerTime

contains

  function galaxyGasMajorMergerTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorGalaxyGasMajorMergerTime} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyGasMajorMergerTime)                :: self
    type(inputParameters                              ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyGasMajorMergerTime()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyGasMajorMergerTimeConstructorParameters

  function galaxyGasMajorMergerTimeConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorGalaxyGasMajorMergerTime} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGalaxyGasMajorMergerTime) :: self
    
    !![
    <addMetaProperty component="basic" name="galaxyGasMajorMergerTime" id="self%galaxyGasMajorMergerTimeID" rank="1" isCreator="no"/>
    !!]
    return
  end function galaxyGasMajorMergerTimeConstructorInternal

  integer function galaxyGasMajorMergerTimeElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyGasMajorMergerTime), intent(inout) :: self

    galaxyGasMajorMergerTimeElementCount=1
    return
  end function galaxyGasMajorMergerTimeElementCount

  function galaxyGasMajorMergerTimeExtract(self,node,instance) result(timeMajorMergers)
    !!{
    Implement a galaxyGasMajorMergerTime output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                                               , dimension(:,:), allocatable :: timeMajorMergers
    class           (nodePropertyExtractorGalaxyGasMajorMergerTime), intent(inout)               :: self
    type            (treeNode                                     ), intent(inout)               :: node
    type            (multiCounter                                 ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic                           )                , pointer     :: basic
    double precision                                               , dimension(:  ), allocatable :: times
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: times

    basic => node %basic                    (                            )
    times =  basic%floatRank1MetaPropertyGet(self%galaxyGasMajorMergerTimeID)
    allocate(timeMajorMergers(size(times),1))
    timeMajorMergers(:,1)=times
    return
  end function galaxyGasMajorMergerTimeExtract
  
  subroutine galaxyGasMajorMergerTimeNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily galaxyGasMajorMergerTime} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyGasMajorMergerTime), intent(inout)                             :: self
    type (varying_string                               ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('galaxyGasMajorMergerTime')
    return
  end subroutine galaxyGasMajorMergerTimeNames

  subroutine galaxyGasMajorMergerTimeDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily galaxyGasMajorMergerTime} properties.
    !!}
    implicit none
    class(nodePropertyExtractorGalaxyGasMajorMergerTime), intent(inout)                             :: self
    type (varying_string                               ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Times of gas-mass-based galaxy major mergers.')
    return
  end subroutine galaxyGasMajorMergerTimeDescriptions

  function galaxyGasMajorMergerTimeUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily galaxyGasMajorMergerTime} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    double precision                                               , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorGalaxyGasMajorMergerTime), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=gigaYear
    return
  end function galaxyGasMajorMergerTimeUnitsInSI
