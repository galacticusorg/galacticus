!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  <nodePropertyExtractor name="nodePropertyExtractorGalaxyMajorMergerTime">
   <description>
     A node property extractor which extracts the time of the last major merger for each galaxy.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorGalaxyMajorMergerTime
     !!{
     A property extractor which extracts the time of the last major merger for each galaxy.
     !!}
     private
     integer :: galaxyMajorMergerTimeID
   contains
     procedure :: extract     => galaxyMajorMergerTimeExtract
     procedure :: name        => galaxyMajorMergerTimeName
     procedure :: description => galaxyMajorMergerTimeDescription
     procedure :: unitsInSI   => galaxyMajorMergerTimeUnitsInSI
  end type nodePropertyExtractorGalaxyMajorMergerTime

  interface nodePropertyExtractorGalaxyMajorMergerTime
     !!{
     Constructors for the ``galaxyMajorMergerTime'' output extractor class.
     !!}
     module procedure galaxyMajorMergerTimeConstructorParameters
     module procedure galaxyMajorMergerTimeConstructorInternal
  end interface nodePropertyExtractorGalaxyMajorMergerTime

contains

  function galaxyMajorMergerTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``galaxyMajorMergerTime'' property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorGalaxyMajorMergerTime)                :: self
    type(inputParameters                           ), intent(inout) :: parameters

    self=nodePropertyExtractorGalaxyMajorMergerTime()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function galaxyMajorMergerTimeConstructorParameters

  function galaxyMajorMergerTimeConstructorInternal() result(self)
    !!{
    Internal constructor for the ``galaxyMajorMergerTime'' output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorGalaxyMajorMergerTime) :: self
    
    !![
    <addMetaProperty component="basic" name="galaxyMajorMergerTime" id="self%galaxyMajorMergerTimeID" rank="1" isCreator="no"/>
    !!]
    return
  end function galaxyMajorMergerTimeConstructorInternal

  function galaxyMajorMergerTimeExtract(self,node,instance) result(timeMajorMergers)
    !!{
    Implement a galaxyMajorMergerTime output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                                            , dimension(:) , allocatable :: timeMajorMergers
    class           (nodePropertyExtractorGalaxyMajorMergerTime), intent(inout)              :: self
    type            (treeNode                                  ), intent(inout)              :: node
    type            (multiCounter                              ), intent(inout), optional    :: instance
    class           (nodeComponentBasic                        )               , pointer     :: basic
    !$GLC attributes unused :: instance

    basic            => node %basic                    (                            )
    timeMajorMergers =  basic%floatRank1MetaPropertyGet(self%galaxyMajorMergerTimeID)
    return
  end function galaxyMajorMergerTimeExtract
  
  function galaxyMajorMergerTimeName(self)
    !!{
    Return the names of the {\normalfont \ttfamily galaxyMajorMergerTime} properties.
    !!}
    implicit none
    type (varying_string                            )                :: galaxyMajorMergerTimeName
    class(nodePropertyExtractorGalaxyMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    galaxyMajorMergerTimeName=var_str('galaxyMajorMergerTime')
    return
  end function galaxyMajorMergerTimeName

  function galaxyMajorMergerTimeDescription(self)
    !!{
    Return the descriptions of the {\normalfont \ttfamily galaxyMajorMergerTime} properties.
    !!}
    implicit none
    type (varying_string                            )                :: galaxyMajorMergerTimeDescription
    class(nodePropertyExtractorGalaxyMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    galaxyMajorMergerTimeDescription=var_str('Time of the last galaxy major merger.')
    return
  end function galaxyMajorMergerTimeDescription

  double precision function galaxyMajorMergerTimeUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily galaxyMajorMergerTime} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorGalaxyMajorMergerTime), intent(inout) :: self
    !$GLC attributes unused :: self

    galaxyMajorMergerTimeUnitsInSI=gigaYear
    return
  end function galaxyMajorMergerTimeUnitsInSI
