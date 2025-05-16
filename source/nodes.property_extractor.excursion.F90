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
  <nodePropertyExtractor name="nodePropertyExtractorExcursion">
   <description>
     A node property extractor which extracts the (infimum of the) excursion corresponding to the mass accretion history for each node.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorExcursion
     !!{
     A property extractor which extracts the (infimum of the) excursion corresponding to the mass accretion history for each node.
     !!}
     private
     integer :: excursionOverdensityID, excursionVarianceID
   contains
     procedure :: elementCount => excursionElementCount
     procedure :: extract      => excursionExtract
     procedure :: names        => excursionNames
     procedure :: descriptions => excursionDescriptions
     procedure :: unitsInSI    => excursionUnitsInSI
  end type nodePropertyExtractorExcursion

  interface nodePropertyExtractorExcursion
     !!{
     Constructors for the \refClass{nodePropertyExtractorExcursion} output extractor class.
     !!}
     module procedure excursionConstructorParameters
     module procedure excursionConstructorInternal
  end interface nodePropertyExtractorExcursion

contains

  function excursionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorExcursion} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(nodePropertyExtractorExcursion)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=nodePropertyExtractorExcursion()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function excursionConstructorParameters

  function excursionConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorExcursion} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorExcursion) :: self
    
    !![
    <addMetaProperty component="basic" name="excursionTime" id="self%excursionOverdensityID" rank="1" isCreator="no"/>
    <addMetaProperty component="basic" name="excursionMass" id="self%excursionVarianceID"    rank="1" isCreator="no"/>
    !!]
    return
  end function excursionConstructorInternal

  integer function excursionElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorExcursion), intent(inout) :: self

    excursionElementCount=2
    return
  end function excursionElementCount

  function excursionExtract(self,node,instance) result(excursion)
    !!{
    Implement a excursion output extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                                , dimension(:,:), allocatable :: excursion
    class           (nodePropertyExtractorExcursion), intent(inout)               :: self
    type            (treeNode                      ), intent(inout)               :: node
    type            (multiCounter                  ), intent(inout) , optional    :: instance
    class           (nodeComponentBasic            )                , pointer     :: basic
    double precision                                , dimension(:  ), allocatable :: overdensities, variances
    !$GLC attributes unused :: instance
    !$GLC attributes initialized :: overdensities, variances
    
    basic         => node %basic                    (                           )
    overdensities =  basic%floatRank1MetaPropertyGet(self%excursionOverdensityID)
    variances     =  basic%floatRank1MetaPropertyGet(self%excursionVarianceID   )
    allocate(excursion(size(overdensities),2))
    excursion(:,1)=overdensities
    excursion(:,2)=variances
    return
  end function excursionExtract
  
  subroutine excursionNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily excursion} properties.
    !!}
    implicit none
    class(nodePropertyExtractorExcursion), intent(inout)                             :: self
    type (varying_string                ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(2))
    names(1)=var_str('haloExcursionOverdensity')
    names(2)=var_str('haloExcursionVariance'   )
    return
  end subroutine excursionNames

  subroutine excursionDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily excursion} properties.
    !!}
    implicit none
    class(nodePropertyExtractorExcursion), intent(inout)                             :: self
    type (varying_string                ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(2))
    descriptions(1)=var_str("The overdensity in the infimum of the halo's excursion.")
    descriptions(2)=var_str("The variance in the infimum of the halo's excursion."   )
    return
  end subroutine excursionDescriptions

  function excursionUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily excursion} properties in the SI system.
    !!}
    implicit none
    double precision                                , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorExcursion), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(2))
    unitsInSI=1.0d0
    return
  end function excursionUnitsInSI
