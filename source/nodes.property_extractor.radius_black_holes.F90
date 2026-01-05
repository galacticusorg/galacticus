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

  use :: Black_Hole_Binary_Separations, only : blackHoleBinarySeparationGrowthRateClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusBlackHoles">
   <description>
     A node property extractor which extracts a list of all super-massive black hole radii and radial migration rates.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorRadiusBlackHoles
     !!{
     A property extractor which extracts a list of all super-massive black hole radii and radial migration rates.
     !!}
     private
     class(blackHoleBinarySeparationGrowthRateClass), pointer :: blackHoleBinarySeparationGrowthRate_ => null()
   contains
     procedure :: elementCount => radiusBlackHolesElementCount
     procedure :: extract      => radiusBlackHolesExtract
     procedure :: names        => radiusBlackHolesNames
     procedure :: descriptions => radiusBlackHolesDescriptions
     procedure :: unitsInSI    => radiusBlackHolesUnitsInSI
  end type nodePropertyExtractorRadiusBlackHoles

  interface nodePropertyExtractorRadiusBlackHoles
     !!{
     Constructors for the \refClass{nodePropertyExtractorRadiusBlackHoles} output extractor class.
     !!}
    module procedure radiusBlackHolesConstructorParameters
  end interface nodePropertyExtractorRadiusBlackHoles

contains

  function radiusBlackHolesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorRadiusBlackHoles} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRadiusBlackHoles   )                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(blackHoleBinarySeparationGrowthRateClass), pointer       :: blackHoleBinarySeparationGrowthRate_

    !![
    <objectBuilder class="blackHoleBinarySeparationGrowthRate" name="blackHoleBinarySeparationGrowthRate_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusBlackHoles(blackHoleBinarySeparationGrowthRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleBinarySeparationGrowthRate_"/>
    !!]
    return
  end function radiusBlackHolesConstructorParameters

  function radiusBlackHolesConstructorInternal(blackHoleBinarySeparationGrowthRate_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorRadiusBlackHoles} node operator class.
    !!}
    implicit none
    type (nodePropertyExtractorRadiusBlackHoles   )                        :: self
    class(blackHoleBinarySeparationGrowthRateClass), intent(in   ), target :: blackHoleBinarySeparationGrowthRate_
    !![
    <constructorAssign variables="*blackHoleBinarySeparationGrowthRate_"/>
    !!]

    return
  end function radiusBlackHolesConstructorInternal

  subroutine radiusBlackHolesDestructor(self)
    !!{
    Destructor for the critical overdensity radiusBlackHoles set barrier class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusBlackHoles), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleBinarySeparationGrowthRate_"/>
    !!]                                                                                                                                                                                                               
    return
  end subroutine radiusBlackHolesDestructor

  integer function radiusBlackHolesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorRadiusBlackHoles), intent(inout) :: self

    radiusBlackHolesElementCount=2
    return
  end function radiusBlackHolesElementCount

  function radiusBlackHolesExtract(self,node,instance) result(radius)
    !!{
    Implement an output extractor for the masses of all supermassive black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    double precision                                       , dimension(:,:), allocatable :: radius
    class           (nodePropertyExtractorRadiusBlackHoles), intent(inout)               :: self
    type            (treeNode                             ), intent(inout)               :: node
    type            (multiCounter                         ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole               )                , pointer     :: blackHole
    integer                                                                              :: i        , countBlackHoles
    !$GLC attributes unused :: instance

    countBlackHoles=node%blackHoleCount()
    allocate(radius(countBlackHoles,2))
    do i=1,countBlackHoles
       blackHole      => node     %blackHole                                      (instance=i        )
       radius   (i,1) =  blackHole%radialPosition                                 (                  )
       radius   (i,2) =  self     %blackHoleBinarySeparationGrowthRate_%growthRate(         blackHole)
    end do
    return
  end function radiusBlackHolesExtract

  subroutine radiusBlackHolesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily radiusBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorRadiusBlackHoles), intent(inout)                             :: self
    type (varying_string                       ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(2))
    names(1)=var_str('radiusBlackHoles'             )
    names(2)=var_str('radialMigrationRateBlackHoles')
    return
  end subroutine radiusBlackHolesNames

  subroutine radiusBlackHolesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily radiusBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorRadiusBlackHoles), intent(inout)                             :: self
    type (varying_string                       ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(2))
    descriptions(1)=var_str('Radial positions of super-massive black holes in this galaxy.'      )
    descriptions(2)=var_str('Radial migration rates of super-massive black holes in this galaxy.')
    return
  end subroutine radiusBlackHolesDescriptions

  function radiusBlackHolesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily radiusBlackHoles} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec, gigaYear
    implicit none
    double precision                                       , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorRadiusBlackHoles), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(2))
    unitsInSI(1)=megaParsec
    unitsInSI(2)=megaParsec/gigaYear
    return
  end function radiusBlackHolesUnitsInSI
