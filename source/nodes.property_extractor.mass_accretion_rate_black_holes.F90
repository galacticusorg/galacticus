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

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorMassAccretionRateBlackHoles">
   <description>
     A node property extractor which extracts a list of all super-massive black hole mass accretion rates.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorMassAccretionRateBlackHoles
     !!{
     A property extractor which extracts a list of all super-massive black hole mass accretion rates.
     !!}
     private
     class(blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
   contains
     final     ::                 massAccretionRateBlackHolesDestructor
     procedure :: elementCount => massAccretionRateBlackHolesElementCount
     procedure :: extract      => massAccretionRateBlackHolesExtract
     procedure :: names        => massAccretionRateBlackHolesNames
     procedure :: descriptions => massAccretionRateBlackHolesDescriptions
     procedure :: unitsInSI    => massAccretionRateBlackHolesUnitsInSI
  end type nodePropertyExtractorMassAccretionRateBlackHoles

  interface nodePropertyExtractorMassAccretionRateBlackHoles
     !!{
     Constructors for the ``massAccretionRateBlackHoles'' output extractor class.
     !!}
     module procedure massAccretionRateBlackHolesConstructorParameters
     module procedure massAccretionRateBlackHolesConstructorInternal
  end interface nodePropertyExtractorMassAccretionRateBlackHoles

contains

  function massAccretionRateBlackHolesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``massAccretionRateBlackHoles'' property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorMassAccretionRateBlackHoles)                :: self
    type (inputParameters                                 ), intent(inout) :: parameters
    class(blackHoleAccretionRateClass                     ), pointer       :: blackHoleAccretionRate_

    !![
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    !!]
    self=nodePropertyExtractorMassAccretionRateBlackHoles(blackHoleAccretionRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    !!]
    return
  end function massAccretionRateBlackHolesConstructorParameters

  function massAccretionRateBlackHolesConstructorInternal(blackHoleAccretionRate_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily massAccretionRateBlackHoles} node operator class.
    !!}
    implicit none
    type (nodePropertyExtractorMassAccretionRateBlackHoles)                        :: self
    class(blackHoleAccretionRateClass                     ), intent(in   ), target :: blackHoleAccretionRate_
    !![
    <constructorAssign variables="*blackHoleAccretionRate_"/>
    !!]

    return
  end function massAccretionRateBlackHolesConstructorInternal

  subroutine massAccretionRateBlackHolesDestructor(self)
    !!{
    Destructor for the critical overdensity massAccretionRateBlackHoles set barrier class.
    !!}
    implicit none
    type(nodePropertyExtractorMassAccretionRateBlackHoles), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    !!]                                                                                                                                                                                                               
    return
  end subroutine massAccretionRateBlackHolesDestructor

  integer function massAccretionRateBlackHolesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorMassAccretionRateBlackHoles), intent(inout) :: self

    massAccretionRateBlackHolesElementCount=1
    return
  end function massAccretionRateBlackHolesElementCount

  function massAccretionRateBlackHolesExtract(self,node,instance) result(massAccretionRate)
    !!{
    Implement an output extractor for the radiative efficiencies of all supermassive black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    double precision                                                  , dimension(:,:), allocatable :: massAccretionRate
    class           (nodePropertyExtractorMassAccretionRateBlackHoles), intent(inout)               :: self
    type            (treeNode                                        ), intent(inout)               :: node
    type            (multiCounter                                    ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole                          )                , pointer     :: blackHole
    integer                                                                                         :: i                                  , countBlackHoles
    double precision                                                                                :: rateMassAccretionSpheroid          , rateMassAccretionHotHalo, &
        &                                                                                              rateMassAccretionNuclearStarCluster
    !$GLC attributes unused :: instance

    countBlackHoles=node%blackHoleCount()
    allocate(massAccretionRate(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole                =>  node%blackHole(instance=i)
       call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo,rateMassAccretionNuclearStarCluster)
       massAccretionRate  (i,1) =  +rateMassAccretionSpheroid           &
            &                      +rateMassAccretionHotHalo            &
            &                      +rateMassAccretionNuclearStarCluster  
    end do
    return
  end function massAccretionRateBlackHolesExtract

  subroutine massAccretionRateBlackHolesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily massAccretionRateBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorMassAccretionRateBlackHoles), intent(inout)                             :: self
    type (varying_string                                  ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('massAccretionRateBlackHoles')
    return
  end subroutine massAccretionRateBlackHolesNames

  subroutine massAccretionRateBlackHolesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily massAccretionRateBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorMassAccretionRateBlackHoles), intent(inout)                             :: self
    type (varying_string                                  ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Mass accretion rates of super-massive black holes in this galaxy.')
    return
  end subroutine massAccretionRateBlackHolesDescriptions

  function massAccretionRateBlackHolesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily massAccretionRateBlackHoles} properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    double precision                                                  , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorMassAccretionRateBlackHoles), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=massSolar/gigaYear
    return
  end function massAccretionRateBlackHolesUnitsInSI
