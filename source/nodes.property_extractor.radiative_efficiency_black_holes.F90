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

  use :: Accretion_Disks           , only : accretionDisksClass
  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiativeEfficiencyBlackHoles">
   <description>
     A node property extractor which extracts a list of all super-massive black hole radiative efficiencies.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorRadiativeEfficiencyBlackHoles
     !!{
     A property extractor which extracts a list of all super-massive black hole radiative efficiencies.
     !!}
     private
     class(accretionDisksClass        ), pointer :: accretionDisks_         => null()
     class(blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
   contains
     final     ::                 radiativeEfficiencyBlackHolesDestructor
     procedure :: elementCount => radiativeEfficiencyBlackHolesElementCount
     procedure :: extract      => radiativeEfficiencyBlackHolesExtract
     procedure :: names        => radiativeEfficiencyBlackHolesNames
     procedure :: descriptions => radiativeEfficiencyBlackHolesDescriptions
     procedure :: unitsInSI    => radiativeEfficiencyBlackHolesUnitsInSI
  end type nodePropertyExtractorRadiativeEfficiencyBlackHoles

  interface nodePropertyExtractorRadiativeEfficiencyBlackHoles
     !!{
     Constructors for the {\normalfont \ttfamily radiativeEfficiencyBlackHoles} output extractor class.
     !!}
     module procedure radiativeEfficiencyBlackHolesConstructorParameters
     module procedure radiativeEfficiencyBlackHolesConstructorInternal
  end interface nodePropertyExtractorRadiativeEfficiencyBlackHoles

contains

  function radiativeEfficiencyBlackHolesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiativeEfficiencyBlackHoles} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRadiativeEfficiencyBlackHoles)                :: self
    type (inputParameters                                   ), intent(inout) :: parameters
    class(accretionDisksClass                               ), pointer       :: accretionDisks_
    class(blackHoleAccretionRateClass                       ), pointer       :: blackHoleAccretionRate_

    !![
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiativeEfficiencyBlackHoles(blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function radiativeEfficiencyBlackHolesConstructorParameters

  function radiativeEfficiencyBlackHolesConstructorInternal(blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily radiativeEfficiencyBlackHoles} node operator class.
    !!}
    implicit none
    type (nodePropertyExtractorRadiativeEfficiencyBlackHoles)                        :: self
    class(accretionDisksClass                               ), intent(in   ), target :: accretionDisks_
    class(blackHoleAccretionRateClass                       ), intent(in   ), target :: blackHoleAccretionRate_
    !![
    <constructorAssign variables="*blackHoleAccretionRate_, *accretionDisks_"/>
    !!]
    
    return
  end function radiativeEfficiencyBlackHolesConstructorInternal

  subroutine radiativeEfficiencyBlackHolesDestructor(self)
    !!{
    Destructor for the critical overdensity radiativeEfficiencyBlackHoles set barrier class.
    !!}
    implicit none
    type(nodePropertyExtractorRadiativeEfficiencyBlackHoles), intent(inout) :: self
    
    !![
    <objectDestructor name="self%accretionDisks_"        />
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    !!]                                                                                                                                                                                                               
    return
  end subroutine radiativeEfficiencyBlackHolesDestructor
  
  integer function radiativeEfficiencyBlackHolesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorRadiativeEfficiencyBlackHoles), intent(inout) :: self

    radiativeEfficiencyBlackHolesElementCount=1
    return
  end function radiativeEfficiencyBlackHolesElementCount

  function radiativeEfficiencyBlackHolesExtract(self,node,instance) result(radiativeEfficiency)
    !!{
    Implement an output extractor for the radiative efficiencies of all supermassive black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    double precision                                                    , dimension(:,:), allocatable :: radiativeEfficiency
    class           (nodePropertyExtractorRadiativeEfficiencyBlackHoles), intent(inout)               :: self
    type            (treeNode                                          ), intent(inout)               :: node
    type            (multiCounter                                      ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole                            )                , pointer     :: blackHole
    integer                                                                                           :: i                                  , countBlackHoles
    double precision                                                                                  :: rateMassAccretionSpheroid          , rateMassAccretionHotHalo, &
        &                                                                                                rateMassAccretionNuclearStarCluster    
    !$GLC attributes unused :: instance

    countBlackHoles=node%blackHoleCount()
    allocate(radiativeEfficiency(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole                => node     %blackHole(instance=i)
       call  self%blackHoleAccretionRate_%rateAccretion(blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo,rateMassAccretionNuclearStarCluster)
       radiativeEfficiency(i,1) =  self%accretionDisks_%efficiencyRadiative(blackHole,rateMassAccretionSpheroid+rateMassAccretionHotHalo+rateMassAccretionNuclearStarCluster)
    end do
    return
  end function radiativeEfficiencyBlackHolesExtract
  
  subroutine radiativeEfficiencyBlackHolesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily radiativeEfficiencyBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorRadiativeEfficiencyBlackHoles), intent(inout)                             :: self
    type (varying_string                                    ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('radiativeEfficiencyBlackHoles')
    return
  end subroutine radiativeEfficiencyBlackHolesNames

  subroutine radiativeEfficiencyBlackHolesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily radiativeEfficiencyBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorRadiativeEfficiencyBlackHoles), intent(inout)                             :: self
    type (varying_string                                    ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Radiative efficiencies of super-massive black holes in this galaxy.')
    return
  end subroutine radiativeEfficiencyBlackHolesDescriptions

  function radiativeEfficiencyBlackHolesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily radiativeEfficiencyBlackHoles} properties in the SI system.
    !!}
    implicit none
    double precision                                                    , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorRadiativeEfficiencyBlackHoles), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=1.0d0
    return
  end function radiativeEfficiencyBlackHolesUnitsInSI
