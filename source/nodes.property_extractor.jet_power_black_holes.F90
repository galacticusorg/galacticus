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
  <nodePropertyExtractor name="nodePropertyExtractorJetPowerBlackHoles">
   <description>
     A node property extractor which extracts a list of all super-massive jet powers.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorList) :: nodePropertyExtractorJetPowerBlackHoles
     !!{
     A property extractor which extracts a list of all super-massive black hole jet powers.
     !!}
     private
     class(accretionDisksClass        ), pointer :: accretionDisks_         => null()
     class(blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
   contains
     final     ::                 jetPowerBlackHolesDestructor
     procedure :: elementCount => jetPowerBlackHolesElementCount
     procedure :: extract      => jetPowerBlackHolesExtract
     procedure :: names        => jetPowerBlackHolesNames
     procedure :: descriptions => jetPowerBlackHolesDescriptions
     procedure :: unitsInSI    => jetPowerBlackHolesUnitsInSI
  end type nodePropertyExtractorJetPowerBlackHoles

  interface nodePropertyExtractorJetPowerBlackHoles
     !!{
     Constructors for the {\normalfont \ttfamily jetPowerBlackHoles} output extractor class.
     !!}
     module procedure jetPowerBlackHolesConstructorParameters
     module procedure jetPowerBlackHolesConstructorInternal
  end interface nodePropertyExtractorJetPowerBlackHoles

contains

  function jetPowerBlackHolesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily jetPowerBlackHoles} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorJetPowerBlackHoles)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(accretionDisksClass                    ), pointer       :: accretionDisks_
    class(blackHoleAccretionRateClass            ), pointer       :: blackHoleAccretionRate_

    !![
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    !!]
    self=nodePropertyExtractorJetPowerBlackHoles(blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function jetPowerBlackHolesConstructorParameters

  function jetPowerBlackHolesConstructorInternal(blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily jetPowerBlackHoles} node operator class.
    !!}
    implicit none
    type (nodePropertyExtractorJetPowerBlackHoles)                        :: self
    class(accretionDisksClass                    ), intent(in   ), target :: accretionDisks_
    class(blackHoleAccretionRateClass            ), intent(in   ), target :: blackHoleAccretionRate_
    !![
    <constructorAssign variables="*blackHoleAccretionRate_, *accretionDisks_"/>
    !!]

    return
  end function jetPowerBlackHolesConstructorInternal

  subroutine jetPowerBlackHolesDestructor(self)
    !!{
    Destructor for the critical overdensity jetPowerBlackHoles set barrier class.
    !!}
    implicit none
    type(nodePropertyExtractorJetPowerBlackHoles), intent(inout) :: self

    !![
    <objectDestructor name="self%accretionDisks_"        />
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    !!]                                                                                                                                                                                                               
    return
  end subroutine jetPowerBlackHolesDestructor

  integer function jetPowerBlackHolesElementCount(self)
    !!{
    Return a count of the number of properties extracted.
    !!}
    implicit none
    class(nodePropertyExtractorJetPowerBlackHoles), intent(inout) :: self

    jetPowerBlackHolesElementCount=1
    return
  end function jetPowerBlackHolesElementCount

  function jetPowerBlackHolesExtract(self,node,instance) result(radiativeEfficiency)
    !!{
    Implement an output extractor for the radiative efficiencies of all supermassive black holes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole
    implicit none
    double precision                                         , dimension(:,:), allocatable :: radiativeEfficiency
    class           (nodePropertyExtractorJetPowerBlackHoles), intent(inout)               :: self
    type            (treeNode                               ), intent(inout)               :: node
    type            (multiCounter                           ), intent(inout) , optional    :: instance
    class           (nodeComponentBlackHole                 )                , pointer     :: blackHole
    integer                                                                                :: i                                  , countBlackHoles
    double precision                                                                       :: rateMassAccretionSpheroid          , rateMassAccretionHotHalo, &
        &                                                                                     rateMassAccretionNuclearStarCluster
    !$GLC attributes unused :: instance

    countBlackHoles=node%blackHoleCount()
    allocate(radiativeEfficiency(countBlackHoles,1))
    do i=1,countBlackHoles
       blackHole                => node     %blackHole(instance=i)
       call  self%blackHoleAccretionRate_%rateAccretion(blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo,rateMassAccretionNuclearStarCluster)
       radiativeEfficiency(i,1) =  self%accretionDisks_%powerJet(blackHole,rateMassAccretionSpheroid+rateMassAccretionHotHalo+rateMassAccretionNuclearStarCluster)
    end do
    return
  end function jetPowerBlackHolesExtract

  subroutine jetPowerBlackHolesNames(self,names)
    !!{
    Return the names of the {\normalfont \ttfamily jetPowerBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorJetPowerBlackHoles), intent(inout)                             :: self
    type (varying_string                         ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self

    allocate(names(1))
    names(1)=var_str('powerJetBlackHoles')
    return
  end subroutine jetPowerBlackHolesNames

  subroutine jetPowerBlackHolesDescriptions(self,descriptions)
    !!{
    Return the descriptions of the {\normalfont \ttfamily jetPowerBlackHoles} properties.
    !!}
    implicit none
    class(nodePropertyExtractorJetPowerBlackHoles), intent(inout)                             :: self
    type (varying_string                         ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self

    allocate(descriptions(1))
    descriptions(1)=var_str('Jet power of super-massive black holes in this galaxy [M☉ (km/s)² Gyr¯¹].')
    return
  end subroutine jetPowerBlackHolesDescriptions

  function jetPowerBlackHolesUnitsInSI(self) result(unitsInSI)
    !!{
    Return the units of the {\normalfont \ttfamily jetPowerBlackHoles} properties in the SI system.
    !!}
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Astronomical, only : massSolar, gigaYear
    implicit none
    double precision                                         , dimension(:) , allocatable :: unitsInSI
    class           (nodePropertyExtractorJetPowerBlackHoles), intent(inout)              :: self
    !$GLC attributes unused :: self

    allocate(unitsInSI(1))
    unitsInSI(1)=+massSolar    &
         &       *kilo     **2 &
         &       /gigaYear
    return
  end function jetPowerBlackHolesUnitsInSI
