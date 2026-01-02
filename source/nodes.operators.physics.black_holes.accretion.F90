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
  Implements a node operator class that implements accretion onto black holes.
  !!}

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use :: Accretion_Disks           , only : accretionDisksClass

  !![
  <nodeOperator name="nodeOperatorBlackHolesAccretion">
   <description>A node operator class that implements accretion onto black holes.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBlackHolesAccretion
     !!{
     A node operator class that implements accretion onto black holes.
     !!}
     private
     class(blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
     class(accretionDisksClass        ), pointer :: accretionDisks_         => null()
   contains
     final     ::                          blackHolesAccretionDestructor
     procedure :: differentialEvolution => blackHolesAccretionDifferentialEvolution
  end type nodeOperatorBlackHolesAccretion
  
  interface nodeOperatorBlackHolesAccretion
     !!{
     Constructors for the \refClass{nodeOperatorBlackHolesAccretion} node operator class.
     !!}
     module procedure blackHolesAccretionConstructorParameters
     module procedure blackHolesAccretionConstructorInternal
  end interface nodeOperatorBlackHolesAccretion
  
contains

  function blackHolesAccretionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorBlackHolesAccretion} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBlackHolesAccretion)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(blackHoleAccretionRateClass    ), pointer       :: blackHoleAccretionRate_
    class(accretionDisksClass            ), pointer       :: accretionDisks_

    !![
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    !!]
    self=nodeOperatorBlackHolesAccretion(blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function blackHolesAccretionConstructorParameters

  function blackHolesAccretionConstructorInternal(blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorBlackHolesAccretion} node operator class.
    !!}
    implicit none
    type (nodeOperatorBlackHolesAccretion)                        :: self
    class(blackHoleAccretionRateClass    ), intent(in   ), target :: blackHoleAccretionRate_
    class(accretionDisksClass            ), intent(in   ), target :: accretionDisks_
    !![
    <constructorAssign variables="*blackHoleAccretionRate_, *accretionDisks_"/>
    !!]
    
    return
  end function blackHolesAccretionConstructorInternal

  subroutine blackHolesAccretionDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorBlackHolesAccretion} node operator class.
    !!}
    implicit none
    type(nodeOperatorBlackHolesAccretion), intent(inout) :: self

    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    <objectDestructor name="self%accretionDisks_"        />
    !!]
    return
  end subroutine blackHolesAccretionDestructor

  subroutine blackHolesAccretionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Account for accretion onto black holes.
    !!}
    use :: Galacticus_Nodes            , only : nodeComponentBlackHole, nodeComponentSpheroid, nodeComponentHotHalo, nodeComponentNSC
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class           (nodeOperatorBlackHolesAccretion), intent(inout), target  :: self
    type            (treeNode                       ), intent(inout), target  :: node
    logical                                          , intent(inout)          :: interrupt
    procedure       (interruptTask                  ), intent(inout), pointer :: functionInterrupt
    integer                                          , intent(in   )          :: propertyType
    class           (nodeComponentBlackHole         )               , pointer :: blackHole
    class           (nodeComponentSpheroid          )               , pointer :: spheroid
    class           (nodeComponentHotHalo           )               , pointer :: hotHalo
    class           (nodeComponentNSC               )               , pointer :: nuclearStarCluster
    integer                                                                   :: countBlackHole       , indexBlackHole
    double precision                                                          :: rateAccretionSpheroid, rateAccretionHotHalo           , &
         &                                                                       efficiencyRadiative  , efficiencyJet                  , &
         &                                                                       rateAccretion        , rateSpinUp                     , &
         &                                                                       rateAccretionReduced , rateAccretionNuclearStarCluster
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! If there are no black holes in this node, we have nothing to do - return immediately.
    countBlackHole=node%blackHoleCount()
    if (countBlackHole < 1) return
    ! Get spheroid, hot halo, and nuclear star cluster components so that accreted mass can be removed from them.
    spheroid           => node%spheroid()
    hotHalo            => node%hotHalo ()
    nuclearStarCluster => node%NSC     ()
    ! Iterate over all black holes in the node.
    do indexBlackHole=1,countBlackHole
       ! Find the accretion rate onto this black hole.
       blackHole => node%blackHole(instance=indexBlackHole)
       call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
       rateAccretion=+rateAccretionSpheroid           &
            &        +rateAccretionHotHalo            &
            &        +rateAccretionNuclearStarCluster
       ! Finish if there is no accretion.
       if (rateAccretion <= 0.0d0) cycle
       ! Find the radiative and jet efficiencies - these will be subtracted from the black hole mass growth rate.
       efficiencyRadiative=+self%accretionDisks_%efficiencyRadiative(blackHole,rateAccretion)
       efficiencyJet      =+self%accretionDisks_%powerJet           (blackHole,rateAccretion) &
            &              /                                                   rateAccretion  &
            &              /                     speedLight**2/kilo**2
       ! Find the rate of increase in mass of the black hole.
       rateAccretionReduced=rateAccretion*(1.0d0-efficiencyRadiative-efficiencyJet)
       ! Skip to the next black hole if this one has non-positive mass and a negative accretion rate.
       if (blackHole%mass() <= 0.0d0 .and. rateAccretionReduced < 0.0d0) cycle
       ! Find the spin-up rate for this black hole.
       rateSpinUp=+self%accretionDisks_%rateSpinUp(blackHole,rateAccretion)
       ! Accumulate rates of mass accretion/loss to the relevant components.
       call blackHole         %       massRate(+rateAccretionReduced                                       )
       call spheroid          %massGasSinkRate(-rateAccretionSpheroid                                      )
       call hotHalo           %   massSinkRate(-rateAccretionHotHalo           ,interrupt,functionInterrupt)
       call nuclearStarCluster%massGasSinkRate(-rateAccretionNuclearStarCluster                            )
       ! Set spin-up rate due to accretion.
       call blackHole         %       spinRate(+rateSpinUp                                                 )
     end do
    return
  end subroutine blackHolesAccretionDifferentialEvolution
  
