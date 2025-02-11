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

  !!{
  Implements a black hole CGM heating class using the accretion disk jet power.
  !!}

  use :: Black_Hole_Accretion_Rates, only : blackHoleAccretionRateClass
  use :: Accretion_Disks           , only : accretionDisksClass

  !![
  <blackHoleCGMHeating name="blackHoleCGMHeatingJetPower">
   <description>
    A black hole CGM heating class using the accretion disk jet power.
   </description>
  </blackHoleCGMHeating>
  !!]
  type, extends(blackHoleCGMHeatingClass) :: blackHoleCGMHeatingJetPower
     !!{
     A black hole CGM heating class using the accretion disk jet power.
     !!}
     private
     class           (blackHoleAccretionRateClass), pointer :: blackHoleAccretionRate_ => null()
     class           (accretionDisksClass        ), pointer :: accretionDisks_         => null()
     double precision                                       :: efficiencyRadioMode
   contains
     final     ::                jetPowerDestructor
     procedure :: heatingRate => jetPowerHeatingRate
  end type blackHoleCGMHeatingJetPower
  
  interface blackHoleCGMHeatingJetPower
     !!{
     Constructors for the {\normalfont \ttfamily jetPower} black hole winds class.
     !!}
     module procedure jetPowerConstructorParameters
     module procedure jetPowerConstructorInternal
  end interface blackHoleCGMHeatingJetPower

contains

  function jetPowerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily jetPower} black hole winds class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleCGMHeatingJetPower)                :: self
    type            (inputParameters            ), intent(inout) :: parameters
    class           (blackHoleAccretionRateClass), pointer       :: blackHoleAccretionRate_
    class           (accretionDisksClass        ), pointer       :: accretionDisks_
    double precision                                             :: efficiencyRadioMode
    
    !![
    <inputParameter>
      <name>efficiencyRadioMode</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Efficiency with which radio-mode feedback is coupled to the \gls{cgm}.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="blackHoleAccretionRate" name="blackHoleAccretionRate_" source="parameters"/>
    <objectBuilder class="accretionDisks"         name="accretionDisks_"         source="parameters"/>
    !!]
    self=blackHoleCGMHeatingJetPower(efficiencyRadioMode,blackHoleAccretionRate_,accretionDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="blackHoleAccretionRate_"/>
    <objectDestructor name="accretionDisks_"        />
    !!]
    return
  end function jetPowerConstructorParameters

  function jetPowerConstructorInternal(efficiencyRadioMode,blackHoleAccretionRate_,accretionDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily jetPower} node operator class.
    !!}
    implicit none
    type            (blackHoleCGMHeatingJetPower)                        :: self
    class           (blackHoleAccretionRateClass), target, intent(in   ) :: blackHoleAccretionRate_
    class           (accretionDisksClass        ), target, intent(in   ) :: accretionDisks_
    double precision                                     , intent(in   ) :: efficiencyRadioMode
    !![
    <constructorAssign variables="efficiencyRadioMode, *blackHoleAccretionRate_, *accretionDisks_"/>
    !!]

    return
  end function jetPowerConstructorInternal

  subroutine jetPowerDestructor(self)
    !!{
    Destructor for the jetPower black hole CGM heating class.
    !!}
    implicit none
    type(blackHoleCGMHeatingJetPower), intent(inout) :: self
    
    !![
    <objectDestructor name="self%blackHoleAccretionRate_"/>
    <objectDestructor name="self%accretionDisks_"        />
    !!]
    return
  end subroutine jetPowerDestructor
  
  double precision function jetPowerHeatingRate(self,blackHole) result(rateHeating)
    !!{
    Compute the heating rate of the CGM based on the accretion disk jet power.
    !!}
    implicit none
    class           (blackHoleCGMHeatingJetPower), intent(inout) :: self
    class           (nodeComponentBlackHole     ), intent(inout) :: blackHole
    double precision                                             :: rateAccretionSpheroid          , rateAccretionHotHalo, &
         &                                                          rateAccretionNuclearStarCluster, rateAccretion

    ! Compute the jet power and CGM heating rate.
    call self%blackHoleAccretionRate_%rateAccretion(blackHole,rateAccretionSpheroid,rateAccretionHotHalo,rateAccretionNuclearStarCluster)
    rateAccretion=+rateAccretionSpheroid           &
         &        +rateAccretionHotHalo            &
         &        +rateAccretionNuclearStarCluster
    rateHeating  =+self                %efficiencyRadioMode                          &
         &        *self%accretionDisks_%powerJet           (blackHole,rateAccretion)
    return
  end function jetPowerHeatingRate
