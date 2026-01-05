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
  Implements the a black hole accretion rate model that tracks the growth of the spheroid.
  !!}

  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass

  !![
  <blackHoleAccretionRate name="blackHoleAccretionRateSpheroidTracking">
   <description>
    A black hole accretion rate calculation that tracks the growth of the spheroid.
   </description>
  </blackHoleAccretionRate>
  !!]
  type, extends(blackHoleAccretionRateClass) :: blackHoleAccretionRateSpheroidTracking
     !!{
     The spheroidTracking black hole accretion rate calculation.      
     !!}
     private
     class           (starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_  => null()
     double precision                                           :: growthRatioToStellarSpheroid
   contains
     final     ::                  spheroidTrackingDestructor
     procedure :: rateAccretion => spheroidTrackingRateAccretion
  end type blackHoleAccretionRateSpheroidTracking

  interface blackHoleAccretionRateSpheroidTracking
     !!{
     Constructors for the \refClass{blackHoleAccretionRateSpheroidTracking} black hole accretion rate class.
     !!}
     module procedure spheroidTrackingConstructorParameters
     module procedure spheroidTrackingConstructorInternal
  end interface blackHoleAccretionRateSpheroidTracking

contains

  function spheroidTrackingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleAccretionRateSpheroidTracking} black hole accretion rate class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleAccretionRateSpheroidTracking)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (starFormationRateSpheroidsClass       ), pointer       :: starFormationRateSpheroids_
    double precision                                                        :: growthRatioToStellarSpheroid

    !![
    <inputParameter>
      <name>growthRatioToStellarSpheroid</name>
      <defaultValue>1.0d-3</defaultValue>
      <description>The ratio of the rates of black hole growth and spheroid stellar mass growth.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    !!]
    self=blackHoleAccretionRateSpheroidTracking(growthRatioToStellarSpheroid,starFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateSpheroids_"/>
    !!]
    return
  end function spheroidTrackingConstructorParameters

  function spheroidTrackingConstructorInternal(growthRatioToStellarSpheroid,starFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the \refClass{blackHoleAccretionRateSpheroidTracking} node operator class.
    !!}
    implicit none
    type            (blackHoleAccretionRateSpheroidTracking)                        :: self
    class           (starFormationRateSpheroidsClass       ), intent(in   ), target :: starFormationRateSpheroids_
    double precision                                        , intent(in   )         :: growthRatioToStellarSpheroid
    !![
    <constructorAssign variables="growthRatioToStellarSpheroid, *starFormationRateSpheroids_"/>
    !!]

    return
  end function spheroidTrackingConstructorInternal

  subroutine spheroidTrackingDestructor(self)
    !!{
    Destructor for the critical overdensity spheroidTracking set barrier class.
    !!}
    implicit none
    type(blackHoleAccretionRateSpheroidTracking), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateSpheroids_"/>
    !!]                                                                                                                                                                                                               
    return
  end subroutine spheroidTrackingDestructor

  subroutine spheroidTrackingRateAccretion(self,blackHole,rateMassAccretionSpheroid,rateMassAccretionHotHalo,rateMassAccretionNuclearStarCluster)
    !!{
    Compute the accretion rate onto a black hole.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class           (blackHoleAccretionRateSpheroidTracking), intent(inout) :: self
    class           (nodeComponentBlackHole                ), intent(inout) :: blackHole
    double precision                                        , intent(  out) :: rateMassAccretionSpheroid          , rateMassAccretionHotHalo, &
         &                                                                     rateMassAccretionNuclearStarCluster
    type            (treeNode                              ), pointer       :: node

    node                                =>  blackHole                            %host                        (    )
    rateMassAccretionSpheroid           =  +self                                 %growthRatioToStellarSpheroid       &
         &                                 *self     %starFormationRateSpheroids_%rate                        (node)
    rateMassAccretionHotHalo            =  +0.0d0
    rateMassAccretionNuclearStarCluster =  +0.0d0
    return
  end subroutine spheroidTrackingRateAccretion
