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
  Implements calculations of satellite merging times using the \cite{boylan-kolchin_dynamical_2008} method.
  !!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesBoylanKolchin2008">
   <description>
    A satellite merging timescale class which computes merging timescales using the dynamical friction calibration of
    \cite{boylan-kolchin_dynamical_2008}.
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesBoylanKolchin2008
     !!{
     A class implementing the \cite{boylan-kolchin_dynamical_2008} method for satellite merging timescales.
     !!}
     private
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     double precision                                     :: timescaleMultiplier
   contains
     final     ::                     boylanKolchin2008Destructor
     procedure :: timeUntilMerging => boylanKolchin2008TimeUntilMerging
  end type satelliteMergingTimescalesBoylanKolchin2008

  interface satelliteMergingTimescalesBoylanKolchin2008
     !!{
     Constructors for the \refClass{satelliteMergingTimescalesBoylanKolchin2008} satellite merging timescale class.
     !!}
     module procedure boylanKolchin2008ConstructorParameters
     module procedure boylanKolchin2008ConstructorInternal
  end interface satelliteMergingTimescalesBoylanKolchin2008

contains

  function boylanKolchin2008ConstructorParameters(parameters) result(self)
    !!{
    A constructor for the {\normalfont \ttfamily boylanKolchin2008} satellite merging timescale class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteMergingTimescalesBoylanKolchin2008)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                   ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                  ), pointer       :: darkMatterProfileDMO_
    double precision                                                             :: timescaleMultiplier

    !![
    <inputParameter>
      <name>timescaleMultiplier</name>
      <defaultValue>0.75d0</defaultValue>
      <description>A multiplier for the merging timescale in dynamical friction timescale calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=satelliteMergingTimescalesBoylanKolchin2008(timescaleMultiplier,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function boylanKolchin2008ConstructorParameters

  function boylanKolchin2008ConstructorInternal(timescaleMultiplier,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily boylanKolchin2008} satellite merging timescale class.
    !!}
    implicit none
    type            (satelliteMergingTimescalesBoylanKolchin2008)                        :: self
    double precision                                             , intent(in   )         :: timescaleMultiplier
    class           (darkMatterHaloScaleClass                   ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                  ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="timescaleMultiplier, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function boylanKolchin2008ConstructorInternal

  subroutine boylanKolchin2008Destructor(self)
    !!{
    Destructor for the \refClass{satelliteMergingTimescalesBoylanKolchin2008} satellite merging timescale class.
    !!}
    implicit none
    type(satelliteMergingTimescalesBoylanKolchin2008), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine boylanKolchin2008Destructor

  double precision function boylanKolchin2008TimeUntilMerging(self,node,orbit)
    !!{
    Return the timescale for merging satellites using the \cite{boylan-kolchin_dynamical_2008} method.
    !!}
    use :: Error             , only : Error_Report
    use :: Galacticus_Nodes  , only : nodeComponentBasic                              , treeNode
    use :: Mass_Distributions, only : massDistributionClass
    use :: Kepler_Orbits     , only : keplerOrbit
    use :: Satellite_Orbits  , only : Satellite_Orbit_Equivalent_Circular_Orbit_Radius, errorCodeNoEquivalentOrbit, errorCodeOrbitUnbound, errorCodeSuccess
    implicit none
    class           (satelliteMergingTimescalesBoylanKolchin2008), intent(inout) :: self
    type            (treeNode                                   ), intent(inout) :: node
    type            (keplerOrbit                                ), intent(inout) :: orbit
    type            (treeNode                                   ), pointer       :: nodeHost
    class           (nodeComponentBasic                         ), pointer       :: basicHost                            , basic
    class           (massDistributionClass                      ), pointer       :: massDistribution_
    logical                                                      , parameter     :: acceptUnboundOrbits          =.false.
    double precision                                             , parameter     :: expArgumentMaximum           =100.0d0
    double precision                                             , parameter     :: A                            =0.216d0, b                 =1.3d0, &  !   Fitting parameters from eqn. (6) of Boylan-Kolchin et al.
         &                                                                          c                            =1.9d0  , d                 =1.0d0
    double precision                                                             :: equivalentCircularOrbitRadius        , massRatio               , &
         &                                                                          orbitalCircularity                   , radialScale             , &
         &                                                                          velocityScale                        , expArgument
    integer                                                                      :: errorCode
    !$GLC attributes unused :: self

    ! Find the host node.
    if (node%isSatellite()) then
       nodeHost => node%parent
    else
       nodeHost => node%parent%firstChild
    end if
    ! Get velocity scale.
    velocityScale=self%darkMatterHaloScale_%velocityVirial(nodeHost)
    radialScale  =self%darkMatterHaloScale_%radiusVirial  (nodeHost)
    ! Get the equivalent circular orbit.
    equivalentCircularOrbitRadius=Satellite_Orbit_Equivalent_Circular_Orbit_Radius(nodeHost,orbit,self%darkMatterHaloScale_,errorCode)
    ! Check error codes.
    select case (errorCode)
    case (errorCodeOrbitUnbound     )
       orbitalCircularity               =1.0d0
       boylanKolchin2008TimeUntilMerging=satelliteMergeTimeInfinite
       return
    case (errorCodeNoEquivalentOrbit)
       ! Circularity is not defined. Assume instantaneous merging.
       orbitalCircularity               =1.0d0
       boylanKolchin2008TimeUntilMerging=0.0d0
       return
    case (errorCodeSuccess          )
       ! Compute orbital circularity.
       massDistribution_  =>  self             %darkMatterProfileDMO_%get            (nodeHost                     )
       orbitalCircularity =  +orbit                                  %angularMomentum(                             ) &
            &                /massDistribution_                      %rotationCurve  (equivalentCircularOrbitRadius) &
            &                /equivalentCircularOrbitRadius
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    case default
       orbitalCircularity=0.0d0
       call Error_Report('unrecognized error code'//{introspection:location})
    end select
    ! Compute mass ratio (mass in host [not including satellite if the node is already a satellite] divided by mass in satellite).
    basic     =>  node     %basic()
    basicHost =>  nodeHost %basic()
    if (node%isSatellite()) then
       massRatio=+basicHost%mass () &
            &    /basic    %mass () &
            &    -1.0d0
    else
       massRatio=+basicHost%mass () &
            &    /basic    %mass ()
    end if
    if (massRatio <= 0.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       boylanKolchin2008TimeUntilMerging=0.0d0
    else
       ! Compute dynamical friction timescale.
       expArgument=min(expArgumentMaximum,c*orbitalCircularity)
       boylanKolchin2008TimeUntilMerging=+self%timescaleMultiplier                               &
            &                            *self%darkMatterHaloScale_%timescaleDynamical(nodeHost) &
            &                            *A                                                      &
            &                            *          massRatio**b                                 &
            &                            /log(1.0d0+massRatio   )                                &
            &                            *exp(expArgument)                                       &
            &                            *(equivalentCircularOrbitRadius/radialScale)**d
    end if
    return
  end function boylanKolchin2008TimeUntilMerging
