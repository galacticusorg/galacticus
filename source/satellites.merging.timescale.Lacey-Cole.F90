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
  Implements a satellite merging timescale class which uses the \cite{lacey_merger_1993} method.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesLaceyCole1993">
   <description>
    A satellite merging timescale class which computes merging timescales using the dynamical friction calculation of
    \cite{lacey_merger_1993}. Timescales are multiplied by the value of the {\normalfont \ttfamily [timescaleMultiplier]} input
    parameter.
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesLaceyCole1993
     !!{
     A class implementing the \cite{lacey_merger_1993} method for satellite merging timescales.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: timescaleMultiplier
   contains
     !![
     <methods>
       <method description="Return the mass-dependent part of the time (in Gyr) until the satellite will merge with its host." method="timeUntilMergingMassDependence" />
     </methods>
     !!]
     final     ::                                   laceyCole1993Destructor
     procedure :: timeUntilMerging               => laceyCole1993TimeUntilMerging
     procedure :: timeUntilMergingMassDependence => laceyCole1993TimeUntilMergingMassDependence
  end type satelliteMergingTimescalesLaceyCole1993

  interface satelliteMergingTimescalesLaceyCole1993
     !!{
     Constructors for the \cite{lacey_merger_1993} merging timescale class.
     !!}
     module procedure laceyCole1993ConstructorParameters
     module procedure laceyCole1993ConstructorInternal
  end interface satelliteMergingTimescalesLaceyCole1993

contains

  function laceyCole1993ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \cite{lacey_merger_1993} merging timescale class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteMergingTimescalesLaceyCole1993)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    double precision                                                         :: timescaleMultiplier

    !![
    <inputParameter>
      <name>timescaleMultiplier</name>
      <defaultValue>0.75d0</defaultValue>
      <description>A multiplier for the merging timescale in dynamical friction timescale calculations.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=satelliteMergingTimescalesLaceyCole1993(timescaleMultiplier,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function laceyCole1993ConstructorParameters

  function laceyCole1993ConstructorInternal(timescaleMultiplier,darkMatterHaloScale_) result(self)
    !!{
    Constructor for the \cite{lacey_merger_1993} merging timescale class.
    !!}
    implicit none
    type            (satelliteMergingTimescalesLaceyCole1993)                        :: self
    double precision                                         , intent(in   )         :: timescaleMultiplier
    class           (darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="timescaleMultiplier, *darkMatterHaloScale_"/>
    !!]

    return
  end function laceyCole1993ConstructorInternal

  subroutine laceyCole1993Destructor(self)
    !!{
    Destructor for the \refClass{satelliteMergingTimescalesLaceyCole1993} satellite merging timescale class.
    !!}
    implicit none
    type(satelliteMergingTimescalesLaceyCole1993), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine laceyCole1993Destructor

  double precision function laceyCole1993TimeUntilMerging(self,node,orbit)
    !!{
    Return the timescale for merging satellites using the \cite{lacey_merger_1993} method.
    !!}
    use :: Kepler_Orbits, only : keplerOrbit
    implicit none
    class           (satelliteMergingTimescalesLaceyCole1993), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    type            (keplerOrbit                            ), intent(inout) :: orbit
    type            (treeNode                               ), pointer       :: nodeHost
    double precision                                                         :: equivalentCircularOrbitRadius, orbitalCircularity, &
         &                                                                      radialScale                  , velocityScale
    
    ! Find the host node.
    if (node%isSatellite()) then
       nodeHost => node%parent
    else
       nodeHost => node%parent%firstChild
    end if
    ! Get velocity scale.
    velocityScale=self%darkMatterHaloScale_%velocityVirial(nodeHost)
    radialScale  =self%darkMatterHaloScale_%radiusVirial  (nodeHost)
    ! Compute radius of orbit with same energy.
    equivalentCircularOrbitRadius=exp(orbit%energy()/velocityScale**2+0.5d0)
    ! Compute orbital circularity.
    orbitalCircularity=orbit%angularMomentum()/velocityScale/radialScale/equivalentCircularOrbitRadius
    ! Compute the merging timescale.
    laceyCole1993TimeUntilMerging                         &
         & =orbitalCircularity           **0.78d0         &
         & *equivalentCircularOrbitRadius**2              &
         & *self%timeUntilMergingMassDependence(node)
    return
  end function laceyCole1993TimeUntilMerging

  double precision function laceyCole1993TimeUntilMergingMassDependence(self,node)
    !!{
    Return the mass-dependent part of the timescale for merging satellites using the \cite{lacey_merger_1993} method.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (satelliteMergingTimescalesLaceyCole1993), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    type            (treeNode                               ), pointer       :: nodeHost
    class           (nodeComponentBasic                     ), pointer       :: basicHost                 , basic
    double precision                                         , parameter     :: inverseTwoB1=1.169335453d0            !  1/2/B(1).
    double precision                                                         :: massRatio
    !$GLC attributes unused :: self

    ! Find the host node.
    nodeHost => node%parent
    ! Compute mass ratio.
    basic     => node    %basic()
    basicHost => nodeHost%basic()
    if (node%isSatellite()) then
       ! Node is already a satellite in its host - compute the mass ratio directly.
       massRatio=+basicHost%mass () &
            &    /basic    %mass ()
    else
       ! Node is not yet a satellite in its host - correct the host mass to what it will be after the node becomes a satellite in the
       ! host.
       massRatio=+basicHost%mass () &
            &    /basic    %mass () &
            &    +1.0d0
    end if
    ! Check for a greater than unity mass ratio.
    if (massRatio <= 1.0d0) then
       ! Assume zero merging time as the satellite is as massive as the host.
       laceyCole1993TimeUntilMergingMassDependence=0.0d0
    else
       ! Compute dynamical friction timescale.
       laceyCole1993TimeUntilMergingMassDependence=+self%timescaleMultiplier                               &
            &                                      *self%darkMatterHaloScale_%timescaleDynamical(nodeHost) &
            &                                      *inverseTwoB1                                           &
            &                                      *    massRatio                                          &
            &                                      /log(massRatio)
    end if
    return
  end function laceyCole1993TimeUntilMergingMassDependence
