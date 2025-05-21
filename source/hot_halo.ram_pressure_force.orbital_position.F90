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
  Implements a model of the ram pressure stripping force from hot halos based on orbital position within the host halo.
  !!}

  !![
  <hotHaloRamPressureForce name="hotHaloRamPressureForceOrbitalPosition">
   <description>
    A hot halo ram pressure force class which computes the force based on the current orbital position within the host
    halo. Specifically, the ram pressure force is    
    \begin{equation}
    \mathcal{F}_\mathrm{ram, hot, host} = \rho_\mathrm{hot, host}(r) v^2(r),
    \end{equation}
    where $\rho_\mathrm{hot, host}(r)$ is the hot halo density profile of the node's host halo, $v(r)$ is the orbital velocity
    of the node in that host, and $r$ is the instantaneous radius of the node's orbit.
   </description>
  </hotHaloRamPressureForce>
  !!]
  type, extends(hotHaloRamPressureForceClass) :: hotHaloRamPressureForceOrbitalPosition
     !!{
     Implementation of a hot halo ram pressure force class based on orbital position within the host halo.
     !!}
     private
   contains
     procedure :: force => orbitalPositionForce
  end type hotHaloRamPressureForceOrbitalPosition

  interface hotHaloRamPressureForceOrbitalPosition
     !!{
     Constructors for the \refClass{hotHaloRamPressureForceOrbitalPosition} hot halo ram pressure force class.
     !!}
     module procedure orbitalPositionConstructorParameters
     module procedure orbitalPositionConstructorInternal
  end interface hotHaloRamPressureForceOrbitalPosition

contains

  function orbitalPositionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hotHaloRamPressureForceOrbitalPosition} hot halo ram pressure force class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(hotHaloRamPressureForceOrbitalPosition)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=hotHaloRamPressureForceOrbitalPosition()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function orbitalPositionConstructorParameters

  function orbitalPositionConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{hotHaloRamPressureForceOrbitalPosition} hot halo ram pressure force class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Error_Report             , Component_List
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none
    type (hotHaloRamPressureForceOrbitalPosition) :: self

    ! Ensure that required methods are supported.
    if     (                                                                                                                         &
         &  .not.                                                                                                                    &
         &       (                                                                                                                   &
         &        defaultSatelliteComponent%positionIsGettable().and.                                                                &
         &        defaultSatelliteComponent%velocityIsGettable()                                                                     &
         &  )                                                                                                                        &
         & ) call Error_Report                                                                                                       &
         &        (                                                                                                                  &
         &         'this method requires that position, and velocity properties must all be gettable for the satellite component.'// &
         &         Component_List(                                                                                                   &
         &                        'satellite'                                                                                     ,  &
         &                        defaultSatelliteComponent%positionAttributeMatch(requireGettable=.true.).intersection.             &
         &                        defaultSatelliteComponent%velocityAttributeMatch(requireGettable=.true.)                           &
         &                       )                                                                                                // &
         &         {introspection:location}                                                                                          &
         &        )
    
    return
  end function orbitalPositionConstructorInternal

  double precision function orbitalPositionForce(self,node)
    !!{
    Return a ram pressure force due to the hot halo based on orbital position within the host halo.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentSatellite
    use :: Vectors                   , only : Vector_Magnitude
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Coordinates               , only : coordinateSpherical   , assignment(=)
    use :: Galactic_Structure_Options, only : componentTypeHotHalo  , massTypeGaseous
    implicit none
    class           (hotHaloRamPressureForceOrbitalPosition), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    class           (nodeComponentSatellite                ), pointer       :: satellite
    type            (treeNode                              ), pointer       :: nodeHost
    class           (massDistributionClass                 ), pointer       :: massDistribution_
    type            (coordinateSpherical                   )                :: coordinates
    double precision                                                        :: radiusOrbital, velocityOrbital

    ! Find the host node.
    nodeHost             =>  node%parent
    ! Get the satellite component.
    satellite            =>  node%satellite()
    ! Compute orbital position and velocity.
    radiusOrbital        =  +Vector_Magnitude(satellite%position())
    velocityOrbital      =  +Vector_Magnitude(satellite%velocity())
    ! Find the ram pressure force this orbital radius.
    coordinates          =  [radiusOrbital,0.0d0,0.0d0]
    massDistribution_    =>  nodeHost         %massDistribution(componentTypeHotHalo,massTypeGaseous)
    orbitalPositionForce =  +massDistribution_%density         (coordinates                         )    &
         &                  *velocityOrbital                                                         **2
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function orbitalPositionForce
