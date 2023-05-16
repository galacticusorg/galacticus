!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

  use :: Galactic_Structure, only : galacticStructureClass

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
     class(galacticStructureClass), pointer :: galacticStructure_=> null()
   contains
     final     ::          orbitalPositionDestructor
     procedure :: force => orbitalPositionForce
  end type hotHaloRamPressureForceOrbitalPosition

  interface hotHaloRamPressureForceOrbitalPosition
     !!{
     Constructors for the {\normalfont \ttfamily orbitalPosition} hot halo ram pressure force class.
     !!}
     module procedure orbitalPositionConstructorParameters
     module procedure orbitalPositionConstructorInternal
  end interface hotHaloRamPressureForceOrbitalPosition

contains

  function orbitalPositionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily orbitalPosition} hot halo ram pressure force class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloRamPressureForceOrbitalPosition)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(galacticStructureClass                ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=hotHaloRamPressureForceOrbitalPosition(galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function orbitalPositionConstructorParameters

  function orbitalPositionConstructorInternal(galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily orbitalPosition} hot halo ram pressure force class.
    !!}
    implicit none
    type (hotHaloRamPressureForceOrbitalPosition)                        :: self
    class(galacticStructureClass                ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function orbitalPositionConstructorInternal

  subroutine orbitalPositionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily orbitalPosition} hot halo ram pressure force class.
    !!}
    implicit none
    type(hotHaloRamPressureForceOrbitalPosition), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine orbitalPositionDestructor

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
    ! Find the ram pressure force at pericenter.
    coordinates          =  [radiusOrbital,0.0d0,0.0d0]
    massDistribution_    =>  nodeHost         %massDistribution(                                                                       )
    orbitalPositionForce =  +massDistribution_%density         (coordinates,componentType=componentTypeHotHalo,massType=massTypeGaseous)    &
         &                  *velocityOrbital                                                                                            **2
    return
  end function orbitalPositionForce
