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

  !!{RST
  Implements a class for the timescale of ram pressure stripping of hot halos in which the timescale is estimated from the ram pressure acceleration.
  !!}

  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleClass
  use :: Hot_Halo_Ram_Pressure_Forces, only : hotHaloRamPressureForceClass

  !![
  <hotHaloRamPressureTimescale name="hotHaloRamPressureTimescaleRamPressureAcceleration" docformat="rst">
   <description>
   A hot halo ram pressure timescale class which computes the ram pressure stripping timescale from the acceleration imparted by the ram pressure force. Following :cite:t:`roediger_ram_2007` this is approximated as:

   .. math::

      a_\mathrm{ram pressure} = P_\mathrm{ram pressure}/\Sigma,

   where :math:`P_\mathrm{ram pressure}` is the ram pressure force per unit area, and :math:`\Sigma` is the surface density of gas. The associated timescale to accelerate gas over a distance :math:`r_\mathrm{outer}` (the current outer radius of the hot halo) is then:

   .. math::

      \tau_\mathrm{ram pressure} = \sqrt{2 r_\mathrm{outer} \Sigma_\mathrm{outer} / P_\mathrm{ram pressure}}.
   </description>
  </hotHaloRamPressureTimescale>
  !!]
  type, extends(hotHaloRamPressureTimescaleClass) :: hotHaloRamPressureTimescaleRamPressureAcceleration
     !!{RST
     Implementation of a hot halo ram pressure timescale class in which the timescale is estimated from the ram pressure acceleration.
     !!}
     private
     class(darkMatterHaloScaleClass    ), pointer :: darkMatterHaloScale_     => null()
     class(hotHaloRamPressureForceClass), pointer :: hotHaloRamPressureForce_ => null()
   contains
     final     ::              ramPressureAccelerationDestructor
     procedure :: timescale => ramPressureAccelerationTimescale
  end type hotHaloRamPressureTimescaleRamPressureAcceleration

  interface hotHaloRamPressureTimescaleRamPressureAcceleration
     !!{RST
     Constructors for the :galacticus-class:`hotHaloRamPressureTimescaleRamPressureAcceleration` hot halo ram pressure timescale class.
     !!}
     module procedure ramPressureAccelerationConstructorParameters
     module procedure ramPressureAccelerationConstructorInternal
  end interface hotHaloRamPressureTimescaleRamPressureAcceleration

contains

  function ramPressureAccelerationConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`hotHaloRamPressureTimescaleRamPressureAcceleration` hot halo ram pressure timescale class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloRamPressureTimescaleRamPressureAcceleration)                :: self
    type (inputParameters                                   ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass                          ), pointer       :: darkMatterHaloScale_
    class(hotHaloRamPressureForceClass                      ), pointer       :: hotHaloRamPressureForce_

    !![
    <objectBuilder class="darkMatterHaloScale"     name="darkMatterHaloScale_"     source="parameters"/>
    <objectBuilder class="hotHaloRamPressureForce" name="hotHaloRamPressureForce_" source="parameters"/>
    !!]
    self=hotHaloRamPressureTimescaleRamPressureAcceleration(darkMatterHaloScale_,hotHaloRamPressureForce_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"    />
    <objectDestructor name="hotHaloRamPressureForce_"/>
    !!]
    return
  end function ramPressureAccelerationConstructorParameters

  function ramPressureAccelerationConstructorInternal(darkMatterHaloScale_,hotHaloRamPressureForce_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`hotHaloRamPressureTimescaleRamPressureAcceleration` hot halo ram pressure timescale class.
    !!}
    implicit none
    type (hotHaloRamPressureTimescaleRamPressureAcceleration)                        :: self
    class(darkMatterHaloScaleClass                          ), intent(in   ), target :: darkMatterHaloScale_
    class(hotHaloRamPressureForceClass                      ), intent(in   ), target :: hotHaloRamPressureForce_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *hotHaloRamPressureForce_"/>
    !!]

    return
  end function ramPressureAccelerationConstructorInternal

  subroutine ramPressureAccelerationDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`hotHaloRamPressureTimescaleRamPressureAcceleration` hot halo ram pressure timescale class.
    !!}
    implicit none
    type(hotHaloRamPressureTimescaleRamPressureAcceleration), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"    />
    <objectDestructor name="self%hotHaloRamPressureForce_"/>
    !!]
    return
  end subroutine ramPressureAccelerationDestructor

  double precision function ramPressureAccelerationTimescale(self,node)
    !!{RST
    Computes the hot halo ram pressure stripping timescale, based on the acceleration due to ram pressure forces. This timescale is approximated as :cite:p:`roediger_ram_2007` :math:`\tau \approx \sqrt{2 r_\mathrm{outer} \Sigma_\mathrm{outer} / P_\mathrm{ ram}}`, where :math:`r_\mathrm{outer}` is the current outer radius of the hot halo, :math:`\Sigma_\mathrm{outer}` is the surface density at that radius, and :math:`P_\mathrm{ram}` is the ram pressure force (per unit area). The surface density is approximated as :math:`\Sigma_\mathrm{outer} \approx r_\mathrm{outer} \rho_\mathrm{outer}`, where :math:`\rho_\mathrm{outer}` is the density at the outer radius.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentHotHalo , treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear             , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Coordinates                     , only : coordinateSpherical  , assignment(=)
    use :: Galactic_Structure_Options      , only : componentTypeHotHalo , massTypeGaseous
    implicit none
    class           (hotHaloRamPressureTimescaleRamPressureAcceleration), intent(inout) :: self
    type            (treeNode                                          ), intent(inout) :: node
    class           (nodeComponentHotHalo                              ), pointer       :: hotHalo
    type            (treeNode                                          ), pointer       :: nodeHost
    double precision                                                    , parameter     :: timescaleInfinite       =huge(1.0d0)
    double precision                                                    , parameter     :: velocityStrippingMaximum=     1.0d1
    class           (massDistributionClass                             ), pointer       :: massDistribution_
    type            (coordinateSpherical                               )                :: coordinates
    double precision                                                                    :: radiusOuter                         , densityOuter       , &
         &                                                                                 forceRamPressure                    , surfaceDensityOuter

    ! Evaluate surface density and ram pressure force.
    massDistribution_   =>  node                                      %massDistribution(componentTypeHotHalo,massTypeGaseous)
    hotHalo             =>  node                                      %hotHalo         (                                    )
    radiusOuter         =   hotHalo                                   %outerRadius     (                                    )
    coordinates         =  [radiusOuter,0.0d0,0.0d0]
    densityOuter        =   massDistribution_                         %density         (coordinates                         )
    forceRamPressure    =   self             %hotHaloRamPressureForce_%force           (node                                )
    surfaceDensityOuter =  +radiusOuter  &
         &                 *densityOuter
    !![
    <objectDestructor name="massDistribution_"/>
    !!]          
    ! Exit with infinite timescale for zero ram pressure force.
    if (forceRamPressure <= 0.0d0) then
       ramPressureAccelerationTimescale=timescaleInfinite
       return
    end if
    ! For zero density or radius, return the halo dynamical time (the timescale should be irrelevant in such cases anyway).
    if (radiusOuter <= 0.0d0 .or. surfaceDensityOuter <= 0.0d0) then
       ramPressureAccelerationTimescale=self%darkMatterHaloScale_%timescaleDynamical(node)
       return
    end if
    ! Find the hosting node.
    nodeHost => node%parent
    ! Evaluate the timescale.
    ramPressureAccelerationTimescale=+(                                                        &
         &                             +megaParsec                                             &
         &                             /kilo                                                   &
         &                             /gigaYear                                               &
         &                            )                                                        &
         &                           *max(                                                     &
         &                                +sqrt(                                               &
         &                                      +2.0d0                                         &
         &                                      *radiusOuter                                   &
         &                                      *surfaceDensityOuter                           &
         &                                      /forceRamPressure                              &
         &                                     )                                             , &
         &                                +radiusOuter                                         &
         &                                /velocityStrippingMaximum                            &
         &                                /self%darkMatterHaloScale_%velocityVirial(nodeHost)  &
         &                               )
    return
  end function ramPressureAccelerationTimescale
