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
  An implementation of virial orbits which assumes an isotropic distribution of infall directions and tangential velocities
  applied to another virial orbit class.
  !!}

  !![
  <virialOrbit name="virialOrbitIsotropic">
   <description>Virial orbits which assumes an isotropic distribution of infall directions and tangential velocities applied to another virial orbit class.</description>
  </virialOrbit>
  !!]
  type, extends(virialOrbitClass) :: virialOrbitIsotropic
     !!{
     A virial orbit class which assumes an isotropic distribution of infall directions and tangential velocities applied to
     another virial orbit class
     !!}
     private
     class(virialOrbitClass), pointer :: virialOrbit_ => null()
   contains
     final     ::                                    isotropicDestructor
     procedure :: orbit                           => isotropicOrbit
     procedure :: densityContrastDefinition       => isotropicDensityContrastDefinition
     procedure :: velocityTangentialMagnitudeMean => isotropicVelocityTangentialMagnitudeMean
     procedure :: velocityTangentialVectorMean    => isotropicVelocityTangentialVectorMean
     procedure :: angularMomentumMagnitudeMean    => isotropicAngularMomentumMagnitudeMean
     procedure :: angularMomentumVectorMean       => isotropicAngularMomentumVectorMean
     procedure :: velocityTotalRootMeanSquared    => isotropicVelocityTotalRootMeanSquared
     procedure :: energyMean                      => isotropicEnergyMean
     procedure :: isAngularlyResolved             => isotropicIsAngularlyResolved
  end type virialOrbitIsotropic

  interface virialOrbitIsotropic
     !!{
     Constructors for the \refClass{virialOrbitIsotropic} virial orbit class.
     !!}
     module procedure isotropicConstructorParameters
     module procedure isotropicConstructorInternal
  end interface virialOrbitIsotropic

contains

  function isotropicConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{virialOrbitIsotropic} satellite virial orbit class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (virialOrbitIsotropic)                :: self
    type (inputParameters     ), intent(inout) :: parameters
    class(virialOrbitClass    ), pointer       :: virialOrbit_

    !![
    <objectBuilder class="virialOrbit"  name="virialOrbit_" source="parameters"/>
    !!]
    self=virialOrbitIsotropic(virialOrbit_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialOrbit_"/>
    !!]
    return
  end function isotropicConstructorParameters

  function isotropicConstructorInternal(virialOrbit_) result(self)
    !!{
    Internal constructor for the \refClass{virialOrbitIsotropic} virial orbits class.
    !!}
    implicit none
    type (virialOrbitIsotropic)                        :: self
    class(virialOrbitClass    ), intent(in   ), target :: virialOrbit_
    !![
    <constructorAssign variables="*virialOrbit_"/>
    !!]

    return
  end function isotropicConstructorInternal

  subroutine isotropicDestructor(self)
    !!{
    Destructor for the \refClass{virialOrbitIsotropic} virial orbits class.
    !!}
    implicit none
    type(virialOrbitIsotropic), intent(inout) :: self

    !![
    <objectDestructor name="self%virialOrbit_"/>
    !!]
    return
  end subroutine isotropicDestructor

  function isotropicOrbit(self,node,host,acceptUnboundOrbits)
    !!{
    Return isotropic orbital parameters for a satellite.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type   (keplerOrbit         )                        :: isotropicOrbit
    class  (virialOrbitIsotropic), intent(inout), target :: self
    type   (treeNode            ), intent(inout)         :: host               , node
    logical                      , intent(in   )         :: acceptUnboundOrbits

    ! Get the underlying orbit.
    isotropicOrbit=self%virialOrbit_%orbit(node,host,acceptUnboundOrbits)
    ! Sample from an isotropic distribution.
    call isotropicOrbit%phiSet    (     2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()       )
    call isotropicOrbit%thetaSet  (acos(2.0d0   *node%hostTree%randomNumberGenerator_%uniformSample()-1.0d0))
    call isotropicOrbit%epsilonSet(     2.0d0*Pi*node%hostTree%randomNumberGenerator_%uniformSample()       )
    return
  end function isotropicOrbit

  function isotropicDensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of virial orbits.
    !!}
    implicit none
    class(virialDensityContrastClass), pointer       :: isotropicDensityContrastDefinition
    class(virialOrbitIsotropic      ), intent(inout) :: self

    isotropicDensityContrastDefinition => self%virialOrbit_%densityContrastDefinition()
    return
  end function isotropicDensityContrastDefinition

  double precision function isotropicVelocityTangentialMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the tangential velocity.
    !!}
    implicit none
    class(virialOrbitIsotropic), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node, host

    isotropicVelocityTangentialMagnitudeMean=self%virialOrbit_%velocityTangentialMagnitudeMean(node,host)
    return
  end function isotropicVelocityTangentialMagnitudeMean

  function isotropicVelocityTangentialVectorMean(self,node,host)
    !!{
    Return the mean of the vector tangential velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                      , dimension(3)  :: isotropicVelocityTangentialVectorMean
    class           (virialOrbitIsotropic), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node                                 , host
    !$GLC attributes unused :: self, node, host

    ! Since the tangential velocity is assumed to be isotropically distributed the mean of the vector tangential velocity is zero.
    isotropicVelocityTangentialVectorMean=0.0d0
    return
  end function isotropicVelocityTangentialVectorMean

  double precision function isotropicAngularMomentumMagnitudeMean(self,node,host)
    !!{
    Return the mean magnitude of the angular momentum.
    !!}
    implicit none
    class(virialOrbitIsotropic), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node, host

    isotropicAngularMomentumMagnitudeMean=self%virialOrbit_%angularMomentumMagnitudeMean(node,host)
    return
  end function isotropicAngularMomentumMagnitudeMean

  function isotropicAngularMomentumVectorMean(self,node,host)
    !!{
    Return the mean of the vector tangential velocity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision                      , dimension(3)  :: isotropicAngularMomentumVectorMean
    class           (virialOrbitIsotropic), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node                              , host
    !$GLC attributes unused :: self, node, host

    ! Since the tangential velocity is assumed to be isotropically distributed the mean of the vector angular momentum is zero.
    isotropicAngularMomentumVectorMean=0.0d0
    return
  end function isotropicAngularMomentumVectorMean

  double precision function isotropicVelocityTotalRootMeanSquared(self,node,host)
    !!{
    Return the root mean squared of the total velocity.
    !!}
    implicit none
    class(virialOrbitIsotropic), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node, host

    isotropicVelocityTotalRootMeanSquared=self%virialOrbit_%velocityTotalRootMeanSquared(node,host)
    return
  end function isotropicVelocityTotalRootMeanSquared

  double precision function isotropicEnergyMean(self,node,host)
    !!{
    Return the mean of the total energy.
    !!}
    implicit none
    class(virialOrbitIsotropic), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node, host

    isotropicEnergyMean=self%virialOrbit_%energyMean(node,host)
    return
  end function isotropicEnergyMean

  logical function isotropicIsAngularlyResolved(self)
    !!{
    Return true indicating that orbits are angularly-resolved.
    !!}
    implicit none
    class(virialOrbitIsotropic), intent(inout) :: self
    !$GLC attributes unused :: self

    isotropicIsAngularlyResolved=.true.
    return
  end function isotropicIsAngularlyResolved
  
