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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a black hole binary separation growth class which follows a modified version of :cite:t:`volonteri_assembly_2003`,
  with dynamical friction computed from the Chandrasekhar integrals of the mass distributions of the host node.
  !!}

  use :: Satellite_Dynamical_Friction_Gaseous, only : dynamicalFrictionGaseousClass

  !![
  <blackHoleBinarySeparationGrowthRate name="blackHoleBinarySeparationGrowthRateChandrasekharIntegral" docformat="rst">
   <description>
   A black hole binary separation growth class identical to :galacticus-class:`blackHoleBinarySeparationGrowthRateStandard`
   (including the hardening due to scattering of stars, and the emission of gravitational waves), except that the acceleration
   due to dynamical friction is computed from the Chandrasekhar integrals of the mass distributions of the host node, in the same
   way as for the dynamical friction experienced by satellites. The acceleration is

   .. math::

      \dot{v}_\mathrm{DF} = 4 \pi \mathrm{G}^2 M_\bullet \mathbf{I} \cdot \hat{\mathbf{v}},

   where :math:`\mathbf{I}` is the Chandrasekhar integral of the full mass distribution of the node (the dark matter halo, the
   galactic components, and the circumgalactic medium), evaluated at the position of the black hole, which is assumed to move on a
   circular orbit, co-rotating with the disk, at the circular velocity, :math:`V_\mathrm{c}`, with velocity
   :math:`\mathbf{v}`. Each component contributes according to its own kinematics: collisionless components through the
   :cite:t:`chandrasekhar_dynamical_1943` formula, with their own velocity dispersions, and collisional (gaseous) components, such
   as the interstellar and circumgalactic media, through the class given by ``[dynamicalFrictionGaseous]``. Since the drag due to
   a disk is computed from the velocity of the perturber relative to the rotating disk, a black hole co-rotating with the disk
   experiences little drag from it.

   A single Coulomb logarithm, :math:`\ln\Lambda = \frac{1}{2}\ln(1+\Lambda^2)` with

   .. math::

      \Lambda =  {a \sigma_\mathrm{spheroid}^2  \over \mathrm{G}(M_{\bullet, 1}+M_{\bullet, 2})},

   is used for all components, where :math:`\sigma_\mathrm{spheroid}` is the velocity dispersion of the spheroid, computed as
   described for :galacticus-class:`blackHoleBinarySeparationGrowthRateStandard`.
   </description>
  </blackHoleBinarySeparationGrowthRate>
  !!]
  type, extends(blackHoleBinarySeparationGrowthRateStandard) :: blackHoleBinarySeparationGrowthRateChandrasekharIntegral
     !!{RST
     A black hole binary separation growth class which follows a modified version of :cite:t:`volonteri_assembly_2003`, with
     dynamical friction computed from the Chandrasekhar integrals of the mass distributions of the host node.
     !!}
     private
     class(dynamicalFrictionGaseousClass), pointer :: dynamicalFrictionGaseous_ => null()
   contains
     final     ::                                  chandrasekharIntegralDestructor
     procedure :: accelerationDynamicalFriction => chandrasekharIntegralAccelerationDynamicalFriction
  end type blackHoleBinarySeparationGrowthRateChandrasekharIntegral

  interface blackHoleBinarySeparationGrowthRateChandrasekharIntegral
     !!{RST
     Constructors for the :galacticus-class:`blackHoleBinarySeparationGrowthRateChandrasekharIntegral` black hole binary
     separation growth rate class.
     !!}
     module procedure chandrasekharIntegralConstructorParameters
     module procedure chandrasekharIntegralConstructorInternal
  end interface blackHoleBinarySeparationGrowthRateChandrasekharIntegral

contains

  function chandrasekharIntegralConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`blackHoleBinarySeparationGrowthRateChandrasekharIntegral` black hole binary separation
    growth rate class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (blackHoleBinarySeparationGrowthRateChandrasekharIntegral)                :: self
    type   (inputParameters                                         ), intent(inout) :: parameters
    class  (darkMatterHaloScaleClass                                ), pointer       :: darkMatterHaloScale_
    class  (dynamicalFrictionGaseousClass                           ), pointer       :: dynamicalFrictionGaseous_
    logical                                                                          :: stellarDensityChangeBinaryMotion, computeVelocityDispersion

    !![
    <inputParameter docformat="rst">
      <name>stellarDensityChangeBinaryMotion</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, account for the change in stellar density caused by the black hole binary's motion through the stellar background
      when computing the hardening rate of the binary via stellar scattering.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>computeVelocityDispersion</name>
      <defaultValue>.false.</defaultValue>
      <description>
      Specifies whether or not the velocity dispersion of the spheroid should be computed using Jeans equation in black hole binary
      hardening calculations, and in the Coulomb logarithm for dynamical friction. If ``false``, then the velocity dispersion is
      assumed to equal the characteristic velocity of the spheroid.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    <objectBuilder class="dynamicalFrictionGaseous" name="dynamicalFrictionGaseous_" source="parameters"/>
    !!]
    self=blackHoleBinarySeparationGrowthRateChandrasekharIntegral(stellarDensityChangeBinaryMotion,computeVelocityDispersion,darkMatterHaloScale_,dynamicalFrictionGaseous_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"     />
    <objectDestructor name="dynamicalFrictionGaseous_"/>
    !!]
    return
  end function chandrasekharIntegralConstructorParameters

  function chandrasekharIntegralConstructorInternal(stellarDensityChangeBinaryMotion,computeVelocityDispersion,darkMatterHaloScale_,dynamicalFrictionGaseous_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`blackHoleBinarySeparationGrowthRateChandrasekharIntegral` black hole binary
    separation growth rate class.
    !!}
    implicit none
    type   (blackHoleBinarySeparationGrowthRateChandrasekharIntegral)                        :: self
    class  (darkMatterHaloScaleClass                                ), intent(in   ), target :: darkMatterHaloScale_
    class  (dynamicalFrictionGaseousClass                           ), intent(in   ), target :: dynamicalFrictionGaseous_
    logical                                                          , intent(in   )         :: stellarDensityChangeBinaryMotion, computeVelocityDispersion
    !![
    <constructorAssign variables="stellarDensityChangeBinaryMotion,computeVelocityDispersion,*darkMatterHaloScale_,*dynamicalFrictionGaseous_"/>
    !!]

    return
  end function chandrasekharIntegralConstructorInternal

  subroutine chandrasekharIntegralDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`blackHoleBinarySeparationGrowthRateChandrasekharIntegral` black hole binary separation
    growth rate class.
    !!}
    implicit none
    type(blackHoleBinarySeparationGrowthRateChandrasekharIntegral), intent(inout) :: self

    !![
    <objectDestructor name="self%dynamicalFrictionGaseous_"/>
    !!]
    return
  end subroutine chandrasekharIntegralDestructor

  double precision function chandrasekharIntegralAccelerationDynamicalFriction(self,node,blackHole,blackHoleCentral,coordinates,velocityRotation,velocityDispersionSpheroid,velocityDispersionDarkMatter) result(accelerationDynamicalFriction)
    !!{RST
    Return the acceleration (in units of km/s/Gyr, and negative for a decelerating black hole) due to dynamical friction of a black
    hole in a binary, computed from the Chandrasekhar integrals of the mass distributions of the host node.
    !!}
    use :: Coordinates                     , only : coordinate               , coordinateCartesian           , assignment(=)
    use :: Galactic_Structure_Options      , only : componentTypeBlackHole   , massTypeBlackHole
    use :: Galacticus_Nodes                , only : treeNode
    use :: Mass_Distributions              , only : massDistributionBlackHole, massDistributionClass
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr        , gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (blackHoleBinarySeparationGrowthRateChandrasekharIntegral), intent(inout) :: self
    type            (treeNode                                                ), intent(inout) :: node
    class           (nodeComponentBlackHole                                  ), intent(inout) :: blackHole                   , blackHoleCentral
    class           (coordinate                                              ), intent(in   ) :: coordinates
    double precision                                                          , intent(in   ) :: velocityRotation            , velocityDispersionSpheroid, &
         &                                                                                       velocityDispersionDarkMatter
    class           (massDistributionClass                                   ), pointer       :: massDistribution_
    type            (massDistributionBlackHole                               ), pointer       :: massDistributionPerturber_
    double precision                                                          , dimension(3)  :: integral_                   , velocity_
    type            (coordinateCartesian                                     )                :: integral                    , velocity
    double precision                                                                          :: argumentCoulomb             , logarithmCoulomb
    !$GLC attributes unused :: velocityDispersionDarkMatter

    accelerationDynamicalFriction=0.0d0
    if (velocityRotation <= 0.0d0) return
    ! Compute the Coulomb logarithm, as in the standard class but using the velocity dispersion of the spheroid for all components.
    argumentCoulomb =+(                                    &
         &              +blackHole       %radialPosition() &
         &              *velocityDispersionSpheroid**2     &
         &             )                                   &
         &            /(                                   &
         &              +gravitationalConstant_internal    &
         &              *(                                 &
         &                +blackHole       %mass        () &
         &                +blackHoleCentral%mass        () &
         &               )                                 &
         &             )
    logarithmCoulomb=+0.5d0                         &
         &           *log(1.0d0+argumentCoulomb**2)
    ! The black hole moves on a circular orbit at the circular velocity, co-rotating with the disk. At the position (r,0,0) the
    ! disk rotates in the -y direction.
    velocity_=[0.0d0,-velocityRotation,0.0d0]
    velocity =        velocity_
    ! The perturber is the (point mass) black hole.
    allocate(massDistributionPerturber_)
    !![
    <referenceConstruct object="massDistributionPerturber_" constructor="massDistributionBlackHole(mass=blackHole%mass(),componentType=componentTypeBlackHole,massType=massTypeBlackHole)"/>
    !!]
    ! Evaluate the Chandrasekhar integral of the full mass distribution of the node.
    massDistribution_ => node%massDistribution()
    integral          =  massDistribution_%chandrasekharIntegral(massDistribution_,massDistributionPerturber_,blackHole%mass(),coordinates,velocity,logarithmCoulomb,self%dynamicalFrictionGaseous_)
    integral_         =  integral
    ! Find the acceleration along the direction of motion.
    accelerationDynamicalFriction=+4.0d0                                   &
         &                        *Pi                                      &
         &                        *gravitationalConstant_internal**2       &
         &                        *blackHole%mass()                        &
         &                        *dot_product(integral_,velocity_       ) &
         &                        /                      velocityRotation  &
         &                        /MpcPerKmPerSToGyr
    !![
    <objectDestructor name="massDistribution_"         />
    <objectDestructor name="massDistributionPerturber_"/>
    !!]
    return
  end function chandrasekharIntegralAccelerationDynamicalFriction
