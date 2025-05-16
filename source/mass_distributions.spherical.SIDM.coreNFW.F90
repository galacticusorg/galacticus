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
  An implementation of a cored-NFW mass distribution to approximate the effects of SIDM based on the model of \cite{jiang_semi-analytic_2023}.
  !!}

  !![
  <massDistribution name="massDistributionSphericalSIDMCoreNFW">
   <description>
     A mass distribution class implementing a cored-NFW dark matter halo profile to approximate the effects of SIDM based
     on the model of Jiang et al. (2022). The profile is defined by the enclosed mass, with \citep{jiang_semi-analytic_2023}:
     \begin{equation}
       M(r) = M_\mathrm{NFW}(r) \mathrm{tanh}\left(\frac{r}{r_\mathrm{c}}\right),
     \end{equation}
     where $r_\mathrm{c} = \alpha r_1$ is a characteristic core size related to the interaction radius $r_1$ by a constant factor
     $\alpha = ${\normalfont \ttfamily [factorRadiusCore]}.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalSIDM) :: massDistributionSphericalSIDMCoreNFW
     !!{
     Implementation of a cored-NFW mass distribution to approximate the effects of SIDM based on the model of \cite{jiang_semi-analytic_2023}.
     !!}
     private
     double precision :: factorRadiusCore
   contains
     !![
     <methods>
       <method method="radiusCore" description="Computes the core radius of halo."/>
     </methods>
     !!]
     final     ::                           sphericalSIDMCoreNFWDestructor
     procedure :: radiusCore             => sphericalSIDMCoreNFWRadiusCore
     procedure :: density                => sphericalSIDMCoreNFWDensity
     procedure :: densityGradientRadial  => sphericalSIDMCoreNFWDensityGradientRadial
     procedure :: massEnclosedBySphere   => sphericalSIDMCoreNFWMassEnclosedBySphere
  end type massDistributionSphericalSIDMCoreNFW

  interface massDistributionSphericalSIDMCoreNFW
     !!{
     Constructors for the \refClass{massDistributionSphericalSIDMCoreNFW} mass distribution class.
     !!}
     module procedure sphericalSIDMCoreNFWConstructorParameters
     module procedure sphericalSIDMCoreNFWConstructorInternal
  end interface massDistributionSphericalSIDMCoreNFW

contains

  function sphericalSIDMCoreNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalSIDMCoreNFW} mass distribution class which takes a parameter set as input.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalSIDMCoreNFW)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (massDistributionClass               ), pointer       :: massDistribution_
    class           (darkMatterParticleClass             ), pointer       :: darkMatterParticle_
    double precision                                                      :: factorRadiusCore   , timeAge
    type            (varying_string                      )                :: componentType      , massType, &
         &                                                                   nonAnalyticSolver

    !![
    <inputParameter>
      <name>timeAge</name>
      <source>parameters</source>
      <description>The age of the halo (in Gyr).</description>
    </inputParameter>
    <inputParameter>
      <name>factorRadiusCore</name>
      <defaultValue>0.45d0</defaultValue>
      <defaultSource>\cite{jiang_semi-analytic_2023}</defaultSource>
      <source>parameters</source>
      <description>The factor $\alpha$ appearing in the definition of the core radius, $r_\mathrm{c}=\alpha r_1$ where $r_1$ is the radius at which an SIDM particle has had, on average, 1 interaction.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="massDistribution"   name="massDistribution_"   source="parameters"/>
    <objectBuilder class="darkMatterParticle" name="darkMatterParticle_" source="parameters"/>
    !!]
    select type (massDistribution_)
    class is (massDistributionNFW)
       self=massDistributionSphericalSIDMCoreNFW(factorRadiusCore,timeAge,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),massDistribution_,darkMatterParticle_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('an NFW mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"/>
    <objectDestructor name="massDistribution_"  />
    !!]
    return
  end function sphericalSIDMCoreNFWConstructorParameters

  function sphericalSIDMCoreNFWConstructorInternal(factorRadiusCore,timeAge,nonAnalyticSolver,massDistribution_,darkMatterParticle_,componentType,massType) result(self)
    !!{
    Internal constructor for the \refClass{massDistributionSphericalSIDMCoreNFW} mass distribution class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    type            (massDistributionSphericalSIDMCoreNFW)                          :: self
    class           (massDistributionNFW                 ), intent(in   ), target   :: massDistribution_
    class           (darkMatterParticleClass             ), intent(in   ), target   :: darkMatterParticle_
    type            (enumerationNonAnalyticSolversType   ), intent(in   )           :: nonAnalyticSolver
    double precision                                      , intent(in   )           :: factorRadiusCore   , timeAge
    type            (enumerationComponentTypeType        ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType             ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="factorRadiusCore, timeAge, nonAnalyticSolver, componentType, massType, *massDistribution_, *darkMatterParticle_"/>
    !!]

    ! Validate the dark matter particle type.
    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! This is as expected.
    class default
       call Error_Report('this class expects a self-interacting dark matter particle'//{introspection:location})
    end select
    ! Initialize state.
    self%dimensionless=.false.
    return
  end function sphericalSIDMCoreNFWConstructorInternal

  subroutine sphericalSIDMCoreNFWDestructor(self)
    !!{
    Destructor for the abstract \refClass{massDistributionSphericalSIDMCoreNFW} class.
    !!}
    implicit none
    type(massDistributionSphericalSIDMCoreNFW), intent(inout) :: self
    
    !![
    <objectDestructor name="self%massDistribution_"  />
    <objectDestructor name="self%darkMatterParticle_"/>
    !!]
    return
  end subroutine sphericalSIDMCoreNFWDestructor

  double precision function sphericalSIDMCoreNFWRadiusCore(self) result(radiusCore)
    !!{
    Returns the core radius (in Mpc) of the ``coreNFW'' approximation to the self-interacting dark matter profile.
    !!}
    implicit none
    class(massDistributionSphericalSIDMCoreNFW), intent(inout) :: self

    radiusCore=+self%factorRadiusCore    &
         &     *self%radiusInteraction()
    return
  end function sphericalSIDMCoreNFWRadiusCore

  double precision function sphericalSIDMCoreNFWDensity(self,coordinates) result(density)
    !!{
    Compute the density at the specified {\normalfont \ttfamily coordinates} for the {\normalfont \ttfamily sidmCoreNFW}
    mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalSIDMCoreNFW), intent(inout) :: self
    class           (coordinate                          ), intent(in   ) :: coordinates
    double precision                                      , parameter     :: radiusFractionalLarge=10.0d0
    double precision                                                      :: radiusFractional            , radiusCore

    radiusCore      =+self       %radiusCore()
    radiusFractional=+coordinates%rSpherical() &
         &            /           radiusCore
    if (radiusFractional < radiusFractionalLarge) then
       ! Use the full solution for sufficiently small radii.
       density=+self%massDistribution_%density     (coordinates             ) &
            &  *tanh(                                                         &
            &        +radiusFractional                                        &
            &       )                                                         &
            &  +self%massDistribution_%massEnclosedBySphere(coordinates%rSpherical()) &
            &  /4.0d0                                                         &
            &  /Pi                                                            &
            &  /      radiusFractional**2                                     &
            &  /      radiusCore      **3                                     &
            &  /cosh(                                                         &
            &        +radiusFractional                                        &
            &       )**2
    else
       ! For large fractional radii avoid floating point overflow by approximating cosh(x) ~ 1/2/exp(-x).
       density=+self%massDistribution_%density      (coordinates            ) &
            &  *tanh(                                                         &
            &        +radiusFractional                                        &
            &       )                                                         &
            &  +self%massDistribution_%massEnclosedBySphere(coordinates%rSpherical()) &
            &  /4.0d0                                                         &
            &  /Pi                                                            &
            &  /      radiusFractional**2                                     &
            &  /      radiusCore      **3                                     &
            &  *4.0d0                                                         &
            &  *exp(                                                          &
            &       -2.0d0                                                    &
            &       * radiusFractional                                        &
            &      )
    end if
    return
  end function sphericalSIDMCoreNFWDensity

  double precision function sphericalSIDMCoreNFWDensityGradientRadial(self,coordinates,logarithmic) result(densityGradient)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a truncated spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalSIDMCoreNFW), intent(inout), target   :: self
    class           (coordinate                          ), intent(in   )           :: coordinates
    logical                                               , intent(in   ), optional :: logarithmic
    double precision                                                                :: radiusCore , massEnclosedNFW   , &
         &                                                                             densityNFW , densityLogSlopeNFW, &
         &                                                                             radius
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]

    radius            =+coordinates           %rSpherical     (                                           )
    radiusCore        =+self                  %radiusCore     (                                           )
    massEnclosedNFW   =+self%massDistribution_%massEnclosedBySphere   (coordinates%rSpherical()                   )
    densityNFW        =+self%massDistribution_%density        (coordinates                                )
    densityLogSlopeNFW=+self%massDistribution_%densityGradientRadial(coordinates             ,logarithmic=.true.)
    densityGradient   =+(                          &
         &               -2.0d0                    &
         &               *massEnclosedNFW          &
         &               *(                        &
         &                 +radiusCore             &
         &                 +radius                 &
         &                 *tanh(                  &
         &                       +radius           &
         &                       /radiusCore       &
         &                      )                  &
         &                )                        &
         &               +radiusCore*radius        &
         &               *(                        &
         &                 +4.0d0                  &
         &                 *Pi                     &
         &                 *radius**2              &
         &                 *densityNFW             &
         &                 +2.0d0                  &
         &                 *Pi                     &
         &                 *radius**2              &
         &                 *(                      &
         &                   +2.0d0                &
         &                   *densityNFW           &
         &                   +(                    &
         &                     +radiusCore         &
         &                     *sinh(              &
         &                           +2.0d0        &
         &                           *radius       &
         &                           /radiusCore   &
         &                         )               &
         &                     *densityLogSlopeNFW &
         &                     *densityNFW         &
         &                    )                    &
         &                 /radius                 &
         &                )                        &
         &               )                         &
         &              )                          &
         &             /(                          &
         &               +radiusCore               &
         &               *(                        &
         &                 +massEnclosedNFW        &
         &                 +2.0d0                  &
         &                 *Pi                     &
         &                 *radiusCore             &
         &                 *radius    **2          &
         &                 *sinh(                  &
         &                       +2.0d0            &
         &                       *radius           &
         &                       /radiusCore       &
         &                      )                  &
         &                 *densityNFW             &
         &                )                        &
         &              )
    if (.not.logarithmic_)                                    &
         & densityGradient=+     densityGradient              &
         &                 *self%density        (coordinates) &
         &                 /     radius
    return
  end function sphericalSIDMCoreNFWDensityGradientRadial
  
  double precision function sphericalSIDMCoreNFWMassEnclosedBySphere(self,radius) result(mass)
    !!{   
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for the {\normalfont \ttfamily sidmCoreNFW}
    mass distribution.
    !!}
    implicit none
    class           (massDistributionSphericalSIDMCoreNFW), intent(inout), target :: self
    double precision                                      , intent(in   )         :: radius

    mass   =+self%massDistribution_%massEnclosedBySphere(radius) &
         &  *tanh(                                               &
         &        +     radius                                   &
         &        /self%radiusCore()                             &
         &       )
    return
  end function sphericalSIDMCoreNFWMassEnclosedBySphere
