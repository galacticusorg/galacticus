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

  !+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

  !!{
  Implements a black hole binary separation growth class which follows a modified version of \cite{volonteri_assembly_2003},
  including terms for dynamical friction, hardening due to scattering of stars and emission of gravitational waves.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <blackHoleBinarySeparationGrowthRate name="blackHoleBinarySeparationGrowthRateStandard">
   <description>
    A black hole binary separation growth class that computes the separation growth rate of the binaries following a modified
    version of \cite {volonteri_assembly_2003} which include terms for dynamical friction, hardening due to scattering of stars
    and gravitational wave emission.
    \begin{equation}
    \dot{a} = \hbox{min} \left( - \frac{\mathrm{G}\rho _{*}a^2 H}{\sigma}, +\frac{2 \dot{v}_\mathrm{DF} a}{v_c} \right) -
    \frac{256 G^3 M_{\bullet, 1} M_{\bullet, 2} (M_{\bullet, 1} +M_{\bullet, 2})}{5 c^5 a^3}
    \end{equation}
    where $a$ is the black hole binary separation, $H$ is a dimensionless hardening parameter $H\approx 15$ in the limit of a
    very hard, equal mass binary, $\rho _\star$ is the density of stars, $\dot{v}_\mathrm{DF}$ is the acceleration (negative)
    due to dynamical friction, $v_\mathrm{c}$ is the circular velocity, $\sigma$ is the velocity dispersion of stars. Here the
    first factor represents hardening due to strong scattering of stars, the second results from dynamical friction with
    distant stars, gas and dark matter and the last results from the emission of gravitational waves
    \cite{peters_gravitational_1964}.
  
    The acceleration due to dynamical friction is computed using Chandrasekhar's formula:
    \begin{equation}
     \dot{v}_\mathrm{DF}=- {2 \pi \mathrm{G}^2 M_\bullet \over V_\mathrm{C}^2} \sum_{i} \rho_i \log(1+\Lambda_i^2) \left[
     \hbox{erf}(X_i)-\left\{ {2 X_i \over \sqrt{\pi}} \exp\left(-X_i^2\right) \right\} \right],
    \end{equation}
    where the sum is taken over the spheroid (gaseous plus stellar mass) and dark matter halo components\footnote{The disk is
    ignored as the black hole is assumed to be orbiting in a circular orbit in the disk.}. Here,
    \begin{equation}
    \Lambda_i =  {a \sigma^2  \over \mathrm{G}(M_{\bullet, 1}+M_{\bullet, 2})},
    \end{equation}
    is the Coulomb logarithm and
    \begin{equation}
    X_i = V_\mathrm{c} / \sqrt{2} \sigma.
    \end{equation}
    In all of the above equations, the velocity dispersion $\sigma_i$ is computed from the spherical Jeans equation assuming an
    isotropic velocity dispersion if {\normalfont \ttfamily [computeVelocityDispersion]}$=${\normalfont \ttfamily
    true}. Otherwise, $\sigma_i$ is set to the halo virial velocity for dark matter and to the spheroid characteristic velocity
    for the spheroid.
    
    In calculating the rate of hardening due to scattering of stars, the stellar density is reduced by a factor
    \citep{volonteri_assembly_2003}
    \begin{equation}
    f_\rho = \hbox{min}\left\{ \left[ { 4 a \sigma_\mathrm{spheroid}^2 \over 3 \mathrm{G} (M_{\bullet, 1}+M_{\bullet, 2})}
    \log\left({\mathrm{G} M_{\bullet, 2} \over 4 \sigma_\mathrm{spheroid}^2 a }\right) \right]^2 , 1 \right\},
    \end{equation}
    if {\normalfont \ttfamily [stellarDensityChangeBinaryMotion]}$=${\normalfont \ttfamily true} to account for the ejection of
    stars from the loss cone.
   </description>
  </blackHoleBinarySeparationGrowthRate>
  !!]
  type, extends(blackHoleBinarySeparationGrowthRateClass) :: blackHoleBinarySeparationGrowthRateStandard
     !!{
     A black hole binary separation growth class which follows a modified version of \cite{volonteri_assembly_2003}, including
     terms for dynamical friction, hardening due to scattering of stars and emission of gravitational waves.
     !!}
     private
     logical                                    :: stellarDensityChangeBinaryMotion          , computeVelocityDispersion
     class  (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_             => null()
   contains
     final     ::               standardDestructor
     procedure :: growthRate => standardGrowthRate
  end type blackHoleBinarySeparationGrowthRateStandard

  interface blackHoleBinarySeparationGrowthRateStandard
     !!{
     Constructors for the \refClass{blackHoleBinarySeparationGrowthRateStandard} black hole binary recoil class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface blackHoleBinarySeparationGrowthRateStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{blackHoleBinarySeparationGrowthRateStandard} black hole binary separation growth rate class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (blackHoleBinarySeparationGrowthRateStandard)                :: self
    type   (inputParameters                            ), intent(inout) :: parameters
    class  (darkMatterHaloScaleClass                   ), pointer       :: darkMatterHaloScale_
    logical                                                             :: stellarDensityChangeBinaryMotion, computeVelocityDispersion

    !![
    <inputParameter>
      <name>stellarDensityChangeBinaryMotion</name>
      <defaultValue>.true.</defaultValue>
      <description>The change in density due to the black hole's motion.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>computeVelocityDispersion</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the velocity dispersion of dark matter and stars should be computed using Jeans equation
         in black hole binary hardening calculations. If {\normalfont \ttfamily false}, then the velocity dispersions are assumed to equal
         the characteristic velocity of dark matter and spheroid.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=blackHoleBinarySeparationGrowthRateStandard(stellarDensityChangeBinaryMotion,computeVelocityDispersion,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(stellarDensityChangeBinaryMotion,computeVelocityDispersion,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{blackHoleBinarySeparationGrowthRateStandard} black hole binary separation growth class.
    !!}
    implicit none
    type   (blackHoleBinarySeparationGrowthRateStandard)                        :: self
    class  (darkMatterHaloScaleClass                   ), intent(in   ), target :: darkMatterHaloScale_
    logical                                             , intent(in   )         :: stellarDensityChangeBinaryMotion, computeVelocityDispersion
    !![
    <constructorAssign variables="stellarDensityChangeBinaryMotion,computeVelocityDispersion,*darkMatterHaloScale_"/>
    !!]

    return
  end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{blackHoleBinarySeparationGrowthRateStandard} black hole binary separation growth class.
    !!}
    implicit none
    type(blackHoleBinarySeparationGrowthRateStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine standardDestructor

  double precision function standardGrowthRate(self,blackHole)
    !!{
    Returns an initial separation growth rate for a binary black holes that follows a modified version of
    \cite{volonteri_assembly_2003}.
    !!}
    use :: Coordinates                     , only : coordinateCylindrical , assignment(=)
    use :: Display                         , only : displayIndent         , displayMessage                , displayUnindent
    use :: Galactic_Structure_Options      , only : componentTypeDarkHalo , componentTypeSpheroid         , massTypeDark   , massTypeGalactic,&
          &                                         massTypeStellar
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentBlackHole, nodeComponentSpheroid         , treeNode
    use :: Mass_Distributions              , only : massDistributionClass , kinematicsDistributionClass
    use :: Numerical_Constants_Astronomical, only : MpcPerKmPerSToGyr     , gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Prefixes    , only : milli
    implicit none
    class           (blackHoleBinarySeparationGrowthRateStandard), intent(inout) :: self
    class           (nodeComponentBlackHole                     ), intent(inout) :: blackHole
    type            (treeNode                                   ), pointer       :: node
    class           (nodeComponentBlackHole                     ), pointer       :: blackHoleCentral
    class           (nodeComponentSpheroid                      ), pointer       :: spheroid
    class           (massDistributionClass                      ), pointer       :: massDistributionSpheroidStellar_             , massDistributionDarkMatterHalo_      , &
         &                                                                          massDistributionGalactic_                    , massDistribution_
    class           (kinematicsDistributionClass                ), pointer       :: kinematicsDistributionSpheroidStellar_       , kinematicsDistributionDarkMatterHalo_
    double precision                                             , parameter     :: hardeningRateDimensionless            =15.0d0
    double precision                                             , parameter     :: outerRadiusMultiplier                 =10.0d0
    double precision                                             , parameter     :: dynamicalFrictionMinimumRadius        = 0.1d0
    double precision                                                             :: coulombLogarithmDarkMatter                   , coulombLogarithmSpheroid             , &
         &                                                                          densityDarkMatter                            , densitySpheroid                      , &
         &                                                                          densityStellar                               , dynamicalFrictionAcceleration        , &
         &                                                                          dynamicalFrictionXDarkMatter                 , dynamicalFrictionXSpheroid           , &
         &                                                                          radiusHardBinary                             , rateGravitationalWaves               , &
         &                                                                          rateScattering                               , rateScatteringDynamicalFriction      , &
         &                                                                          rateScatteringStars                          , rotationCurveGradient                , &
         &                                                                          stellarDensityFractionRemaining              , velocityDispersionDarkMatter         , &
         &                                                                          velocityDispersionSpheroid                   , velocityRotation
    type            (coordinateCylindrical                      )                :: coordinates
    character       (len=24                                     )                :: message

    ! Get the host node.
    node             => blackHole%host()
    ! Get the primary (central) black hole of this node.
    blackHoleCentral => node%blackHole(instance=1)
    ! Return a zero separation growth rate if the black hole has non-positive radial position or either black hole (active or
    ! central) has negative mass.
    if     (                                        &
         &   blackHole   %radialPosition() <= 0.0d0 &
         &  .or.                                    &
         &   blackHole   %mass          () <= 0.0d0 &
         &  .or.                                    &
         &   blackHoleCentral%mass      () <= 0.0d0 &
         & ) then
       standardGrowthRate=0.0d0
       return
    end if
    ! Get required mass distributions.
    massDistribution_                => node%massDistribution(                                                             )
    massDistributionSpheroidStellar_ => node%massDistribution(componentType=componentTypeSpheroid,massType=massTypeStellar )
    massDistributionDarkMatterHalo_  => node%massDistribution(componentType=componentTypeDarkHalo,massType=massTypeDark    )
    massDistributionGalactic_        => node%massDistribution(                                    massType=massTypeGalactic)
    ! Get the spheroid component.
    spheroid                         => node%spheroid        (                                                             )
    ! Set coordinates of the black hole.
    coordinates                      =  [blackHole%radialPosition(),0.0d0,0.0d0]
   ! Compute the velocity dispersion of stars and dark matter.
    if (self%computeVelocityDispersion) then
       kinematicsDistributionSpheroidStellar_ => massDistributionSpheroidStellar_      %kinematicsDistribution(                                               )
       kinematicsDistributionDarkMatterHalo_  => massDistributionDarkMatterHalo_       %kinematicsDistribution(                                               )
       velocityDispersionSpheroid             =  kinematicsDistributionSpheroidStellar_%velocityDispersion1D  (coordinates,massDistribution_,massDistribution_)
       velocityDispersionDarkMatter           =  kinematicsDistributionDarkMatterHalo_ %velocityDispersion1D  (coordinates,massDistribution_,massDistribution_)
       !![
       <objectDestructor name="kinematicsDistributionSpheroidStellar_"/>
       <objectDestructor name="kinematicsDistributionDarkMatterHalo_" />
       !!]
    else
       velocityDispersionSpheroid            =  spheroid                     %velocity      (    )
       velocityDispersionDarkMatter          =  self    %darkMatterHaloScale_%velocityVirial(node)
    end if
    ! Compute the separation growth rate due to emission of gravitational waves.
    rateGravitationalWaves=-(                                        &
         &                   +256.0d0                                &
         &                   *gravitationalConstant_internal**3      &
         &                   *  blackHole       %mass          ()    &
         &                   *  blackHoleCentral%mass          ()    &
         &                   *(                                      &
         &                     +blackHole       %mass          ()    &
         &                     +blackHoleCentral%mass          ()    &
         &                    )                                      &
         &                  )                                        &
         &                 /(                                        &
         &                   +5.0d0                                  &
         &                   *(speedLight*milli)**5                  &
         &                   *  blackHole       %radialPosition()**3 &
         &                  )                                        &
         &                 /MpcPerKmPerSToGyr
    ! Compute the hard binary radius, where shrinking of the binary switches from dynamical friction to hardening due to strong
    ! scattering of individual stars.
    if (velocityDispersionSpheroid > 0.0d0) then
       radiusHardBinary=+(                                &
            &             +gravitationalConstant_internal &
            &             *(                              &
            &               +blackHole       %mass()      &
            &               +blackHoleCentral%mass()      &
            &              )                              &
            &            )                                &
            &           /(                                &
            &             +4.0d0                          &
            &             *velocityDispersionSpheroid**2  &
            &            )
    else
       radiusHardBinary=0.0d0
    end if
    ! First, check if the change in stellar density due to the binary's inward motion is to be computed.
    stellarDensityFractionRemaining=1.0d0
    ! If it does change, we first compute the fraction of that change according to Volonteri et al. (2003) else we set it as the
    ! normal density.
    if (self%stellarDensityChangeBinaryMotion .and. velocityDispersionSpheroid > 0.0d0)                                         &
         &  stellarDensityFractionRemaining=(+                                                blackHole       %radialPosition() &
         &                                   *                                                velocityDispersionSpheroid**2     &
         &                                   *(4.0d0/3.0d0)/(gravitationalConstant_internal*(+blackHoleCentral%mass          () &
         &                                                                                   +blackHole       %mass          () &
         &                                                                                  )                                   &
         &                                   *log                                           (                                   &
         &                                                   gravitationalConstant_internal*  blackHole       %mass          () &
         &                                                                                   /4.0d0                             &
         &                                                                                   /velocityDispersionSpheroid**2     &
         &                                                                                   /blackHole       %radialPosition() &
         &                                                                                  )                                   &
         &                                                  )                                                                   &
         &                                  )**2
    ! Limit the density fraction to unity.
    stellarDensityFractionRemaining=min(stellarDensityFractionRemaining,1.0d0)
    ! Compute the stellar density, accounting for any loss.
    densityStellar=+massDistributionSpheroidStellar_%density(coordinates) &
         &         *stellarDensityFractionRemaining
    ! Compute the hardening rate due to strong scattering of individual stars.
    if (velocityDispersionSpheroid > 0.0d0) then
       rateScatteringStars=-gravitationalConstant_internal &
            &              *densityStellar                 &
            &              *blackHole%radialPosition()**2  &
            &              *hardeningRateDimensionless     &
            &              /velocityDispersionSpheroid     &
            &              /MpcPerKmPerSToGyr
    else
       rateScatteringStars=0.0d0
    end if
    ! Check if the binary has sufficiently large separation that we should compute the rate of hardening due to dynamical friction.
    if (blackHole%radialPosition() > dynamicalFrictionMinimumRadius*radiusHardBinary) then
       ! Compute the total density, including dark matter.
       densitySpheroid  =massDistributionGalactic_      %density(coordinates)
       densityDarkMatter=massDistributionDarkMatterHalo_%density(coordinates)
       ! Compute the Coulomb logarithms for dynamical friction.
       coulombLogarithmSpheroid  = (                                     &
            &                          blackHole       %radialPosition() &
            &                       *  velocityDispersionSpheroid**2     &
            &                      )                                     &
            &                     /(                                     &
            &                        gravitationalConstant_internal      &
            &                       *(                                   &
            &                          blackHole       %mass          () &
            &                         +blackHoleCentral%mass          () &
            &                        )                                   &
            &                      )
       coulombLogarithmDarkMatter= (                                     &
            &                          blackHole       %radialPosition() &
            &                       *  velocityDispersionDarkMatter**2   &
            &                      )                                     &
            &                     /(                                     &
            &                        gravitationalConstant_internal      &
            &                       *(                                   &
            &                          blackHole       %mass          () &
            &                         +blackHoleCentral%mass          () &
            &                        )                                   &
            &                      )
       ! Compute the rotation curve of the galaxy and the additional contribution from the active black hole. Add them in
       ! quadrature to get an estimate of the actual orbital speed of the black hole binary.
       ! Precompute the "X" term appearing in the dynamical friction formula.
       velocityRotation=massDistribution_%rotationCurve(blackHole%radialPosition())
       if (velocityDispersionSpheroid > 0.0d0) then
          dynamicalFrictionXSpheroid=+velocityRotation                                                       &
               &                     /sqrt(2.0d0)                                                            &
               &                     /velocityDispersionSpheroid
       else
          dynamicalFrictionXSpheroid=+0.0d0
       end if
       dynamicalFrictionXDarkMatter =+velocityRotation                                                       &
            &                        /sqrt(2.0d0)                                                            &
            &                        /velocityDispersionDarkMatter
       ! Compute the acceleration due to dynamical friction.
       dynamicalFrictionAcceleration=-4.0d0                                                                  &
            &                        *Pi                                                                     &
            &                        *gravitationalConstant_internal**2                                      &
            &                        *blackHole%mass()                                                       &
            &                        *0.5d0                                                                  &
            &                        /velocityRotation**2                                                    &
            &                        /MpcPerKmPerSToGyr                                                      &
            &                        *(                                                                      &
            &                           densitySpheroid                                                      &
            &                          *log(1.0d0+coulombLogarithmSpheroid**2)                               &
            &                          *(                                                                    &
            &                             erf(dynamicalFrictionXSpheroid)                                    &
            &                            -(                                                                  &
            &                               2.0d0                                                            &
            &                              *dynamicalFrictionXSpheroid                                       &
            &                              /sqrt(Pi)                                                         &
            &                              *exp(-dynamicalFrictionXSpheroid**2)                              &
            &                             )                                                                  &
            &                           )                                                                    &
            &                          +densityDarkMatter                                                    &
            &                          *log(1.0d0+coulombLogarithmDarkMatter**2)                             &
            &                          *(                                                                    &
            &                             erf(dynamicalFrictionXDarkMatter)                                  &
            &                            -(                                                                  &
            &                               2.0d0                                                            &
            &                              *dynamicalFrictionXDarkMatter                                     &
            &                              /sqrt(Pi)                                                         &
            &                              *exp(-dynamicalFrictionXDarkMatter**2)                            &
            &                             )                                                                  &
            &                           )                                                                    &
            &                         )
       ! Compute the radial inflow velocity due to dynamical friction.
       rotationCurveGradient          =(                                                                     &
            &                           +                  velocityRotation                                  &
            &                           +0.5d0                                                               &
            &                           *                                        blackHole%radialPosition()  &
            &                           /                  velocityRotation                                  &
            &                           *massDistribution_%rotationCurveGradient(blackHole%radialPosition()) &
            &                          )
       if (rotationCurveGradient == 0.0d0) then
          call displayIndent('dynamical friction calculation report')
          write (message,'(a,i12)  ') 'nodeIndex = ',node%index()
          write (message,'(a,e12.6)') '  V (r)     = ',                  velocityRotation
          call displayMessage(trim(message))
          write (message,'(a,e12.6)') '     r     = ',                                                blackHole%radialPosition()
          call displayMessage(trim(message))
          write (message,'(a,e12.6)') ' dV²(r)/dr = ',massDistribution_%rotationCurveGradient        (blackHole%radialPosition())
          call displayMessage(trim(message))
          write (message,'(a,e12.6)') '    a_{df} = ',                  dynamicalFrictionAcceleration
          call displayMessage(trim(message))
          call displayUnindent('done')
          call Error_Report('rotation curve gradient is zero'//{introspection:location})
       end if
       rateScatteringDynamicalFriction=+2.0d0                         &
            &                          *dynamicalFrictionAcceleration &
            &                          *blackHole%radialPosition()    &
            &                          /rotationCurveGradient
       ! Take the most negative rate from dynamical friction or scattering of stars.
       rateScattering=min(rateScatteringDynamicalFriction,rateScatteringStars)
    else
       ! For sufficiently small radii assume that scattering of individual stars always dominates over dynamical friction.
       rateScattering=rateScatteringStars
    end if
    ! Sum the two contributions to the radial growth rate.
    standardGrowthRate=rateScattering+rateGravitationalWaves
    ! Clean up.
    !![
    <objectDestructor name="massDistribution_"               />
    <objectDestructor name="massDistributionSpheroidStellar_"/>
    <objectDestructor name="massDistributionDarkMatterHalo_" />
    <objectDestructor name="massDistributionGalactic_"       />
    !!]
    return
  end function standardGrowthRate
