!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  An implementation of a cored-NFW dark matter halo profile to approximate the effects of SIDM based on the model of Jiang et al. (2022).
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOSIDMCoreNFW">
   <description>Cored-NFW dark matter halo profile to approximate the effects of SIDM based on the model of Jiang et al. (2022).</description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOSIDM) :: darkMatterProfileDMOSIDMCoreNFW
     !!{
     A dark matter halo profile class implementing a cored-NFW dark matter halo profile to approximate the effects of SIDM based
     on the model of Jiang et al. (2022). The profile is defined by the enclosed mass, with (Jiang et al. 2022):
     \begin{equation}
       M(r) = M_\mathrm{NFW}(r) \mathrm{tanh}\left(\frac{r}{r_\mathrm{c}}\right),
     \end{equation}
     where $r_\mathrm{c} = \alpha r_1$ is a characteristic core size related to the interaction radius $r_1$ by a constant factor
     $\alpha = ${\normalfont \ttfamily [factorRadiusCore]}.
     !!}
     private
     double precision :: factorRadiusCore
   contains
     !![
     <methods>
       <method method="radiusCore"       description="Computes the core radius of halo."/>
       <method method="calculationReset" description="Reset memoized calculations."     />
     </methods>
     !!]
     final     ::                                      sidmCoreNFWDestructor
     procedure :: autoHook                          => sidmCoreNFWAutoHook
     procedure :: calculationReset                  => sidmCoreNFWCalculationReset
     procedure :: radiusCore                        => sidmCoreNFWRadiusCore
     procedure :: density                           => sidmCoreNFWDensity
     procedure :: densityLogSlope                   => sidmCoreNFWDensityLogSlope
     procedure :: radiusEnclosingDensity            => sidmCoreNFWRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => sidmCoreNFWRadiusEnclosingMass
     procedure :: radialMoment                      => sidmCoreNFWRadialMoment
     procedure :: enclosedMass                      => sidmCoreNFWEnclosedMass
     procedure :: potential                         => sidmCoreNFWPotential
     procedure :: circularVelocity                  => sidmCoreNFWCircularVelocity
     procedure :: circularVelocityMaximum           => sidmCoreNFWCircularVelocityMaximum
     procedure :: radiusCircularVelocityMaximum     => sidmCoreNFWRadiusCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => sidmCoreNFWRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => sidmCoreNFWRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => sidmCoreNFWRotationNormalization
     procedure :: energy                            => sidmCoreNFWEnergy
     procedure :: kSpace                            => sidmCoreNFWKSpace
     procedure :: freefallRadius                    => sidmCoreNFWFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => sidmCoreNFWFreefallRadiusIncreaseRate
  end type darkMatterProfileDMOSIDMCoreNFW

  interface darkMatterProfileDMOSIDMCoreNFW
     !!{
     Constructors for the {\normalfont \ttfamily sidmCoreNFW} dark matter halo profile class.
     !!}
     module procedure sidmCoreNFWConstructorParameters
     module procedure sidmCoreNFWConstructorInternal
  end interface darkMatterProfileDMOSIDMCoreNFW

contains

  function sidmCoreNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sidmCoreNFW} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOSIDMCoreNFW)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    class           (darkMatterParticleClass        ), pointer       :: darkMatterParticle_
    double precision                                                 :: factorRadiusCore

    !![
    <inputParameter>
      <name>factorRadiusCore</name>
      <defaultValue>0.45d0</defaultValue>
      <defaultSource>Jiang et al. (2022)</defaultSource>
      <source>parameters</source>
      <description>The factor $\alpha$ appearing in the definition of the core radius, $r_\mathrm{c}=\alpha r_1$ where $r_1$ is the radius at which an SIDM particle has had, on average, 1 interaction.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !!]
    self=darkMatterProfileDMOSIDMCoreNFW(factorRadiusCore,darkMatterHaloScale_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function sidmCoreNFWConstructorParameters

  function sidmCoreNFWConstructorInternal(factorRadiusCore,darkMatterHaloScale_,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sidmCoreNFW} dark matter profile class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    type            (darkMatterProfileDMOSIDMCoreNFW)                        :: self
    class           (darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterParticleClass        ), intent(in   ), target :: darkMatterParticle_
    double precision                                 , intent(in   )         :: factorRadiusCore
    !![
    <constructorAssign variables="factorRadiusCore, *darkMatterHaloScale_, *darkMatterParticle_"/>
    !!]

    ! Validate the dark matter particle type.
    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! This is as expected.
    class default
       call Error_Report('transfer function expects a self-interacting dark matter particle'//{introspection:location})
    end select
    ! Construct an NFW profile.
    allocate(darkMatterProfileDMONFW :: self%darkMatterProfileDMO_)
    select type (darkMatterProfileDMO_ => self%darkMatterProfileDMO_)
    type is (darkMatterProfileDMONFW)
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="darkMatterProfileDMO_"  object="darkMatterProfileDMO_">
	 <constructor>
	   darkMatterProfileDMONFW(                                                           &amp;
	    &amp;                  velocityDispersionUseSeriesExpansion=.true.              , &amp;
	    &amp;                  darkMatterHaloScale_                =darkMatterHaloScale_  &amp;
	    &amp;                 )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    self%genericLastUniqueID =-1_kind_int8
    self%uniqueIDPreviousSIDM=-1_kind_int8
    return
  end function sidmCoreNFWConstructorInternal

  subroutine sidmCoreNFWAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self

    call calculationResetEvent%attach(self,sidmCoreNFWCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOSIDMCoreNFW')
    return
  end subroutine sidmCoreNFWAutoHook

  subroutine sidmCoreNFWDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily sidmCoreNFW} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterParticle_"  />
    !!]
    if (calculationResetEvent%isAttached(self,sidmCoreNFWCalculationReset)) call calculationResetEvent%detach(self,sidmCoreNFWCalculationReset)
    return
  end subroutine sidmCoreNFWDestructor

  subroutine sidmCoreNFWCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    self%genericLastUniqueID                         =node%uniqueID()
    self%uniqueIDPreviousSIDM                        =node%uniqueID()
    self%radiusInteractivePrevious                   =-1.0d0
    self%genericEnclosedMassRadiusMinimum            =+huge(0.0d0)
    self%genericEnclosedMassRadiusMaximum            =-huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMinimum=+huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMaximum=-huge(0.0d0)
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine sidmCoreNFWCalculationReset

  double precision function sidmCoreNFWRadiusCore(self,node)
    !!{
    Returns the core radius (in Mpc) of the ``coreNFW'' approximation to the self-interacting dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    sidmCoreNFWRadiusCore=+self%factorRadiusCore        &
         &                *self%radiusInteraction(node)
    return
  end function sidmCoreNFWRadiusCore

  double precision function sidmCoreNFWDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius
    double precision                                 , parameter     :: radiusFractionalLarge=10.0d0
    double precision                                                 :: radiusFractional            , radiusCore

    radiusCore      =+self%radiusCore(node)
    radiusFractional=+     radius           &
         &            /    radiusCore
    if (radiusFractional < radiusFractionalLarge) then
       ! Use the full solution for sufficiently small radii.
       sidmCoreNFWDensity=+self%darkMatterProfileDMO_%density(node,radius)      &
            &             *tanh(                                                &
            &                   +radiusFractional                               &
            &                  )                                                &
            &             +self%darkMatterProfileDMO_%enclosedMass(node,radius) &
            &             /4.0d0                                                &
            &             /Pi                                                   &
            &             /      radiusFractional**2                            &
            &             /      radiusCore      **3                            &
            &             /cosh(                                                &
            &                   +radiusFractional                               &
            &             )**2
    else
       ! For large fractional radii avoid floating point overflow by approximating cosh(x) ~ 1/2/exp(-x).
       sidmCoreNFWDensity=+self%darkMatterProfileDMO_%density(node,radius)      &
            &             *tanh(                                                &
            &                   +radiusFractional                               &
            &                  )                                                &
            &             +self%darkMatterProfileDMO_%enclosedMass(node,radius) &
            &             /4.0d0                                                &
            &             /Pi                                                   &
            &             /      radiusFractional**2                            &
            &             /      radiusCore      **3                            &
            &             *4.0d0                                                &
            &             *exp(                                                 &
            &                  -2.0d0                                           &
            &                  * radiusFractional                               &
            &                 )
    end if
    return
  end function sidmCoreNFWDensity

  double precision function sidmCoreNFWDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius
    double precision                                                 :: radiusCore, massEnclosedNFW   , &
         &                                                              densityNFW, densityLogSlopeNFW

    radiusCore                =+self%radiusCore                           (node       )
    massEnclosedNFW           =+self%darkMatterProfileDMO_%enclosedMass   (node,radius)
    densityNFW                =+self%darkMatterProfileDMO_%density        (node,radius)
    densityLogSlopeNFW        =+self%darkMatterProfileDMO_%densityLogSlope(node,radius)
    sidmCoreNFWDensityLogSlope=+(                                     &
         &                       -2.0d0                               &
         &                       *massEnclosedNFW                     &
         &                       *(                                   &
         &                         +radiusCore                        &
         &                         +radius                            &
         &                         *tanh(                             &
         &                               +radius                      &
         &                               /radiusCore                  &
         &                              )                             &
         &                        )                                   &
         &                       +radiusCore*radius                   &
         &                       *(                                   &
         &                         +4.0d0                             &
         &                         *Pi                                &
         &                         *radius**2                         &
         &                         *densityNFW                        &
         &                         +2.0d0                             &
         &                         *Pi                                &
         &                         *radius**2                         &
         &                         *(                                 &
         &                           +2.0d0                           &
         &                           *densityNFW                      &
         &                           +(                               &
         &                             +radiusCore                    &
         &                             *sinh(                         &
         &                                   +2.0d0                   &
         &                                   *radius                  &
         &                                   /radiusCore              &
         &                                 )                          &
         &                             *densityLogSlopeNFW*densityNFW &
         &                            )                               &
         &                         /radius                            &
         &                        )                                   &
         &                       )                                    &
         &                      )                                     &
         &                     /(                                     &
         &                       +radiusCore                          &
         &                       *(                                   &
         &                         +massEnclosedNFW                   &
         &                         +2.0d0                             &
         &                         *Pi                                &
         &                         *radiusCore                        &
         &                         *radius    **2                     &
         &                         *sinh(                             &
         &                              +2.0d0                        &
         &                              *radius                       &
         &                              /radiusCore                   &
         &                        )                                   &
         &                       *densityNFW                          &
         &                      )                                     &
         &                     )
    return
  end function sidmCoreNFWDensityLogSlope

  double precision function sidmCoreNFWEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    sidmCoreNFWEnclosedMass=+      self%darkMatterProfileDMO_%enclosedMass(node,radius) &
         &                  *tanh(                                                      &
         &                        +                                             radius  &
         &                        /self%radiusCore                        (node       ) &
         &                       )
    return
  end function sidmCoreNFWEnclosedMass

  double precision function sidmCoreNFWRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout), target :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: density
    
    sidmCoreNFWRadiusEnclosingDensity=self%radiusEnclosingDensityNumerical(node,density)
    return
  end function sidmCoreNFWRadiusEnclosingDensity

  double precision function sidmCoreNFWRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout), target :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: mass

    sidmCoreNFWRadiusEnclosingMass=self%radiusEnclosingMassNumerical(node,mass)
    return
  end function sidmCoreNFWRadiusEnclosingMass

  double precision function sidmCoreNFWRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout)           :: self
    type            (treeNode                       ), intent(inout)           :: node
    double precision                                 , intent(in   )           :: moment
    double precision                                 , intent(in   ), optional :: radiusMinimum, radiusMaximum

    sidmCoreNFWRadialMoment=self%radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum)
    return
  end function sidmCoreNFWRadialMoment

  double precision function sidmCoreNFWPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW  ), intent(inout)           :: self
    type            (treeNode                         ), intent(inout), target   :: node
    double precision                                   , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status

    sidmCoreNFWPotential=self%potentialNumerical(node,radius,status)
    return
  end function sidmCoreNFWPotential

  double precision function sidmCoreNFWCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    sidmCoreNFWCircularVelocity=self%circularVelocityNumerical(node,radius)
    return
  end function sidmCoreNFWCircularVelocity

  double precision function sidmCoreNFWCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    sidmCoreNFWCircularVelocityMaximum=self%circularVelocityMaximumNumerical(node)
    return
  end function sidmCoreNFWCircularVelocityMaximum

  double precision function sidmCoreNFWRadiusCircularVelocityMaximum(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity is acheived in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    sidmCoreNFWRadiusCircularVelocityMaximum=self%radiusCircularVelocityMaximumNumerical(node)
    return
  end function sidmCoreNFWRadiusCircularVelocityMaximum

  double precision function sidmCoreNFWRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    sidmCoreNFWRadialVelocityDispersion=self%radialVelocityDispersionNumerical(node,radius)
    return
  end function sidmCoreNFWRadialVelocityDispersion

  double precision function sidmCoreNFWRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: specificAngularMomentum

    sidmCoreNFWRadiusFromSpecificAngularMomentum=self%radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum)
    return
  end function sidmCoreNFWRadiusFromSpecificAngularMomentum

  double precision function sidmCoreNFWRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    sidmCoreNFWRotationNormalization=self%rotationNormalizationNumerical(node)
    return
  end function sidmCoreNFWRotationNormalization

  double precision function sidmCoreNFWEnergy(self,node)
    !!{
    Return the energy of a sidmCoreNFW halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node

    sidmCoreNFWEnergy=self%energyNumerical(node)
    return
  end function sidmCoreNFWEnergy

  double precision function sidmCoreNFWKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the sidmCoreNFW density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout)         :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: waveNumber

    sidmCoreNFWKSpace=self%kSpaceNumerical(node,waveNumber)
    return
  end function sidmCoreNFWKSpace

  double precision function sidmCoreNFWFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the sidmCoreNFW density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout), target :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: time

    sidmCoreNFWFreefallRadius=self%freefallRadiusNumerical(node,time)
    return
  end function sidmCoreNFWFreefallRadius

  double precision function sidmCoreNFWFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the sidmCoreNFW density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMOSIDMCoreNFW), intent(inout), target :: self
    type            (treeNode                       ), intent(inout), target :: node
    double precision                                 , intent(in   )         :: time

    sidmCoreNFWFreefallRadiusIncreaseRate=self%freefallRadiusIncreaseRateNumerical(node,time)
    return
  end function sidmCoreNFWFreefallRadiusIncreaseRate
