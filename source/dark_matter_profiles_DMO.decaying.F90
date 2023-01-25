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
  An implementation of decaying dark matter halo profiles.
  !!}
  use :: Dark_Matter_Particles       , only : darkMatterParticleClass
  use :: Statistics_Distributions    , only : distributionFunction1DNonCentralChi

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMODecaying">
   <description>Decaying dark matter halo profiles.</description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMODecaying
     !!{
     A dark matter halo profile class implementing decaying dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_              => null()
     class           (darkMatterParticleClass  ), pointer :: darkMatterParticle_                => null()
     double precision                                     :: lifetime_, massSplitting_
     logical                                              :: massLoss_
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset"/>
     </methods>
     !!]
     final     ::                                      decayingDestructor
     procedure :: autoHook                          => decayingAutoHook
     procedure :: calculationReset                  => decayingCalculationReset
     procedure :: density                           => decayingDensity
     procedure :: densityLogSlope                   => decayingDensityLogSlope
     procedure :: radiusEnclosingDensity            => decayingRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => decayingRadiusEnclosingMass
     procedure :: radialMoment                      => decayingRadialMoment
     procedure :: enclosedMass                      => decayingEnclosedMass
     procedure :: potential                         => decayingPotential
     procedure :: circularVelocity                  => decayingCircularVelocity
     procedure :: radiusCircularVelocityMaximum     => decayingRadiusCircularVelocityMaximum
     procedure :: circularVelocityMaximum           => decayingCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => decayingRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => decayingRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => decayingRotationNormalization
     procedure :: energy                            => decayingEnergy
     procedure :: kSpace                            => decayingKSpace
     procedure :: freefallRadius                    => decayingFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => decayingFreefallRadiusIncreaseRate
     procedure :: decayingFactor                    => decayingDecayingFactor
     procedure :: escapeFraction                    => decayingEscapeFraction
  end type darkMatterProfileDMODecaying

  interface darkMatterProfileDMODecaying
     !!{
     Constructors for the {\normalfont \ttfamily decaying} dark matter halo profile class.
     !!}
     module procedure decayingConstructorParameters
     module procedure decayingConstructorInternal
  end interface darkMatterProfileDMODecaying

contains

  function decayingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily decaying} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMODecaying)                 :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass    ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    class           (darkMatterParticleClass      ), pointer       :: darkMatterParticle_
    double precision                                               :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, &
         &                                                            toleranceRelativePotential

    !![
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"   name="darkMatterProfileDMO_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMODecaying(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,toleranceRelativePotential,darkMatterParticle_,darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"  />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="darkMatterParticle_"  />
    !!]
    return
  end function decayingConstructorParameters

  function decayingConstructorInternal(toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,toleranceRelativePotential,darkMatterParticle_,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily decaying} dark matter profile class.
    !!}
    use :: Error                , only : Error_Report
    use :: Dark_Matter_Particles, only : darkMatterParticleDecayingDarkMatter
    implicit none
    type            (darkMatterProfileDMODecaying)                        :: self
    class           (darkMatterProfileDMOClass   ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass    ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterParticleClass     ), intent(in   ), target :: darkMatterParticle_
    double precision                                   , intent(in   )    :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, &
         &                                                                   toleranceRelativePotential
   !![
    <constructorAssign variables="toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, toleranceRelativePotential, *darkMatterParticle_, *darkMatterProfileDMO_,*darkMatterHaloScale_"/>
    !!]
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime_ = darkMatterParticle_%lifetime()
       self%massSplitting_ = darkMatterParticle_%massSplitting()
       self%massLoss_ = darkMatterParticle_%massLoss()
    class default
       ! No decays.
       self%lifetime_=-1.0d0
       self%massSplitting_=0.0d0
       self%massLoss_=.false.
    end select
    ! Initialize.
    self%genericLastUniqueID=-1_kind_int8
   return
  end function decayingConstructorInternal

  subroutine decayingAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self

    call calculationResetEvent%attach(self,decayingCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine decayingAutoHook

  subroutine decayingDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily decaying} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMODecaying), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterParticle_"/>
    !!]
    return
  end subroutine decayingDestructor

  subroutine decayingCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    ! Reset calculations for this profile.
    self%genericLastUniqueID                         =node%uniqueID()
    self%genericEnclosedMassRadiusMinimum            =+huge(0.0d0)
    self%genericEnclosedMassRadiusMaximum            =-huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMinimum=+huge(0.0d0)
    self%genericVelocityDispersionRadialRadiusMaximum=-huge(0.0d0)
    if (allocated(self%genericVelocityDispersionRadialVelocity)) deallocate(self%genericVelocityDispersionRadialVelocity)
    if (allocated(self%genericVelocityDispersionRadialRadius  )) deallocate(self%genericVelocityDispersionRadialRadius  )
    if (allocated(self%genericEnclosedMassMass                )) deallocate(self%genericEnclosedMassMass                )
    if (allocated(self%genericEnclosedMassRadius              )) deallocate(self%genericEnclosedMassRadius              )
    return
  end subroutine decayingCalculationReset

  subroutine decayingDecayingFactor(self, node, radius, factor)
    !!{
    Return the change in mass factor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    double precision                              , intent(  out) :: factor
    class           (nodeComponentBasic          ), pointer       :: basic
    double precision                                              :: fe
    basic             => node%basic()
    if (self%massLoss_) then
      call self%escapeFraction(radius, node, fe)
      factor = +1.0d0 - (fe - self%massSplitting_ * (1 - fe)) &
        &    *(+1.0d0 - exp(-basic%time() / self%lifetime_))
    else
      factor = +1.0d0
    end if
    return
  end subroutine decayingDecayingFactor

  subroutine decayingEscapeFraction(self, radius, node, fraction)
    !!{
    Return the fraction of the particles that escape at a given radius.
    !!}
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class            (darkMatterProfileDMODecaying       ), intent(inout) :: self
    type            (treeNode                    ), intent(inout), target :: node
    type             (distributionFunction1DNonCentralChi)                :: nonCentralChi
    double precision                                      , intent(in   ) :: radius
    double precision                                      , intent(  out) :: fraction
    double precision                                                      :: sigma, vEsc
    sigma = self%darkMatterProfileDMO_%radialVelocityDispersion(node, radius)
    if (sigma > 0.0d0) then
      vEsc = +2.0d0 * abs(self%darkMatterProfileDMO_%potential(node, radius) - self%darkMatterProfileDMO_%potential(node, +1.0d3 * self%darkMatterHaloScale_%radiusVirial(node))) !! Check sign
      nonCentralChi = distributionFunction1DNonCentralChi(       &
                    & (self%massSplitting_/(speedLight/kilo))**2 &
                    & / sigma**2)
      fraction = +1.0d0 - nonCentralChi%cumulative(vEsc / sigma**2)
    else
      fraction = +1.0d0
    end if
    return
  end subroutine decayingEscapeFraction


  double precision function decayingDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    double precision                                              :: factor
    
    call self%decayingFactor(node, radius, factor)
    decayingDensity=+self%darkMatterProfileDMO_%density(node,radius) * factor
    return
  end function decayingDensity
  
  double precision function decayingDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    
    decayingDensityLogSlope=self%densityLogSlopeNumerical(node, radius)
    return
  end function decayingDensityLogSlope

  
  double precision function decayingRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: density
    
    decayingRadiusEnclosingDensity = self%radiusEnclosingDensityNumerical(node, density)
    return
  end function decayingRadiusEnclosingDensity

  double precision function decayingRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: mass

    decayingRadiusEnclosingMass=self%radiusEnclosingMassNumerical(node, mass)
    return
  end function decayingRadiusEnclosingMass

  double precision function decayingRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum, radiusMaximum
    
    decayingRadialMoment=self%radialMomentNumerical(node, moment, radiusMinimum, radiusMaximum)
    return
  end function decayingRadialMoment

  double precision function decayingEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    
    decayingEnclosedMass=+self%enclosedMassNumerical(node, radius)
    return
  end function decayingEnclosedMass

  double precision function decayingPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout)           :: self
    type            (treeNode                    ), intent(inout), target   :: node
    double precision                              , intent(in   )           :: radius
    type        (enumerationStructureErrorCodeType), intent(  out), optional:: status
    
    decayingPotential=self%potentialNumerical(node, radius, status)
    return
  end function decayingPotential

  double precision function decayingCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    
    decayingCircularVelocity=self%circularVelocityNumerical(node, radius)
    return
  end function decayingCircularVelocity

  double precision function decayingCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    
    decayingCircularVelocityMaximum=self%circularVelocityMaximumNumerical(node)
    return
  end function decayingCircularVelocityMaximum

  double precision function decayingRadiusCircularVelocityMaximum(self,node)
    !!{
    Returns the radius (in Mpc) at which the maximum circular velocity is acheived in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    decayingRadiusCircularVelocityMaximum=self%radiusCircularVelocityMaximumNumerical(node)
    return
  end function decayingRadiusCircularVelocityMaximum

  ! Multiply by sqrt of factor
  double precision function decayingRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    
    decayingRadialVelocityDispersion=self%radialVelocityDispersionNumerical(node, radius)
    return
  end function decayingRadialVelocityDispersion

  double precision function decayingRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: specificAngularMomentum
    
    decayingRadiusFromSpecificAngularMomentum=self%radiusFromSpecificAngularMomentumNumerical(node, specificAngularMomentum)
    return
  end function decayingRadiusFromSpecificAngularMomentum

  double precision function decayingRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    
    decayingRotationNormalization=self%rotationNormalizationNumerical(node)
    return
  end function decayingRotationNormalization

  double precision function decayingEnergy(self,node)
    !!{
    Return the energy of a decaying halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    
    decayingEnergy=self%energyNumerical(node)
    return
  end function decayingEnergy

  double precision function decayingKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the decaying density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: waveNumber
    
    decayingKSpace=self%kSpaceNumerical(node, waveNumber)
    return
  end function decayingKSpace

  ! multiply by sqrt factor
  double precision function decayingFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the decaying density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: time
    
    decayingFreefallRadius=self%freefallRadiusNumerical(node, time)
    return
  end function decayingFreefallRadius

  double precision function decayingFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the rate of increase of the freefall radius in the decaying density profile at the specified {\normalfont
    \ttfamily time} (given in Gyr).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: time

    decayingFreefallRadiusIncreaseRate=self%freefallRadiusIncreaseRateNumerical(node, time)
    return
  end function decayingFreefallRadiusIncreaseRate
