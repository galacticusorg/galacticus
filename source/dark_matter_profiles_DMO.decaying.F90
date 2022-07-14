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

  use :: Dark_Matter_Profiles_Generic, only : enumerationNonAnalyticSolversEncode, enumerationNonAnalyticSolversIsValid, nonAnalyticSolversFallThrough

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
     integer                                              :: nonAnalyticSolver
     double precision                                     :: lifetime  , massSplitting
   contains
     final     ::                                      decayingDestructor
     procedure :: density                           => decayingDensity
     procedure :: densityLogSlope                   => decayingDensityLogSlope
     procedure :: radiusEnclosingDensity            => decayingRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => decayingRadiusEnclosingMass
     procedure :: radialMoment                      => decayingRadialMoment
     procedure :: enclosedMass                      => decayingEnclosedMass
     procedure :: potential                         => decayingPotential
     procedure :: circularVelocity                  => decayingCircularVelocity
     procedure :: circularVelocityMaximum           => decayingCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => decayingRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => decayingRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => decayingRotationNormalization
     procedure :: energy                            => decayingEnergy
     procedure :: kSpace                            => decayingKSpace
     procedure :: freefallRadius                    => decayingFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => decayingFreefallRadiusIncreaseRate
     procedure :: decayingFactor                    => decayingDecayingFactor
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
    type            (varying_string               )                :: nonAnalyticSolver

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"   name="darkMatterProfileDMO_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMODecaying(darkMatterParticle_, enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"  />
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="darkMatterParticle_"  />
    !!]
    return
  end function decayingConstructorParameters

  function decayingConstructorInternal(lifetime, massSplitting, nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily decaying} dark matter profile class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterProfileDMODecaying)                        :: self
    class           (darkMatterProfileDMOClass   ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass    ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterParticleClass      ), intent(in   ), target :: darkMatterParticle_
    double precision                              , intent(in   )         :: lifetime, massSplitting
    integer                                       , intent(in   )         :: nonAnalyticSolver
    !![
    <constructorAssign variables="*darkMatterParticle_, nonAnalyticSolver,*darkMatterProfileDMO_,*darkMatterHaloScale_"/>
    !!]
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime = darkMatterParticle_%lifetime()
       self%massSplitting = darkMatterParticle_%massSplitting()
    class default
       ! No decays.
       self%lifetime=-1.0d0
       self%massSplitting=0.0d0
    end select
    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    return
  end function decayingConstructorInternal

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

  subroutine decayingDecayingFactor(self, factor)
    !!{
    Return the change in mass factor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    double precision                              , intent(  out) :: factor
    class           (nodeComponentBasic          ), pointer       :: basic

    basic  => node%basic()
    factor = +1.0d0 - self%massSplitting * (+1.0d0 - exp(-basic%time() / self%lifetime))
    return
  end subroutine decayingDecayingFactor

  double precision function decayingDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileDMODecaying), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    double precision                                               :: factor
    
    call self%decayingDecayingFactor(factor)
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
    
    decayingDensityLogSlope=+self%darkMatterProfileDMO_%densityLogSlope(node,radius)
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
    double precision                                                      :: factor, oldDensity
    
    call self%decayingDecayingFactor(factor)
    oldDensity = density / factor
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingRadiusEnclosingDensity=self%darkMatterProfileDMO_%radiusEnclosingDensity         (node, oldDensity)
    else
       decayingRadiusEnclosingDensity=self                      %radiusEnclosingDensityNumerical(node, oldDensity)
    end if
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
    double precision                                                      :: factor
    class           (nodeComponentBasic          ), pointer               :: basic
    double precision                                                      :: factor, oldDensity
    
    call self%decayingDecayingFactor(factor)
    oldMass = mass / factor
    basic  => node%basic()
    factor = +1.0d0 - self%massSplitting * (+1.0d0 - exp(-basic%time() / self%lifetime))
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingRadiusEnclosingMass=self%darkMatterProfileDMO_%radiusEnclosingMass(node, oldMass)
    else
      decayingRadiusEnclosingMass=self                      %radiusEnclosingMassNumerical(node, oldMass)
    end if
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
    double precision                                                        :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
      decayingRadialMoment=self%darkMatterProfileDMO_%radialMoment         (node,moment,radiusMinimum,radiusMaximum) * factor
    else
      decayingRadialMoment=self                      %radialMomentNumerical(node,moment,radiusMinimum,radiusMaximum) * factor
    end if
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
    double precision                                              :: factor
    
    call self%decayingDecayingFactor(factor)
    decayingEnclosedMass=+self%darkMatterProfileDMO_%enclosedMass(node,radius) * factor
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
    integer                                       , intent(  out), optional :: status
    double precision                                                        :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingPotential=self%darkMatterProfileDMO_%potential         (node,radius,status) * factor
    else
       decayingPotential=self                      %potentialNumerical(node,radius,status) * factor
    end if
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
    double precision                                              :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
      decayingCircularVelocity=self%darkMatterProfileDMO_%circularVelocity         (node,radius) * sqrt(factor)
    else
      decayingCircularVelocity=self                      %circularVelocityNumerical(node,radius) * sqrt(factor)
    end if
    return
  end function decayingCircularVelocity

  double precision function decayingCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    double precision                                    :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
      decayingCircularVelocityMaximum=self%darkMatterProfileDMO_%circularVelocityMaximum         (node) * sqrt(factor)
    else
      decayingCircularVelocityMaximum=self                      %circularVelocityMaximumNumerical(node) * sqrt(factor)
    end if
    return
  end function decayingCircularVelocityMaximum

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
    double precision                                              :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingRadialVelocityDispersion=self%darkMatterProfileDMO_%radialVelocityDispersion(node,radius) * sqrt(factor)
    else
       decayingRadialVelocityDispersion=self%radialVelocityDispersionNumerical(node,radius) * sqrt(factor)
    end if
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
    double precision                                              :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingRadiusFromSpecificAngularMomentum=self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum(node,specificAngularMomentum / (sqrt(factor)))
    else
       decayingRadiusFromSpecificAngularMomentum=self                      %radiusFromSpecificAngularMomentumNumerical(node,specificAngularMomentum / (sqrt(factor)))
    end if
    return
  end function decayingRadiusFromSpecificAngularMomentum

  double precision function decayingRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    double precision                                   :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingRotationNormalization=self%darkMatterProfileDMO_%rotationNormalization(node)
    else
       decayingRotationNormalization=self                      %rotationNormalizationNumerical(node)
    end if
    return
  end function decayingRotationNormalization

  double precision function decayingEnergy(self,node)
    !!{
    Return the energy of a decaying halo density profile.
    !!}
    implicit none
    class(darkMatterProfileDMODecaying), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    double precision                                   :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingEnergy=self%darkMatterProfileDMO_%energy         (node) * factor**2
    else
       decayingEnergy=self                      %energyNumerical(node) * factor**2
    end if
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
    double precision                                                      :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingKSpace=self%darkMatterProfileDMO_%kSpace         (node,waveNumber) * factor
    else
       decayingKSpace=self                      %kSpaceNumerical(node,waveNumber) * factor
    end if
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
    double precision                                                      :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingFreefallRadius=self%darkMatterProfileDMO_%freefallRadius         (node,time * sqrt(factor))
    else
       decayingFreefallRadius=self                      %freefallRadiusNumerical(node,time * sqrt(factor))
    end if
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
    double precision                                                      :: factor
    
    call self%decayingDecayingFactor(factor)
    if (self%nonAnalyticSolver == nonAnalyticSolversFallThrough) then
       decayingFreefallRadiusIncreaseRate=self%darkMatterProfileDMO_%freefallRadiusIncreaseRate         (node,time * sqrt(factor))
    else
       decayingFreefallRadiusIncreaseRate=self                      %freefallRadiusIncreaseRateNumerical(node,time * sqrt(factor))
    end if
    return
  end function decayingFreefallRadiusIncreaseRate