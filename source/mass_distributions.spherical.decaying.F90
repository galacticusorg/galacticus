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
  An implementation of decaying dark matter halo profiles.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <massDistribution name="massDistributionSphericalDecaying">
   <description>
     Decaying dark matter halo profiles.
    </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSphericalDecorator) :: massDistributionSphericalDecaying
     !!{
     A dark matter halo profile class implementing decaying dark matter halos.
     !!}
     private
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_     => null()
     double precision                                   :: lifetime_                        , massSplitting  , &
          &                                                velocityKick                     , potentialEscape, &
          &                                                radiusUndepleted                 , time           , &
          &                                                radiusEscape
     logical                                            :: potentialEscapeComputed
   contains
     !![
     <methods>
       <method description="Compute the density reduction factor." method="decayingFactor"/>
     </methods>
     !!]
     final     ::                         sphericalDecayingDestructor
     procedure :: decayingFactor       => sphericalDecayingDecayingFactor
     procedure :: density              => sphericalDecayingDensity
     procedure :: massEnclosedBySphere => sphericalDecayingMassEnclosedBySphere
     procedure :: radiusEnclosingMass  => sphericalDecayingRadiusEnclosingMass
     procedure :: useUndecorated       => sphericalDecayingUseUndecorated
  end type massDistributionSphericalDecaying

  interface massDistributionSphericalDecaying
     !!{
     Constructors for the \refClass{massDistributionSphericalDecaying} mass distribution class.
     !!}
     module procedure sphericalDecayingConstructorParameters
     module procedure sphericalDecayingConstructorInternal
  end interface massDistributionSphericalDecaying

contains

  function sphericalDecayingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalDecaying} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionSphericalDecaying)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (massDistributionClass            ), pointer       :: massDistribution_
    class           (darkMatterParticleClass          ), pointer       :: darkMatterParticle_
    type            (varying_string                   )                :: massType                              , componentType
    double precision                                                   :: toleranceRelativePotential            , radiusEscape                       , &
         &                                                                time
    logical                                                            :: tolerateVelocityMaximumFailure        , toleratePotentialIntegrationFailure, &
         &                                                                tolerateEnclosedMassIntegrationFailure

    !![
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
    <inputParameter>
      <name>toleranceRelativePotential</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the gravitational potential.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerateEnclosedMassIntegrationFailure</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, tolerate failures to find the mass enclosed as a function of radius.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerateVelocityMaximumFailure</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, tolerate failures to find the radius of the maximum circular velocity.</description>
    </inputParameter>
    <inputParameter>
      <name>toleratePotentialIntegrationFailure</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, tolerate failures to compute the potential.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusEscape</name>
      <source>parameters</source>
      <description>The radius beyond which a particle is assumed to have escaped the potential.</description>
    </inputParameter>
    <inputParameter>
      <name>time</name>
      <source>parameters</source>
      <description>The time at which decays should be evaluated.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle" name="darkMatterParticle_" source="parameters"/>
    <objectBuilder class="massDistribution"   name="massDistribution_"   source="parameters"/>
      !!]
    select type (massDistribution_)
    class is (massDistributionSpherical)
       self=massDistributionSphericalDecaying(toleranceRelativePotential,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,tolerateEnclosedMassIntegrationFailure,radiusEscape,time,darkMatterParticle_,massDistribution_,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    class default
       call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
    end select
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="massDistribution_"  />
    <objectDestructor name="darkMatterParticle_"/>
    !!]
    return
  end function sphericalDecayingConstructorParameters
  
  function sphericalDecayingConstructorInternal(toleranceRelativePotential,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,tolerateEnclosedMassIntegrationFailure,radiusEscape,time,darkMatterParticle_,massDistribution_,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionSphericalDecaying} mass distribution class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleDecayingDarkMatter
    implicit none
    type            (massDistributionSphericalDecaying)                          :: self
    class           (massDistributionSpherical        ), intent(in   ), target   :: massDistribution_
    class           (darkMatterParticleClass          ), intent(in   ), target   :: darkMatterParticle_
    double precision                                   , intent(in   )           :: toleranceRelativePotential            , radiusEscape                       , &
         &                                                                          time
    logical                                            , intent(in   )           :: tolerateVelocityMaximumFailure        , toleratePotentialIntegrationFailure, &
         &                                                                          tolerateEnclosedMassIntegrationFailure
    type            (enumerationComponentTypeType     ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType          ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="toleranceRelativePotential, tolerateVelocityMaximumFailure, toleratePotentialIntegrationFailure, tolerateEnclosedMassIntegrationFailure, radiusEscape, time, *darkMatterParticle_, *massDistribution_"/>
    !!]
 
    ! In models with decays, tolerate failures in integration of the density profile (as this can become almost fully disrupted).
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime_                             =+darkMatterParticle_%lifetime     ()
       self%massSplitting                         =+darkMatterParticle_%massSplitting()
       self%velocityKick                          =+darkMatterParticle_%velocityKick ()
       self%tolerateEnclosedMassIntegrationFailure=.true.
    class default
       call Error_Report('expected a member of the `darkMatterParticleDecayingDarkMatter` class'//{introspection:location})
    end select
    ! Initialize.
    self%potentialEscapeComputed=.false.
    self%dimensionless          =.false.
    self%radiusUndepleted       =+0.0d0
    return
  end function sphericalDecayingConstructorInternal

  subroutine sphericalDecayingDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionSphericalDecaying} mass distribution class.
    !!}
    implicit none
    type(massDistributionSphericalDecaying), intent(inout) :: self

    !![
    <objectDestructor name="self%massDistribution_"  />
    <objectDestructor name="self%darkMatterParticle_"/>
    !!]
    return
  end subroutine sphericalDecayingDestructor

  logical function sphericalDecayingUseUndecorated(self) result(useUndecorated)
    !!{
    Determines whether to use the undecorated solution.
    !!}
    implicit none
    class(massDistributionSphericalDecaying), intent(inout) :: self

    useUndecorated=.false.
    return
  end function sphericalDecayingUseUndecorated
  
  double precision function sphericalDecayingDecayingFactor(self,radius) result(factor)
    !!{
    Return the remaining mass fraction in the profile.
    !!}
    use :: Coordinates         , only : coordinateSpherical               , assignment(=)
    use :: Decaying_Dark_Matter, only : decayingDarkMatterFractionRetained
    use :: Error               , only : Error_Report
    implicit none
    class           (massDistributionSphericalDecaying), intent(inout) :: self
    double precision                                   , intent(in   ) :: radius
    double precision                                   , parameter     :: depletionNegligble     =1.0d-2
    class           (kinematicsDistributionClass      ), pointer       :: kinematicsDistribution_
    double precision                                                   :: velocityDispersion            , velocityEscape , &
         &                                                                fractionRetained              , fractionDecayed, &
         &                                                                potentialDifference
    type            (coordinateSpherical              )                :: coordinates

    ! For radii within the undepleted region of the halo, the depletion factor is, by definition, 1.
    if (radius < self%radiusUndepleted) then
       factor=1.0d0
       return
    end if
    ! Outside of the undepleted region, compute the depletion factor directly.
    !! Find the escape velocity.
    if (radius < self%radiusEscape) then
       if (.not.self%potentialEscapeComputed) then
          coordinates                 =[self%radiusEscape,0.0d0,0.0d0]
          self%potentialEscape        =self%massDistribution_%potential(coordinates)
          self%potentialEscapeComputed=.true.
       end if
       coordinates        =[radius,0.0d0,0.0d0]
       potentialDifference=+self                  %potentialEscape              &
            &              -self%massDistribution_%potential      (coordinates)
       if (potentialDifference > 0.0d0) then
          velocityEscape=+sqrt(                     &
               &               +2.0d0               &
               &               *potentialDifference &
               &              )
       else
          velocityEscape=+0.0d0
       end if
    else
       velocityEscape=+0.0d0
    end if
    kinematicsDistribution_ => self%massDistribution_%kinematicsDistribution()
    coordinates             =  [radius,0.0d0,0.0d0]
    velocityDispersion      =  kinematicsDistribution_%velocityDispersion1D(coordinates,self%massDistribution_,self%massDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    if (velocityDispersion > 0.0d0) then
       fractionRetained   =+decayingDarkMatterFractionRetained(velocityDispersion,velocityEscape,self%velocityKick)
    else
       if (self%velocityKick > velocityEscape) then
          fractionRetained=+0.0d0
       else
          fractionRetained=+1.0d0
       end if
    end if
    fractionDecayed =  +1.0d0               &
         &             -exp(                &
         &                  -self%    time  &
         &                  /self%lifetime_ &
         &                 )
    factor          =  +(+1.0d0-     fractionDecayed ) & ! { Fraction of particles undecayed - have 100% of their original mass.
         &             +             fractionDecayed   & ! ⎧ Fraction of particles decayed...
         &             *             fractionRetained  & ! ⎨  ...but not escaped...
         &             *(+1.0d0-self%massSplitting   )   ! ⎩  ...have 1-ε of their original mass.
    ! Check for negligible depletion, and update the undepleted radius to the largest such radius yet found.
    if (factor > 1.0d0-depletionNegligble .and. radius > self%radiusUndepleted) self%radiusUndepleted=radius
    return
  end function sphericalDecayingDecayingFactor

  double precision function sphericalDecayingDensity(self,coordinates) result(density)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a scaled spherical mass distribution.
    !!}
    implicit none
    class(massDistributionSphericalDecaying), intent(inout) :: self
    class(coordinate                       ), intent(in   ) :: coordinates

    density=+self%massDistribution_%density       (coordinates             ) &
         &  *self                  %decayingFactor(coordinates%rSpherical())
    return
  end function sphericalDecayingDensity

  double precision function sphericalDecayingMassEnclosedBySphere(self,radius) result(mass)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a decaying mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecaying), intent(inout), target :: self
    double precision                                   , intent(in   )         :: radius

    if (radius <= 0.0d0) then
       mass=0.0d0
    else if (radius <= self%radiusUndepleted) then
       ! Within the undepleted radius the mass is unchanged.
       mass=+self%massDistribution_%massEnclosedBySphere         (radius)
    else
       mass=+self                  %massEnclosedBySphereNumerical(radius)
    end if
    return
  end function sphericalDecayingMassEnclosedBySphere
  
  double precision function sphericalDecayingRadiusEnclosingMass(self,mass,massFractional) result(radius)
    !!{
    Computes the radius enclosing a given mass or mass fraction for decaying spherical mass distributions.
    !!}
    implicit none
    class           (massDistributionSphericalDecaying), intent(inout), target   :: self
    double precision                                   , intent(in   ), optional :: mass, massFractional

    radius=self%radiusEnclosingMassNumerical(mass,massFractional)
    return
  end function sphericalDecayingRadiusEnclosingMass
