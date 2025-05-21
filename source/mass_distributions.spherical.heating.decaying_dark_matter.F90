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
  Implements a decaying dark matter mass distribution heating class.
  !!}

  !![
  <massDistributionHeating name="massDistributionHeatingDecayingDarkMatter">
    <description>
      Implements heating from decays and response to mass loss. The mass loss heating is parameterized as:
      \begin{equation}
      \epsilon(r) = \gamma \frac{\mathrm{G}\Delta M(r)}{r},
      \end{equation}
      where $\gamma=${\normalfont \ttfamily [gamma]} sets the magnitude of the heating.
    </description>
  </massDistributionHeating>
  !!]
  type, extends(massDistributionHeatingClass) :: massDistributionHeatingDecayingDarkMatter
     !!{
     Implementation of a decaying dark matter mass distribution heating class.
     !!}
     private
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_    => null()
     logical                                            :: includeKickHeating
     double precision                                   :: lifetime                        , massSplitting                 , &
          &                                                gamma                           , velocityKick
     logical                                            :: energySpecificComputed          , energySpecificGradientComputed, &
          &                                                factorsComputed                 , potentialEscapeComputed       , &
          &                                                massEnclosedComputed
     double precision                                   :: radius                          , massEnclosed                  , &
          &                                                fractionRetained                , energyRetained                , &
          &                                                velocityDispersion              , velocityEscape                , &
          &                                                massLossEnergy                  , fractionDecayed               , &
          &                                                energySpecificGradient          , energySpecific                , &
          &                                                potentialEscape                 , radiusEscape                  , &
          &                                                time
   contains
     !![
     <methods>
       <method description="Compute memoized factors." method="computeFactors"/>
     </methods>
     !!]
     final     ::                                   decayingDarkMatterDestructor
     procedure :: specificEnergy                 => decayingDarkMatterSpecificEnergy
     procedure :: specificEnergyGradient         => decayingDarkMatterSpecificEnergyGradient
     procedure :: specificEnergyIsEveryWhereZero => decayingDarkMatterSpecificEnergyIsEverywhereZero
     procedure :: computeFactors                 => decayingDarkMatterComputeFactors
  end type massDistributionHeatingDecayingDarkMatter

  interface massDistributionHeatingDecayingDarkMatter
     !!{
     Constructors for the \refClass{massDistributionHeatingDecayingDarkMatter} mass distribution class.
     !!}
     module procedure decayingDarkMatterConstructorParameters
     module procedure decayingDarkMatterConstructorInternal
  end interface massDistributionHeatingDecayingDarkMatter

contains

  function decayingDarkMatterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingDecayingDarkMatter} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (massDistributionHeatingDecayingDarkMatter)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (darkMatterParticleClass                  ), pointer       :: darkMatterParticle_
    double precision                                                           :: radiusEscape       , time, &
         &                                                                        gamma
    logical                                                                    :: includeKickHeating

    !![
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
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <description>Parameter controlling the magnitude of heating due to mass loss.</description>
      <defaultValue>0.5d0</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>includeKickHeating</name>
      <source>parameters</source>
      <description>Parameter controlling whether heating due to velocity kicks is to be included.</description>
      <defaultValue>.true.</defaultValue>
    </inputParameter>
    <objectBuilder class="darkMatterParticle" name="darkMatterParticle_" source="parameters"/>
    !!]
    self=massDistributionHeatingDecayingDarkMatter(radiusEscape,time,gamma,includeKickHeating,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"/>
    !!]
    return
  end function decayingDarkMatterConstructorParameters
  
  function decayingDarkMatterConstructorInternal(radiusEscape,time,gamma,includeKickHeating,darkMatterParticle_) result(self)
    !!{
    Constructor for the \refClass{massDistributionHeatingDecayingDarkMatter} heating class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleDecayingDarkMatter
    use :: Error                , only : Error_Report
    implicit none
    type            (massDistributionHeatingDecayingDarkMatter)                        :: self
    class           (darkMatterParticleClass                  ), intent(in   ), target :: darkMatterParticle_
    double precision                                           , intent(in   )         :: time               , radiusEscape, &
         &                                                                                gamma
    logical                                                    , intent(in   )         :: includeKickHeating
    !![
    <constructorAssign variables="radiusEscape, time, gamma, includeKickHeating, *darkMatterParticle_"/>
    !!]
 
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime     = darkMatterParticle_%lifetime     ()
       self%massSplitting= darkMatterParticle_%massSplitting()
       self%velocityKick = darkMatterParticle_%velocityKick ()
    class default
       call Error_Report('expected a member of the `darkMatterParticleDecayingDarkMatter` class'//{introspection:location})
    end select
    ! Validate.
    if (gamma < 0.0d0) call Error_Report('`gamma` ≥ 0 is required'//{introspection:location})
    ! Initialize memoized calculations.
    self%radius                        =-huge(0.0d0)
    self%factorsComputed               =.false.
    self%massEnclosedComputed          =.false.
    self%potentialEscapeComputed       =.false.
    self%energySpecificComputed        =.false.
    self%energySpecificGradientComputed=.false.
    return
  end function decayingDarkMatterConstructorInternal

  subroutine decayingDarkMatterDestructor(self)
    !!{
    Destructor for the \refClass{massDistributionHeatingDecayingDarkMatter} dark matter profile heating class.
    !!}
    implicit none
    type(massDistributionHeatingDecayingDarkMatter), intent(inout) :: self
    
    !![
    <objectDestructor name="self%darkMatterParticle_"/>
    !!]
    return
  end subroutine decayingDarkMatterDestructor

  subroutine decayingDarkMatterComputeFactors(self,radius,massDistribution_,gradientRequired)
    !!{
    Compute various factors for specific energy of heating.
    !!}
    use :: Coordinates                     , only : coordinateSpherical               , assignment(=)
    use :: Decaying_Dark_Matter            , only : decayingDarkMatterFractionRetained, decayingDarkMatterEnergyRetained, decayingDarkMatterFractionRetainedDerivatives, decayingDarkMatterEnergyRetainedDerivatives
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionHeatingDecayingDarkMatter), intent(inout) :: self
    double precision                                           , intent(in   ) :: radius
    class           (massDistributionClass                    ), intent(inout) :: massDistribution_
    logical                                                    , intent(in   ) :: gradientRequired
    class           (kinematicsDistributionClass              ), pointer       :: kinematicsDistribution_
    double precision                                                           :: density                                  , densityLogGradient                , &
         &                                                                        fractionDerivativeVelocityEscapeScaleFree, fractionDerivativeVelocityKick    , &
         &                                                                        energyDerivativeVelocityEscapeScaleFree  , energyDerivativeVelocityKick      , &
         &                                                                        fractionDerivativeVelocityEscapeScaleFree, fractionDerivativeVelocityKick    , &
         &                                                                        energyDerivativeVelocityEscapeScaleFree  , energyDerivativeVelocityKick      , &
         &                                                                        fractionDerivativeVelocityDispersion     , energyDerivativeVelocityDispersion, &
         &                                                                        velocityDispersionGradient               , velocityEscapeGradient            , &
         &                                                                        potentialDifference                      , massLossGradient
    type            (coordinateSpherical                      )                :: coordinates
    
    if (radius /= self%radius) then
       self%radius                        =-huge(0.0d0)
       self%factorsComputed               =.false.
       self%massEnclosedComputed          =.false.
       self%potentialEscapeComputed       =.false.
       self%energySpecificComputed        =.false.
       self%energySpecificGradientComputed=.false.       
    end if
    if (.not.self%potentialEscapeComputed) then
       ! Compute the gravitational potential corresponding to escape (i.e. at some large radius).
       coordinates                 =[self%radiusEscape,0.0d0,0.0d0]
       self%potentialEscape        =massDistribution_%potential(coordinates)
       self%potentialEscapeComputed=.true.
    end if
    if (.not.self%massEnclosedComputed) then
       ! Enclosed mass is needed if mass loss is to be accounted for, or if computing gradients.
       if (self%gamma > 0.0d0 .or. gradientRequired) then
          self%massEnclosed        =massDistribution_%massEnclosedBySphere(radius)
          self%massEnclosedComputed=.true.
       end if
    end if
    if (.not.self%factorsComputed) then
       ! Find the fraction of particles that have decayed by this time.
       self%fractionDecayed =  +1.0d0              &
            &                  -exp(               &
            &                       -self%    time &
            &                       /self%lifetime &
            &                      )
       ! Compute the change in energy due to mass loss (assuming all decayed particles are lost).
       if (self%gamma > 0.0d0) then
          if (radius <= 0.0d0) then
             self%massLossEnergy=+0.0d0
          else
             self%massLossEnergy=+self%gamma                          &
                  &              *     gravitationalConstant_internal &
                  &              *self%massEnclosed                   &
                  &              /     radius
          end if
       else
          self%massLossEnergy=+0.0d0
       end if
       ! Compute the velocity dispersion.
       kinematicsDistribution_ => massDistribution_%kinematicsDistribution()
       coordinates             =  [radius,0.0d0,0.0d0]
       self%velocityDispersion =kinematicsDistribution_%velocityDispersion1D(coordinates,massDistribution_,massDistribution_)
       !![
       <objectDestructor name="kinematicsDistribution_"/>
       !!]
       ! Find the escape velocity.
       if (radius < self%radiusEscape) then
          coordinates        = [radius,0.0d0,0.0d0]
          potentialDifference=+self             %potentialEscape              &
               &              -massDistribution_%potential      (coordinates)
          if (potentialDifference > 0.0d0) then
             self%velocityEscape=+sqrt(                     &
                  &                    +2.0d0               &
                  &                    *potentialDifference &
                  &                   )
          else
             self%velocityEscape=+0.0d0
          end if
       else
          self%velocityEscape   =+0.0d0
       end if
       ! Store computation point.
       self%radius         =radius
       self%factorsComputed=.true.
    end if
    if (.not.self%energySpecificComputed) then
       ! Compute the fraction of particles, and kick energy retained.
       if (self%velocityDispersion > 0.0d0) then
          self%fractionRetained   =+decayingDarkMatterFractionRetained(self%velocityDispersion,self%velocityEscape,self%velocityKick)
          self%energyRetained     =+decayingDarkMatterEnergyRetained  (self%velocityDispersion,self%velocityEscape,self%velocityKick)
       else
          if (self%velocityKick > self%velocityEscape) then
             self%fractionRetained=+0.0d0
             self%energyRetained  =+0.0d0
          else
             self%fractionRetained=+1.0d0
             self%energyRetained  =+decayingDarkMatterEnergyRetained  (self%velocityDispersion,self%velocityEscape,self%velocityKick)
          end if
       end if
       ! If heating is not to be included, set the energy retained to zero.
       if (.not.self%includeKickHeating) self%energyRetained=0.0d0
       ! Compute the specific heating energy.
       self%energySpecific=+(                                                    &
            &                +self%energyRetained                                &
            &                +(                                                  &
            &                  +(1.0d0-self%fractionRetained)                    &
            &                  +       self%fractionRetained *self%massSplitting &
            &                 )                                                  &
            &                *self%massLossEnergy                                &
            &               )                                                    &
            &              *self%fractionDecayed
       self%energySpecificComputed=.true.
    end if
    if (.not.self%energySpecificGradientComputed.and.gradientRequired) then
       ! Compute properties of the dark matter density profile.
       coordinates       = [radius,0.0d0,0.0d0]
       density           =+massDistribution_%density              (coordinates                   )
       densityLogGradient=+massDistribution_%densityGradientRadial(coordinates,logarithmic=.true.)
       ! Compute the change in energy due to mass loss (assuming all decayed particles are lost).
       if (self%gamma > 0.0d0) then
          massLossGradient=+self%gamma                             &
               &           *gravitationalConstant_internal         &
               &           *(                                      &
               &             -         self%massEnclosed/radius    &
               &             +4.0d0*Pi*density          *radius**2 &
               &            )                                      &
               &           /radius
       else
          massLossGradient=+0.0d0
       end if
       ! Compute the fraction of particles, and kick energy retained derivatives with respect to the scale-free escape and kick velocities..
       call decayingDarkMatterFractionRetainedDerivatives(self%velocityDispersion,self%velocityEscape,self%velocityKick,fractionDerivativeVelocityEscapeScaleFree,fractionDerivativeVelocityKick)
       call decayingDarkMatterEnergyRetainedDerivatives  (self%velocityDispersion,self%velocityEscape,self%velocityKick,  energyDerivativeVelocityEscapeScaleFree,  energyDerivativeVelocityKick)
       ! Convert derivatives to dimensionful form.
       fractionDerivativeVelocityEscapeScaleFree=fractionDerivativeVelocityEscapeScaleFree/self%velocityDispersion
       fractionDerivativeVelocityKick           =fractionDerivativeVelocityKick           /self%velocityDispersion
       energyDerivativeVelocityEscapeScaleFree  =energyDerivativeVelocityEscapeScaleFree  /self%velocityDispersion
       energyDerivativeVelocityKick             =energyDerivativeVelocityKick             /self%velocityDispersion
       ! Construct the derivatives with respect to velocity dispersion.
       fractionDerivativeVelocityDispersion     =-      (     fractionDerivativeVelocityEscapeScaleFree*self%velocityEscape+fractionDerivativeVelocityKick*self%velocityKick)/self%velocityDispersion
       energyDerivativeVelocityDispersion       =-      (       energyDerivativeVelocityEscapeScaleFree*self%velocityEscape+  energyDerivativeVelocityKick*self%velocityKick)/self%velocityDispersion &
            &                                    +2.0d0*   self%energyRetained                                                                                               /self%velocityDispersion
       ! Compute the gradient in velocity dispersion.
       velocityDispersionGradient=+(                                        &
            &                       -     gravitationalConstant_internal    &
            &                       *self%massEnclosed                      &
            &                       /     radius                        **2 &
            &                       -self%velocityDispersion            **2 &
            &                       *     densityLogGradient                &
            &                       /     radius                            &
            &                      )                                        &
            &                     /2.0d0                                    &
            &                     /self%velocityDispersion
       ! Compute the gradient in escape velocity.
       if (self%velocityEscape > 0.0d0) then
          velocityEscapeGradient=-     gravitationalConstant_internal &
               &                 *self%massEnclosed                   &
               &                 /self%velocityEscape                 &
               &                 /     radius**2
       else
          velocityEscapeGradient=+0.0d0
       end if
       ! If heating is not to be included, set the energy retained derivatives to zero.
       if (.not.self%includeKickHeating) then
          energyDerivativeVelocityDispersion     =0.0d0
          energyDerivativeVelocityEscapeScaleFree=0.0d0
       end if
       ! Compute the specific heating gradient.
       self%energySpecificGradient=+(                                                                                                                                                                &
            &                        +((1.0d0-self%fractionRetained)+self%massSplitting*self%fractionRetained)*massLossGradient                                                                      &
            &                        +(energyDerivativeVelocityDispersion     +(-1.0d0+self%massSplitting)*fractionDerivativeVelocityDispersion     *self%massLossEnergy)*velocityDispersionGradient &
            &                        +(energyDerivativeVelocityEscapeScaleFree+(-1.0d0+self%massSplitting)*fractionDerivativeVelocityEscapeScaleFree*self%massLossEnergy)*velocityEscapeGradient     &
            &                       )                                                                                                                                                                &
            &                      *self%fractionDecayed
       self%energySpecificGradientComputed=.true.
    end if
    return
  end subroutine decayingDarkMatterComputeFactors

  double precision function decayingDarkMatterSpecificEnergy(self,radius,massDistribution_) result(energySpecific)
    !!{
    Compute the specific energy in a decaying dark matter-heated mass distribution.
    !!}
    implicit none
    class           (massDistributionHeatingDecayingDarkMatter), intent(inout) :: self
    double precision                                           , intent(in   ) :: radius
    class           (massDistributionClass                    ), intent(inout) :: massDistribution_
  
    call self%computeFactors(radius,massDistribution_,gradientRequired=.false.)
    energySpecific=self%energySpecific
    return
  end function decayingDarkMatterSpecificEnergy

  double precision function decayingDarkMatterSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating.
    !!}
    use :: Coordinates, only : coordinateSpherical, assignment(=)
    implicit none
    class           (massDistributionHeatingDecayingDarkMatter), intent(inout) :: self
    double precision                                           , intent(in   ) :: radius
    class           (massDistributionClass                    ), intent(inout) :: massDistribution_
    
    call self%computeFactors(radius,massDistribution_,gradientRequired=.true.)
    energySpecificGradient=self%energySpecificGradient
    return
  end function decayingDarkMatterSpecificEnergyGradient

  logical function decayingDarkMatterSpecificEnergyIsEverywhereZero(self) result(energySpecificIsEverywhereZero)
    !!{
    Returns true if the specific energy is everywhere zero.
    !!}
    implicit none
    class(massDistributionHeatingDecayingDarkMatter), intent(inout) :: self

    energySpecificIsEverywhereZero=self%lifetime <= 0.0d0 .or. self%massSplitting <= 0.0d0 .or. (.not.self%includeKickHeating .and. self%gamma == 0.0d0)
    return
  end function decayingDarkMatterSpecificEnergyIsEverywhereZero
