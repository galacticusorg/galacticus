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
  A dark matter halo profile heating class which accounts for heating from decays.
  !!}

  use :: Kind_Numbers         , only : kind_int8
  use :: Dark_Matter_Particles, only : darkMatterParticleClass
  use :: Decaying_Dark_Matter , only : decayingDarkMatterFractionRetained, decayingDarkMatterEnergyRetained, decayingDarkMatterFractionRetainedDerivatives, decayingDarkMatterEnergyRetainedDerivatives

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingDDMv4">
   <description>
    Implements heating from decays and response to mass loss.
   </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingDDMv4
     !!{
     A dark matter profile heating class which accounts for heating due to decays.
     !!}
     private
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_  => null()
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     logical                                             :: heating_                      , massLoss_
     double precision                                    :: lifetime_                     , massSplitting_,&
          &                                                 gamma_                        , velocityKick
  contains
     procedure :: specificEnergy                 => DDMv4SpecificEnergy
     procedure :: specificEnergyGradient         => DDMv4SpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero => DDMv4SpecificEnergyIsEverywhereZero
  end type darkMatterProfileHeatingDDMv4

  interface darkMatterProfileHeatingDDMv4
     !!{
     Constructors for the {\normalfont \ttfamily DDMv4} dark matter profile heating class.
     !!}
     module procedure DDMv4ConstructorParameters
     module procedure DDMv4ConstructorInternal
  end interface darkMatterProfileHeatingDDMv4

contains

  function DDMv4ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily DDM} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (darkMatterProfileHeatingDDMv4), target        :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(darkMatterParticleClass      ), pointer       :: darkMatterParticle_
    class(darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
         
    !![
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileHeatingDDMv4(darkMatterParticle_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function DDMv4ConstructorParameters

  function DDMv4ConstructorInternal(darkMatterParticle_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily DDM} dark matter profile heating scales class.
    !!}
    use :: Dark_Matter_Particles       , only : darkMatterParticleDecayingDarkMatter
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    type (darkMatterProfileHeatingDDMv4)                        :: self
    class(darkMatterParticleClass      ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterHaloScaleClass     ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterParticle_,*darkMatterHaloScale_"/>
    !!]
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime_     = darkMatterParticle_%lifetime     ()
       self%massSplitting_= darkMatterParticle_%massSplitting()
       self%heating_      = darkMatterParticle_%heating      ()
       self%massLoss_     = darkMatterParticle_%massLoss     ()
       self%gamma_        = darkMatterParticle_%gamma        ()
       self%velocityKick  =+self               %massSplitting_  &
            &              *speedLight                          &
            &              /kilo
    class default
       ! No decays.
       self%lifetime_     =-1.0d0
       self%massSplitting_=+0.0d0
       self%heating_      =.false.
       self%massLoss_     =.false.
       self%gamma_        =+1.0d0
       self%velocityKick  =+0.0d0
    end select
    return
  end function DDMv4ConstructorInternal

  double precision function DDMv4SpecificEnergy(self,node,radius,darkMatterProfileDMO_) result(energySpecific)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class           (darkMatterProfileHeatingDDMv4), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic           ), pointer       :: basic
    double precision                               , parameter     :: fractionRadiusVirialMaximum=1.0d3
    double precision                                               :: massLossEnergy                   , velocityDispersion , &
         &                                                            velocityEscape                   , potentialDifference, &
         &                                                            fractionRetained                 , energyRetained     , &
         &                                                            fractionDecayed

    ! Find the fraction of particles that have decayed by this time.
    basic           =>  node%basic()
    fractionDecayed =  +1.0d0                  &
         &             -exp(                   &
         &                  -basic%time     () &
         &                  /self %lifetime_   &
         &                 )
    ! Compute the change in energy due to mass loss (assuming all decayed particles are lost).
    if (self%massLoss_) then
       massLossEnergy=+self%gamma_                                     &
            &         *gravitationalConstantGalacticus                 &
            &         *darkMatterProfileDMO_%enclosedMass(node,radius) &
            &         /radius
    else
       massLossEnergy=+0.0d0
    end if
    ! Compute the velocity dispersion.
    velocityDispersion=darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    ! Find the escape velocity.
    if (radius < fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) then
       potentialDifference=+darkMatterProfileDMO_%potential(node,fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) &
            &              -darkMatterProfileDMO_%potential(node,                                                      radius            )
       if (potentialDifference > 0.0d0) then
          velocityEscape=+sqrt(                     &
               &               +2.0d0               &
               &               *potentialDifference &
               &              )
       else
          velocityEscape=+0.0d0
       end if
    else
       velocityEscape   =+0.0d0
    end if
    ! Compute the fraction of particles, and kick energy retained.
    if (velocityDispersion > 0.0d0) then
       fractionRetained   =+decayingDarkMatterFractionRetained(velocityDispersion,velocityEscape,self%velocityKick)
       energyRetained     =+decayingDarkMatterEnergyRetained  (velocityDispersion,velocityEscape,self%velocityKick)
    else
       if (self%velocityKick > velocityEscape) then
          fractionRetained=+0.0d0
          energyRetained  =+0.0d0
       else
          fractionRetained=+1.0d0
          energyRetained  =+decayingDarkMatterEnergyRetained  (velocityDispersion,velocityEscape,self%velocityKick)
       end if
    end if
    ! If heating is not to be included, set the energy retained to zero.
    if (.not.self%heating_) energyRetained=0.0d0
    ! Compute the specific heating energy.
    energySpecific=+(                                                &
         &           +energyRetained                                 &
         &           +(                                              &
         &             +(1.0d0-fractionRetained)                     &
         &             +       fractionRetained *self%massSplitting_ &
         &            )                                              &
         &           *massLossEnergy                                 &
         &          )                                                &
         &         *fractionDecayed
    return
  end function DDMv4SpecificEnergy

  double precision function DDMv4SpecificEnergyGradient(self,node,radius,darkMatterProfileDMO_) result(energySpecificGradient)
    !!{
    Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Galacticus_Nodes                , only : nodeComponentBasic
    implicit none
    class           (darkMatterProfileHeatingDDMv4), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic           ), pointer       :: basic
    double precision                               , parameter     :: fractionRadiusVirialMaximum=1.0d3
    double precision                                               :: massLossEnergy, massLossGradient,&
   & velocityDispersion,velocityDispersionGradient,&
   & velocityEscape,  velocityEscapeGradient,&
   & fractionRetained, energyRetained,&
   & fractionDerivativeVelocityEscapeScaleFree, energyDerivativeVelocityEscapeScaleFree,&
   & fractionDerivativeVelocityKick,energyDerivativeVelocityKick,&
   & fractionDerivativeVelocityDispersion,energyDerivativeVelocityDispersion,&
   & massEnclosed, potentialDifference,&
   & density, densityLogGradient,&
   & fractionDecayed
    
    ! Find the fraction of particles that have decayed by this time.
    basic           =>  node%basic()
    fractionDecayed =  +1.0d0                  &
         &             -exp(                   &
         &                  -basic%time     () &
         &                  /self %lifetime_   &
         &                 )
    ! Compute properties of the dark matter density profile.
    massEnclosed      =+darkMatterProfileDMO_%enclosedMass   (node,radius)
    density           =+darkMatterProfileDMO_%density        (node,radius)
    densityLogGradient=+darkMatterProfileDMO_%densityLogSlope(node,radius)
    ! Compute the change in energy due to mass loss (assuming all decayed particles are lost).
    if (self%massLoss_) then
       massLossGradient=self%gamma_                        &
            &           *gravitationalConstantGalacticus   &
            &           *(                                 &
            &             -         massEnclosed/radius    &
            &             +4.0d0*Pi*density     *radius**2 &
            &            )                                 &
            &           /radius
       massLossEnergy  =+self%gamma_                       &
            &           *gravitationalConstantGalacticus   &
            &           *massEnclosed                      &
            &           /radius
    else
       massLossEnergy  =+0.0d0
       massLossGradient=+0.0d0
    end if
    ! Compute the velocity dispersion.
    velocityDispersion=darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    ! Find the escape velocity.
    if (radius < fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) then
       potentialDifference=+darkMatterProfileDMO_%potential(node,fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) &
            &              -darkMatterProfileDMO_%potential(node,                                                      radius            )
       if (potentialDifference > 0.0d0) then
          velocityEscape=+sqrt(                     &
               &               +2.0d0               &
               &               *potentialDifference &
               &              )
       else
          velocityEscape=+0.0d0
       end if
    else
       velocityEscape   =+0.0d0
    end if
    ! Compute the fraction of particles, and kick energy retained, along with their derivatives with respect to the scale-free escape and kick velocities..
    fractionRetained=decayingDarkMatterFractionRetained(velocityDispersion,velocityEscape,self%velocityKick)
    energyRetained  =decayingDarkMatterEnergyRetained  (velocityDispersion,velocityEscape,self%velocityKick)
    call decayingDarkMatterFractionRetainedDerivatives(velocityDispersion,velocityEscape,self%velocityKick,fractionDerivativeVelocityEscapeScaleFree,fractionDerivativeVelocityKick)
    call decayingDarkMatterEnergyRetainedDerivatives  (velocityDispersion,velocityEscape,self%velocityKick,  energyDerivativeVelocityEscapeScaleFree,  energyDerivativeVelocityKick)
    ! Convert derivatives to dimensionful form.
    fractionDerivativeVelocityEscapeScaleFree=fractionDerivativeVelocityEscapeScaleFree/velocityDispersion
    fractionDerivativeVelocityKick           =fractionDerivativeVelocityKick           /velocityDispersion
    energyDerivativeVelocityEscapeScaleFree  =energyDerivativeVelocityEscapeScaleFree  /velocityDispersion
    energyDerivativeVelocityKick             =energyDerivativeVelocityKick             /velocityDispersion
    ! Construct the derivatives with respect to velocity dispersion.
    fractionDerivativeVelocityDispersion     =-      (fractionDerivativeVelocityEscapeScaleFree*velocityEscape+fractionDerivativeVelocityKick*self%velocityKick)/velocityDispersion
    energyDerivativeVelocityDispersion       =-      (  energyDerivativeVelocityEscapeScaleFree*velocityEscape+  energyDerivativeVelocityKick*self%velocityKick)/velocityDispersion &
         &                                    +2.0d0*   energyRetained                                                                                          /velocityDispersion
    ! Compute the gradient in velocity disperson.
    velocityDispersionGradient=+(                                    &
         &                       -gravitationalConstantGalacticus    &
         &                       *massEnclosed                       &
         &                       /radius                         **2 &
         &                       -velocityDispersion             **2 &
         &                       *densityLogGradient                 &
         &                       /radius                             &
         &                      )                                    &
         &                     /2.0d0                                &
         &                     /velocityDispersion
    ! Compute the gradient in escape velocity.
    if (velocityEscape > 0.0d0) then
       velocityEscapeGradient=-gravitationalConstantGalacticus &
            &                 *massEnclosed                    &
            &                 /velocityEscape                  &
            &                 *radius**2
    else
       velocityEscapeGradient=+0.0d0
    end if
    ! If heating is not to be included, set the energy retained derivatives to zero.
    if (.not.self%heating_) then
       energyDerivativeVelocityDispersion     =0.0d0
       energyDerivativeVelocityEscapeScaleFree=0.0d0
    end if
    ! Compute the specific heating gradient.
    energySpecificGradient=+(                                                                                                                                                            &
         &                   +((1.0d0-fractionRetained)+self%massSplitting_*fractionRetained)*massLossGradient                                                                           &
         &                   +(energyDerivativeVelocityDispersion     +(-1.0d0+self%massSplitting_)*fractionDerivativeVelocityDispersion     *massLossEnergy)*velocityDispersionGradient &
         &                   +(energyDerivativeVelocityEscapeScaleFree+(-1.0d0+self%massSplitting_)*fractionDerivativeVelocityEscapeScaleFree*massLossEnergy)*velocityEscapeGradient     &
         &                  )                                                                                                                                                            &
         &                 *fractionDecayed
    return
  end function DDMv4SpecificEnergyGradient

  logical function DDMv4SpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !!{
    Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileHeatingDDMv4), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    !$GLC attributes unused :: darkMatterProfileDMO_

    DDMv4SpecificEnergyIsEverywhereZero=(self%lifetime_ <= 0.0d0) .or. (self%massSplitting_ <= 0.0d0) .or. (.not. (self%heating_ .or. self%massLoss_))       
    return
  end function DDMv4SpecificEnergyIsEverywhereZero
