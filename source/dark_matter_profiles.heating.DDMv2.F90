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

  use :: Kind_Numbers            , only : kind_int8
  use :: Dark_Matter_Particles   , only : darkMatterParticleClass

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingDDMv2">
   <description>
    Implements heating from decays and response to mass loss.
   </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingDDMv2
     !!{
     A dark matter profile heating class which accounts for heating due to decays.
     !!}
     private
     class           (darkMatterParticleClass), pointer     :: darkMatterParticle_     => null()
     class           (darkMatterHaloScaleClass     ), pointer             :: darkMatterHaloScale_ => null()
     logical                                                :: heating_, massLoss_
     double precision                                       :: lifetime_, massSplitting_, gamma_
  contains
     procedure :: specificEnergy                  => DDMv2SpecificEnergy
     procedure :: specificEnergyGradient          => DDMv2SpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero  => DDMv2SpecificEnergyIsEverywhereZero
  end type darkMatterProfileHeatingDDMv2

  interface darkMatterProfileHeatingDDMv2
     !!{
     Constructors for the {\normalfont \ttfamily DDMv2} dark matter profile heating class.
     !!}
     module procedure DDMv2ConstructorParameters
     module procedure DDMv2ConstructorInternal
  end interface darkMatterProfileHeatingDDMv2

contains

  function DDMv2ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily DDM} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileHeatingDDMv2), target              :: self
    type            (inputParameters              ), intent(inout)       :: parameters
    class           (darkMatterParticleClass      ), pointer             :: darkMatterParticle_
    class           (darkMatterHaloScaleClass     ), pointer             :: darkMatterHaloScale_
         
    !![
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=darkMatterProfileHeatingDDMv2(darkMatterParticle_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"  />
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function DDMv2ConstructorParameters

  function DDMv2ConstructorInternal(darkMatterParticle_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily DDM} dark matter profile heating scales class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleDecayingDarkMatter
    implicit none
    type            (darkMatterProfileHeatingDDMv2)                        :: self
    class           (darkMatterParticleClass      ), intent(in   ), target :: darkMatterParticle_
    class           (darkMatterHaloScaleClass    ), intent(in   ), target  :: darkMatterHaloScale_
    !![
    <constructorAssign variables="*darkMatterParticle_,*darkMatterHaloScale_"/>
    !!]
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime_ = darkMatterParticle_%lifetime()
       self%massSplitting_ = darkMatterParticle_%massSplitting()
       self%heating_ = darkMatterParticle_%heating()
       self%massLoss_ = darkMatterParticle_%massLoss()
       self%gamma_ = darkMatterParticle_%gamma()
    class default
       ! No decays.
       self%lifetime_=-1.0d0
       self%massSplitting_=0.0d0
       self%heating_ = .false.
       self%massLoss_ = .false.
       self%gamma_ =+1.0d0
    end select
    return
  end function DDMv2ConstructorInternal

  double precision function DDMv2SpecificEnergy(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Statistics_Distributions    , only : distributionFunction1DNonCentralChiDegree3
    implicit none
    class           (darkMatterProfileHeatingDDMv2), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic           ), pointer       :: basic
    double precision                               , parameter     :: fractionRadiusVirialMaximum=1.0d3
    type            (distributionFunction1DNonCentralChiDegree3)   :: distributionSpeed
    double precision                                               :: heatingEnergy, massLossEnergy, decayingEscapeFraction, velocityDispersion, velocityEscape, potentialDifference

    basic             => node%basic()
    heatingEnergy = 0.0d0
    massLossEnergy = 0.0d0
    if (self%heating_) heatingEnergy = +0.5d0*(                                                     &
                        &              +1.0d0 - exp(-basic%time() / self%lifetime_)                 &
                        &                     )                                                     &
                        &              *self%massSplitting_**2                                      &
                        &              *(speedLight/kilo)**2
    if (self%massLoss_) massLossEnergy =  self%gamma_*self%massSplitting_                          &
                            &             *darkMatterProfileDMO_%enclosedMass(node, radius)        &
                            &             *gravitationalConstantGalacticus/radius                  &
                            &             *(+1.0d0 - exp(-basic%time() / self%lifetime_))
if (heatingEnergy > 0.0d0) then
    velocityDispersion=darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
    if (radius < fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) then
       potentialDifference=+darkMatterProfileDMO_%potential(node,fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) &
            &              -darkMatterProfileDMO_%potential(node,                                                      radius            )
       if (potentialDifference > 0.0d0) then
          velocityEscape=+sqrt(                                                                                                                  &
               &               +2.0d0                                                                                                            &
               &               *potentialDifference                                                                                                                &
               &              )
       else
          velocityEscape=0.0d0
       end if
    else
       velocityEscape=+0.0d0
    end if
    if (velocityDispersion > 0.0d0) then
       distributionSpeed     =distributionFunction1DNonCentralChiDegree3(                                                                   &
            &                                                            +(                                                                &
            &              +self%massSplitting_                                                                                            &
            &              *speedLight                                                                                                     &
            &              /kilo                                                           &
            &              /velocityDispersion                                                                                             &
            &                                                             )**2                                                             &
            &                                                           )
       decayingEscapeFraction=+1.0d0                           &
            &                 -distributionSpeed%cumulative(  &
            &                                               ( &
            &                                 +velocityEscape &
            &                             /velocityDispersion &
            &                                            )**2 &
            &                                              )
    else
       if (self%massSplitting_*speedLight/kilo > velocityEscape) then
          decayingEscapeFraction = +1.0d0
       else
          decayingEscapeFraction = +0.0d0
       end if
    end if
 else
    decayingEscapeFraction=0.0d0
    end if
    DDMv2SpecificEnergy = (1.0d0-decayingEscapeFraction)*heatingEnergy + (decayingEscapeFraction+self%massSplitting_*(1.0d0-decayingEscapeFraction))*massLossEnergy
    return
  end function DDMv2SpecificEnergy

  double precision function DDMv2SpecificEnergyGradient(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Statistics_Distributions    , only : distributionFunction1DNonCentralChiDegree3
    implicit none
    class           (darkMatterProfileHeatingDDMv2), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic           ), pointer       :: basic
    double precision                               , parameter     :: fractionRadiusVirialMaximum=1.0d3
    double precision                               , parameter     :: radiusEpsilon=1.0d-6
    type            (distributionFunction1DNonCentralChiDegree3)   :: distributionSpeed
    double precision                                               :: heatGradient, massLossGradient, decayingEscapeFractionGradient, decayingEscapeFraction, decayingEscapeFractionUpper, velocityDispersion, velocityDispersionUpper, velocityEscape, velocityEscapeUpper, radiusUpper, massLossEnergy, heatingEnergy, potentialDifference

    basic          => node%basic()
    massLossEnergy = +0.0d0
    heatingEnergy  = +0.0d0
    if (self%massLoss_) then
       massLossGradient=self%gamma_*self%massSplitting_                                                              &
            &        *(darkMatterProfileDMO_%enclosedMass(node, radius)*gravitationalConstantGalacticus/radius&
            &        +gravitationalConstantGalacticus*4.0d0*Pi*darkMatterProfileDMO_%density(node, radius)*radius**2)/radius &
            &        *(+1.0d0 - exp(-basic%time() / self%lifetime_))
       massLossEnergy =  self%gamma_*self%massSplitting_ &
            & *darkMatterProfileDMO_%enclosedMass(node, radius)        &
            & *gravitationalConstantGalacticus/radius                  &
            &             *(+1.0d0 - exp(-basic%time() /self%lifetime_))
    else
       massLossGradient=+0.0d0
    end if
    if (self%heating_) then
       heatingEnergy = +0.5d0*( &
            &              +1.0d0 - exp(-basic%time() / self%lifetime_)                 &
            &                     ) &
            &              *self%massSplitting_**2 &
            &              *(speedLight/kilo)**2
    end if
    if (heatingEnergy > 0.0d0) then
       velocityDispersion=darkMatterProfileDMO_%radialVelocityDispersion(node,radius)
       if (velocityDispersion > 0.0d0) then
          if (radius < fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) then
             potentialDifference=+darkMatterProfileDMO_%potential(node,fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) &
                  &              -darkMatterProfileDMO_%potential(node,                                                      radius            )
             if (potentialDifference > 0.0d0) then
                velocityEscape=+sqrt(                 &
                     &               +2.0d0            &
                     &               *potentialDifference                &
                     &              )
             else
                velocityEscape=+0.0d0
             end if
          else
             velocityEscape=+0.0d0
          end if
          distributionSpeed     =distributionFunction1DNonCentralChiDegree3( &
               &                                                            +( &
               &              +self%massSplitting_                             &
               &              *speedLight                                      &
               &                                                              /kilo &
               &              /velocityDispersion                                   &
               &                                                             )**2   &
               &                                                           )
          decayingEscapeFraction=+1.0d0                           &
               &                 -distributionSpeed%cumulative(  &
               &                                               ( &
               &                                 +velocityEscape &
               &                             /velocityDispersion &
               &                                            )**2 &
               &                                              )
       else
          decayingEscapeFraction = +1.0d0
       end if
       radiusUpper = radius+radiusEpsilon
       velocityDispersionUpper=darkMatterProfileDMO_%radialVelocityDispersion(node,radiusUpper)
       if (velocityDispersionUpper > 0.0d0) then
          if (radiusUpper < fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) then
             potentialDifference=+darkMatterProfileDMO_%potential(node,fractionRadiusVirialMaximum*self%darkMatterHaloScale_%radiusVirial(node)) &
                  &              -darkMatterProfileDMO_%potential(node,                                                      radiusUpper       )
             if (potentialDifference > 0.0d0) then
                velocityEscapeUpper=+sqrt(                 &
                  &               +2.0d0            &
                  &               *potentialDifference                &
                  &              )
             else
                velocityEscapeUpper=+0.0d0
             end if
          else
             velocityEscapeUpper=+0.0d0
          end if
          distributionSpeed     =distributionFunction1DNonCentralChiDegree3( &
               &                                                            +( &
               &              +self%massSplitting_                             &
               &              *speedLight                                      &
               &                                                              /kilo &
               &              /velocityDispersionUpper                                   &
               &                                                             )**2   &
               &                                                           )
          decayingEscapeFractionUpper=1.0d0                           &
               &                 -distributionSpeed%cumulative(  &
               &                                               ( &
               &                                 +velocityEscapeUpper &
               &                             /velocityDispersionUpper &
               &                                            )**2 &
               &                                              )
       else
          decayingEscapeFractionUpper=+1.0d0
       end if
       decayingEscapeFractionGradient = (decayingEscapeFractionUpper-decayingEscapeFractionUpper-decayingEscapeFraction)/radiusEpsilon
    else
       decayingEscapeFraction=0.0d0
       decayingEscapeFractionGradient=0.0d0
    end if
    heatGradient=+0.0d0
    DDMv2SpecificEnergyGradient = (decayingEscapeFraction-self%massSplitting_*(1.0d0-decayingEscapeFraction))*massLossGradient + decayingEscapeFractionGradient*(1.0d0-self%massSplitting_)*massLossEnergy - decayingEscapeFractionGradient*heatingEnergy
    return
  end function DDMv2SpecificEnergyGradient

  logical function DDMv2SpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !!{
    Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileHeatingDDMv2), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    !$GLC attributes unused :: darkMatterProfileDMO_

    DDMv2SpecificEnergyIsEverywhereZero=(self%lifetime_ <= 0.0d0) .or. (self%massSplitting_ <= 0.0d0) .or. (.not. (self%heating_ .or. self%massLoss_))       
    return
  end function DDMv2SpecificEnergyIsEverywhereZero
