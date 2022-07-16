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
     logical                                                :: heating, massLoss
     double precision                                       :: lifetime_, massSplitting_, gamma
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
    logical                                                              :: heating   , massLoss
    double precision                                                     :: gamma
         
    !![
    <inputParameter>
      <name>heating</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>Indicates if heating from decays should be considered.</description>
    </inputParameter>
    <inputParameter>
      <name>massLoss</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>Indicates if mass loss effects should be considered.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <defaultValue>+1.0d0</defaultValue>
      <source>parameters</source>
      <description>Free parameter in mass loss heating term.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    !!]
    self=darkMatterProfileHeatingDDMv2(darkMatterParticle_, heating, massLoss, gamma)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"  />
    !!]
    return
  end function DDMv2ConstructorParameters

  function DDMv2ConstructorInternal(darkMatterParticle_, heating, massLoss, gamma) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily DDM} dark matter profile heating scales class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleDecayingDarkMatter
    implicit none
    type            (darkMatterProfileHeatingDDMv2)                        :: self
    class           (darkMatterParticleClass      ), intent(in   ), target :: darkMatterParticle_
    logical                                        , intent(in   )         :: heating , massLoss
    double precision                               , intent(in   )         :: gamma
    !![
    <constructorAssign variables="*darkMatterParticle_, heating, massLoss, gamma"/>
    !!]
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleDecayingDarkMatter)
       self%lifetime_ = darkMatterParticle_%lifetime()
       self%massSplitting_ = darkMatterParticle_%massSplitting()
    class default
       ! No decays.
       self%lifetime_=-1.0d0
       self%massSplitting_=0.0d0
    end select
    return
  end function DDMv2ConstructorInternal

  double precision function DDMv2SpecificEnergy(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Physical, only : speedLight
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class           (darkMatterProfileHeatingDDMv2), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    double precision                             , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass  ), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic         ), pointer       :: basic
    double precision :: heatingEnergy, massLossEnergy

    basic             => node%basic()
    if (self%heating) then
      heatingEnergy = +0.5d0*(                                    &
           &          +1.0d0 - exp(-basic%time() / self%lifetime_)&
           &                 )                                    &
           &          *self%massSplitting_**2                     &
           &              *(speedLight/kilo)**2
    else
      heatingEnergy = 0.0d0
    end if
    if (self%massLoss) then
      massLossEnergy = (+0.5d0*self%massSplitting_**2*(speedLight/kilo)**2&
           &            +self%gamma*self%massSplitting_                   &
           &            *darkMatterProfileDMO_%potential(node, radius)    &
           &           )*(+1.0d0 - exp(-basic%time() / self%lifetime_))
    else
      massLossEnergy = 0.0d0
    end if
    DDMv2SpecificEnergy = heatingEnergy + massLossEnergy
    return
  end function DDMv2SpecificEnergy

  double precision function DDMv2SpecificEnergyGradient(self,node,radius,darkMatterProfileDMO_)
    !!{
    Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (darkMatterProfileHeatingDDMv2), intent(inout) :: self
    type            (treeNode                     ), intent(inout) :: node
    double precision                               , intent(in   ) :: radius
    class           (darkMatterProfileDMOClass    ), intent(inout) :: darkMatterProfileDMO_
    class           (nodeComponentBasic           ), pointer       :: basic
    
    basic             => node%basic()
    if (self%massLoss) then
      DDMv2SpecificEnergyGradient=self%gamma*self%massSplitting_                                                                &
             &                    *(-darkMatterProfileDMO_%potential(node, radius)/radius                                       &
             &                    +gravitationalConstantGalacticus*4.0d0*Pi*darkMatterProfileDMO_%density(node, radius)*radius) &
             &                    *(+1.0d0 - exp(-basic%time() / self%lifetime_))
    else
      DDMv2SpecificEnergyGradient=+0.0d0
    end if
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

    DDMv2SpecificEnergyIsEverywhereZero=(self%lifetime_ <= 0.0d0) .or. (self%massSplitting_ <= 0.0d0) .or. (.not. (self%heating .or. self%massLoss))
    return
  end function DDMv2SpecificEnergyIsEverywhereZero
