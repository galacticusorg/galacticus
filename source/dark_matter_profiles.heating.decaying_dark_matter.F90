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

  !!{
  A dark matter halo profile heating class which accounts for heating from decaying dark matter.
  !!}
  
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Dark_Matter_Particles  , only : darkMatterParticleClass

  !![
  <darkMatterProfileHeating name="darkMatterProfileHeatingDecayingDarkMatter">
   <description>
    A dark matter profile heating class that constructs \refClass{massDistributionHeatingDecayingDarkMatter} objects to compute heating due to
    decaying dark matter.
   </description>
  </darkMatterProfileHeating>
  !!]
  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingDecayingDarkMatter
     !!{
     A dark matter profile heating class which accounts for heating due to decaying dark matter.
     !!}
     private
     class           (darkMatterParticleClass  ), pointer :: darkMatterParticle_  => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_ => null()
     double precision                                     :: gamma
     logical                                              :: includeKickHeating
   contains
     final     ::        decayingDarkMatterDestructor
     procedure :: get => decayingDarkMatterGet
  end type darkMatterProfileHeatingDecayingDarkMatter

  interface darkMatterProfileHeatingDecayingDarkMatter
     !!{
     Constructors for the \refClass{darkMatterProfileHeatingDecayingDarkMatter} dark matter profile heating class.
     !!}
     module procedure decayingDarkMatterConstructorParameters
     module procedure decayingDarkMatterConstructorInternal
  end interface darkMatterProfileHeatingDecayingDarkMatter

contains

  function decayingDarkMatterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileHeatingDecayingDarkMatter} dark matter profile heating scales class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileHeatingDecayingDarkMatter), target        :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                  ), pointer       :: darkMatterHaloScale_
    class           (darkMatterParticleClass                   ), pointer       :: darkMatterParticle_
    double precision                                                            :: gamma
    logical                                                                     :: includeKickHeating

    !![
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
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileHeatingDecayingDarkMatter(gamma,includeKickHeating,darkMatterParticle_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function decayingDarkMatterConstructorParameters

  function decayingDarkMatterConstructorInternal(gamma,includeKickHeating,darkMatterParticle_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileHeatingDecayingDarkMatter} dark matter profile heating scales class.
    !!}
    implicit none
    type (darkMatterProfileHeatingDecayingDarkMatter)                        :: self
    class(darkMatterHaloScaleClass                  ), intent(in   ), target :: darkMatterHaloScale_
    class(darkMatterParticleClass                   ), intent(in   ), target :: darkMatterParticle_
    double precision                                 , intent(in   )         :: gamma
    logical                                          , intent(in   )         :: includeKickHeating
    !![
    <constructorAssign variables="gamma, includeKickHeating, *darkMatterParticle_, *darkMatterHaloScale_"/>
    !!]
    
    return
  end function decayingDarkMatterConstructorInternal

  subroutine decayingDarkMatterDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileHeatingDecayingDarkMatter} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileHeatingDecayingDarkMatter), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterParticle_"  />
    !!]
    return
  end subroutine decayingDarkMatterDestructor

  function decayingDarkMatterGet(self,node) result(massDistributionHeating_)
    !!{
    Return the dark matter mass distribution heating for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: Mass_Distributions, only : massDistributionHeatingDecayingDarkMatter
    implicit none
    class           (massDistributionHeatingClass              ), pointer       :: massDistributionHeating_
    class           (darkMatterProfileHeatingDecayingDarkMatter), intent(inout) :: self
    type            (treeNode                                  ), intent(inout) :: node
    class           (nodeComponentBasic                        ), pointer       :: basic
    double precision                                            , parameter     :: factorRadiusEscape       =1000.0d0
    double precision                                                            :: radiusEscape

    ! Create the mass distribution.
    allocate(massDistributionHeatingDecayingDarkMatter :: massDistributionHeating_)
    select type(massDistributionHeating_)
    type is (massDistributionHeatingDecayingDarkMatter)
       basic        =>  node                     %basic             (    )
       radiusEscape =  +                          factorRadiusEscape       &
            &          *self%darkMatterHaloScale_%radiusVirial      (node)
       !![
       <referenceConstruct object="massDistributionHeating_">
	 <constructor>
           massDistributionHeatingDecayingDarkMatter(                                                 &amp;
            &amp;                                    radiusEscape       =      radiusEscape         , &amp;
            &amp;                                    time               =basic%time               (), &amp;
            &amp;                                    gamma              =self %gamma                , &amp;
            &amp;                                    includeKickHeating =self %includeKickHeating   , &amp;
            &amp;                                    darkMatterParticle_=self %darkMatterParticle_    &amp;
            &amp;                                   )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function decayingDarkMatterGet

