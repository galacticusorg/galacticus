!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Dark_Matter_Particles  , only : darkMatterParticleClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMODecaying">
   <description>
     A dark matter profile DMO class which builds \refClass{massDistributionSphericalDecaying} objects to account for dark matter
     particle decays in some other dark matter profile.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMODecaying
     !!{
     A dark matter halo profile class implementing decaying dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_                      => null()
     class           (darkMatterParticleClass  ), pointer :: darkMatterParticle_                        => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_                       => null()
     double precision                                     :: toleranceRelativePotential                          , toleranceRelativeVelocityDispersion, &
          &                                                  toleranceRelativeVelocityDispersionMaximum
     logical                                              :: tolerateVelocityMaximumFailure                      , toleratePotentialIntegrationFailure, &
          &                                                  tolerateEnclosedMassIntegrationFailure
   contains
     final     ::        decayingDestructor
     procedure :: get => decayingGet
  end type darkMatterProfileDMODecaying

  interface darkMatterProfileDMODecaying
     !!{
     Constructors for the \refClass{darkMatterProfileDMODecaying} dark matter halo profile class.
     !!}
     module procedure decayingConstructorParameters
     module procedure decayingConstructorInternal
  end interface darkMatterProfileDMODecaying

contains

  function decayingConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily decaying} dark matter halo profile class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (darkMatterProfileDMODecaying)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass   ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_
    class           (darkMatterParticleClass     ), pointer       :: darkMatterParticle_
    double precision                                              :: toleranceRelativePotential                , toleranceRelativeVelocityDispersion, &
         &                                                           toleranceRelativeVelocityDispersionMaximum
    logical                                                       :: tolerateVelocityMaximumFailure            , toleratePotentialIntegrationFailure, &
         &                                                           tolerateEnclosedMassIntegrationFailure

    !![
    <inputParameter>
      <name>toleranceRelativeVelocityDispersion</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the velocity dispersion.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersionMaximum</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum relative tolerance to use in numerical solutions for the velocity dispersion.</description>
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
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=darkMatterProfileDMODecaying(toleranceRelativePotential,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,tolerateEnclosedMassIntegrationFailure,darkMatterParticle_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterParticle_"  />
    !!]
    return
  end function decayingConstructorParameters

  function decayingConstructorInternal(toleranceRelativePotential,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,tolerateEnclosedMassIntegrationFailure,darkMatterParticle_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMODecaying} dark matter profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversIsValid
    use :: Error             , only : Error_Report
    implicit none
    type            (darkMatterProfileDMODecaying)                        :: self
    class           (darkMatterProfileDMOClass   ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass    ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterParticleClass     ), intent(in   ), target :: darkMatterParticle_
    double precision                              , intent(in   )         :: toleranceRelativePotential                , toleranceRelativeVelocityDispersion, &
         &                                                                   toleranceRelativeVelocityDispersionMaximum
    logical                                       , intent(in   )         :: tolerateVelocityMaximumFailure            , toleratePotentialIntegrationFailure, &
         &                                                                   tolerateEnclosedMassIntegrationFailure
    !![
    <constructorAssign variables="toleranceRelativePotential, toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, tolerateVelocityMaximumFailure, toleratePotentialIntegrationFailure, tolerateEnclosedMassIntegrationFailure, *darkMatterParticle_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function decayingConstructorInternal

  subroutine decayingDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMODecaying} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMODecaying), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterParticle_"  />
    !!]
    return
  end subroutine decayingDestructor

  function decayingGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo            , massTypeDark                       , weightByMass
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Mass_Distributions        , only : massDistributionSphericalDecaying, kinematicsDistributionCollisionless, massDistributionSpherical
    implicit none
    class           (massDistributionClass              ), pointer                 :: massDistribution_
    type            (kinematicsDistributionCollisionless), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMODecaying       ), intent(inout)           :: self
    type            (treeNode                           ), intent(inout)           :: node
    type            (enumerationWeightByType            ), intent(in   ), optional :: weightBy
    integer                                              , intent(in   ), optional :: weightIndex
    double precision                                     , parameter               :: factorRadiusEscape       =1000.0d0
    class           (nodeComponentBasic                 ), pointer                 :: basic
    class           (massDistributionClass              ), pointer                 :: massDistributionDecorated
    double precision                                                               :: radiusEscape
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalDecaying :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalDecaying)
       massDistributionDecorated => self%darkMatterProfileDMO_%get(node,weightBy,weightIndex)
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          basic        =>  node                     %basic             (    )
          radiusEscape =  +                          factorRadiusEscape       &
               &          *self%darkMatterHaloScale_%radiusVirial      (node)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalDecaying(                                                                                       &amp;
	       &amp;                            toleranceRelativePotential            =self %toleranceRelativePotential              , &amp;
	       &amp;                            tolerateVelocityMaximumFailure        =self %tolerateVelocityMaximumFailure          , &amp;
	       &amp;                            toleratePotentialIntegrationFailure   =self %toleratePotentialIntegrationFailure     , &amp;
	       &amp;                            tolerateEnclosedMassIntegrationFailure=self %tolerateEnclosedMassIntegrationFailure  , &amp;
	       &amp;                            radiusEscape                          =      radiusEscape                            , &amp;
	       &amp;                            time                                  =basic%time                                  (), &amp;
	       &amp;                            darkMatterParticle_                   =self %darkMatterParticle_                     , &amp;
               &amp;                            massDistribution_                     =      massDistributionDecorated               , &amp;
               &amp;                            componentType                         =      componentTypeDarkHalo                   , &amp;
               &amp;                            massType                              =      massTypeDark                              &amp;
               &amp;                           )
	    </constructor>
	  </referenceConstruct>
	  <objectDestructor name="massDistributionDecorated"/>
          !!]
       class default
          call Error_Report('expected a spherical mass distribution'//{introspection:location})
       end select
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionCollisionless(                                                                                            &amp;
         &amp;                              toleranceRelativeVelocityDispersion       =self%toleranceRelativeVelocityDispersion       , &amp; 
         &amp;                              toleranceRelativeVelocityDispersionMaximum=self%toleranceRelativeVelocityDispersionMaximum  &amp; 
	 &amp;                             )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function decayingGet
