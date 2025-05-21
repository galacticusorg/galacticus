!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of dark matter halo profiles for self-interacting dark matter following the ``isothermal'' model of Jiang et
  al. (2022), including the effects of a baryonic potential.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <darkMatterProfile name="darkMatterProfileSIDMIsothermal">
    <description>
      A dark matter halo profile class that builds \refClass{massDistributionSphericalSIDMIsothermalBaryons} objects for
      isothermal SIDM profiles containing baryons.
    </description>
  </darkMatterProfile>
  !!]
  type, extends(darkMatterProfileClass) :: darkMatterProfileSIDMIsothermal
     !!{
     A dark matter halo profile class implementing profiles for self-interacting dark matter following the ``isothermal'' model of Jiang et al. (2022).
     !!}
     private
     class(darkMatterProfileClass ), pointer :: darkMatterProfile_  => null()
     class(darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
   contains
     final     ::        sidmIsothermalDestructor
     procedure :: get => sidmIsothermalGet
  end type darkMatterProfileSIDMIsothermal

  interface darkMatterProfileSIDMIsothermal
     !!{
     Constructors for the \refClass{darkMatterProfileSIDMIsothermal} dark matter halo profile class.
     !!}
     module procedure sidmIsothermalConstructorParameters
     module procedure sidmIsothermalConstructorInternal
  end interface darkMatterProfileSIDMIsothermal

contains

  function sidmIsothermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileSIDMIsothermal} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (darkMatterProfileSIDMIsothermal)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(darkMatterParticleClass        ), pointer       :: darkMatterParticle_
    class(darkMatterProfileClass         ), pointer       :: darkMatterProfile_

    !![
    <objectBuilder class="darkMatterParticle" name="darkMatterParticle_" source="parameters"/>
    <objectBuilder class="darkMatterProfile"  name="darkMatterProfile_"  source="parameters"/>
    !!]
    self=darkMatterProfileSIDMIsothermal(darkMatterProfile_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"/>
    <objectDestructor name="darkMatterProfile_" />
    !!]
    return
  end function sidmIsothermalConstructorParameters

  function sidmIsothermalConstructorInternal(darkMatterProfile_,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileSIDMIsothermal} dark matter profile class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    type (darkMatterProfileSIDMIsothermal)                        :: self
    class(darkMatterParticleClass        ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterProfileClass         ), intent(in   ), target :: darkMatterProfile_
    !![
    <constructorAssign variables="*darkMatterProfile_, *darkMatterParticle_"/>
    !!]

    ! Validate the dark matter particle type.
    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! This is as expected.
    class default
       call Error_Report('SIDM isothermal dark matter profile expects a self-interacting dark matter particle'//{introspection:location})
    end select
    return
  end function sidmIsothermalConstructorInternal

  subroutine sidmIsothermalDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileSIDMIsothermal} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileSIDMIsothermal), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_"/>
    <objectDestructor name="self%darkMatterProfile_" />
    !!]
    return
  end subroutine sidmIsothermalDestructor

  function sidmIsothermalGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                         , massTypeDark                        , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalSIDMIsothermalBaryons, kinematicsDistributionSIDMIsothermal, nonAnalyticSolversNumerical, massDistributionSpherical, &
         &                                    sphericalSIDMIsothermalBaryonsInitializor
    implicit none
    class    (massDistributionClass                    ), pointer                 :: massDistribution_
    type     (kinematicsDistributionSIDMIsothermal     ), pointer                 :: kinematicsDistribution_
    class    (darkMatterProfileSIDMIsothermal          ), intent(inout), target   :: self
    type     (treeNode                                 ), intent(inout), target   :: node
    type     (enumerationWeightByType                  ), intent(in   ), optional :: weightBy
    integer                                             , intent(in   ), optional :: weightIndex
    class    (massDistributionClass                    ), pointer                 :: massDistributionBaryonic, massDistributionDecorated
    class    (nodeComponentBasic                       ), pointer                 :: basic
    procedure(sphericalSIDMIsothermalBaryonsInitializor), pointer                 :: initializationFunction
    class    (*                                        ), pointer                 :: initializationSelf      , initializationArgument
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalSIDMIsothermalBaryons :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalSIDMIsothermalBaryons)
       massDistributionDecorated => self%darkMatterProfile_%get  (node,weightBy,weightIndex)
       massDistributionBaryonic  => null                         (                         )
       basic                     => node                   %basic(                         )
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          initializationFunction => sidmIsothermalInitialize
          initializationSelf     => self
          initializationArgument => node
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalSIDMIsothermalBaryons(                                                              &amp;
	      &amp;                                          timeAge                 =basic%time                       (), &amp;
	      &amp;                                          nonAnalyticSolver       =      nonAnalyticSolversNumerical  , &amp;
	      &amp;                                          massDistribution_       =      massDistributionDecorated    , &amp;
	      &amp;                                          massDistributionBaryonic=      massDistributionBaryonic     , &amp;
	      &amp;                                          darkMatterParticle_     =self %darkMatterParticle_          , &amp;
	      &amp;                                          initializationFunction  =      initializationFunction       , &amp;
	      &amp;                                          initializationSelf      =      initializationSelf           , &amp;
	      &amp;                                          initializationArgument  =      initializationArgument       , &amp;
              &amp;                                          componentType           =      componentTypeDarkHalo        , &amp;
              &amp;                                          massType                =      massTypeDark                   &amp;
              &amp;                                         )
	    </constructor>
          </referenceConstruct>
          !!]
       class default
          call Error_Report('expected a spherical mass distribution'//{introspection:location})
       end select
       !![
       <objectDestructor name="massDistributionDecorated"/>
       !!]
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionSIDMIsothermal( &amp;
	 &amp;                              )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function sidmIsothermalGet

  subroutine sidmIsothermalInitialize(self,node,massDistributionBaryonic)
    !!{
    Initialize the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : massTypeBaryonic
    implicit none
    class(*                    ), intent(inout), target  :: self                    , node
    class(massDistributionClass), intent(  out), pointer :: massDistributionBaryonic

    select type (self)
    type is (darkMatterProfileSIDMIsothermal)
       select type (node)
       type is (treeNode)
          massDistributionBaryonic => node%massDistribution(massType=massTypeBaryonic)
        class default
          call Error_Report('unexpected class'//{introspection:location})
       end select
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end subroutine sidmIsothermalInitialize
