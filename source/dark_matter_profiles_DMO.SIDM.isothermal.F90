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
  An implementation of dark matter halo profiles for self-interacting dark matter following the ``isothermal'' model of \cite{jiang_semi-analytic_2023}.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOSIDMIsothermal">
    <description>
      Dark matter halo profiles for self-interacting dark matter following the ``isothermal'' model of
      \cite{jiang_semi-analytic_2023} are built via the \refClass{massDistributionSphericalSIDMIsothermal} class.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOSIDMIsothermal
     !!{
     A dark matter halo profile class implementing profiles for self-interacting dark matter following the ``isothermal'' model of \cite{jiang_semi-analytic_2023}.
     !!}
     private
     class(darkMatterParticleClass  ), pointer :: darkMatterParticle_   => null()
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
   contains
     final     ::        sidmIsothermalDestructor
     procedure :: get => sidmIsothermalGet    
  end type darkMatterProfileDMOSIDMIsothermal

  interface darkMatterProfileDMOSIDMIsothermal
     !!{
     Constructors for the {\normalfont \ttfamily sidmIsothermal} dark matter halo profile class.
     !!}
     module procedure sidmIsothermalConstructorParameters
     module procedure sidmIsothermalConstructorInternal
  end interface darkMatterProfileDMOSIDMIsothermal

contains

  function sidmIsothermalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sidmIsothermal} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMOSIDMIsothermal)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(darkMatterParticleClass           ), pointer       :: darkMatterParticle_
    class(darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_

    !![
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOSIDMIsothermal(darkMatterProfileDMO_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"  />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function sidmIsothermalConstructorParameters

  function sidmIsothermalConstructorInternal(darkMatterProfileDMO_,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sidmIsothermal} dark matter profile class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    type (darkMatterProfileDMOSIDMIsothermal)                        :: self
    class(darkMatterParticleClass           ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="*darkMatterProfileDMO_, *darkMatterParticle_"/>
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
    Destructor for the {\normalfont \ttfamily sidmIsothermal} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOSIDMIsothermal), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterParticle_"  />
    !!]
    return
  end subroutine sidmIsothermalDestructor

  function sidmIsothermalGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                  , massTypeDark                        , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalSIDMIsothermal, kinematicsDistributionSIDMIsothermal, nonAnalyticSolversNumerical, massDistributionSpherical
    implicit none
    class           (massDistributionClass               ), pointer                 :: massDistribution_
    type            (kinematicsDistributionSIDMIsothermal), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOSIDMIsothermal  ), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    type            (enumerationWeightByType             ), intent(in   ), optional :: weightBy
    integer                                               , intent(in   ), optional :: weightIndex
    class           (massDistributionClass               ), pointer                 :: massDistributionDecorated
    class           (nodeComponentBasic                  ), pointer                 :: basic
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalSIDMIsothermal :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalSIDMIsothermal)
       massDistributionDecorated => self%darkMatterProfileDMO_%get  (node,weightBy,weightIndex)
       basic                     => node                      %basic(                         )
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalSIDMIsothermal(                                                         &amp;
	      &amp;                                   timeAge            =basic%time                       (), &amp;
	      &amp;                                   nonAnalyticSolver  =      nonAnalyticSolversNumerical  , &amp;
	      &amp;                                   massDistribution_  =      massDistributionDecorated    , &amp;
	      &amp;                                   darkMatterParticle_=self %darkMatterParticle_          , &amp;
              &amp;                                   componentType      =      componentTypeDarkHalo        , &amp;
              &amp;                                   massType           =      massTypeDark                   &amp;
              &amp;                                  )
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
