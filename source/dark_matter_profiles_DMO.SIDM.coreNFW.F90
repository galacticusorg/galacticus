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
  An implementation of a cored-NFW dark matter halo profile to approximate the effects of SIDM based on the model of \cite{jiang_semi-analytic_2023}.
  !!}

  use :: Dark_Matter_Particles  , only : darkMatterParticleClass
  use :: Dark_Matter_Halo_Scales, only : darkmatterHaloScaleClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOSIDMCoreNFW">
    <description>
      Cored-NFW dark matter halo profiles to approximate the effects of SIDM based on the model of \cite{jiang_semi-analytic_2023}
      are built via \refClass{massDistributionSphericalSIDMCoreNFW} objects.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOSIDMCoreNFW
     !!{
     A dark matter halo profile class implementing a cored-NFW dark matter halo profile to approximate the effects of SIDM based
     on the model of Jiang et al. (2022). The profile is defined by the enclosed mass, with (Jiang et al. 2022):
     \begin{equation}
       M(r) = M_\mathrm{NFW}(r) \mathrm{tanh}\left(\frac{r}{r_\mathrm{c}}\right),
     \end{equation}
     where $r_\mathrm{c} = \alpha r_1$ is a characteristic core size related to the interaction radius $r_1$ by a constant factor
     $\alpha = ${\normalfont \ttfamily [factorRadiusCore]}.
     !!}
     private
     class           (darkMatterParticleClass  ), pointer :: darkMatterParticle_   => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     double precision                                     :: factorRadiusCore
   contains
     final     ::        sidmCoreNFWDestructor
     procedure :: get => sidmCoreNFWGet
  end type darkMatterProfileDMOSIDMCoreNFW

  interface darkMatterProfileDMOSIDMCoreNFW
     !!{
     Constructors for the \refClass{darkMatterProfileDMOSIDMCoreNFW} dark matter halo profile class.
     !!}
     module procedure sidmCoreNFWConstructorParameters
     module procedure sidmCoreNFWConstructorInternal
  end interface darkMatterProfileDMOSIDMCoreNFW

contains

  function sidmCoreNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOSIDMCoreNFW} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOSIDMCoreNFW)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (darkMatterParticleClass        ), pointer       :: darkMatterParticle_
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    double precision                                                 :: factorRadiusCore

    !![
    <inputParameter>
      <name>factorRadiusCore</name>
      <defaultValue>0.45d0</defaultValue>
      <defaultSource>Jiang et al. (2022)</defaultSource>
      <source>parameters</source>
      <description>The factor $\alpha$ appearing in the definition of the core radius, $r_\mathrm{c}=\alpha r_1$ where $r_1$ is the radius at which an SIDM particle has had, on average, 1 interaction.</description>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOSIDMCoreNFW(factorRadiusCore,darkMatterHaloScale_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function sidmCoreNFWConstructorParameters

  function sidmCoreNFWConstructorInternal(factorRadiusCore,darkMatterHaloScale_,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOSIDMCoreNFW} dark matter profile class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    type            (darkMatterProfileDMOSIDMCoreNFW)                        :: self
    class           (darkMatterParticleClass        ), intent(in   ), target :: darkMatterParticle_
    class           (darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                 , intent(in   )         :: factorRadiusCore
    !![
    <constructorAssign variables="factorRadiusCore, *darkMatterHaloScale_, *darkMatterParticle_"/>
    !!]

    ! Validate the dark matter particle type.
    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! This is as expected.
    class default
       call Error_Report('transfer function expects a self-interacting dark matter particle'//{introspection:location})
    end select
    ! Construct an NFW profile.
    allocate(darkMatterProfileDMONFW :: self%darkMatterProfileDMO_)
    select type (darkMatterProfileDMO_ => self%darkMatterProfileDMO_)
    type is (darkMatterProfileDMONFW)
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="darkMatterProfileDMO_"  object="darkMatterProfileDMO_">
	 <constructor>
	   darkMatterProfileDMONFW(                                                           &amp;
	    &amp;                  velocityDispersionUseSeriesExpansion=.true.              , &amp;
	    &amp;                  darkMatterHaloScale_                =darkMatterHaloScale_  &amp;
	    &amp;                 )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function sidmCoreNFWConstructorInternal

  subroutine sidmCoreNFWDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOSIDMCoreNFW} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOSIDMCoreNFW), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterParticle_"  />
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine sidmCoreNFWDestructor

  function sidmCoreNFWGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo               , massTypeDark                       , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalSIDMCoreNFW, kinematicsDistributionCollisionless, massDistributionNFW, nonAnalyticSolversNumerical
    implicit none
    class           (massDistributionClass              ), pointer                 :: massDistribution_
    type            (kinematicsDistributionCollisionless), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOSIDMCoreNFW    ), intent(inout)           :: self
    type            (treeNode                           ), intent(inout)           :: node
    type            (enumerationWeightByType            ), intent(in   ), optional :: weightBy
    integer                                              , intent(in   ), optional :: weightIndex
    class           (massDistributionClass              ), pointer                 :: massDistributionDecorated
    class           (nodeComponentBasic                 ), pointer                 :: basic
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalSIDMCoreNFW :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalSIDMCoreNFW)
       massDistributionDecorated => self%darkMatterProfileDMO_%get  (node,weightBy,weightIndex)
       basic                     => node                      %basic(                         )
       select type (massDistributionDecorated)
       class is (massDistributionNFW)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalSIDMCoreNFW(                                                         &amp;
	      &amp;                                factorRadiusCore   =self %factorRadiusCore             , &amp;
	      &amp;                                timeAge            =basic%time                       (), &amp;
	      &amp;                                nonAnalyticSolver  =      nonAnalyticSolversNumerical  , &amp;
	      &amp;                                massDistribution_  =      massDistributionDecorated    , &amp;
	      &amp;                                darkMatterParticle_=self %darkMatterParticle_          , &amp;
              &amp;                                componentType      =      componentTypeDarkHalo        , &amp;
              &amp;                                massType           =      massTypeDark                   &amp;
              &amp;                               )
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
        kinematicsDistributionCollisionless( &amp;
	 &amp;                             )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function sidmCoreNFWGet
