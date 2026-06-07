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
  An implementation of dark matter halo profiles for self-interacting dark matter following the ``SIDM\_parametric'' model of \cite{yang_parametric_2024}.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass
  use :: Dark_Matter_Halo_Scales   , only : darkMatterHaloScaleClass

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOSIDMParametric">
    <description>
      Dark matter halo profiles for self-interacting dark matter following the ``SIDM\_parametric'' model of
      \cite{yang_parametric_2024} are built via the \refClass{massDistributionSIDMParametricProfile} class.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOSIDMParametric
     !!{
     A dark matter halo profile class implementing profiles for self-interacting dark matter following the ``SIDM\_parametric'' model
     of \cite{yang_parametric_2024}.
     !!}
     private
     class(darkMatterParticleClass  ), pointer :: darkMatterParticle_   => null()
     class(darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     integer                                   :: RhosSIDMID, RsSIDMID, RcSIDMID
     double precision                          :: beta
   contains
     final     ::        sidmParametricDestructor
     procedure :: get => sidmParametricGet    
  end type darkMatterProfileDMOSIDMParametric

  interface darkMatterProfileDMOSIDMParametric
     !!{
     Constructors for the {\normalfont \ttfamily sidmParametric} dark matter halo profile class.
     !!}
     module procedure sidmParametricConstructorParameters
     module procedure sidmParametricConstructorInternal
  end interface darkMatterProfileDMOSIDMParametric

contains

  function sidmParametricConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily sidmParametric} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (darkMatterProfileDMOSIDMParametric)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(darkMatterParticleClass           ), pointer       :: darkMatterParticle_
    class(darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    double precision                                         :: beta

    !![
    <inputParameter>
      <name>beta</name>
      <defaultValue>4.0d0</defaultValue>
      <description>The value $\beta$ in a SIDMParametric-model mass distribution.</description>
      <source>parameters</source>
    </inputParameter>
    !!]

    !![
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOSIDMParametric(beta,darkMatterParticle_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"  />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function sidmParametricConstructorParameters

  function sidmParametricConstructorInternal(beta,darkMatterParticle_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily sidmParametric} dark matter profile class.
    !!}
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter
    use :: Error                , only : Error_Report
    implicit none
    type (darkMatterProfileDMOSIDMParametric)                        :: self
    class(darkMatterParticleClass           ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    double precision                         , intent(in   )         :: beta

    !![
    <constructorAssign variables="beta, *darkMatterParticle_, *darkMatterHaloScale_"/>
    !!]

    ! Validate the dark matter particle type.
    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       ! This is as expected.
    class default
       call Error_Report('SIDM parametric dark matter profile expects a self-interacting dark matter particle'//{introspection:location})
    end select    

    !![
    <addMetaProperty component="darkMatterProfile" name="RhosSIDM" id="self%RhosSIDMID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="RsSIDM"   id="self%RsSIDMID"   isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="RcSIDM"   id="self%RcSIDMID"   isEvolvable="yes" isCreator="no"/>
    !!]

    return
  end function sidmParametricConstructorInternal

  subroutine sidmParametricDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily sidmParametric} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOSIDMParametric), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_"  />
    <objectDestructor name="self%darkMatterHaloScale_"  />
    !!]
    return
  end subroutine sidmParametricDestructor

  function sidmParametricGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic                     , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                  , massTypeDark                        , weightByMass
    use :: Mass_Distributions        , only : massDistributionSIDMParametricProfile, massDistributionNFW, kinematicsDistributionCollisionlessTabulated, kinematicsDistributionNFW, kinematicsDistributionClass
    implicit none
    class           (massDistributionClass               ), pointer                 :: massDistribution_
    class            (kinematicsDistributionClass ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOSIDMParametric  ), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    type            (enumerationWeightByType             ), intent(in   ), optional :: weightBy
    integer                                               , intent(in   ), optional :: weightIndex
    class           (nodeComponentBasic                  ), pointer                 :: basic
    class           (nodeComponentDarkMatterProfile),       pointer                 :: darkMatterProfile

    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    ! Use NFW for zero core radius.
    if (darkMatterProfile%floatRank0MetaPropertyGet(self%RcSIDMID) == 0.0d0) then
       allocate(massDistributionNFW :: massDistribution_)
       select type(massDistribution_)
       type is (massDistributionNFW)
          !![
	  <referenceConstruct object="massDistribution_">
            <constructor>
              massDistributionNFW(                                                         &amp;
	      &amp;               densityNormalization= darkMatterProfile%floatRank0MetaPropertyGet(self%RhosSIDMID), &amp;
              &amp;               scaleLength         = darkMatterProfile%floatRank0MetaPropertyGet(self%RsSIDMID), &amp;
              &amp;               componentType       =      componentTypeDarkHalo        , &amp;
              &amp;               massType            =      massTypeDark                   &amp;
              &amp;              )
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
       allocate(kinematicsDistributionNFW :: kinematicsDistribution_)
       select type (kinematicsDistribution_)
       type is (kinematicsDistributionNFW)
          !![
	  <referenceConstruct object="kinematicsDistribution_">
	    <constructor>
              kinematicsDistributionNFW(useSeriesApproximation=.true.)
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
    else
       ! Use the SIDM parametric profile.
       allocate(massDistributionSIDMParametricProfile :: massDistribution_)
       select type(massDistribution_)
       type is (massDistributionSIDMParametricProfile)
          !![
	  <referenceConstruct object="massDistribution_">
            <constructor>
              massDistributionSIDMParametricProfile(                                                         &amp;
	      &amp;                                   beta                = self%beta, &amp;
              &amp;                                   densityNormalization= darkMatterProfile%floatRank0MetaPropertyGet(self%RhosSIDMID), &amp;
              &amp;                                   radiusScale         = darkMatterProfile%floatRank0MetaPropertyGet(self%RsSIDMID), &amp;
	      &amp;                                   radiusCore          = darkMatterProfile%floatRank0MetaPropertyGet(self%RcSIDMID), &amp;
              &amp;                                   componentType       =      componentTypeDarkHalo        , &amp;
              &amp;                                   massType            =      massTypeDark                   &amp;
              &amp;                                  )
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
       allocate(kinematicsDistributionCollisionlessTabulated :: kinematicsDistribution_)
       select type (kinematicsDistribution_)
       type is (kinematicsDistributionCollisionlessTabulated)
          !![
	  <referenceConstruct object="kinematicsDistribution_">
	    <constructor>
              kinematicsDistributionCollisionlessTabulated( &amp;
	      &amp;                              )
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
    end if
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function sidmParametricGet
