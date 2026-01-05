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
  An implementation of non-dark-matter-only dark matter halo profiles which are unchanged from their dark-matter-only counterpart.
  !!}

  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <darkMatterProfile name="darkMatterProfileDarkMatterOnly">
   <description>An implementation of non-dark-matter-only dark matter halo profiles which are unchanged from their dark-matter-only counterpart.</description>
  </darkMatterProfile>
  !!]
  type, extends(darkMatterProfileClass) :: darkMatterProfileDarkMatterOnly
     !!{
     A dark matter halo profile class implementing non-dark-matter-only dark matter halo profiles which are unchanged from their dark-matter-only counterpart.
     !!}
     private
     class           (cosmologyParametersClass ), pointer :: cosmologyParameters_                           => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_                          => null()
     double precision                                     :: darkMatterFraction
     logical                                              :: chandrasekharIntegralComputeVelocityDispersion
   contains
     final     ::        darkMatterOnlyDestructor
     procedure :: get => darkMatterOnlyGet
  end type darkMatterProfileDarkMatterOnly

  interface darkMatterProfileDarkMatterOnly
     !!{
     Constructors for the \refClass{darkMatterProfileDarkMatterOnly} non-dark-matter-only dark matter halo profile class.
     !!}
     module procedure darkMatterOnlyConstructorParameters
     module procedure darkMatterOnlyConstructorInternal
  end interface darkMatterProfileDarkMatterOnly

contains

  function darkMatterOnlyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDarkMatterOnly} non-dark-matter-only dark matter halo profile class which takes
    a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (darkMatterProfileDarkMatterOnly)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    class  (cosmologyParametersClass       ), pointer       :: cosmologyParameters_
    class  (darkMatterProfileDMOClass      ), pointer       :: darkMatterProfileDMO_
    logical                                                 :: chandrasekharIntegralComputeVelocityDispersion

    !![
    <inputParameter>
      <name>chandrasekharIntegralComputeVelocityDispersion</name>
      <defaultValue>.true.</defaultValue>
      <description>
        If true, the Chandrasekhar integral is computed using the velocity dispersion, $\sigma_mathrm{r}(r)$. Otherwise, the
        velocity dispersion is approximated as $V_\mathrm{c}(r)/\sqrt{2}$.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=darkMatterProfileDarkMatterOnly(chandrasekharIntegralComputeVelocityDispersion,cosmologyParameters_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function darkMatterOnlyConstructorParameters

  function darkMatterOnlyConstructorInternal(chandrasekharIntegralComputeVelocityDispersion,cosmologyParameters_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDarkMatterOnly} dark matter profile class.
    !!}
    implicit none
    type   (darkMatterProfileDarkMatterOnly)                        :: self
    class  (cosmologyParametersClass       ), intent(in   ), target :: cosmologyParameters_
    class  (darkMatterProfileDMOClass      ), intent(in   ), target :: darkMatterProfileDMO_
    logical                                 , intent(in   )         :: chandrasekharIntegralComputeVelocityDispersion
    !![
    <constructorAssign variables="chandrasekharIntegralComputeVelocityDispersion, *cosmologyParameters_, *darkMatterProfileDMO_"/>
    !!]

    ! Evaluate the dark matter fraction.
    self%darkMatterFraction=+1.0d0                                   &
         &                  -self%cosmologyParameters_%OmegaBaryon() &
         &                  /self%cosmologyParameters_%OmegaMatter()
    return
  end function darkMatterOnlyConstructorInternal

  subroutine darkMatterOnlyDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDarkMatterOnly} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDarkMatterOnly), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine darkMatterOnlyDestructor

  function darkMatterOnlyGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : weightByMass
    use :: Mass_Distributions        , only : massDistributionSpherical, massDistributionSphericalScaler, kinematicsDistributionSphericalScaler, kinematicsDistributionClass
    use :: Error                     , only : Error_Report
    implicit none
    class  (massDistributionClass                ), pointer                 :: massDistribution_
    type   (kinematicsDistributionSphericalScaler), pointer                 :: kinematicsDistribution_
    class  (kinematicsDistributionClass          ), pointer                 :: kinematicsDistributionDMO
    class  (darkMatterProfileDarkMatterOnly      ), intent(inout), target   :: self
    type   (treeNode                             ), intent(inout), target   :: node
    type   (enumerationWeightByType              ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    class  (massDistributionClass                ), pointer                 :: massDistributionDMO
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Get the dark matter-only mass distribution.
    massDistributionDMO       => self               %darkMatterProfileDMO_%get                   (node,weightBy,weightIndex)
    kinematicsDistributionDMO => massDistributionDMO                      %kinematicsDistribution(                         )
    if (.not.associated(massDistributionDMO)) return
    select type (massDistributionDMO)
    class is (massDistributionSpherical)
       ! Create the mass distribution.
       allocate(massDistributionSphericalScaler :: massDistribution_)
       select type(massDistribution_)
       type is (massDistributionSphericalScaler)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
	      massDistributionSphericalScaler(                                                                                                    &amp;
	        &amp;                         factorScalingLength                           =     1.0d0                                         , &amp;
	        &amp;                         factorScalingMass                             =self%darkMatterFraction                            , &amp;
	        &amp;                         massDistribution_                             =     massDistributionDMO                           , &amp;
	        &amp;                         chandrasekharIntegralComputeVelocityDispersion=self%chandrasekharIntegralComputeVelocityDispersion  &amp;	      
	        &amp;                        )
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
       allocate(kinematicsDistribution_)
       !![
       <referenceConstruct object="kinematicsDistribution_">
	 <constructor>
           kinematicsDistributionSphericalScaler(                                                        &amp;
	        &amp;                            factorScalingLength    =     1.0d0                    , &amp;
	        &amp;                            factorScalingMass      =self%darkMatterFraction       , &amp;
	        &amp;                            kinematicsDistribution_=     kinematicsDistributionDMO  &amp;
                &amp;                           )
	 </constructor>
       </referenceConstruct>
       !!]
       call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
       !![
       <objectDestructor name="kinematicsDistribution_"/>
       !!]
    class default
       call Error_Report('a spherical mass distribution is required'//{introspection:location})
    end select
    !![
    <objectDestructor name="massDistributionDMO"      />
    <objectDestructor name="kinematicsDistributionDMO"/>
    !!]
    return
  end function darkMatterOnlyGet
