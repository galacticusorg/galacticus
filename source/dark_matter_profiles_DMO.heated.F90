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
  An implementation of heated dark matter halo profiles.
  !!}

  use :: Mass_Distributions, only : enumerationNonAnalyticSolversType

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOHeated">
   <description>
     A dark matter profile DMO class which builds \refClass{massDistributionSphericalHeated} objects to account for heating of
     some other dark matter profile.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOHeated
     !!{
     A dark matter halo profile class implementing heated dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_                  => null()
     class           (darkMatterProfileHeatingClass    ), pointer :: darkMatterProfileHeating_              => null()
     double precision                                             :: toleranceRelativeVelocityDispersion             , toleranceRelativeVelocityDispersionMaximum, &
          &                                                          fractionRadiusFinalSmall                        , toleranceRelativePotential
     type            (enumerationNonAnalyticSolversType)          :: nonAnalyticSolver
     logical                                                      :: velocityDispersionApproximate                   , tolerateVelocityMaximumFailure            , &
          &                                                          tolerateEnclosedMassIntegrationFailure          , tolerateVelocityDispersionFailure         , &
          &                                                          toleratePotentialIntegrationFailure
   contains
     final     ::        heatedDestructor
     procedure :: get => heatedGet
  end type darkMatterProfileDMOHeated

  interface darkMatterProfileDMOHeated
     !!{
     Constructors for the \refClass{darkMatterProfileDMOHeated} dark matter halo profile class.
     !!}
     module procedure heatedConstructorParameters
     module procedure heatedConstructorInternal
  end interface darkMatterProfileDMOHeated

contains

  function heatedConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily heated} dark matter halo profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    use :: Input_Parameters  , only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOHeated    )                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass     ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterProfileHeatingClass ), pointer       :: darkMatterProfileHeating_
    type            (varying_string                )                :: nonAnalyticSolver
    logical                                                         :: velocityDispersionApproximate         , tolerateVelocityMaximumFailure            , &
         &                                                             tolerateEnclosedMassIntegrationFailure, tolerateVelocityDispersionFailure         , &
         &                                                             toleratePotentialIntegrationFailure
    double precision                                                :: toleranceRelativeVelocityDispersion   , toleranceRelativeVelocityDispersionMaximum, &
         &                                                             fractionRadiusFinalSmall              , toleranceRelativePotential

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityDispersionApproximate</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, radial velocity dispersion is computed using an approximate method in which we assume that $\sigma_\mathrm{r}^2(r) \rightarrow \sigma_\mathrm{r}^2(r) - (2/3) \epsilon(r)$, where $\epsilon(r)$ is the specific heating energy. If {\normalfont \ttfamily false} then radial velocity dispersion is computed by numerically solving the Jeans equation.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerateEnclosedMassIntegrationFailure</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, tolerate failures to find the mass enclosed as a function of radius.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerateVelocityDispersionFailure</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, tolerate failures to compute the velocity dispersion.</description>
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
    <inputParameter>
      <name>toleranceRelativeVelocityDispersion</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersionMaximum</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum relative tolerance to use in numerical solutions for the velocity dispersion in dark-matter-only density profiles.</description>
    </inputParameter>
    <inputParameter>
      <name>fractionRadiusFinalSmall</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The initial radius is limited to be no smaller than this fraction of the final radius. This can help avoid problems in profiles that are extremely close to being disrupted.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativePotential</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum allowed relative tolerance to use in numerical solutions for the gravitational potential in dark-matter-only density profiles before aborting.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerateVelocityMaximumFailure</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, tolerate failures to find the radius of the peak in the rotation curve.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOHeated(enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),velocityDispersionApproximate,tolerateEnclosedMassIntegrationFailure,tolerateVelocityDispersionFailure,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,fractionRadiusFinalSmall,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,toleranceRelativePotential,darkMatterProfileDMO_,darkMatterProfileHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"    />
    <objectDestructor name="darkMatterProfileHeating_"/>
    !!]
    return
  end function heatedConstructorParameters

  function heatedConstructorInternal(nonAnalyticSolver,velocityDispersionApproximate,tolerateEnclosedMassIntegrationFailure,tolerateVelocityDispersionFailure,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,fractionRadiusFinalSmall,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,toleranceRelativePotential,darkMatterProfileDMO_,darkMatterProfileHeating_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOHeated} dark matter profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversIsValid
    use :: Error             , only : Error_Report
    implicit none
    type            (darkMatterProfileDMOHeated       )                        :: self
    class           (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterProfileHeatingClass    ), intent(in   ), target :: darkMatterProfileHeating_
    type            (enumerationNonAnalyticSolversType), intent(in   )         :: nonAnalyticSolver
    logical                                            , intent(in   )         :: velocityDispersionApproximate            , tolerateVelocityMaximumFailure                   , &
         &                                                                        toleratePotentialIntegrationFailure      , tolerateEnclosedMassIntegrationFailure           , &
         &                                                                        tolerateVelocityDispersionFailure
    double precision                                   , intent(in   )         :: toleranceRelativeVelocityDispersion      , toleranceRelativeVelocityDispersionMaximum       , &
         &                                                                        fractionRadiusFinalSmall                 , toleranceRelativePotential
    double precision                                   , parameter             :: toleranceAbsolute                  =0.0d0, toleranceRelative                         =1.0d-6
    !![
    <constructorAssign variables="nonAnalyticSolver, velocityDispersionApproximate, tolerateVelocityMaximumFailure, toleratePotentialIntegrationFailure, tolerateEnclosedMassIntegrationFailure, tolerateVelocityDispersionFailure, fractionRadiusFinalSmall, toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, toleranceRelativePotential, *darkMatterProfileDMO_, *darkMatterProfileHeating_"/>
    !!]

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    return
  end function heatedConstructorInternal

  subroutine heatedDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOHeated} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOHeated), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"    />
    <objectDestructor name="self%darkMatterProfileHeating_"/>
    !!]
    return
  end subroutine heatedDestructor

  function heatedGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo          , massTypeDark                , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalHeated, kinematicsDistributionHeated, massDistributionSpherical, massDistributionHeatingClass
    implicit none
    class           (massDistributionClass       ), pointer                 :: massDistribution_
    type            (kinematicsDistributionHeated), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOHeated  ), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    type            (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                                       , intent(in   ), optional :: weightIndex
    class           (massDistributionClass       ), pointer                 :: massDistributionDecorated
    class           (massDistributionHeatingClass), pointer                 :: massDistributionHeating_
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalHeated :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalHeated)
       massDistributionDecorated => self%darkMatterProfileDMO_    %get(node,weightBy,weightIndex)
       massDistributionHeating_  => self%darkMatterProfileHeating_%get(node                     )
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalHeated(                                                                                    &amp;
               &amp;                          nonAnalyticSolver                     =self%nonAnalyticSolver                     , &amp;
	       &amp;                          tolerateVelocityMaximumFailure        =self%tolerateVelocityMaximumFailure        , &amp;
	       &amp;                          tolerateEnclosedMassIntegrationFailure=self%tolerateEnclosedMassIntegrationFailure, &amp;
	       &amp;                          toleratePotentialIntegrationFailure   =self%toleratePotentialIntegrationFailure   , &amp;
	       &amp;                          fractionRadiusFinalSmall              =self%fractionRadiusFinalSmall              , &amp;
	       &amp;                          toleranceRelativePotential            =self%toleranceRelativePotential            , &amp;
               &amp;                          massDistribution_                     =     massDistributionDecorated             , &amp;
               &amp;                          massDistributionHeating_              =     massDistributionHeating_              , &amp;
               &amp;                          componentType                         =     componentTypeDarkHalo                 , &amp;
               &amp;                          massType                              =     massTypeDark                            &amp;
               &amp;                         )
	    </constructor>
	  </referenceConstruct>
	  <objectDestructor name="massDistributionDecorated"/>
	  <objectDestructor name="massDistributionHeating_" />
          !!]
       class default
          call Error_Report('expected a spherical mass distribution'//{introspection:location})
       end select
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionHeated(                                                                                            &amp;
         &amp;                       nonAnalyticSolver                         =self%nonAnalyticSolver                         , &amp; 
         &amp;                       velocityDispersionApproximate             =self%velocityDispersionApproximate             , &amp; 
         &amp;                       toleranceRelativeVelocityDispersion       =self%toleranceRelativeVelocityDispersion       , &amp; 
         &amp;                       toleranceRelativeVelocityDispersionMaximum=self%toleranceRelativeVelocityDispersionMaximum  &amp; 
	 &amp;                      )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function heatedGet
