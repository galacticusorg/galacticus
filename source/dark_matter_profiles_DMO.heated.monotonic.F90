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
  An implementation of heated dark matter halo profiles based on the energy ordering of shells.
  !!}

  use :: Mass_Distributions, only : enumerationNonAnalyticSolversType

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOHeatedMonotonic">
   <description>
     A dark matter profile DMO class in which builds \refClass{massDistributionSphericalHeatedMonotonic} objects to account for
     heating of some other dark matter profile.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOHeatedMonotonic
     !!{
     A dark matter halo profile class implementing heated dark matter halos.
     !!}
     private
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_                => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_               => null()
     class           (darkMatterProfileHeatingClass    ), pointer :: darkMatterProfileHeating_           => null()
     type            (enumerationNonAnalyticSolversType)          :: nonAnalyticSolver
     double precision                                             :: toleranceRelativeVelocityDispersion          , toleranceRelativeVelocityDispersionMaximum
   contains
     final     ::        heatedMonotonicDestructor
     procedure :: get => heatedMonotonicGet
  end type darkMatterProfileDMOHeatedMonotonic

  interface darkMatterProfileDMOHeatedMonotonic
     !!{
     Constructors for the \refClass{darkMatterProfileDMOHeatedMonotonic} dark matter halo profile class.
     !!}
     module procedure heatedMonotonicConstructorParameters
     module procedure heatedMonotonicConstructorInternal
  end interface darkMatterProfileDMOHeatedMonotonic

contains

  function heatedMonotonicConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily heatedMonotonic} dark matter halo profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    use :: Input_Parameters  , only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOHeatedMonotonic)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass          ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileHeatingClass      ), pointer       :: darkMatterProfileHeating_
    type            (varying_string                     )                :: nonAnalyticSolver
    double precision                                                     :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
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
    <objectBuilder class="darkMatterProfileDMO"     name="darkMatterProfileDMO_"     source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOHeatedMonotonic(enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"    />
    <objectDestructor name="darkMatterHaloScale_"     />
    <objectDestructor name="darkMatterProfileHeating_"/>
    !!]
    return
  end function heatedMonotonicConstructorParameters

  function heatedMonotonicConstructorInternal(nonAnalyticSolver,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,darkMatterProfileDMO_,darkMatterHaloScale_,darkMatterProfileHeating_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOHeatedMonotonic} dark matter profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversIsValid
    use :: Error             , only : Error_Report
    implicit none
    type            (darkMatterProfileDMOHeatedMonotonic)                        :: self
    class           (darkMatterProfileDMOClass          ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileHeatingClass      ), intent(in   ), target :: darkMatterProfileHeating_
    type            (enumerationNonAnalyticSolversType  ), intent(in   )         :: nonAnalyticSolver
    double precision                                     , intent(in   )         :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="nonAnalyticSolver, toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, *darkMatterProfileDMO_, *darkMatterHaloScale_, *darkMatterProfileHeating_"/>
    !!]

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})   
    return
  end function heatedMonotonicConstructorInternal

  subroutine heatedMonotonicDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOHeatedMonotonic} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOHeatedMonotonic), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"     />
    <objectDestructor name="self%darkMatterHaloScale_"      />
    <objectDestructor name="self%darkMatterProfileHeating_" />
    !!]
    return
  end subroutine heatedMonotonicDestructor

  function heatedMonotonicGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                   , massTypeDark                       , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalHeatedMonotonic, kinematicsDistributionCollisionless, massDistributionSpherical, massDistributionHeatingClass
    implicit none
    class           (massDistributionClass              ), pointer                 :: massDistribution_
    type            (kinematicsDistributionCollisionless), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOHeatedMonotonic), intent(inout)           :: self
    type            (treeNode                           ), intent(inout)           :: node
    type            (enumerationWeightByType            ), intent(in   ), optional :: weightBy
    integer                                              , intent(in   ), optional :: weightIndex
    class           (massDistributionClass              ), pointer                 :: massDistributionDecorated
    class           (massDistributionHeatingClass       ), pointer                 :: massDistributionHeating_
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalHeatedMonotonic :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalHeatedMonotonic)
       massDistributionDecorated => self%darkMatterProfileDMO_    %get(node,weightBy,weightIndex)
       massDistributionHeating_  => self%darkMatterProfileHeating_%get(node                     )
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalHeatedMonotonic(                                                                                    &amp;
               &amp;                                   radiusVirial            =self%darkMatterHaloScale_%radiusVirial             (node), &amp;
               &amp;                                   nonAnalyticSolver       =self                     %nonAnalyticSolver              , &amp;
               &amp;                                   massDistribution_       =                          massDistributionDecorated      , &amp;
               &amp;                                   massDistributionHeating_=                          massDistributionHeating_       , &amp;
               &amp;                                   componentType           =                          componentTypeDarkHalo          , &amp;
               &amp;                                   massType                =                          massTypeDark                     &amp;
               &amp;                                  )
	    </constructor>
	  </referenceConstruct>
          !!]
       class default
          call Error_Report('expected a spherical mass distribution'//{introspection:location})
       end select
       !![
       <objectDestructor name="massDistributionDecorated"/>
       <objectDestructor name="massDistributionHeating_" />
       !!]
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
  end function heatedMonotonicGet
