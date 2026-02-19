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
  An implementation of truncated dark matter halo profiles.
  !!}

  use :: Mass_Distributions     , only : enumerationNonAnalyticSolversType
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOTruncated">
    <description>
      Truncated dark matter halo profiles are built via the \refClass{massDistributionSphericalTruncated} class.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOTruncated
     !!{
     A dark matter halo profile class implementing truncated dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_           => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_            => null()
     double precision                                             :: radiusFractionalTruncateMinimum          , radiusFractionalTruncateMaximum
     type            (enumerationNonAnalyticSolversType)          :: nonAnalyticSolver
   contains
     final     ::        truncatedDestructor
     procedure :: get => truncatedGet
  end type darkMatterProfileDMOTruncated

  interface darkMatterProfileDMOTruncated
     !!{
     Constructors for the \refClass{darkMatterProfileDMOTruncated} dark matter halo profile class.
     !!}
     module procedure truncatedConstructorParameters
     module procedure truncatedConstructorInternal
  end interface darkMatterProfileDMOTruncated

contains

  function truncatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOTruncated} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters  , only : inputParameters
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    implicit none
    type            (darkMatterProfileDMOTruncated)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass    ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass     ), pointer       :: darkMatterHaloScale_
    type            (varying_string               )                :: nonAnalyticSolver
    double precision                                               :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusFractionalTruncateMinimum</name>
      <defaultValue>2.0d0</defaultValue>
      <source>parameters</source>
      <description>The minimum radius (in units of the virial radius) to begin truncating the density profile.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusFractionalTruncateMaximum</name>
      <defaultValue>4.0d0</defaultValue>
      <source>parameters</source>
      <description>The maximum radius (in units of the virial radius) to finish truncating the density profile.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO"   name="darkMatterProfileDMO_"   source="parameters"/>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=darkMatterProfileDMOTruncated(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"  />
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function truncatedConstructorParameters

  function truncatedConstructorInternal(radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOTruncated} dark matter profile class.
    !!}
    use :: Error             , only : Error_Report
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversIsValid
    implicit none
    type            (darkMatterProfileDMOTruncated    )                        :: self
    class           (darkMatterProfileDMOClass        ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                   , intent(in   )         :: radiusFractionalTruncateMinimum, radiusFractionalTruncateMaximum
    type            (enumerationNonAnalyticSolversType), intent(in   )         :: nonAnalyticSolver
    !![
    <constructorAssign variables="radiusFractionalTruncateMinimum,radiusFractionalTruncateMaximum,nonAnalyticSolver,*darkMatterProfileDMO_,*darkMatterHaloScale_"/>
    !!]

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    return
  end function truncatedConstructorInternal

  subroutine truncatedDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOTruncated} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOTruncated), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine truncatedDestructor

  function truncatedGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo             , massTypeDark                   , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalTruncated, kinematicsDistributionTruncated, massDistributionSpherical
    implicit none
    class           (massDistributionClass          ), pointer                 :: massDistribution_
    type            (kinematicsDistributionTruncated), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOTruncated  ), intent(inout)           :: self
    type            (treeNode                       ), intent(inout)           :: node
    type            (enumerationWeightByType        ), intent(in   ), optional :: weightBy
    integer                                          , intent(in   ), optional :: weightIndex
    double precision                                                           :: radiusVirial
    class           (massDistributionClass          ), pointer                 :: massDistributionDecorated
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalTruncated :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalTruncated)
       radiusVirial              =  self%darkMatterHaloScale_ %radiusVirial(node                     )
       massDistributionDecorated => self%darkMatterProfileDMO_%get         (node,weightBy,weightIndex)
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalTruncated(                                                                         &amp;
	      &amp;                              radiusTruncateMinimum=self%radiusFractionalTruncateMinimum*radiusVirial, &amp;
	      &amp;                              radiusTruncateMaximum=self%radiusFractionalTruncateMaximum*radiusVirial, &amp;
	      &amp;                              nonAnalyticSolver    =self%nonAnalyticSolver                           , &amp;
	      &amp;                              massDistribution_    =     massDistributionDecorated                   , &amp;
              &amp;                              componentType        =     componentTypeDarkHalo                       , &amp;
              &amp;                              massType             =     massTypeDark                                  &amp;
              &amp;                             )
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
        kinematicsDistributionTruncated( &amp;
	 &amp;                         )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function truncatedGet
