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
  An accelerator class for dark matter halo profiles.
  !!}

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOAccelerator">
   <description>
    An accelerator wrapper for dark-matter-only profile classes that caches recently computed profile quantities (enclosed mass, density, potential, etc.) and returns interpolated estimates when the requested radius is close to a cached value, significantly speeding up calculations that repeatedly query the profile at nearby radii.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOAccelerator
     !!{
     An accelerator class for the dark matter halo profile class.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer      :: darkMatterProfileDMO_ => null()
     double precision                                          :: toleranceRelative              , factorRadiusMaximum
   contains
     final     ::        acceleratorDestructor
     procedure :: get => acceleratorGet
  end type darkMatterProfileDMOAccelerator

  interface darkMatterProfileDMOAccelerator
     !!{
     Constructors for the \refClass{darkMatterProfileDMOAccelerator} dark matter halo profile class.
     !!}
     module procedure acceleratorConstructorParameters
     module procedure acceleratorConstructorInternal
  end interface darkMatterProfileDMOAccelerator

contains

  function acceleratorConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOAccelerator} dark matter halo profile class which takes a parameter set as input.
    !!}
    implicit none
    type   (darkMatterProfileDMOAccelerator)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    class  (darkMatterProfileDMOClass      ), pointer       :: darkMatterProfileDMO_
    double precision                                        :: toleranceRelative    , factorRadiusMaximum

    !![
    <inputParameter>
      <name>toleranceRelative</name>
      <defaultValue>1.0d-2</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance within which a cached (accelerated) profile estimate is accepted; if the requested radius differs from the cached value by more than this fractional amount, the full profile calculation is recomputed.</description>
    </inputParameter>
    <inputParameter>
      <name>factorRadiusMaximum</name>
      <defaultValue>3.0d0</defaultValue>
      <source>parameters</source>
      <description>The maximum radial extrapolation factor allowed when using cached profile values for the dark-matter-only profile; if the requested radius exceeds the cached radius by more than this factor, the full calculation is recomputed.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" parameterName="darkMatterProfileDMO" source="parameters"/>
    !!]
    self=darkMatterProfileDMOAccelerator(toleranceRelative,factorRadiusMaximum,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function acceleratorConstructorParameters

  function acceleratorConstructorInternal(toleranceRelative,factorRadiusMaximum,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOAccelerator} dark matter profile class.
    !!}
    implicit none
    type            (darkMatterProfileDMOAccelerator)                        :: self
    class           (darkMatterProfileDMOClass      ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                 , intent(in)            :: toleranceRelative    , factorRadiusMaximum
    !![
    <constructorAssign variables="toleranceRelative, factorRadiusMaximum, *darkMatterProfileDMO_"/>
    !!]
    
    return
  end function acceleratorConstructorInternal

  subroutine acceleratorDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOAccelerator} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOAccelerator), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine acceleratorDestructor

  function acceleratorGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given \mono{node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo               , massTypeDark                       , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalAccelerator, kinematicsDistributionCollisionless, massDistributionSpherical, nonAnalyticSolversNumerical
    implicit none
    class           (massDistributionClass              ), pointer                 :: massDistribution_
    type            (kinematicsDistributionCollisionless), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOAccelerator    ), intent(inout)           :: self
    type            (treeNode                           ), intent(inout)           :: node
    type            (enumerationWeightByType            ), intent(in   ), optional :: weightBy
    integer                                              , intent(in   ), optional :: weightIndex
    class           (massDistributionClass              ), pointer                 :: massDistributionDecorated
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalAccelerator :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalAccelerator)
       massDistributionDecorated => self%darkMatterProfileDMO_%get(node,weightBy,weightIndex)
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalAccelerator(                                                      &amp;
	      &amp;                                toleranceRelative  =self%toleranceRelative          , &amp;
	      &amp;                                factorRadiusMaximum=self%factorRadiusMaximum        , &amp;
	      &amp;                                massDistribution_  =     massDistributionDecorated  , &amp;
              &amp;                                nonAnalyticSolver  =     nonAnalyticSolversNumerical, &amp;
              &amp;                                componentType      =     componentTypeDarkHalo      , &amp;
              &amp;                                massType           =     massTypeDark                 &amp;
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
  end function acceleratorGet
