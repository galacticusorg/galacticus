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
  An implementation of dark matter halo profiles with finite resolution (to mimic the effects of resolution in N-body
  simulations for example).
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOFiniteResolutionNFW">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionSphericalFiniteResolutionNFW} objects to mimic a finite
    resolution to an NFW density profile.
   </description>

  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOFiniteResolution) :: darkMatterProfileDMOFiniteResolutionNFW
     !!{
     A dark matter halo profile class implementing finiteResolutionNFW dark matter halos.
     !!}
     private
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::        finiteResolutionNFWDestructor
     procedure :: get => finiteResolutionNFWGet
  end type darkMatterProfileDMOFiniteResolutionNFW

  interface darkMatterProfileDMOFiniteResolutionNFW
     !!{
     Constructors for the \refClass{darkMatterProfileDMOFiniteResolutionNFW} dark matter halo profile class.
     !!}
     module procedure finiteResolutionNFWConstructorParameters
     module procedure finiteResolutionNFWConstructorInternal
  end interface darkMatterProfileDMOFiniteResolutionNFW

contains

  function finiteResolutionNFWConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily finiteResolutionNFW} dark matter halo profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    use :: Input_Parameters  , only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOFiniteResolutionNFW)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    double precision                                                         :: lengthResolution    , massResolution
    type            (varying_string                         )                :: nonAnalyticSolver
    logical                                                                  :: resolutionIsComoving

    !![
    <inputParameter>
      <name>lengthResolution</name>
      <source>parameters</source>
      <description>The resolution length, $\Delta x$.</description>
    </inputParameter>
    <inputParameter>
      <name>massResolution</name>
      <source>parameters</source>
      <description>The resolution mass, $\Delta M$.</description>
    </inputParameter>
    <inputParameter>
      <name>resolutionIsComoving</name>
      <source>parameters</source>
      <description>If true, the resolution length is assumed to be fixed in comoving coordinates, otherwise in physical coordinates.</description>
    </inputParameter>
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !!]
    self=darkMatterProfileDMOFiniteResolutionNFW(lengthResolution,massResolution,resolutionIsComoving,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterHaloScale_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="cosmologyFunctions_"  />
    !!]
    return
  end function finiteResolutionNFWConstructorParameters

  function finiteResolutionNFWConstructorInternal(lengthResolution,massResolution,resolutionIsComoving,nonAnalyticSolver,darkMatterHaloScale_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOFiniteResolutionNFW} dark matter profile class.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    implicit none
    type            (darkMatterProfileDMOFiniteResolutionNFW)                        :: self
    class           (darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    double precision                                         , intent(in   )         :: lengthResolution    , massResolution
    type            (enumerationNonAnalyticSolversType      ), intent(in   )         :: nonAnalyticSolver
    logical                                                  , intent(in   )         :: resolutionIsComoving
    !![
    <constructorAssign variables="lengthResolution, massResolution, resolutionIsComoving, nonAnalyticSolver, *darkMatterHaloScale_, *cosmologyFunctions_"/>
    !!]

    allocate(darkMatterProfileDMONFW :: self%darkMatterProfileDMO_)
    select type (darkMatterProfileDMO_ => self%darkMatterProfileDMO_)
    type is (darkMatterProfileDMONFW)
       !![
       <referenceConstruct owner="self" isResult="yes" object="darkMatterProfileDMO_" nameAssociated="darkMatterProfileDMO_">
         <constructor>
           darkMatterProfileDMONFW(                                                                &amp;
              &amp;                velocityDispersionUseSeriesExpansion=.false.                  , &amp;
              &amp;                darkMatterHaloScale_                =self%darkMatterHaloScale_  &amp;
              &amp;               )
         </constructor>
       </referenceConstruct>
       !!]
    end select    
    return
  end function finiteResolutionNFWConstructorInternal
  
  subroutine finiteResolutionNFWDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOFiniteResolutionNFW} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOFiniteResolutionNFW), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine finiteResolutionNFWDestructor

  function finiteResolutionNFWGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic                          , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                       , massTypeDark                             , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalFiniteResolutionNFW, kinematicsDistributionFiniteResolutionNFW, massDistributionSpherical
    implicit none
    class  (massDistributionClass                    ), pointer                 :: massDistribution_
    type   (kinematicsDistributionFiniteResolutionNFW), pointer                 :: kinematicsDistribution_
    class  (darkMatterProfileDMOFiniteResolutionNFW  ), intent(inout)           :: self
    type   (treeNode                                 ), intent(inout)           :: node
    type   (enumerationWeightByType                  ), intent(in   ), optional :: weightBy
    integer                                           , intent(in   ), optional :: weightIndex
    class  (nodeComponentBasic                       ), pointer                 :: basic
    class  (nodeComponentDarkMatterProfile           ), pointer                 :: darkMatterProfile
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalFiniteResolutionNFW :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalFiniteResolutionNFW)
       basic => node%basic()
       darkMatterProfile => node%darkMatterProfile()
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionSphericalFiniteResolutionNFW(                                                                                         &amp;
	   &amp;                                        lengthResolution =self                                  %lengthResolutionPhysical(node), &amp;
	   &amp;                                        radiusScale      =darkMatterProfile                     %scale                   (    ), &amp;
	   &amp;                                        radiusVirial     =self             %darkMatterHaloScale_%radiusVirial            (node), &amp;
	   &amp;                                        mass             =basic                                 %mass                    (    ), &amp;
           &amp;                                        componentType    =                                       componentTypeDarkHalo         , &amp;
           &amp;                                        massType         =                                       massTypeDark                    &amp;
           &amp;                                       )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionFiniteResolutionNFW( &amp;
	 &amp;                                   )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function finiteResolutionNFWGet
