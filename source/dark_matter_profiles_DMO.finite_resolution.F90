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
  An implementation of dark matter halo profiles with finite resolution (to mimic the effects of resolution in N-body
  simulations for example).
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Mass_Distributions , only : enumerationNonAnalyticSolversType

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOFiniteResolution">
    <description>
      A dark matter profile DMO class which builds \refClass{massDistributionSphericalFiniteResolution} objects to mimic the
      effects of finite resolution in an N-body simulation.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOFiniteResolution
     !!{
     A dark matter halo profile class implementing finiteResolution dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_    => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_      => null()
     type            (enumerationNonAnalyticSolversType)          :: nonAnalyticSolver
     logical                                                      :: resolutionIsComoving
     double precision                                             :: lengthResolution                  , massResolution, &
          &                                                          lengthResolutionPrevious
     integer         (kind_int8                        )          :: lastUniqueID
   contains
     !![
     <methods>
       <method description="Return the resolution length in physical units." method="lengthResolutionPhysical"/>
       <method description="Reset memoized calculations."                    method="calculationReset"        />
     </methods>
     !!]
     final     ::                             finiteResolutionDestructor
     procedure :: autoHook                 => finiteResolutionAutoHook
     procedure :: calculationReset         => finiteResolutionCalculationReset
     procedure :: get                      => finiteResolutionGet    
     procedure :: lengthResolutionPhysical => finiteResolutionLengthResolutionPhysical
  end type darkMatterProfileDMOFiniteResolution

 interface darkMatterProfileDMOFiniteResolution
     !!{
     Constructors for the \refClass{darkMatterProfileDMOFiniteResolution} dark matter halo profile class.
     !!}
     module procedure finiteResolutionConstructorParameters
     module procedure finiteResolutionConstructorInternal
  end interface darkMatterProfileDMOFiniteResolution

  double precision, parameter :: radiusLengthResolutionRatioMaximum=100.0d0

contains

  function finiteResolutionConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily finiteResolution} dark matter halo profile class.
    !!}
    use :: Input_Parameters  , only : inputParameter                     , inputParameters
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    implicit none
    type            (darkMatterProfileDMOFiniteResolution)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass           ), pointer       :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass             ), pointer       :: cosmologyFunctions_
    double precision                                                      :: lengthResolution       , massResolution
    type            (varying_string                      )                :: nonAnalyticSolver
    logical                                                               :: resolutionIsComoving

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
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"   source="parameters"/>
    !!]
    self=darkMatterProfileDMOFiniteResolution(lengthResolution,massResolution,resolutionIsComoving,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterProfileDMO_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="cosmologyFunctions_"  />
    !!]
    return
  end function finiteResolutionConstructorParameters

  function finiteResolutionConstructorInternal(lengthResolution,massResolution,resolutionIsComoving,nonAnalyticSolver,darkMatterProfileDMO_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOFiniteResolution} dark matter profile class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (darkMatterProfileDMOFiniteResolution)                        :: self
    class           (darkMatterProfileDMOClass           ), intent(in   ), target :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass             ), intent(in   ), target :: cosmologyFunctions_
    double precision                                      , intent(in   )         :: lengthResolution       , massResolution
    type            (enumerationNonAnalyticSolversType   ), intent(in   )         :: nonAnalyticSolver
    logical                                               , intent(in   )         :: resolutionIsComoving
    !![
    <constructorAssign variables="lengthResolution, massResolution, resolutionIsComoving, nonAnalyticSolver, *darkMatterProfileDMO_, *cosmologyFunctions_"/>
    !!]
    
    self%lastUniqueID=-huge(1_kind_int8)
    return
  end function finiteResolutionConstructorInternal

  subroutine finiteResolutionAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    
    call calculationResetEvent%attach(self,finiteResolutionCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOFiniteResolution')
    return
  end subroutine finiteResolutionAutoHook
  
  subroutine finiteResolutionDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOFiniteResolution} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOFiniteResolution), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%cosmologyFunctions_"  />
    !!]
    if (calculationResetEvent%isAttached(self,finiteResolutionCalculationReset)) call calculationResetEvent%detach(self,finiteResolutionCalculationReset)
    return
  end subroutine finiteResolutionDestructor

  subroutine finiteResolutionCalculationReset(self,node,uniqueID)
    !!{
    Reset the dark matter profile calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type   (treeNode                            ), intent(inout) :: node
    integer(kind_int8                           ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%lastUniqueID            =uniqueID
    self%lengthResolutionPrevious=-huge(0.0d0)
    return
  end subroutine finiteResolutionCalculationReset

  function finiteResolutionGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                    , massTypeDark                       , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalFiniteResolution, kinematicsDistributionCollisionless, massDistributionSpherical
    implicit none
    class           (massDistributionClass               ), pointer                 :: massDistribution_
    type            (kinematicsDistributionCollisionless ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOFiniteResolution), intent(inout)           :: self
    type            (treeNode                            ), intent(inout)           :: node
    type            (enumerationWeightByType             ), intent(in   ), optional :: weightBy
    integer                                               , intent(in   ), optional :: weightIndex
    class           (massDistributionClass               ), pointer                 :: massDistributionDecorated
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalFiniteResolution :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalFiniteResolution)
       massDistributionDecorated => self%darkMatterProfileDMO_%get(node,weightBy,weightIndex)
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalFiniteResolution(                                                        &amp;
	      &amp;                                     lengthResolution =self%lengthResolutionPhysical (node), &amp;
	      &amp;                                     nonAnalyticSolver=self%nonAnalyticSolver              , &amp;
	      &amp;                                     massDistribution_=     massDistributionDecorated      , &amp;
              &amp;                                     componentType    =     componentTypeDarkHalo          , &amp;
              &amp;                                     massType         =     massTypeDark                     &amp;
              &amp;                                    )
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
  end function finiteResolutionGet

  double precision function finiteResolutionLengthResolutionPhysical(self,node)
    !!{
    Return the resolution length in physical units.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(darkMatterProfileDMOFiniteResolution), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    class(nodeComponentBasic                  ), pointer       :: basic
    class(massDistributionClass               ), pointer       :: massDistribution_

    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    if (self%lengthResolutionPrevious < 0.0d0) then
       self%lengthResolutionPrevious=self%lengthResolution
       if (self%resolutionIsComoving) then
          basic                            =>                                           node %basic()
          self%lengthResolutionPrevious    =  +self%lengthResolutionPrevious                           &
               &                              *self%cosmologyFunctions_%expansionFactor(basic%time ())
       end if
       if (self%massResolution > 0.0d0) then
          massDistribution_             => self%darkMatterProfileDMO_%get(node)
          self%lengthResolutionPrevious =  max(                                                                 &
               &                               self             %lengthResolutionPrevious                     , &
               &                               massDistribution_%radiusEnclosingMass     (self%massResolution)  &
               &                              )
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       end if
    end if
    finiteResolutionLengthResolutionPhysical=self%lengthResolutionPrevious
    return
  end function finiteResolutionLengthResolutionPhysical
