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
  Implementation of a simple solver for galactic structure (self-gravity of baryons is ignored).
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass

  !![
  <galacticStructureSolver name="galacticStructureSolverSimple">
   <description>
    A galactic structure solver class that determines the sizes of galactic components by assuming that their self-gravity is
    negligible (i.e. that the gravitational potential well is dominated by dark matter) and that, therefore, baryons do not
    modify the dark matter density profile. The radius of a given \gls{component} is then found by solving
    \begin{equation}
     j = \sqrt{\G M_\mathrm{DM}(r) r},
    \end{equation}
    where $j$ is the specific angular momentum of the \gls{component} (at whatever point in the profile is to be solved for),
    $r$ is radius and $M(r)$ is the mass of dark matter within radius $r$. The parameter {\normalfont \ttfamily
    [useFormationHalo]} controls whether the structure of the galaxy will be solved for using the properties of its present
    \gls{node} or those of its \gls{node} at the time of \gls{node} formation (which requires that ``node formation'' has been
    suitably defined and implemented by a component).
   </description>
  </galacticStructureSolver>
  !!]
  type, extends(galacticStructureSolverClass) :: galacticStructureSolverSimple
     !!{
     Implementation of a simple solver for galactic structure (self-gravity of baryons is ignored).
     !!}
     private
     class  (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     logical                                     :: useFormationHalo               , solveForInactiveProperties
   contains
     final     ::             simpleDestructor
     procedure :: solve    => simpleSolve
     procedure :: revert   => simpleRevert
     procedure :: autoHook => simpleAutoHook
  end type galacticStructureSolverSimple

  interface galacticStructureSolverSimple
     !!{
     Constructors for the \refClass{galacticStructureSolverSimple} galactic structure solver class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface galacticStructureSolverSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticStructureSolverSimple} galactic structure solver class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (galacticStructureSolverSimple)                :: self
    type   (inputParameters              ), intent(inout) :: parameters
    class  (darkMatterProfileDMOClass    ), pointer       :: darkMatterProfileDMO_
    logical                                               :: useFormationHalo     , solveForInactiveProperties

    !![
    <inputParameter>
      <name>useFormationHalo</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the ``formation halo'' should be used when solving for the radii of galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>solveForInactiveProperties</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, galactic structure is solved for during evaluation of inactive property integrals. Otherwise, structure is not solved for during this phase---this should only be used if the inactive property integrands \emph{do not} depend on galactic structure.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=galacticStructureSolverSimple(useFormationHalo,solveForInactiveProperties,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(useFormationHalo,solveForInactiveProperties,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{galacticStructureSolverSimple} galactic structure solver class.
    !!}
    implicit none
    type   (galacticStructureSolverSimple)                        :: self
    logical                               , intent(in   )         :: useFormationHalo     , solveForInactiveProperties
    class  (darkMatterProfileDMOClass    ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="useFormationHalo, solveForInactiveProperties, *darkMatterProfileDMO_"/>
    !!]

    return
  end function simpleConstructorInternal

  subroutine simpleAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent  , openMPThreadBindingAtLevel, postEvolveEvent, preDerivativeEvent, &
          &                     satelliteMergerEvent, dependencyDirectionAfter  , dependencyRegEx
    implicit none
    class(galacticStructureSolverSimple), intent(inout) :: self
    type (dependencyRegEx              ), dimension(1)  :: dependencies

    dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^nodeComponent')
    call   preDerivativeEvent%attach(self,simpleSolvePreDeriativeHook,openMPThreadBindingAtLevel,label='structureSolverSimple',dependencies=dependencies)
    call      postEvolveEvent%attach(self,simpleSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverSimple',dependencies=dependencies)
    call satelliteMergerEvent%attach(self,simpleSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverSimple',dependencies=dependencies)
    call   nodePromotionEvent%attach(self,simpleSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverSimple',dependencies=dependencies)
    return
  end subroutine simpleAutoHook

  subroutine simpleDestructor(self)
    !!{
    Destructor for the \refClass{galacticStructureSolverSimple} galactic structure solver class.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent, postEvolveEvent, preDerivativeEvent, satelliteMergerEvent
    implicit none
    type(galacticStructureSolverSimple), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    if (  preDerivativeEvent%isAttached(self,simpleSolvePreDeriativeHook)) call   preDerivativeEvent%detach(self,simpleSolvePreDeriativeHook)
    if (     postEvolveEvent%isAttached(self,simpleSolveHook            )) call      postEvolveEvent%detach(self,simpleSolveHook            )
    if (satelliteMergerEvent%isAttached(self,simpleSolveHook            )) call satelliteMergerEvent%detach(self,simpleSolveHook            )
    if (  nodePromotionEvent%isAttached(self,simpleSolveHook            )) call   nodePromotionEvent%detach(self,simpleSolveHook            )
    return
  end subroutine simpleDestructor

  subroutine simpleSolveHook(self,node)
    !!{
    Hookable wrapper around the solver.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node

    select type (self)
    type is (galacticStructureSolverSimple)
       call self%solve(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine simpleSolveHook

  subroutine simpleSolvePreDeriativeHook(self,node,propertyType)
    !!{
    Hookable wrapper around the solver.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : propertyTypeInactive
    implicit none
    class  (*       ), intent(inout)         :: self
    type   (treeNode), intent(inout), target :: node
    integer          , intent(in   )         :: propertyType
    !$GLC attributes unused :: propertyType

    select type (self)
    type is (galacticStructureSolverSimple)
       if (propertyType /= propertyTypeInactive .or. self%solveForInactiveProperties) call self%solve(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine simpleSolvePreDeriativeHook

  subroutine simpleSolve(self,node,plausibilityOnly)
    !!{
    Solve for the structure of galactic components.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Error              , only : Error_Report
    !![
    <include directive="radiusSolverTask" type="moduleUse">
    !!]
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    !![
    </include>
    <include directive="radiusSolverPlausibility" type="moduleUse">
    !!]
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (galacticStructureSolverSimple), intent(inout)           :: self
    type            (treeNode                     ), intent(inout), target   :: node
    logical                                        , intent(in   ), optional :: plausibilityOnly
    logical                                        , parameter               :: specificAngularMomentumRequired=.true.
    procedure       (solverGet                    ), pointer                 :: radiusGet                             , velocityGet
    procedure       (solverSet                    ), pointer                 :: radiusSet                             , velocitySet
    type            (treeNode                     ), pointer                 :: haloNode
    logical                                                                  :: componentActive
    double precision                                                         :: specificAngularMomentum
    !![
    <optionalArgument name="plausibilityOnly" defaultsTo=".false."/>
    !!]

    ! Check that the galaxy is physical plausible. In this simple solver, we don't act on this.
    node%isPhysicallyPlausible=.true.
    node%isSolvable           =.true.
    !![
    <include directive="radiusSolverPlausibility" type="functionCall" functionType="void">
     <functionArgs>node</functionArgs>
    !!]
    include 'galactic_structure.radius_solver.plausible.inc'
    !![
    </include>
    !!]
    if (node%isPhysicallyPlausible .and. .not.plausibilityOnly_) then
       ! Determine which node to use for halo properties.
       if (self%useFormationHalo) then
          if (.not.associated(node%formationNode)) call Error_Report('no formation node exists'//{introspection:location})
          haloNode => node%formationNode
       else
          haloNode => node
       end if
       ! Solve for each component.
       call Calculations_Reset(node)
       !![
       <include directive="radiusSolverTask" type="functionCall" functionType="void">
        <functionArgs>node,componentActive,specificAngularMomentumRequired,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet</functionArgs>
        <onReturn>if (componentActive) call radiusSolve(node,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet)</onReturn>
       !!]
       include 'galactic_structure.radius_solver.tasks.inc'
       !![
       </include>
       !!]
    end if
    return

  contains

    subroutine radiusSolve(node,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet)
      !!{
      Solve for the equilibrium radius of the given component.
      !!}
      use :: Mass_Distributions, only : massDistributionClass
      implicit none
      type            (treeNode             ), intent(inout)          :: node
      double precision                       , intent(in   )          :: specificAngularMomentum
      procedure       (solverGet            ), intent(in   ), pointer :: radiusGet              , velocityGet
      procedure       (solverSet            ), intent(in   ), pointer :: radiusSet              , velocitySet
      class           (massDistributionClass)               , pointer :: massDistribution_
      double precision                                                :: radius                 , velocity
      !$GLC attributes unused :: radiusGet, velocityGet

      ! Return immediately if the specific angular momentum is zero.
      if (specificAngularMomentum <= 0.0d0) return
      massDistribution_ => self%darkMatterProfileDMO_%get(haloNode)
      ! Find the radius in the dark matter profile with the required specific angular momentum
      radius  =massDistribution_%radiusFromSpecificAngularMomentum(specificAngularMomentum)
      ! Find the velocity at this radius.
      velocity=massDistribution_%rotationCurve                    (radius                 )
      ! Set the component size to new radius and velocity.
      call radiusSet  (node,radius  )
      call velocitySet(node,velocity)
      !![
      <objectDestructor name="massDistribution_"/>
      !!]
      return
    end subroutine radiusSolve

  end subroutine simpleSolve

  subroutine simpleRevert(self,node)
    !!{
    Revert radii for the simple galactic structure solve. Not necessary for this algorithm.
    !!}
    implicit none
    class(galacticStructureSolverSimple), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    return
  end subroutine simpleRevert
