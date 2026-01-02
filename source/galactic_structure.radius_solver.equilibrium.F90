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
  Implementation of an ``equilibrium'' solver for galactic structure.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass

  !![
  <galacticStructureSolver name="galacticStructureSolverEquilibrium">
   <description>An ``equilibrium'' solver for galactic structure.</description>
  </galacticStructureSolver>
  !!]
  type, extends(galacticStructureSolverClass) :: galacticStructureSolverEquilibrium
     !!{
     Implementation of an ``equilibrium'' solver for galactic structure.
     !!}
     private
     logical                                              :: includeBaryonGravity                , useFormationHalo         , &
          &                                                  solveForInactiveProperties          , convergenceFailureIsFatal
     double precision                                     :: solutionTolerance
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_       => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_      => null()
   contains
     final     ::             equilibriumDestructor
     procedure :: solve    => equilibriumSolve
     procedure :: revert   => equilibriumRevert
     procedure :: autoHook => equilibriumAutoHook
  end type galacticStructureSolverEquilibrium

  interface galacticStructureSolverEquilibrium
     !!{
     Constructors for the \refClass{galacticStructureSolverEquilibrium} galactic structure solver class.
     !!}
     module procedure equilibriumConstructorParameters
     module procedure equilibriumConstructorInternal
  end interface galacticStructureSolverEquilibrium

  ! Module-scope variables used to communicate current state of radius solver.
  integer                                               :: countComponentsActive        , countIterations
  double precision                                      :: fitMeasure
  type            (treeNode), pointer                   :: node_
  double precision          , allocatable, dimension(:) :: radiusStored                , velocityStored
  logical                                               :: revertStructure     =.false.
  !$omp threadprivate(countIterations,countComponentsActive,fitMeasure,node_,radiusStored,velocityStored,revertStructure)

contains

  function equilibriumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticStructureSolverEquilibrium} galactic structure solver class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticStructureSolverEquilibrium)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
    logical                                                             :: useFormationHalo          , includeBaryonGravity     , &
         &                                                                 solveForInactiveProperties, convergenceFailureIsFatal
    double precision                                                    :: solutionTolerance

    !![
    <inputParameter>
      <name>convergenceFailureIsFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, failure to achieve convergence in radii results in a fatal error.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>includeBaryonGravity</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether or not gravity from baryons is included when solving for sizes of galactic components.</description>
      <source>parameters</source>
    </inputParameter>
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
    <inputParameter>
      <name>solutionTolerance</name>
      <defaultValue>1.0d-2</defaultValue>
      <description>Maximum allowed mean fractional error in the radii of all components when seeking equilibrium solutions for galactic structure.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=galacticStructureSolverEquilibrium(convergenceFailureIsFatal,useFormationHalo,includeBaryonGravity,solutionTolerance,solveForInactiveProperties,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function equilibriumConstructorParameters

  function equilibriumConstructorInternal(convergenceFailureIsFatal,useFormationHalo,includeBaryonGravity,solutionTolerance,solveForInactiveProperties,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the \refClass{galacticStructureSolverEquilibrium} galactic structure solver class.
    !!}
    implicit none
    type            (galacticStructureSolverEquilibrium)                        :: self
    logical                                             , intent(in   )         :: useFormationHalo          , includeBaryonGravity     , &
         &                                                                         solveForInactiveProperties, convergenceFailureIsFatal
    double precision                                    , intent(in   )         :: solutionTolerance
    class           (darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="convergenceFailureIsFatal, useFormationHalo, includeBaryonGravity, solutionTolerance, solveForInactiveProperties, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function equilibriumConstructorInternal

  subroutine equilibriumAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : dependencyDirectionAfter, dependencyRegEx   , nodePromotionEvent  , openMPThreadBindingAtLevel, &
          &                     postEvolveEvent         , preDerivativeEvent, satelliteMergerEvent
    implicit none
    class(galacticStructureSolverEquilibrium), intent(inout) :: self
    type (dependencyRegEx                   ), dimension(1)  :: dependencies

    dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^nodeComponent')
    call   preDerivativeEvent%attach(self,equilibriumSolvePreDeriativeHook,openMPThreadBindingAtLevel,label='structureSolverEquilibrium'                          )
    call      postEvolveEvent%attach(self,equilibriumSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverEquilibrium',dependencies=dependencies)
    call satelliteMergerEvent%attach(self,equilibriumSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverEquilibrium',dependencies=dependencies)
    call   nodePromotionEvent%attach(self,equilibriumSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverEquilibrium',dependencies=dependencies)
    return
  end subroutine equilibriumAutoHook

  subroutine equilibriumDestructor(self)
    !!{
    Destructor for the \refClass{galacticStructureSolverEquilibrium} galactic structure solver class.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent, postEvolveEvent, preDerivativeEvent, satelliteMergerEvent
    implicit none
    type(galacticStructureSolverEquilibrium), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    if (  preDerivativeEvent%isAttached(self,equilibriumSolvePreDeriativeHook)) call   preDerivativeEvent%detach(self,equilibriumSolvePreDeriativeHook)
    if (     postEvolveEvent%isAttached(self,equilibriumSolveHook            )) call      postEvolveEvent%detach(self,equilibriumSolveHook            )
    if (satelliteMergerEvent%isAttached(self,equilibriumSolveHook            )) call satelliteMergerEvent%detach(self,equilibriumSolveHook            )
    if (  nodePromotionEvent%isAttached(self,equilibriumSolveHook            )) call   nodePromotionEvent%detach(self,equilibriumSolveHook            )
    return
  end subroutine equilibriumDestructor

  subroutine equilibriumSolveHook(self,node)
    !!{
    Hookable wrapper around the solver.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node

    select type (self)
    type is (galacticStructureSolverEquilibrium)
       call self%solve(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine equilibriumSolveHook

  subroutine equilibriumSolvePreDeriativeHook(self,node,propertyType)
    !!{
    Hookable wrapper around the solver for pre-derivative events.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : propertyTypeInactive, treeNode
    implicit none
    class  (*       ), intent(inout)         :: self
    type   (treeNode), intent(inout), target :: node
    integer          , intent(in   )         :: propertyType

    select type (self)
    type is (galacticStructureSolverEquilibrium)
       call self%solve(node,plausibilityOnly=propertyType == propertyTypeInactive .and. .not.self%solveForInactiveProperties)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine equilibriumSolvePreDeriativeHook

  subroutine equilibriumSolve(self,node,plausibilityOnly)
    !!{
    Solve for the structure of galactic components.
    !!}
    use :: Calculations_Resets       , only : Calculations_Reset
    use :: Display                   , only : displayMessage
    use :: Error                     , only : Error_Report                , Warn
    use :: Galactic_Structure_Options, only : enumerationComponentTypeType
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    implicit none
    class           (galacticStructureSolverEquilibrium), intent(inout)           :: self
    type            (treeNode                          ), intent(inout), target   :: node
    logical                                             , intent(in   ), optional :: plausibilityOnly
    logical                                             , parameter               :: specificAngularMomentumRequired=.true.
    integer                                             , parameter               :: iterationMaximum               =  100
    procedure       (solverGet                         ), pointer                 :: radiusGet                             , velocityGet
    procedure       (solverSet                         ), pointer                 :: radiusSet                             , velocitySet
    logical                                                                       :: componentActive
    double precision                                                              :: specificAngularMomentum
    type            (enumerationComponentTypeType      )                          :: component
    !![
    <optionalArgument name="plausibilityOnly" defaultsTo=".false."/>
    !!]
    
    ! Check that the galaxy is physical plausible. In this equilibrium solver, we don't act on this.
    node%isPhysicallyPlausible=.true.
    node%isSolvable           =.true.
    include 'galactic_structure.radius_solver.plausible.inc'
    if (node%isPhysicallyPlausible .and. .not.plausibilityOnly_) then
       ! Initialize the solver state.
       countIterations=0
       fitMeasure    =2.0d0*self%solutionTolerance
       ! Determine which node to use for halo properties.
       if (self%useFormationHalo) then
          if (.not.associated(node%formationNode)) call Error_Report('no formation node exists'//{introspection:location})
          node_ => node%formationNode
       else
          node_ => node
       end if
       ! Begin iteration to find a converged solution.
       do while (countIterations <= 1 .or. ( fitMeasure > self%solutionTolerance .and. countIterations < iterationMaximum ) )
          countIterations      =countIterations+1
          countComponentsActive=0
          if (countIterations > 1) fitMeasure=0.0d0
          include 'galactic_structure.radius_solver.tasks.inc'
          ! Check that we have some active components.
          if (countComponentsActive == 0) then
             fitMeasure=0.0d0
             exit
          else
             ! Normalize the fit measure by the number of active components.
             fitMeasure=fitMeasure/dble(countComponentsActive)
          end if
       end do
       ! Check that we found a converged solution.
       if (fitMeasure > self%solutionTolerance) then
          if (self%convergenceFailureIsFatal) then
             call displayMessage('dumping node for which radii are currently being sought')
             call node%serializeASCII()
             call Error_Report('failed to find converged solution'//{introspection:location})
          else
             call Warn('failed to find converged solution for galactic structure radii')
          end if
       end if
    end if
    ! Unset structure reversion flag.
    revertStructure=.false.
    return

  contains

    subroutine radiusSolve(node,component,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet)
      !!{
      Solve for the equilibrium radius of the given component.
      !!}
      use :: Display                         , only : displayVerbosity              , displayVerbositySet, verbosityLevelStandard
      use :: Galactic_Structure_Options      , only : massTypeBaryonic              , radiusLarge        , massTypeDark          , componentTypeDarkHalo
      use :: Mass_Distributions              , only : massDistributionClass
      use :: Error                           , only : Error_Report
      use :: ISO_Varying_String              , only : varying_string
      use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
      use :: String_Handling                 , only : operator(//)
      implicit none
      type            (treeNode                    ), intent(inout)                     :: node
      type            (enumerationComponentTypeType), intent(in   )                     :: component
      double precision                              , intent(in   )                     :: specificAngularMomentum
      procedure       (solverGet                   ), intent(in   ) , pointer           :: radiusGet                         , velocityGet
      procedure       (solverSet                   ), intent(in   ) , pointer           :: radiusSet                         , velocitySet
      double precision                              , dimension(:,:), allocatable, save :: radiusHistory
      !$omp threadprivate(radiusHistory)
      double precision                              , dimension(:,:), allocatable       :: radiusHistoryTemporary
      double precision                              , dimension(:  ), allocatable       :: radiusStoredTmp                   , velocityStoredTmp
      class           (massDistributionClass       ), pointer                           :: massDistribution_
      integer                                       , parameter                         :: storeIncrement                 =10
      integer                                       , parameter                         :: iterationsForBisectionMinimum  =10
      integer                                       , parameter                         :: activeComponentMaximumIncrement= 2
      integer                                                                           :: activeComponentMaximumCurrent
      character       (len=14                      )                                    :: label
      type            (varying_string              ), save                              :: message
      !$omp threadprivate(message)
      double precision                                                                  :: baryonicVelocitySquared           , darkMatterMassFinal, &
           &                                                                               darkMatterVelocitySquared         , velocity           , &
           &                                                                               radius                            , radiusNew          , &
           &                                                                               specificAngularMomentumMaximum
      !$GLC attributes unused :: component

      ! Count the number of active components.
      countComponentsActive=countComponentsActive+1
      if (countIterations == 1) then
         ! If structure is to be reverted, do so now.
         if (revertStructure.and.allocated(radiusStored)) then
            call          radiusSet(node,  radiusStored(countComponentsActive))
            call        velocitySet(node,velocityStored(countComponentsActive))
            call Calculations_Reset(node                                      )
         end if
         ! On first iteration, see if we have a previous radius set for this component.
         radius=radiusGet(node)
         if (radius <= 0.0d0) then
            ! No previous radius was set, so make a simple estimate of sizes of all components ignoring equilibrium contraction and self-gravity.
            massDistribution_ => self%darkMatterProfileDMO_%get(node)
            ! First check that there is a solution within a reasonable radius.
            specificAngularMomentumMaximum=+massDistribution_%rotationCurve(radiusLarge) &
                 &                         *                                radiusLarge
            if (specificAngularMomentumMaximum < specificAngularMomentum) then
               ! No solution exists even within a very large radius. Use a simple estimate of the virial radius.
               radius=self%darkMatterHaloScale_%radiusVirial                      (node_                        )
            else
               ! Find the radius in the dark matter profile with the required specific angular momentum.
               radius=massDistribution_%radiusFromSpecificAngularMomentum(specificAngularMomentum)
            end if
            ! Find the velocity at this radius.
            velocity=massDistribution_%rotationCurve(radius)
            !![
            <objectDestructor name="massDistribution_"/>
            !!]
         else
            ! A previous radius was set, so use it, and the previous circular velocity, as the initial guess.
            velocity=velocityGet(node)
         end if
         ! If structure is not being reverted, store the new values of radius and velocity.
         if (.not.revertStructure .and. countIterations == 1) then
            ! Store these quantities.
            if (.not.allocated(radiusStored)) then
               allocate(  radiusStored(storeIncrement))
               allocate(velocityStored(storeIncrement))
               radiusStored  =0.0d0
               velocityStored=0.0d0
            else if (countComponentsActive > size(radiusStored)) then
               call move_alloc(  radiusStored,       radiusStoredTmp                )
               call move_alloc(velocityStored,     velocityStoredTmp                )
               allocate       (  radiusStored(size(  radiusStoredTmp)+storeIncrement))
               allocate       (velocityStored(size(velocityStoredTmp)+storeIncrement))
               radiusStored  (                        1:size(  radiusStoredTmp))=  radiusStoredTmp
               velocityStored(                        1:size(velocityStoredTmp))=velocityStoredTmp
               radiusStored  (size(  radiusStoredTmp)+1:size(  radiusStored   ))=0.0d0
               velocityStored(size(velocityStoredTmp)+1:size(velocityStored   ))=0.0d0
               deallocate(  radiusStoredTmp)
               deallocate(velocityStoredTmp)
            end if
            radiusStored  (countComponentsActive)=radius
            velocityStored(countComponentsActive)=velocity
         end if
      else
         ! On subsequent iterations do the full calculation providing component has non-zero specific angular momentum.
         if (specificAngularMomentum <= 0.0d0) return
         ! Get current radius of the component.
         radius=radiusGet(node)
         ! Find the enclosed mass in the dark matter halo.
         massDistribution_   => node             %massDistribution    (componentTypeDarkHalo,massTypeDark)
         darkMatterMassFinal =  massDistribution_%massEnclosedBySphere(radius                            )
         !![
	 <objectDestructor name="massDistribution_"/>
         !!]
         ! Compute dark matter contribution to rotation curve.
         darkMatterVelocitySquared=gravitationalConstant_internal*darkMatterMassFinal/radius
         ! Compute baryonic contribution to rotation curve.
         if (self%includeBaryonGravity) then
            massDistribution_       => node             %massDistribution(massType=massTypeBaryonic)
            baryonicVelocitySquared =  massDistribution_%rotationCurve   (         radius          )**2
            !![
	    <objectDestructor name="massDistribution_"/>
	    !!]
         else
            baryonicVelocitySquared=0.0d0
         end if
         ! Compute new estimate of velocity.
         velocity=sqrt(darkMatterVelocitySquared+baryonicVelocitySquared)
         ! Compute new estimate of radius.
         if (radius > 0.0d0) then
            radiusNew=sqrt(specificAngularMomentum/velocity*radius)
         else
            radiusNew=     specificAngularMomentum/velocity
         endif
         ! Ensure that the radius history array is sufficiently sized.
         if (.not.allocated(radiusHistory)) then
            allocate(radiusHistory(2,countComponentsActive+activeComponentMaximumIncrement))
            radiusHistory=-1.0d0
         else if (size(radiusHistory,dim=2) < countComponentsActive) then
            activeComponentMaximumCurrent=size(radiusHistory,dim=2)
            call Move_Alloc(radiusHistory,radiusHistoryTemporary)
            allocate(radiusHistory(2,countComponentsActive+activeComponentMaximumIncrement))
            radiusHistory(:,                              1:                      activeComponentMaximumCurrent  )=radiusHistoryTemporary
            radiusHistory(:,activeComponentMaximumCurrent+1:countComponentsActive+activeComponentMaximumIncrement)=-1.0d0
            deallocate(radiusHistoryTemporary)
         end if
         ! Detect oscillations in the radius solver. Only do this after a few bisection iterations have passed as we don't want to
         ! declare a true oscillation until the solver has had time to "burn in".
         if (countIterations == 1) radiusHistory=-1.0d0
         if     (                                                                                                                                                                            &
              &             countIterations                                                                                                                >  iterationsForBisectionMinimum  &
              &  .and. all( radiusHistory(:,countComponentsActive)                                                                                         >= 0.0d0                        ) &
              &  .and.     (radiusHistory(2,countComponentsActive)-radiusHistory(1,countComponentsActive))*(radiusHistory(1,countComponentsActive)-radius) <  0.0d0                          &
              & ) then
            ! An oscillation has been detected - attempt to break out of it. The following heuristic has been found to work quite
            ! well - we bisect previous solutions in the oscillating sequence in a variety of different ways
            ! (arithmetic/geometric and using the current+previous or two previous solutions), alternating the bisection method
            ! sequentially. There's no guarantee that this will work in every situation however.
            select case (mod(countIterations,4))
            case (0)
               radius=sqrt  (radius                                *radiusHistory(1,countComponentsActive))
            case (1)
               radius=0.5d0*(radius                                +radiusHistory(1,countComponentsActive))
            case (2)
               radius=sqrt  (radiusHistory(1,countComponentsActive)*radiusHistory(2,countComponentsActive))
            case (3)
               radius=0.5d0*(radiusHistory(1,countComponentsActive)+radiusHistory(2,countComponentsActive))
            end select
            radiusHistory(:,countComponentsActive)=-1.0d0
         end if
         radiusHistory(2,countComponentsActive)=radiusHistory(1,countComponentsActive)
         radiusHistory(1,countComponentsActive)=radius
         ! Compute a fit measure.
         if (radius > 0.0d0 .and. radiusNew > 0.0d0) fitMeasure=fitMeasure+abs(log(radiusNew/radius))
         ! Set radius to new radius.
         radius=radiusNew
         ! Catch unphysical states.
         if (radius <= 0.0d0) then
            if (displayVerbosity() < verbosityLevelStandard) call displayVerbositySet(verbosityLevelStandard)
            call node%serializeASCII()
            message='radius has reached zero for node '
            message=message//node%index()//' - report follows:'       //char(10)
            write (label,'(e12.6)') specificAngularMomentum
            message=message//'  specific angular momentum:    '//label//char(10)
            write (label,'(e12.6)') velocity
            message=message//'  rotation velocity:            '//label//char(10)
            write (label,'(e12.6)') sqrt(darkMatterVelocitySquared)
            message=message//'   -> dark matter contribution: '//label//char(10)
            write (label,'(e12.6)') sqrt(baryonicVelocitySquared  )
            message=message//'   -> baryonic contribution:    '//label
            call Error_Report(message//{introspection:location})
         end if
      end if
      ! Set the component size to new radius and velocity.
      call          radiusSet(node,radius  )
      call        velocitySet(node,velocity)
      call Calculations_Reset(node         )
      return
    end subroutine radiusSolve

  end subroutine equilibriumSolve

  subroutine equilibriumRevert(self,node)
    !!{
    Revert radii for the equilibrium galactic structure solve.
    !!}
    implicit none
    class(galacticStructureSolverEquilibrium), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Simply record that reversion should be performed on the next call to the solver.
    revertStructure=.true.
    return
  end subroutine equilibriumRevert
