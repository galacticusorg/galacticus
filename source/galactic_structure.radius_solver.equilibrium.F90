!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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

  use :: Dark_Matter_Profiles    , only : darkMatterProfileClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Galactic_Structure      , only : galacticStructureClass

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
     class           (darkMatterProfileClass   ), pointer :: darkMatterProfile_         => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_      => null()
     class           (galacticStructureClass   ), pointer :: galacticStructure_         => null()
   contains
     final     ::             equilibriumDestructor
     procedure :: solve    => equilibriumSolve
     procedure :: revert   => equilibriumRevert
     procedure :: autoHook => equilibriumAutoHook
  end type galacticStructureSolverEquilibrium

  interface galacticStructureSolverEquilibrium
     !!{
     Constructors for the {\normalfont \ttfamily equilibrium} galactic structure solver class.
     !!}
     module procedure equilibriumConstructorParameters
     module procedure equilibriumConstructorInternal
  end interface galacticStructureSolverEquilibrium

  ! Module-scope variables used to communicate current state of radius solver.
  integer                                               :: equilibriumActiveComponentCount        , equilibriumIterationCount
  double precision                                      :: equilibriumFitMeasure
  type            (treeNode), pointer                   :: equilibriumHaloNode
  double precision          , allocatable, dimension(:) :: equilibriumRadiusStored                , equilibriumVelocityStored
  logical                                               :: equilibriumRevertStructure     =.false.
  !$omp threadprivate(equilibriumIterationCount,equilibriumActiveComponentCount,equilibriumFitMeasure,equilibriumHaloNode,equilibriumRadiusStored,equilibriumVelocityStored,equilibriumRevertStructure)

contains

  function equilibriumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily equilibrium} galactic structure solver class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticStructureSolverEquilibrium)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (darkMatterProfileClass            ), pointer       :: darkMatterProfile_
    class           (darkMatterProfileDMOClass         ), pointer       :: darkMatterProfileDMO_
    class           (galacticStructureClass            ), pointer       :: galacticStructure_
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
      <description>Specifies whether or not gravity from baryons is included when solving for sizes of galactic components in equilibriumally contracted dark matter halos.</description>
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
    <objectBuilder class="darkMatterProfile"    name="darkMatterProfile_"    source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="galacticStructure"    name="galacticStructure_"    source="parameters"/>
    !!]
    self=galacticStructureSolverEquilibrium(convergenceFailureIsFatal,useFormationHalo,includeBaryonGravity,solutionTolerance,solveForInactiveProperties,darkMatterProfile_,darkMatterProfileDMO_,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfile_"   />
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="galacticStructure_"   />
    !!]
    return
  end function equilibriumConstructorParameters

  function equilibriumConstructorInternal(convergenceFailureIsFatal,useFormationHalo,includeBaryonGravity,solutionTolerance,solveForInactiveProperties,darkMatterProfile_,darkMatterProfileDMO_,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily equilibrium} galactic structure solver class.
    !!}
    implicit none
    type            (galacticStructureSolverEquilibrium)                        :: self
    logical                                             , intent(in   )         :: useFormationHalo          , includeBaryonGravity     , &
         &                                                                         solveForInactiveProperties, convergenceFailureIsFatal
    double precision                                    , intent(in   )         :: solutionTolerance
    class           (darkMatterProfileClass            ), intent(in   ), target :: darkMatterProfile_
    class           (darkMatterProfileDMOClass         ), intent(in   ), target :: darkMatterProfileDMO_
    class           (galacticStructureClass            ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="convergenceFailureIsFatal, useFormationHalo, includeBaryonGravity, solutionTolerance, solveForInactiveProperties, *darkMatterProfile_, *darkMatterProfileDMO_, *galacticStructure_"/>
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
    call   preDerivativeEvent%attach(self,equilibriumSolvePreDeriativeHook,openMPThreadBindingAtLevel                                                             )
    call      postEvolveEvent%attach(self,equilibriumSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverEquilibrium',dependencies=dependencies)
    call satelliteMergerEvent%attach(self,equilibriumSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverEquilibrium',dependencies=dependencies)
    call   nodePromotionEvent%attach(self,equilibriumSolveHook            ,openMPThreadBindingAtLevel,label='structureSolverEquilibrium',dependencies=dependencies)
    return
  end subroutine equilibriumAutoHook

  subroutine equilibriumDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily equilibrium} galactic structure solver class.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent, postEvolveEvent, preDerivativeEvent, satelliteMergerEvent
    implicit none
    type(galacticStructureSolverEquilibrium), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfile_"   />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%galacticStructure_"   />
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
       if (propertyType /= propertyTypeInactive .or. self%solveForInactiveProperties) call self%solve(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine equilibriumSolvePreDeriativeHook

  subroutine equilibriumSolve(self,node)
    !!{
    Solve for the structure of galactic components.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Display            , only : displayMessage
    use :: Error              , only : Error_Report      , Warn
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    implicit none
    class           (galacticStructureSolverEquilibrium), intent(inout)         :: self
    type            (treeNode                          ), intent(inout), target :: node
    logical                                             , parameter             :: specificAngularMomentumRequired=.true.
    integer                                             , parameter             :: iterationMaximum               =  100
    procedure       (solverGet                         ), pointer               :: radiusGet                             , velocityGet
    procedure       (solverSet                         ), pointer               :: radiusSet                             , velocitySet
    logical                                                                     :: componentActive
    double precision                                                            :: specificAngularMomentum

    ! Check that the galaxy is physical plausible. In this equilibrium solver, we don't act on this.
    node%isPhysicallyPlausible=.true.
    node%isSolvable           =.true.
    include 'galactic_structure.radius_solver.plausible.inc'
    if (node%isPhysicallyPlausible) then
       ! Initialize the solver state.
       equilibriumIterationCount=0
       equilibriumFitMeasure    =2.0d0*self%solutionTolerance
       ! Determine which node to use for halo properties.
       if (self%useFormationHalo) then
          if (.not.associated(node%formationNode)) call Error_Report('no formation node exists'//{introspection:location})
          equilibriumHaloNode => node%formationNode
       else
          equilibriumHaloNode => node
       end if
       ! Begin iteration to find a converged solution.
       do while (equilibriumIterationCount <= 2 .or. ( equilibriumFitMeasure > self%solutionTolerance .and. equilibriumIterationCount < iterationMaximum ) )
          call Calculations_Reset(node)
          equilibriumIterationCount      =equilibriumIterationCount+1
          equilibriumActiveComponentCount=0
          if (equilibriumIterationCount > 1) equilibriumFitMeasure=0.0d0
          include 'galactic_structure.radius_solver.tasks.inc'
          ! Check that we have some active components.
          if (equilibriumActiveComponentCount == 0) then
             equilibriumFitMeasure=0.0d0
             exit
          else
             ! Normalize the fit measure by the number of active components.
             equilibriumFitMeasure=equilibriumFitMeasure/dble(equilibriumActiveComponentCount)
          end if
       end do
       ! Check that we found a converged solution.
       if (equilibriumFitMeasure > self%solutionTolerance) then
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
    equilibriumRevertStructure=.false.
    return

  contains

    subroutine radiusSolve(node,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet)
      !!{
      Solve for the equilibrium radius of the given component.
      !!}
      use :: Display                         , only : displayVerbosity               , displayVerbositySet, verbosityLevelStandard
      use :: Galactic_Structure_Options      , only : massTypeBaryonic
      use :: Error                           , only : Error_Report
      use :: ISO_Varying_String              , only : varying_string
      use :: Memory_Management               , only : allocateArray                  , deallocateArray
      use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
      use :: String_Handling                 , only : operator(//)
      implicit none
      type            (treeNode          ), intent(inout)                     :: node
      double precision                    , intent(in   )                     :: specificAngularMomentum
      procedure       (solverGet         ), intent(in   ) , pointer           :: radiusGet                           , velocityGet
      procedure       (solverSet         ), intent(in   ) , pointer           :: radiusSet                           , velocitySet
      double precision                    , dimension(:,:), allocatable, save :: radiusHistory
      !$omp threadprivate(radiusHistory)
      double precision                    , dimension(:,:), allocatable       :: radiusHistoryTemporary
      double precision                    , dimension(:  ), allocatable       :: equilibriumRadiusStoredTmp          , equilibriumVelocityStoredTmp
      integer                             , dimension(1  ), parameter         :: storeIncrement                 =[10]
      integer                             , parameter                         :: iterationsForBisectionMinimum  = 10
      integer                             , parameter                         :: activeComponentMaximumIncrement=  2
      integer                                                                 :: activeComponentMaximumCurrent
      character       (len=14            )                                    :: label
      type            (varying_string    )                                    :: message
      double precision                                                        :: baryonicVelocitySquared             , darkMatterMassFinal         , &
           &                                                                     darkMatterVelocitySquared           , velocity                    , &
           &                                                                     radius                              , radiusNew

      ! Count the number of active comonents.
      equilibriumActiveComponentCount=equilibriumActiveComponentCount+1
      if (equilibriumIterationCount == 1) then
         ! If structure is to be reverted, do so now.
         if (equilibriumRevertStructure.and.allocated(equilibriumRadiusStored)) then
            call   radiusSet(node,  equilibriumRadiusStored(equilibriumActiveComponentCount))
            call velocitySet(node,equilibriumVelocityStored(equilibriumActiveComponentCount))
         end if
         ! On first iteration, see if we have a previous radius set for this component.
         radius=radiusGet(node)
         if (radius <= 0.0d0) then
            ! No previous radius was set, so make a simple estimate of sizes of all components ignoring equilibrium contraction and self-gravity.
            ! Find the radius in the dark matter profile with the required specific angular momentum
            radius  =self%darkMatterProfileDMO_%radiusFromSpecificAngularMomentum(equilibriumHaloNode,specificAngularMomentum)
            ! Find the velocity at this radius.
            velocity=self%darkMatterProfileDMO_%circularVelocity                 (equilibriumHaloNode,radius                 )
         else
            ! A previous radius was set, so use it, and the previous circular velocity, as the initial guess.
            velocity=velocityGet(node)
         end if
         ! If structure is not being reverted, store the new values of radius and velocity.
         if (.not.equilibriumRevertStructure .and. equilibriumIterationCount == 1) then
            ! Store these quantities.
            if (.not.allocated(equilibriumRadiusStored)) then
               call allocateArray(  equilibriumRadiusStored,storeIncrement)
               call allocateArray(equilibriumVelocityStored,storeIncrement)
               equilibriumRadiusStored  =0.0d0
               equilibriumVelocityStored=0.0d0
            else if (equilibriumActiveComponentCount > size(equilibriumRadiusStored)) then
               call move_alloc     (  equilibriumRadiusStored,        equilibriumRadiusStoredTmp                )
               call move_alloc     (equilibriumVelocityStored,      equilibriumVelocityStoredTmp                )
               call allocateArray  (  equilibriumRadiusStored,shape(  equilibriumRadiusStoredTmp)+storeIncrement)
               call allocateArray  (equilibriumVelocityStored,shape(equilibriumVelocityStoredTmp)+storeIncrement)
               equilibriumRadiusStored  (                        1:size(  equilibriumRadiusStoredTmp))=  equilibriumRadiusStoredTmp
               equilibriumVelocityStored(                        1:size(equilibriumVelocityStoredTmp))=equilibriumVelocityStoredTmp
               equilibriumRadiusStored  (size(  equilibriumRadiusStoredTmp)+1:size(  equilibriumRadiusStored   ))=0.0d0
               equilibriumVelocityStored(size(equilibriumVelocityStoredTmp)+1:size(equilibriumVelocityStored   ))=0.0d0
               call deallocateArray(  equilibriumRadiusStoredTmp)
               call deallocateArray(equilibriumVelocityStoredTmp)
            end if
            equilibriumRadiusStored  (equilibriumActiveComponentCount)=radius
            equilibriumVelocityStored(equilibriumActiveComponentCount)=velocity
         end if
       else
         ! On subsequent iterations do the full calculation providing component has non-zero specific angular momentum.
         if (specificAngularMomentum <= 0.0d0) return
         ! Get current radius of the component.
         radius                   =radiusGet(node)
         ! Find the enclosed mass in the dark matter halo.
         darkMatterMassFinal      =self%darkMatterProfile_%enclosedMass(equilibriumHaloNode,radius)
         ! Compute dark matter contribution to rotation curve.
         darkMatterVelocitySquared=gravitationalConstantGalacticus*darkMatterMassFinal/radius
         ! Compute baryonic contribution to rotation curve.
         if (self%includeBaryonGravity) then
            baryonicVelocitySquared=self%galacticStructure_%velocityRotation(node,radius,massType=massTypeBaryonic)**2
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
            call allocateArray(radiusHistory,[2,equilibriumActiveComponentCount+activeComponentMaximumIncrement])
            radiusHistory=-1.0d0
         else if (size(radiusHistory,dim=2) < equilibriumActiveComponentCount) then
            activeComponentMaximumCurrent=size(radiusHistory,dim=2)
            call Move_Alloc(radiusHistory,radiusHistoryTemporary)
            call allocateArray(radiusHistory,[2,equilibriumActiveComponentCount+activeComponentMaximumIncrement])
            radiusHistory(:,                              1:                                activeComponentMaximumCurrent  )=radiusHistoryTemporary
            radiusHistory(:,activeComponentMaximumCurrent+1:equilibriumActiveComponentCount+activeComponentMaximumIncrement)=-1.0d0
            call deallocateArray(radiusHistoryTemporary)
         end if
         ! Detect oscillations in the radius solver. Only do this after a few bisection iterations have passed as we don't want to
         ! declare a true oscillation until the solver has had time to "burn in".
         if (equilibriumIterationCount == 1) radiusHistory=-1.0d0
         if     (                                                                                                                                                                                                          &
              &             equilibriumIterationCount                                                                                                                                    >  iterationsForBisectionMinimum  &
              &  .and. all( radiusHistory(:,equilibriumActiveComponentCount)                                                                                                             >= 0.0d0                        ) &
              &  .and.     (radiusHistory(2,equilibriumActiveComponentCount)-radiusHistory(1,equilibriumActiveComponentCount))*(radiusHistory(1,equilibriumActiveComponentCount)-radius) <  0.0d0                          &
              & ) then
            ! An oscillation has been detected - attempt to break out of it. The following heuristic has been found to work quite
            ! well - we bisect previous solutions in the oscillating sequence in a variety of different ways
            ! (arithmetic/geometric and using the current+previous or two previous solutions), alternating the bisection method
            ! sequentially. There's no guarantee that this will work in every situation however.
            select case (mod(equilibriumIterationCount,4))
            case (0)
               radius=sqrt  (radius                                          *radiusHistory(1,equilibriumActiveComponentCount))
            case (1)
               radius=0.5d0*(radius                                          +radiusHistory(1,equilibriumActiveComponentCount))
            case (2)
               radius=sqrt  (radiusHistory(1,equilibriumActiveComponentCount)*radiusHistory(2,equilibriumActiveComponentCount))
            case (3)
               radius=0.5d0*(radiusHistory(1,equilibriumActiveComponentCount)+radiusHistory(2,equilibriumActiveComponentCount))
            end select
            radiusHistory(:,equilibriumActiveComponentCount)=-1.0d0
         end if
         radiusHistory(2,equilibriumActiveComponentCount)=radiusHistory(1,equilibriumActiveComponentCount)
         radiusHistory(1,equilibriumActiveComponentCount)=radius
         ! Compute a fit measure.
         if (radius > 0.0d0 .and. radiusNew > 0.0d0) equilibriumFitMeasure=equilibriumFitMeasure+abs(log(radiusNew/radius))
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
      call   radiusSet(node,radius  )
      call velocitySet(node,velocity)
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
    equilibriumRevertStructure=.true.
    return
  end subroutine equilibriumRevert
