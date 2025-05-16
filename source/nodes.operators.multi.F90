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
  Implements a multi node operator class.
  !!}

  type, public :: multiProcessList
     class(nodeOperatorClass), pointer :: process_ => null()
     type (multiProcessList ), pointer :: next     => null()
  end type multiProcessList

  !![
  <nodeOperator name="nodeOperatorMulti">
   <description>A multi node operator property process class.</description>
   <linkedList type="multiProcessList" variable="processes" next="next" object="process_" objectType="nodeOperatorClass"/>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorMulti
     !!{
     A multi node operator output process class, which applies multiple node operators.
     !!}
     private
     type(multiProcessList), pointer :: processes => null()
   contains
     !![
     <methods>
	<method method="isActive" description="Return true if the operators are active for the given {\normalfont \ttfamily node}."/>
     </methods>
     !!]
     final     ::                                        multiDestructor
     procedure :: nodeTreeInitialize                  => multiNodeTreeInitialize
     procedure :: nodeInitialize                      => multiNodeInitialize
     procedure :: nodesMerge                          => multiNodesMerge
     procedure :: nodePromote                         => multiNodePromote
     procedure :: galaxiesMerge                       => multiGalaxiesMerge
     procedure :: differentialEvolutionPre            => multiDifferentialEvolutionPre
     procedure :: differentialEvolution               => multiDifferentialEvolution
     procedure :: differentialEvolutionScales         => multiDifferentialEvolutionScales
     procedure :: differentialEvolutionAnalytics      => multiDifferentialEvolutionAnalytics
     procedure :: predeterminedSolveAnalytics         => multiPredeterminedSolveAnalytics
     procedure :: differentialEvolutionSolveAnalytics => multiDifferentialEvolutionSolveAnalytics
     procedure :: differentialEvolutionInactives      => multiDifferentialEvolutionInactives
     procedure :: differentialEvolutionStepFinalState => multiDifferentialEvolutionStepFinalState
     procedure :: differentialEvolutionPost           => multiDifferentialEvolutionPost
     procedure :: differentialEvolutionPostStep       => multiDifferentialEvolutionPostStep
     procedure :: isActive                            => multiIsActive
  end type nodeOperatorMulti

  interface nodeOperatorMulti
     !!{
     Constructors for the \refClass{nodeOperatorMulti} node operator class.
     !!}
     module procedure multiConstructorParameters
     module procedure multiConstructorInternal
  end interface nodeOperatorMulti

contains

  function multiConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorMulti} node operator property process class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodeOperatorMulti)                :: self
    type   (inputParameters  ), intent(inout) :: parameters
    type   (multiProcessList ), pointer       :: process_
    integer                                   :: i

    self    %processes => null()
    process_           => null()
    do i=1,parameters%copiesCount('nodeOperator',zeroIfNotPresent=.true.)
       if (associated(process_)) then
          allocate(process_%next)
          process_ => process_%next
       else
          allocate(self%processes)
          process_ => self%processes
       end if
       !![
       <objectBuilder class="nodeOperator" name="process_%process_" source="parameters" copy="i" />
       !!]
    end do
    !![
    <inputParametersValidate source="parameters" multiParameters="nodeOperator"/>
    !!]
    return
  end function multiConstructorParameters

  function multiConstructorInternal(processes) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorMulti} output process property process class.
    !!}
    implicit none
    type(nodeOperatorMulti)                         :: self
    type(multiProcessList ), target , intent(in   ) :: processes
    type(multiProcessList ), pointer                :: process_

    self    %processes => processes
    process_           => processes
    do while (associated(process_))
       !![
       <referenceCountIncrement owner="process_" object="process_"/>
       !!]
       process_ => process_%next
    end do
    return
  end function multiConstructorInternal

  subroutine multiDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorMulti} output process property process class.
    !!}
    implicit none
    type(nodeOperatorMulti), intent(inout) :: self
    type(multiProcessList ), pointer       :: process_, processNext

    if (associated(self%processes)) then
       process_ => self%processes
       do while (associated(process_))
          processNext => process_%next
          !![
          <objectDestructor name="process_%process_"/>
          !!]
          deallocate(process_)
          process_ => processNext
       end do
    end if
    return
  end subroutine multiDestructor

  subroutine multiNodeTreeInitialize(self,node)
    !!{
    Perform node tree initialization.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout), target  :: self
    type (treeNode         ), intent(inout), target  :: node
    type (multiProcessList )               , pointer :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%nodeTreeInitialize(node)
       process_ => process_%next
    end do
    return
  end subroutine multiNodeTreeInitialize

  subroutine multiNodeInitialize(self,node)
    !!{
    Perform node initialization.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout), target  :: self
    type (treeNode         ), intent(inout), target  :: node
    type (multiProcessList )               , pointer :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%nodeInitialize(node)
       process_ => process_%next
    end do
    return
  end subroutine multiNodeInitialize

  subroutine multiNodesMerge(self,node)
    !!{
    Act on a merger between galaxies.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%nodesMerge(node)
       process_ => process_%next
    end do
    return
  end subroutine multiNodesMerge

  subroutine multiNodePromote(self,node)
    !!{
    Act on a node promotion event.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%nodePromote(node)
       process_ => process_%next
    end do
    return
  end subroutine multiNodePromote

  subroutine multiGalaxiesMerge(self,node)
    !!{
    Act on a merger between galaxies.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%galaxiesMerge(node)
       process_ => process_%next
    end do
    return
  end subroutine multiGalaxiesMerge

  subroutine multiDifferentialEvolutionPre(self,node)
    !!{
    Act on a node before differential evolution.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolutionPre(node)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolutionPre

  subroutine multiDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scales prior to differential evolution.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolutionScales(node)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolutionScales

  subroutine multiDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark (meta-)properties as analytically solved for the ODE solver prior to differential evolution.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolutionAnalytics(node)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolutionAnalytics

  subroutine multiDifferentialEvolutionInactives(self,node)
    !!{
    Mark (meta-)properties as inactive for the ODE solver prior to differential evolution.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolutionInactives(node)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolutionInactives

  subroutine multiDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Act on a node during differential evolution.
    !!}
    implicit none
    class    (nodeOperatorMulti), intent(inout), target  :: self
    type     (treeNode         ), intent(inout), target  :: node
    logical                     , intent(inout)          :: interrupt
    procedure(interruptTask    ), intent(inout), pointer :: functionInterrupt
    integer                     , intent(in   )          :: propertyType
    type     (multiProcessList )               , pointer :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolution(node,interrupt,functionInterrupt,propertyType)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolution

  subroutine multiDifferentialEvolutionSolveAnalytics(self,node,time)
    !!{
    Set the values of analytically-solvable properties of a node during differential evolution.
    !!}
    implicit none
    class           (nodeOperatorMulti), intent(inout) :: self
    type            (treeNode         ), intent(inout) :: node
    double precision                   , intent(in   ) :: time
    type            (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolutionSolveAnalytics(node,time)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolutionSolveAnalytics

  subroutine multiPredeterminedSolveAnalytics(self,node,time)
    !!{
    Set the pre-determined values of analytically-solvable properties of a node.
    !!}
    implicit none
    class           (nodeOperatorMulti), intent(inout) :: self
    type            (treeNode         ), intent(inout) :: node
    double precision                   , intent(in   ) :: time
    type            (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%predeterminedSolveAnalytics(node,time)
       process_ => process_%next
    end do
    return
  end subroutine multiPredeterminedSolveAnalytics

  subroutine multiDifferentialEvolutionStepFinalState(self,node)
    !!{
    Act on a node after a differential evolution ODE step.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolutionStepFinalState(node)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolutionStepFinalState

  subroutine multiDifferentialEvolutionPost(self,node)
    !!{
    Act on a node after differential evolution.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    type (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolutionPost(node)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolutionPost

  subroutine multiDifferentialEvolutionPostStep(self,node,status)
    !!{
    Act on a node after a differential evolution step.
    !!}
    implicit none
    class  (nodeOperatorMulti), intent(inout) :: self
    type   (treeNode         ), intent(inout) :: node
    integer                   , intent(inout) :: status
    type   (multiProcessList ), pointer       :: process_

    if (.not.self%isActive(node)) return
    process_ => self%processes
    do while (associated(process_))
       call process_%process_%differentialEvolutionPostStep(node,status)
       process_ => process_%next
    end do
    return
  end subroutine multiDifferentialEvolutionPostStep

  logical function multiIsActive(self,node) result(isActive)
    !!{
    Return true if the operators are active for the given {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(nodeOperatorMulti), intent(inout) :: self
    type (treeNode         ), intent(inout) :: node
    !$GLC attributes unused :: self, node
    
    isActive=.true.
    return
  end function multiIsActive
