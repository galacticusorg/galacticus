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
  Contains a merger tree evolution timestep class which limits the step to the next satellite merger.
  !!}

  use :: Nodes_Operators, only : nodeOperatorClass

  !![
  <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepSatellite">
   <description>  
    A merger tree evolution timestepping class which enforces the following for satellite \glspl{node}. If the satellite's merge
    target has been advanced to at least a time of $t_\mathrm{required} = t_\mathrm{satellite} + \Delta t_\mathrm{merge} - \delta
    t_\mathrm{merge,maximum}$ then
    \begin{equation}
     \Delta t \le \Delta t_\mathrm{merge},
    \end{equation}
    where $t_\mathrm{satellite}$ is the current time for the satellite \gls{node}, $\Delta t_\mathrm{merge}$ is the time until the
    satellite is due to merge and $\delta t_\mathrm{merge,maximum}$ is the maximum allowed time difference between merging
    galaxies. This ensures that the satellite is not evolved past the time at which it is due to merge. If this criterion is the
    limiting criteria for $\Delta t$ then the merging of the satellite will be triggered at the end of the timestep.
  
    If the merge target has not been advanced to at least $t_\mathrm{required}$ then instead
    \begin{equation}
     \Delta t \le \hbox{max}(\Delta t_\mathrm{merge}-\delta t_\mathrm{merge,maximum}/2,0),
    \end{equation}
    is asserted to ensure that the satellite does not reach the time of merging until its merge target is sufficiently close (within
    $\delta t_\mathrm{merge,maximum}$) of the time of merging.
   </description>
  </mergerTreeEvolveTimestep>
  !!]
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepSatellite
     !!{
     Implementation of a merger tree evolution timestep class which limits the step to the next satellite merger.
     !!}
     private
     class           (nodeOperatorClass), pointer :: nodeOperator_             => null()
     double precision                             :: timeOffsetMaximumAbsolute          , timeOffsetMaximumRelative
     logical                                      :: limitTimesteps
   contains
     final     ::                 satelliteDestructor
     procedure :: timeEvolveTo => satelliteTimeEvolveTo
  end type mergerTreeEvolveTimestepSatellite

  interface mergerTreeEvolveTimestepSatellite
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepSatellite} merger tree evolution timestep class.
     !!}
     module procedure satelliteConstructorParameters
     module procedure satelliteConstructorInternal
  end interface mergerTreeEvolveTimestepSatellite

contains

  function satelliteConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepSatellite} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeEvolveTimestepSatellite)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (nodeOperatorClass                ), pointer       :: nodeOperator_
    double precision                                                   :: timeOffsetMaximumAbsolute, timeOffsetMaximumRelative

    !![
    <inputParameter>
      <name>timeOffsetMaximumAbsolute</name>
      <defaultValue>0.010d0</defaultValue>
      <description>The maximum absolute time difference (in Gyr) allowed between merging pairs of galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeOffsetMaximumRelative</name>
      <defaultValue>0.001d0</defaultValue>
      <description>The maximum time difference (relative to the cosmic time at the merger epoch) allowed between merging pairs of galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="nodeOperator" name="nodeOperator_" source="parameters"/>
    !!]
    self=mergerTreeEvolveTimestepSatellite(timeOffsetMaximumAbsolute,timeOffsetMaximumRelative,nodeOperator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodeOperator_"/>
    !!]
    return
  end function satelliteConstructorParameters

  function satelliteConstructorInternal(timeOffsetMaximumAbsolute,timeOffsetMaximumRelative,nodeOperator_) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepSatellite} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : defaultSatelliteComponent
    implicit none
    type            (mergerTreeEvolveTimestepSatellite)                        :: self
    class           (nodeOperatorClass                ), intent(in   ), target :: nodeOperator_
    double precision                                   , intent(in   )         :: timeOffsetMaximumAbsolute, timeOffsetMaximumRelative
    !![
    <constructorAssign variables="timeOffsetMaximumAbsolute, timeOffsetMaximumRelative, *nodeOperator_"/>
    !!]

    self%limitTimesteps=defaultSatelliteComponent%timeOfMergingIsGettable()
    return
  end function satelliteConstructorInternal

  subroutine satelliteDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolveTimestepSatellite} erger tree evolution timestep class.
    !!}
    implicit none
    type(mergerTreeEvolveTimestepSatellite), intent(inout) :: self

    !![
    <objectDestructor name="self%nodeOperator_"/>
    !!]
    return
  end subroutine satelliteDestructor

  double precision function satelliteTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType)
    !!{
    Determine a suitable timestep for {\normalfont \ttfamily node} such that it does not exceed the time of the next satellite merger.
    !!}
    use :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use :: Galacticus_Nodes      , only : nodeComponentBasic   , nodeComponentSatellite, treeNode
    use :: ISO_Varying_String    , only : varying_string
    implicit none
    class           (mergerTreeEvolveTimestepSatellite), intent(inout), target            :: self
    double precision                                   , intent(in   )                    :: timeEnd
    type            (treeNode                         ), intent(inout), target            :: node
    procedure       (timestepTask                     ), intent(  out), pointer           :: task
    class           (*                                ), intent(  out), pointer           :: taskSelf
    logical                                            , intent(in   )                    :: report
    type            (treeNode                         ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                   ), intent(  out)         , optional :: lockType
    type            (treeNode                         )               , pointer           :: nodeHost
    class           (nodeComponentBasic               )               , pointer           :: basicHost             , basic
    class           (nodeComponentSatellite           )               , pointer           :: satellite
    double precision                                                                      :: mergeTargetTimeMinimum, mergeTargetTimeOffsetMaximum, &
         &                                                                                   timeOfMerging
    !$GLC attributes unused :: timeEnd

    ! By default set a huge timestep so that this class has no effect.
    satelliteTimeEvolveTo           =  huge(0.0d0)
    task                            => null(     )
    taskSelf                        => null(     )
    if (present(lockNode)) lockNode => null(     )
    if (present(lockType)) lockType =  ""
    ! If not limiting timesteps, return.
    if (.not.self%limitTimesteps  ) return
    ! If node is not a satellite, return.
    if (.not.node%isSatellite   ()) return
    ! Find the time of merging.
    basic         => node     %basic        ()
    satellite     => node     %satellite    ()
    timeOfMerging =  satellite%timeOfMerging()
    ! If time is infinite, no merging will ever occur, so return.
    if (timeOfMerging == huge(0.0d0)) return
    ! Compute the minimum time to which the node we will merge with must have been evolved before merging is allowed.
    mergeTargetTimeOffsetMaximum=min(                                 &
         &                           +self%timeOffsetMaximumAbsolute, &
         &                           +self%timeOffsetMaximumRelative  &
         &                           *timeOfMerging                   &
         &                          )
    mergeTargetTimeMinimum=+timeOfMerging-mergeTargetTimeOffsetMaximum
    ! Find the node to merge with.
    nodeHost  => node    %mergesWith()
    basicHost => nodeHost%basic     ()
    if (basicHost%time() < mergeTargetTimeMinimum .and. associated(nodeHost%parent)) then
       ! Do not set an end of timestep task in this case - we want to wait for the merge target to catch up before triggering a
       ! merger.
       satelliteTimeEvolveTo           =  max(timeOfMerging-basic%time()-0.5d0*mergeTargetTimeOffsetMaximum,0.0d0)+basic%time()
       if (present(lockNode)) lockNode => nodeHost
       if (present(lockType)) lockType =  "satellite (host)"
       if (        report   ) call Evolve_To_Time_Report("satellite (host): ",satelliteTimeEvolveTo)
    else
       ! Set return value if our timestep is smaller than current one.
       if (present(lockNode)) lockNode => nodeHost
       if (present(lockType)) lockType =  "satellite (self)"
       satelliteTimeEvolveTo=timeOfMerging
       ! Trigger a merger event only if the target node has no children. If it has children, we need to wait for them to be
       ! evolved before merging.
       if (.not.associated(nodeHost%firstChild)) then
          task     => satelliteMergerProcess
          taskSelf => self
       end if
    end if
    if (report) call Evolve_To_Time_Report("satellite (self): ",satelliteTimeEvolveTo,nodeHost%index())
    return
  end function satelliteTimeEvolveTo

  subroutine satelliteMergerProcess(self,tree,node,deadlockStatus)
    !!{
    Process a satellite node which has undergone a merger with its host node.
    !!}
    use :: Display                            , only : displayMessage               , displayVerbosity, verbosityLevelInfo
    use :: Error                              , only : Error_Report
    use :: ISO_Varying_String                 , only : varying_string
    use :: Merger_Trees_Evolve_Deadlock_Status, only : deadlockStatusIsNotDeadlocked
    use :: Satellite_Promotion                , only : Satellite_Move_To_New_Host
    use :: String_Handling                    , only : operator(//)
    implicit none
    class(*                            ), intent(inout)          :: self
    type (mergerTree                   ), intent(in   )          :: tree
    type (treeNode                     ), intent(inout), pointer :: node
    type (enumerationDeadlockStatusType), intent(inout)          :: deadlockStatus
    type (treeNode                     )               , pointer :: mergee        , mergeeNext       , &
         &                                                          nodeSatellite , nodeSatelliteNext
    type (varying_string               )                         :: message
    !$GLC attributes unused :: tree

    ! Report if necessary.
    if (displayVerbosity() >= verbosityLevelInfo) then
       message='Satellite node ['
       message=message//node%index()//'] is being merged'
       call displayMessage(message)
    end if
    ! Perform any node operations on the galaxy merger, and then trigger the satellite merger event.
    select type (self)
    class is (mergerTreeEvolveTimestepSatellite)
       call self%nodeOperator_%galaxiesMerge(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    !![
    <eventHook name="satelliteMerger">
     <callWith>node</callWith>
    </eventHook>
    !!]
    ! Any mergees of the merging node must become mergees of its merge target.
    mergee => node%firstMergee
    do while (associated(mergee))
       mergeeNext => mergee%siblingMergee
       call mergee%removeFromMergee()
       mergee%siblingMergee             => node      %mergeTarget%firstMergee
       node  %mergeTarget  %firstMergee => mergee
       mergee%mergeTarget               => node      %mergeTarget
       mergee                           => mergeeNext
    end do
    ! Move the sub-sub-halo into its host's host.
    nodeSatellite => node%firstSatellite
    do while (associated(nodeSatellite))
       nodeSatelliteNext => nodeSatellite%sibling
       call Satellite_Move_To_New_Host(nodeSatellite,node%parent)
       nodeSatellite     => nodeSatelliteNext
    end do
    ! Finally remove the satellite node from the host and merge targets and destroy it.
    call node%removeFromHost  ()
    call node%removeFromMergee()
    call node%destroy         ()
    deallocate(node)
    node => null()
    ! The tree was changed, so mark that it is not deadlocked.
    deadlockStatus=deadlockStatusIsNotDeadlocked
    return
  end subroutine satelliteMergerProcess
