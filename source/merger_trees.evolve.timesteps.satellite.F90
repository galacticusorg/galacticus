!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a time-stepping criterion for merger tree evolution which stops evolution when a merger is
!% about to happen.

module Merger_Tree_Timesteps_Satellite
  implicit none
  private
  public :: Merger_Tree_Timestep_Satellite

  ! Flag indicating whether this module is initialized.
  logical          :: mergerTimestepsInitialized          =.false.                                       
  
  ! Flag indicating if this module is limiting timesteps.
  logical          :: limitTimesteps                                                                     
  
  ! The largest time difference allowed between satellite and merge target at the time or merging.
  double precision :: mergeTargetTimeOffsetMaximumAbsolute        , mergeTargetTimeOffsetMaximumRelative 
  
contains

  !# <timeStepsTask>
  !#  <unitName>Merger_Tree_Timestep_Satellite</unitName>
  !# </timeStepsTask>
  subroutine Merger_Tree_Timestep_Satellite(thisNode,timeStep,End_Of_Timestep_Task,report,lockNode,lockType)
    !% Determines the timestep to go to the time at which the node merges.
    use Galacticus_Nodes
    use Evolve_To_Time_Reports
    use Merger_Trees_Evolve_Timesteps_Template
    use Input_Parameters
    use ISO_Varying_String
    use String_Handling
    use ISO_Varying_String
    implicit none
    type            (treeNode                     ), intent(inout)          , pointer :: thisNode                                                
    procedure       (End_Of_Timestep_Task_Template), intent(inout)          , pointer :: End_Of_Timestep_Task                                    
    double precision                               , intent(inout)                    :: timeStep                                                
    logical                                        , intent(in   )                    :: report                                                  
    type            (treeNode                     ), intent(inout), optional, pointer :: lockNode                                                
    type            (varying_string               ), intent(inout), optional          :: lockType                                                
    type            (treeNode                     )                         , pointer :: hostNode                                                
    class           (nodeComponentBasic           )                         , pointer :: hostBasicComponent    , thisBasicComponent              
    class           (nodeComponentSatellite       )                         , pointer :: thisSatelliteComponent                                  
    double precision                                                                  :: mergeTargetTimeMinimum, mergeTargetTimeOffsetMaximum, & 
         &                                                                               timeStepAllowed       , timeUntilMerging                
    
    ! Initialize the module.
    if (.not.mergerTimestepsInitialized) then
       !$omp critical (Merger_Tree_Timestep_Satellite_Initialize)
       if (.not.mergerTimestepsInitialized) then
          ! Check that the merge time property exists.
          limitTimesteps=defaultSatelliteComponent%mergeTimeIsGettable()
          
          ! Get parameters controlling time maximum allowed time difference between galaxies at merging.
          !@ <inputParameter>
          !@   <name>mergeTargetTimeOffsetMaximumAbsolute</name>
          !@   <defaultValue>0.01</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The maximum absolute time difference (in Gyr) allowed between merging pairs of galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>timeStepping</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergeTargetTimeOffsetMaximumAbsolute',mergeTargetTimeOffsetMaximumAbsolute,defaultValue=0.010d0)
          !@ <inputParameter>
          !@   <name>mergeTargetTimeOffsetMaximumRelative</name>
          !@   <defaultValue>0.001</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@      The maximum time difference (relative to the cosmic time at the merger epoch) allowed between merging pairs of galaxies.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>timeStepping</group>
          !@ </inputParameter>
          call Get_Input_Parameter('mergeTargetTimeOffsetMaximumRelative',mergeTargetTimeOffsetMaximumRelative,defaultValue=0.001d0)
          ! Flag that the module is initialized.
          mergerTimestepsInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Timestep_Satellite_Initialize)
    end if

    ! Exit if we are not limiting timesteps.
    if (.not.limitTimesteps) return

    ! Get the satellite component.
    thisSatelliteComponent => thisNode%satellite()
    
    ! Get the time until this node merges.
    timeUntilMerging=thisSatelliteComponent%mergeTime()

    ! If time is negative, implies this is not a satellite, so return.
    if (timeUntilMerging < 0.0d0) return

    ! Get the basic component.
    thisBasicComponent => thisNode%basic()

    ! Compute the minimum time to which the node we will merge with must have been evolved before merging is allowed.
    mergeTargetTimeOffsetMaximum=min(mergeTargetTimeOffsetMaximumAbsolute,(thisBasicComponent%time()+timeUntilMerging)&
         &*mergeTargetTimeOffsetMaximumRelative)
    mergeTargetTimeMinimum=thisBasicComponent%time()+timeUntilMerging-mergeTargetTimeOffsetMaximum

    ! Find the node to merge with.
    hostNode           => thisNode%mergesWith()
    hostBasicComponent => hostNode%basic     ()
    if (hostBasicComponent%time() < mergeTargetTimeMinimum) then
       timeStepAllowed=max(timeUntilMerging-0.5d0*mergeTargetTimeOffsetMaximum,0.0d0)
       
       ! Set return value if our timestep is smaller than current one. Do not set an end of timestep task in this case - we want
       ! to wait for the merge target to catch up before triggering a merger.
       if (timeStepAllowed <= timeStep) then
          if (present(lockNode)) lockNode => hostNode
          if (present(lockType)) lockType =  "satellite (host)"
          timeStep=timeStepAllowed
          End_Of_Timestep_Task => null()
       end if
       if (report) call Evolve_To_Time_Report("satellite (host): ",timeStep)
    else
       ! Set return value if our timestep is smaller than current one.
       if (timeUntilMerging <= timeStep) then
          if (present(lockNode)) lockNode => hostNode
          if (present(lockType)) lockType =  "satellite (self)"
          timeStep=timeUntilMerging
          ! Trigger a merger event only if the target node has no children. If it has children, we need to wait for them to be
          ! evolved before merging.
          if (.not.associated(hostNode%firstChild)) then
             End_Of_Timestep_Task => Satellite_Merger_Process
          else
             End_Of_Timestep_Task => null()
          end if
       end if
       if (report) call Evolve_To_Time_Report("satellite (self): ",timeStep,hostNode%index())
    end if
    return
  end subroutine Merger_Tree_Timestep_Satellite

  subroutine Satellite_Merger_Process(thisTree,thisNode,deadlockStatus)
    !% Process a satellite node which has undergone a merger with its host node.
    use Merger_Trees
    use Galacticus_Nodes
    use Merger_Trees_Evolve_Deadlock_Status
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    !# <include directive="satelliteMergerTask" type="moduleUse">
    include 'merger_trees.evolve.timesteps.satellite.moduleUse.inc'
    !# </include>
    implicit none
    type   (mergerTree    ), intent(in   )          :: thisTree       
    type   (treeNode      ), intent(inout), pointer :: thisNode       
    integer                , intent(inout)          :: deadlockStatus 
    type   (varying_string)                         :: message        
    
    ! Report if necessary.
    ! Report if necessary.
    if (Galacticus_Verbosity_Level() >= verbosityInfo) then
       message='Satellite node ['
       message=message//thisNode%index()//'] is being merged'
       call Galacticus_Display_Message(message)
    end if

    ! Allow arbitrary routines to process the merger.
    !# <include directive="satelliteMergerTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode</functionArgs>
    include 'merger_trees.evolve.timesteps.satellite.inc'
    !# </include>

    ! Finally remove the satellite node from the host and merge targets and destroy it.
    call thisNode%removeFromHost  ()
    call thisNode%removeFromMergee()
    call thisNode%destroy         ()
    deallocate(thisNode)
    thisNode => null()

    ! The tree was changed, so mark that it is not deadlocked.
    deadlockStatus=isNotDeadlocked
    return
  end subroutine Satellite_Merger_Process

end module Merger_Tree_Timesteps_Satellite
