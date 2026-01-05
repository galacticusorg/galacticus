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
  Implements a node operator class that triggers orphanizing of satellites that reach the end of their history.
  !!}

  !![
  <nodeOperator name="nodeOperatorSatelliteOrphanize">
   <description>A node operator class that triggers orphanizing of satellites that reach the end of their history.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteOrphanize
     !!{
     A node operator class that triggers orphanizing of satellites that reach the end of their history.
     !!}
     private
   contains
     final     ::                          satelliteOrphanizeDestructor
     procedure :: autoHook              => satelliteOrphanizeAutoHook
     procedure :: differentialEvolution => satelliteOrphanizeDifferentialEvolution
  end type nodeOperatorSatelliteOrphanize
  
  interface nodeOperatorSatelliteOrphanize
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteOrphanize} node operator class.
     !!}
     module procedure satelliteOrphanizeConstructorParameters
  end interface nodeOperatorSatelliteOrphanize

contains

  function satelliteOrphanizeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteOrphanize} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorSatelliteOrphanize)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=nodeOperatorSatelliteOrphanize()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteOrphanizeConstructorParameters
  
  subroutine satelliteOrphanizeAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent, branchJumpPostProcessEvent, interTreePostProcessEvent, openMPThreadBindingAtLevel
    implicit none
    class(nodeOperatorSatelliteOrphanize), intent(inout) :: self
    
    call   satelliteHostChangeEvent%attach(self,satelliteHostChange  ,openMPThreadBindingAtLevel,label='satelliteOrphanize')
    call branchJumpPostProcessEvent%attach(self,branchJumpPostProcess,openMPThreadBindingAtLevel,label='satelliteOrphanize')
    call  interTreePostProcessEvent%attach(self,branchJumpPostProcess,openMPThreadBindingAtLevel,label='satelliteOrphanize')
   return
  end subroutine satelliteOrphanizeAutoHook

  subroutine satelliteOrphanizeDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteOrphanize} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent, branchJumpPostProcessEvent, interTreePostProcessEvent
    implicit none
    type(nodeOperatorSatelliteOrphanize), intent(inout) :: self

    if (  satelliteHostChangeEvent%isAttached(self,satelliteHostChange)) call   satelliteHostChangeEvent%detach(self,satelliteHostChange)
    if (branchJumpPostProcessEvent%isAttached(self,satelliteHostChange)) call branchJumpPostProcessEvent%detach(self,satelliteHostChange)
    if ( interTreePostProcessEvent%isAttached(self,satelliteHostChange)) call  interTreePostProcessEvent%detach(self,satelliteHostChange)
    return
  end subroutine satelliteOrphanizeDestructor

  subroutine satelliteOrphanizeDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Trigger orphanization of a satellite.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite, &
          &                         propertyInactive
    use :: Histories       , only : history
    implicit none
    class    (nodeOperatorSatelliteOrphanize), intent(inout), target  :: self
    type     (treeNode                      ), intent(inout), target  :: node
    logical                                  , intent(inout)          :: interrupt
    procedure(interruptTask                 ), intent(inout), pointer :: functionInterrupt
    integer                                  , intent(in   )          :: propertyType
    class    (nodeComponentBasic            )               , pointer :: basic
    class    (nodeComponentSatellite        )               , pointer :: satellite
    type     (history               )                                 :: massBoundHistory
    logical                                                           :: exceedsHistoryTime

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    basic            => node     %basic           ()
    satellite        => node     %satellite       ()
    massBoundHistory =  satellite%boundMassHistory()
    if (massBoundHistory%exists()) then
       exceedsHistoryTime= basic           %time(                           ) &
            &             >=                                                  &
            &              massBoundHistory%time(size(massBoundHistory%time))
    else
       exceedsHistoryTime=.true.
    end if
    ! Test for orphanization.
    if     (                                         &
         &   .not.          satellite%isOrphan   ()  &
         &  .and.                                    &
         &                  node     %isSatellite()  &
         &  .and.                                    &
         &        associated(node    %mergeTarget  ) &
         &  .and.                                    &
         &   exceedsHistoryTime                      &
         & ) then
       ! Merging criterion met - trigger an interrupt.
       interrupt         =  .true.
       functionInterrupt => orphanizePerform
    end if
    return
  end subroutine satelliteOrphanizeDifferentialEvolution

  subroutine orphanizePerform(node,timeEnd)
    !!{
    Perform the orphanization the satellite.
    !!}
    use :: Display           , only : displayMessage    , displayVerbosity      , verbosityLevelInfo
    use :: Galacticus_Nodes  , only : nodeComponentBasic, nodeComponentSatellite, treeNode
    use :: ISO_Varying_String, only : operator(//)      , var_str               , varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    type            (treeNode              ), intent(inout), target      :: node
    double precision                        , intent(in   ), optional    :: timeEnd
    type            (treeNode              )               , pointer     :: nodeHost
    class           (nodeComponentBasic    )               , pointer     :: basic    , basicHost
    class           (nodeComponentSatellite)               , pointer     :: satellite
    !$GLC attributes unused :: timeEnd

    nodeHost  => node%mergeTarget
    satellite => node%satellite  ()
    call satellite%isOrphanSet(.true.)
    ! For satellite merge targets, step up through parents until an isolated host is found.
    do while (nodeHost%isSatellite())
       nodeHost => nodeHost%parent
    end do
    basic     => node    %basic()
    basicHost => nodeHost%basic()
    ! Trace the merge target progenitors back until one is found which exists at the time of the orphaned satellite.
    do while (basicHost%time() > basic%time())
       ! For satellite merge targets, step up through parents until an isolated host is found.
       do while (nodeHost%isSatellite())
          nodeHost => nodeHost%parent
       end do
       ! If a progenitor exists, move to it.
       if (associated(nodeHost%firstChild)) then
          nodeHost  => nodeHost%firstChild
          basicHost => nodeHost%basic     ()
       else
          ! No further progenitors exist, so stop here.
          exit
       end if
    end do
    ! Report.
    if (displayVerbosity() >= verbosityLevelInfo) then
       block
         type(varying_string) :: message
         message=var_str('Satellite node [')//node%index()//'] is being orphanized'
         if (associated(node%parent,nodeHost)) then
            message=message//' - remains in same host ['//nodeHost%index()//']'
         else
            message=message//' - moves from host ['//node%parent%index()//'] to host ['//nodeHost%index()//']'
         end if
         call displayMessage(message)
       end block
    end if
    ! Move to the new host. (If the new host is the same as the current host, do nothing.)
    if (.not.associated(node%parent,nodeHost)) then
       if (associated(node%parent)) call node%removeFromHost()
       node    %sibling        => nodeHost%firstSatellite
       node    %parent         => nodeHost
       nodeHost%firstSatellite => node
    end if
    return
  end subroutine orphanizePerform

  subroutine satelliteHostChange(self,node)
    !!{
    Handle cases where a satellite switches host node.
    !!}
    use :: Error             , only : Error_Report
    use :: Display           , only : displayMessage        , displayVerbosity, verbosityLevelInfo
    use :: Galacticus_Nodes  , only : nodeComponentSatellite, treeNode        , treeNodeLinkedList
    use :: ISO_Varying_String, only : operator(//)          , var_str         , varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentSatellite)               , pointer :: satellite
    type (treeNodeLinkedList    )               , pointer :: nodeStack, nodeNext, &
         &                                                   nodeNew
    type (treeNode              )               , pointer :: nodeWork , mergee

    select type (self)
    class is (nodeOperatorSatelliteOrphanize)
       satellite => node%satellite()
       if (satellite%isOrphan().and.associated(node%mergeTarget)) then
          if (displayVerbosity() >= verbosityLevelInfo) then
             block
               type (varying_string) :: message
               message=var_str('Satellite node [')//node%index()//'] will be orphanized due to host change'
               call displayMessage(message)
             end block
          end if
          ! Initialize a stack of nodes to allow us to process all mergees.
          allocate(nodeStack)
          nodeStack%node => node
          nodeStack%next => null()
          ! Process the stack.
          do while (associated(nodeStack))
             ! Pop a node from the stack.
             nodeWork => nodeStack%node
             nodeNext => nodeStack%next
             deallocate(nodeStack)
             nodeStack => nodeNext
             ! Push any mergees onto the stack.
             mergee => nodeWork%firstMergee
             do while (associated(mergee))
                ! Only push orphaned nodes onto the stack.
                satellite => mergee%satellite()
                if (satellite%isOrphan()) then
                   allocate(nodeNew)
                   nodeNew  %node => mergee
                   nodeNew  %next => nodeStack
                   nodeStack      => nodeNew
                end if
                mergee => mergee%siblingMergee
             end do
             ! Process the node.
             call orphanizePerform(nodeWork)
          end do
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteHostChange

  subroutine branchJumpPostProcess(self,node)
    !!{    
    For inter-tree node transfers and branch jumps, ensure that any orphaned mergees of the transferred node are transferred over
    to the new branch.
    !!}
    use :: Error             , only : Error_Report
    use :: Display           , only : displayMessage        , displayVerbosity, verbosityLevelInfo
    use :: Galacticus_Nodes  , only : nodeComponentSatellite, treeNode        , treeNodeLinkedList
    use :: ISO_Varying_String, only : operator(//)          , var_str         , varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), pointer :: node
    type (treeNode              )               , pointer :: mergee         , nodeWork
    type (treeNodeLinkedList    )               , pointer :: nodeStack      , nodeNext, &
         &                                                   nodeNew
    class(nodeComponentSatellite)               , pointer :: satelliteMergee

    select type (self)
    class is (nodeOperatorSatelliteOrphanize)
       nodeStack => null()
       mergee => node%firstMergee
       do while (associated(mergee))
          satelliteMergee => mergee%satellite()
          if (satelliteMergee%isOrphan()) then
             if (displayVerbosity() >= verbosityLevelInfo) then
                block
                  type(varying_string) :: message
                  message=var_str('Satellite node [')//mergee%index()//'] will be orphanized due to event'
                  call displayMessage(message)
                end block
             end if
             allocate(nodeNew)
             nodeNew  %node => mergee
             nodeNew  %next => nodeStack
             nodeStack      => nodeNew
          end if
          mergee => mergee%siblingMergee
       end do
       ! Process the stack.
       do while (associated(nodeStack))
          ! Pop a node from the stack.
          nodeWork => nodeStack%node
          nodeNext => nodeStack%next
          deallocate(nodeStack)
          nodeStack => nodeNext
          ! Push any mergees onto the stack.
          mergee => nodeWork%firstMergee
          do while (associated(mergee))
             ! Only push orphaned nodes onto the stack
             satelliteMergee => mergee%satellite()
             if (satelliteMergee%isOrphan()) then
                allocate(nodeNew)
                nodeNew  %node => mergee
                nodeNew  %next => nodeStack
                nodeStack      => nodeNew
             end if
             mergee => mergee%siblingMergee
          end do
          ! Process the node.
          call orphanizePerform(nodeWork)
       end do
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine branchJumpPostProcess
