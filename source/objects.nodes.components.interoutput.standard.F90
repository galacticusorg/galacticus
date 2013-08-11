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

!% Contains a module which implements the standard indices component.

module Node_Component_Inter_Output_Standard
  !% Implements the standard indices component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Inter_Output_Standard_Rate_Compute    , Node_Component_Inter_Output_Standard_Reset, &
       &    Node_Component_Inter_Output_Standard_Satellite_Merger

  !# <component>
  !#  <class>interOutput</class>
  !#  <name>standard</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>diskStarFormationRate</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar/gigaYear" comment="Disk star formation rate averaged over time between current and previous outputs."/>
  !#   </property>
  !#   <property>
  !#     <name>spheroidStarFormationRate</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar/gigaYear" comment="Spheroid star formation rate averaged over time between current and previous outputs."/>
  !#     <classDefault>0.0d0</classDefault>
  !#   </property>
  !#  </properties>
  !# </component>

contains


  !# <rateComputeTask>
  !#  <unitName>Node_Component_Inter_Output_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Inter_Output_Standard_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the exponential disk node mass rate of change.
    use Galacticus_Output_Times
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: thisNode
    logical                                       , intent(inout)          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentInterOutput    )               , pointer :: thisInterOutput
    class           (nodeComponentDisk           )               , pointer :: thisDisk
    class           (nodeComponentSpheroid       )               , pointer :: thisSpheroid
    class           (nodeComponentBasic          )               , pointer :: thisBasic
    double precision                                                       :: diskStarFormationRate, spheroidStarFormationRate, &
         &                                                                    timeCurrent          , timeOutputNext           , &
         &                                                                    timeOutputPrevious

    ! Return immediately if the standard inter-output component is not active.
    if (.not.defaultInteroutputComponent%standardIsActive()) return
    ! Get the disk and check that it is of our class.
    thisInterOutput => thisNode%interOutput()
    ! Get disk and spheroid star formation rates.
    thisDisk                  => thisNode    %disk             ()
    thisSpheroid              => thisNode    %spheroid         ()
    diskStarFormationRate     =  thisDisk    %starFormationRate()
    spheroidStarFormationRate =  thisSpheroid%starFormationRate()
    ! Find the time interval between previous and next outputs.
    thisBasic          => thisNode %basic()
    timeCurrent        =  thisBasic%time ()
    timeOutputPrevious =  Galacticus_Previous_Output_Time(timeCurrent)
    timeOutputNext     =  Galacticus_Next_Output_Time    (timeCurrent)
    ! Return if there is no next output.
    if (timeOutputNext     < 0.0d0) return
    ! Set previous time to zero if there is no previous output.
    if (timeOutputPrevious < 0.0d0) timeOutputPrevious=0.0d0
    ! Accumulate rates.
    call thisInterOutput%    diskStarFormationRateRate(    diskStarFormationRate/(timeOutputNext-timeOutputPrevious),interrupt,interruptProcedure)
    call thisInterOutput%spheroidStarFormationRateRate(spheroidStarFormationRate/(timeOutputNext-timeOutputPrevious),interrupt,interruptProcedure)
    return
  end subroutine Node_Component_Inter_Output_Standard_Rate_Compute

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Inter_Output_Standard_Reset</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Inter_Output_Standard_Reset(thisNode,iOutput,treeIndex,nodePassesFilter)
    !% Reset interoutput accumulated quantities.
    use Kind_Numbers
    implicit none
    type   (treeNode                ), intent(inout), pointer :: thisNode
    integer                          , intent(in   )          :: iOutput
    integer(kind=kind_int8          ), intent(in   )          :: treeIndex
    logical                          , intent(in   )          :: nodePassesFilter
    class  (nodeComponentInterOutput)               , pointer :: thisInterOutput

    ! Get the interoutput component and check it is of our class.
    thisInterOutput => thisNode%interOutput()
    select type (thisInterOutput)
       class is (nodeComponentInterOutputStandard)
       call thisInterOutput%    diskStarFormationRateSet(0.0d0)
       call thisInterOutput%spheroidStarFormationRateSet(0.0d0)
    end select
    return
  end subroutine Node_Component_Inter_Output_Standard_Reset

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Inter_Output_Standard_Satellite_Merger</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Inter_Output_Standard_Satellite_Merger(thisNode)
    !% Remove any inter-output quantities associated with {\tt thisNode} and add them to the merge target.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    implicit none
    type (treeNode                ), intent(inout), pointer :: thisNode
    type (treeNode                )               , pointer :: hostNode
    class(nodeComponentInterOutput)               , pointer :: hostInterOutput, thisInterOutput

    ! Get the inter-output component.
    thisInterOutput => thisNode%interOutput()
    ! Ensure that it is of the standard class.
    select type (thisInterOutput)
    class is (nodeComponentInterOutputStandard)
       ! Find the node to merge with.
       hostNode        => thisNode%mergesWith ()
       hostInterOutput => hostNode%interOutput()
       ! Move the star formation rates from secondary to primary.
       select case (thisMergerStarsMoveTo)
       case (movesToDisk    )
          call hostInterOutput%    diskStarFormationRateSet(                                             &
               &                                             hostInterOutput%    diskStarFormationRate() &
               &                                            +thisInterOutput%    diskStarFormationRate() &
               &                                           )
          call hostInterOutput%spheroidStarFormationRateSet(                                             &
               &                                             hostInterOutput%spheroidStarFormationRate() &
               &                                            +thisInterOutput%spheroidStarFormationRate() &
               &                                           )
       case (movesToSpheroid)
       case default
          call Galacticus_Error_Report('Node_Component_Inter_Output_Standard_Satellite_Merger','unrecognized movesTo descriptor')
       end select
       ! Zero rates in the secondary,
       call thisInterOutput%    diskStarFormationRateSet(                                             &
            &                                             0.0d0                                       &
            &                                           )
       call thisInterOutput%spheroidStarFormationRateSet(                                             &
            &                                             0.0d0                                       &
            &                                           )
       ! Move star formation rates within the host if necessary.
       select case (thisHostStarsMoveTo)
       case (movesToDisk)
          call hostInterOutput%    diskStarFormationRateSet(                                             &
               &                                             hostInterOutput%    diskStarFormationRate() &
               &                                            +hostInterOutput%spheroidStarFormationRate() &
               &                                           )
          call hostInterOutput%spheroidStarFormationRateSet(                                             &
               &                                             0.0d0                                       &
               &                                           )
       case (movesToSpheroid)
          call hostInterOutput%spheroidStarFormationRateSet(                                             &
               &                                             hostInterOutput%spheroidStarFormationRate() &
               &                                            +hostInterOutput%    diskStarFormationRate() &
               &                                           )
          call hostInterOutput%    diskStarFormationRateSet(                                             &
               &                                             0.0d0                                       &
               &                                           )
       case (doesNotMove)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('Node_Component_Inter_Output_Standard_Satellite_Merger','unrecognized movesTo descriptor')
       end select
    end select
    return
  end subroutine Node_Component_Inter_Output_Standard_Satellite_Merger

end module Node_Component_Inter_Output_Standard
