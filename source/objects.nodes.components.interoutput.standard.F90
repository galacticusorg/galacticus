!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements the standard indices component.

module Node_Component_Inter_Output_Standard
  !% Implements the standard indices component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Inter_Output_Standard_Rate_Compute     , Node_Component_Inter_Output_Standard_Reset    , &
       &    Node_Component_Inter_Output_Standard_Satellite_Merging, Node_Component_Inter_Output_Standard_Scale_Set

  !# <component>
  !#  <class>interOutput</class>
  !#  <name>standard</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>diskStarFormationRate</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar/gigaYear" comment="Disk star formation rate averaged over time between current and previous outputs."/>
  !#   </property>
  !#   <property>
  !#     <name>spheroidStarFormationRate</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar/gigaYear" comment="Spheroid star formation rate averaged over time between current and previous outputs."/>
  !#     <classDefault>0.0d0</classDefault>
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Inter_Output_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Inter_Output_Standard_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    implicit none
    type            (treeNode                ), intent(inout), pointer :: node
    class           (nodeComponentInterOutput)               , pointer :: interOutput
    class           (nodeComponentDisk       )               , pointer :: disk
    class           (nodeComponentSpheroid   )               , pointer :: spheroid
    double precision                          , parameter              :: massMinimum    =1.0d0
    double precision                          , parameter              :: timeScale      =1.0d0
    double precision                                                   :: mass

    ! Get the interoutput component.
    interOutput => node%interOutput()
    ! Check if component is of standard class.
    select type (interOutput)
    class is (nodeComponentInterOutputStandard)
       ! Get disk and spheroid components.
       disk     => node%disk    ()
       spheroid => node%spheroid()       
       ! Set scale for masses.
       mass   = disk%massGas    ()+spheroid%massGas    () &
            &  +disk%massStellar()+spheroid%massStellar()
       call interOutput%    diskStarFormationRateScale(max(mass,massMinimum)/timeScale)
       call interOutput%spheroidStarFormationRateScale(max(mass,massMinimum)/timeScale)
    end select
    return
  end subroutine Node_Component_Inter_Output_Standard_Scale_Set

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Inter_Output_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Inter_Output_Standard_Rate_Compute(node,odeConverged,interrupt,interruptProcedure)
    !% Compute the exponential disk node mass rate of change.
    use Galacticus_Output_Times
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: node
    logical                                       , intent(in   )          :: odeConverged
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentInterOutput    )               , pointer :: interOutput
    class           (nodeComponentDisk           )               , pointer :: disk
    class           (nodeComponentSpheroid       )               , pointer :: spheroid
    class           (nodeComponentBasic          )               , pointer :: basic
    double precision                                                       :: diskStarFormationRate, spheroidStarFormationRate, &
         &                                                                    timeCurrent          , timeOutputNext           , &
         &                                                                    timeOutputPrevious
    !GCC$ attributes unused :: odeConverged
    
    ! Return immediately if the standard inter-output component is not active.
    if (.not.defaultInteroutputComponent%standardIsActive()) return
    ! Get the disk and check that it is of our class.
    interOutput => node%interOutput()
    ! Get disk and spheroid star formation rates.
    disk                      => node    %disk             ()
    spheroid                  => node    %spheroid         ()
    diskStarFormationRate     =  disk    %starFormationRate()
    spheroidStarFormationRate =  spheroid%starFormationRate()
    ! Find the time interval between previous and next outputs.
    basic              => node %basic()
    timeCurrent        =  basic%time ()
    timeOutputPrevious =  Galacticus_Previous_Output_Time(timeCurrent)
    timeOutputNext     =  Galacticus_Next_Output_Time    (timeCurrent)
    ! Return if there is no next output.
    if (timeOutputNext     < 0.0d0) return
    ! Set previous time to zero if there is no previous output.
    if (timeOutputPrevious < 0.0d0) timeOutputPrevious=0.0d0
    ! Accumulate rates.
    call interOutput%    diskStarFormationRateRate(    diskStarFormationRate/(timeOutputNext-timeOutputPrevious),interrupt,interruptProcedure)
    call interOutput%spheroidStarFormationRateRate(spheroidStarFormationRate/(timeOutputNext-timeOutputPrevious),interrupt,interruptProcedure)
    return
  end subroutine Node_Component_Inter_Output_Standard_Rate_Compute

  !# <mergerTreeExtraOutputTask>
  !#  <unitName>Node_Component_Inter_Output_Standard_Reset</unitName>
  !# </mergerTreeExtraOutputTask>
  subroutine Node_Component_Inter_Output_Standard_Reset(node,iOutput,treeIndex,nodePassesFilter)
    !% Reset interoutput accumulated quantities.
    use, intrinsic :: ISO_C_Binding
    use Kind_Numbers
    implicit none
    type   (treeNode                ), intent(inout), pointer :: node
    integer(c_size_t                ), intent(in   )          :: iOutput
    integer(kind=kind_int8          ), intent(in   )          :: treeIndex
    logical                          , intent(in   )          :: nodePassesFilter
    class  (nodeComponentInterOutput)               , pointer :: interOutput
    !GCC$ attributes unused :: iOutput, nodePassesFilter, treeIndex
    
    ! Get the interoutput component and check it is of our class.
    interOutput => node%interOutput()
    select type (interOutput)
       class is (nodeComponentInterOutputStandard)
       call interOutput%    diskStarFormationRateSet(0.0d0)
       call interOutput%spheroidStarFormationRateSet(0.0d0)
    end select
    return
  end subroutine Node_Component_Inter_Output_Standard_Reset

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Inter_Output_Standard_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Inter_Output_Standard_Satellite_Merging(node)
    !% Remove any inter-output quantities associated with {\normalfont \ttfamily node} and add them to the merge target.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    implicit none
    type (treeNode                ), intent(inout), pointer :: node
    type (treeNode                )               , pointer :: nodeHost
    class(nodeComponentInterOutput)               , pointer :: interOutputHost, interOutput

    ! Get the inter-output component.
    interOutput => node%interOutput()
    ! Ensure that it is of the standard class.
    select type (interOutput)
    class is (nodeComponentInterOutputStandard)
       ! Find the node to merge with.
       nodeHost        => node%mergesWith ()
       interOutputHost => nodeHost%interOutput()
       ! Move the star formation rates from secondary to primary.
       select case (thisMergerStarsMoveTo)
       case (movesToDisk    )
          call interOutputHost%    diskStarFormationRateSet(                                             &
               &                                             interOutputHost%    diskStarFormationRate() &
               &                                            +interOutput    %    diskStarFormationRate() &
               &                                           )
          call interOutputHost%spheroidStarFormationRateSet(                                             &
               &                                             interOutputHost%spheroidStarFormationRate() &
               &                                            +interOutput    %spheroidStarFormationRate() &
               &                                           )
       case (movesToSpheroid)
       case default
          call Galacticus_Error_Report('Node_Component_Inter_Output_Standard_Satellite_Merging','unrecognized movesTo descriptor')
       end select
       ! Zero rates in the secondary,
       call interOutput%    diskStarFormationRateSet(                                             &
            &                                         0.0d0                                       &
            &                                       )
       call interOutput%spheroidStarFormationRateSet(                                             &
            &                                         0.0d0                                       &
            &                                       )
       ! Move star formation rates within the host if necessary.
       select case (thisHostStarsMoveTo)
       case (movesToDisk)
          call interOutputHost%    diskStarFormationRateSet(                                             &
               &                                             interOutputHost%    diskStarFormationRate() &
               &                                            +interOutputHost%spheroidStarFormationRate() &
               &                                           )
          call interOutputHost%spheroidStarFormationRateSet(                                             &
               &                                             0.0d0                                       &
               &                                           )
       case (movesToSpheroid)
          call interOutputHost%spheroidStarFormationRateSet(                                             &
               &                                             interOutputHost%spheroidStarFormationRate() &
               &                                            +interOutputHost%    diskStarFormationRate() &
               &                                           )
          call interOutputHost%    diskStarFormationRateSet(                                             &
               &                                             0.0d0                                       &
               &                                           )
       case (doesNotMove)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('Node_Component_Inter_Output_Standard_Satellite_Merging','unrecognized movesTo descriptor')
       end select
    end select
    return
  end subroutine Node_Component_Inter_Output_Standard_Satellite_Merging

end module Node_Component_Inter_Output_Standard
