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

!% Contains a module which implements the standard galaxy age statistics component.

module Node_Component_Age_Statistics_Standard
  !% Implements the standard galaxy age statistics component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Age_Statistics_Standard_Scale_Set       , Node_Component_Age_Statistics_Standard_Rate_Compute, &
       &    Node_Component_Age_Statistics_Standard_Satellite_Merger

  !# <component>
  !#  <class>ageStatistics</class>
  !#  <name>standard</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>diskTimeWeightedIntegratedSFR</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar*gigaYear" comment="Time-weighted integral over disk star formation rate."/>
  !#   </property>
  !#   <property>
  !#     <name>spheroidTimeWeightedIntegratedSFR</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar*gigaYear" comment="Time-weighted integral over spheroid star formation rate."/>
  !#   </property>
  !#   <property>
  !#     <name>diskIntegratedSFR</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar" comment="Integral over disk star formation rate."/>
  !#   </property>
  !#   <property>
  !#     <name>spheroidIntegratedSFR</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar" comment="Integral over spheroid star formation rate."/>
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Age_Statistics_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Age_Statistics_Standard_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type            (treeNode                  ), intent(inout), pointer :: thisNode
    class           (nodeComponentAgeStatistics)               , pointer :: thisAgeStatistics
    class           (nodeComponentDisk         )               , pointer :: thisDisk
    class           (nodeComponentSpheroid     )               , pointer :: thisSpheroid
    double precision                            , parameter              :: massMinimum    =1.0d0
    double precision                            , parameter              :: timeScale      =1.0d-3
    double precision                                                     :: mass

    ! Get the age statistics component.
    thisAgeStatistics => thisNode%ageStatistics()
    ! Check if component is of standard class.
    select type (thisAgeStatistics)
    class is (nodeComponentAgeStatisticsStandard)
       ! Get disk and spheroid components.
       thisDisk     => thisNode%disk    ()
       thisSpheroid => thisNode%spheroid()       
       ! Set scale for masses.
       mass   = thisDisk%massGas    ()+thisSpheroid%massGas    () &
            &  +thisDisk%massStellar()+thisSpheroid%massStellar()
       call thisAgeStatistics%    diskTimeWeightedIntegratedSFRScale(max(mass,massMinimum)*timeScale)
       call thisAgeStatistics%spheroidTimeWeightedIntegratedSFRScale(max(mass,massMinimum)*timeScale)
       call thisAgeStatistics%                diskIntegratedSFRScale(max(mass,massMinimum)          )
       call thisAgeStatistics%            spheroidIntegratedSFRScale(max(mass,massMinimum)          )
    end select
    return
  end subroutine Node_Component_Age_Statistics_Standard_Scale_Set

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Age_Statistics_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Age_Statistics_Standard_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the exponential disk node mass rate of change.
    use Galacticus_Nodes
    use Galacticus_Output_Times
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: thisNode
    logical                                       , intent(inout)          :: interrupt
    procedure       (Interrupt_Procedure_Template), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentAgeStatistics  )               , pointer :: thisAgeStatistics
    class           (nodeComponentDisk           )               , pointer :: thisDisk
    class           (nodeComponentSpheroid       )               , pointer :: thisSpheroid
    class           (nodeComponentBasic          )               , pointer :: thisBasic
    double precision                                                       :: time                     , diskStarFormationRate, &
         &                                                                    spheroidStarFormationRate

    ! Return immediately if the standard age statistics component is not active.
    if (.not.defaultAgeStatisticsComponent%standardIsActive()) return
    ! Get the age statistics component.
    thisAgeStatistics => thisNode%ageStatistics()
    ! Get the star formation rates.
    thisDisk                  => thisNode    %disk             ()
    thisSpheroid              => thisNode    %spheroid         ()
    diskStarFormationRate     =  thisDisk    %starFormationRate()
    spheroidStarFormationRate =  thisSpheroid%starFormationRate()
    ! Find the current cosmic time.
    thisBasic => thisNode %basic() 
    time      =  thisBasic%time ()
    ! Accumulate rates.
    call thisAgeStatistics%    diskTimeWeightedIntegratedSFRRate(    diskStarFormationRate*time,interrupt,interruptProcedure)
    call thisAgeStatistics%spheroidTimeWeightedIntegratedSFRRate(spheroidStarFormationRate*time,interrupt,interruptProcedure)
    call thisAgeStatistics%                diskIntegratedSFRRate(    diskStarFormationRate     ,interrupt,interruptProcedure)
    call thisAgeStatistics%            spheroidIntegratedSFRRate(spheroidStarFormationRate     ,interrupt,interruptProcedure)
    return
  end subroutine Node_Component_Age_Statistics_Standard_Rate_Compute

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Age_Statistics_Standard_Satellite_Merger</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Age_Statistics_Standard_Satellite_Merger(thisNode)
    !% Remove any age statistics quantities associated with {\tt thisNode} and add them to the merge target.
    use Satellite_Merging_Mass_Movements_Descriptors
    use Galacticus_Error
    implicit none
    type (treeNode                  ), intent(inout), pointer :: thisNode
    type (treeNode                  )               , pointer :: hostNode
    class(nodeComponentAgeStatistics)               , pointer :: thisAgeStatistics, hostAgeStatistics

    ! Get the inter-output component.
    thisAgeStatistics => thisNode%ageStatistics()
    ! Ensure that it is of the standard class.
    select type (thisAgeStatistics)
    class is (nodeComponentAgeStatisticsStandard)
       ! Find the node to merge with.
       hostNode          => thisNode%mergesWith   ()
       hostAgeStatistics => hostNode%ageStatistics()
       ! Move the star formation rates from secondary to primary.
       select case (thisMergerStarsMoveTo)
       case (movesToDisk    )
          call hostAgeStatistics%    diskTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       hostAgeStatistics%    diskTimeWeightedIntegratedSFR() &
               &                                                      +thisAgeStatistics%    diskTimeWeightedIntegratedSFR() &
               &                                                     )
          call hostAgeStatistics%spheroidTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       hostAgeStatistics%spheroidTimeWeightedIntegratedSFR() &
               &                                                      +thisAgeStatistics%spheroidTimeWeightedIntegratedSFR() &
               &                                                     )
          call hostAgeStatistics%                diskIntegratedSFRSet(                                                       &
               &                                                       hostAgeStatistics%                diskIntegratedSFR() &
               &                                                      +thisAgeStatistics%                diskIntegratedSFR() &
               &                                                     )
          call hostAgeStatistics%            spheroidIntegratedSFRSet(                                                       &
               &                                                       hostAgeStatistics%            spheroidIntegratedSFR() &
               &                                                      +thisAgeStatistics%            spheroidIntegratedSFR() &
               &                                                     )
       case (movesToSpheroid)
       case default
          call Galacticus_Error_Report('Node_Component_Age_Statistics_Standard_Satellite_Merger','unrecognized movesTo descriptor')
       end select
       ! Zero rates in the secondary,
       call thisAgeStatistics%    diskTimeWeightedIntegratedSFRSet(                                                          &
            &                                                       0.0d0                                                    &
            &                                                     )
       call thisAgeStatistics%spheroidTimeWeightedIntegratedSFRSet(                                                          &
            &                                                       0.0d0                                                    &
            &                                                     )
       call thisAgeStatistics%                diskIntegratedSFRSet(                                                          &
            &                                                       0.0d0                                                    &
            &                                                     )
       call thisAgeStatistics%            spheroidIntegratedSFRSet(                                                          &
            &                                                       0.0d0                                                    &
            &                                                     )
       ! Move star formation rates within the host if necessary.
       select case (thisHostStarsMoveTo)
       case (movesToDisk)
          call hostAgeStatistics%    diskTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       hostAgeStatistics%    diskTimeWeightedIntegratedSFR() &
               &                                                      +hostAgeStatistics%spheroidTimeWeightedIntegratedSFR() &
               &                                                     )
          call hostAgeStatistics%spheroidTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       0.0d0                                                 &
               &                                                     )
          call hostAgeStatistics%                diskIntegratedSFRSet(                                                       &
               &                                                       hostAgeStatistics%                diskIntegratedSFR() &
               &                                                      +hostAgeStatistics%            spheroidIntegratedSFR() &
               &                                                     )
          call hostAgeStatistics%            spheroidIntegratedSFRSet(                                                       &
               &                                                       0.0d0                                                 &
               &                                                     )
       case (movesToSpheroid)
          call hostAgeStatistics%spheroidTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       hostAgeStatistics%spheroidTimeWeightedIntegratedSFR() &
               &                                                      +hostAgeStatistics%    diskTimeWeightedIntegratedSFR() &
               &                                                     )
          call hostAgeStatistics%    diskTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       0.0d0                                                 &
               &                                                     )
          call hostAgeStatistics%            spheroidIntegratedSFRSet(                                                       &
               &                                                       hostAgeStatistics%            spheroidIntegratedSFR() &
               &                                                      +hostAgeStatistics%                diskIntegratedSFR() &
               &                                                     )
          call hostAgeStatistics%                diskIntegratedSFRSet(                                                       &
               &                                                       0.0d0                                                 &
               &                                                     )
       case (doesNotMove)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('Node_Component_Age_Statistics_Standard_Satellite_Merger','unrecognized movesTo descriptor')
       end select       
    end select
    return
  end subroutine Node_Component_Age_Statistics_Standard_Satellite_Merger

end module Node_Component_Age_Statistics_Standard
