!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements the standard galaxy age statistics component.

module Node_Component_Age_Statistics_Standard
  !% Implements the standard galaxy age statistics component.
  implicit none
  private
  public :: Node_Component_Age_Statistics_Standard_Scale_Set        , Node_Component_Age_Statistics_Standard_Rate_Compute, &
       &    Node_Component_Age_Statistics_Standard_Satellite_Merging, Node_Component_Age_Statistics_Standard_Inactive    , &
       &    Node_Component_Age_Statistics_Standard_Initialize

  !# <component>
  !#  <class>ageStatistics</class>
  !#  <name>standard</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>diskTimeWeightedIntegratedSFR</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar*gigaYear" comment="Time-weighted integral over disk star formation rate."/>
  !#   </property>
  !#   <property>
  !#     <name>spheroidTimeWeightedIntegratedSFR</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar*gigaYear" comment="Time-weighted integral over spheroid star formation rate."/>
  !#   </property>
  !#   <property>
  !#     <name>diskIntegratedSFR</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar" comment="Integral over disk star formation rate."/>
  !#   </property>
  !#   <property>
  !#     <name>spheroidIntegratedSFR</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <classDefault>0.0d0</classDefault>
  !#     <output unitsInSI="massSolar" comment="Integral over spheroid star formation rate."/>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Record of whether variables in this component are inactive.
  logical :: ageStatisticsStandardIsInactive
  
contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Age_Statistics_Standard_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Age_Statistics_Standard_Initialize(globalParameters_)
    !% Initializes the tree node standard disk methods module.
    use Input_Parameters
    use Galacticus_Nodes, only : defaultAgeStatisticsComponent
    implicit none
    type(inputParameters), intent(inout) :: globalParameters_

    if (defaultAgeStatisticsComponent%standardIsActive()) then
       !# <inputParameter>
       !#   <name>ageStatisticsStandardIsInactive</name>
       !#   <cardinality>1</cardinality>
       !#   <defaultValue>.false.</defaultValue>
       !#   <description>Specifies whether or not the variables of the standard age statistics component are inactive (i.e. do not appear in any ODE being solved).</description>
       !#   <source>globalParameters_</source>
       !#   <type>boolean</type>
       !# </inputParameter>
    end if
    return
  end subroutine Node_Component_Age_Statistics_Standard_Initialize

  !# <inactiveSetTask>
  !#  <unitName>Node_Component_Age_Statistics_Standard_Inactive</unitName>
  !# </inactiveSetTask>
  subroutine Node_Component_Age_Statistics_Standard_Inactive(node)
    !% Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    use Stellar_Luminosities_Structure
    use Galacticus_Nodes              , only : treeNode, nodeComponentAgeStatistics, nodeComponentAgeStatisticsStandard
    implicit none
    type (treeNode                  ), intent(inout), pointer :: node
    class(nodeComponentAgeStatistics)               , pointer :: ageStatistics
    
    ! Get the age statistics component.
    ageStatistics => node%ageStatistics()
    ! Check if an standard disk component exists.
    select type (ageStatistics)
    class is (nodeComponentAgeStatisticsStandard)
       if (ageStatisticsStandardIsInactive) then
          call ageStatistics%    diskTimeWeightedIntegratedSFRInactive()
          call ageStatistics%spheroidTimeWeightedIntegratedSFRInactive()
          call ageStatistics%                diskIntegratedSFRInactive()
          call ageStatistics%            spheroidIntegratedSFRInactive()
       end if
    end select
    return
  end subroutine Node_Component_Age_Statistics_Standard_Inactive
  
  !# <scaleSetTask>
  !#  <unitName>Node_Component_Age_Statistics_Standard_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Age_Statistics_Standard_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use Galacticus_Nodes, only : treeNode, nodeComponentAgeStatistics, nodeComponentAgeStatisticsStandard, nodeComponentDisk, nodeComponentSpheroid
    implicit none
    type            (treeNode                  ), intent(inout), pointer :: node
    class           (nodeComponentAgeStatistics)               , pointer :: ageStatistics
    class           (nodeComponentDisk         )               , pointer :: disk
    class           (nodeComponentSpheroid     )               , pointer :: spheroid
    double precision                            , parameter              :: massMinimum    =1.0d0
    double precision                            , parameter              :: timeScale      =1.0d-3
    double precision                                                     :: mass

    ! Get the age statistics component.
    ageStatistics => node%ageStatistics()
    ! Check if component is of standard class.
    select type (ageStatistics)
    class is (nodeComponentAgeStatisticsStandard)
       ! Get disk and spheroid components.
       disk     => node%disk    ()
       spheroid => node%spheroid()       
       ! Set scale for masses.
       mass   = disk%massGas    ()+spheroid%massGas    () &
            &  +disk%massStellar()+spheroid%massStellar()
       call ageStatistics%    diskTimeWeightedIntegratedSFRScale(max(mass,massMinimum)*timeScale)
       call ageStatistics%spheroidTimeWeightedIntegratedSFRScale(max(mass,massMinimum)*timeScale)
       call ageStatistics%                diskIntegratedSFRScale(max(mass,massMinimum)          )
       call ageStatistics%            spheroidIntegratedSFRScale(max(mass,massMinimum)          )
    end select
    return
  end subroutine Node_Component_Age_Statistics_Standard_Scale_Set

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Age_Statistics_Standard_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Age_Statistics_Standard_Rate_Compute(node,odeConverged,interrupt,interruptProcedure,propertyType)
    !% Compute the exponential disk node mass rate of change.
    use Galacticus_Nodes, only : treeNode             , nodeComponentAgeStatistics, nodeComponentAgeStatisticsStandard, nodeComponentDisk            , &
         &                       nodeComponentSpheroid, nodeComponentBasic        , interruptTask                     , defaultAgeStatisticsComponent, &
         &                       propertyTypeActive   , propertyTypeInactive      , propertyTypeAll
    use Output_Times
    implicit none
    type            (treeNode                    ), intent(inout), pointer :: node
    logical                                       , intent(in   )          :: odeConverged
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: interruptProcedure
    integer                                       , intent(in   )          :: propertyType
    class           (nodeComponentAgeStatistics  )               , pointer :: ageStatistics
    class           (nodeComponentDisk           )               , pointer :: disk
    class           (nodeComponentSpheroid       )               , pointer :: spheroid
    class           (nodeComponentBasic          )               , pointer :: basic
    double precision                                                       :: time                     , diskStarFormationRate, &
         &                                                                    spheroidStarFormationRate
    logical                                                                :: isGeneric
    !GCC$ attributes unused :: odeConverged
    
    ! Return immediately if the standard age statistics component is not active.
    if (.not.defaultAgeStatisticsComponent%standardIsActive()) return
    ! Get the age statistics component.
    ageStatistics => node%ageStatistics()
    select type (ageStatistics)
    type is (nodeComponentAgeStatistics)
       isGeneric=.true.
    class default
       isGeneric=.false.
    end select
    ! Return immediately if the wrong property type is required.
    if     (                                                                                     &
         &   (                                                                                   &
         &     (propertyType == propertyTypeActive   .and.      ageStatisticsStandardIsInactive) &
         &    .or.                                                                               &
         &     (propertyType == propertyTypeInactive .and. .not.ageStatisticsStandardIsInactive) &
         &    .or.                                                                               &
         &      propertyType == propertyTypeAll                                                  &
         &   )                                                                                   &
         &  .and.                                                                                &
         &                                                 .not.isGeneric                        &
         & ) return
    ! Get the star formation rates.
    disk                      => node    %disk             ()
    spheroid                  => node    %spheroid         ()
    diskStarFormationRate     =  disk    %starFormationRate()
    spheroidStarFormationRate =  spheroid%starFormationRate()
    ! Find the current cosmic time.
    basic => node %basic() 
    time  =  basic%time ()
    ! Accumulate rates.
    call ageStatistics%    diskTimeWeightedIntegratedSFRRate(    diskStarFormationRate*time,interrupt,interruptProcedure)
    call ageStatistics%spheroidTimeWeightedIntegratedSFRRate(spheroidStarFormationRate*time,interrupt,interruptProcedure)
    call ageStatistics%                diskIntegratedSFRRate(    diskStarFormationRate     ,interrupt,interruptProcedure)
    call ageStatistics%            spheroidIntegratedSFRRate(spheroidStarFormationRate     ,interrupt,interruptProcedure)
    return
  end subroutine Node_Component_Age_Statistics_Standard_Rate_Compute

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Age_Statistics_Standard_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Remnant_Compute</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Age_Statistics_Standard_Satellite_Merging(node)
    !% Remove any age statistics quantities associated with {\normalfont \ttfamily node} and add them to the merge target.
    use Satellite_Merging_Remnant_Properties
    use Satellite_Merging_Mass_Movements
    use Galacticus_Error
    use Galacticus_Nodes                    , only : treeNode, nodeComponentAgeStatistics, nodeComponentAgeStatisticsStandard
    implicit none
    type   (treeNode                  ), intent(inout), pointer :: node
    type   (treeNode                  )               , pointer :: nodeHost
    class  (nodeComponentAgeStatistics)               , pointer :: ageStatistics, ageStatisticsHost
 
    ! Get the inter-output component.
    ageStatistics => node%ageStatistics()
    ! Ensure that it is of the standard class.
    select type (ageStatistics)
    class is (nodeComponentAgeStatisticsStandard)
       ! Find the node to merge with.
       nodeHost          => node    %mergesWith   (                 )
       ageStatisticsHost => nodeHost%ageStatistics(autoCreate=.true.)
       ! Move the star formation rates from secondary to primary.
       select case (destinationStarsSatellite)
       case (destinationMergerDisk    )
          call ageStatisticsHost%    diskTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       ageStatisticsHost%    diskTimeWeightedIntegratedSFR() &
               &                                                      +ageStatistics    %    diskTimeWeightedIntegratedSFR() &
               &                                                     )
          call ageStatisticsHost%spheroidTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       ageStatisticsHost%spheroidTimeWeightedIntegratedSFR() &
               &                                                      +ageStatistics    %spheroidTimeWeightedIntegratedSFR() &
               &                                                     )
          call ageStatisticsHost%                diskIntegratedSFRSet(                                                       &
               &                                                       ageStatisticsHost%                diskIntegratedSFR() &
               &                                                      +ageStatistics    %                diskIntegratedSFR() &
               &                                                     )
          call ageStatisticsHost%            spheroidIntegratedSFRSet(                                                       &
               &                                                       ageStatisticsHost%            spheroidIntegratedSFR() &
               &                                                      +ageStatistics    %            spheroidIntegratedSFR() &
               &                                                     )
       case (destinationMergerSpheroid)
       case default
          call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       ! Zero rates in the secondary,
       call ageStatistics%    diskTimeWeightedIntegratedSFRSet(                                                          &
            &                                                   0.0d0                                                    &
            &                                                 )
       call ageStatistics%spheroidTimeWeightedIntegratedSFRSet(                                                          &
            &                                                   0.0d0                                                    &
            &                                                 )
       call ageStatistics%                diskIntegratedSFRSet(                                                          &
            &                                                   0.0d0                                                    &
            &                                                 )
       call ageStatistics%            spheroidIntegratedSFRSet(                                                          &
            &                                                   0.0d0                                                    &
            &                                                 )
       ! Move star formation rates within the host if necessary.
       select case (destinationStarsHost)
       case (destinationMergerDisk)
          call ageStatisticsHost%    diskTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       ageStatisticsHost%    diskTimeWeightedIntegratedSFR() &
               &                                                      +ageStatisticsHost%spheroidTimeWeightedIntegratedSFR() &
               &                                                     )
          call ageStatisticsHost%spheroidTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       0.0d0                                                 &
               &                                                     )
          call ageStatisticsHost%                diskIntegratedSFRSet(                                                       &
               &                                                       ageStatisticsHost%                diskIntegratedSFR() &
               &                                                      +ageStatisticsHost%            spheroidIntegratedSFR() &
               &                                                     )
          call ageStatisticsHost%            spheroidIntegratedSFRSet(                                                       &
               &                                                       0.0d0                                                 &
               &                                                     )
       case (destinationMergerSpheroid)
          call ageStatisticsHost%spheroidTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       ageStatisticsHost%spheroidTimeWeightedIntegratedSFR() &
               &                                                      +ageStatisticsHost%    diskTimeWeightedIntegratedSFR() &
               &                                                     )
          call ageStatisticsHost%    diskTimeWeightedIntegratedSFRSet(                                                       &
               &                                                       0.0d0                                                 &
               &                                                     )
          call ageStatisticsHost%            spheroidIntegratedSFRSet(                                                       &
               &                                                       ageStatisticsHost%            spheroidIntegratedSFR() &
               &                                                      +ageStatisticsHost%                diskIntegratedSFR() &
               &                                                     )
          call ageStatisticsHost%                diskIntegratedSFRSet(                                                       &
               &                                                       0.0d0                                                 &
               &                                                     )
       case (destinationMergerUnmoved)
          ! Do nothing.
       case default
          call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select       
    end select
    return
  end subroutine Node_Component_Age_Statistics_Standard_Satellite_Merging

end module Node_Component_Age_Statistics_Standard
