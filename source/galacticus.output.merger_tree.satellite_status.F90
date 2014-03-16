!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which handles outputting of satellite status (i.e. whether orphaned or not) to the \glc\ output file.

module Galacticus_Output_Trees_Satellite_Status
  !% Handles outputting of satellite status (i.e. whether orphaned or not) to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Satellite_Status, Galacticus_Output_Tree_Satellite_Status_Property_Count,&
       & Galacticus_Output_Tree_Satellite_Status_Names

  ! Number of host properties.
  integer :: satelliteStatusPropertyCount

  ! Flag indicating whether or not satellite host information is to be output.
  logical :: outputSatelliteStatus

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputSatelliteStatusInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Satellite_Status_Initialize
    !% Initializes the module by determining whether or not satellite host data should be output.
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Nodes
    use ISO_Varying_String
    implicit none

    if (.not.outputSatelliteStatusInitialized) then
       !$omp critical(Galacticus_Output_Tree_Satellite_Status_Initialize)
       if (.not.outputSatelliteStatusInitialized) then
          !@ <inputParameter>
          !@   <name>outputSatelliteStatus</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not satellite status (i.e. whether the satellite is orphaned or not) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputSatelliteStatus',outputSatelliteStatus,defaultValue=.false.)
          ! Count number of properties to output and check that required properties are gettable.
          satelliteStatusPropertyCount=0
          if (outputSatelliteStatus) then
             satelliteStatusPropertyCount=satelliteStatusPropertyCount+1
             if     (                                                                                                                    &
                  &  .not.                                                                                                               &
                  &       (                                                                                                              &
                  &        defaultSatelliteComponent% boundMassHistoryIsGettable()                                                       &
                  &       )                                                                                                              &
                  & ) call Galacticus_Error_Report                                                                                       &
                  &        (                                                                                                             &
                  &         'Galacticus_Output_Tree_Satellite_Status_Initialize'                                                      ,  &
                  &         'this method requires that boundMassHistory property must be gettable for the satellite component.'       // &
                  &         Galacticus_Component_List(                                                                                   &
                  &                                   'satellite'                                                                     ,  &
                  &                                   defaultSatelliteComponent%boundMassHistoryAttributeMatch(requireGettable=.true.)   &
                  &                                  )                                                                                   &
                  &        )
           end if
          ! Flag that module is now initialized.
          outputSatelliteStatusInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Satellite_Status_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Status_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Status_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Status</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Satellite_Status_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of satellite orbital extremum properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: time
    integer                                 , intent(inout)          :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout)          :: doublePropertyComments , doublePropertyNames   , &
         &                                                              integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout)          :: doublePropertyUnitsSI  , integerPropertyUnitsSI

    ! Initialize the module.
    call Galacticus_Output_Tree_Satellite_Status_Initialize

    ! Return property names if we are outputting satellite host data.
    if (outputSatelliteStatus) then
       integerProperty=integerProperty+1
       !@ <outputProperty>
       !@   <name>satelliteStatus</name>
       !@   <datatype>integer</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Satellite status flag (0=not a satellite; 1=satellite with halo; 2=orphaned satellite).</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       integerPropertyNames   (integerProperty)='satelliteStatus'
       integerPropertyComments(integerProperty)="Satellite status flag (0=not a satellite; 1=satellite with halo; 2=orphaned satellite)."
       integerPropertyUnitsSI (integerProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Status_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Status_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Status</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Satellite_Status_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of satellite host properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Satellite_Status_Initialize

    ! Increment property count.
    integerPropertyCount=integerPropertyCount+satelliteStatusPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Satellite_Status_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Status</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Status</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Satellite_Status(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store satellite host halo properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    use Histories
use iso_varying_string
    implicit none
    double precision                        , intent(in   )          :: time
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(inout)          :: doubleBufferCount     , doubleProperty , integerBufferCount, &
         &                                                              integerProperty
    integer         (kind=kind_int8        ), intent(inout)          :: integerBuffer    (:,:)
    double precision                        , intent(inout)          :: doubleBuffer     (:,:)
    class           (nodeComponentBasic    )               , pointer :: thisBasic
    class           (nodeComponentSatellite)               , pointer :: thisSatellite
    type            (history               )                         :: boundMassHistory
    integer         (kind=kind_int8        )                         :: status
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Satellite_Status_Initialize

    ! Store property data if we are outputting satellite orbital pericenter data.
    if (outputSatelliteStatus) then
       ! Test for satellite.
       if (thisNode%isSatellite()) then
          ! Is a satellite. Determine if we have halo information or not.
          thisSatellite    => thisNode     %satellite       ()
          thisBasic        => thisNode     %basic           ()
          boundMassHistory =  thisSatellite%boundMassHistory()
          if (boundMassHistory%exists()) then
             if (thisBasic%time() > boundMassHistory%time(size(boundMassHistory%time))) then
                status=2
             else
                status=1
             end if
          else
             status=2
          end if
       else
          ! Not a satellite.
          status=0
       end if
       ! Store the orbital properties.
       integerProperty=integerProperty+1
       integerBuffer(integerBufferCount,integerProperty)=status
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Status

end module Galacticus_Output_Trees_Satellite_Status
