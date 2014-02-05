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

!% Contains a module which handles outputting of satellite host properties to the \glc\ output file.

module Galacticus_Output_Trees_Satellite_Host
  !% Handles outputting of satellite host data to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Satellite_Host, Galacticus_Output_Tree_Satellite_Host_Property_Count,&
       & Galacticus_Output_Tree_Satellite_Host_Names

  ! Number of host properties.
  integer :: satelliteHostPropertyCount

  ! Flag indicating whether or not satellite host information is to be output.
  logical :: outputSatelliteHostData

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputSatelliteHostDataInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Satellite_Host_Initialize
    !% Initializes the module by determining whether or not satellite host data should be output.
    use Input_Parameters
    implicit none

    if (.not.outputSatelliteHostDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Satellite_Host_Initialize)
       if (.not.outputSatelliteHostDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputSatelliteHostData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not satellite host data (node mass) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputSatelliteHostData',outputSatelliteHostData,defaultValue=.false.)
          ! Count number of properties to output.
          satelliteHostPropertyCount=0
          if (outputSatelliteHostData) satelliteHostPropertyCount=satelliteHostPropertyCount+1
          ! Flag that module is now initialized.
          outputSatelliteHostDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Satellite_Host_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Host_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Host_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Host</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Satellite_Host_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
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
    call Galacticus_Output_Tree_Satellite_Host_Initialize

    ! Return property names if we are outputting satellite host data.
    if (outputSatelliteHostData) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>satelliteHostMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Mass of the satellite's host halo [Msun].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='satelliteHostMass'
       doublePropertyComments(doubleProperty)="Mass of the satellite's host halo [Msun]."
       doublePropertyUnitsSI (doubleProperty)=massSolar
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Host_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Host_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Host</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Satellite_Host_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of satellite host properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Satellite_Host_Initialize

    ! Increment property count.
    doublePropertyCount=doublePropertyCount+satelliteHostPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Satellite_Host_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Host</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Host</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Satellite_Host(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store satellite host halo properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    use Kepler_Orbits
    use Satellite_Orbits
    implicit none
    double precision                    , intent(in   )          :: time
    type            (treeNode          ), intent(inout), pointer :: thisNode
    integer                             , intent(inout)          :: doubleBufferCount     , doubleProperty , integerBufferCount, &
         &                                                          integerProperty
    integer         (kind=kind_int8    ), intent(inout)          :: integerBuffer    (:,:)
    double precision                    , intent(inout)          :: doubleBuffer     (:,:)
    type            (treeNode          )               , pointer :: hostNode
    class           (nodeComponentBasic)               , pointer :: hostBasic
    double precision                                             :: hostMass

    ! Initialize the module.
    call Galacticus_Output_Tree_Satellite_Host_Initialize

    ! Store property data if we are outputting satellite orbital pericenter data.
    if (outputSatelliteHostData) then
       ! Test for satellite.
       if (thisNode%isSatellite()) then
          ! Find the host node.
          hostNode  => thisNode %parent
          ! Get the basic component of the host halo.
          hostBasic => hostNode %basic ()
          ! Extract the mass.
          hostMass  =  hostBasic%mass  ()
       else
          hostMass  =  -1.0d0
       end if
       ! Store the orbital properties.
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=hostMass
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Host

end module Galacticus_Output_Trees_Satellite_Host
