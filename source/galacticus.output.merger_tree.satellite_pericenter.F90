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

!% Contains a module which handles outputting of satellite orbital pericenter data to the \glc\ output file.

module Galacticus_Output_Trees_Satellite_Pericenter
  !% Handles outputting of satellite orbital pericenter data to the \glc\ output file.
  implicit none
  private
  public :: Galacticus_Output_Tree_Satellite_Pericenter, Galacticus_Output_Tree_Satellite_Pericenter_Property_Count,&
       & Galacticus_Output_Tree_Satellite_Pericenter_Names

  ! Number of orbital properties.
  integer, parameter :: satellitePericenterPropertyCount        =2

  ! Flag indicating whether or not satellite orbital pericenter information is to be output.
  logical            :: outputSatellitePericenterData

  ! Flag indicating whether or not this module has been initialized.
  logical            :: outputSatellitePericenterDataInitialized=.false.

contains

  subroutine Galacticus_Output_Tree_Satellite_Pericenter_Initialize
    !% Initializes the module by determining whether or not satellite pericenter data should be output.
    use Input_Parameters
    implicit none

    if (.not.outputSatellitePericenterDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Satellite_Pericenter_Initialize)
       if (.not.outputSatellitePericenterDataInitialized) then
          !@ <inputParameter>
          !@   <name>outputSatellitePericenterData</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not satellite orbital pericenter data (radius, velocity) should be included in the output.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('outputSatellitePericenterData',outputSatellitePericenterData,defaultValue=.false.)

          ! Flag that module is now initialized.
          outputSatellitePericenterDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Satellite_Pericenter_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Pericenter_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Pericenter_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Pericenter</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Satellite_Pericenter_Names(thisNode,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set the names of satellite orbital pericenter properties to be written to the \glc\ output file.
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
    call Galacticus_Output_Tree_Satellite_Pericenter_Initialize

    ! Return property names if we are outputting satellite orbital pericenter data.
    if (outputSatellitePericenterData) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>satellitePericenterRadius</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Pericenteric radius of satellite orbit [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='satellitePericenterRadius'
       doublePropertyComments(doubleProperty)='Pericenteric radius of satellite orbit [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>satellitePericenterVelocity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Pericenteric velocity of satellite orbit [km/s].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='satellitePericenterVelocity'
       doublePropertyComments(doubleProperty)='Pericenteric velocity of satellite orbit [km/s].'
       doublePropertyUnitsSI (doubleProperty)=kilo
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Pericenter_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Pericenter_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Pericenter</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Satellite_Pericenter_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of satellite orbital pericenter properties to be written to the \glc\ output file.
    use Galacticus_Nodes
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: time
    integer                   , intent(inout)          :: doublePropertyCount, integerPropertyCount

    ! Initialize the module.
    call Galacticus_Output_Tree_Satellite_Pericenter_Initialize

    ! Increment property count if we are outputting satellite orbital pericenter data.
    if (outputSatellitePericenterData) doublePropertyCount=doublePropertyCount+satellitePericenterPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Satellite_Pericenter_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Satellite_Pericenter</unitName>
  !#  <sortName>Galacticus_Output_Tree_Satellite_Pericenter</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Satellite_Pericenter(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store satellite orbital pericenter properties in the \glc\ output file buffers.
    use Galacticus_Nodes
    use Kind_Numbers
    use Kepler_Orbits
    use Satellite_Orbits
    implicit none
    double precision                        , intent(in   )          :: time
    type            (treeNode              ), intent(inout), pointer :: thisNode
    integer                                 , intent(inout)          :: doubleBufferCount          , doubleProperty , integerBufferCount, &
         &                                                              integerProperty
    integer         (kind=kind_int8        ), intent(inout)          :: integerBuffer         (:,:)
    double precision                        , intent(inout)          :: doubleBuffer          (:,:)
    type            (treeNode              )               , pointer :: hostNode
    class           (nodeComponentSatellite)               , pointer :: thisSatelliteComponent
    type            (keplerOrbit           )                         :: thisOrbit
    double precision                                                 :: orbitalRadius              , orbitalVelocity

    ! Initialize the module.
    call Galacticus_Output_Tree_Satellite_Pericenter_Initialize

    ! Store property data if we are outputting satellite orbital pericenter data.
    if (outputSatellitePericenterData) then
       ! Test for satellite.
       if (thisNode%isSatellite()) then
          ! Find the host node.
          hostNode      => thisNode%parent
          ! Get the satellite component.
          thisSatelliteComponent => thisNode%satellite()
          ! Get the orbit for this node.
          thisOrbit=thisSatelliteComponent%virialOrbit()
          ! Get the orbital radius and velocity at pericenter.
          call Satellite_Orbit_Pericenter_Phase_Space_Coordinates(hostNode,thisOrbit,orbitalRadius,orbitalVelocity)
       else
          orbitalRadius  =0.0d0
          orbitalVelocity=0.0d0
       end if
       ! Store the orbital properties.
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=orbitalRadius
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=orbitalVelocity
    end if
    return
  end subroutine Galacticus_Output_Tree_Satellite_Pericenter

end module Galacticus_Output_Trees_Satellite_Pericenter
