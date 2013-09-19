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

!% Contains a module that provides and object that implements satellite merging timescales.

module Satellite_Merging_Timescales
  !% Provides and object that implements satellite merging timescales.
  use ISO_Varying_String
  use Galacticus_Nodes
  use Kepler_Orbits
  !# <include directive="satelliteMergingTimescales" type="functionModules" >
  include 'satelliteMergingTimescales.functionModules.inc'
  !# </include>
  private
  public :: Satellite_Merging_Timescale_State_Store   , Satellite_Merging_Timescale_State_Retrieve, &
       &    Satellite_Merging_Timescale_State_Snapshot

  !# <include directive="satelliteMergingTimescales" type="function" >
  !#  <description>Object providing merging timescales for satellites.</description>
  !#  <default>jiang2008</default>
  !#  <defaultThreadPrivate>yes</defaultThreadPrivate>
  !#  <method name="timeUntilMerging" >
  !#   <description>Return the time (in Gyr) until the satellite will merge with its host given the current orbit.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode   ), intent(inout), pointer :: thisNode</argument>
  !#   <argument>type(keplerOrbit), intent(inout)          :: thisOrbit</argument>
  !#  </method>
  !#  <method name="stateStore" >
  !#   <description>Store the state of the object to file.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <modules>FGSL</modules>
  !#   <argument>integer           , intent(in   ) :: stateFile    </argument>
  !#   <argument>type   (fgsl_file), intent(in   ) :: fgslStateFile</argument>
  !#   <code>! Do nothing.</code>
  !#  </method>
  !#  <method name="stateRestore" >
  !#   <description>Restore the state of the object from file.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <modules>FGSL</modules>
  !#   <argument>integer           , intent(in   ) :: stateFile    </argument>
  !#   <argument>type   (fgsl_file), intent(in   ) :: fgslStateFile</argument>
  !#   <code>! Do nothing.</code>
  !#  </method>
  !#  <method name="stateSnapshot" >
  !#   <description>Stores a snapshot of the object state.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <code>! Do nothing.</code>
  !#  </method>
  include 'satelliteMergingTimescales.type.inc'
  !# </include>

  !# <galacticusStateStoreTask>
  !#  <unitName>Satellite_Merging_Timescale_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Satellite_Merging_Timescale_State_Store(stateFile,fgslStateFile)
    !% Store the state to file.
    implicit none
    integer                                 , intent(in   ) :: stateFile
    type   (fgsl_file                      ), intent(in   ) :: fgslStateFile
    class  (satelliteMergingTimescalesClass), pointer       :: thisSatelliteMergingTimescales

    thisSatelliteMergingTimescales => satelliteMergingTimescales()
    call thisSatelliteMergingTimescales%stateStore(stateFile,fgslStateFile)
    return
  end subroutine Satellite_Merging_Timescale_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Satellite_Merging_Timescale_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Satellite_Merging_Timescale_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the state from file.
    implicit none
    integer                                 , intent(in   ) :: stateFile
    type   (fgsl_file                      ), intent(in   ) :: fgslStateFile
    class  (satelliteMergingTimescalesClass), pointer       :: thisSatelliteMergingTimescales

    thisSatelliteMergingTimescales => satelliteMergingTimescales()
    call thisSatelliteMergingTimescales%stateRestore(stateFile,fgslStateFile)
    return
  end subroutine Satellite_Merging_Timescale_State_Retrieve

  !# <galacticusStateSnapshotTask>
  !#  <unitName>Satellite_Merging_Timescale_State_Snapshot</unitName>
  !# </galacticusStateSnapshotTask>
  subroutine Satellite_Merging_Timescale_State_Snapshot()
    !% Retrieve the state from file.
    implicit none
    class(satelliteMergingTimescalesClass), pointer :: thisSatelliteMergingTimescales

    thisSatelliteMergingTimescales => satelliteMergingTimescales()
    call thisSatelliteMergingTimescales%stateSnapshot()
    return
  end subroutine Satellite_Merging_Timescale_State_Snapshot

end module Satellite_Merging_Timescales
