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

!% Contains a module which outputs data on all merging events which occur.

module Satellites_Merging_Output
  !% Implement outputting of all merging events.
  implicit none
  private
  public :: Satellite_Merging_Output

  ! Module initialization state.
  logical :: initialized           =.false.

  ! Options controlling output.
  logical :: outputSatelliteMergers        , outputSatelliteMergersMainBranchOnly
  
contains
  
  !# <satelliteMergerTask>
  !#  <unitName>Satellite_Merging_Output</unitName>
  !# </satelliteMergerTask>
  subroutine Satellite_Merging_Output(node)
    !% Outputs properties of merging nodes.
    use Galacticus_Nodes
    use Galacticus_HDF5
    use IO_HDF5
    use Input_Parameters
    implicit none
    type (treeNode          ), intent(inout), pointer :: node
    type (treeNode          )               , pointer :: host
    class(nodeComponentBasic)               , pointer :: basic       , basicHost
    type (hdf5Object        )                         :: mergersGroup
    
    ! Initialize the module if necessary.
    !$omp critical (Satellite_Merging_Output_Initialize)
    if (.not.initialized) then
       ! Read controlling parameters.
       !@ <inputParameter>
       !@   <name>outputSatelliteMergers</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>Specifies whether satellite merger information should be output.</description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>???</group>
       !@ </inputParameter>
       call Get_Input_Parameter("outputSatelliteMergers",outputSatelliteMergers,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>outputSatelliteMergersMainBranchOnly</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>Specifies whether satellite merger information should be output only for the main branch host halo.</description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>???</group>
       !@ </inputParameter>
         call Get_Input_Parameter("outputSatelliteMergersMainBranchOnly",outputSatelliteMergersMainBranchOnly,defaultValue=.true.)
     initialized=.true.
    end if
    !$omp end critical (Satellite_Merging_Output_Initialize)
    ! Exit if merger data is not to be output.
    if (.not.outputSatelliteMergers) return
    ! Find the node to merge with.
    host => node%mergesWith()
    ! Exit if the host is not on the main branch and we want to output only main branch mergers.
    if (outputSatelliteMergersMainBranchOnly.and..not.host%isOnMainBranch()) return
    ! Get the basic components of both node and host.
    basic     => node%basic()
    basicHost => host%basic()
    ! Open the group to which satellite mergers should be written.
    !$omp critical (HDF5_Access)
    mergersGroup=galacticusOutputFile%openGroup("satelliteMergers","Satellite mergers data.")    
    ! Append to the datasets.
    call mergersGroup%writeDataset([node              %index()],"indexSatellite","Index of the satellite."               ,appendTo=.true.)
    call mergersGroup%writeDataset([host              %index()],"indexHost"     ,"Index of the satellite."               ,appendTo=.true.)
    call mergersGroup%writeDataset([node     %hostTree%index  ],"indexTree"     ,"Index of the tree."                    ,appendTo=.true.)
    call mergersGroup%writeDataset([basic             %time ()],"time"          ,"Time at which the merger occurs [Gyr].",appendTo=.true.)
    call mergersGroup%writeDataset([basic             %mass ()],"massSatellite" ,"Mass of the satellite halo [M☉]."      ,appendTo=.true.)
    call mergersGroup%writeDataset([basicHost         %mass ()],"massHost"      ,"Mass of the host halo [M☉]."           ,appendTo=.true.)
    ! Close the group.
    call mergersGroup%close()
    !$omp end critical (HDF5_Access)
    return
  end subroutine Satellite_Merging_Output

end module Satellites_Merging_Output
