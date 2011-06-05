!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements reading of merger trees from an HDF5 file.

module Merger_Tree_Read
  !% Implements reading of merger trees from an HDF5 file.
  use Merger_Trees
  use Tree_Nodes
  use ISO_Varying_String
  use IO_HDF5
  use HDF5
  use Kind_Numbers
  private
  public :: Merger_Tree_Read_Initialize

  ! The name of the file from which to read merger trees and its internal object.
  type(hdf5Object)                                   :: mergerTreeFile
  type(varying_string)                               :: mergerTreeReadFileName

  ! The halo trees group.
  type(hdf5Object)                                   :: haloTreesGroup

  ! The particles group.
  type(hdf5Object)                                   :: particlesGroup
  integer                                            :: particleEpochType
  integer,          parameter                        :: particleEpochTypeTime           =0
  integer,          parameter                        :: particleEpochTypeExpansionFactor=1
  integer,          parameter                        :: particleEpochTypeRedshift       =2
  character(len=16)                                  :: particleEpochDatasetName

  ! Cosmological parameters group.
  type(hdf5Object)                                   :: cosmologicalParametersGroup

  ! Tree indexing group.
  type(hdf5Object)                                   :: treeIndexGroup

  ! Units group.
  type(hdf5Object)                                   :: unitsGroup

  ! Simulation properties group.
  type(hdf5Object)                                   :: simulationGroup

  ! Index of the next merger tree to read.
  integer                                            :: nextTreeToRead=0

  ! Index of the first tree to process.
  integer(kind=kind_int8)                            :: mergerTreeReadBeginAt

  ! Scaling factors to convert read data to Galacticus internal units.
  double precision                                   :: unitConversionMass         , unitConversionLength,     &
       &                                                unitConversionVelocity     , unitConversionTime
  integer                                            :: scaleFactorExponentMass    , scaleFactorExponentLength,&
       &                                                scaleFactorExponentVelocity, scaleFactorExponentTime

  ! Indexing information for merger tree halos.
  integer,                 allocatable, dimension(:) :: mergerTreeFirstNodeIndex,mergerTreeNodeCount
  integer(kind=kind_int8), allocatable, dimension(:) :: mergerTreeIndex

  ! Volume weight factor for trees.
  double precision                                   :: treeVolumeWeightUniform,treeVolumeWeightCurrent
  double precision,        allocatable, dimension(:) :: treeVolumeWeight
  !$omp threadprivate(treeVolumeWeightCurrent)

  ! Flag indicating whether or not subhalo masses are included in halo masses.
  logical                                            :: haloMassesIncludeSubhalos

  ! Flags indicating whether or not to preset subhalo properties.
  logical                                            :: mergerTreeReadPresetMergerTimes
  logical                                            :: mergerTreeReadPresetMergerNodes
  logical                                            :: mergerTreeReadPresetSubhaloMasses
  logical                                            :: mergerTreeReadPresetPositions

  ! Buffer to hold additional merger trees.
  type(mergerTree),        allocatable, dimension(:) :: mergerTreesQueued
  integer                                            :: mergerTreeQueuePosition

  ! Labels used for node properties.
  integer,                 parameter                 :: nodeIsUnreachable=-1
  integer,                 parameter                 :: nodeIsReachable  = 0

  ! Internal list of output times and the relative tolerance used to "snap" nodes to output times.
  integer                                            :: outputTimesCount
  double precision                                   :: mergerTreeReadOutputTimeSnapTolerance
  double precision,        allocatable, dimension(:) :: outputTimes

  ! Type used to store raw data.
  type nodeData
     !% Structure used to store raw data read from merger tree files.
     integer(kind=kind_int8)                         :: nodeIndex,hostIndex,descendentIndex,particleIndexStart,particleIndexCount&
          &,isolatedNodeIndex,mergesWithIndex,treeIndex
     double precision                                :: nodeMass,nodeTime
     double precision,       dimension(3)            :: position,velocity
     logical                                         :: isSubhalo,childIsSubhalo
     type(nodeData),         pointer                 :: descendentNode,parentNode,hostNode
     type(treeNode),         pointer                 :: node
  end type nodeData

contains

  !# <mergerTreeConstructMethod>
  !#  <unitName>Merger_Tree_Read_Initialize</unitName>
  !# </mergerTreeConstructMethod>
  subroutine Merger_Tree_Read_Initialize(mergerTreeConstructMethod,Merger_Tree_Construct)
    !% Initializes the merger tree reading module.
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Display
    use Galacticus_Output_Times
    use Numerical_Comparison
    use Cosmological_Parameters
    use CDM_Power_Spectrum
    use Numerical_Constants_Astronomical
    use Memory_Management
    implicit none
    type(varying_string),          intent(in)    :: mergerTreeConstructMethod
    procedure(),          pointer, intent(inout) :: Merger_Tree_Construct
    integer                                      :: treesAreSelfContained,hubbleExponent,haloMassesIncludeSubhalosInteger,iOutput
    double precision                             :: cosmologicalParameter,lengthSimulationBox
    character(len=14)                            :: valueString
    type(varying_string)                         :: message

    ! Check if our method is to be used.
    if (mergerTreeConstructMethod == 'read') then
       ! Assign pointer to our merger tree construction subroutine.
       Merger_Tree_Construct => Merger_Tree_Read_Do
       ! Read parameters for halo mass sampling.
       !@ <inputParameter>
       !@   <name>mergerTreeReadFileName</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum mass of merger tree base halos to consider when building merger trees, in units of $M_\odot$.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadFileName',mergerTreeReadFileName)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetMergerTimes</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether merging times for subhalos should be preset when reading merger trees from a file.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetMergerTimes',mergerTreeReadPresetMergerTimes,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetMergerNodes</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether the target nodes for mergers should be preset (i.e. determined from descendent nodes). If they are not, merging will be with each satellite's host node.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetMergerNodes',mergerTreeReadPresetMergerNodes,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetSubhaloMasses</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether subhalo mass should be preset when reading merger trees from a file.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetSubhaloMasses',mergerTreeReadPresetSubhaloMasses,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetPositions</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether node positions should be preset when reading merger trees from a file.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetPositions',mergerTreeReadPresetPositions,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadBeginAt</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies the index of the tree to begin at. (Use -1 to always begin with the first tree.)
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadBeginAt',mergerTreeReadBeginAt,defaultValue=-1)
       !@ <inputParameter>
       !@   <name>mergerTreeReadOutputTimeSnapTolerance</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The relative tolerance required to ``snap'' a node time to the closest output time.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadOutputTimeSnapTolerance',mergerTreeReadOutputTimeSnapTolerance,defaultValue=0.0d0)

       ! Validate input parameters.
       if (mergerTreeReadPresetMergerNodes.and..not.mergerTreeReadPresetMergerTimes) call Galacticus_Error_Report("Merger_Tree_Read_Initialize","presetting of merger target nodes requires that merger times also be preset")

       ! Get array of output times.
       outputTimesCount=Galacticus_Output_Time_Count()
       call Alloc_Array(outputTimes,[outputTimesCount])
       do iOutput=1,outputTimesCount
          outputTimes(iOutput)=Galacticus_Output_Time(iOutput)
       end do

       ! Read basic data from the merger tree file.
       ! Open the file.
       call mergerTreeFile%openFile(char(mergerTreeReadFileName),readOnly=.true.)
       ! Open the merger trees group.
       haloTreesGroup=IO_HDF5_Open_Group(mergerTreeFile,"haloTrees")

       ! Check that we can handle the type of tree in this file.
       if (haloTreesGroup%hasAttribute("haloMassesIncludeSubhalos")) then
          call haloTreesGroup%readAttribute("haloMassesIncludeSubhalos",haloMassesIncludeSubhalosInteger)
          haloMassesIncludeSubhalos=(haloMassesIncludeSubhalosInteger == 1)
       else
          call Galacticus_Error_Report('Merger_Tree_Read_Initialize','required attribute "haloMassesIncludeSubhalos" not present')
       end if

       ! Determine if subhalo masses have been included in halo masses.
       if (haloTreesGroup%hasAttribute("treesAreSelfContained")) then
          call haloTreesGroup%readAttribute("treesAreSelfContained",treesAreSelfContained)
          if (treesAreSelfContained == 0) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','only self-contained trees are supported')
       end if

       ! Compute unit conversion factors.
       unitsGroup=IO_HDF5_Open_Group(mergerTreeFile,"units")
       call unitsGroup%readAttribute("massUnitsInSI"              ,unitConversionMass       )
       call unitsGroup%readAttribute("massScaleFactorExponent"    ,scaleFactorExponentMass  )
       call unitsGroup%readAttribute("massHubbleExponent"         ,hubbleExponent           )
       unitConversionMass    =unitConversionMass    *(Little_H_0()**hubbleExponent)/massSolar
       call unitsGroup%readAttribute("lengthUnitsInSI"            ,unitConversionLength     )
       call unitsGroup%readAttribute("lengthScaleFactorExponent"  ,scaleFactorExponentLength)
       call unitsGroup%readAttribute("lengthHubbleExponent"       ,hubbleExponent           )
       unitConversionLength  =unitConversionLength  *(Little_H_0()**hubbleExponent)/megaParsec
       call unitsGroup%readAttribute("velocityUnitsInSI"          ,unitConversionVelocity     )
       call unitsGroup%readAttribute("velocityScaleFactorExponent",scaleFactorExponentVelocity)
       call unitsGroup%readAttribute("velocityHubbleExponent"     ,hubbleExponent           )
       unitConversionVelocity=unitConversionVelocity*(Little_H_0()**hubbleExponent)/kilo
       if (unitsGroup%hasAttribute("timeUnitsInSI")) then
          call unitsGroup%readAttribute("timeUnitsInSI"           ,unitConversionTime       )
          call unitsGroup%readAttribute("timeScaleFactorExponent" ,scaleFactorExponentTime  )
          call unitsGroup%readAttribute("timeHubbleExponent"      ,hubbleExponent           )
          unitConversionTime =unitConversionTime    *(Little_H_0()**hubbleExponent)/gigayear
          if (scaleFactorExponentTime /= 0) call Galacticus_Error_Report("Merger_Tree_Read_Do","expect no scaling of time units with expansion factor")
       end if
       call unitsGroup%close()

       ! Get the volume of the simulation.
       simulationGroup=IO_HDF5_Open_Group(mergerTreeFile,"simulation")
       if (simulationGroup%hasAttribute("boxSize")) then
          call simulationGroup%readAttribute("boxSize",lengthSimulationBox)
          lengthSimulationBox=lengthSimulationBox*unitConversionLength
          treeVolumeWeightUniform=1.0d0/lengthSimulationBox**3
       else
          call Galacticus_Error_Report('Merger_Tree_Read_Initialize','the boxSize attribute of the simulation group is required')
       end if
       call simulationGroup%close()

       ! Check that cosmological parameters are consistent with the internal ones.
       cosmologicalParametersGroup=IO_HDF5_Open_Group(mergerTreeFile,"cosmology")
       if (cosmologicalParametersGroup%hasAttribute("Omega0")) then
          call cosmologicalParametersGroup%readAttribute("Omega0",cosmologicalParameter)
          if (Values_Differ(cosmologicalParameter,Omega_0(),absTol=0.001d0)) then
             message='Omega_0 in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') Omega_0()
             message=message//trim(valueString)//']'
             call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
          end if
       end if
       if (cosmologicalParametersGroup%hasAttribute("OmegaBaryon")) then
          call cosmologicalParametersGroup%readAttribute("OmegaBaryon",cosmologicalParameter)
          if (Values_Differ(cosmologicalParameter,Omega_b(),absTol=0.001d0)) then
             message='Omega_b in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') Omega_b()
             message=message//trim(valueString)//']'
             call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
          end if
       end if
       if (cosmologicalParametersGroup%hasAttribute("OmegaLambda")) then
          call cosmologicalParametersGroup%readAttribute("OmegaLambda",cosmologicalParameter)
          if (Values_Differ(cosmologicalParameter,Omega_DE(),absTol=0.001d0)) then
             message='Omega_DE in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') Omega_DE()
             message=message//trim(valueString)//']'
             call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
          end if
       end if
       if (cosmologicalParametersGroup%hasAttribute("HubbleParam")) then
          call cosmologicalParametersGroup%readAttribute("HubbleParam",cosmologicalParameter)
          if (Values_Differ(cosmologicalParameter,Little_H_0(),absTol=0.00001d0)) then
             message='H_0 in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') Little_H_0()
             message=message//trim(valueString)//']'
             call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
          end if
       end if
       if (cosmologicalParametersGroup%hasAttribute("sigma_8")) then
          call cosmologicalParametersGroup%readAttribute("sigma_8",cosmologicalParameter)
          if (Values_Differ(cosmologicalParameter,sigma_8(),absTol=0.00001d0)) then
             message='sigma_0 in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') sigma_8()
             message=message//trim(valueString)//'] - may not matter if sigma_8 is not other functions'
             call Galacticus_Display_Message(message)
          end if
       end if

       ! Read the indexing data for merger tree halos.
       treeIndexGroup=IO_HDF5_Open_Group(mergerTreeFile,"treeIndex")
       call treeIndexGroup%readDataset("firstNode"    ,mergerTreeFirstNodeIndex)
       call treeIndexGroup%readDataset("numberOfNodes",mergerTreeNodeCount     )
       call treeIndexGroup%readDataset("treeIndex"    ,mergerTreeIndex         )
       if (treeIndexGroup%hasDataset("treeWeight")) then
          call treeIndexGroup%readDataset("treeWeight",treeVolumeWeight)
          treeVolumeWeight=treeVolumeWeight/unitConversionLength**3
       end if
       call treeIndexGroup%close()

       ! Check that position information is present if required.
       if (mergerTreeReadPresetPositions) then
          if (.not.(haloTreesGroup%hasDataset("position").and.haloTreesGroup%hasDataset("velocity"))) call&
               & Galacticus_Error_Report("Merger_Tree_Read_Initialize","presetting positions requires that both position and&
               & velocity datasets be present in merger tree file")
       end if

       ! Reset first node indices to Fortran array standard.
       mergerTreeFirstNodeIndex=mergerTreeFirstNodeIndex+1

       ! Open the group of particle data.
       if (mergerTreeFile%hasGroup("particles")) then
          particlesGroup=IO_HDF5_Open_Group(mergerTreeFile,"particles")
          if (particlesGroup%hasDataset("time")) then
             particleEpochType=particleEpochTypeTime
             particleEpochDatasetName="time"
          else if (particlesGroup%hasDataset("expansionFactor")) then
             particleEpochType=particleEpochTypeExpansionFactor
             particleEpochDatasetName="expansionFactor"
          else if (particlesGroup%hasDataset("redshift")) then
             particleEpochType=particleEpochTypeRedshift
             particleEpochDatasetName="redshift"
          else
             call Galacticus_Error_Report("Merger_Tree_Read_Do","particles group must have one of time, redshift or expansionFactor datasets")
          end if
       end if
    end if
    return
  end subroutine Merger_Tree_Read_Initialize

  subroutine Merger_Tree_Read_Do(thisTree,skipTree)
    !% Read a merger tree from file.
    use Cosmology_Functions
    use Galacticus_Error
    use String_Handling
    use Memory_Management
    use Kind_Numbers
    use Histories
    implicit none
    type(mergerTree),        intent(inout)                     :: thisTree
    logical,                 intent(in)                        :: skipTree
    double precision,        allocatable, dimension(:)         :: historyTime,historyMass
    double precision,        allocatable, dimension(:,:)       :: position,velocity
    type(nodeData),          allocatable, dimension(:), target :: nodes
    type(treeNodeList),      allocatable, dimension(:)         :: thisNodeList
    logical,                 allocatable, dimension(:)         :: childIsSubhalo
    integer(kind=HSIZE_T),                dimension(1)         :: nodeCount,firstNodeIndex
    integer                                                    :: isolatedNodeCount,primaryRootIndex,iExtraTree
    integer(kind=kind_int8)                                    :: iNode,historyCountMaximum
    logical                                                    :: haveTree

    !$omp critical(mergerTreeReadTree)
    ! Do we have some queued trees to process?
    if (allocated(mergerTreesQueued)) then
       ! We do have queued trees - simply return the next one.
       mergerTreeQueuePosition=mergerTreeQueuePosition+1
       thisTree=mergerTreesQueued(mergerTreeQueuePosition)
       ! If the last tree in the queue has been reached, destroy the queue.
       if (mergerTreeQueuePosition == size(mergerTreesQueued)) then
          call Memory_Usage_Record(sizeof(mergerTreesQueued),addRemove=-1)
          deallocate(mergerTreesQueued)
       end if
    else
       ! We do not have any queued trees. Therefore, increment the tree to read index.
       nextTreeToRead=nextTreeToRead+1
       ! Keep incrementing the tree index until we find the first tree to process (if we haven't done so already). Also skip trees
       ! that contain 1 or fewer nodes and these are unprocessable.
       do while ( (       mergerTreeReadBeginAt            > 0                     &
            &       .and. mergerTreeIndex(nextTreeToRead) /= mergerTreeReadBeginAt &
            &     )                                                                &
            &     .or. mergerTreeNodeCount(nextTreeToRead) <= 1                    &
            &   )
          nextTreeToRead=nextTreeToRead+1
          ! If the end of the list has been reached, exit.
          if (nextTreeToRead > size(mergerTreeFirstNodeIndex)) exit
       end do
       ! Flag that we've now found the first merger tree to process.
       mergerTreeReadBeginAt=-1
       if (nextTreeToRead > size(mergerTreeFirstNodeIndex))  then
          ! All trees have been read.
          ! Close the halo trees group.
          call haloTreesGroup%close()
          ! Close the file.
          call mergerTreeFile%close()
          ! Flag that we do not have a tree.
          haveTree=.false.
       else
          ! Flag that we do have a tree.
          haveTree=.true.
       end if

       ! Continue only if we have a tree.
       if (haveTree) then
          ! If the tree is to be skipped, do not read it.
          if (skipTree) then
             ! Simply allocate a base node to indicate that the tree exists.
             call thisTree%createNode(thisTree%baseNode)
          else
             ! Set tree properties.
             ! treeIndex
             thisTree%index=mergerTreeIndex(nextTreeToRead)
             ! volumeWeight
             if (allocated(treeVolumeWeight)) then
                ! Use the tree-specific weight if available.
                treeVolumeWeightCurrent=treeVolumeWeight(nextTreeToRead)
             else
                ! Otherwise use the same weight for all trees.
                treeVolumeWeightCurrent=treeVolumeWeightUniform
             end if
             thisTree%volumeWeight=treeVolumeWeightCurrent
             
             ! Get the start index and extent of the datasets for this tree.
             firstNodeIndex(1)=mergertreeFirstNodeIndex(nextTreeToRead)
             nodeCount     (1)=mergerTreeNodeCount     (nextTreeToRead)
             
             ! Allocate temporary arrays to hold the merger tree data.
             allocate(nodes(nodeCount(1)))
             call Memory_Usage_Record(sizeof(nodes))
             
             ! Read data from the file.
             call Read_Node_Data(nodes,firstNodeIndex,nodeCount)
             
             ! Convert masses to Galacticus internal units.
             nodes%nodeMass=nodes%nodeMass*unitConversionMass
             if (scaleFactorExponentMass /= 0) then
                do iNode=1,nodeCount(1)
                   nodes(iNode)%nodeMass=nodes(iNode)%nodeMass*Expansion_Factor(nodes(iNode)%nodeTime)**scaleFactorExponentMass
                end do
             end if
             
             ! Read optional datasets.
             call Read_Particle_Data(nodes,firstNodeIndex,nodeCount,position,velocity)
             
             ! Identify subhalos.
             nodes%isSubhalo=nodes%nodeIndex /= nodes%hostIndex
             
             ! Build pointers to descendent nodes.
             call Build_Descendent_Pointers(nodes)
             
             ! Find cases where something that was a subhalo stops being a subhalo.
             call Enforce_Subhalo_Status(nodes)
             
             ! If necessary, add masses of subhalos to host halos.
             if (.not.haloMassesIncludeSubhalos) then
                do iNode=1,nodeCount(1)
                   if (nodes(iNode)%hostNode%nodeIndex /= nodes(iNode)%nodeIndex) nodes(iNode)%hostNode%nodeMass= &
                        &                                                         nodes(iNode)%hostNode%nodeMass+nodes(iNode)%nodeMass
                end do
             end if
             
             ! Associate parent pointers with the descendent host.
             call Build_Parent_Pointers(nodes)
             
             ! Create an array of standard nodes.
             call Create_Node_Array(thisTree,nodes,thisNodeList,isolatedNodeCount,childIsSubhalo)
             
             ! Assign nodes to individual trees.
             call Assign_Nodes_To_Trees(nodes,iExtraTree,primaryRootIndex)
             
             ! Assign parent pointers and properties.
             call Build_Isolated_Parent_Pointers(thisTree,nodes,thisNodeList,primaryRootIndex,iExtraTree)
             
             ! Now build child and sibling links.
             call Build_Child_and_Sibling_Links(nodes,thisNodeList,childIsSubhalo)
             
             ! Check that all required properties exist.
             if (mergerTreeReadPresetPositions) then
                ! Position and velocity methods are required.
                if (.not.(associated(Tree_Node_Position_Set).and.associated(Tree_Node_Velocity_Set))) call&
                     & Galacticus_Error_Report('Merger_Tree_Read_Do','presetting positions requires a component that supports&
                     & position and velocity setting')
             end if
             if (mergerTreeReadPresetMergerTimes) then
                ! Time of merging property is required.
                if (.not.associated(Tree_Node_Satellite_Time_of_Merging_Set)) call Galacticus_Error_Report('Merger_Tree_Read_Do',&
                     & 'presetting merging times requires a component that supports setting of merging times')
             end if
             
             ! Assign isolated node indices to subhalos.
             call Assign_Isolated_Node_Indices(nodes,thisNodeList)
             
             ! Scan subhalos to determine when and how they merge.
             call Scan_For_Mergers(nodes,thisNodeList,historyCountMaximum)
             
             ! Search for any nodes which were flagged as merging with another node and assign appropriate pointers.
             call Assign_Mergers(nodes,thisNodeList)
             
             ! Allocate arrays for history building.
             if (historyCountMaximum > 0) then
                if (allocated(position)) call Dealloc_Array(position)
                if (allocated(velocity)) call Dealloc_Array(velocity)
                call Alloc_Array(historyTime,[int(historyCountMaximum)])
                if (mergerTreeReadPresetSubhaloMasses) call Alloc_Array(historyMass,[  int(historyCountMaximum)])
                if (mergerTreeReadPresetPositions    ) call Alloc_Array(position   ,[3,int(historyCountMaximum)])
                if (mergerTreeReadPresetPositions    ) call Alloc_Array(velocity   ,[3,int(historyCountMaximum)])
             end if
             
             ! Build subhalo mass histories if required.
             call Build_Subhalo_Mass_Histories(nodes,thisNodeList,historyTime,historyMass,position,velocity)
             
             ! Deallocate history building arrays.
             if (allocated(historyTime)) call Dealloc_Array(historyTime)
             if (allocated(historyMass)) call Dealloc_Array(historyMass)
             if (allocated(position   )) call Dealloc_Array(position   )
             if (allocated(velocity   )) call Dealloc_Array(velocity   )
             
             ! Deallocate the temporary arrays.
             call Memory_Usage_Record(sizeof(nodes       ),addRemove=-1)
             deallocate(nodes)
             call Memory_Usage_Record(sizeof(thisNodeList),addRemove=-1)
             deallocate(thisNodeList)
          end if
       end if
    end if
    !$omp end critical(mergerTreeReadTree)
    return
  end subroutine Merger_Tree_Read_Do

  subroutine Read_Node_Data(nodes,firstNodeIndex,nodeCount)
    !% Read node data from an HDF5 file.
    use Galacticus_Error
    use Cosmology_Functions
    use Arrays_Search
    use Numerical_Comparison
    implicit none
    type(nodeData),        intent(inout), dimension(:) :: nodes
    integer(kind=HSIZE_T), intent(in),    dimension(1) :: nodeCount,firstNodeIndex
    integer(kind=kind_int8)                            :: iNode
    integer                                            :: iOutput

    ! nodeIndex
    call haloTreesGroup%readDatasetStatic("nodeIndex"      ,nodes%nodeIndex      ,firstNodeIndex,nodeCount)
    ! hostIndex
    call haloTreesGroup%readDatasetStatic("hostIndex"      ,nodes%hostIndex      ,firstNodeIndex,nodeCount)
    ! parentNode
    call haloTreesGroup%readDatasetStatic("descendentIndex",nodes%descendentIndex,firstNodeIndex,nodeCount)
    ! nodeMass
    call haloTreesGroup%readDatasetStatic("nodeMass"       ,nodes%nodeMass       ,firstNodeIndex,nodeCount)
    ! nodeTime
    if      (haloTreesGroup%hasDataset("time"           )) then
       ! Time is present, so read it.
       call haloTreesGroup%readDatasetStatic("time"           ,nodes%nodeTime,firstNodeIndex,nodeCount)
       ! Convert to Galacticus internal units.
       nodes%nodeTime=nodes%nodeTime*unitConversionTime
    else if (haloTreesGroup%hasDataset("expansionFactor")) then
       ! Expansion factor is present, read it instead.
       call haloTreesGroup%readDatasetStatic("expansionFactor",nodes%nodeTime,firstNodeIndex,nodeCount)
       ! Convert expansion factors to times.
       do iNode=1,nodeCount(1)
          nodes(iNode)%nodeTime=Cosmology_Age(nodes(iNode)%nodeTime)
       end do
    else if (haloTreesGroup%hasDataset("redshift"       )) then
       ! Redshift is present, read it instead.
       call haloTreesGroup%readDatasetStatic("redshift"       ,nodes%nodeTime,firstNodeIndex,nodeCount)
       ! Convert redshifts to times.
       do iNode=1,nodeCount(1)
          nodes(iNode)%nodeTime=Cosmology_Age(Expansion_Factor_from_Redshift(nodes(iNode)%nodeTime))
       end do
    else
       call Galacticus_Error_Report("Merger_Tree_Read_Do","one of time, redshift or expansionFactor data sets must be present in haloTrees group")
    end if

    ! Snap node times to output times if a tolerance has been specified.
    if (mergerTreeReadOutputTimeSnapTolerance > 0.0d0) then
       ! Loop over all nodes.
       do iNode=1,nodeCount(1)
          ! Find closest output time to the node time.
          iOutput=Search_Array_For_Closest(outputTimes,nodes(iNode)%nodeTime)
          ! Test if this time is sufficiently close that we should snap the node time to it.
          if (Values_Agree(nodes(iNode)%nodeTime,outputTimes(iOutput),relTol=mergerTreeReadOutputTimeSnapTolerance)) nodes(iNode)%nodeTime=outputTimes(iOutput)
       end do
    end if

    return
  end subroutine Read_Node_Data
  
  subroutine Read_Particle_Data(nodes,firstNodeIndex,nodeCount,position,velocity)
    !% Read data on particle positions/velocities.
    use Galacticus_Error
    use Cosmology_Functions
    implicit none
    type(nodeData),        intent(inout), dimension(:) :: nodes
    integer(kind=HSIZE_T), intent(in),    dimension(1) :: nodeCount,firstNodeIndex
    double precision,      intent(inout), dimension(:,:), allocatable :: position,velocity
    integer(kind=kind_int8)                            :: iNode
    
    if (mergerTreeReadPresetPositions) then
       ! position.
       call haloTreesGroup%readDataset("position",position,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
       ! velocity.
       call haloTreesGroup%readDataset("velocity",velocity,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
       ! If a set of most bound particle indices are present, read them.
       if (haloTreesGroup%hasDataset("particleIndexStart").and.haloTreesGroup%hasDataset("particleIndexCount")) then
          call haloTreesGroup%readDatasetStatic("particleIndexStart",nodes%particleIndexStart,firstNodeIndex,nodeCount)
          call haloTreesGroup%readDatasetStatic("particleIndexCount",nodes%particleIndexCount,firstNodeIndex,nodeCount)
       else
          nodes%particleIndexStart=-1
          nodes%particleIndexCount=-1
       end if
       ! Convert to Galacticus internal units.
       position=position*unitConversionLength
       if (scaleFactorExponentLength   /= 0) then
          do iNode=1,nodeCount(1)
             position(:,iNode)=position(:,iNode)*Expansion_Factor(nodes(iNode)%nodeTime)**scaleFactorExponentLength
          end do
       end if
       velocity=velocity*unitConversionVelocity
       if (scaleFactorExponentVelocity /= 0) then
          do iNode=1,nodeCount(1)
             velocity(:,iNode)=velocity(:,iNode)*Expansion_Factor(nodes(iNode)%nodeTime)**scaleFactorExponentVelocity
          end do
       end if
       ! Transfer to nodes.
       forall(iNode=1:nodeCount(1))
          nodes(iNode)%position=position(:,iNode)
          nodes(iNode)%velocity=velocity(:,iNode)
       end forall
    end if
    return
  end subroutine Read_Particle_Data

  subroutine Build_Descendent_Pointers(nodes)
    !% Builds pointers from each node to its descendent node.
    use String_Handling
    use Galacticus_Error
    implicit none
    type(nodeData),      intent(inout), dimension(:), target :: nodes
    integer(kind=kind_int8)                                  :: iNode,jNode
    type(varying_string)                                     :: message
    
    do iNode=1,size(nodes)
       nodes(iNode)%descendentNode => null()
       nodes(iNode)%hostNode       => null()
       do jNode=1,size(nodes)
          if (nodes(iNode)%descendentIndex >= 0 .and. nodes(jNode)%nodeIndex == nodes(iNode)%descendentIndex) &
               & nodes(iNode)%descendentNode => nodes(jNode)
          if (                                       nodes(jNode)%nodeIndex == nodes(iNode)%hostIndex      ) &
               & nodes(iNode)%hostNode       => nodes(jNode)
       end do
       if (nodes(iNode)%descendentIndex >= 0 .and. .not.associated(nodes(iNode)%descendentNode)) then
          message='failed to find descendent node: '
          message=message//nodes(iNode)%descendentIndex//' of '//nodes(iNode)%nodeIndex
          call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
       end if
       if (                                       .not.associated(nodes(iNode)%hostNode      )) then
          message='failed to find host node: '
          message=message//nodes(iNode)%hostIndex//' of '//nodes(iNode)%nodeIndex
          call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
       end if
    end do
    return
  end subroutine Build_Descendent_Pointers

  subroutine Enforce_Subhalo_Status(nodes)
    !% Ensure that any node which was once a subhalo remains a subhalo.
    use Galacticus_Error
    use String_Handling
    implicit none
    type(nodeData), intent(inout), dimension(:), target  :: nodes
    type(nodeData),                              pointer :: descendentNode
    integer(kind=kind_int8)                              :: iNode
    logical                                              :: failed
    type(varying_string)                                 :: message

    do iNode=1,size(nodes)
       if (nodes(iNode)%isSubhalo) then
          descendentNode => nodes(iNode)%descendentNode
          do while (associated(descendentNode))
             ! Is this node isolated?
             if (.not.descendentNode%isSubhalo) then
                ! Check if there is any isolated node which descends into this node.
                if (.not.any(nodes%descendentIndex == descendentNode%nodeIndex .and. nodes%nodeIndex == nodes%hostIndex)) then
                   ! Node is isolated, has no isolated node that descends into it. Therefore, our current node is not allowed to be a subhalo.
                   nodes(iNode)%isSubhalo=.false.
                   nodes(iNode)%hostNode => nodes(iNode)
                   nodes(iNode)%hostIndex=nodes(iNode)%nodeIndex                   
                end if
             end if             
             descendentNode => descendentNode%descendentNode
          end do
       end if
    end do
    ! Check that subhalo enforcement was successful.
    failed=.false.
    do iNode=1,size(nodes)
       ! Find nodes which have no isolated node descending into them.
       if (any(nodes%descendentIndex == nodes(iNode)%nodeIndex) .and. .not.any(nodes%descendentIndex == nodes(iNode)%nodeIndex .and. nodes%nodeIndex == nodes%hostIndex)) then
          ! Such nodes must be subhalos. If they are not, report an error
          if (.not.nodes(iNode)%isSubhalo) then
             if (failed) then
                message=message//', '
             else
                message='failed to enforce persistent subhalo status for node ['
             end if
             message=message//nodes(iNode)%nodeIndex
             failed=.true.
          end if
       end if
    end do
    if (failed) then
       message=message//']'
       call Galacticus_Error_Report('Enforce_Subhalo_Status',message)
    end if
    
    return
  end subroutine Enforce_Subhalo_Status

  subroutine Build_Parent_Pointers(nodes)
    !% Build pointers to node parents.
    use Galacticus_Error
    use String_Handling
    implicit none
    type(nodeData), intent(inout), dimension(:), target  :: nodes
    type(nodeData),                              pointer :: parentNode
    integer(kind=kind_int8)                              :: iNode
    type(varying_string)                                 :: message

    do iNode=1,size(nodes)
       if (associated(nodes(iNode)%descendentNode)) then
          if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
             ! Find an isolated parent node, by repeatedly jumping from host to host.
             parentNode => nodes(iNode)%descendentNode%hostNode
             do while (parentNode%isSubhalo)
              if (associated(parentNode,parentNode%hostNode)) then
                   message='node ['
                   message=message//parentNode%nodeIndex//'] flagged as subhalo is self-hosting - exiting to avoid infinite loop'
                   call Galacticus_Error_Report('Build_Parent_Pointers',message)
                end if
                parentNode => parentNode%hostNode
             end do
             nodes(iNode)%parentNode => parentNode
          else
             nodes(iNode)%parentNode => null()
          end if
       else
          nodes(iNode)%parentNode => null()
       end if
    end do
    ! Check for self-parents.
    do iNode=1,size(nodes)
       if (associated(nodes(iNode)%parentNode)) then
          if (nodes(iNode)%nodeIndex == nodes(iNode)%parentNode%nodeIndex) then
             message='node ['
             message=message//nodes(iNode)%nodeIndex//'] is its own parent - exiting to avoid infinite loop'
             call Galacticus_Error_Report('Build_Parent_Pointers',message)
          end if
       end if
    end do
    return
  end subroutine Build_Parent_Pointers

  subroutine Create_Node_Array(thisTree,nodes,nodeList,isolatedNodeCount,childIsSubhalo)
    !% Create an array of standard nodes and associated structures.
    use Memory_Management
    use Galacticus_Error
    implicit none
    type(mergerTree),   intent(inout)                            :: thisTree
    type(nodeData),     intent(inout), dimension(:)              :: nodes
    type(treeNodeList), intent(inout), dimension(:), allocatable :: nodeList
    logical,            intent(inout), dimension(:), allocatable :: childIsSubhalo
    integer,            intent(out)                              :: isolatedNodeCount
    integer(kind=kind_int8)                                      :: iNode
    integer                                                      :: iIsolatedNode

    ! Determine how many nodes are isolated (i.e. not subhalos).
    isolatedNodeCount=count(.not.nodes%isSubhalo)
    
    ! Allocate nodes.
    allocate(nodeList(isolatedNodeCount))
    call Memory_Usage_Record(sizeof(nodeList))
    call Alloc_Array(childIsSubhalo,[isolatedNodeCount])
    
    ! Create the nodes.
    iIsolatedNode          =0
    nodes%isolatedNodeIndex=nodeIsUnreachable
    do iNode=1,size(nodes)
       if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
          iIsolatedNode=iIsolatedNode+1
          ! Store a record of where this node goes in the isolated node list.
          nodes(iNode)%isolatedNodeIndex=iIsolatedNode
          call thisTree%createNode(nodeList(iIsolatedNode)%node)
          call nodeList(iIsolatedNode)%node%indexSet(nodes(iNode)%nodeIndex)
          nodes(iNode)%node => nodeList(iIsolatedNode)%node
       end if
    end do
    return
  end subroutine Create_Node_Array

  subroutine Assign_Nodes_To_Trees(nodes,iExtraTree,primaryRootIndex)
    !% Determine how many tree roots exist in this tree and assign nodes to these trees.
    use Memory_Management
    implicit none
    type(nodeData), intent(inout), dimension(:), target  :: nodes
    integer,        intent(out)                          :: iExtraTree,primaryRootIndex
    type(nodeData),                              pointer :: thisNode
    integer                                              :: iNode,treeRootCount
    double precision                                     :: baseNodeMass
    
    ! Determine how many tree roots we have.
    treeRootCount=0
    baseNodeMass=-1.0d0
    do iNode=1,size(nodes)
       if (.not.nodes(iNode)%isSubhalo .and. .not.associated(nodes(iNode)%descendentNode)) then
          treeRootCount=treeRootCount+1
          if (nodes(iNode)%nodeMass > baseNodeMass) then
             primaryRootIndex=iNode
             baseNodeMass=nodes(iNode)%nodeMass
          end if
       end if
    end do
    
    ! Scan through all nodes and find which of these trees they belong to.
    nodes%treeIndex=-1       
    do iNode=1,size(nodes)
       thisNode => nodes(iNode)
       do while (associated(thisNode))
          nodes(iNode)%treeIndex=thisNode%nodeIndex
          thisNode => thisNode%hostNode%descendentNode
       end do
    end do
    
    ! Allocate a tree queue for any extra roots.
    if (treeRootCount > 1) then
       if (allocated(mergerTreesQueued)) then
          call Memory_Usage_Record(sizeof(mergerTreesQueued),addRemove=-1)
          deallocate(mergerTreesQueued)
       end if
       allocate(mergerTreesQueued(treeRootCount-1))
       call Memory_Usage_Record(sizeof(mergerTreesQueued))
       iExtraTree=0
       mergerTreeQueuePosition=0
    end if
    return
  end subroutine Assign_Nodes_To_Trees

  subroutine Build_Isolated_Parent_Pointers(thisTree,nodes,nodeList,primaryRootIndex,iExtraTree)
    !% Create parent pointer links between isolated nodes and assign times and masses to those nodes.
    implicit none
    type(mergerTree),   intent(inout)               :: thisTree
    type(nodeData),     intent(inout), dimension(:) :: nodes
    type(treeNodeList), intent(inout), dimension(:) :: nodeList
    integer,            intent(inout)               :: iExtraTree,primaryRootIndex
    integer                                         :: iNode,iIsolatedNode
    
    iIsolatedNode=0
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
          iIsolatedNode=iIsolatedNode+1
          if (associated(nodes(iNode)%parentNode)) then
             nodeList(iIsolatedNode)%node%parentNode  => nodes(iNode)%parentNode%node
          else
             nodeList(iIsolatedNode)%node%parentNode  => null()
             if (iNode == primaryRootIndex) then
                thisTree                     %baseNode    => nodeList(iIsolatedNode)%node
                thisTree                     %index       =  nodes(iNode)%nodeIndex
                thisTree                     %volumeWeight=  treeVolumeWeightCurrent
             else
                iExtraTree=iExtraTree+1
                mergerTreesQueued(iExtraTree)%baseNode    => nodeList(iIsolatedNode)%node
                mergerTreesQueued(iExtraTree)%index       =  nodes(iNode)%nodeIndex
                mergerTreesQueued(iExtraTree)%volumeWeight=  treeVolumeWeightCurrent
             end if
          end if
          call Tree_Node_Mass_Set(nodeList(iIsolatedNode)%node,nodes(iNode)%nodeMass)
          call Tree_Node_Time_Set(nodeList(iIsolatedNode)%node,nodes(iNode)%nodeTime)
       end if
    end do
    return
  end subroutine Build_Isolated_Parent_Pointers

  subroutine Build_Child_and_Sibling_Links(nodes,nodeList,childIsSubhalo)
    !% Build child and sibling links between nodes.
    use Memory_Management
    implicit none
    type(nodeData),     intent(inout), dimension(:)              :: nodes
    type(treeNodeList), intent(inout), dimension(:)              :: nodeList
    logical,            intent(inout), dimension(:), allocatable :: childIsSubhalo
    integer                                                      :: iNode,iIsolatedNode
    logical                                                      :: descendsToSubhalo 

    iIsolatedNode=0
    childIsSubhalo=.false.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
          iIsolatedNode=iIsolatedNode+1
          ! Check if the node has a parent.
          if (associated(nodeList(iIsolatedNode)%node%parentNode)) then
             ! Determine if this node definitely descends to a subhalo - in which case it can never be the primary progenitor.
             descendsToSubhalo=nodes(iNode)%descendentIndex /= nodeList(iIsolatedNode)%node%parentNode%index()                   
             ! It does, so set the child pointer of the parent appropriately.
             if (associated(nodeList(iIsolatedNode)%node%parentNode%childNode)) then
                ! A child is already associated. Check if current node does not descend to a subhalo and is more massive.
                if (.not.descendsToSubhalo                                                                                                             &
                     & .and. (                                                                                                                         &
                     &        childIsSubhalo(nodes(iNode)%parentNode%isolatedNodeIndex)                                                                &
                     &         .or.                                                                                                                    &
                     &        Tree_Node_Mass(nodeList(iIsolatedNode)%node) > Tree_Node_Mass(nodeList(iIsolatedNode)%node%parentNode%childNode) &
                     &       )                                                                                                                         &
                     & ) then
                   ! It is, so make this the main progenitor.
                   nodeList(iIsolatedNode)%node%siblingNode                      => nodeList(iIsolatedNode)%node%parentNode%childNode
                   nodeList(iIsolatedNode)%node%parentNode%childNode             => nodeList(iIsolatedNode)%node
                   ! Record that the main child is now not a subhalo.
                   childIsSubhalo(nodes(iNode)%parentNode%isolatedNodeIndex)=.false.
                else
                   ! It is not, so add after the main child.
                   nodeList(iIsolatedNode)%node%siblingNode                      => nodeList(iIsolatedNode)%node%parentNode%childNode%siblingNode
                   nodeList(iIsolatedNode)%node%parentNode%childNode%siblingNode => nodeList(iIsolatedNode)%node
                end if
             else
                ! No child is currently associated. Simply point to the current node.
                nodeList(iIsolatedNode)%node%parentNode%childNode => nodeList(iIsolatedNode)%node
                ! Record whether or not this child is a known subhalo or not.
                childIsSubhalo(nodes(iNode)%parentNode%isolatedNodeIndex)=descendsToSubhalo
             end if
          end if
       end if
    end do
    call Dealloc_Array(childIsSubhalo)
    return
  end subroutine Build_Child_and_Sibling_Links
  
  subroutine Assign_Isolated_Node_Indices(nodes,nodeList)
    !% Assign to each node the number of the corresponding isolated node.
    implicit none
    type(nodeData),     intent(inout), dimension(:) :: nodes
    type(treeNodeList), intent(inout), dimension(:) :: nodeList
    type(nodeData),     pointer                     :: thisNode
    integer                                         :: iNode,iIsolatedNode
    logical                                         :: endOfBranch

    iIsolatedNode=0
    do iNode=1,size(nodes)
       if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
          iIsolatedNode=iIsolatedNode+1
          ! Find the subset with descendents.
          if (associated(nodes(iNode)%descendentNode)) then
             ! Skip halos that are primary progenitors.
             if (.not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                ! Select the subset which have a subhalo as a descendent.
                if (nodes(iNode)%descendentNode%isSubhalo) then
                   ! Trace descendents until merging or final time.
                   thisNode   => nodes(iNode)%descendentNode
                   endOfBranch=  .false.
                   do while (.not.endOfBranch)
                      ! Record that this node was reachable via descendents of an isolated node.
                      if (thisNode%isolatedNodeIndex == nodeIsUnreachable) thisNode%isolatedNodeIndex=nodeIsReachable
                      if (.not.associated(thisNode%descendentNode)) then
                         ! If there is no descendent then the end of the branch has been reached.                               
                         endOfBranch=.true.
                      else
                         ! Step to the next descendent.
                         thisNode => thisNode%descendentNode
                      end if
                   end do
                end if
             end if
          end if
       end if
    end do
    return
  end subroutine Assign_Isolated_Node_Indices

  subroutine Scan_For_Mergers(nodes,nodeList,historyCountMaximum)
    !% Scan for and record mergers between nodes.
    implicit none
    type(nodeData),          intent(inout), dimension(:) :: nodes
    type(treeNodeList),      intent(inout), dimension(:) :: nodeList
    integer(kind=kind_int8), intent(out)                 :: historyCountMaximum
    type(nodeData),          pointer                     :: thisNode
    type(treeNode),          pointer                     :: firstProgenitor
    double precision,        parameter                   :: timeUntilMergingInfinite=1.0d30
    integer                                              :: iNode,jNode,iIsolatedNode,historyCount
    logical                                              :: branchMerges,branchTipReached,endOfBranch,nodeWillMerge
    double precision                                     :: timeSubhaloMerges
    
    iIsolatedNode        = 0
    historyCountMaximum  = 0
    nodes%mergesWithIndex=-1 
    do iNode=1,size(nodes)
       if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
          iIsolatedNode=iIsolatedNode+1
          ! Find the subset with descendents.
          if (associated(nodes(iNode)%descendentNode)) then
             ! Flag indicating if this is a node for which a merging time should be set.
             nodeWillMerge=.false.
             ! Select the subset which have a subhalo as a descendent.
             if (nodes(iNode)%descendentNode%isSubhalo) then
                ! Trace descendents until merging or final time.
                thisNode        => nodes(iNode)%descendentNode
                endOfBranch     =.false.
                branchTipReached=.false.
                branchMerges    =.false.
                historyCount    =0
                do while (.not.endOfBranch)
                   ! Record which isolated node this node belongs to.
                   thisNode%isolatedNodeIndex=iIsolatedNode
                   ! Increment the history count for this branch.
                   historyCount=historyCount+1
                   ! Test the branch.
                   if (.not.associated(thisNode%descendentNode)) then
                      ! No descendent, indicating tip of branch has been reached
                      branchTipReached        =.true.
                      endOfBranch             =.true.
                      historyCount            =historyCount+max(0,thisNode%particleIndexCount)
                   else if (.not.thisNode%descendentNode%isSubhalo) then
                      ! Descendent is not a subhalo, treat as a merging event.  
                      branchMerges            =.true.
                      endOfBranch             =.true.
                      thisNode%mergesWithIndex=thisNode%descendentNode%nodeIndex
                      historyCount            =historyCount+max(0,thisNode%particleIndexCount)
                      thisNode                => thisNode%descendentNode
                   else
                      ! Merges with another subhalo.
                      do jNode=1,size(nodes)
                         if (                     nodes(jNode)%nodeIndex         /= thisNode%nodeIndex        &
                              & .and.             nodes(jNode)%isolatedNodeIndex /= nodeIsUnreachable         &
                              & .and.  associated(nodes(jNode)%descendentNode                         )) then
                            if (associated(nodes(jNode)%descendentNode, thisNode%descendentNode)       &
                                 & .and.   nodes(jNode)%nodeMass      > thisNode%nodeMass       ) then
                               ! Another node mergers into current node's descendent subhalo and is more massive than current
                               ! node. Therefore, class this as a subhalo-subhalo merger.
                               branchMerges            =.true.                                  
                               endOfBranch             =.true.
                               thisNode%mergesWithIndex=nodes(jNode)%descendentNode%nodeIndex
                               historyCount            =historyCount+max(0,thisNode%particleIndexCount)
                               thisNode                => thisNode%descendentNode
                               exit
                            end if
                         end if
                      end do
                      ! Step to the next descendent.
                      if (.not.endOfBranch) thisNode => thisNode%descendentNode
                   end if
                end do
                ! Only set a merging time if this node is not the primary progenitor of its parent.
                if (.not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                   ! Record the largest history.
                   historyCountMaximum=max(historyCountMaximum,historyCount)
                   ! Set an appropriate merging time for this subhalo.
                   if (branchTipReached) timeSubhaloMerges=timeUntilMergingInfinite ! Subhalo never merges, so set merging time to effective infinity.
                   if (branchMerges    ) timeSubhaloMerges=thisNode%nodeTime
                   ! Flag that this node will merge.
                   nodeWillMerge=.true.
                end if
             else if (.not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                ! Descendent is not a subhalo but this node is not the primary progenitor. Assume instantaneous merging.
                timeSubhaloMerges=Tree_Node_Time(nodeList(iIsolatedNode)%node)
                ! Flag that this node will merge.
                nodeWillMerge=.true.
                ! Record the node with which the merger occurs.
                nodes(iNode)%mergesWithIndex=nodes(iNode)%descendentNode%nodeIndex
                ! Ensure the history arrays will be large enough to hold data for this node.
                historyCountMaximum=max(historyCountMaximum,max(0,nodes(iNode)%particleIndexCount))
             end if
             ! Set a merging time if this node will merge.
             if (nodeWillMerge.and.mergerTreeReadPresetMergerTimes) then
                ! Store the time of merging for this node and all of its primary progenitors.
                firstProgenitor => nodeList(iIsolatedNode)%node
                do while (associated(firstProgenitor))
                   call Tree_Node_Satellite_Time_of_Merging_Set(firstProgenitor,timeSubhaloMerges)
                   firstProgenitor => firstProgenitor%childNode
                end do
             end if
          end if
          ! Set position and velocity if required.
          if (mergerTreeReadPresetPositions) then
             call Tree_Node_Position_Set(nodeList(iIsolatedNode)%node,nodes(iNode)%position)
             call Tree_Node_Velocity_Set(nodeList(iIsolatedNode)%node,nodes(iNode)%velocity)
          end if
       end if
    end do
    return
  end subroutine Scan_For_Mergers

  subroutine Assign_Mergers(nodes,nodeList)
    !% Assign pointers to merge targets.
    use Galacticus_Error
    use String_Handling
    implicit none
    type(nodeData),     intent(inout), dimension(:) :: nodes
    type(treeNodeList), intent(inout), dimension(:) :: nodeList
    type(treeNode),     pointer                     :: rootNode
    integer                                         :: iNode,jNode
    type(varying_string)                            :: message
 
    if (mergerTreeReadPresetMergerNodes) then
       do iNode=1,size(nodes)
          ! Check if this node was flagged as merging with another node.
          if (nodes(iNode)%mergesWithIndex /= -1) then
             ! Search for the node that it merges with.
             do jNode=1,size(nodes)
                if (nodes(jNode)%nodeIndex == nodes(iNode)%mergesWithIndex) then
                   if (nodes(jNode)%isolatedNodeIndex <= 0) then
                      ! This node does not belong to any isolated halo - this should not happen.
                      message='subhalo-subhalo ['
                      message=message//nodes(iNode)%nodeIndex//":"//nodes(jNode)%nodeIndex
                      message=message//'] merger in which subhalo has no isolated node progenitor - this should not happen'
                      call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
                   else
                      ! Find the root node into which this node will descend.
                      rootNode => nodeList(nodes(iNode)%isolatedNodeIndex)%node
                      do while (associated(rootNode%parentNode))
                         rootNode => rootNode%parentNode
                      end do
                      ! Set a pointer between the isolated nodes corresponding to these subhalos if and only if they descend
                      ! into the same root node.
                      if (rootNode%index() == nodes(jNode)%treeIndex) then
                         ! Set pointer from merging node (a.k.a. the "mergee") to node that will be merged with.
                         nodeList(nodes(iNode)%isolatedNodeIndex)%node%mergeNode => nodeList(nodes(jNode)%isolatedNodeIndex)%node
                         ! Make a backward pointer from the merge target to the mergee. Check if the target already has mergees associated with it.
                         if (associated(nodeList(nodes(jNode)%isolatedNodeIndex)%node%mergeeNode)) then
                            ! It does: unlink them and attached to the "nextMergee" pointer of the current mergee.
                            nodeList(nodes(iNode)%isolatedNodeIndex)%node%nextMergee => nodeList(nodes(jNode)%isolatedNodeIndex)%node%mergeeNode
                         else
                            ! It does not: simply nullify the next mergee pointer of the mergee.
                            nodeList(nodes(iNode)%isolatedNodeIndex)%node%nextMergee => null()
                         end if
                         ! Append the mergee as the first mergee on the target node.
                         nodeList(nodes(jNode)%isolatedNodeIndex)%node%mergeeNode => nodeList(nodes(iNode)%isolatedNodeIndex)%node
                      end if
                   end if
                   exit
                end if
             end do
          end if
       end do
    end if
    return
  end subroutine Assign_Mergers

  subroutine Build_Subhalo_Mass_Histories(nodes,nodeList,historyTime,historyMass,position,velocity)
    !% Build and attached bound mass histories to subhalos.
    use Galacticus_Error
    use String_Handling
    use Cosmology_Functions
    use Histories
    implicit none
    type(nodeData),     intent(inout), dimension(:),  target  :: nodes
    type(treeNodeList), intent(inout), dimension(:)           :: nodeList
    double precision,   intent(inout), dimension(:)           :: historyTime,historyMass
    double precision,   intent(inout), dimension(:,:)         :: position,velocity
    type(nodeData),                                   pointer :: thisNode
    type(treeNode),                                   pointer :: firstProgenitor
    integer                                                   :: iIsolatedNode,iNode,jNode,iAxis,iTime,historyCount
    logical                                                   :: endOfBranch
    double precision                                          :: expansionFactor
    type(varying_string)                                      :: message
    type(history)                                             :: subhaloHistory
    
    if (mergerTreeReadPresetSubhaloMasses.or.mergerTreeReadPresetPositions) then
       ! Check that preset subhalo masses are supported.
       if (mergerTreeReadPresetSubhaloMasses.and..not.associated(Tree_Node_Bound_Mass_History)) &
            & call Galacticus_Error_Report('Merger_Tree_Read_Do','presetting subhalo masses requires a component that supports node bound mass histories')
       iIsolatedNode=0
       historyBuildNodeLoop: do iNode=1,size(nodes)
          historyBuildIsolatedSelect: if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
             iIsolatedNode=iIsolatedNode+1
             ! Find the subset with descendents.
             historyBuildHasDescendentSelect: if (associated(nodes(iNode)%descendentNode)) then
                ! Set a pointer to the current node - this will be updated if any descendents are traced.
                thisNode => nodes(iNode)
                ! Set initial number of times in the history to zero.
                historyCount=0
                ! Select the subset which have a subhalo as a descendent.
                historyBuildSubhaloSelect: if (nodes(iNode)%descendentNode%isSubhalo) then
                  ! Trace descendents until merging or final time.
                   thisNode    => nodes(iNode)%descendentNode
                   endOfBranch =.false.
                   historyBuildBranchWalk: do while (.not.endOfBranch)
                      ! Increment the history count for this branch.
                      historyCount=historyCount+1
                      ! Store the history.
                      historyTime(historyCount)=thisNode%nodeTime
                      if (mergerTreeReadPresetSubhaloMasses) historyMass(  historyCount)=thisNode%nodeMass
                      if (mergerTreeReadPresetPositions    ) position   (:,historyCount)=thisNode%position
                      if (mergerTreeReadPresetPositions    ) velocity   (:,historyCount)=thisNode%velocity
                      ! Test the branch.
                      if (.not.associated(thisNode%descendentNode).or..not.thisNode%descendentNode%isSubhalo) then
                         ! End of branch reached.
                         endOfBranch=.true.
                      else
                         ! Check if merges with another subhalo.
                         do jNode=1,size(nodes)
                            if (                    nodes(jNode)%nodeIndex         /= thisNode%nodeIndex        &
                                 & .and.            nodes(jNode)%isolatedNodeIndex /= nodeIsUnreachable         &
                                 & .and. associated(nodes(jNode)%descendentNode                         )) then
                               if (associated(nodes(jNode)%descendentNode,thisNode%descendentNode) .and. nodes(jNode)%nodeMass > thisNode%nodeMass) then
                                  ! Subhalo-subhalo merger.
                                  endOfBranch =.true.
                                  exit
                               end if
                            end if
                         end do
                         ! Step to the next descendent.
                         if (.not.endOfBranch) thisNode => thisNode%descendentNode
                      end if
                   end do historyBuildBranchWalk
                   ! Set the mass history for this node.
                   if (mergerTreeReadPresetSubhaloMasses) then
                      call subhaloHistory%destroy()
                      call subhaloHistory%create(1,int(historyCount))
                      subhaloHistory%time(:  )=historyTime(1:historyCount)
                      subhaloHistory%data(:,1)=historyMass(1:historyCount)
                      call nodeList(iIsolatedNode)%node%earliestProgenitor(firstProgenitor)
                      call Tree_Node_Bound_Mass_History_Set(firstProgenitor,subhaloHistory)
                   end if
                end if historyBuildSubhaloSelect

                ! Set the position history for this node.
                if (mergerTreeReadPresetPositions.and..not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                   ! Check if particle data is available for this node.
                   if (thisNode%particleIndexStart >= 0) then
                      ! Check that arrays are large enough to hold particle data. They should be. If they are not, it's a
                      ! bug.
                      if (historyCount+thisNode%particleIndexCount > size(historyTime)) then
                         message='history arrays are too small to hold data for node '
                         message=message//nodeList(iIsolatedNode)%node%index()//': ['//historyCount//'+'//thisNode%particleIndexCount//']='//(historyCount+thisNode%particleIndexCount)//'>'//size(historyTime)
                         call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
                      end if
                      ! Read position data.
                      call particlesGroup%readDatasetStatic("position",position   (:,historyCount+1:historyCount&
                           &+thisNode%particleIndexCount),[1_kind_int8,thisNode%particleIndexStart+1],[3_kind_int8&
                           &,thisNode%particleIndexCount])
                      ! Read velocity data.
                      call particlesGroup%readDatasetStatic("velocity",velocity   (:,historyCount+1:historyCount&
                           &+thisNode%particleIndexCount),[1_kind_int8,thisNode%particleIndexStart+1],[3_kind_int8&
                           &,thisNode%particleIndexCount])
                      ! Read epoch data.
                      call particlesGroup%readDatasetStatic(trim(particleEpochDatasetName),historyTime(historyCount&
                           &+1:historyCount +thisNode%particleIndexCount),[thisNode%particleIndexStart+1]&
                           &,[thisNode%particleIndexCount])
                      do iTime=historyCount+1,historyCount+thisNode%particleIndexCount
                         ! Get the cosmic age and expansion factor.
                         select case (particleEpochType)
                         case (particleEpochTypeTime           )
                            expansionFactor=Expansion_Factor(historyTime(iTime))
                         case (particleEpochTypeExpansionFactor)
                            expansionFactor=historyTime(iTime)
                            historyTime(iTime)=Cosmology_Age(expansionFactor)
                         case (particleEpochTypeRedshift       )
                            expansionFactor=Expansion_Factor_from_Redshift(historyTime(iTime))
                            historyTime(iTime)=Cosmology_Age(expansionFactor)
                         end select
                         ! Convert position units.
                         position(:,iTime)=position(:,iTime)*unitConversionLength
                         if (scaleFactorExponentLength /= 0) then
                            do iAxis=1,3
                               position(iAxis,iTime)=position(iAxis,iTime)*expansionFactor**scaleFactorExponentLength
                            end do
                         end if
                         ! Convert velocity units.
                         velocity(:,iTime)=velocity(:,iTime)*unitConversionVelocity
                         if (scaleFactorExponentLength /= 0) then
                            do iAxis=1,3
                               velocity(iAxis,iTime)=velocity(iAxis,iTime)*expansionFactor**scaleFactorExponentVelocity
                            end do
                         end if
                      end do
                      ! Increment the history count for this node.
                      historyCount=historyCount+thisNode%particleIndexCount
                   end if
                   if (historyCount > 0) then
                      call subhaloHistory%destroy()
                      call subhaloHistory%create(6,int(historyCount))
                      subhaloHistory%time(:    )=          historyTime(    1:historyCount)
                      subhaloHistory%data(:,1:3)=transpose(position   (1:3,1:historyCount))
                      subhaloHistory%data(:,4:6)=transpose(velocity   (1:3,1:historyCount))
                      call Tree_Node_Position_6D_History_Set(nodeList(iIsolatedNode)%node,subhaloHistory)
                   end if
                end if
                
             end if historyBuildHasDescendentSelect
          end if historyBuildIsolatedSelect
       end do historyBuildNodeLoop
    end if
    return
  end subroutine Build_Subhalo_Mass_Histories
  
  subroutine Dump_Tree(nodes,highlightNodes)
    !% Dumps the tree structure to a file in a format suitable for processing with \href{http://www.graphviz.org/}{\sc dot}.
    use File_Utilities
    implicit none
    type(nodeData),          intent(in), dimension(:)           :: nodes
    integer(kind=kind_int8), intent(in), dimension(:), optional :: highlightNodes
    integer                                                     :: iNode,fileUnit
    character(len=20)                                           :: color,style

    ! Open an output file and write the GraphViz opening.
    fileUnit=File_Units_Get()
    open(fileUnit,file='mergerTreeConstructReadTree.gv',status='unknown',form='formatted')
    write (fileUnit,*) 'digraph Tree {'

    ! Loop over all nodes.
    do iNode=1,size(nodes)
       ! Write each node, setting the node shape to a box for subhalos and a circle for halos. Node label consists of the node
       ! index plus the redshift, separated by a colon.
       ! Determine node color.
       if (present(highlightNodes)) then
          if (any(highlightNodes == nodes(iNode)%nodeIndex)) then
             color='green'
             style='filled'
          else
             color='black'
             style='solid'
          end if
       else
          color='black'
          style='solid'
       end if
       if (nodes(iNode)%isSubhalo) then
          write (fileUnit,'(a,i16.16,a,i16.16,a,f4.2,a,a,a,a,a)') '"',nodes(iNode)%nodeIndex,'" [shape=box   , label="',nodes(iNode)%nodeIndex,':',nodes(iNode)%nodeTime,'", color=',trim(color),', style=',trim(style),'];'
          ! If a host node is given, add a link to it as a red line.
          if (associated(nodes(iNode)%hostNode)) write (fileUnit,'(a,i16.16,a,i16.16,a)') '"',nodes(iNode)%nodeIndex,'" -> "',nodes(iNode)%hostNode%nodeIndex,'" [color=red];'
       else
          write (fileUnit,'(a,i16.16,a,i16.16,a,f4.2,a,a,a,a,a)') '"',nodes(iNode)%nodeIndex,'" [shape=circle, label="',nodes(iNode)%nodeIndex,':',nodes(iNode)%nodeTime,'", color=',trim(color),', style=',trim(style),'];'
       endif
       ! Make a link to the descendent node using a black line.
       if (associated(nodes(iNode)%descendentNode)) write (fileUnit,'(a,i16.16,a,i16.16,a)') '"',nodes(iNode)%nodeIndex,'" -> "',nodes(iNode)%descendentNode%nodeIndex,'" ;'
    end do

    ! Close the file.
    write (fileUnit,*) '}'
    close(fileUnit)
    return
  end subroutine Dump_Tree

end module Merger_Tree_Read
