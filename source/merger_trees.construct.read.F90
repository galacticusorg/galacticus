!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
  use FGSL
  implicit none
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

  ! Size of the simulation box.
  double precision                                   :: lengthSimulationBox

  ! Flag indicating whether to skip the main tree (and consider only queued trees).
  logical                                            :: mergerTreeReadSkipMainTree

  ! Flag indicating whether mismatches in cosmological parameters should be considered fatal.
  logical                                            :: mergerTreeReadMismatchIsFatal

  ! Flag indicating whether or not subhalo masses are included in halo masses.
  logical                                            :: haloMassesIncludeSubhalos

  ! Flag indicating whether or not positions are periodic.
  logical                                            :: simulationIsPeriodic

  ! Flag indicating whether or not velocities include Hubble flow.
  logical                                            :: velocitiesIncludeHubbleFlow

  ! Flags indicating whether or not to preset subhalo properties.
  logical                                            :: mergerTreeReadPresetMergerTimes
  logical                                            :: mergerTreeReadPresetMergerNodes
  logical                                            :: mergerTreeReadPresetSubhaloMasses
  logical                                            :: mergerTreeReadPresetPositions
  logical                                            :: mergerTreeReadPresetScaleRadii
  double precision                                   :: mergerTreeReadPresetScaleRadiiMinimumMass
  logical                                            :: mergerTreeReadPresetSpins
  logical                                            :: mergerTreeReadPresetOrbits,mergerTreeReadPresetOrbitsBoundOnly&
       &,mergerTreeReadPresetOrbitsSetAll, mergerTreeReadPresetOrbitsAssertAllSet

  ! Option controlling fatality of missing host node condition.
  logical                                            :: mergerTreeReadMissingHostsAreFatal

  ! Option controlling whether tree indices should always be set to the corresponding node index.
  logical                                            :: mergerTreeReadTreeIndexToRootNodeIndex

  ! Buffer to hold additional merger trees.
  type(mergerTree),        allocatable, dimension(:) :: mergerTreesQueued

  ! Labels used for node properties.
  integer,                 parameter                 :: nodeIsUnreachable=-1
  integer,                 parameter                 :: nodeIsReachable  = 0

  ! Labels for switches.
  integer,                 parameter                 :: integerTrue   = 1
  integer,                 parameter                 :: integerFalse  = 0
  integer,                 parameter                 :: integerUnknown=-1

  ! Internal list of output times and the relative tolerance used to "snap" nodes to output times.
  integer                                            :: outputTimesCount
  double precision                                   :: mergerTreeReadOutputTimeSnapTolerance
  double precision,        allocatable, dimension(:) :: outputTimes

  ! Node used in root finding.
  logical                                            :: haveScaleRadii
  type(treeNode),          pointer                   :: activeNode
  double precision                                   :: halfMassRadius
  !$omp threadprivate(activeNode,halfMassRadius)

  ! Sorted node index list.
  integer(kind=kind_int8), allocatable, dimension(:) :: nodeLocations    ,descendentLocations
  integer(kind=kind_int8), allocatable, dimension(:) :: nodeIndicesSorted,descendentIndicesSorted

  ! Type used to store raw data.
  type nodeData
     !% Structure used to store raw data read from merger tree files.
     integer(kind=kind_int8)                         :: nodeIndex,hostIndex,descendentIndex,particleIndexStart,particleIndexCount&
          &,isolatedNodeIndex,mergesWithIndex
     double precision                                :: nodeMass,nodeTime,halfMassRadius,scaleRadius,angularMomentum
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
    integer                                      :: treesAreSelfContained,hubbleExponent,haloMassesIncludeSubhalosInteger&
         &,simulationIsPeriodicInteger,velocitiesIncludeHubbleFlowInteger,iOutput,treesHaveSubhalos
    double precision                             :: cosmologicalParameter
    character(len=14)                            :: valueString
    type(varying_string)                         :: message
    double precision                             :: localLittleH0,localOmegaMatter,localOmegaDE,localOmegaBaryon,localSigma8

    ! Check if our method is to be used.
    if (mergerTreeConstructMethod == 'read') then
       ! Assign pointer to our merger tree construction subroutine.
       Merger_Tree_Construct => Merger_Tree_Read_Do
       ! Read parameters for halo mass sampling.
       !@ <inputParameter>
       !@   <name>mergerTreeReadFileName</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the file from which merger tree data should be read when using the {\tt [mergerTreeConstructMethod]}$=${\tt read} tree construction method.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadFileName',mergerTreeReadFileName)
       !@ <inputParameter>
       !@   <name>mergerTreeReadSkipMainTree</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>false</defaultValue>
       !@   <description>
       !@     Specifies whether the main tree is to be skipped (and only queued trees are processed). Useful primarily for debugging purposes.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadSkipMainTree',mergerTreeReadSkipMainTree,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadMismatchIsFatal</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether mismatches in cosmological parameter values between \glc\ and the merger tree file should be considered fatal.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadMismatchIsFatal',mergerTreeReadMismatchIsFatal,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetMergerTimes</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether merging times for subhalos should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetMergerTimes',mergerTreeReadPresetMergerTimes,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetMergerNodes</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether the target nodes for mergers should be preset (i.e. determined from descendent nodes). If they are not, merging will be with each satellite's host node.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetMergerNodes',mergerTreeReadPresetMergerNodes,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetSubhaloMasses</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether subhalo mass should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetSubhaloMasses',mergerTreeReadPresetSubhaloMasses,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetPositions</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether node positions should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetPositions',mergerTreeReadPresetPositions,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetScaleRadii</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether node scale radii should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetScaleRadii',mergerTreeReadPresetScaleRadii,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetScaleRadiiMinimumMass</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>0</defaultValue>
       !@   <description>
       !@     The minimum halo mass for which scale radii should be preset (if {\tt [mergerTreeReadPresetScaleRadii]}$=${\tt true}).
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetScaleRadiiMinimumMass',mergerTreeReadPresetScaleRadiiMinimumMass,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetSpins</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether node spins should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetSpins',mergerTreeReadPresetSpins,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetOrbits</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether node orbits should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetOrbits',mergerTreeReadPresetOrbits,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetOrbitsSetAll</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Forces all orbits to be set. If the computed orbit does not cross the virial radius, then select one at random instead.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetOrbitsSetAll',mergerTreeReadPresetOrbitsSetAll,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetOrbitsAssertAllSet</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Asserts that all virial orbits must be preset. If any can not be set, \glc\ will stop.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetOrbitsAssertAllSet',mergerTreeReadPresetOrbitsAssertAllSet,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetOrbitsBoundOnly</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether only bound node orbits should be set.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetOrbitsBoundOnly',mergerTreeReadPresetOrbitsBoundOnly,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadBeginAt</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>-1</defaultValue>
       !@   <description>
       !@     Specifies the index of the tree to begin at. (Use -1 to always begin with the first tree.)
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadBeginAt',mergerTreeReadBeginAt,defaultValue=-1_kind_int8)
       !@ <inputParameter>
       !@   <name>mergerTreeReadOutputTimeSnapTolerance</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>0</defaultValue>
       !@   <description>
       !@     The relative tolerance required to ``snap'' a node time to the closest output time.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadOutputTimeSnapTolerance',mergerTreeReadOutputTimeSnapTolerance,defaultValue=0.0d0)
       !@ <inputParameter>
       !@   <name>mergerTreeReadMissingHostsAreFatal</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether nodes with missing host nodes should be considered to be fatal---see \S\ref{sec:MergerTreeFileProcessing}.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadMissingHostsAreFatal',mergerTreeReadMissingHostsAreFatal,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadTreeIndexToRootNodeIndex</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>false</defaultValue>
       !@   <description>
       !@     Specifies whether tree indices should always be set to the index of their root node.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadTreeIndexToRootNodeIndex',mergerTreeReadTreeIndexToRootNodeIndex,defaultValue=.false.)

       ! Get array of output times.
       outputTimesCount=Galacticus_Output_Time_Count()
       call Alloc_Array(outputTimes,[outputTimesCount])
       do iOutput=1,outputTimesCount
          outputTimes(iOutput)=Galacticus_Output_Time(iOutput)
       end do

       ! Get cosmological parameters. We do this in advance to avoid HDF5 thread conflicts.
       localLittleH0   =Little_H_0()
       localOmegaMatter=Omega_Matter()
       localOmegaDE    =Omega_DE()
       localOmegaBaryon=Omega_b()
       localSigma8     =sigma_8()

       ! Read basic data from the merger tree file.
       !$omp critical(HDF5_Access)
       ! Open the file.
       call mergerTreeFile%openFile(char(mergerTreeReadFileName),readOnly=.true.)
       ! Open the merger trees group.
       haloTreesGroup=mergerTreeFile%openGroup("haloTrees")

       ! Determine if trees have subhalos.
       if (haloTreesGroup%hasAttribute("treesHaveSubhalos")) then
          call haloTreesGroup%readAttribute("treesHaveSubhalos",treesHaveSubhalos,allowPseudoScalar=.true.)
       else
          treesHaveSubhalos=integerUnknown
       end if

       ! Validate input parameters.
       if (mergerTreeReadPresetMergerNodes.and..not.mergerTreeReadPresetMergerTimes) then
          message="presetting of merger target nodes requires that merger times also be preset;"//char(10)
          message=message//" try setting [mergerTreeReadPresetMergerTimes]=true."//char(10)
          if (treesHaveSubhalos /= integerTrue) then
             message=message//" Note: presetting merger target nodes and merger times is usually only a good idea when subhalo information is present in the merger trees"
             if (treesHaveSubhalos == integerFalse) then
                message=message//" (which the current trees do not)"
             else
                message=message//" (subhalo presence in the current trees was not specified)"
             end if
          end if
          call Galacticus_Error_Report("Merger_Tree_Read_Initialize",message)
       end if

       ! Perform sanity checks if subhalos are not included.
       if (treesHaveSubhalos == integerFalse) then
          if (mergerTreeReadPresetMergerTimes  ) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','cannot preset merger times as no subhalos are present; try setting [mergerTreeReadPresetMergerTimes]=false')
          if (mergerTreeReadPresetMergerNodes  ) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','cannot preset merger nodes as no subhalos are present; try setting [mergerTreeReadPresetMergerNodes]=false')
          if (mergerTreeReadPresetSubhaloMasses) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','cannot preset subhalo masses as no subhalos are present; try setting [mergerTreeReadPresetSubhaloMasses]=false')
       end if

       ! Check that we can handle the type of tree in this file.
       if (haloTreesGroup%hasAttribute("haloMassesIncludeSubhalos")) then
          call haloTreesGroup%readAttribute("haloMassesIncludeSubhalos",haloMassesIncludeSubhalosInteger,allowPseudoScalar=.true.)
          haloMassesIncludeSubhalos=(haloMassesIncludeSubhalosInteger == 1)
       else
          call Galacticus_Error_Report('Merger_Tree_Read_Initialize','required attribute "haloMassesIncludeSubhalos" not present')
       end if

       ! Determine if subhalo masses have been included in halo masses.
       if (haloTreesGroup%hasAttribute("treesAreSelfContained")) then
          call haloTreesGroup%readAttribute("treesAreSelfContained",treesAreSelfContained,allowPseudoScalar=.true.)
          if (treesAreSelfContained == 0) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','only self-contained trees are supported')
       end if

       ! Compute unit conversion factors.
       unitsGroup=mergerTreeFile%openGroup("units")
       call unitsGroup%readAttribute("massUnitsInSI"              ,unitConversionMass       ,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("massScaleFactorExponent"    ,scaleFactorExponentMass  ,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("massHubbleExponent"         ,hubbleExponent           ,allowPseudoScalar=.true.)
       unitConversionMass    =unitConversionMass    *(localLittleH0**hubbleExponent)/massSolar
       if (unitConversionMass <= 0.0d0) call Galacticus_Error_Report('Merger_Tree_Read_Do','non-positive units for mass')
       if ((mergerTreeReadPresetPositions.or.mergerTreeReadPresetOrbits).and..not.(unitsGroup%hasAttribute("lengthUnitsInSI").and.unitsGroup%hasAttribute("velocityUnitsInSI"))) call Galacticus_Error_Report('Merger_Tree_Read_Do','length and velocity units must be given if positions/velocities or orbits are to be preset')
       if (unitsGroup%hasAttribute("lengthUnitsInSI")) then
          call unitsGroup%readAttribute("lengthUnitsInSI"            ,unitConversionLength     ,allowPseudoScalar=.true.)
          call unitsGroup%readAttribute("lengthScaleFactorExponent"  ,scaleFactorExponentLength,allowPseudoScalar=.true.)
          call unitsGroup%readAttribute("lengthHubbleExponent"       ,hubbleExponent           ,allowPseudoScalar=.true.)
          unitConversionLength  =unitConversionLength  *(localLittleH0**hubbleExponent)/megaParsec
          if (unitConversionLength <= 0.0d0) call Galacticus_Error_Report('Merger_Tree_Read_Do','non-positive units for length')
       end if
       if (unitsGroup%hasAttribute("velocityUnitsInSI")) then
          call unitsGroup%readAttribute("velocityUnitsInSI"          ,unitConversionVelocity     ,allowPseudoScalar=.true.)
          call unitsGroup%readAttribute("velocityScaleFactorExponent",scaleFactorExponentVelocity,allowPseudoScalar=.true.)
          call unitsGroup%readAttribute("velocityHubbleExponent"     ,hubbleExponent             ,allowPseudoScalar=.true.)
          unitConversionVelocity=unitConversionVelocity*(localLittleH0**hubbleExponent)/kilo
          if (unitConversionVelocity <= 0.0d0) call Galacticus_Error_Report('Merger_Tree_Read_Do','non-positive units for velocity')
       end if
       if (unitsGroup%hasAttribute("timeUnitsInSI")) then
          call unitsGroup%readAttribute("timeUnitsInSI"           ,unitConversionTime       ,allowPseudoScalar=.true.)
          call unitsGroup%readAttribute("timeScaleFactorExponent" ,scaleFactorExponentTime  ,allowPseudoScalar=.true.)
          call unitsGroup%readAttribute("timeHubbleExponent"      ,hubbleExponent           ,allowPseudoScalar=.true.)
          unitConversionTime =unitConversionTime    *(localLittleH0**hubbleExponent)/gigayear
          if (scaleFactorExponentTime /= 0) call Galacticus_Error_Report("Merger_Tree_Read_Do","expect no scaling of time units with expansion factor")
       end if
       call unitsGroup%close()

       ! Determine if velocities include Hubble flow.
       if (haloTreesGroup%hasAttribute("velocitiesIncludeHubbleFlow")) then
          call haloTreesGroup%readAttribute("velocitiesIncludeHubbleFlow",velocitiesIncludeHubbleFlowInteger,allowPseudoScalar=.true.)
          velocitiesIncludeHubbleFlow=(velocitiesIncludeHubbleFlowInteger == 1)
       else
          velocitiesIncludeHubbleFlow=.true.
       end if

       ! Determine if positions are periodic.
       if (haloTreesGroup%hasAttribute("positionsArePeriodic")) then
          call haloTreesGroup%readAttribute("positionsArePeriodic",simulationIsPeriodicInteger,allowPseudoScalar=.true.)
          simulationIsPeriodic=(simulationIsPeriodicInteger == 1)
       else
          simulationIsPeriodic=.false.
       end if

       ! Get the volume of the simulation.
       if (mergerTreeFile%hasGroup("simulation")) then
          simulationGroup=mergerTreeFile%openGroup("simulation")
          if (simulationGroup%hasAttribute("boxSize")) then
             call simulationGroup%readAttribute("boxSize",lengthSimulationBox,allowPseudoScalar=.true.)
             lengthSimulationBox=lengthSimulationBox*unitConversionLength
             treeVolumeWeightUniform=1.0d0/lengthSimulationBox**3
          else
             call Galacticus_Error_Report('Merger_Tree_Read_Initialize','the boxSize attribute of the simulation group is required')
          end if
          call simulationGroup%close()
       else
          treeVolumeWeightUniform=1.0d0
          if (simulationIsPeriodic) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','the boxSize attribute of the simulation group is required')
       end if

       ! Check that cosmological parameters are consistent with the internal ones.
       cosmologicalParametersGroup=mergerTreeFile%openGroup("cosmology")
       if (cosmologicalParametersGroup%hasAttribute("OmegaMatter")) then
          call cosmologicalParametersGroup%readAttribute("OmegaMatter",cosmologicalParameter,allowPseudoScalar=.true.)
          if (Values_Differ(cosmologicalParameter,localOmegaMatter,absTol=0.001d0)) then
             message='Omega_Matter in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') localOmegaMatter
             message=message//trim(valueString)//']'
             if (mergerTreeReadMismatchIsFatal) then
                call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
             else
                call Galacticus_Display_Message(message,verbosityWarn)
             end if
          end if
       else if (cosmologicalParametersGroup%hasAttribute("Omega0")) then
          ! <expiry>
          !  <label>Warning regarding use of Omega0 in merger tree files.</label>
          !  <date>25-Aug-2012</date>
          ! </expiry>
          call Galacticus_Display_Message('WARNING: Use of "Omega0" in merger tree files is deprecated - use OmegaMatter instead',verbosityWarn)
          call cosmologicalParametersGroup%readAttribute("Omega0",cosmologicalParameter,allowPseudoScalar=.true.)
          if (Values_Differ(cosmologicalParameter,localOmegaMatter,absTol=0.001d0)) then
             message='Omega_Matter in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') localOmegaMatter
             message=message//trim(valueString)//']'
             if (mergerTreeReadMismatchIsFatal) then
                call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
             else
                call Galacticus_Display_Message(message,verbosityWarn)
             end if
          end if
       end if
       if (cosmologicalParametersGroup%hasAttribute("OmegaBaryon")) then
          call cosmologicalParametersGroup%readAttribute("OmegaBaryon",cosmologicalParameter,allowPseudoScalar=.true.)
          if (Values_Differ(cosmologicalParameter,localOmegaBaryon,absTol=0.001d0)) then
             message='Omega_b in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') localOmegaBaryon
             message=message//trim(valueString)//']'
             if (mergerTreeReadMismatchIsFatal) then
                call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
             else
                call Galacticus_Display_Message(message,verbosityWarn)
             end if
          end if
       end if
       if (cosmologicalParametersGroup%hasAttribute("OmegaLambda")) then
          call cosmologicalParametersGroup%readAttribute("OmegaLambda",cosmologicalParameter,allowPseudoScalar=.true.)
          if (Values_Differ(cosmologicalParameter,localOmegaDE,absTol=0.001d0)) then
             message='Omega_DE in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') localOmegaDE
             message=message//trim(valueString)//']'
             if (mergerTreeReadMismatchIsFatal) then
                call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
             else
                call Galacticus_Display_Message(message,verbosityWarn)
             end if
          end if
       end if
       if (cosmologicalParametersGroup%hasAttribute("HubbleParam")) then
          call cosmologicalParametersGroup%readAttribute("HubbleParam",cosmologicalParameter,allowPseudoScalar=.true.)
          if (Values_Differ(cosmologicalParameter,localLittleH0,absTol=0.00001d0)) then
             message='H_0 in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') localLittleH0
             message=message//trim(valueString)//']'
             if (mergerTreeReadMismatchIsFatal) then
                call Galacticus_Error_Report('Merger_Tree_Read_Initialize',message)
             else
                call Galacticus_Display_Message(message,verbosityWarn)
             end if
          end if
       end if
       if (cosmologicalParametersGroup%hasAttribute("sigma_8")) then
          call cosmologicalParametersGroup%readAttribute("sigma_8",cosmologicalParameter,allowPseudoScalar=.true.)
          if (Values_Differ(cosmologicalParameter,localSigma8,absTol=0.00001d0)) then
             message='sigma_8 in merger tree file ['
             write (valueString,'(e14.8)') cosmologicalParameter
             message=message//trim(valueString)//'] differs from the internal value ['
             write (valueString,'(e14.8)') localSigma8
             message=message//trim(valueString)//'] - may not matter if sigma_8 is not use in other functions'
             call Galacticus_Display_Message(message)
          end if
       end if
       call cosmologicalParametersGroup%close()

       if (mergerTreeFile%hasGroup("treeIndex")) then
          treeIndexGroup=mergerTreeFile%openGroup("treeIndex")
          call treeIndexGroup%readDataset("firstNode"    ,mergerTreeFirstNodeIndex)
          call treeIndexGroup%readDataset("numberOfNodes",mergerTreeNodeCount     )
          call treeIndexGroup%readDataset("treeIndex"    ,mergerTreeIndex         )
          if (treeIndexGroup%hasDataset("treeWeight")) then
             call treeIndexGroup%readDataset("treeWeight",treeVolumeWeight)
             treeVolumeWeight=treeVolumeWeight/unitConversionLength**3
          end if
          call treeIndexGroup%close()
       else
          call Galacticus_Error_Report('Merger_Tree_Read_Initialize','merger tree file must contain the treeIndex group')
       end if

       ! Check that position information is present if required.
       if (mergerTreeReadPresetPositions.or.mergerTreeReadPresetOrbits) then
          if (.not.(haloTreesGroup%hasDataset("position").and.haloTreesGroup%hasDataset("velocity"))) call&
               & Galacticus_Error_Report("Merger_Tree_Read_Initialize","presetting positions requires that both position and&
               & velocity datasets be present in merger tree file")
       end if

       ! Check that half-mass radius information is present if required.
       if (mergerTreeReadPresetScaleRadii) then
          if     (                                                                                                                                                   &
               &  .not.(                                                                                                                                             &
               &        haloTreesGroup%hasDataset("halfMassRadius")                                                                                                  &
               &         .or.                                                                                                                                        &
               &        haloTreesGroup%hasDataset("scaleRadius"   )                                                                                                  &
               &       )                                                                                                                                             &
               & ) call Galacticus_Error_Report(                                                                                                                     &
               &                                "Merger_Tree_Read_Initialize",                                                                                       &  
               &                                "presetting scale radii requires that at least one of halfMassRadius or scaleRadius datasets be present in merger"// &
               &                                "tree file; try setting"//char(10)//"  [mergerTreeReadPresetScaleRadii]=false"                                       &
               &                               )
       end if

       ! Check that angular momentum information is present if required.
       if (mergerTreeReadPresetSpins) then
          if (.not.haloTreesGroup%hasDataset("angularMomentum")) call&
               & Galacticus_Error_Report("Merger_Tree_Read_Initialize","presetting spins requires that the angularMomentum dataset be present in merger tree file; try setting"//char(10)//"  [mergerTreeReadPresetSpins]=false")
       end if

       ! Reset first node indices to Fortran array standard.
       mergerTreeFirstNodeIndex=mergerTreeFirstNodeIndex+1

       ! Open the group of particle data.
       if (mergerTreeFile%hasGroup("particles")) then
          particlesGroup=mergerTreeFile%openGroup("particles")
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
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Merger_Tree_Read_Initialize

  subroutine Merger_Tree_Read_Do(thisTree,skipTree)
    !% Read a merger tree from file.
    use Galacticus_State
    use Cosmology_Functions
    use Galacticus_Error
    use String_Handling
    use Memory_Management
    use Kind_Numbers
    use Histories
    use ISO_Varying_String
    use Merger_Tree_Read_State
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
    logical                                                    :: haveTree,finished
    type(varying_string)                                       :: message

    !$omp critical(mergerTreeReadTree)
    finished=.false.
    do while (.not.finished)
       ! Do we have some queued trees to process?
       if (allocated(mergerTreesQueued)) then
          ! Retrieve stored internal state if possible.
          call Galacticus_State_Retrieve
          ! Take a snapshot of the internal state and store it.
          call Galacticus_State_Snapshot
          message='Storing state for subtree #'
          message=message//(mergerTreeQueuePosition+1)//" of tree #"//nextTreeToRead
          call Galacticus_State_Store(message)
          ! We do have queued trees - simply return the next one.
          mergerTreeQueuePosition=mergerTreeQueuePosition+1
          thisTree=mergerTreesQueued(mergerTreeQueuePosition)
          ! If the last tree in the queue has been reached, destroy the queue.
          if (mergerTreeQueuePosition == size(mergerTreesQueued)) then
             call Memory_Usage_Record(sizeof(mergerTreesQueued),addRemove=-1)
             deallocate(mergerTreesQueued)
             mergerTreeQueuePosition=0
          end if
          finished=.true.
       else
          ! We do not have any queued trees. Therefore, increment the tree to read index.
          nextTreeToRead=nextTreeToRead+1
          ! Keep incrementing the tree index until we find the first tree to process (if we haven't done so already). Also skip trees
          ! that contain 1 or fewer nodes and these are unprocessable.
          if (nextTreeToRead <= size(mergerTreeFirstNodeIndex)) then
             do while ( (       mergerTreeReadBeginAt            > 0                     &
                  &       .and. mergerTreeIndex(nextTreeToRead) /= mergerTreeReadBeginAt &
                  &     )                                                                &
                  &     .or. mergerTreeNodeCount(nextTreeToRead) <= 1                    &
                  &   )
                nextTreeToRead=nextTreeToRead+1
                ! If the end of the list has been reached, exit.
                if (nextTreeToRead > size(mergerTreeFirstNodeIndex)) exit
             end do
          end if

          ! Flag that we've now found the first merger tree to process.
          mergerTreeReadBeginAt=-1
          if (nextTreeToRead > size(mergerTreeFirstNodeIndex))  then
             ! All trees have been read.
             !$omp critical(HDF5_Access)
             ! Close the halo trees group.
             if (haloTreesGroup%isOpen()) call haloTreesGroup%close()
             ! Test if the merger tree file is still open.
             if (mergerTreeFile%isOpen()) then
                ! Close the particles group.
                if (mergerTreeFile%hasGroup("particles")) then
                   if (particlesGroup%isOpen()) call particlesGroup%close()
                end if
                ! Close the file.
                call mergerTreeFile%close()
             end if
             !$omp end critical(HDF5_Access)
             ! Flag that we do not have a tree.
             haveTree=.false.
          else
             ! Flag that we do have a tree.
             haveTree=.true.
          end if

          ! Continue only if we have a tree.
          if (haveTree) then
             ! Retrieve stored internal state if possible.
             call Galacticus_State_Retrieve
             ! Take a snapshot of the internal state and store it.
             call Galacticus_State_Snapshot
             message='Storing state for tree #'
             message=message//nextTreeToRead
             call Galacticus_State_Store(message)

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

                ! Sort node indices.
                call Create_Node_Indices(nodes,nodeCount)

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
                if (mergerTreeReadPresetPositions.or.mergerTreeReadPresetOrbits) then
                   ! Position and velocity methods are required.
                   if (.not.(associated(Tree_Node_Position_Set).and.associated(Tree_Node_Velocity_Set))) call&
                        & Galacticus_Error_Report('Merger_Tree_Read_Do','presetting positions or orbits requires a component that supports&
                        & position and velocity setting (e.g. set [treeNodeMethodPosition]=preset); alternatively setting [mergerTreeReadPresetPositions]=false and [mergerTreeReadPresetOrbits]=false will remove the need to store positions and velocities')
                end if
                if (mergerTreeReadPresetMergerTimes) then
                   ! Time of merging property is required.
                   if (.not.associated(Tree_Node_Satellite_Time_of_Merging_Set)) call Galacticus_Error_Report('Merger_Tree_Read_Do',&
                        & 'presetting merging times requires a component that supports setting of merging times')
                end if
                if (mergerTreeReadPresetScaleRadii) then
                   ! Scale radius property is required.
                   if (.not.associated(Tree_Node_Dark_Matter_Profile_Scale_Set)) call Galacticus_Error_Report('Merger_Tree_Read_Do',&
                        & 'presetting scale radii requires a component that supports setting of scale radii')
                end if
                if (mergerTreeReadPresetSpins      ) then
                   ! Spin property is required.
                   if (.not.associated(Tree_Node_Spin_Set                     )) call Galacticus_Error_Report('Merger_Tree_Read_Do',&
                        & 'presetting spins requires a component that supports setting of spins')
                end if
                if (mergerTreeReadPresetOrbits     ) then
                   ! Orbit property is required.
                   if (.not.associated(Tree_Node_Satellite_Virial_Orbit_Set   )) call Galacticus_Error_Report('Merger_Tree_Read_Do',&
                        & 'presetting orbits requires a component that supports setting of orbits (e.g. [treeNodeMethodSatelliteOrbit]=preset); alternatively, set [mergerTreeReadPresetOrbits]=false to prevent attempts to set orbits')
                end if

                ! Assign scale radii.
                if (mergerTreeReadPresetScaleRadii) call Assign_Scale_Radii    (nodes,thisNodeList)

                ! Assign spin parameters.
                if (mergerTreeReadPresetSpins     ) call Assign_Spin_Parameters(nodes,thisNodeList)

                ! Assign isolated node indices to subhalos.
                call Assign_Isolated_Node_Indices(nodes,thisNodeList)

                ! Ensure that isolated nodes with progenitors that descend into subhalos have valid primary progenitors.
                call Validate_Isolated_Halos(nodes)

                ! Scan subhalos to determine when and how they merge.
                call Scan_For_Mergers(nodes,thisNodeList,historyCountMaximum)

                ! Search for any nodes which were flagged as merging with another node and assign appropriate pointers.
                call Assign_Mergers(nodes,thisNodeList)

                ! Allocate arrays for history building.
                if (historyCountMaximum > 0) then
                   if (allocated(position)) call Dealloc_Array(position)
                   if (allocated(velocity)) call Dealloc_Array(velocity)
                   call Alloc_Array(historyTime,[int(historyCountMaximum)])
                   if (mergerTreeReadPresetSubhaloMasses                              ) call Alloc_Array(historyMass,[  int(historyCountMaximum)])
                   if (mergerTreeReadPresetPositions    .or.mergerTreeReadPresetOrbits) call Alloc_Array(position   ,[3,int(historyCountMaximum)])
                   if (mergerTreeReadPresetPositions    .or.mergerTreeReadPresetOrbits) call Alloc_Array(velocity   ,[3,int(historyCountMaximum)])
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

                ! Destroy sorted node indices.
                call Destroy_Node_Indices()
             end if
          end if

          ! If the main gtree is to be skipped, destroy it and don't set the finished flag so that we'll
          ! make another pass through the main tree getting loop and pick up a queued tree instead.
          if (mergerTreeReadSkipMainTree) then
             call thisTree%destroy()
          else
             finished=.true.
          end if

       end if
    end do
    !$omp end critical(mergerTreeReadTree)
    return
  end subroutine Merger_Tree_Read_Do

  subroutine Read_Node_Data(nodes,firstNodeIndex,nodeCount)
    !% Read node data from an HDF5 file.
    use Galacticus_Error
    use Cosmology_Functions
    use Arrays_Search
    use Numerical_Comparison
    use Vectors
    use Memory_Management
    implicit none
    type(nodeData),          intent(inout), dimension(:)   :: nodes
    integer(kind=HSIZE_T),   intent(in),    dimension(1)   :: nodeCount,firstNodeIndex
    double precision,        allocatable,   dimension(:,:) :: angularMomentum
    integer(kind=kind_int8)                                :: iNode
    integer                                                :: iOutput,scaleFactorExponentAngularMomentum

    !$omp critical(HDF5_Access)
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
    ! Scale or half-mass radius.
    if (mergerTreeReadPresetScaleRadii) then
       haveScaleRadii=haloTreesGroup%hasDataset("scaleRadius")
       if (haveScaleRadii) then
          call haloTreesGroup%readDatasetStatic("scaleRadius",nodes%scaleRadius,firstNodeIndex,nodeCount)
          nodes%scaleRadius=nodes%scaleRadius*unitConversionLength
          if (scaleFactorExponentLength /= 0) then
             do iNode=1,nodeCount(1)
                nodes(iNode)%scaleRadius=nodes(iNode)%scaleRadius*Expansion_Factor(nodes(iNode)%nodeTime)**scaleFactorExponentLength
             end do
          end if
       else
          call haloTreesGroup%readDatasetStatic("halfMassRadius",nodes%halfMassRadius,firstNodeIndex,nodeCount)
          nodes%halfMassRadius=nodes%halfMassRadius*unitConversionLength
          if (scaleFactorExponentLength /= 0) then
             do iNode=1,nodeCount(1)
                nodes(iNode)%halfMassRadius=nodes(iNode)%halfMassRadius*Expansion_Factor(nodes(iNode)%nodeTime)**scaleFactorExponentLength
             end do
          end if
       end if
    end if
    ! Halo spin.
    if (mergerTreeReadPresetSpins     ) then
       call haloTreesGroup%readDataset("angularMomentum",angularMomentum,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
       angularMomentum=angularMomentum*unitConversionLength*unitConversionVelocity*unitConversionMass
       scaleFactorExponentAngularMomentum=scaleFactorExponentLength+scaleFactorExponentVelocity+scaleFactorExponentMass
       if (scaleFactorExponentAngularMomentum /= 0) then
          do iNode=1,nodeCount(1)
             angularMomentum(:,iNode)=angularMomentum(:,iNode)*Expansion_Factor(nodes(iNode)%nodeTime)**scaleFactorExponentAngularMomentum
          end do
       end if
       ! Transfer to nodes.
       forall(iNode=1:nodeCount(1))
          nodes(iNode)%angularMomentum=Vector_Magnitude(angularMomentum(:,iNode))
       end forall
       call Dealloc_Array(angularMomentum)
    end if
    !$omp end critical(HDF5_Access)

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

  subroutine Create_Node_Indices(nodes,nodeCount)
    !% Create a sorted list of node indices with an index into the original array.
    use Galacticus_Display
    use String_Handling
    use Memory_Management
    use Sort
    implicit none
    type(nodeData),          intent(inout), dimension(:)   :: nodes
    integer(kind=HSIZE_T),   intent(in),    dimension(1)   :: nodeCount
    integer                                                :: iNode
    type(varying_string)                                   :: message

    ! Build a sorted list of node indices with an index into the original arrays.
    call Alloc_Array(nodeLocations          ,int(nodeCount))
    call Alloc_Array(nodeIndicesSorted      ,int(nodeCount))
    call Alloc_Array(descendentLocations    ,int(nodeCount))
    call Alloc_Array(descendentIndicesSorted,int(nodeCount))
    nodeLocations      =Sort_Index_Do(nodes%nodeIndex      )
    descendentLocations=Sort_Index_Do(nodes%descendentIndex)
    forall (iNode=1:int(nodeCount(1)))
       nodeIndicesSorted      (iNode)=nodes(nodeLocations      (iNode))%nodeIndex
       descendentIndicesSorted(iNode)=nodes(descendentLocations(iNode))%descendentIndex
    end forall
    if (Galacticus_Verbosity_Level() >= verbosityWarn) then
       do iNode=2,int(nodeCount(1))
          if (nodeIndicesSorted(iNode) == nodeIndicesSorted(iNode-1)) then
             if (Galacticus_Verbosity_Level() >= verbosityWarn) then
                message="duplicate node index found ["
                message=message//nodeIndicesSorted(iNode)//']'
                call Galacticus_Display_Message(message,verbosityWarn)
             end if
          end if
       end do
    end if
    return
  end subroutine Create_Node_Indices

  function Node_Location(nodeIndex)
    !% Return the location in the original array of the given {\tt nodeIndex}.
    use Arrays_Search
    implicit none
    integer(kind=kind_int8)             :: Node_Location
    integer(kind=kind_int8), intent(in) :: nodeIndex
    integer(kind=kind_int8)             :: iNode

    iNode=Search_Array(nodeIndicesSorted,nodeIndex)
    if (iNode < 1 .or. iNode > size(nodeIndicesSorted)) then
       Node_Location=1
    else
       Node_Location=nodeLocations(iNode)
    end if
    return
  end function Node_Location

  integer function Descendent_Node_Sort_Index(descendentIndex)
    !% Return the sort index of the given {\tt descendentIndex}.
    use Arrays_Search
    implicit none
    integer(kind=kind_int8), intent(in) :: descendentIndex

    Descendent_Node_Sort_Index=Search_Array(descendentIndicesSorted,descendentIndex)
    return
  end function Descendent_Node_Sort_Index

  subroutine Destroy_Node_Indices()
    !% Destroy the sorted list of node indices.
    use Memory_Management
    implicit none

    if (allocated(nodeLocations          )) call Dealloc_Array(nodeLocations          )
    if (allocated(nodeIndicesSorted      )) call Dealloc_Array(nodeIndicesSorted      )
    if (allocated(descendentLocations    )) call Dealloc_Array(descendentLocations    )
    if (allocated(descendentIndicesSorted)) call Dealloc_Array(descendentIndicesSorted)
    return
  end subroutine Destroy_Node_Indices

  subroutine Read_Particle_Data(nodes,firstNodeIndex,nodeCount,position,velocity)
    !% Read data on particle positions/velocities.
    use Galacticus_Error
    use Cosmology_Functions
    implicit none
    type(nodeData),        intent(inout), dimension(:)                :: nodes
    integer(kind=HSIZE_T), intent(in),    dimension(1)                :: nodeCount,firstNodeIndex
    double precision,      intent(inout), dimension(:,:), allocatable :: position,velocity
    integer(kind=HSIZE_T)                                             :: iNode

    ! Initial particle data to null values.
    nodes%particleIndexStart=-1_kind_int8
    nodes%particleIndexCount=-1_kind_int8

    if (mergerTreeReadPresetPositions.or.mergerTreeReadPresetOrbits) then
       !$omp critical(HDF5_Access)
       ! position.
       call haloTreesGroup%readDataset("position",position,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
       ! velocity.
       call haloTreesGroup%readDataset("velocity",velocity,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
       ! If a set of most bound particle indices are present, read them.
       if (haloTreesGroup%hasDataset("particleIndexStart").and.haloTreesGroup%hasDataset("particleIndexCount")) then
          call haloTreesGroup%readDatasetStatic("particleIndexStart",nodes%particleIndexStart,firstNodeIndex,nodeCount)
          call haloTreesGroup%readDatasetStatic("particleIndexCount",nodes%particleIndexCount,firstNodeIndex,nodeCount)
       end if
       !$omp end critical(HDF5_Access)
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
    use Galacticus_Display
    implicit none
    type(nodeData),          intent(inout), dimension(:), target :: nodes
    integer(kind=kind_int8)                                      :: iNode,nodeLocation
    type(varying_string)                                         :: message

    do iNode=1,size(nodes)
       if (nodes(iNode)%descendentIndex >= 0) then
          nodeLocation=Node_Location(nodes(iNode)%descendentIndex)
          if (nodes(nodeLocation)%nodeIndex /= nodes(iNode)%descendentIndex) then
             message='failed to find descendent node: '
             message=message//nodes(iNode)%descendentIndex//' of '//nodes(iNode)%nodeIndex
             call Galacticus_Error_Report('Build_Descendent_Pointers',message)
          end if
          nodes(iNode)%descendentNode => nodes(nodeLocation)
       else
          nodes(iNode)%descendentNode => null()
       end if
       if (nodes(iNode)%hostIndex       >= 0) then
          nodeLocation=Node_Location(nodes(iNode)%hostIndex      )
          if (nodes(nodeLocation)%nodeIndex /= nodes(iNode)%hostIndex       ) then
             message='failed to find host node: '
             message=message//nodes(iNode)%hostIndex//' of '//nodes(iNode)%nodeIndex
             if (mergerTreeReadMissingHostsAreFatal) then
                call Galacticus_Error_Report('Build_Descendent_Pointers',message)
             else
                message=message//" - resetting this node to be an isolated node"
                call Galacticus_Display_Message(message,verbosity=verbosityInfo)
                nodes(iNode)%hostIndex =  nodes(iNode)%nodeIndex
                nodes(iNode)%hostNode  => nodes(iNode)
                nodes(iNode)%isSubhalo =  .false.
             end if
          else
             nodes(iNode)%hostNode    => nodes(nodeLocation)
          end if
       else
          nodes(iNode)%hostNode       => null()
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
    integer(kind=kind_int8)                              :: iNode,descendentLocation
    logical                                              :: failed,isolatedProgenitorExists,descendentsFound
    integer                                              :: descendentIndex
    type(varying_string)                                 :: message

    do iNode=1,size(nodes)
       if (nodes(iNode)%isSubhalo) then
          descendentNode => nodes(iNode)%descendentNode
          do while (associated(descendentNode))
             ! Is this node isolated?
             if (.not.descendentNode%isSubhalo) then
                ! Check if there is any isolated node which descends into this node.
                descendentIndex=Descendent_Node_Sort_Index(descendentNode%nodeIndex)
                if (descendentIndex > 0 .and. descendentIndex <= size(nodes)) then
                   descendentLocation=descendentLocations(descendentIndex)
                   descendentsFound=(nodes(descendentLocation)%descendentIndex == descendentNode%nodeIndex)
                else
                   descendentsFound=.false.
                end if
                if (descendentsFound) then
                   isolatedProgenitorExists=.false.
                   do while (nodes(descendentLocation)%descendentIndex == descendentNode%nodeIndex .and. descendentIndex > 0 .and. .not.isolatedProgenitorExists)
                      isolatedProgenitorExists=(nodes(descendentLocation)%nodeIndex == nodes(descendentLocation)%hostIndex) 
                      descendentIndex   =descendentIndex-1
                      descendentLocation=descendentLocations(descendentIndex)
                   end do
                end if
                if (.not.descendentsFound .or. .not.isolatedProgenitorExists) then
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
       descendentIndex=Descendent_Node_Sort_Index(nodes(iNode)%nodeIndex)
       if (descendentIndex > 0 .and. descendentIndex <= size(nodes)) then
          descendentLocation=descendentLocations(descendentIndex)
          descendentsFound=(nodes(descendentLocation)%descendentIndex == nodes(iNode)%nodeIndex)
       else
          descendentsFound=.false.
       end if
       if (descendentsFound) then
          isolatedProgenitorExists=.false.
          do while (nodes(descendentLocation)%descendentIndex == nodes(iNode)%nodeIndex .and. descendentIndex > 0 .and. .not.isolatedProgenitorExists)
             isolatedProgenitorExists=(nodes(descendentLocation)%nodeIndex == nodes(descendentLocation)%hostIndex) 
             descendentIndex   =descendentIndex-1
             descendentLocation=descendentLocations(descendentIndex)
          end do
       end if
       if (descendentsFound .and. .not.isolatedProgenitorExists) then
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
    use Merger_Tree_Read_State
    implicit none
    type(nodeData), intent(inout), dimension(:), target  :: nodes
    integer,        intent(out)                          :: iExtraTree,primaryRootIndex
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

    ! Allocate a tree queue for any extra roots.
    if (treeRootCount > 1) then
       if (allocated(mergerTreesQueued)) then
          call Memory_Usage_Record(sizeof(mergerTreesQueued),addRemove=-1)
          deallocate(mergerTreesQueued)
       end if
       allocate(mergerTreesQueued(treeRootCount-1))
       call Memory_Usage_Record(sizeof(mergerTreesQueued))
       iExtraTree=0
    end if
    return
  end subroutine Assign_Nodes_To_Trees

  subroutine Build_Isolated_Parent_Pointers(thisTree,nodes,nodeList,primaryRootIndex,iExtraTree)
    !% Create parent pointer links between isolated nodes and assign times and masses to those nodes.
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Error
    implicit none
    type(mergerTree),   intent(inout)               :: thisTree
    type(nodeData),     intent(inout), dimension(:) :: nodes
    type(treeNodeList), intent(inout), dimension(:) :: nodeList
    integer,            intent(inout)               :: iExtraTree,primaryRootIndex
    integer                                         :: iNode,iIsolatedNode
    type(varying_string)                            :: message
    character(len=12)                               :: label

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
                if (mergerTreeReadTreeIndexToRootNodeIndex) thisTree%index = nodes(iNode)%nodeIndex
                thisTree                     %volumeWeight=  treeVolumeWeightCurrent
             else
                iExtraTree=iExtraTree+1
                ! Catch attempts to queue the root node of the tree.
                if (nodes(iNode)%nodeIndex == thisTree%index .and. .not.mergerTreeReadTreeIndexToRootNodeIndex) then
                   message=         "the root node of a tree is being placed on the tree queue,"       //char(10)
                   message=message//"but [mergerTreeReadTreeIndexToRootNodeIndex]=false."              //char(10)
                   message=message//"This will result in duplicated tree indices."                     //char(10)
                   message=message//"Recommendation: set [mergerTreeReadTreeIndexToRootNodeIndex]=true"//char(10)
                   call Galacticus_Error_Report('Build_Isolated_Parent_Pointers',message)
                end if
                mergerTreesQueued(iExtraTree)%baseNode    => nodeList(iIsolatedNode)%node
                mergerTreesQueued(iExtraTree)%index       =  nodes(iNode)%nodeIndex
                mergerTreesQueued(iExtraTree)%volumeWeight=  treeVolumeWeightCurrent
                mergerTreesQueued(iExtraTree)%initialized =  .false.
             end if
          end if
          if (nodes(iNode)%nodeMass <= 0.0d0) then
             write (label,'(e12.6)') nodes(iNode)%nodeMass
             message='non-positive mass ['//label//'] found for node '
             message=message//nodeList(iIsolatedNode)%node%index()
             call Galacticus_Error_Report('Build_Isolated_Parent_Pointers',message)
          end if
          if (nodes(iNode)%nodeTime <= 0.0d0) then
             write (label,'(e12.6)') nodes(iNode)%nodeTime
             message='non-positive time ['//label//'] found for node '
             message=message//nodeList(iIsolatedNode)%node%index()
             call Galacticus_Error_Report('Build_Isolated_Parent_Pointers',message)
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

  subroutine Assign_Scale_Radii(nodes,nodeList)
    !% Assign scale radii to nodes.
    use Root_Finder
    use Dark_Matter_Halo_Scales
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Display
    use Galacticus_Error
    implicit none
    type(nodeData),          intent(inout), dimension(:) :: nodes
    type(treeNodeList),      intent(inout), dimension(:) :: nodeList
    type(fgsl_function),     save                        :: rootFunction
    type(fgsl_root_fsolver), save                        :: rootFunctionSolver
    double precision,        parameter                   :: toleranceAbsolute=1.0d-9,toleranceRelative=1.0d-9,scaleRadiusMaximumAllowed=100.0d0
    logical,                 save                        :: excessiveScaleRadiiReported=.false.
    type(c_ptr)                                          :: parameterPointer
    integer                                              :: iNode,iIsolatedNode
    double precision                                     :: radiusMinimum,radiusMaximum,radiusScale
    logical                                              :: excessiveScaleRadii,excessiveHalfMassRadii
    character(len=50)                                    :: message

    iIsolatedNode=0
    excessiveScaleRadii   =.false.
    excessiveHalfMassRadii=.false.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
          iIsolatedNode=iIsolatedNode+1

          ! Check if the node is sufficiently massive.
          if (Tree_Node_Mass(nodeList(iIsolatedNode)%node) >= mergerTreeReadPresetScaleRadiiMinimumMass) then

             ! Check if we have scale radii read directly from file.
             if (haveScaleRadii) then
                ! We do, so simply use them to set the scale radii in tree nodes.
                call Tree_Node_Dark_Matter_Profile_Scale_Set(nodeList(iIsolatedNode)%node,nodes(iNode)%scaleRadius)

             else
                ! We do not have scale radii read directly. Instead, compute them from half-mass radii.

                ! Set the active node and target half mass radius.
                activeNode => nodeList(iIsolatedNode)%node
                halfMassRadius=nodes(iNode)%halfMassRadius

                ! Find minimum and maximum scale radii such that the target half-mass radius is bracketed.
                radiusMinimum=halfMassRadius
                radiusMaximum=halfMassRadius
                do while (Half_Mass_Radius_Root(radiusMinimum,parameterPointer) <= 0.0d0)
                   radiusMinimum=0.5d0*radiusMinimum
                end do
                do while (Half_Mass_Radius_Root(radiusMaximum,parameterPointer) >= 0.0d0)
                   if (radiusMaximum > scaleRadiusMaximumAllowed*Dark_Matter_Halo_Virial_Radius(activeNode)) then
                      write (message,'(a,i10,a)') 'node ',activeNode%index(),' has unphysical scale-radius:'
                      call Galacticus_Display_Message(trim(message))
                      write (message,'(a,e12.6,a)') ' half mass radius = ',halfMassRadius                            ,' Mpc'
                      call Galacticus_Display_Message(trim(message))
                      write (message,'(a,e12.6,a)') '    virial radius = ',Dark_Matter_Halo_Virial_Radius(activeNode),' Mpc'
                      call Galacticus_Display_Message(trim(message))
                      call Galacticus_Error_Report('Assign_Scale_Radii','scale radius greatly exceeds virial radius - this seems unphysical')
                   end if
                   radiusMaximum=2.0d0*radiusMaximum
                end do

                ! Solve for the scale radius.
                radiusScale=Root_Find(radiusMinimum,radiusMaximum,Half_Mass_Radius_Root,parameterPointer &
                     &,rootFunction,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
                call Tree_Node_Dark_Matter_Profile_Scale_Set(nodeList(iIsolatedNode)%node,radiusScale)

                ! Check for scale radii exceeding the virial radius.
                if (radiusScale    > Dark_Matter_Halo_Virial_Radius(activeNode)) excessiveScaleRadii   =.true.

                ! Check for half-mass radii exceeding the virial radius.
                if (halfMassRadius > Dark_Matter_Halo_Virial_Radius(activeNode)) excessiveHalfMassRadii=.true.
             end if

          end if
       end if
    end do

    ! Report warning on excessive scale radii if not already done.
    if (excessiveScaleRadii.and..not.excessiveScaleRadiiReported) then
       excessiveScaleRadiiReported=.true.
       call Galacticus_Display_Message('Assign_Scale_Radii(): warning - some scale radii exceed the corresponding virial radii - suggests&
            & inconsistent definitions of halo mass/radius',verbosityWarn)
    end if

    ! Exit on excessive half mass radii.
    if (excessiveHalfMassRadii) call Galacticus_Error_Report('Assign_Scale_Radii','some half mass radii exceed corresponding virial radii')

    return
  end subroutine Assign_Scale_Radii

  subroutine Assign_Spin_Parameters(nodes,nodeList)
    !% Assign spin parameters to nodes.
    use Numerical_Constants_Physical
    use Dark_Matter_Profiles
    implicit none
    type(nodeData),          intent(inout), dimension(:) :: nodes
    type(treeNodeList),      intent(inout), dimension(:) :: nodeList
    integer                                              :: iNode,iIsolatedNode
    double precision                                     :: spin

    iIsolatedNode=0
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
          iIsolatedNode=iIsolatedNode+1
          ! Compute the spin parameter.
          spin=                                         nodes(iNode)%angularMomentum           &
               & *dsqrt(dabs(Dark_Matter_Profile_Energy(nodeList(iIsolatedNode)%node)))        &
               & /gravitationalConstantGalacticus                                              &
               & /Tree_Node_Mass(                       nodeList(iIsolatedNode)%node  )**2.5d0
          ! Assign to the node.
          call Tree_Node_Spin_Set(nodeList(iIsolatedNode)%node,spin)
       end if
    end do
    return
  end subroutine Assign_Spin_Parameters

  function Half_Mass_Radius_Root(radius,parameterPointer) bind(c)
    !% Function used to find scale radius of dark matter halos given their half-mass radius.
    use, intrinsic :: ISO_C_Binding
    use Dark_Matter_Profiles
    implicit none
    real(c_double)        :: Half_Mass_Radius_Root
    real(c_double), value :: radius
    type(c_ptr)           :: parameterPointer

    ! Set scale radius to current guess.
    call Tree_Node_Dark_Matter_Profile_Scale_Set(activeNode,radius)
    ! Compute difference between mass fraction enclosed at half mass radius and one half.
    Half_Mass_Radius_Root=Dark_Matter_Profile_Enclosed_Mass(activeNode,halfMassRadius)/Tree_Node_Mass(activeNode)-0.50d0
    return
  end function Half_Mass_Radius_Root

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
    end do
    return
  end subroutine Assign_Isolated_Node_Indices

  subroutine Scan_For_Mergers(nodes,nodeList,historyCountMaximum)
    !% Scan for and record mergers between nodes.
    use Vectors
    use Virial_Orbits
    use Kepler_Orbits_Structure
    use Tree_Nodes
    use Cosmology_Functions
    use Dark_Matter_Halo_Scales
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Error
    use Numerical_Comparison
    implicit none
    type(nodeData),          intent(inout), dimension(:) :: nodes
    type(treeNodeList),      intent(inout), dimension(:) :: nodeList
    integer(kind=kind_int8), intent(out)                 :: historyCountMaximum
    type(nodeData),          pointer                     :: thisNode
    type(treeNode),          pointer                     :: firstProgenitor,satelliteNode,hostNode
    double precision,        parameter                   :: timeUntilMergingInfinite=1.0d30
    double precision,                       dimension(3) :: satellitePosition,hostPosition,relativePosition
    double precision,                       dimension(3) :: satelliteVelocity,hostVelocity,relativeVelocity
    logical,                 parameter                   :: acceptUnboundOrbits=.false.
    type(keplerOrbit)                                    :: thisOrbit
    integer                                              :: iNode,iIsolatedNode,descendentIndex
    integer(kind=kind_int8)                              :: historyCount,descendentLocation
    logical                                              :: branchMerges,branchTipReached,endOfBranch,nodeWillMerge,descendentsFound
    double precision                                     :: timeSubhaloMerges,radiusPericenter,radiusApocenter,radiusVirial
    type(varying_string)                                 :: message

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
             ! Select the subset which have a subhalo.
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
                      branchTipReached            =.true.
                      endOfBranch                 =.true.
                      historyCount                =historyCount+max(0_kind_int8,thisNode%particleIndexCount)
                   else if (.not.thisNode%descendentNode%isSubhalo) then
                      ! Descendent is not a subhalo, treat as a merging event.  
                      branchMerges                =.true.
                      endOfBranch                 =.true.
                      nodes(iNode)%mergesWithIndex=thisNode%descendentNode%nodeIndex
                      historyCount                =historyCount+max(0_kind_int8,thisNode%particleIndexCount)
                      thisNode                    => thisNode%descendentNode
                   else
                      ! Merges with another subhalo.
                      descendentIndex=Descendent_Node_Sort_Index(thisNode%descendentNode%nodeIndex)
                      if (descendentIndex > 0 .and. descendentIndex <= size(nodes)) then
                         descendentLocation=descendentLocations(descendentIndex)
                         descendentsFound=(nodes(descendentLocation)%descendentNode%nodeIndex == thisNode%descendentNode%nodeIndex)
                      else
                         descendentsFound=.false.
                      end if
                      if (descendentsFound) then
                         do while (nodes(descendentLocation)%descendentNode%nodeIndex == thisNode%descendentNode%nodeIndex .and. descendentIndex > 0)
                            if     (                                                                                     &
                                 &                    nodes(descendentLocation)%nodeIndex         /= thisNode%nodeIndex  &
                                 &  .and.             nodes(descendentLocation)%isolatedNodeIndex /= nodeIsUnreachable   &
                                 &  .and.  associated(nodes(descendentLocation)%descendentNode                         ) &
                                 &  .and.             nodes(descendentLocation)%nodeMass           > thisNode%nodeMass   &
                                 & ) then
                               ! Another node mergers into current node's descendent subhalo and is more massive than current
                               ! node. Therefore, class this as a subhalo-subhalo merger.
                               branchMerges                =.true.                                  
                               endOfBranch                 =.true.
                               nodes(iNode)%mergesWithIndex=nodes(descendentLocation)%descendentNode%nodeIndex
                               historyCount                =historyCount+max(0_kind_int8,thisNode%particleIndexCount)
                               thisNode                    => thisNode%descendentNode
                               exit
                            end if
                            do while (descendentIndex > 0)
                               descendentIndex   =descendentIndex-1
                               if (descendentIndex <= 0) exit
                               descendentLocation=descendentLocations(descendentIndex)
                               if (associated(nodes(descendentLocation)%descendentNode)) exit
                            end do
                            if (descendentIndex == 0) exit
                         end do
                      end if
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
                historyCountMaximum=max(historyCountMaximum,max(0_kind_int8,nodes(iNode)%particleIndexCount))
             end if
             ! Set a merging time and/or orbit if this node will merge.
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
             ! Detect if the node parent has no isolated child - in which case one will have been made for it using a
             ! direct copy of itself. Note that this is an ugly solution - once trees can handle nodes with no primary
             ! progenitor (but with secondary progenitors) a cleaner test could be used here.
             if (associated(nodeList(iIsolatedNode)%node%childNode)) then
                if     (                                                                                                                                                       &
                     &                                    nodeList(iIsolatedNode)%node%index() ==                nodeList(iIsolatedNode)%node%childNode%index()                &
                     &  .and. Values_Agree(Tree_Node_Time(nodeList(iIsolatedNode)%node)        ,  Tree_Node_Time(nodeList(iIsolatedNode)%node%childNode)       ,relTol=2.0d-6) &
                     & ) then
                   ! Set the position and velocity of the pseudo-primary progenitor here also.
                   call Tree_Node_Position_Set(nodeList(iIsolatedNode)%node%childNode,nodes(iNode)%position)
                   call Tree_Node_Velocity_Set(nodeList(iIsolatedNode)%node%childNode,nodes(iNode)%velocity)
                end if
             end if
          end if
       end if
    end do

    ! Set orbits.
    if (mergerTreeReadPresetOrbits) then
       iIsolatedNode=0
       do iNode=1,size(nodes)
          if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
             iIsolatedNode=iIsolatedNode+1
             ! Set the orbit for this halo.
             satelliteNode => nodeList(iIsolatedNode)%node
             if (associated(satelliteNode%parentNode).and..not.satelliteNode%isPrimaryProgenitor()) then
                ! Find the orbital partner.
                hostNode => satelliteNode%parentNode%childNode
                ! Get position and velocity.
                call Tree_Node_Position(satelliteNode,satellitePosition)
                call Tree_Node_Velocity(satelliteNode,satelliteVelocity)
                call Tree_Node_Position(hostNode     ,hostPosition     )
                call Tree_Node_Velocity(hostNode     ,hostVelocity     )
                ! Find relative position and velocity.
                relativePosition=satellitePosition-hostPosition
                relativeVelocity=satelliteVelocity-hostVelocity
                ! Catch zero separation halos.
                if (Vector_Magnitude(relativePosition) == 0.0d0) then
                   message='merging halos ['
                   message=message//satelliteNode%index()//' & '//hostNode%index()//'] have zero separation'
                   call Galacticus_Error_Report('Scan_For_Mergers',message)
                end if
                ! Account for periodicity.
                if (simulationIsPeriodic) then
                   relativePosition=mod(relativePosition+0.5d0*lengthSimulationBox,lengthSimulationBox)-0.5d0*lengthSimulationBox
                   relativePosition=mod(relativePosition-0.5d0*lengthSimulationBox,lengthSimulationBox)+0.5d0*lengthSimulationBox
                end if
                ! Account for Hubble flow.
                if (.not.velocitiesIncludeHubbleFlow) relativeVelocity=relativeVelocity+relativePosition*Hubble_Parameter(tCosmological=Tree_Node_Time(satelliteNode))
                ! Create the orbit.
                call thisOrbit%reset()
                call thisOrbit%massesSet            (                                     &
                     &                               Tree_Node_Mass(satelliteNode      ), &
                     &                               Tree_Node_Mass(hostNode%parentNode)  &
                     &                              )
                call thisOrbit%radiusSet            (                                                                   Vector_Magnitude(relativePosition))
                call thisOrbit%velocityRadialSet    (                 Dot_Product   (relativeVelocity,relativePosition)/Vector_Magnitude(relativePosition))
                call thisOrbit%velocityTangentialSet(Vector_Magnitude(Vector_Product(relativeVelocity,relativePosition)/Vector_Magnitude(relativePosition)))
                ! Propagate to the virial radius.
                radiusPericenter=thisOrbit%radiusPericenter()
                radiusApocenter =thisOrbit%radiusApocenter ()
                radiusVirial    =Dark_Matter_Halo_Virial_Radius(hostNode%parentNode)
                ! Check if the orbit intersects the virial radius.
                if     (                                                                          &
                     &    radiusVirial >= radiusPericenter                                        &
                     &  .and.                                                                     &
                     &   (radiusVirial <= radiusApocenter          .or. .not.thisOrbit%isBound()) &
                     &  .and.                                                                     &
                     &   (.not.mergerTreeReadPresetOrbitsBoundOnly .or.      thisOrbit%isBound()) &
                     & ) then
                   call thisOrbit%propagate(radiusVirial,infalling=.true.)
                   ! Set the orbit.
                   call Tree_Node_Satellite_Virial_Orbit_Set(satelliteNode,thisOrbit)
                else if (mergerTreeReadPresetOrbitsSetAll) then
                   ! The given orbit does not cross the virial radius. Since all orbits must be set, choose an orbit at random.
                   thisOrbit=Virial_Orbital_Parameters(satelliteNode,hostNode,acceptUnboundOrbits)
                   call Tree_Node_Satellite_Virial_Orbit_Set(satelliteNode,thisOrbit)
                else if (mergerTreeReadPresetOrbitsAssertAllSet) then
                   message='virial orbit could not be set for node '
                   message=message//satelliteNode%index()//char(10)
                   message=message//' -> set [mergerTreeReadPresetOrbitsAssertAllSet]=false to ignore this problem'//char(10)
                   message=message//'    (this may lead to other problems)'
                   call Galacticus_Error_Report('Scan_For_Mergers',message)
                end if
             end if
          end if
       end do
    end if

    return
  end subroutine Scan_For_Mergers

  subroutine Assign_Mergers(nodes,nodeList)
    !% Assign pointers to merge targets.
    use Galacticus_Error
    use String_Handling
    implicit none
    type(nodeData),     intent(inout), dimension(:) :: nodes
    type(treeNodeList), intent(inout), dimension(:) :: nodeList
    type(treeNode),     pointer                     :: rootNode,mergeRootNode
    integer                                         :: iNode
    integer(kind=kind_int8)                         :: jNode
    type(varying_string)                            :: message

    if (mergerTreeReadPresetMergerNodes) then
       do iNode=1,size(nodes)
          ! Check if this node was flagged as merging with another node.
          if (nodes(iNode)%mergesWithIndex /= -1) then
             ! Search for the node that it merges with.
             jNode=Node_Location(nodes(iNode)%mergesWithIndex)
             if (nodes(jNode)%isolatedNodeIndex <= 0) then
                ! This node does not belong to any isolated halo - this should not happen.
                message='subhalo-subhalo ['
                message=message//nodes(iNode)%nodeIndex//":"//nodes(jNode)%nodeIndex
                message=message//'] merger in which subhalo has no isolated node progenitor - this should not happen'
                call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
             else
                ! Find the root nodes into which these nodes will descend.
                rootNode => nodeList(nodes(iNode)%isolatedNodeIndex)%node
                do while (associated(rootNode%parentNode))
                   rootNode => rootNode%parentNode
                end do
                mergeRootNode => nodeList(nodes(jNode)%isolatedNodeIndex)%node
                do while (associated(mergeRootNode%parentNode))
                   mergeRootNode => mergeRootNode%parentNode
                end do
                ! Set a pointer between the isolated nodes corresponding to these subhalos if and only if they descend
                ! into the same root node.
                if (rootNode%index() == mergeRootNode%index()) then
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
    integer(kind=kind_int8)                                   :: iTime,historyCount,descendentLocation
    integer                                                   :: iIsolatedNode,iNode,iAxis,descendentIndex
    logical                                                   :: endOfBranch,descendentsFound
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
                ! Select the subset which have a subhalo as a descendent and are not the primary progenitor.
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
                         descendentIndex=Descendent_Node_Sort_Index(thisNode%descendentNode%nodeIndex)
                         if (descendentIndex > 0 .and. descendentIndex <= size(nodes)) then
                            descendentLocation=descendentLocations(descendentIndex)
                            descendentsFound=(nodes(descendentLocation)%descendentNode%nodeIndex == thisNode%descendentNode%nodeIndex)
                         else
                            descendentsFound=.false.
                         end if
                         if (descendentsFound) then
                            do while (nodes(descendentLocation)%descendentNode%nodeIndex == thisNode%descendentNode%nodeIndex .and. descendentIndex > 0)
                               if     (                                                                                     &
                                    &                    nodes(descendentLocation)%nodeIndex         /= thisNode%nodeIndex  &
                                    &  .and.             nodes(descendentLocation)%isolatedNodeIndex /= nodeIsUnreachable   &
                                    &  .and.  associated(nodes(descendentLocation)%descendentNode                         ) &
                                    &  .and.             nodes(descendentLocation)%nodeMass           > thisNode%nodeMass   &
                                    & ) then
                                  ! Subhalo-subhalo merger.
                                  endOfBranch =.true.
                                  exit
                               end if
                               do while (descendentIndex > 0)
                                  descendentIndex   =descendentIndex-1
                                  if (descendentIndex <= 0) exit
                                  descendentLocation=descendentLocations(descendentIndex)
                                  if (associated(nodes(descendentLocation)%descendentNode)) exit
                               end do
                               if (descendentIndex == 0) exit
                            end do
                         end if
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
                      !$omp critical(HDF5_Access)
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
                      !$omp end critical(HDF5_Access)
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

  subroutine Validate_Isolated_Halos(nodes)
    !% Ensure that nodes have valid primary progenitors.
    implicit none
    type(nodeData), intent(inout), dimension(:) :: nodes
    type(treeNode), pointer                     :: newNode
    integer                                     :: iNode,iIsolatedNode

    ! Search for cases where a node has no progenitors which do not descend into subhalos.
    do iNode=1,size(nodes)
       iIsolatedNode=0
       ! Select isolated nodes.
       if (nodes(iNode)%nodeIndex == nodes(iNode)%hostNode%nodeIndex) then
          iIsolatedNode=iIsolatedNode+1
          ! Select nodes with parents.
          if (associated(nodes(iNode)%node%parentNode)) then
             ! Select nodes with subhalo descendents which are also the primary progenitor of their parent.
             if (nodes(iNode)%descendentNode%isSubhalo.and.associated(nodes(iNode)%node%parentNode%childNode,nodes(iNode)%node)) then
                ! Insert a copy of the parent node as its own primary progenitor. This avoids current node being promoted into its
                ! parent even though it is intended to descend into a subhalo. The copy is shifted to a very slightly earlier
                ! time to avoid having two identical halos existing simultaneously (which can be problematic if outputting
                ! quantities which use the node index as a label in dataset names for example).
                allocate(newNode)
                call newNode%copy(nodes(iNode)%node%parentNode)
                newNode%siblingNode                    => nodes(iNode)%node
                newNode%parentNode                     => nodes(iNode)%node%parentNode
                newNode%childNode                      => null()
                nodes(iNode)%node%parentNode%childNode => newNode
                call Tree_Node_Time_Set(newNode,Tree_Node_Time(newNode)*(1.0d0-1.0d-6))
             end if
          end if
       end if
    end do
    return
  end subroutine Validate_Isolated_Halos

  subroutine Dump_Tree(nodes,highlightNodes)
    !% Dumps the tree structure to a file in a format suitable for processing with \href{http://www.graphviz.org/}{\sc dot}.
    implicit none
    type(nodeData),          intent(in), dimension(:)           :: nodes
    integer(kind=kind_int8), intent(in), dimension(:), optional :: highlightNodes
    integer                                                     :: iNode,fileUnit
    character(len=20)                                           :: color,style

    ! Open an output file and write the GraphViz opening.
    open(newunit=fileUnit,file='mergerTreeConstructReadTree.gv',status='unknown',form='formatted')
    write (fileUnit,*) 'digraph Tree {'

    ! Loop over all nodes.
    do iNode=1,size(nodes)
       ! Write each node, setting the node shape to a box for subhalos and a circle for halos. Node label consists of the node
       ! index plus the time, separated by a colon.
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
          write (fileUnit,'(a,i16.16,a,i16.16,a,f5.2,a,a,a,a,a)') '"',nodes(iNode)%nodeIndex,'" [shape=box   , label="',nodes(iNode)%nodeIndex,':',nodes(iNode)%nodeTime,'", color=',trim(color),', style=',trim(style),'];'
          ! If a host node is given, add a link to it as a red line.
          if (associated(nodes(iNode)%hostNode)) write (fileUnit,'(a,i16.16,a,i16.16,a)') '"',nodes(iNode)%nodeIndex,'" -> "',nodes(iNode)%hostNode%nodeIndex,'" [color=red];'
       else
          write (fileUnit,'(a,i16.16,a,i16.16,a,f5.2,a,a,a,a,a)') '"',nodes(iNode)%nodeIndex,'" [shape=circle, label="',nodes(iNode)%nodeIndex,':',nodes(iNode)%nodeTime,'", color=',trim(color),', style=',trim(style),'];'
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
