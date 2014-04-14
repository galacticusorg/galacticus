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

!% Contains a module which implements reading of merger trees from a file.

module Merger_Tree_Read
  !% Implements reading of merger trees from a file.
  use Galacticus_Nodes
  use ISO_Varying_String
  use Kind_Numbers
  use Merger_Tree_Read_Importers
  implicit none
  private
  public :: Merger_Tree_Read_Initialize

  ! The name of the file from which to read merger trees and its internal object.
  type            (varying_string                )                                                :: mergerTreeReadFileName                                                                       
  class           (mergerTreeImporterClass       )                                      , pointer :: defaultImporter                                                                              
  
  ! Index of the next merger tree to read.
  integer                                                                                         :: nextTreeToRead                           =0                                                  
  
  ! Index of the first tree to process.
  integer         (kind=kind_int8                )                                                :: mergerTreeReadBeginAt                                                                        
  
  ! Volume weight factor for trees.
  double precision                                                                                :: treeVolumeWeightCurrent                                                                      
  !$omp threadprivate(treeVolumeWeightCurrent)

  ! Flag indicating whether branch jumps are allowed.
  logical                                                                                         :: mergerTreeReadAllowBranchJumps                                                               
  
  ! Flag indicating whether subhalo promotions are allowed.
  logical                                                                                         :: mergerTreeReadAllowSubhaloPromotions                                                         
  
  ! Flags indicating whether or not to preset subhalo properties.
  logical                                                                                         :: mergerTreeReadPresetMergerTimes                                                              
  logical                                                                                         :: mergerTreeReadPresetMergerNodes                                                              
  logical                                                                                         :: mergerTreeReadPresetSubhaloMasses                                                            
  logical                                                                                         :: mergerTreeReadPresetPositions                                                                
  logical                                                                                         :: mergerTreeReadPresetScaleRadii                    , mergerTreeReadPresetScaleRadiiFailureIsFatal
  double precision                                                                                :: mergerTreeReadPresetScaleRadiiMinimumMass         , mergerTreeReadPresetScaleRadiiConcentrationMinimum, &
       &                                                                                             mergerTreeReadPresetScaleRadiiConcentrationMaximum
  logical                                                                                         :: mergerTreeReadPresetSpins                         , mergerTreeReadPresetSpins3D                                                         
  logical                                                                                         :: mergerTreeReadPresetOrbits                        , mergerTreeReadPresetOrbitsAssertAllSet            , & 
       &                                                                                             mergerTreeReadPresetOrbitsBoundOnly               , mergerTreeReadPresetOrbitsSetAll          
  logical                                                                                         :: mergerTreeReadPresetParticleCounts                , mergerTreeReadPresetVelocityMaxima                , &
       &                                                                                             mergerTreeReadPresetVelocityDispersions

  ! Option controlling fatality of missing host node condition.
  logical                                                                                         :: mergerTreeReadMissingHostsAreFatal                                                           
  
  ! Option controlling whether tree indices should always be set to the corresponding node index.
  logical                                                                                         :: mergerTreeReadTreeIndexToRootNodeIndex                                                       
  
  ! Labels used for node properties.
  integer                                                                    , parameter          :: nodeIsUnreachable                        =-1                                                 
  integer                                                                    , parameter          :: nodeIsReachable                          =0                                                  
  
  ! Internal list of output times and the relative tolerance used to "snap" nodes to output times.
  integer                                                                                         :: outputTimesCount                                                                             
  double precision                                                                                :: mergerTreeReadOutputTimeSnapTolerance                                                        
  double precision                                , allocatable, dimension(:)                     :: outputTimes                                                                                  
  
  ! Node used in root finding.
  class           (nodeComponentDarkMatterProfile)                                      , pointer :: activeDarkMatterProfileComponent                                                             
  class           (nodeComponentBasic            )                                      , pointer :: activeBasicComponent                                                                         
  type            (treeNode                      )                                      , pointer :: activeNode                                                                                   
  double precision                                                                                :: halfMassRadius                                                                               
  !$omp threadprivate(activeDarkMatterProfileComponent,activeBasicComponent,activeNode,halfMassRadius)
  ! Sorted node index list.
  integer         (kind=kind_int8                ), allocatable, dimension(:)                     :: descendentLocations                              , nodeLocations                             
  integer         (kind=kind_int8                ), allocatable, dimension(:)                     :: descendentIndicesSorted                          , nodeIndicesSorted                         
  
  ! Effective infinity for merging times.
  double precision                                                           , parameter          :: timeUntilMergingInfinite                 =1.0d30                                             
  
  ! Record of warnings issued.
  logical                                                                                         :: warningNestedHierarchyIssued             =.false.                                            
  
  ! Iterator object for iterating over progenitor nodes.
  type :: progenitorIterator
     integer(kind=kind_int8) :: progenitorIndex , progenitorLocation, targetIndex 
     logical                 :: progenitorsFound                                  
   contains
     !@ <objectMethods>
     !@   <object>progenitorIterator</object>
     !@   <objectMethod>
     !@     <method>descendentSet</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless class(nodeData)\textgreater} node\argin, \textcolor{red}{\textless class(nodeData)(:)\textgreater} nodes\argin</arguments>
     !@     <description>Set the target descendent node and initialize the iterator.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>next</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless class(nodeData)(:)\textgreater} nodes\argin</arguments>
     !@     <description>Move to the next progenitor. Returns true if the next progenitor exists, false otherwise.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>index</method>
     !@     <type>\textcolor{red}{\textless integer(kind=kind\_int8)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless class(nodeData)(:)\textgreater} nodes\argin</arguments>
     !@     <description>Return the index of the current progenitor.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>current</method>
     !@     <type>\textcolor{red}{\textless *class(nodeData)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless class(nodeData)(:)\textgreater} nodes\argin</arguments>
     !@     <description>Return a pointer to the current progenitor.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>exist</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if any progenitors exist, false otherwise.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: descendentSet=>progenitorIteratorDescendentSet 
     procedure :: next         =>progenitorIteratorNext          
     procedure :: index        =>progenitorIteratorIndex         
     procedure :: current      =>progenitorIteratorCurrent       
     procedure :: exist        =>progenitorIteratorExist         
  end type progenitorIterator

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
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Boolean
    use Memory_Management
    use Cosmology_Parameters
    implicit none
    type     (varying_string          ), intent(in   )          :: mergerTreeConstructMethod 
    procedure(Merger_Tree_Read_Do     ), intent(inout), pointer :: Merger_Tree_Construct     
    integer                                                     :: iOutput                   
    type     (varying_string          )                         :: message                   
    
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
       !@   <name>mergerTreeReadPresetScaleRadiiFailureIsFatal</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether failure to set a node scale radii should be regarded as a fatal error. (If not, a fallback method to set scale radius is used in such cases.)
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetScaleRadiiFailureIsFatal',mergerTreeReadPresetScaleRadiiFailureIsFatal,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetScaleRadiiConcentrationMinimum</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>3</defaultValue>
       !@   <description>
       !@     The lowest concentration ($c=r_{\rm vir}/r_{\rm s}$) allowed when setting scale radii, $r_{\rm s}$.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetScaleRadiiConcentrationMinimum',mergerTreeReadPresetScaleRadiiConcentrationMinimum,defaultValue=3.0d0)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetScaleRadiiConcentrationMaximum</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>60</defaultValue>
       !@   <description>
       !@     The largest concentration ($c=r_{\rm vir}/r_{\rm s}$) allowed when setting scale radii, $r_{\rm s}$.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetScaleRadiiConcentrationMaximum',mergerTreeReadPresetScaleRadiiConcentrationMaximum,defaultValue=60.0d0)
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
       !@   <name>mergerTreeReadPresetSpins3D</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether node 3-D spin vectors should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetSpins3D',mergerTreeReadPresetSpins3D,defaultValue=.false.)
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
       !@   <name>mergerTreeReadPresetParticleCounts</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>false</defaultValue>
       !@   <description>
       !@     Specifies whether node particle counts should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetParticleCounts',mergerTreeReadPresetParticleCounts,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetVelocityMaxima</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>false</defaultValue>
       !@   <description>
       !@     Specifies whether node rotation curve velocity maxima should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetVelocityMaxima',mergerTreeReadPresetVelocityMaxima,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadPresetVelocityDispersions</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>false</defaultValue>
       !@   <description>
       !@     Specifies whether node velocity dispersions should be preset when reading merger trees from a file.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadPresetVelocityDispersions',mergerTreeReadPresetVelocityDispersions,defaultValue=.false.)
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
       !@ <inputParameter>
       !@   <name>mergerTreeReadAllowBranchJumps</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether nodes are allowed to jump between branches.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadAllowBranchJumps',mergerTreeReadAllowBranchJumps,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>mergerTreeReadAllowSubhaloPromotions</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>true</defaultValue>
       !@   <description>
       !@     Specifies whether subhalos are permitted to be promoted to being isolated halos.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadAllowSubhaloPromotions',mergerTreeReadAllowSubhaloPromotions,defaultValue=.true.)

       ! Get array of output times.
       outputTimesCount=Galacticus_Output_Time_Count()
       call Alloc_Array(outputTimes,[outputTimesCount])
       do iOutput=1,outputTimesCount
          outputTimes(iOutput)=Galacticus_Output_Time(iOutput)
       end do

       ! Open the file.
       defaultImporter => mergerTreeImporter()
       call defaultImporter%open(mergerTreeReadFileName)

       ! Validate input parameters.
       if (mergerTreeReadPresetMergerNodes.and..not.mergerTreeReadPresetMergerTimes) then
          message="presetting of merger target nodes requires that merger times also be preset;"//char(10)
          message=message//" try setting [mergerTreeReadPresetMergerTimes]=true."//char(10)
          if (defaultimporter%treesHaveSubhalos() /= booleanTrue) then
             message=message//" Note: presetting merger target nodes and merger times is usually only a good idea when subhalo information is present in the merger trees"
             if (defaultimporter%treesHaveSubhalos() == booleanFalse) then
                message=message//" (which the current trees do not)"
             else
                message=message//" (subhalo presence in the current trees was not specified)"
             end if
          end if
          call Galacticus_Error_Report("Merger_Tree_Read_Initialize",message)
       end if

       ! Warn if subhalo promotions are allowed, but branch jumps are not.
       if (mergerTreeReadAllowSubhaloPromotions.and..not.mergerTreeReadAllowBranchJumps) then
          message='WARNING: allowing subhalo promotions while not allowed branch jumps can lead to deadlocking of trees.'//char(10)
          message=message//'Be sure that your trees have no promotion events for subhalos which survive beyond the end of their original branch'//char(10)
          message=message//'For example, without branch jumping, "a" in the following tree (in which "==>" indicates subhalo host) is stuck in "1", so it cannot evolve to become "b" causing the subhalo promotion "b" to "c" to be unreachable, resulting in a deadlock of the tree:'//char(10)//char(10)
          message=message//' ---   -----'//char(10)
          message=message//' |a|==>| 1 |'//char(10)
          message=message//' ---   -----'//char(10)
          message=message//'  |         '//char(10)
          message=message//' ---   -----'//char(10)
          message=message//' |b|==>| 2 |'//char(10)
          message=message//' ---   -----'//char(10)
          message=message//'  |      |  '//char(10)
          message=message//' ---   -----'//char(10)
          message=message//' |c|   | 3 |'//char(10)
          message=message//' ---   -----'//char(10)
          call Galacticus_Display_Message(message,verbosityWarn)
       end if

       ! Perform sanity checks if subhalos are not included.
       if (defaultimporter%treesHaveSubhalos() == booleanFalse) then
          if (mergerTreeReadPresetMergerTimes  ) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','cannot preset merger times as no subhalos are present; try setting [mergerTreeReadPresetMergerTimes]=false')
          if (mergerTreeReadPresetMergerNodes  ) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','cannot preset merger nodes as no subhalos are present; try setting [mergerTreeReadPresetMergerNodes]=false')
          if (mergerTreeReadPresetSubhaloMasses) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','cannot preset subhalo masses as no subhalos are present; try setting [mergerTreeReadPresetSubhaloMasses]=false')
       end if
       ! Determine if subhalo masses have been included in halo masses.   
       if (defaultImporter%treesAreSelfContained() == booleanFalse) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','only self-contained trees are supported')   
       ! Check that position information is present if required.
       if     (                                                                                   &
            &   (                                                                                 &
            &     mergerTreeReadPresetPositions                                                   &
            &    .or.                                                                             &
            &     mergerTreeReadPresetOrbits                                                      &
            &   )                                                                                 &
            &  .and.                                                                              &
            &   .not.defaultImporter%positionsAvailable(positions=.true.,velocities=.true.)       &
            & )                                                                                   &
            & call Galacticus_Error_Report(                                                       &
            &                              "Merger_Tree_Read_Initialize"                        , &
            &                              "presetting positions requires that both position and  &
            &                               velocity datasets be present in merger tree file"     &
            &                             )
       ! Check that half-mass radius information is present if required.
       if     (                                                                                                                                                    &
            &   mergerTreeReadPresetScaleRadii                                                                                                                     &
            &  .and.                                                                                                                                               &
            &   .not.defaultImporter%scaleRadiiAvailable()                                                                                                         &
            & )  call Galacticus_Error_Report(                                                                                                                     &
            &                                "Merger_Tree_Read_Initialize"                                                                                       , &
            &                                "presetting scale radii requires that at least one of halfMassRadius or scaleRadius datasets be present in merger"//  &
            &                                "tree file; try setting"//char(10)//"  [mergerTreeReadPresetScaleRadii]=false"                                        &
            &                               )
       ! Check that particle count information is present if required.
       if (                                                                                                                                   &
            &   mergerTreeReadPresetParticleCounts                                                                                            &
            &  .and.                                                                                                                          &
            &   .not.defaultImporter%particleCountAvailable()                                                                                 &
            & )  call Galacticus_Error_Report(                                                                                                &
            &                                "Merger_Tree_Read_Initialize"                                                                  , &
            &                                "presetting particle counts requires that the particleCount dataset be present in the merger"//  &
            &                                "tree file; try setting"//char(10)//"  [mergerTreeReadPresetParticleCount]=false"                &
            &                               )
       ! Check that velocity maximum information is present if required.
       if (                                                                                                                                     &
            &   mergerTreeReadPresetVelocityMaxima                                                                                              &
            &  .and.                                                                                                                            &
            &   .not.defaultImporter%velocityMaximumAvailable()                                                                                 &
            & )  call Galacticus_Error_Report(                                                                                                  &
            &                                "Merger_Tree_Read_Initialize"                                                                    , &
            &                                "presetting velocity maxima requires that the velocityMaximum dataset be present in the merger"//  &
            &                                "tree file; try setting"//char(10)//"  [mergerTreeReadPresetVelocityMaxima]=false"                 &
            &                               )
       ! Check that velocity dispersion information is present if required.
       if (                                                                                                                                             &
            &   mergerTreeReadPresetVelocityDispersions                                                                                                 &
            &  .and.                                                                                                                                    &
            &   .not.defaultImporter%velocityDispersionAvailable()                                                                                      &
            & )  call Galacticus_Error_Report(                                                                                                          &
            &                                "Merger_Tree_Read_Initialize"                                                                            , &
            &                                "presetting velocity dispersions requires that the velocityDispersion dataset be present in the merger"//  &
            &                                "tree file; try setting"//char(10)//"  [mergerTreeReadPresetVelocityDispersions]=false"                    &
            &                               )
       ! Check that angular momentum information is present if required.
       if     (                                                                                                        &
            &   mergerTreeReadPresetSpins                                                                              &
            &  .and.                                                                                                   &
            &   .not.                                                                                                  &
            &        (                                                                                                 &
            &          defaultImporter%          spinAvailable()                                                       &
            &         .or.                                                                                             &
            &          defaultImporter%angularMomentaAvailable()                                                       &
            &        )                                                                                                 &
            & )                                                                                                        &
            & call Galacticus_Error_Report(                                                                            &
            &                              "Merger_Tree_Read_Initialize"                                            ,  &
            &                              "presetting spins requires that the spin or angularMomentum datasets        &
            &                               be present in merger tree file; try setting"                            // &
            &                              char(10)                                                                 // &
            &                              " [mergerTreeReadPresetSpins]=false"                                        &
            &                             )
       if     (                                                                                                        &
            &   mergerTreeReadPresetSpins3D                                                                            &
            &  .and.                                                                                                   &
            &   .not.                                                                                                  &
            &        (                                                                                                 &
            &          defaultImporter%          spin3DAvailable()                                                     &
            &         .or.                                                                                             &
            &          defaultImporter%angularMomenta3DAvailable()                                                     &
            &        )                                                                                                 &
            & )                                                                                                        &
            & call Galacticus_Error_Report(                                                                            &
            &                              "Merger_Tree_Read_Initialize"                                            ,  &
            &                              "presetting spin vectors requires that the spin or angularMomentum vector   &
            &                               datasets be present in merger tree file; try setting"                   // &
            &                              char(10)                                                                 // &
            &                              " [mergerTreeReadPresetSpins3D]=false"                                      &
            &                             )
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
    use Arrays_Search
    use Array_Utilities
    use Numerical_Comparison
use omp_lib
    implicit none
    type            (mergerTree                    )                             , intent(inout), target :: thisTree                               
    logical                                                                      , intent(in   )         :: skipTree                               
    double precision                                , allocatable, dimension(:  )                        :: historyMass        , historyTime       
    double precision                                , allocatable, dimension(:,:)                        :: position           , velocity          
    class           (nodeData                      ), allocatable, dimension(:  )               , target :: nodes                                  
    type            (treeNodeList                  ), allocatable, dimension(:  )                        :: thisNodeList                           
    logical                                         , allocatable, dimension(:  )                        :: childIsSubhalo                         
    integer                                                                                              :: iOutput            , isolatedNodeCount 
    integer         (kind=kind_int8                )                                                     :: historyCountMaximum, iNode             
    logical                                                                                              :: haveTree                               
    type            (varying_string                )                                                     :: message                                
    
    !$omp critical(mergerTreeReadTree)
    ! Increment the tree to read index.
    nextTreeToRead=nextTreeToRead+1
    ! Keep incrementing the tree index until we find the first tree to process (if we haven't done so already). Also skip trees
    ! that contain 1 or fewer nodes and these are unprocessable.
    if (nextTreeToRead <= defaultImporter%treeCount()) then
       do while (                                                                      &
            &     (                                                                    &
            &       mergerTreeReadBeginAt                     > 0                      &
            &      .and.                                                               &
            &       defaultImporter%treeIndex(nextTreeToRead) /= mergerTreeReadBeginAt &
            &     )                                                                    &
            &    .or.                                                                  &
            &       defaultImporter%nodeCount(nextTreeToRead) <= 1                     &
            &   )
          nextTreeToRead=nextTreeToRead+1
          ! If the end of the list has been reached, exit.
          if (nextTreeToRead > defaultImporter%treeCount()) exit
       end do
    end if

    ! Flag that we've now found the first merger tree to process.
    mergerTreeReadBeginAt=-1

    if (nextTreeToRead > defaultImporter%treeCount())  then
       ! All trees have been read.
       call defaultImporter%close()
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
          thisTree%baseNode => treeNode(hostTree=thisTree)
       else
          ! Set tree properties.
          ! treeIndex
          thisTree%index=defaultImporter%treeIndex(nextTreeToRead)
          ! volumeWeight
          treeVolumeWeightCurrent=defaultImporter%treeWeight(nextTreeToRead)
          thisTree%volumeWeight=treeVolumeWeightCurrent
          ! Initialize no events.
          thisTree%event => null()
          ! Read data from the file.
          call defaultImporter%import(                                                                                                                      &
               &                      nextTreeToRead                                                                                                      , &
               &                      nodes                                                                                                               , &
               &                      requireScaleRadii         = mergerTreeReadPresetScaleRadii                                                          , &
               &                      requireParticleCounts     = mergerTreeReadPresetParticleCounts                                                      , &
               &                      requireVelocityMaxima     = mergerTreeReadPresetVelocityMaxima                                                      , &
               &                      requireVelocityDispersions= mergerTreeReadPresetVelocityDispersions                                                 , &
               &                      requireAngularMomenta     =(mergerTreeReadPresetSpins              .and.defaultImporter%angularMomentaAvailable  ()), &
               &                      requireAngularMomenta3D   =(mergerTreeReadPresetSpins3D            .and.defaultImporter%angularMomenta3DAvailable()), &
               &                      requireSpin               =(mergerTreeReadPresetSpins              .and.defaultImporter%spinAvailable            ()), &
               &                      requireSpin3D             =(mergerTreeReadPresetSpins3D            .and.defaultImporter%spin3DAvailable          ()), &
               &                      requirePositions          =(mergerTreeReadPresetPositions          .or. mergerTreeReadPresetOrbits                 )  &
               &                     )

          ! Snap node times to output times if a tolerance has been specified.
          if (mergerTreeReadOutputTimeSnapTolerance > 0.0d0) then
             ! Loop over all nodes.
             do iNode=1,size(nodes)
                ! Find closest output time to the node time.
                iOutput=Search_Array_For_Closest(outputTimes,nodes(iNode)%nodeTime)
                ! Test if this time is sufficiently close that we should snap the node time to it.
                if (Values_Agree(nodes(iNode)%nodeTime,outputTimes(iOutput),relTol=mergerTreeReadOutputTimeSnapTolerance)) &
                     & nodes(iNode)%nodeTime=outputTimes(iOutput)
             end do
          end if
          
          ! Sort node indices.
          call Create_Node_Indices(nodes)

          ! Identify subhalos.
          nodes%isSubhalo=nodes%nodeIndex /= nodes%hostIndex

          ! Build pointers to descendent nodes.
          call Build_Descendent_Pointers(nodes)

          ! Find cases where something that was a subhalo stops being a subhalo and prevent them if necessary.
          call Enforce_Subhalo_Status(nodes)

          ! If necessary, add masses of subhalos to host halos.
          if (.not.defaultImporter%massesIncludeSubhalos()) then
             do iNode=1,size(nodes)
                if (nodes(iNode)%host%nodeIndex /= nodes(iNode)%nodeIndex) nodes(iNode)%host%nodeMass= &
                     &                                                     nodes(iNode)%host%nodeMass+nodes(iNode)%nodeMass
             end do
          end if

          ! Associate parent pointers with the descendent host.
          call Build_Parent_Pointers(nodes)

          ! Create an array of standard nodes.
          call Create_Node_Array(thisTree,nodes,thisNodeList,isolatedNodeCount,childIsSubhalo)

          ! Assign parent pointers and properties.
          call Build_Isolated_Parent_Pointers(thisTree,nodes,thisNodeList)

          ! Now build child and sibling links.
          call Build_Child_and_Sibling_Links(nodes,thisNodeList,childIsSubhalo)

          ! Check that all required properties exist.
          if (mergerTreeReadPresetPositions.or.mergerTreeReadPresetOrbits) then
             ! Position and velocity methods are required.
             if     (                                                                                                                                                                   &
                  &  .not.(                                                                                                                                                             &
                  &         defaultPositionComponent%positionIsSettable()                                                                                                               &
                  &        .and.                                                                                                                                                        &
                  &         defaultPositionComponent%velocityIsSettable()                                                                                                               &
                  &       )                                                                                                                                                             &
                  & )                                                                                                                                                                   &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting positions or orbits requires a component that supports position and velocity setting (e.g. set [treeNodeMethodPosition]=preset);'              // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'darkMatterProfile'                                                                                                               , &
                  &                                  defaultPositionComponent        %     positionAttributeMatch(requireSettable=.true.)                                               &
                  &                                 .intersection.                                                                                                                      &
                  &                                  defaultPositionComponent        %     velocityAttributeMatch(requireSettable=.true.)                                               &
                  &                                )                                                                                                                                 // &
                  &       char(10)                                                                                                                                                   // &
                  &       'alternatively setting [mergerTreeReadPresetPositions]=false and [mergerTreeReadPresetOrbits]=false will remove the need to store positions and velocities'   &
                  & )
          end if
          if (mergerTreeReadPresetMergerTimes) then
             ! Time of merging property is required.
             if (.not.defaultSatelliteComponent%timeOfMergingIsSettable      ())                                                                                                        &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting merging times requires a component that supports setting of merging times.'                                                                    // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'satellite'                                                                                                                       , &
                  &                                  defaultSatelliteComponent       %timeOfMergingAttributeMatch(requireSettable=.true.)                                               &
                  &                                 )                                                                                                                                   &
                  &      )
          end if
          if (mergerTreeReadPresetScaleRadii) then
             ! Scale radius property is required.
             if (.not.defaultDarkMatterProfileComponent%scaleIsSettable      ())                                                                                                        &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting scale radii requires a component that supports setting of scale radii.'                                                                        // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'darkMatterProfile'                                                                                                               , &
                  &                                 defaultDarkMatterProfileComponent%        scaleAttributeMatch(requireSettable=.true.)                                               &
                  &                                )                                                                                                                                    &
                  &      )
          end if
          if (mergerTreeReadPresetParticleCounts) then
             ! Particle count property is required.
             if (.not.defaultNBodyComponent%particleCountIsSettable          ())                                                                                                        &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting particle counts requires an nBody component that supports setting of particle count.'                                                          // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'nBody'                                                                                                                           , &
                  &                                 defaultNBodyComponent            %particleCountAttributeMatch(requireSettable=.true.)                                               &
                  &                                )                                                                                                                                    &
                  &      )
          end if
          if (mergerTreeReadPresetVelocityMaxima) then
             ! Velocity maximum property is required.
             if (.not.defaultNBodyComponent%velocityMaximumIsSettable        ())                                                                                                        &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting velocity maxima requires an nBody component that supports setting of velocity maxima.'                                                         // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'nBody'                                                                                                                           , &
                  &                                 defaultNBodyComponent            %velocityMaximumAttributeMatch(requireSettable=.true.)                                             &
                  &                                )                                                                                                                                    &
                  &      )
          end if
          if (mergerTreeReadPresetVelocityDispersions) then
             ! Velocity dispersion property is required.
             if (.not.defaultNBodyComponent%velocityDispersionIsSettable     ())                                                                                                        &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting velocity dispersions requires an nBody component that supports setting of velocity dispersions.'                                               // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'nBody'                                                                                                                           , &
                  &                                 defaultNBodyComponent            %velocityDispersionAttributeMatch(requireSettable=.true.)                                          &
                  &                                )                                                                                                                                    &
                  &      )
          end if
          if (mergerTreeReadPresetSpins      ) then
             ! Spin property is required.
             if (.not.defaultSpinComponent             %spinIsSettable       ())                                                                                                        &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting spins requires a component that supports setting of spins.'                                                                                    // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'spin'                                                                                                                            , &
                  &                                 defaultSpinComponent             %         spinAttributeMatch(requireSettable=.true.)                                               &
                  &                                )                                                                                                                                    &
                  &      )
          end if
          if (mergerTreeReadPresetSpins3D    ) then
             ! Spin property is required.
             if (.not.defaultSpinComponent             %spinVectorIsSettable ())                                                                                                        &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting spin vectors requires a component that supports setting of spin vectors.'                                                                      // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'spinVector'                                                                                                                      , &
                  &                                 defaultSpinComponent             %     spinVectorAttributeMatch(requireSettable=.true.)                                             &
                  &                                )                                                                                                                                    &
                  &      )
          end if
          if (mergerTreeReadPresetOrbits     ) then
             ! Orbit property is required.
             if (.not.defaultSatelliteComponent        %virialOrbitIsSettable())                                                                                                        &
                  & call Galacticus_Error_Report                                                                                                                                        &
                  &      (                                                                                                                                                              &
                  &       'Merger_Tree_Read_Do'                                                                                                                                       , &
                  &       'presetting orbits requires a component that supports setting of orbits (e.g. [treeNodeMethodSatelliteOrbit]=preset);'                                     // &
                  &       Galacticus_Component_List(                                                                                                                                    &
                  &                                 'satellite'                                                                                                                       , &
                  &                                 defaultSatelliteComponent        %  virialOrbitAttributeMatch(requireSettable=.true.)                                               &
                  &                                )                                                                                                                                 // &
                  &       char(10)                                                                                                                                                   // &
                  &       'Alternatively, set [mergerTreeReadPresetOrbits]=false to prevent attempts to set orbits)')
          end if

          ! Assign scale radii.
          if (mergerTreeReadPresetScaleRadii                               ) call Assign_Scale_Radii         (nodes,thisNodeList)

          ! Assign particle counts.
          if (mergerTreeReadPresetParticleCounts                           ) call Assign_Particle_Counts     (nodes,thisNodeList)

          ! Assign velocity maxima.
          if (mergerTreeReadPresetVelocityMaxima                           ) call Assign_Velocity_Maxima     (nodes,thisNodeList)

          ! Assign velocity dispersions.
          if (mergerTreeReadPresetVelocityDispersions                      ) call Assign_Velocity_Dispersions(nodes,thisNodeList)

          ! Assign spin parameters.
          if (mergerTreeReadPresetSpins     .or.mergerTreeReadPresetSpins3D) call Assign_Spin_Parameters     (nodes,thisNodeList)

          ! Assign isolated node indices to subhalos.
          call Assign_Isolated_Node_Indices(nodes)

          ! Ensure that isolated nodes with progenitors that descend into subhalos have valid primary progenitors.
          call Validate_Isolated_Halos(nodes)

          ! Scan subhalos to determine when and how they merge.
          call Scan_For_Mergers(nodes,thisNodeList,historyCountMaximum)

          ! Search for any nodes which were flagged as merging with another node and assign appropriate pointers.
          call Assign_Mergers(nodes,thisNodeList)

          ! Find cases where something that was a subhalo stops being a subhalo and add events to handle.
          call Scan_for_Subhalo_Promotions(nodes,thisNodeList)

          ! Search for subhalos which move between branches/trees.
          call Scan_for_Branch_Jumps(nodes,thisNodeList)

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

          ! Assign new uniqueIDs to any cloned nodes inserted into the trees.
          call Assign_UniqueIDs_To_Clones(thisNodeList)

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
    !$omp end critical(mergerTreeReadTree)
    return
  end subroutine Merger_Tree_Read_Do

  subroutine Create_Node_Indices(nodes)
    !% Create a sorted list of node indices with an index into the original array.
    use Galacticus_Display
    use String_Handling
    use Memory_Management
    use Sort
    implicit none
    class  (nodeData      ), dimension(:), intent(inout) :: nodes   
    integer                                              :: iNode   
    type   (varying_string)                              :: message 
    
    ! Build a sorted list of node indices with an index into the original arrays.
    call Alloc_Array(nodeLocations          ,shape(nodes))
    call Alloc_Array(nodeIndicesSorted      ,shape(nodes))
    call Alloc_Array(descendentLocations    ,shape(nodes))
    call Alloc_Array(descendentIndicesSorted,shape(nodes))
    nodeLocations      =Sort_Index_Do(nodes%nodeIndex      )
    descendentLocations=Sort_Index_Do(nodes%descendentIndex)
    forall (iNode=1:size(nodes))
       nodeIndicesSorted      (iNode)=nodes(nodeLocations      (iNode))%nodeIndex
       descendentIndicesSorted(iNode)=nodes(descendentLocations(iNode))%descendentIndex
    end forall
    if (Galacticus_Verbosity_Level() >= verbosityWarn) then
       do iNode=2,size(nodes)
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
    integer(kind=kind_int8)                :: Node_Location 
    integer(kind=kind_int8), intent(in   ) :: nodeIndex     
    integer(kind=kind_int8)                :: iNode         
    
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
    integer(kind=kind_int8), intent(in   ) :: descendentIndex 
    
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

  subroutine Build_Descendent_Pointers(nodes)
    !% Builds pointers from each node to its descendent node.
    use String_Handling
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    class  (nodeData      ), dimension(:), intent(inout), target :: nodes                 
    integer(kind=kind_int8)                                      :: iNode  , nodeLocation 
    type   (varying_string)                                      :: message               
    
    do iNode=1,size(nodes)
       if (nodes(iNode)%descendentIndex >= 0) then
          nodeLocation=Node_Location(nodes(iNode)%descendentIndex)
          if (nodes(nodeLocation)%nodeIndex /= nodes(iNode)%descendentIndex) then
             message='failed to find descendent node: '
             message=message//nodes(iNode)%descendentIndex//' of '//nodes(iNode)%nodeIndex
             call Galacticus_Error_Report('Build_Descendent_Pointers',message)
          end if
          nodes(iNode)%descendent => nodes(nodeLocation)
       else
          nodes(iNode)%descendent => null()
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
                nodes(iNode)%host      => nodes(iNode)
                nodes(iNode)%isSubhalo =  .false.
             end if
          else
             nodes(iNode)%host    => nodes(nodeLocation)
          end if
       else
          call Galacticus_Error_Report('Build_Descendent_Pointers','negative values are not allowed for hostIndex - if node is self-hosting [i.e. not a subhalo] set hostIndex=nodeIndex')
       end if
    end do
    return
  end subroutine Build_Descendent_Pointers

  subroutine Enforce_Subhalo_Status(nodes)
    !% Ensure that any node which was once a subhalo remains a subhalo.
    use Galacticus_Error
    use String_Handling
    implicit none
    class  (nodeData          ), dimension(:), intent(inout), target :: nodes                                    
    class  (nodeData          ), pointer                             :: descendentNode, progenitorNode           
    integer(kind=kind_int8    )                                      :: iNode                                    
    logical                                                          :: failed        , isolatedProgenitorExists 
    type   (varying_string    )                                      :: message                                  
    type   (progenitorIterator)                                      :: progenitors                              
    
    ! Return immediately if subhalo promotions are allowed.
    if (mergerTreeReadAllowSubhaloPromotions) return
    ! Subhalo promotions are not allowed, so enforce subhalo status.
    do iNode=1,size(nodes)
       if (nodes(iNode)%isSubhalo) then
          descendentNode => nodes(iNode)%descendent
          do while (associated(descendentNode))
             ! Is this node isolated?
             if (.not.descendentNode%isSubhalo) then
                ! Check if there is any isolated node which descends into this node.
                isolatedProgenitorExists=.false.
                call progenitors%descendentSet(descendentNode,nodes)
                do while (progenitors%next(nodes) .and. .not.isolatedProgenitorExists)
                   progenitorNode => progenitors%current(nodes)
                   isolatedProgenitorExists=(progenitorNode%nodeIndex == progenitorNode%hostIndex)
                end do
                if (.not.progenitors%exist() .or. .not.isolatedProgenitorExists) then
                   ! Node is isolated, has no isolated node that descends into it. Therefore, our current node is not allowed to
                   ! be a subhalo.
                   nodes(iNode)%isSubhalo=.false.
                   nodes(iNode)%host => nodes(iNode)
                   nodes(iNode)%hostIndex=nodes(iNode)%nodeIndex
                end if
             end if
             descendentNode => descendentNode%descendent
          end do
       end if
    end do
    ! Check that subhalo enforcement was successful.
    failed=.false.
    do iNode=1,size(nodes)
       ! Find nodes which have no isolated node descending into them.       
       isolatedProgenitorExists=.false.
       call progenitors%descendentSet(nodes(iNode),nodes)
       do while (progenitors%next(nodes) .and. .not.isolatedProgenitorExists)
          progenitorNode => progenitors%current(nodes)
          isolatedProgenitorExists=(progenitorNode%nodeIndex == progenitorNode%hostIndex)
       end do
       if (progenitors%exist() .and. .not.isolatedProgenitorExists) then
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

  subroutine Scan_for_Subhalo_Promotions(nodes,nodeList)
    !% Scan for cases where a subhalo stops being a subhalo and so must be promoted.
    use Node_Subhalo_Promotions
    implicit none
    class  (nodeData          )         , dimension(:), intent(inout), target :: nodes                                    
    type   (treeNodeList      )         , dimension(:), intent(inout)         :: nodeList                                 
    class  (nodeData          ), pointer                                      :: descendentNode          , progenitorNode 
    type   (nodeEvent         ), pointer                                      :: newEvent                , pairEvent      
    type   (treeNode          ), pointer                                      :: promotionNode           , thisNode       
    integer(kind=kind_int8    )                                               :: iNode                                    
    logical                                                                   :: isolatedProgenitorExists                 
    type   (progenitorIterator)                                               :: progenitors                              
    
    ! Return immediately if subhalo promotion is not allowed.
    if (.not.mergerTreeReadAllowSubhaloPromotions) return
    ! Find subhalos to be promoted.
    do iNode=1,size(nodes)
       if (nodes(iNode)%isSubhalo.and.associated(nodes(iNode)%descendent)) then
          descendentNode => nodes(iNode)%descendent
          ! Is this node isolated?
          if (.not.descendentNode%isSubhalo) then
             ! Check if there is any isolated node which descends into this node.
             isolatedProgenitorExists=.false.
             call progenitors%descendentSet(descendentNode,nodes)
             do while (progenitors%next(nodes) .and. .not.isolatedProgenitorExists)
                progenitorNode => progenitors%current(nodes)
                isolatedProgenitorExists=(progenitorNode%nodeIndex == progenitorNode%hostIndex)
             end do
             if (.not.progenitors%exist() .or. .not.isolatedProgenitorExists) then
                ! Node is isolated, has no isolated node that descends into it. Therefore, our subhalo must be promoted to
                ! become an isolated halo again.
                thisNode       => nodeList(nodes(inode)  %isolatedNodeIndex)%node
                promotionNode  => nodeList(descendentNode%isolatedNodeIndex)%node
                newEvent       => thisNode     %createEvent()
                newEvent %time =  descendentNode%nodeTime
                newEvent %node => promotionNode
                newEvent %task => Node_Subhalo_Promotion
                pairEvent      => promotionNode%createEvent()
                pairEvent%time =  descendentNode%nodeTime
                pairEvent%node => thisNode
                pairEvent%task => null()
                pairEvent%ID   =  newEvent%ID
             end if
          end if
       end if
    end do
    return
  end subroutine Scan_for_Subhalo_Promotions

  subroutine Build_Parent_Pointers(nodes)
    !% Build pointers to node parents.
    use Galacticus_Error
    use String_Handling
    implicit none
    class  (nodeData      ), dimension(:), intent(inout), target :: nodes      
    class  (nodeData      ), pointer                             :: parentNode 
    integer(kind=kind_int8)                                      :: iNode      
    type   (varying_string)                                      :: message    
    
    do iNode=1,size(nodes)
       if (associated(nodes(iNode)%descendent)) then
          if (nodes(iNode)%nodeIndex == nodes(iNode)%host%nodeIndex) then
             ! Find an isolated parent node, by repeatedly jumping from host to host.
             parentNode => nodes(iNode)%descendent%host
             do while (parentNode%isSubhalo)
                if (associated(parentNode,parentNode%host)) then
                   message='node ['
                   message=message//parentNode%nodeIndex//'] flagged as subhalo is self-hosting - exiting to avoid infinite loop'
                   call Galacticus_Error_Report('Build_Parent_Pointers',message)
                end if
                parentNode => parentNode%host
             end do
             nodes(iNode)%parent => parentNode
          else
             nodes(iNode)%parent => null()
          end if
       else
          nodes(iNode)%parent => null()
       end if
    end do
    ! Check for self-parents.
    do iNode=1,size(nodes)
       if (associated(nodes(iNode)%parent)) then
          if (nodes(iNode)%nodeIndex == nodes(iNode)%parent%nodeIndex) then
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
    implicit none
    type   (mergerTree        )                           , intent(inout) :: thisTree                                 
    class  (nodeData          )             , dimension(:), intent(inout) :: nodes                                    
    type   (treeNodeList      ), allocatable, dimension(:), intent(inout) :: nodeList                                 
    logical                    , allocatable, dimension(:), intent(inout) :: childIsSubhalo                           
    integer                                               , intent(  out) :: isolatedNodeCount                        
    integer(kind=kind_int8    )                                           :: iNode                                    
    integer                                                               :: iIsolatedNode    , initialSatelliteCount 
    logical                                                               :: createNode                               
    type   (progenitorIterator)                                           :: progenitors                              
    
    ! Determine how many nodes are isolated (i.e. not subhalos).
    isolatedNodeCount=count(.not.nodes%isSubhalo)

    ! Scan here for nodes that are subhalos and have no progenitor. These objects must be
    ! created as satellites within the tree.
    initialSatelliteCount=0
    do iNode=1,size(nodes)
       if (nodes(iNode)%isSubhalo) then
          call progenitors%descendentSet(nodes(iNode),nodes)
          if (.not.progenitors%exist()) initialSatelliteCount=initialSatelliteCount+1
       end if
    end do

    ! Allocate nodes.
    allocate(nodeList(isolatedNodeCount+initialSatelliteCount))
    call Memory_Usage_Record(sizeof(nodeList))
    call Alloc_Array(childIsSubhalo,[isolatedNodeCount+initialSatelliteCount])

    ! Create the nodes.
    iIsolatedNode          =0
    nodes%isolatedNodeIndex=nodeIsUnreachable
    do iNode=1,size(nodes)
       createNode=.false.
       if (nodes(iNode)%nodeIndex == nodes(iNode)%host%nodeIndex) then
          createNode=.true.
       else if (nodes(iNode)%isSubhalo) then
          call progenitors%descendentSet(nodes(iNode),nodes)
          if (.not.progenitors%exist()) createNode=.true.
       end if
       if (createNode) then
          iIsolatedNode=iIsolatedNode+1
          ! Store a record of where this node goes in the isolated node list.
          nodes(iNode)%isolatedNodeIndex=iIsolatedNode
          nodeList(iIsolatedNode)%node => treeNode(hostTree=thisTree)
          call nodeList(iIsolatedNode)%node%indexSet(nodes(iNode)%nodeIndex)
          nodes(iNode)%node => nodeList(iIsolatedNode)%node
       end if
    end do
    return
  end subroutine Create_Node_Array

  subroutine Build_Isolated_Parent_Pointers(thisTree,nodes,nodeList)
    !% Create parent pointer links between isolated nodes and assign times and masses to those nodes.
    use String_Handling
    use Galacticus_Error
    implicit none
    type     (mergerTree        )                       , intent(inout) , target::                    thisTree 
    class    (nodeData          )         , dimension(:), intent(inout) ::      nodes                          
    type     (treeNodeList      )         , dimension(:), intent(inout) ::      nodeList                       
    class    (nodeComponentBasic), pointer                              ::      nodeBasicComponent             
    type     (mergerTree        ), pointer                              ::      currentTree                    
    class    (nodeData          ), pointer                              ::      parentNode                     
    integer                                                             ::      iNode                          
    integer  (kind=kind_int8    )                                       ::      iIsolatedNode                  
    type     (varying_string    )                                       ::      message                        
    character(len=12            )                                       ::      label                          
    logical                                                             ::      assignLastIsolatedTime         
    
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node (or an initial satellite).
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Check for an isolated node.
          if (nodes(iNode)%nodeIndex == nodes(iNode)%host%nodeIndex) then
             assignLastIsolatedTime=.false.
             if (associated(nodes(iNode)%parent)) then
                nodeList(iIsolatedNode)%node%parent  => nodes(iNode)%parent%node
             else
                nodeList(iIsolatedNode)%node%parent  => null()
                ! Find a tree to attach this base node to. Begin with the original tree passed to us.
                currentTree => thisTree
                ! Check if its baseNode is already assigned.
                do while (associated(currentTree%baseNode))
                   ! While it is, create the next tree (unless it already exists), then step to it.
                   if (.not.associated(currentTree%nextTree)) allocate(currentTree%nextTree)
                   currentTree => currentTree%nextTree
                end do
                ! Assign this node as the base node of the current tree.
                currentTree   %baseNode    => nodeList(iIsolatedNode)%node
                if (mergerTreeReadTreeIndexToRootNodeIndex) then
                   currentTree%index       =  nodes   (iNode        )%nodeIndex
                else
                   currentTree%index       =  thisTree               %index
                end if
                currentTree   %volumeWeight=  treeVolumeWeightCurrent
                currentTree   %initialized =  .false.
                currentTree   %event       => null()
             end if
          else
             ! Node is not isolated, so must be an initial satellite.
             assignLastIsolatedTime=.true.
             if (associated(nodes(iNode)%host)) then
                parentNode => nodes(iNode)%host
                do while (parentNode%isSubhalo)
                   parentNode => parentNode%host
                end do
                nodeList(iIsolatedNode)%node%parent => parentNode%node
             else
                call Galacticus_Error_Report('Build_Isolated_Parent_Pointers','initial satellite has no parent defined')
             end if
          end if
          ! Assign mass and time. For the case of satellites we also assign the time at
          ! which the satellite was last isolated. Since we do not know this, we simply
          ! set it equal to the current time (which is, obviously, an upper limit).
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
          nodeBasicComponent => nodeList(iIsolatedNode)%node%basic(autoCreate=.true.)
          call        nodeBasicComponent%massSet            (nodes(iNode)%nodeMass)
          call        nodeBasicComponent%timeSet            (nodes(iNode)%nodeTime)
          if (assignLastIsolatedTime) &
               & call nodeBasicComponent%timeLastIsolatedSet(nodes(iNode)%nodeTime)
       end if
    end do
    return
  end subroutine Build_Isolated_Parent_Pointers

  subroutine Build_Child_and_Sibling_Links(nodes,nodeList,childIsSubhalo)
    !% Build child and sibling links between nodes.
    use Memory_Management
    implicit none
    class  (nodeData          )             , dimension(:), intent(inout) :: nodes                                     
    type   (treeNodeList      )             , dimension(:), intent(inout) :: nodeList                                  
    logical                    , allocatable, dimension(:), intent(inout) :: childIsSubhalo                            
    class  (nodeComponentBasic), pointer                                  :: nodeBasicComponent, primaryBasicComponent 
    integer                                                               :: iNode                                     
    integer(kind=kind_int8    )                                           :: iIsolatedNode                             
    logical                                                               :: descendsToSubhalo                         
    
    childIsSubhalo=.false.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Check for an isolated node.
          if (nodes(iNode)%nodeIndex == nodes(iNode)%host%nodeIndex) then
             ! Check if the node has a parent.
             if (associated(nodeList(iIsolatedNode)%node%parent)) then
                ! Determine if this node definitely descends to a subhalo - in which case it can never be the primary progenitor.
                descendsToSubhalo=nodes(iNode)%descendentIndex /= nodeList(iIsolatedNode)%node%parent%index()
                ! It does, so set the child pointer of the parent appropriately.
                if (associated(nodeList(iIsolatedNode)%node%parent%firstChild)) then
                   ! A child is already associated. Check if current node does not descend to a subhalo and is more massive.
                   nodeBasicComponent    => nodeList(iIsolatedNode)%node%basic()
                   primaryBasicComponent => nodeList(iIsolatedNode)%node%parent%firstChild%basic()
                   if (.not.descendsToSubhalo                                             &
                        & .and. (                                                         &
                        &        childIsSubhalo(nodes(iNode)%parent%isolatedNodeIndex)    &
                        &         .or.                                                    &
                        &        nodeBasicComponent%mass() > primaryBasicComponent%mass() &
                        &       )                                                         &
                        & ) then
                      ! It is, so make this the main progenitor.
                      nodeList(iIsolatedNode)%node%sibling           => nodeList(iIsolatedNode)%node%parent%firstChild
                      nodeList(iIsolatedNode)%node%parent%firstChild => nodeList(iIsolatedNode)%node
                      ! Record that the main child is now not a subhalo.
                      childIsSubhalo(nodes(iNode)%parent%isolatedNodeIndex)=.false.
                   else
                      ! It is not, so add after the main child.
                      nodeList(iIsolatedNode)%node%sibling                   => nodeList(iIsolatedNode)%node%parent%firstChild%sibling
                      nodeList(iIsolatedNode)%node%parent%firstChild%sibling => nodeList(iIsolatedNode)%node
                   end if
                else
                   ! No child is currently associated. Simply point to the current node.
                   nodeList(iIsolatedNode)%node%parent%firstChild => nodeList(iIsolatedNode)%node
                   ! Record whether or not this child is a known subhalo or not.
                   childIsSubhalo(nodes(iNode)%parent%isolatedNodeIndex)=descendsToSubhalo
                end if
             end if
          else
             ! Node must be an initial satellite.
             if (associated(nodeList(iIsolatedNode)%node%parent%firstSatellite)) then
                ! The parent halo already has some satellites. Add this one to the list.
                nodeList(iIsolatedNode)%node%sibling               => nodeList(iIsolatedNode)%node%parent%firstSatellite
                nodeList(iIsolatedNode)%node%parent%firstSatellite => nodeList(iIsolatedNode)%node
             else
                ! The parent halo does not yet have any satellites. Simply add this one as the first.
                nodeList(iIsolatedNode)%node%sibling               => null()
                nodeList(iIsolatedNode)%node%parent%firstSatellite => nodeList(iIsolatedNode)%node
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
    use Dark_Matter_Profiles_Concentration
    use Galacticus_Display
    use Galacticus_Error
    use Input_Parameters
    implicit none
    class           (nodeData                           )                 , dimension(:), intent(inout) :: nodes                                                                 
    type            (treeNodeList                       )                 , dimension(:), intent(inout) :: nodeList                                                              
    double precision                                     , parameter                                    :: scaleRadiusMaximumAllowed                =100.0d0, toleranceAbsolute  =1.0d-9, & 
         &                                                                                                 toleranceRelative                        =1.0d-9                                 
    logical                                                         , save                              :: excessiveScaleRadiiReported              =.false.                                
    class           (nodeComponentBasic                 ), pointer                                      :: thisBasicComponent                                                    
    class           (nodeComponentDarkMatterProfile     ), pointer                                      :: thisDarkMatterProfileComponent                                        
    integer                                                                                             :: iNode                                            , status                     , &
         &                                                                                                 messageVerbosity
    integer         (kind=kind_int8                     )                                               :: iIsolatedNode                                                         
    double precision                                                                                    :: radiusScale                                                           
    logical                                                                                             :: excessiveHalfMassRadii                           , excessiveScaleRadii        , &
         &                                                                                                 useFallbackScaleMethod
    type            (rootFinder                         )           , save                              :: finder                                                                
    !$omp threadprivate(finder)
    class           (darkMatterProfileConcentrationClass), pointer  , save                              :: fallbackConcentration
    logical                                                         , save                              :: functionInitialized                      =.false.
    type            (varying_string                     )                                               :: mergerTreeReadConcentrationFallbackMethod        , message
    character       (len=16                             )                                               :: label

    ! Initialize if necessary.
    if (.not.functionInitialized) then
       !$omp critical(Assign_Scale_Radii_Initialize)
       if (.not.functionInitialized) then
          ! Construct the fallback concentration method.
          if (Input_Parameter_Is_Present('mergerTreeReadConcentrationFallbackMethod')) then
             !@ <inputParameter>
             !@   <name>mergerTreeReadConcentrationFallbackMethod</name>
             !@   <defaultValue>{\tt [darkMatterProfileConcentrationMethod]}</defaultValue>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The method to be used for setting node scale radii when reading merger trees from file and the node mass falls below the reliability threshold.
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>1</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('mergerTreeReadConcentrationFallbackMethod',mergerTreeReadConcentrationFallbackMethod)
             fallbackConcentration => darkMatterProfileConcentration(char(mergerTreeReadConcentrationFallbackMethod))
          else
             fallbackConcentration => darkMatterProfileConcentration()
          end if
          ! Record that we are now initialized.
          functionInitialized=.true.
       end if
       !$omp end critical(Assign_Scale_Radii_Initialize)
    end if
    ! Initialize our root finder.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(Half_Mass_Radius_Root              )
       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
    end if
    ! Find the scale radius.
    excessiveScaleRadii   =.false.
    excessiveHalfMassRadii=.false.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Assume that we need to use a fallback method to set halo scale radius.
          useFallbackScaleMethod=.true.
          ! Check if the node is sufficiently massive.
          thisBasicComponent             => nodeList(iIsolatedNode)%node%basic            (                 )
          thisDarkMatterProfileComponent => nodeList(iIsolatedNode)%node%darkMatterProfile(autoCreate=.true.)          
          if (thisBasicComponent%mass() >= mergerTreeReadPresetScaleRadiiMinimumMass) then
             ! Check if we have scale radii read directly from file.
             if     (                                                                &
                  &     nodes(iNode)%scaleRadius                                     &
                  &   >                                                              &
                  &     0.0d0                                                        &
                  &  .and.                                                           &
                  &     nodes(iNode)%scaleRadius                                     &
                  &   <                                                              &
                  &     Dark_Matter_Halo_Virial_Radius(nodeList(iIsolatedNode)%node) &
                  &    /mergerTreeReadPresetScaleRadiiConcentrationMinimum           &
                  &  .and.                                                           &
                  &     nodes(iNode)%scaleRadius                                     &
                  &   >                                                              &
                  &     Dark_Matter_Halo_Virial_Radius(nodeList(iIsolatedNode)%node) &
                  &    /mergerTreeReadPresetScaleRadiiConcentrationMaximum           &
                  & ) then
                ! We do, so simply use them to set the scale radii in tree nodes.
                call thisDarkMatterProfileComponent%scaleSet(nodes(iNode)%scaleRadius)
                useFallbackScaleMethod=.false.
             else if (nodes(iNode)%halfMassRadius > 0.0d0) then
                ! We do not have scale radii read directly. Instead, compute them from half-mass radii.
                ! Set the active node and target half mass radius.
                activeNode                       => nodeList(iIsolatedNode)%node
                activeDarkMatterProfileComponent => activeNode%darkMatterProfile()
                activeBasicComponent             => activeNode%basic            ()
                halfMassRadius                   =  nodes(iNode)%halfMassRadius
                ! Solve for the scale radius.
                call finder%rangeExpand    (                                                                                   &
                     &                      rangeExpandDownward          =0.5d0                                              , &
                     &                      rangeExpandUpward            =2.0d0                                              , &
                     &                      rangeDownwardLimit           = Dark_Matter_Halo_Virial_Radius(activeNode)          &
                     &                                                    /mergerTreeReadPresetScaleRadiiConcentrationMaximum, &
                     &                      rangeUpwardLimit             = Dark_Matter_Halo_Virial_Radius(activeNode)          &
                     &                                                    /mergerTreeReadPresetScaleRadiiConcentrationMinimum, &
                     &                      rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive                      , &
                     &                      rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative                      , &
                     &                      rangeExpandType              =rangeExpandMultiplicative                            &
                     &                     )
                radiusScale=finder%find(rootGuess=halfMassRadius,status=status)
                if (status == errorStatusSuccess) then
                   call thisDarkMatterProfileComponent%scaleSet(radiusScale)
                   ! Check for scale radii exceeding the virial radius.
                   if (radiusScale    > Dark_Matter_Halo_Virial_Radius(activeNode)) excessiveScaleRadii   =.true.
                   ! Check for half-mass radii exceeding the virial radius.
                   if (halfMassRadius > Dark_Matter_Halo_Virial_Radius(activeNode)) excessiveHalfMassRadii=.true.
                   useFallbackScaleMethod=.false.
                else
                   if (mergerTreeReadPresetScaleRadiiFailureIsFatal) then
                      messageVerbosity=Galacticus_Verbosity_Level()
                   else
                      messageVerbosity=verbosityWarn
                   end if
                   call Galacticus_Display_Indent  ("failed to find scale radius consistent with specified half-mass radius",messageVerbosity)
                   write (label,'(i16)') activeNode%hostTree%index
                   message="      tree index: "//trim(label)
                   call Galacticus_Display_Message(message,messageVerbosity)
                    write (label,'(i16)') activeNode%index()
                   message="      node index: "//trim(label)
                   call Galacticus_Display_Message(message,messageVerbosity)
                   write (label,'(e12.6)') Dark_Matter_Halo_Virial_Radius(activeNode)
                   message="   virial radius: "//trim(label)
                   call Galacticus_Display_Message(message,messageVerbosity)
                   write (label,'(e12.6)') halfMassRadius
                   message="half-mass radius: "//trim(label)
                   call Galacticus_Display_Message(message,messageVerbosity)
                   call Galacticus_Display_Unindent("",messageVerbosity)
                   if (mergerTreeReadPresetScaleRadiiFailureIsFatal) then
                      call Galacticus_Error_Report('Assign_Scale_Radii','aborting')
                   else 
                     useFallbackScaleMethod=.true.
                   end if
                end if
             end if
          end if
          if (useFallbackScaleMethod) then
             ! The node mass is below the reliability threshold, or no scale information is available. Set the scale radius using
             ! the fallback concentration method.
             radiusScale=                                                              &
                  &  Dark_Matter_Halo_Virial_Radius     (nodeList(iIsolatedNode)%node) &
                  & /fallBackConcentration%concentration(nodeList(iIsolatedNode)%node)
             call thisDarkMatterProfileComponent%scaleSet(radiusScale)
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
    use Galacticus_Error
    implicit none
    class           (nodeData          )         , dimension(:), intent(inout) :: nodes              
    type            (treeNodeList      )         , dimension(:), intent(inout) :: nodeList           
    class           (nodeComponentBasic), pointer                              :: thisBasicComponent 
    class           (nodeComponentSpin ), pointer                              :: thisSpinComponent  
    integer                                                                    :: iNode              
    integer         (kind=kind_int8    )                                       :: iIsolatedNode      
    double precision                                                           :: spin             , spinNormalization        
    double precision                             , dimension(3)                :: spin3D               
    
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Get basic and spin components.
          thisBasicComponent => nodeList(iIsolatedNode)%node%basic(                 )
          thisSpinComponent  => nodeList(iIsolatedNode)%node%spin (autoCreate=.true.)
          if (mergerTreeReadPresetSpins  ) then
             if      (defaultImporter%          spinAvailable()) then
                ! If spins are available directly, use them.
                call thisSpinComponent%spinSet(nodes(iNode)%spin)
             else if (defaultImporter%angularMomentaAvailable()) then
                ! Compute the spin parameter normalization to convert from angular momentum.
                spinNormalization= sqrt(abs(Dark_Matter_Profile_Energy(nodeList(iIsolatedNode)%node))) &
                     &            /gravitationalConstantGalacticus                                     &
                     &            /thisBasicComponent%mass()**2.5d0
                spin  = spinNormalization              &
                     & *nodes(iNode)%angularMomentum
                call thisSpinComponent%spinSet(spin)
             else
                call Galacticus_Error_Report('Assign_Spin_Parameters','no method exists to set spins')
             end if
          end if
          if (mergerTreeReadPresetSpins3D) then
             if      (defaultImporter%          spin3DAvailable()) then
                ! If spins are available directly, use them.
                call thisSpinComponent%spinVectorSet(nodes(iNode)%spin3D)
             else if (defaultImporter%angularMomenta3DAvailable()) then
                ! Compute the spin parameter normalization to convert from angular momentum.
                spin3D= spinNormalization              &
                     & *nodes(iNode)%angularMomentum3D
                call thisSpinComponent%spinVectorSet(spin3D)
             else
                call Galacticus_Error_Report('Assign_Spin_Parameters','no method exists to set vector spins')
             end if
          end if
       end if
    end do
    return
  end subroutine Assign_Spin_Parameters

  subroutine Assign_Particle_Counts(nodes,nodeList)
    !% Assign particle counts to nodes.
    implicit none
    class           (nodeData          )         , dimension(:), intent(inout) :: nodes              
    type            (treeNodeList      )         , dimension(:), intent(inout) :: nodeList           
    class           (nodeComponentNBody), pointer                              :: thisNBody
    integer                                                                    :: iNode              
    integer         (kind=kind_int8    )                                       :: iIsolatedNode      
    
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Get N-body component.
          thisNBody => nodeList(iIsolatedNode)%node%nBody(autoCreate=.true.)
          ! Assign the particle count.
          call thisNBody%particleCountSet(nodes(iNode)%particleCount)
       end if
    end do
    return
  end subroutine Assign_Particle_Counts

  subroutine Assign_Velocity_Maxima(nodes,nodeList)
    !% Assign velocity maxima to nodes.
    implicit none
    class           (nodeData          )         , dimension(:), intent(inout) :: nodes              
    type            (treeNodeList      )         , dimension(:), intent(inout) :: nodeList           
    class           (nodeComponentNBody), pointer                              :: thisNBody
    integer                                                                    :: iNode              
    integer         (kind=kind_int8    )                                       :: iIsolatedNode      
    
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Get N-body component.
          thisNBody => nodeList(iIsolatedNode)%node%nBody(autoCreate=.true.)
          ! Assign the velocity maximum.
          call thisNBody%velocityMaximumSet(nodes(iNode)%velocityMaximum)
       end if
    end do
    return
  end subroutine Assign_Velocity_Maxima

  subroutine Assign_Velocity_Dispersions(nodes,nodeList)
    !% Assign velocity dispersions to nodes.
    implicit none
    class           (nodeData          )         , dimension(:), intent(inout) :: nodes              
    type            (treeNodeList      )         , dimension(:), intent(inout) :: nodeList           
    class           (nodeComponentNBody), pointer                              :: thisNBody
    integer                                                                    :: iNode              
    integer         (kind=kind_int8    )                                       :: iIsolatedNode      
    
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Get N-body component.
          thisNBody => nodeList(iIsolatedNode)%node%nBody(autoCreate=.true.)
          ! Assign the velocity dispersion.
          call thisNBody%velocityDispersionSet(nodes(iNode)%velocityDispersion)
       end if
    end do
    return
  end subroutine Assign_Velocity_Dispersions
  
  double precision function Half_Mass_Radius_Root(radius)
    !% Function used to find scale radius of dark matter halos given their half-mass radius.
    use Dark_Matter_Profiles
    implicit none
    double precision, intent(in   ) :: radius 
    
    ! Set scale radius to current guess.
    call activeDarkMatterProfileComponent%scaleSet(radius)
    ! Compute difference between mass fraction enclosed at half mass radius and one half.
    Half_Mass_Radius_Root=Dark_Matter_Profile_Enclosed_Mass(activeNode,halfMassRadius)/activeBasicComponent%mass()-0.50d0
    return
  end function Half_Mass_Radius_Root

  subroutine Assign_Isolated_Node_Indices(nodes)
    !% Assign to each node the number of the corresponding isolated node.
    implicit none
    class  (nodeData      ), dimension(:), intent(inout) :: nodes         
    class  (nodeData      ), pointer                     :: thisNode      
    integer(kind=kind_int8)                              :: iIsolatedNode 
    integer                                              :: iNode         
    logical                                              :: endOfBranch   
    
    ! First make a copy of the currently assigned isolated node indices. These will be used
    ! later to reference the nodes which are the primary node associated with objects in nodeList.
    nodes%primaryIsolatedNodeIndex=nodes%isolatedNodeIndex
    ! Iterate over nodes.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Find the subset with descendents.
          if (associated(nodes(iNode)%descendent)) then
             ! Select the subset which have a subhalo as a descendent.
             if (nodes(iNode)%descendent%isSubhalo) then
                ! Trace descendents until merging or final time.
                thisNode   => nodes(iNode)%descendent
                endOfBranch=  .false.
                do while (.not.endOfBranch)
                   ! Record that this node was reachable via descendents of an isolated node.
                   if (thisNode%isolatedNodeIndex == nodeIsUnreachable) thisNode%isolatedNodeIndex=nodeIsReachable

                   if (.not.associated(thisNode%descendent)) then
                      ! If there is no descendent then the end of the branch has been reached.
                      endOfBranch=.true.
                   else
                      ! Step to the next descendent.
                      thisNode => thisNode%descendent
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
    use Kepler_Orbits
    use Dark_Matter_Halo_Scales
    use String_Handling
    use Galacticus_Error
    implicit none
    class           (nodeData                       )           , dimension(:), intent(inout), target :: nodes                                                                
    type            (treeNodeList                   )           , dimension(:), intent(inout)         :: nodeList                                                             
    integer         (kind=kind_int8                 )                         , intent(  out)         :: historyCountMaximum                                                  
    class           (nodeData                       ), pointer                                        :: lastSeenNode                       , progenitorNode              , & 
         &                                                                                               thisNode                                                             
    type            (treeNode                       ), pointer                                        :: firstProgenitor                    , hostNode                    , & 
         &                                                                                               orbitalPartner                     , satelliteNode                   
    double precision                                            , dimension(3)                        :: hostPosition                       , relativePosition            , & 
         &                                                                                               satellitePosition                                                    
    double precision                                            , dimension(3)                        :: hostVelocity                       , relativeVelocity            , & 
         &                                                                                               satelliteVelocity                                                    
    logical                                          , parameter                                      :: acceptUnboundOrbits        =.false.                                  
    class           (nodeComponentBasic             ), pointer                                        :: childBasicComponent                , orbitalPartnerBasicComponent, & 
         &                                                                                               satelliteBasicComponent            , thisBasicComponent              
    class           (nodeComponentPosition          ), pointer                                        :: childPositionComponent             , hostPositionComponent       , & 
         &                                                                                               satellitePositionComponent         , thisPositionComponent           
    class           (nodeComponentSatellite         ), pointer                                        :: satelliteSatelliteComponent        , thisSatelliteComponent          
    type            (keplerOrbit                    )                                                 :: thisOrbit                                                            
    integer                                                                                           :: iNode                                                                
    integer         (kind=kind_int8                 )                                                 :: historyCount                       , iIsolatedNode                   
    logical                                                                                           :: branchMerges                       , branchTipReached            , & 
         &                                                                                               endOfBranch                        , isolatedProgenitorExists    , & 
         &                                                                                               nodeWillMerge                                                        
    double precision                                                                                  :: radiusApocenter                    , radiusPericenter            , & 
         &                                                                                               radiusVirial                       , timeSubhaloMerges               
    type            (varying_string                 )                                                 :: message                                                              
    type            (progenitorIterator             )                                                 :: progenitors                                                          
    
    historyCountMaximum  = 0
    nodes%mergesWithIndex=-1
    do iNode=1,size(nodes)
       if (nodes(iNode)%primaryIsolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%primaryIsolatedNodeIndex
          ! Find the subset with descendents.
          if (associated(nodes(iNode)%descendent)) then
             ! Flag indicating if this is a node for which a merging time should be set.
             nodeWillMerge=.false.
             ! Select the subset which have a subhalo, or which are an initial subhalo.
             if (nodes(iNode)%descendent%isSubhalo.or.nodes(iNode)%isSubhalo) then
                ! Trace descendents until merging or final time.
                endOfBranch     =.false.
                branchTipReached=.false.
                branchMerges    =.false.
                historyCount    =0
                if (nodes(iNode)%isSubhalo) then
                   thisNode => nodes(iNode)
                else
                   ! Check for an immediate subhalo-subhalo merger.
                   if (Is_Subhalo_Subhalo_Merger(nodes,nodes(iNode))) then
                      endOfBranch =.true.
                      branchMerges=.true.
                      nodes(iNode)%mergesWithIndex=nodes(iNode)%descendent%nodeIndex
                      historyCount=historyCount+max(0_kind_int8,defaultImporter%subhaloTraceCount(nodes(iNode)))
                   end if
                   lastSeenNode => nodes(iNode)
                   thisNode     => nodes(iNode)%descendent
                end if
                do while (.not.endOfBranch)
                   ! Record which isolated node this node belongs to.
                   thisNode%isolatedNodeIndex=iIsolatedNode
                   ! Increment the history count for this branch.
                   historyCount=historyCount+1
                   ! Test the branch.
                   if (.not.associated(thisNode%descendent)) then
                      ! No descendent, indicating tip of branch has been reached
                      branchTipReached            =.true.
                      endOfBranch                 =.true.
                      historyCount                =historyCount+max(0_kind_int8,defaultImporter%subhaloTraceCount(thisNode))
                   else if (.not.thisNode%descendent%isSubhalo) then
                      ! Descendent is not a subhalo, treat as a merging event or a subhalo promotion.
                      endOfBranch                 =.true.
                      historyCount                =historyCount+max(0_kind_int8,defaultImporter%subhaloTraceCount(thisNode))
                      ! Search for any isolated progenitors of the node's descendent.
                      isolatedProgenitorExists=.false.
                      call progenitors%descendentSet(thisNode%descendent,nodes)
                      do while (progenitors%next(nodes) .and. .not.isolatedProgenitorExists)
                         progenitorNode => progenitors%current(nodes)
                         isolatedProgenitorExists=(progenitorNode%nodeIndex == progenitorNode%hostIndex)
                      end do
                      ! If an isolated progenitor exists, this is a merger event. If not, it is a subhalo promotion (which will be
                      ! handled elsewhere).
                      if (isolatedProgenitorExists) then
                         branchMerges                =.true.
                         nodes(iNode)%mergesWithIndex=thisNode%descendent%nodeIndex
                         lastSeenNode                => thisNode
                         thisNode                    => thisNode%descendent
                      end if
                   else
                      ! Merges with another subhalo.
                      call progenitors%descendentSet(thisNode%descendent,nodes)
                      do while (progenitors%next(nodes))
                         progenitorNode => progenitors%current(nodes)
                         if     (                                                                           &
                              &                    progenitorNode%nodeIndex         /= thisNode%nodeIndex  &
                              &  .and.             progenitorNode%isolatedNodeIndex /= nodeIsUnreachable   &
                              &  .and.  associated(progenitorNode%descendent                             ) &
                              &  .and.             progenitorNode%nodeMass           > thisNode%nodeMass   &
                              & ) then
                            ! Another node merges into current node's descendent subhalo and is more massive than current
                            ! node. Therefore, class this as a subhalo-subhalo merger.
                            branchMerges                =.true.
                            endOfBranch                 =.true.
                            nodes(iNode)%mergesWithIndex=progenitorNode%descendent%nodeIndex
                            historyCount                =historyCount+max(0_kind_int8,defaultImporter%subhaloTraceCount(thisNode))
                            lastSeenNode                => thisNode
                            thisNode                    => thisNode%descendent
                            exit
                         end if
                      end do
                      ! Step to the next descendent.
                      if (.not.endOfBranch) thisNode => thisNode%descendent
                   end if
                end do
                ! Only set a merging time if this node is not the primary progenitor of its parent.
                if (.not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                   ! Record the largest history.
                   historyCountMaximum=max(historyCountMaximum,historyCount)
                   ! Set an appropriate merging time for this subhalo.
                   if      (branchTipReached) then
                      timeSubhaloMerges=timeUntilMergingInfinite ! Subhalo never merges, so set merging time to effective infinity.
                   else if (branchMerges    ) then
                      ! Find the time of merging, accounting for any additional (subresolution) time.
                      timeSubhaloMerges=thisNode%nodeTime
                      call Time_Until_Merging_Subresolution(lastSeenNode,nodes,iNode,timeSubhaloMerges)
                   else
                      ! Neither the branch tip was reached, not does this branch merge. Therefore, this must be a subhalo which is
                      ! promoted to be an isolated halo. Simply set an infinite merging time as we do not wish this node to merge.
                      timeSubhaloMerges=timeUntilMergingInfinite
                   end if
                   ! Flag that this node will merge.
                   nodeWillMerge=.true.
                end if
             else if (.not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                ! Descendent is not a subhalo but this node is not the primary progenitor. Assume instantaneous merging.
                thisBasicComponent => nodeList(iIsolatedNode)%node%basic()
                timeSubhaloMerges=thisBasicComponent%time()
                ! Flag that this node will merge.
                nodeWillMerge=.true.
                ! Record the node with which the merger occurs.
                nodes(iNode)%mergesWithIndex=nodes(iNode)%descendent%nodeIndex
                ! Ensure the history arrays will be large enough to hold data for this node.
                historyCountMaximum=max(historyCountMaximum,max(0_kind_int8,defaultImporter%subhaloTraceCount(nodes(iNode))))
                ! Account for any subresolution merging time.
                call Time_Until_Merging_Subresolution(nodes(iNode),nodes,iNode,timeSubhaloMerges)
             end if
             ! Set a merging time and/or orbit if this node will merge.
             if (nodeWillMerge.and.mergerTreeReadPresetMergerTimes) then
                ! Store the time of merging for this node and all of its primary progenitors.
                firstProgenitor => nodeList(iIsolatedNode)%node
                do while (associated(firstProgenitor))
                   thisSatelliteComponent => firstProgenitor%satellite(autoCreate=.true.)
                   call thisSatelliteComponent%timeOfMergingSet(timeSubhaloMerges)
                   firstProgenitor => firstProgenitor%firstChild
                end do
             end if
          end if
          ! Set position and velocity if required.
          if (mergerTreeReadPresetPositions) then
             thisPositionComponent => nodeList(iIsolatedNode)%node%position(autoCreate=.true.)
             call thisPositionComponent%positionSet(nodes(iNode)%position)
             call thisPositionComponent%velocitySet(nodes(iNode)%velocity)
             ! Detect if the node parent has no isolated child - in which case one will have been made for it using a
             ! direct copy of itself. Note that this is an ugly solution - once trees can handle nodes with no primary
             ! progenitor (but with secondary progenitors) a cleaner test could be used here.
             if (associated(nodeList(iIsolatedNode)%node%firstChild)) then
                thisBasicComponent  => nodeList(iIsolatedNode)%node           %basic()
                childBasicComponent => nodeList(iIsolatedNode)%node%firstChild%basic()
                if (nodeList(iIsolatedNode)%node%uniqueID() == nodeList(iIsolatedNode)%node%firstChild%uniqueID()) then
                   ! Set the position and velocity of the pseudo-primary progenitor here also.
                   childPositionComponent => nodeList(iIsolatedNode)%node%firstChild%position(autoCreate=.true.)
                   call childPositionComponent%positionSet(nodes(iNode)%position)
                   call childPositionComponent%velocitySet(nodes(iNode)%velocity)
                end if
             end if
          end if
       end if
    end do

    ! Set orbits.
    if (mergerTreeReadPresetOrbits) then
       iIsolatedNode=0
       do iNode=1,size(nodes)
         if (nodes(iNode)%primaryIsolatedNodeIndex /= nodeIsUnreachable) then
             iIsolatedNode=nodes(iNode)%primaryIsolatedNodeIndex
            ! Set the orbit for this halo.
             satelliteNode => nodeList(iIsolatedNode)%node
             if (associated(satelliteNode%parent).and..not.satelliteNode%isPrimaryProgenitor()) then
                ! Find the orbital partner.
                hostNode => satelliteNode%parent%firstChild
                ! If the parent node has no progenitors, then we are forced to use the parent node
                ! itself as the orbital partner.
                if (associated(hostNode)) then
                   orbitalPartner => hostNode     %parent
                else
                   hostNode       => satelliteNode%parent
                   orbitalPartner => hostNode
                end if
                ! Get components.
                satelliteBasicComponent      => satelliteNode             %basic    (                 )
                satellitePositionComponent   => satelliteNode             %position (                 )
                satelliteSatelliteComponent  => satelliteNode             %satellite(autoCreate=.true.)
                hostPositionComponent        => hostNode                  %position (                 )
                orbitalPartnerBasicComponent => orbitalPartner            %basic    (                 )
                ! Get position and velocity.
                satellitePosition            =  satellitePositionComponent%position (                 )
                satelliteVelocity            =  satellitePositionComponent%velocity (                 )
                hostPosition                 =       hostPositionComponent%position (                 )
                hostVelocity                 =       hostPositionComponent%velocity (                 )
                ! Find relative position and velocity.
                relativePosition=satellitePosition-hostPosition
                relativeVelocity=satelliteVelocity-hostVelocity
                ! Update position/velocity for periodicity and Hubble flow.
                call Phase_Space_Position_Realize(satelliteBasicComponent%time(),relativePosition,relativeVelocity)
                ! Catch zero separation halos.
                if (Vector_Magnitude(relativePosition) == 0.0d0) then
                   if (mergerTreeReadPresetOrbitsSetAll) then
                      ! The given orbit does not cross the virial radius. Since all orbits must be set, choose an orbit at random.
                      thisOrbit=Virial_Orbital_Parameters(satelliteNode,hostNode,acceptUnboundOrbits)
                      call satelliteSatelliteComponent%virialOrbitSet(thisOrbit)
                   else 
                      message='merging halos ['
                      message=message//satelliteNode%index()//' & '//hostNode%index()//'] have zero separation'
                      call Galacticus_Error_Report('Scan_For_Mergers',message)
                   end if
                else
                   ! Create the orbit.
                   thisOrbit=Orbit_Construct(satelliteBasicComponent%mass(),orbitalPartnerBasicComponent%mass(),relativePosition,relativeVelocity)
                   ! Propagate to the virial radius.
                   radiusPericenter=thisOrbit%radiusPericenter()
                   radiusApocenter =thisOrbit%radiusApocenter ()
                   radiusVirial    =Dark_Matter_Halo_Virial_Radius(orbitalPartner)
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
                      call satelliteSatelliteComponent%virialOrbitSet(thisOrbit)
                   else if (mergerTreeReadPresetOrbitsSetAll) then
                      ! The given orbit does not cross the virial radius. Since all orbits must be set, choose an orbit at random.
                      thisOrbit=Virial_Orbital_Parameters(satelliteNode,hostNode,acceptUnboundOrbits)
                      call satelliteSatelliteComponent%virialOrbitSet(thisOrbit)
                   else if (mergerTreeReadPresetOrbitsAssertAllSet) then
                      message='virial orbit could not be set for node '
                      message=message//satelliteNode%index()//char(10)
                      message=message//' -> set [mergerTreeReadPresetOrbitsAssertAllSet]=false to ignore this problem'//char(10)
                      message=message//'    (this may lead to other problems)'
                      call Galacticus_Error_Report('Scan_For_Mergers',message)
                   end if
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
    class  (nodeData      ), dimension(:), intent(inout) :: nodes    
    type   (treeNodeList  ), dimension(:), intent(inout) :: nodeList 
    integer                                              :: iNode    
    integer(kind=kind_int8)                              :: jNode    
    type   (varying_string)                              :: message  
    
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
                ! Set pointer from merging node (a.k.a. the "mergee") to node that will be merged with.
                nodeList(nodes(iNode)%isolatedNodeIndex)%node%mergeTarget => nodeList(nodes(jNode)%isolatedNodeIndex)%node

                ! Make a backward pointer from the merge target to the mergee. Check if the target already has mergees associated with it.
                if (associated(nodeList(nodes(jNode)%isolatedNodeIndex)%node%firstMergee)) then
                   ! It does: unlink them and attached to the "siblingMergee" pointer of the current mergee.
                   nodeList(nodes(iNode)%isolatedNodeIndex)%node%siblingMergee => nodeList(nodes(jNode)%isolatedNodeIndex)%node%firstMergee
                else
                   ! It does not: simply nullify the next mergee pointer of the mergee.
                   nodeList(nodes(iNode)%isolatedNodeIndex)%node%siblingMergee => null()
                end if
                ! Append the mergee as the first mergee on the target node.
                nodeList(nodes(jNode)%isolatedNodeIndex)%node%firstMergee => nodeList(nodes(iNode)%isolatedNodeIndex)%node
             end if
          end if
       end do
    end if
    return
  end subroutine Assign_Mergers

  subroutine Scan_for_Branch_Jumps(nodes,nodeList)
    !% Search for subhalos which move between branches/trees.
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    implicit none
    class           (nodeData      ), dimension(:), intent(inout), target :: nodes                                                        
    type            (treeNodeList  ), dimension(:), intent(inout)         :: nodeList                                                     
    class           (nodeData      ), pointer                             :: currentHost  , descendentNode, hostDescendent, jumpToHost, & 
         &                                                                   previousNode                                                 
    integer                                                               :: iNode                                                        
    integer         (kind=kind_int8)                                      :: iIsolatedNode                                                
    logical                                                               :: isMergerEvent, subhaloJumps  , wasMergerEvent                
    double precision                                                      :: timeOfJump                                                   
    type            (varying_string)                                      :: message                                                      
    
    ! If branch jumps are not allowed, simply return.
    if (.not.mergerTreeReadAllowBranchJumps) return
    ! Search for subhalos whose descendents live in a different host than that to which their
    ! host descends. These subhalos are jumping between tree branches (or between trees). Add
    ! an event to such nodes to handle the jump.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node (or an initial satellite).
       if (nodes(iNode)%primaryIsolatedNodeIndex /= nodeIsUnreachable) then
          iIsolatedNode=nodes(iNode)%primaryIsolatedNodeIndex
          ! Find those which are a subhalo, or whose descendent is a subhalo.
          descendentNode => null()
          if      (           nodes(iNode)%isSubhalo  ) then
             descendentNode => nodes(iNode)
             previousNode   => nodes(iNode)
          else if (associated(nodes(iNode)%descendent)) then
             if (nodes(iNode)%descendent%isSubhalo) then
                descendentNode => nodes(iNode)%descendent
                previousNode   => nodes(iNode)
             end if
          end if
          ! Check for an immediate subhalo-subhalo merger. If found, nullify the descendent,
          ! so we do not attempt to process this branch.
          if (Is_Subhalo_Subhalo_Merger(nodes,nodes(iNode))) then
             descendentNode => null()
             currentHost => Last_Host_Descendent(nodes(iNode))
             ! Add a jump if the tree ends before the descendent time.
             if (currentHost%nodeTime <= nodes(iNode)%descendent%nodeTime) then
                timeOfJump     =  currentHost%nodeTime
                jumpToHost     => nodes(iNode)%descendent%host
                do while (jumpToHost%isSubhalo)
                   if (.not.warningNestedHierarchyIssued) then
                      message='nested hierarchy detected [node '
                      message=message//nodes(iNode)%descendent%nodeIndex//']'
                      message=message//char(10)//'ignoring as not currently supported'
                      message=message//char(10)//'warning will not be issued again'
                      call Galacticus_Display_Message(message,verbosityWarn)
                      warningNestedHierarchyIssued=.true.
                   end if
                   jumpToHost => jumpToHost%host
                end do
                call Create_Branch_Jump_Event(                                                    &
                     &                        nodeList(iIsolatedNode                      )%node, &
                     &                        nodeList(jumpToHost%primaryIsolatedNodeIndex)%node, &
                     &                        timeOfJump                                          &
                     &                       )
             end if
          end if
          ! If a subhalo was found, follow its descent.
          wasMergerEvent=.false.
          do while (associated(descendentNode))
             subhaloJumps=.false.
             timeOfJump  =-1.0d0
             if (descendentNode%isSubhalo.and.associated(descendentNode%descendent)) then
                ! Determine if this is actually a merger event rather than a branch jump.
                ! Assume it is not a merger initially.
                isMergerEvent=.false.
                if (descendentNode%descendent%isSubhalo) then
                   ! Descendent is a subhalo. Check for subhalo-subhalo merger.
                   if (Is_Subhalo_Subhalo_Merger(nodes,previousNode)) isMergerEvent=.true.
                else
                   ! Descendent is not a subhalo, so this must be a merger event.
                   isMergerEvent=.true.
                end if
                ! If this is a merger event, then check that the current descendent's host has a
                ! descendent that exists beyond the time of the merger. If it does not, then we
                ! still need to allow our node to jump branches (if necessary) as it will not be
                ! able to evolve in the descendentless host.
                wasMergerEvent=isMergerEvent
                if (isMergerEvent) then
                   currentHost => Last_Host_Descendent(descendentNode)
                   if (currenthost%nodeTime <= descendentNode%descendent%nodeTime) then
                      isMergerEvent=.false.
                      timeOfJump=currentHost%nodeTime
                   endif
                end if
                ! Proceed only if this is not a merger event.
                if (.not.isMergerEvent) then
                   ! Does this subhalo's descendent live in the host to which the subhalo's host descends.
                   if (.not.associated(descendentNode%host%descendent)) then
                      ! Host has no descendent, so this must be a branch jump.
                      subhaloJumps=.true.
                   else if (descendentNode%descendent%hostIndex /= descendentNode%host%descendent%hostIndex) then
                      ! Host has a descendent, but it's host is not the same as our descendent's host.
                      subhaloJumps=.true.
                      ! Check that is not simply a case of the subhalo skipping one or more timesteps before
                      ! reappearing in the expected host.
                      hostDescendent => descendentNode%host%descendent%host
                      do while (descendentNode%descendent%host%nodeTime > hostDescendent%nodeTime)
                         if (associated(hostDescendent%descendent)) then
                            hostDescendent => hostDescendent%descendent%host
                         else
                            exit
                         end if
                      end do
                      ! Subhalo reappeared in the expected host. This is not a branch jump.
                      if (descendentNode%descendent%hostIndex == hostDescendent%nodeIndex) subhaloJumps=.false.
                   end if
                else
                   ! Since this is a merger event, we're finished checking this branch.
                   exit
                end if
             end if
             ! If a jump was detected, create an event.
             if (subhaloJumps) then
                if (timeOfJump < 0.0d0)                   &
                     & timeOfJump=descendentNode%nodeTime
                jumpToHost => descendentNode%descendent%host
                ! Find an isolated host.
                do while (jumpToHost%isSubhalo)
                   jumpToHost => jumpToHost%host
                end do
                call Create_Branch_Jump_Event(                                                    &
                     &                        nodeList(iIsolatedNode                      )%node, &
                     &                        nodeList(jumpToHost%primaryIsolatedNodeIndex)%node, &
                     &                        timeOfJump                                          &
                     &                       )
             end if
             ! Move to the descendent.
             previousNode   => descendentNode
             descendentNode => descendentNode%descendent
             ! If the descendent is not a subhalo, then we're finished checking this branch.
             if (associated(descendentNode)) then
                if (.not.descendentNode%isSubhalo) exit
             end if
             ! If this was a merger event, then we're finished checking this branch.
             if (wasMergerEvent) exit
          end do
       end if
    end do
    return
  end subroutine Scan_for_Branch_Jumps

  function Last_Host_Descendent(thisNode) result (currentHost)
    !% Return a pointer to the last descendent that can be reached from {\tt thisNode} when descending through hosts.
    implicit none
    class(nodeData), pointer               :: currentHost 
    class(nodeData), intent(inout), target :: thisNode    
    
    currentHost => thisNode%host
    do while(associated(currentHost%descendent))
       currentHost => currentHost%descendent%host
    end do
    return
  end function Last_Host_Descendent

  subroutine Create_Branch_Jump_Event(thisNode,jumpToHost,timeOfJump)
    !% Create a matched-pair of branch jump events in the given nodes.
    use Node_Branch_Jumps
    implicit none
    type            (treeNode ), intent(inout), pointer :: jumpToHost, thisNode  
    double precision           , intent(in   )          :: timeOfJump            
    type            (nodeEvent)               , pointer :: newEvent  , pairEvent 
    
    newEvent       => thisNode%createEvent()
    newEvent %time =  timeOfJump
    newEvent %node => jumpToHost
    newEvent %task => Node_Branch_Jump
    pairEvent      => jumpToHost%createEvent()
    pairEvent%time =  timeOfJump
    pairEvent%node => thisNode
    pairEvent%task => null()
    pairEvent%ID   =  newEvent%ID
    return
  end subroutine Create_Branch_Jump_Event

  subroutine Build_Subhalo_Mass_Histories(nodes,nodeList,historyTime,historyMass,position,velocity)
    !% Build and attached bound mass histories to subhalos.
    use Galacticus_Error
    use String_Handling
    use Cosmology_Functions
    use Histories
    implicit none
    class           (nodeData               )         , dimension(:  ), intent(inout), target :: nodes                                    
    type            (treeNodeList           )         , dimension(:  ), intent(inout)         :: nodeList                                 
    double precision                                  , dimension(:  ), intent(inout)         :: historyMass              , historyTime   
    double precision                                  , dimension(:,:), intent(inout)         :: position                 , velocity      
    class           (nodeData               ), pointer                                        :: progenitorNode           , thisNode      
    type            (treeNode               ), pointer                                        :: firstProgenitor                          
    class           (nodeComponentSatellite ), pointer                                        :: thisSatelliteComponent                   
    class           (nodeComponentPosition  ), pointer                                        :: thisPositionComponent                    
    class           (cosmologyFunctionsClass), pointer                                        :: cosmologyFunctionsDefault                
    integer         (kind=kind_int8         )                                                 :: historyCount             , iIsolatedNode 
    integer                                                                                   :: iNode                                    
    logical                                                                                   :: endOfBranch                              
    type            (varying_string         )                                                 :: message                                  
    type            (history                )                                                 :: subhaloHistory                           
    type            (progenitorIterator     )                                                 :: progenitors                              
    
    if (mergerTreeReadPresetSubhaloMasses.or.mergerTreeReadPresetPositions) then
       ! Check that preset subhalo masses are supported.
       if (mergerTreeReadPresetSubhaloMasses.and..not.defaultSatelliteComponent%boundMassHistoryIsSettable()) &
            & call Galacticus_Error_Report('Merger_Tree_Read_Do','presetting subhalo masses requires a component that supports setting of node bound mass histories')
       ! Get the default cosmology functions object.
       cosmologyFunctionsDefault => cosmologyFunctions()
       historyBuildNodeLoop: do iNode=1,size(nodes)
          historyBuildIsolatedSelect: if (nodes(iNode)%primaryIsolatedNodeIndex /= nodeIsUnreachable) then
             iIsolatedNode=nodes(iNode)%primaryIsolatedNodeIndex
             ! Find the subset with descendents.
             historyBuildHasDescendentSelect: if (associated(nodes(iNode)%descendent)) then
                ! Set a pointer to the current node - this will be updated if any descendents are traced.
                thisNode => nodes(iNode)
                  ! Set initial number of times in the history to zero.
                historyCount=0
                ! Select the subset which have a subhalo as a descendent and are not the primary progenitor.
                historyBuildSubhaloSelect: if (nodes(iNode)%descendent%isSubhalo) then
                   ! Trace descendents until merging or final time.
                   thisNode    => nodes(iNode)%descendent
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
                      if (.not.associated(thisNode%descendent).or..not.thisNode%descendent%isSubhalo) then
                         ! End of branch reached.
                         endOfBranch=.true.
                      else
                         ! Check if merges with another subhalo.
                         call progenitors%descendentSet(thisNode%descendent,nodes)
                         do while (progenitors%next(nodes))
                            progenitorNode => progenitors%current(nodes)
                            if     (                                                                          &
                                 &                    progenitorNode%nodeIndex         /= thisNode%nodeIndex  &
                                 &  .and.             progenitorNode%isolatedNodeIndex /= nodeIsUnreachable   &
                                 &  .and.  associated(progenitorNode%descendent                             ) &
                                 &  .and.             progenitorNode%nodeMass           > thisNode%nodeMass   &
                                 & ) then
                               ! Subhalo-subhalo merger.
                               endOfBranch =.true.
                               exit
                            end if
                         end do
                         ! Step to the next descendent.
                         if (.not.endOfBranch) thisNode => thisNode%descendent
                      end if
                   end do historyBuildBranchWalk
                   ! Set the mass history for this node.
                   if (mergerTreeReadPresetSubhaloMasses) then
                      call subhaloHistory%destroy()
                      call subhaloHistory%create(1,int(historyCount))
                      subhaloHistory%time(:  )=historyTime(1:historyCount)
                      subhaloHistory%data(:,1)=historyMass(1:historyCount)
                      firstProgenitor        => nodeList(iIsolatedNode)%node%earliestProgenitor()
                      thisSatelliteComponent => firstProgenitor%satellite()
                      call thisSatelliteComponent%boundMassHistorySet(subhaloHistory)
                   end if
                end if historyBuildSubhaloSelect

                ! Set the position history for this node.
                if (mergerTreeReadPresetPositions.and..not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                   ! Check if particle data is available for this node.
                   if (defaultImporter%subhaloTraceCount(thisNode) > 0) then
                      ! Check that arrays are large enough to hold particle data. They should be. If they are not, it's a
                      ! bug.
                      if (historyCount+defaultImporter%subhaloTraceCount(thisNode) > size(historyTime)) then
                         message='history arrays are too small to hold data for node '
                         message=message//nodeList(iIsolatedNode)%node%index()//': ['//historyCount//'+'//defaultImporter%subhaloTraceCount(thisNode)//']='//(historyCount+defaultImporter%subhaloTraceCount(thisNode))//'>'//size(historyTime)
                         call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
                      end if
                      ! Read subhalo position trace data.
                      call defaultImporter%subhaloTrace                                               &
                           & (                                                                        &
                           &  thisNode                                                              , &
                           &  historyTime(  historyCount+1:historyCount+defaultImporter%subhaloTraceCount(thisNode)), &
                           &  position   (:,historyCount+1:historyCount+defaultImporter%subhaloTraceCount(thisNode)), &
                           &  velocity   (:,historyCount+1:historyCount+defaultImporter%subhaloTraceCount(thisNode))  &
                           & )
                      ! Increment the history count for this node.
                      historyCount=historyCount+defaultImporter%subhaloTraceCount(thisNode)
                   end if
                   if (historyCount > 0) then
                      call subhaloHistory%destroy()
                      call subhaloHistory%create(6,int(historyCount))
                      subhaloHistory%time(:    )=          historyTime(    1:historyCount)
                      subhaloHistory%data(:,1:3)=transpose(position   (1:3,1:historyCount))
                      subhaloHistory%data(:,4:6)=transpose(velocity   (1:3,1:historyCount))
                      thisPositionComponent => nodeList(iIsolatedNode)%node%position()
                      call thisPositionComponent%positionHistorySet(subhaloHistory)
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
    class  (nodeData          ), dimension(:), intent(inout) :: nodes                            
    type   (treeNode          ), pointer                     :: newNode          , thisSatellite 
    class  (nodeComponentBasic), pointer                     :: newBasicComponent                
    integer                                                  :: iNode                            
    integer(kind=kind_int8    )                              :: iIsolatedNode                    
    
    ! Search for cases where a node has no progenitors which do not descend into subhalos.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeIsUnreachable .and. .not.nodes(iNode)%isSubhalo) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Select nodes with parents.
          if (associated(nodes(iNode)%node%parent)) then
             ! Select nodes with subhalo descendents which are also the primary progenitor of their parent.
             if (nodes(iNode)%descendent%isSubhalo.and.associated(nodes(iNode)%node%parent%firstChild,nodes(iNode)%node)) then
                ! Insert a copy of the parent node as its own primary progenitor. This avoids current node being promoted into its
                ! parent even though it is intended to descend into a subhalo. The copy is shifted to a very slightly earlier
                ! time to avoid having two identical halos existing simultaneously (which can be problematic if outputting
                ! quantities which use the node index as a label in dataset names for example).
                allocate(newNode)
                call nodes(iNode)%node%parent%copyNodeTo(newNode)
                newNode%sibling                         => nodes(iNode)%node
                newNode%parent                          => nodes(iNode)%node%parent
                newNode%firstChild                      => null()
                nodes(iNode)%node%parent%firstChild     => newNode
                newBasicComponent                       => newNode%basic()
                call newBasicComponent%timeSet(newBasicComponent%time()*(1.0d0-1.0d-6))
                ! Any satellites are now attached to the copy.
                nodes(iNode)%node%parent%firstSatellite => null()
                thisSatellite => newNode%firstSatellite
                do while (associated(thisSatellite))
                   thisSatellite%parent => newNode
                   thisSatellite        => thisSatellite%sibling
                end do
             end if
          end if
       end if
    end do
    return
  end subroutine Validate_Isolated_Halos

  subroutine Assign_UniqueIDs_To_Clones(nodeList)
    !% Assign new uniqueID values to any cloned nodes inserted into the trees.
    implicit none
    type   (treeNodeList), dimension(:), intent(inout) :: nodeList 
    integer                                            :: iNode    
    
    do iNode=1,size(nodeList)
       if (associated(nodeList(iNode)%node%firstChild)) then
          if (nodeList(iNode)%node%uniqueID() == nodeList(iNode)%node%firstChild%uniqueID()) &
               &  call nodeList(iNode)%node%firstChild%uniqueIDSet()
       end if
    end do
    return
  end subroutine Assign_UniqueIDs_To_Clones

  logical function Is_Subhalo_Subhalo_Merger(nodes,thisNode)
    !% Returns true if {\tt thisNode} undergoes a subhalo-subhalo merger.
    implicit none
    class(nodeData          ), dimension(:), intent(inout) :: nodes          
    class(nodeData          )              , intent(in   ) :: thisNode       
    class(nodeData          ), pointer                     :: progenitorNode 
    type (progenitorIterator)                              :: progenitors    
    
    Is_Subhalo_Subhalo_Merger=.false.
    ! Return immediately if there is no descendent. (Since there can be no merger if there is no descendent.)
    if (.not.associated(thisNode%descendent          )) return
    ! Return immediately if descendent is not a subhalo, as this could then not be a subhalo-subhalo merger.
    if (.not.           thisNode%descendent%isSubhalo ) return
    ! Check if thisNode's descendent has any progenitor nodes.
    call progenitors%descendentSet(thisNode%descendent,nodes)
    do while (progenitors%next(nodes))
       progenitorNode => progenitors%current(nodes)
       if     (                                                                          &
            &                    progenitorNode%nodeIndex         /= thisNode%nodeIndex  &
            &  .and.             progenitorNode%isolatedNodeIndex /= nodeIsUnreachable   &
            &  .and.  associated(progenitorNode%descendent                             ) &
            &  .and.             progenitorNode%nodeMass           > thisNode%nodeMass   &
            & ) then
          ! It does, so this is a subhalo-subhalo merger.
          Is_Subhalo_Subhalo_Merger =.true.
          exit
       end if
    end do
    return
  end function Is_Subhalo_Subhalo_Merger

  subroutine Dump_Tree(nodes,highlightNodes,branchRoot)
    !% Dumps the tree structure to a file in a format suitable for processing with \href{http://www.graphviz.org/}{\sc dot}.
    implicit none
    class    (nodeData      ), dimension(:), intent(in   ), target   :: nodes                 
    integer  (kind=kind_int8), dimension(:), intent(in   ), optional :: highlightNodes        
    integer  (kind=kind_int8)              , intent(in   ), optional :: branchRoot            
    class    (nodeData      ), pointer                               :: thisNode              
    integer                                                          :: fileUnit      , iNode 
    character(len=20        )                                        :: color         , style 
    logical                                                          :: outputNode            
    
    ! Open an output file and write the GraphViz opening.
    open(newunit=fileUnit,file='mergerTreeConstructReadTree.gv',status='unknown',form='formatted')
    write (fileUnit,*) 'digraph Tree {'

    ! Loop over all nodes.
    do iNode=1,size(nodes)
       ! Determine if node is in the branch to be output.
       if (present(branchRoot)) then
          outputNode=.false.
          thisNode => nodes(iNode)
          do while (associated(thisNode%descendent))
             if (thisNode%nodeIndex == branchRoot) then
                outputNode=.true.
                exit
             end if
             thisNode => thisNode%descendent
          end do
          outputNode=outputNode.or.(thisNode%nodeIndex == branchRoot)
       else
          outputNode=.true.
       end if
       if (outputNode) then
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
             write (fileUnit,'(a,i16.16,a,i16.16,a,f5.2,a,a,a,a,a,f5.2,a)') '"',nodes(iNode)%nodeIndex,'" [shape=box   , label="',nodes(iNode)%nodeIndex,':',nodes(iNode)%nodeTime,'", color=',trim(color),', style=',trim(style),', z=',nodes(iNode)%nodeTime,'];'
             ! If a host node is given, add a link to it as a red line.
             if (associated(nodes(iNode)%host)) write (fileUnit,'(a,i16.16,a,i16.16,a)') '"',nodes(iNode)%nodeIndex,'" -> "',nodes(iNode)%host%nodeIndex,'" [color=red];'
          else
             write (fileUnit,'(a,i16.16,a,i16.16,a,f5.2,a,a,a,a,a,f5.2,a)') '"',nodes(iNode)%nodeIndex,'" [shape=circle, label="',nodes(iNode)%nodeIndex,':',nodes(iNode)%nodeTime,'", color=',trim(color),', style=',trim(style),', z=',nodes(iNode)%nodeTime,'];'
          endif
          ! Make a link to the descendent node using a black line.
          if (associated(nodes(iNode)%descendent)) write (fileUnit,'(a,i16.16,a,i16.16,a)') '"',nodes(iNode)%nodeIndex,'" -> "',nodes(iNode)%descendent%nodeIndex,'" ;'
       end if
    end do

    ! Close the file.
    write (fileUnit,*) '}'
    close(fileUnit)
    return
  end subroutine Dump_Tree

  subroutine Time_Until_Merging_Subresolution(lastSeenNode,nodes,iNode,timeSubhaloMerges)
    !% Compute the additional time until merging after a subhalo is lost from the tree (presumably due to limited resolution).
    use Vectors
    use Kepler_Orbits
    use Galacticus_Error
    use Satellite_Merging_Timescales
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    implicit none
    class           (nodeData                       )                             , intent(in   ) ::      lastSeenNode                                                                      
    class           (nodeData                       )               , dimension(:), intent(inout) , target::                                      nodes                                     
    integer                                                                       , intent(in   ) ::      iNode                                                                             
    double precision                                                              , intent(inout) ::      timeSubhaloMerges                                                                 
    class           (nodeData                       ), pointer                                    ::      primaryProgenitor                                    , progenitorNode         , & 
         &                                                                                                thisNode                                                                          
    class           (satelliteMergingTimescalesClass), pointer, save                              ::      thisSatelliteMergingTimescales                                                    
    type            (treeNode                       ), pointer                                    ::      hostNode                                             , satelliteNode              
    class           (nodeComponentBasic             ), pointer                                    ::      hostBasicComponent                                   , satelliteBasicComponent    
    double precision                                                , dimension(3)                ::      relativePosition                                     , relativeVelocity           
    logical                                                   , save                              ::      functionInitialized                          =.false.                             
    type            (keplerOrbit                    )                                             ::      thisOrbit                                                                         
    double precision                                                                              ::      primaryProgenitorMass                                , timeUntilMerging           
    type            (varying_string                 )                                             ::      mergerTreeReadSubresolutionMergingMethod                                          
    type            (progenitorIterator             )                                             ::      progenitors                                                                       
    character       (len=42                         )                                             ::      coordinateLabel                                                                   
    type            (varying_string                 )                                             ::      message                                                                           
    
    ! Initialize if necessary.
    if (.not.functionInitialized) then
       !$omp critical(Time_Until_Merging_Subresolution_Initialize)
       if (.not.functionInitialized) then
          ! Construct the satellite merging timescale object.
          !@ <inputParameter>
          !@   <name>mergerTreeReadSubresolutionMergingMethod</name>
          !@   <defaultValue>null</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of a satellite merging timescale method to be used for computing the extra time until merging for subhalos.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeReadSubresolutionMergingMethod',mergerTreeReadSubresolutionMergingMethod,defaultValue='null')
          thisSatelliteMergingTimescales => satelliteMergingTimescales(char(mergerTreeReadSubresolutionMergingMethod))
          ! Record that we are now initialized.
          functionInitialized=.true.
       end if
       !$omp end critical(Time_Until_Merging_Subresolution_Initialize)
    end if
    ! Find the nodes that descendent into our target node's descendent.
    call progenitors%descendentSet(lastSeenNode%descendent,nodes)
    if (progenitors%exist()) then
       primaryProgenitor     => null()
       primaryProgenitorMass =  0.0d0
       do while (progenitors%next(nodes))
          progenitorNode => progenitors%current(nodes)
          if (progenitorNode%nodeIndex /= lastSeenNode%nodeIndex .and. progenitorNode%nodeMass > primaryProgenitorMass) then
             primaryProgenitorMass =  progenitorNode%nodeMass
             primaryProgenitor     => progenitorNode
          end if
       end do
       ! Initialize time until merging to zero.
       timeUntilMerging=0.0d0
       ! If position information is available, compute the subresolution orbit.
       if (mergerTreeReadPresetPositions) then
          ! Find relative position and velocity.
          relativePosition=lastSeenNode%position-primaryProgenitor%position
          relativeVelocity=lastSeenNode%velocity-primaryProgenitor%velocity
          ! Update position/velocity for periodicity and Hubble flow.
          call Phase_Space_Position_Realize(lastSeenNode%nodeTime,relativePosition,relativeVelocity)
          ! Catch zero separation halos.
          if (Vector_Magnitude(relativePosition) == 0.0d0) then
             message='merging halos ['
             message=message//lastSeenNode%nodeIndex//' & '//primaryProgenitor%nodeIndex//'] have zero separation'
             call Galacticus_Display_Indent  (message                          ,verbosityWarn)
             write (coordinateLabel,'("[",e12.6,",",e12.6,",",e12.6,"]")') primaryProgenitor%position
             message="position [primary  ] = "//trim(coordinateLabel)
             call Galacticus_Display_Message (message                          ,verbosityWarn)
             write (coordinateLabel,'("[",e12.6,",",e12.6,",",e12.6,"]")') lastSeenNode     %position
             message="position [satellite] = "//trim(coordinateLabel)
             call Galacticus_Display_Message (message                          ,verbosityWarn)
             write (coordinateLabel,'("[",e12.6,",",e12.6,",",e12.6,"]")') primaryProgenitor%velocity
             message="velocity [primary  ] = "//trim(coordinateLabel)
             call Galacticus_Display_Message (message                          ,verbosityWarn)
             write (coordinateLabel,'("[",e12.6,",",e12.6,",",e12.6,"]")') lastSeenNode     %velocity
             message="velocity [satellite] = "//trim(coordinateLabel)
             call Galacticus_Display_Message (message                          ,verbosityWarn)
             call Galacticus_Display_Unindent('assuming instantaneous merging' ,verbosityWarn)
          else
             ! Create the orbit.
             thisOrbit=Orbit_Construct(lastSeenNode%nodeMass,primaryProgenitor%nodeMass,relativePosition,relativeVelocity)
             ! Construct temporary nodes.
             satelliteNode                => treeNode           (                 )
             hostNode                     => treeNode           (                 )
             satelliteNode%parent         => hostNode
             hostNode     %firstSatellite => satelliteNode
             satelliteBasicComponent      => satelliteNode%basic(autoCreate=.true.)
             hostBasicComponent           => hostNode     %basic(autoCreate=.true.)
             call satelliteBasicComponent%timeSet(lastSeenNode     %nodeTime)
             call      hostBasicComponent%timeSet(primaryProgenitor%nodeTime)
             call satelliteBasicComponent%massSet(lastSeenNode     %nodeMass)
             call      hostBasicComponent%massSet(primaryProgenitor%nodeMass)
             ! Determine the time until merging.
             timeUntilMerging=thisSatelliteMergingTimescales%timeUntilMerging(satelliteNode,thisOrbit)
             ! Clean up.
             call satelliteNode%destroy()
             call hostNode     %destroy()
             deallocate(satelliteNode)
             deallocate(hostNode     )
          end if
       end if
       ! Find the new merging time, and the node with which the merging will occur.
       thisNode          => lastSeenNode%descendent
       timeSubhaloMerges =  timeSubhaloMerges+timeUntilMerging
       do while (associated(thisNode%descendent))
          if (thisNode%descendent%nodeTime > timeSubhaloMerges) then
             nodes(iNode)%mergesWithIndex=thisNode%nodeIndex
             exit
          else
             thisNode => thisNode%descendent
          end if
       end do
       ! Merging time is beyond the end of the tree. Set merging time to infinity.
       if (.not.associated(thisNode%descendent)) then
          timeSubhaloMerges=timeUntilMergingInfinite
          nodes(iNode)%mergesWithIndex=-1
       end if
    else
       call Galacticus_Error_Report('Time_Until_Merging_Subresolution','no descendents found')
    end if
    return
  end subroutine Time_Until_Merging_Subresolution

  function Orbit_Construct(mass1,mass2,position,velocity) result(thisOrbit)
    !% Construct a Keplerian orbit given body masses, positions, and relative velocities.
    use Kepler_Orbits
    use Vectors
    implicit none
    type            (keplerOrbit)                              :: thisOrbit           
    double precision                           , intent(in   ) :: mass1    , mass2    
    double precision             , dimension(3), intent(in   ) :: position , velocity 
    
    call thisOrbit%reset()
    call thisOrbit%massesSet            (       &
         &                               mass1, &
         &                               mass2  &
         &                              )
    call thisOrbit%radiusSet            (                                                   Vector_Magnitude(position))
    call thisOrbit%velocityRadialSet    (                    Dot_Product(velocity,position)/Vector_Magnitude(position))
    call thisOrbit%velocityTangentialSet(Vector_Magnitude(Vector_Product(velocity,position)/Vector_Magnitude(position)))
    return
  end function Orbit_Construct

  subroutine Phase_Space_Position_Realize(time,position,velocity)
    !% Modify relative positions and velocities to account for both any periodicity of the simulated volume, and for Hubble flow.
    use Cosmology_Functions
    use Numerical_Constants_Boolean
    implicit none
    double precision                                       , intent(in   ) :: time                                
    double precision                         , dimension(3), intent(inout) :: position                 , velocity 
    class           (cosmologyFunctionsClass), pointer                     :: cosmologyFunctionsDefault           
    double precision                                                       :: lengthSimulationBox                 
    
    ! Account for periodicity.
    if (defaultImporter%positionsArePeriodic() /= booleanFalse) then
       lengthSimulationBox=defaultImporter%cubeLength(time)
       position=mod(position+0.5d0*lengthSimulationBox,lengthSimulationBox)-0.5d0*lengthSimulationBox
       position=mod(position-0.5d0*lengthSimulationBox,lengthSimulationBox)+0.5d0*lengthSimulationBox
    end if
    ! Account for Hubble flow.
    if (defaultImporter%velocitiesIncludeHubbleFlow() /= booleanTrue) then
       cosmologyFunctionsDefault => cosmologyFunctions()
       velocity=velocity                                                    &
            &  +position                                                    &
            &  *cosmologyFunctionsDefault%hubbleParameterEpochal(time=time)
    end if
    return
  end subroutine Phase_Space_Position_Realize

  subroutine progenitorIteratorDescendentSet(self,node,nodes)
    !% Initialize a progenitor iterator object by storing the index of the taget {\tt node} and finding the location of the first
    !% progenitor (if any).
    implicit none
    class(progenitorIterator)              , intent(inout) :: self  
    class(nodeData          )              , intent(in   ) :: node  
    class(nodeData          ), dimension(:), intent(in   ) :: nodes 
    
    ! Store the index of the target node.
    self%targetIndex     =node%nodeIndex
    ! Assume no progenitors descendants by default.
    self%progenitorsFound=.false.
    ! Find the index of matching nodes in the list sorted by descendent index.
    self%progenitorIndex=Descendent_Node_Sort_Index(node%nodeIndex)
    if (self%progenitorIndex > 0 .and. self%progenitorIndex <= size(nodes)) then
       ! Progenitors may exist, store the location of the first progenitor if found.
       self%progenitorLocation=descendentLocations(self%progenitorIndex)
       if (associated(nodes(self%progenitorLocation)%descendent)) &
            & self%progenitorsFound=(nodes(self%progenitorLocation)%descendent%nodeIndex == node%nodeIndex)
       ! Increment the initial index so that the first call to get the next progenitor can find the first progenitor by
       ! subtracting one from this index.
       self%progenitorIndex=self%progenitorIndex+1
    end if
    return
  end subroutine progenitorIteratorDescendentSet
  
  logical function progenitorIteratorNext(self,nodes)
    !% Move to the next progenitor using a progenitor iterator object, returning true if the next progenitor exists, false if it
    !% does not.
    implicit none
    class(progenitorIterator)              , intent(inout) :: self  
    class(nodeData          ), dimension(:), intent(in   ) :: nodes 
    
    if (self%progenitorsFound) then
       progenitorIteratorNext=.true.
      do while (self%progenitorIndex > 0)
          self%progenitorIndex=self%progenitorIndex-1
          if (self%progenitorIndex <= 0) exit
          self%progenitorLocation=descendentLocations(self%progenitorIndex)
          if (associated(nodes(self%progenitorLocation)%descendent)) exit
       end do
       if (.not.associated(nodes(self%progenitorLocation)%descendent)) then
          progenitorIteratorNext=.false.
          return
       end if
       if (      self%progenitorIndex                          == 0               ) progenitorIteratorNext=.false.
       if (nodes(self%progenitorLocation)%descendent%nodeIndex /= self%targetIndex) progenitorIteratorNext=.false.
    else
       progenitorIteratorNext=.false.
    end if
    return
  end function progenitorIteratorNext
  
  function progenitorIteratorIndex(self,nodes)
    !% Return the node index of the current progenitor in a progenitor iterator object.
    implicit none
    integer(kind=kind_int8    )                              :: progenitorIteratorIndex 
    class  (progenitorIterator)              , intent(in   ) :: self                    
    class  (nodeData          ), dimension(:), intent(in   ) :: nodes                   
    
    progenitorIteratorIndex=nodes(self%progenitorLocation)%nodeIndex
   return
 end function progenitorIteratorIndex
  
 function progenitorIteratorCurrent(self,nodes)
   !% Return a pointer to the current progenitor in a progenitor iterator object.
   implicit none
   class(nodeData          ), pointer                             :: progenitorIteratorCurrent 
   class(progenitorIterator)              , intent(in   )         :: self                      
   class(nodeData          ), dimension(:), intent(in   ), target :: nodes                     
   
   progenitorIteratorCurrent => nodes(self%progenitorLocation)
   return
 end function progenitorIteratorCurrent
  
 logical function progenitorIteratorExist(self)
   !% Return true if progenitors exist, false otherwise.
   implicit none
   class(progenitorIterator), intent(in   ) :: self 
   
   progenitorIteratorExist=self%progenitorsFound
   return
 end function progenitorIteratorExist

end module Merger_Tree_Read
