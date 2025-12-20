!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  !!{
  Implements a merger tree constructor class which constructs trees by reading their definitions from file.
  !!}

  use    :: Cosmology_Functions               , only : cosmologyFunctionsClass
  use    :: Dark_Matter_Halo_Scales           , only : darkMatterHaloScaleClass
  use    :: Dark_Matter_Profile_Scales        , only : darkMatterProfileScaleRadius       , darkMatterProfileScaleRadiusClass
  use    :: Dark_Matter_Profiles_Concentration, only : darkMatterProfileConcentrationClass
  use    :: Dark_Matter_Profiles_DMO          , only : darkMatterProfileDMOClass
  use    :: Galacticus_Nodes                  , only : nodeComponentBasic                 , nodeComponentDarkMatterProfile   , treeNode
  use    :: Halo_Spin_Distributions           , only : haloSpinDistributionClass
  use    :: Kind_Numbers                      , only : kind_int8
  use    :: MPI_Utilities                     , only : mpiCounter
  use    :: Merger_Tree_Read_Importers        , only : mergerTreeImporterClass
  use    :: Merger_Tree_Seeds                 , only : mergerTreeSeedsClass
  use    :: Nodes_Operators                   , only : nodeOperatorClass
  use    :: Numerical_Random_Numbers          , only : randomNumberGeneratorClass
  !$ use :: OMP_Lib                           , only : omp_lock_kind
  use    :: Output_Times                      , only : outputTimes                        , outputTimesClass
  use    :: Satellite_Merging_Timescales      , only : satelliteMergingTimescalesClass
  use    :: Virial_Orbits                     , only : virialOrbitClass

  ! Enumeration of cross-tree event types.
  !![
  <enumeration>
   <name>pushType</name>
   <description>Cross-tree event type enumeration.</description>
   <entry label="branchJump"      />
   <entry label="subhaloPromotion"/>
  </enumeration>
  !!]

  ! Enumeration of node reachability status.
  !![
  <enumeration>
   <name>nodeReachability</name>
   <description>Node reachability status.</description>
   <entry label="unreachable"/>
   <entry label="reachable"  />
  </enumeration>
  !!]

  ! Enumeration of subhalo angular momentum methods.
  !![
  <enumeration>
   <name>subhaloAngularMomentaMethod</name>
   <description>Subhalo angular momentum methods.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="scale"    />
   <entry label="summation"/>
  </enumeration>
  !!]

  !![
  <mergerTreeConstructor name="mergerTreeConstructorRead">
   <description>
    A merger tree constructor class from data imported from a file and processed into a form suitable for \glc\ to evolve. Merger
    trees are inherently complex structures, particularly when the possibility of subhalos are considered. \glc\ is currently designed
    to work with single descendant merger trees, i.e. ones in which the tree structure is entirely defined by specifying which
    \gls{node} a given \gls{node} is physically associated with at a later time. Additionally, \glc\ expects the merger tree file to
    contain information on the host \gls{node}, i.e. the node within which a given node is physically located. In the following, these
    two properties are labeled {\normalfont \ttfamily descendantNode} and {\normalfont \ttfamily hostNode}. \glc\ assumes that nodes
    for which {\normalfont \ttfamily descendantNode}$=${\normalfont \ttfamily hostNode} are isolated halos (i.e. they are their own
    hosts) while other nodes are subhalos (i.e. they are hosted by some other node). An example of a simple tree structure is shown in
    Fig.~\ref{fig:MergerTreeSimple}. The particular structure would be represented by the following list of nodes and node properties
    (a $-1$ indicates that no descendant node exists):
    \begin{center}
    \begin{tabular}{rrr}
    \hline
    {\normalfont \ttfamily node} &amp; {\normalfont \ttfamily descendantNode} &amp; {\normalfont \ttfamily hostNode} \\
    \hline
    1 &amp; -1 &amp; 1 \\
    2 &amp;  1 &amp; 2 \\
    3 &amp;  2 &amp; 3 \\
    4 &amp;  1 &amp; 4 \\
    5 &amp;  4 &amp; 5 \\
    6 &amp; -1 &amp; 1 \\
    7 &amp;  6 &amp; 4 \\
    8 &amp;  7 &amp; 8 \\
    \hline
    \end{tabular}
    \end{center}
    
    \begin{figure}
     \begin{center}
     \includegraphics[width=160mm]{Diagrams/MergerTreeSimple.pdf}
     \end{center}
     \caption{An example of a simple merger tree structure. Colored circles represent nodes in the merger tree. Each node has a unique
       index indicated by the number inside each circle. Black arrows link each node to its descendant node (as specified by the
       {\normalfont \ttfamily descendantNode} property. Where a node is not its own host node it is placed inside its host node.}
     \label{fig:MergerTreeSimple}
    \end{figure}
    
    The following should be noted when constructing merger tree files:
    \begin{itemize}
    \item Note that \glc\ does not require that nodes be placed on a uniform grid of times/redshifts, nor that mass be conserved along
      a branch of the tree. After processing the tree in this way, \glc\ builds additional links which identify the child node of each
      halo and any sibling nodes. These are not required to specify the tree structure but are computationally convenient.
    \item It is acceptable for a node to begin its existence as a subhalo (i.e. to have never had an isolated node progenitor). Such
      nodes will be created as satellites in the merger tree and, providing the selected node components (see
      \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#sec.Components}{here})
      initialize their properties appropriately, will be evolved correctly.
    \item It is acceptable for an isolated node to have progenitors, none of which are a primary progenitor. This can happen if all
      progenitors descend into subhalos in the isolated node. In such cases, \glc\ will create a clone of the isolated node at a very
      slightly earlier time to act as the primary progenitor. This is necessary to allow the tree to be processed correctly, but does
      not affect the evolution of the tree.
    \item \hyperdef{physics}{mergerTreeConstructRead.missingHosts}{} Normally, cases where a node's host node cannot be found in
      the \gls{forest} will cause \glc\ to exit with an error. Setting {\normalfont \ttfamily
      [missingHostsAreFatal]}$=${\normalfont \ttfamily false} will instead circumvent this issue by making any such
      nodes self-hosting (i.e. they become isolated nodes rather than subhalos). Note that this behavior is not a physically
      correct way to treat such cases---it is intended only to allow trees to be processed in cases where the full \gls{forest}
      is not available.
    \item It is acceptable for nodes to jump between branches in a tree, or even to jump between branches in different trees. In the
      latter case, all trees linked by jumping nodes (a so-called ``\gls{forest}'' of connected trees) must be stored as a single
      forest (with multiple root-nodes) in the merger tree file. \glc\ will process this \gls{forest} of trees simultaneously,
      allowing to nodes to move between their branches.
    \item It is acceptable for a subhalo to later become an isolated halo (as can happen due to three-body interactions; see
      \citealt{sales_cosmic_2007}). If {\normalfont \ttfamily [allowSubhaloPromotions]}$=${\normalfont \ttfamily true} then such
      cases will be handled correctly (i.e. the subhalo will be promoted back to being an isolated halo). If the parameter
      {\normalfont \ttfamily [alwaysPromoteMostMassive]}$=${\normalfont \ttfamily true} then the most massive progenitor is treated
      as the primary progenitor, even if that progenitor is a subhalo. Alternatively, if {\normalfont \ttfamily
      [alwaysPromoteMostMassive]}$=${\normalfont \ttfamily false} then a most massive progenitor that is a subhalo is only treated
      as the primary progenitor \emph{if} no isolated progenitors exist (otherwise, the most massive of the isolated progenitors
      is treated as the primary progenitor). If {\normalfont \ttfamily [allowSubhaloPromotions]}$=${\normalfont \ttfamily false}
      then subhalos are not permitted to become isolated halos. In this case, the following logic will be applied to remove all
      such cases from the tree:\\
    
      \noindent\hspace{ 5mm} $\rightarrow$ \parbox[t]{150mm}{For any branch in a tree which at some point is a subhalo:}\\
    
      \noindent\hspace{10mm} $\rightarrow$ \parbox[t]{145mm}{Beginning from the earliest node in the branch that is a subhalo,
        repeatedly step to the next descendant node;}\\
    
      \noindent\hspace{10mm} $\rightarrow$ \parbox[t]{145mm}{If that descendant is \emph{not} a subhalo then:}\\
    
      \noindent\hspace{15mm} $\rightarrow$ \parbox[t]{140mm}{If there is not currently any non-subhalo node which has the present node
        as its descendant then current node is only descendant of a subhalo. Therefore, try to make this node a subhalo, and propose
        the descendant of the host node of the previous node visited in the branch as the new host:}\\
      
      \noindent\hspace{20mm} $\rightarrow$ \parbox[t]{135mm}{If the proposed host exists:}\\
      
      \noindent\hspace{25mm} $\rightarrow$ \parbox[t]{130mm}{If the mass of the current node is less than that of the proposed host:}\\
      
      \noindent\hspace{30mm} $\rightarrow$ \parbox[t]{125mm}{If the proposed hosts exists before the current node, repeatedly step to
        its descendants until one is found which exists at or after the time of the current node. This is the new proposed host.}\\
      
      \noindent\hspace{30mm} $\rightarrow$ \parbox[t]{125mm}{If the proposed host is a subhalo, make it an isolated node.}\\
      
      \noindent\hspace{30mm} $\rightarrow$ \parbox[t]{125mm}{The current node is made a subhalo within the proposed host.}\\
      
      \noindent\hspace{25mm} $\rightarrow$ \parbox[t]{130mm}{Otherwise:}\\
      
      \noindent\hspace{30mm} $\rightarrow$ \parbox[t]{125mm}{The current node remains an isolated node, while the proposed host is
        instead made a subhalo within the current node.}\\
      
      \noindent\hspace{20mm} $\rightarrow$ \parbox[t]{135mm}{Otherwise:}\\
      
      \noindent\hspace{25mm} $\rightarrow$ \parbox[t]{130mm}{The proposed host does not exists, which implies the end of a branch has
        been reached. Therefore, flag the current node as being a subhalo with a host identical to that of the node from which it
        descended.}\\
    \end{itemize}
    
    \textbf{Requirements for \glc\ Input Parameters:} The following requirements must be met for the input parameters to \glc\ when
    using merger trees read from file:
    \begin{itemize}
    \item The cosmological parameters ($\Omega_\mathrm{M}$, $\Omega_\Lambda$, $\Omega_\mathrm{b}$, $H_0$, $\sigma_8$), if defined in
      the file, must be set identically in the \glc\ input file unless you set {\normalfont \ttfamily
        [mismatchIsFatal]}$=${\normalfont \ttfamily false} in which case you'll just be warned about any mismatch;
    \item \glc\ assumes by default that all merger trees exist at the final output time---if this is not the case set {\normalfont
        \ttfamily [allTreesExistAtFinalTime]}$=${\normalfont \ttfamily false}.
    \end{itemize}
    
    \textbf{Dark Matter Scale Radii}: \index{dark matter halo!concentration}\index{dark matter halo!scale radius} If {\normalfont
      \ttfamily [presetScaleRadii]}$=${\normalfont \ttfamily true} and the {\normalfont \ttfamily halfMassRadius}
    dataset is available within the {\normalfont \ttfamily haloTrees} group (see
    \href{https://github.com/galacticusorg/galacticus/wiki/Merger-Tree-File-Format#forest-halos-group}{here}) then the half-mass radii
    of nodes will be used to compute the corresponding scale length of the dark matter halo profile\footnote{The scale radius is found
      by seeking a value which gives the correct half mass radius. It is therefore important that the definition of halo mass
      (specifically the virial overdensity) in \protect\glc\ be the same as was used in computing the input half mass radii.}. This
    requires a dark matter profile scale component which supports setting of the scale length (see
    \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#sec.DarkMatterProfileScale}{here}).
    
    \textbf{Satellite Merger Times:}\index{merger times}\index{satellite!merger times} If {\normalfont \ttfamily
      [presetMergerTimes]}$=${\normalfont \ttfamily true} then merger times for satellites will be computed directly
    from the merger tree data read from file. When a subhalo has an isolated halo as a descendant it is assumed to undergo a merger
    with that isolated halo at that time. Note that this requires a satellite orbit component method which supports setting of merger
    times (e.g. {\normalfont \ttfamily [componentSatellite]}$=${\normalfont \ttfamily preset}).
    
    \textbf{Dark Matter Halo Angular Momenta:}\index{dark matter halo!angular momentum} If {\normalfont \ttfamily
      [presetAngularMomenta]}$=${\normalfont \ttfamily true} and the {\normalfont \ttfamily angularMomentum} dataset is available
    within the {\normalfont \ttfamily haloTrees} group (see
    \href{https://github.com/galacticusorg/galacticus/wiki/Merger-Tree-File-Format#forest-halos-group}{here}) then the angular momenta
    of nodes will be computed and set. This requires a dark matter halo spin component which supports setting of the angular momentum (see
    \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#sec.DarkMatterHaloSpinComponent}{here}).
   </description>
  </mergerTreeConstructor>
  !!]
  type, extends(mergerTreeConstructorClass) :: mergerTreeConstructorRead
     !!{
     A class implementing merger tree construction by reading trees from file.
     !!}
     private
     class           (cosmologyFunctionsClass                   ), pointer                   :: cosmologyFunctions_                   => null()
     class           (mergerTreeImporterClass                   ), pointer                   :: mergerTreeImporter_                   => null()
     class           (darkMatterProfileConcentrationClass       ), pointer                   :: darkMatterProfileConcentration_       => null()
     class           (satelliteMergingTimescalesClass           ), pointer                   :: satelliteMergingTimescales_           => null()
     class           (darkMatterHaloScaleClass                  ), pointer                   :: darkMatterHaloScale_                  => null()
     class           (darkMatterProfileDMOClass                 ), pointer                   :: darkMatterProfileDMO_                 => null()
     class           (haloSpinDistributionClass                 ), pointer                   :: haloSpinDistribution_                 => null()
     class           (virialOrbitClass                          ), pointer                   :: virialOrbit_                          => null()
     class           (outputTimesClass                          ), pointer                   :: outputTimes_                          => null()
     class           (darkMatterProfileScaleRadiusClass         ), pointer                   :: darkMatterProfileScaleRadius_         => null()
     class           (randomNumberGeneratorClass                ), pointer                   :: randomNumberGenerator_                => null()
     class           (nodeOperatorClass                         ), pointer                   :: nodeOperator_                         => null()
     class           (mergerTreeSeedsClass                      ), pointer                   :: mergerTreeSeeds_                      => null()
     integer                                                                                 :: fileCurrent
     type            (varying_string                            ), allocatable, dimension(:) :: fileNames                                       , presetNamedReals                    , &
          &                                                                                     presetNamedIntegers
     integer                                                     , allocatable, dimension(:) :: indexNamedReals                                 , indexNamedIntegers
     logical                                                                                 :: importerOpen                          =  .false.
     integer         (kind_int8                                 )                            :: beginAt
     logical                                                                                 :: foundBeginAt                          =  .false.
     double precision                                                                        :: treeWeightCurrent
     logical                                                                                 :: allowBranchJumps
     logical                                                                                 :: allowSubhaloPromotions                          , alwaysPromoteMostMassive
     integer         (c_size_t                                  )                            :: forestSizeMaximum                               , treeNumberOffset
     logical                                                                                 :: presetMergerTimes
     logical                                                                                 :: presetMergerNodes
     logical                                                                                 :: presetSubhaloMasses
     logical                                                                                 :: presetSubhaloIndices
     logical                                                                                 :: presetPositions
     logical                                                                                 :: presetScaleRadii                                , scaleRadiiFailureIsFatal
     double precision                                                                        :: presetScaleRadiiMinimumMass                     , presetScaleRadiiConcentrationMinimum, &
          &                                                                                     presetScaleRadiiConcentrationMaximum
     logical                                                                                 :: presetAngularMomenta                            , presetAngularMomenta3D              , &
          &                                                                                     presetUnphysicalAngularMomenta
     logical                                                                                 :: presetOrbits                                    , presetOrbitsAssertAllSet            , &
          &                                                                                     presetOrbitsBoundOnly                           , presetOrbitsSetAll
     type            (enumerationSubhaloAngularMomentaMethodType)                            :: subhaloAngularMomentaMethod
     logical                                                                                 :: missingHostsAreFatal
     logical                                                                                 :: treeIndexToRootNodeIndex
     integer         (c_size_t                                  )                            :: outputTimesCount
     double precision                                                                        :: outputTimeSnapTolerance
     double precision                                            , allocatable, dimension(:) :: outputTimes
     integer         (c_size_t                                  ), allocatable, dimension(:) :: descendantLocations                              , nodeLocations
     integer         (kind_int8                                 ), allocatable, dimension(:) :: descendantIndicesSorted                          , nodeIndicesSorted
     !$ integer      (omp_lock_kind                             )                            :: splitForestLock
     integer                                                                                 :: splitForestActiveForest
     integer         (c_size_t                                  )                            :: splitForestNextTree                              , splitForestUniqueID
     integer         (c_size_t                                  ), allocatable, dimension(:) :: splitForestTreeSize                              , splitForestTreeStart               , &
          &                                                                                     splitForestMapIndex
     integer         (kind_int8                                 ), allocatable, dimension(:) :: splitForestPushTo                                , splitForestPullFrom
     type            (enumerationPushTypeType                   ), allocatable, dimension(:) :: splitForestPushType
     double precision                                            , allocatable, dimension(:) :: splitForestPushTime
     logical                                                     , allocatable, dimension(:) :: splitForestIsPrimary                            , splitForestPushDone                , &
          &                                                                                     splitForestPullDone
     logical                                                                                 :: warningNestedHierarchyIssued
     logical                                                                                 :: warningSplitForestNestedHierarchyIssued
     real                                                        , allocatable, dimension(:) :: timingTimes
     type            (varying_string                            ), allocatable, dimension(:) :: timingLabels
   contains
     !![
     <methods>
       <method description="Ensure that any node which was once a subhalo remains a subhalo." method="enforceSubhaloStatus" />
       <method description="Scan for cases where a subhalo stops being a subhalo and so must be promoted." method="scanForSubhaloPromotions" />
       <method description="Create a sorted list of node indices with an index into the original array." method="createNodeIndices" />
       <method description="Return the location in the original array of the given {\normalfont \ttfamily nodeIndex}." method="nodeLocation" />
       <method description="Return the sort index of the given {\normalfont \ttfamily descendantIndex}." method="descendantNodeSortIndex" />
       <method description="Destroy the sorted list of node indices." method="destroyNodeIndices" />
       <method description="Builds pointers from each node to its descendant node." method="buildDescendantPointers" />
       <method description="Create parent pointer links between isolated nodes and assign times and masses to those nodes." method="buildIsolatedParentPointers" />
       <method description="Assign named properties to nodes." method="assignNamedProperties" />
       <method description="Assign scale radii to nodes." method="assignScaleRadii" />
       <method description="Scan for and record mergers between nodes." method="scanForMergers" />
       <method description="Set the virial orbit in a node." method="setOrbit" />
       <method description="Assign angular momenta to nodes." method="assignAngularMomenta" />
       <method description="Assign pointers to merge targets." method="assignMergers" />
       <method description="Search for subhalos which move between branches/trees." method="scanForBranchJumps" />
       <method description="Build and attached bound mass histories to subhalos." method="buildSubhaloMassHistories" />
       <method description="Compute the additional time until merging after a subhalo is lost from the tree (presumably due to limited resolution)." method="timeUntilMergingSubresolution" />
       <method description="Modify relative positions and velocities to account for both any periodicity of the simulated volume, and for Hubble flow." method="phaseSpacePositionRealize" />
       <method description="Record timing data." method="timingRecord" />
       <method description="Report on time taken in various steps of processing merger trees read from file." method="timingReport" />
       <method description="Find initial root node affinities for all nodes." method="rootNodeAffinitiesInitial" />
       <method description="Return true if the given node is on the current ``push-to'' list of nodes for split forests." method="isOnPushList" />
       <method description="Return true if the given node is on the current ``pull-from'' list of nodes for split forests." method="isOnPullList" />
       <method description="Return the index of the given node in the ``push-to'' list of nodes for split forests." method="pushListIndex" />
       <method description="Return the index of the given node in the ``pull-from'' list of nodes for split forests." method="pullListIndex" />
       <method description="Return the number of the given node in the ``pull-from'' list of nodes for split forests." method="pullListCount" />
       <method description="Assign events to nodes if they jump between trees in a forest." method="assignSplitForestEvents" />
       <method description="Returns true if {\normalfont \ttfamily node} undergoes a subhalo-subhalo merger." method="isSubhaloSubhaloMerger" />
       <method description="Create an array of standard nodes and associated structures." method="createNodeArray" />
     </methods>
     !!]
     final     ::                                  readDestructor
     procedure :: construct                     => readConstruct
     procedure :: enforceSubhaloStatus          => readEnforceSubhaloStatus
     procedure :: scanForSubhaloPromotions      => readScanForSubhaloPromotions
     procedure :: createNodeIndices             => readCreateNodeIndices
     procedure :: nodeLocation                  => readNodeLocation
     procedure :: descendantNodeSortIndex       => readDescendantNodeSortIndex
     procedure :: destroyNodeIndices            => readDestroyNodeIndices
     procedure :: buildDescendantPointers       => readBuildDescendantPointers
     procedure :: buildIsolatedParentPointers   => readBuildIsolatedParentPointers
     procedure :: assignNamedProperties         => readAssignNamedProperties
     procedure :: assignScaleRadii              => readAssignScaleRadii
     procedure :: scanForMergers                => readScanForMergers
     procedure :: setOrbit                      => readSetOrbit
     procedure :: assignAngularMomenta          => readAssignAngularMomenta
     procedure :: assignMergers                 => readAssignMergers
     procedure :: scanForBranchJumps            => readScanForBranchJumps
     procedure :: buildSubhaloMassHistories     => readBuildSubhaloMassHistories
     procedure :: timeUntilMergingSubresolution => readTimeUntilMergingSubresolution
     procedure :: phaseSpacePositionRealize     => readPhaseSpacePositionRealize
     procedure :: timingRecord                  => readTimingRecord
     procedure :: timingReport                  => readTimingReport
     procedure :: rootNodeAffinitiesInitial     => readRootNodeAffinitiesInitial
     procedure :: isOnPushList                  => readIsOnPushList
     procedure :: isOnPullList                  => readIsOnPullList
     procedure :: pushListIndex                 => readPushListIndex
     procedure :: pullListIndex                 => readPullListIndex
     procedure :: pullListCount                 => readPullListCount
     procedure :: assignSplitForestEvents       => readAssignSplitForestEvents
     procedure :: isSubhaloSubhaloMerger        => readIsSubhaloSubhaloMerger
     procedure :: createNodeArray               => readCreateNodeArray
  end type mergerTreeConstructorRead

  interface mergerTreeConstructorRead
     !!{
     Constructors for the \refClass{mergerTreeConstructorRead} merger tree constructor class.
     !!}
     module procedure readConstructorParameters
     module procedure readConstructorInternal
  end interface mergerTreeConstructorRead

  ! Iterator object for iterating over progenitor nodes.
  type :: progenitorIterator
     integer(c_size_t                 )          :: progenitorLocation
     integer(kind_int8                )          :: progenitorIndex             , targetIndex
     logical                                     :: progenitorsFound
     type   (mergerTreeConstructorRead), pointer :: constructor        => null()
   contains
     !![
     <methods>
       <method description="Set the target descendant node and initialize the iterator." method="descendantSet" />
       <method description="Move to the next progenitor. Returns true if the next progenitor exists, false otherwise." method="next" />
       <method description="Return the index of the current progenitor." method="index" />
       <method description="Return a pointer to the current progenitor." method="current" />
       <method description="Return true if any progenitors exist, false otherwise." method="exist" />
     </methods>
     !!]
     procedure :: descendantSet => progenitorIteratorDescendantSet
     procedure :: next          => progenitorIteratorNext
     procedure :: index         => progenitorIteratorIndex
     procedure :: current       => progenitorIteratorCurrent
     procedure :: exist         => progenitorIteratorExist
  end type progenitorIterator

  ! Variables used in root-finding.
  class           (nodeComponentDarkMatterProfile     ), pointer :: darkMatterProfile_
  class           (nodeComponentBasic                 ), pointer :: basic_
  type            (treeNode                           ), pointer :: node_
  double precision                                               :: radiusHalfMass_
  class           (mergerTreeConstructorRead          ), pointer :: self_
  !$omp threadprivate(darkMatterProfile_,basic_,node_,radiusHalfMass_,self_)

  ! Counter used to assign unique IDs to split forests.
  type(mpiCounter) :: splitForestUniqueID

  ! Fractional offset in time used for cloned nodes.
  double precision :: fractionOffsetTimeClones=1.0d-9

contains

  function readConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeConstructorRead} merger tree constructor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeConstructorRead          )                            :: self
    type            (inputParameters                    ), intent(inout)             :: parameters
    class           (cosmologyFunctionsClass            ), pointer                   :: cosmologyFunctions_
    class           (mergerTreeImporterClass            ), pointer                   :: mergerTreeImporter_
    class           (darkMatterProfileConcentrationClass), pointer                   :: darkMatterProfileConcentration_
    class           (satelliteMergingTimescalesClass    ), pointer                   :: satelliteMergingTimescales_
    class           (darkMatterHaloScaleClass           ), pointer                   :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass          ), pointer                   :: darkMatterProfileDMO_
    class           (haloSpinDistributionClass          ), pointer                   :: haloSpinDistribution_
    class           (virialOrbitClass                   ), pointer                   :: virialOrbit_
    class           (outputTimesClass                   ), pointer                   :: outputTimes_
    class           (darkMatterProfileScaleRadiusClass  ), pointer                   :: darkMatterProfileScaleRadius_
    class           (randomNumberGeneratorClass         ), pointer                   :: randomNumberGenerator_
    class           (nodeOperatorClass                  ), pointer                   :: nodeOperator_
    class           (mergerTreeSeedsClass               ), pointer                   :: mergerTreeSeeds_
    type            (varying_string                     ), allocatable, dimension(:) :: fileNames                           , presetNamedReals                    , &
         &                                                                              presetNamedIntegers
    integer                                                                          :: fileCount
    integer         (c_size_t                           )                            :: forestSizeMaximum
    integer         (kind_int8                          )                            :: beginAt
    logical                                                                          :: missingHostsAreFatal                , presetMergerTimes                   , &
         &                                                                              presetMergerNodes                   , presetSubhaloMasses                 , &
         &                                                                              presetSubhaloIndices                , presetPositions                     , &
         &                                                                              presetScaleRadii                    , scaleRadiiFailureIsFatal            , &
         &                                                                              presetUnphysicalAngularMomenta      , presetAngularMomenta                , &
         &                                                                              presetAngularMomenta3D              , presetOrbits                        , &
         &                                                                              presetOrbitsSetAll                  , presetOrbitsAssertAllSet            , &
         &                                                                              presetOrbitsBoundOnly               , allowSubhaloPromotions              , &
         &                                                                              treeIndexToRootNodeIndex            , allowBranchJumps                    , &
         &                                                                              alwaysPromoteMostMassive
    type            (varying_string                     )                            :: subhaloAngularMomentaMethod
    double precision                                                                 :: presetScaleRadiiConcentrationMinimum, presetScaleRadiiConcentrationMaximum, &
         &                                                                              presetScaleRadiiMinimumMass         , outputTimeSnapTolerance
    !$GLC attributes initialized :: self

    fileCount=parameters%count('fileNames')
    allocate(fileNames(fileCount))
    !![
    <inputParameter>
      <name>fileNames</name>
      <description>The name of the file(s) from which merger tree data should be read when using the {\normalfont \ttfamily [mergerTreeConstruct]}$=${\normalfont \ttfamily read} tree construction method.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>forestSizeMaximum</name>
      <defaultValue>0_c_size_t</defaultValue>
      <description>The maximum number of nodes allowed in a forest before it will be broken up into trees and processed individually. A value of 0 implies that forests should never be split.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetMergerTimes</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether merging times for subhalos should be preset when reading merger trees from a file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetMergerNodes</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether the target nodes for mergers should be preset (i.e. determined from descendant nodes). If they are not, merging will be with each satellite's host node.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetSubhaloMasses</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether subhalo mass should be preset when reading merger trees from a file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>subhaloAngularMomentaMethod</name>
      <defaultValue>var_str('summation')</defaultValue>
      <description>Specifies how to account for subhalo angular momentum when adding subhalo mass to host halo mass.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetSubhaloIndices</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether subhalo indices should be preset when reading merger trees from a file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetPositions</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether node positions should be preset when reading merger trees from a file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetScaleRadii</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether node scale radii should be preset when reading merger trees from a file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleRadiiFailureIsFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether failure to set a node scale radii should be regarded as a fatal error. (If not, a fallback method to set scale radius is used in such cases.)</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetScaleRadiiConcentrationMinimum</name>
      <defaultValue>3.0d0</defaultValue>
      <description>The lowest concentration ($c=r_\mathrm{vir}/r_\mathrm{s}$) allowed when setting scale radii, $r_\mathrm{s}$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetScaleRadiiConcentrationMaximum</name>
      <defaultValue>60.0d0</defaultValue>
      <description>The largest concentration ($c=r_\mathrm{vir}/r_\mathrm{s}$) allowed when setting scale radii, $r_\mathrm{s}$.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetScaleRadiiMinimumMass</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum halo mass for which scale radii should be preset (if {\normalfont \ttfamily [presetScaleRadii]}$=${\normalfont \ttfamily true}).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetUnphysicalAngularMomenta</name>
      <defaultValue>.false.</defaultValue>
      <description>When reading merger trees from file and presetting halo angular momenta, detect unphysical (&lt;=0) angular momenta and preset them using the selected halo spin method.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetAngularMomenta</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether node angular momenta should be preset when reading merger trees from a file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetAngularMomenta3D</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether node 3-D angular momenta vectors should be preset when reading merger trees from a file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetOrbits</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether node orbits should be preset when reading merger trees from a file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetOrbitsSetAll</name>
      <defaultValue>.true.</defaultValue>
      <description>Forces all orbits to be set. If the computed orbit does not cross the virial radius, then select one at random instead.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetOrbitsAssertAllSet</name>
      <defaultValue>.true.</defaultValue>
      <description>Asserts that all virial orbits must be preset. If any can not be set, \glc\ will stop.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>presetOrbitsBoundOnly</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether only bound node orbits should be set.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>beginAt</name>
      <defaultValue>-1_kind_int8</defaultValue>
      <description>Specifies the index of the tree to begin at. (Use -1 to always begin with the first tree.)</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputTimeSnapTolerance</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The relative tolerance required to ``snap'' a node time to the closest output time.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>missingHostsAreFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether nodes with missing host nodes should be considered to be fatal---see \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Physics.pdf\#physics.mergerTreeConstructRead.missingHosts}{here}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>treeIndexToRootNodeIndex</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether tree indices should always be set to the index of their root node.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>allowBranchJumps</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether nodes are allowed to jump between branches.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>allowSubhaloPromotions</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether subhalos are permitted to be promoted to being isolated halos.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>alwaysPromoteMostMassive</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, the most massive progenitor is always promoted to be the primary progenitor \emph{even if} it is a subhalo. Otherwise, isolated progenitors are given priority over subhalo progenitors, even if they are less massive.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    allocate(presetNamedReals   (parameters%count('presetNamedReals'   ,zeroIfNotPresent=.true.)))
    if (size(presetNamedReals   ) > 0) then
       !![
       <inputParameter>
         <name>presetNamedReals</name>
         <description>Names of real datasets to be additionally read and stored in the nodes of the merger tree when using the {\normalfont \ttfamily [mergerTreeConstruct]}$=${\normalfont \ttfamily read} tree construction method.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    allocate(presetNamedIntegers(parameters%count('presetNamedIntegers',zeroIfNotPresent=.true.)))
    if (size(presetNamedIntegers) > 0) then
       !![
       <inputParameter>
         <name>presetNamedIntegers</name>
         <description>Names of integer datasets to be additionally read and stored in the nodes of the merger tree when using the {\normalfont \ttfamily [mergerTreeConstruct]}$=${\normalfont \ttfamily read} tree construction method.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"                                                         />
    <objectBuilder class="mergerTreeImporter"             name="mergerTreeImporter_"             source="parameters"                                                         />
    <objectBuilder class="darkMatterHaloScale"            name="darkMatterHaloScale_"            source="parameters"                                                         />
    <objectBuilder class="darkMatterProfileDMO"           name="darkMatterProfileDMO_"           source="parameters"                                                         />
    <objectBuilder class="darkMatterProfileConcentration" name="darkMatterProfileConcentration_" source="parameters"                                                         />
    <objectBuilder class="haloSpinDistribution"           name="haloSpinDistribution_"           source="parameters"                                                         />
    <objectBuilder class="virialOrbit"                    name="virialOrbit_"                    source="parameters"                                                         />
    <objectBuilder class="outputTimes"                    name="outputTimes_"                    source="parameters"                                                         />
    <objectBuilder class="darkMatterProfileScaleRadius"   name="darkMatterProfileScaleRadius_"   source="parameters"                                                         />
    <objectBuilder class="mergerTreeSeeds"                name="mergerTreeSeeds_"                source="parameters"                                                         />
    <objectBuilder class="nodeOperator"                   name="nodeOperator_"                   source="parameters" parameterName="nodeOperatorTreeInitializor"              >
     <default>
      <nodeOperatorTreeInitializor value="null"/>
     </default>
    </objectBuilder>
    <objectBuilder class="randomNumberGenerator"          name="randomNumberGenerator_"          source="parameters"                                                         />
    <objectBuilder class="satelliteMergingTimescales"     name="satelliteMergingTimescales_"     source="parameters" parameterName="satelliteMergingTimescalesSubresolution"  >
     <default>
      <satelliteMergingTimescalesSubresolution value="zero"/>
     </default>
    </objectBuilder>
    !!]
    self=mergerTreeConstructorRead(                                                                                                                 &
         &                                                                           fileNames                                                    , &
         &                                                                           outputTimeSnapTolerance                                      , &
         &                                                                           forestSizeMaximum                                            , &
         &                                                                           beginAt                                                      , &
         &                                                                           missingHostsAreFatal                                         , &
         &                                                                           treeIndexToRootNodeIndex                                     , &
         &                         enumerationSubhaloAngularMomentaMethodEncode(char(subhaloAngularMomentaMethod         ),includesPrefix=.false.), &
         &                                                                           allowBranchJumps                                             , &
         &                                                                           allowSubhaloPromotions                                       , &
         &                                                                           alwaysPromoteMostMassive                                     , &
         &                                                                           presetMergerTimes                                            , &
         &                                                                           presetMergerNodes                                            , &
         &                                                                           presetSubhaloMasses                                          , &
         &                                                                           presetSubhaloIndices                                         , &
         &                                                                           presetPositions                                              , &
         &                                                                           presetScaleRadii                                             , &
         &                                                                           presetScaleRadiiConcentrationMinimum                         , &
         &                                                                           presetScaleRadiiConcentrationMaximum                         , &
         &                                                                           presetScaleRadiiMinimumMass                                  , &
         &                                                                           scaleRadiiFailureIsFatal                                     , &
         &                                                                           presetUnphysicalAngularMomenta                               , &
         &                                                                           presetAngularMomenta                                         , &
         &                                                                           presetAngularMomenta3D                                       , &
         &                                                                           presetOrbits                                                 , &
         &                                                                           presetOrbitsSetAll                                           , &
         &                                                                           presetOrbitsAssertAllSet                                     , &
         &                                                                           presetOrbitsBoundOnly                                        , &
         &                                                                           presetNamedReals                                             , &
         &                                                                           presetNamedIntegers                                          , &
         &                                                                           cosmologyFunctions_                                          , &
         &                                                                           mergerTreeImporter_                                          , &
         &                                                                           mergerTreeSeeds_                                             , &
         &                                                                           darkMatterHaloScale_                                         , &
         &                                                                           darkMatterProfileDMO_                                        , &
         &                                                                           darkMatterProfileConcentration_                              , &
         &                                                                           haloSpinDistribution_                                        , &
         &                                                                           satelliteMergingTimescales_                                  , &
         &                                                                           virialOrbit_                                                 , &
         &                                                                           outputTimes_                                                 , &
         &                                                                           darkMatterProfileScaleRadius_                                , &
         &                                                                           nodeOperator_                                                , &
         &                                                                           randomNumberGenerator_                                         &
         &                        )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="mergerTreeImporter_"            />
    <objectDestructor name="darkMatterHaloScale_"           />
    <objectDestructor name="darkMatterProfileDMO_"          />
    <objectDestructor name="darkMatterProfileConcentration_"/>
    <objectDestructor name="mergerTreeSeeds_"               />
    <objectDestructor name="haloSpinDistribution_"          />
    <objectDestructor name="virialOrbit_"                   />
    <objectDestructor name="outputTimes_"                   />
    <objectDestructor name="darkMatterProfileScaleRadius_"  />
    <objectDestructor name="nodeOperator_"                  />
    <objectDestructor name="satelliteMergingTimescales_"    />
    <objectDestructor name="randomNumberGenerator_"         />
    !!]
    return
  end function readConstructorParameters

  function readConstructorInternal(fileNames,outputTimeSnapTolerance,forestSizeMaximum,beginAt,missingHostsAreFatal,treeIndexToRootNodeIndex,subhaloAngularMomentaMethod,allowBranchJumps,allowSubhaloPromotions,alwaysPromoteMostMassive,presetMergerTimes,presetMergerNodes,presetSubhaloMasses,presetSubhaloIndices,presetPositions,presetScaleRadii,presetScaleRadiiConcentrationMinimum,presetScaleRadiiConcentrationMaximum,presetScaleRadiiMinimumMass,scaleRadiiFailureIsFatal,presetUnphysicalAngularMomenta,presetAngularMomenta,presetAngularMomenta3D,presetOrbits,presetOrbitsSetAll,presetOrbitsAssertAllSet,presetOrbitsBoundOnly,presetNamedReals,presetNamedIntegers,cosmologyFunctions_,mergerTreeImporter_,mergerTreeSeeds_,darkMatterHaloScale_,darkMatterProfileDMO_,darkMatterProfileConcentration_,haloSpinDistribution_,satelliteMergingTimescales_,virialOrbit_,outputTimes_,darkMatterProfileScaleRadius_,nodeOperator_,randomNumberGenerator_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeConstructorRead} merger tree constructor class.
    !!}
    use    :: Display                    , only : displayMagenta, displayReset
    use    :: Error                      , only : Error_Report  , Warn
    use    :: Numerical_Constants_Boolean, only : booleanFalse  , booleanTrue
    !$ use :: OMP_Lib                    , only : OMP_Init_Lock
    implicit none
    type            (mergerTreeConstructorRead                 )                              :: self
    class           (cosmologyFunctionsClass                   ), intent(in   ), target       :: cosmologyFunctions_
    class           (mergerTreeImporterClass                   ), intent(in   ), target       :: mergerTreeImporter_
    class           (mergerTreeSeedsClass                      ), intent(in   ), target       :: mergerTreeSeeds_
    class           (darkMatterHaloScaleClass                  ), intent(in   ), target       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass                 ), intent(in   ), target       :: darkMatterProfileDMO_
    class           (darkMatterProfileConcentrationClass       ), intent(in   ), target       :: darkMatterProfileConcentration_
    class           (haloSpinDistributionClass                 ), intent(in   ), target       :: haloSpinDistribution_
    class           (satelliteMergingTimescalesClass           ), intent(in   ), target       :: satelliteMergingTimescales_
    class           (virialOrbitClass                          ), intent(in   ), target       :: virialOrbit_
    class           (outputTimesClass                          ), intent(in   ), target       :: outputTimes_
    class           (darkMatterProfileScaleRadiusClass         ), intent(in   ), target       :: darkMatterProfileScaleRadius_
    class           (randomNumberGeneratorClass                ), intent(in   ), target       :: randomNumberGenerator_
    class           (nodeOperatorClass                         ), intent(in   ), target       :: nodeOperator_
    type            (varying_string                            ), intent(in   ), dimension(:) :: fileNames                           , presetNamedReals                    , &
         &                                                                                       presetNamedIntegers
    integer         (c_size_t                                  ), intent(in   )               :: forestSizeMaximum
    integer         (kind_int8                                 ), intent(in   )               :: beginAt
    logical                                                     , intent(in   )               :: missingHostsAreFatal                , presetMergerTimes                   , &
         &                                                                                       presetMergerNodes                   , presetSubhaloMasses                 , &
         &                                                                                       presetSubhaloIndices                , presetPositions                     , &
         &                                                                                       presetScaleRadii                    , scaleRadiiFailureIsFatal            , &
         &                                                                                       presetUnphysicalAngularMomenta      , presetAngularMomenta                , &
         &                                                                                       presetAngularMomenta3D              , presetOrbits                        , &
         &                                                                                       presetOrbitsSetAll                  , presetOrbitsAssertAllSet            , &
         &                                                                                       presetOrbitsBoundOnly               , allowSubhaloPromotions              , &
         &                                                                                       treeIndexToRootNodeIndex            , allowBranchJumps                    , &
         &                                                                                       alwaysPromoteMostMassive
    type            (enumerationSubhaloAngularMomentaMethodType), intent(in   )               :: subhaloAngularMomentaMethod
    double precision                                            , intent(in   )               :: presetScaleRadiiConcentrationMinimum, presetScaleRadiiConcentrationMaximum, &
         &                                                                                       presetScaleRadiiMinimumMass         , outputTimeSnapTolerance
    integer         (c_size_t                                  )                              :: iOutput                             , i
    type            (varying_string                            )                              :: message
    !![
    <constructorAssign variables="fileNames, outputTimeSnapTolerance, forestSizeMaximum, beginAt, missingHostsAreFatal, treeIndexToRootNodeIndex, subhaloAngularMomentaMethod, allowBranchJumps, allowSubhaloPromotions, alwaysPromoteMostMassive, presetMergerTimes, presetMergerNodes, presetSubhaloMasses, presetSubhaloIndices, presetPositions, presetScaleRadii,  presetScaleRadiiConcentrationMinimum, presetScaleRadiiConcentrationMaximum, presetScaleRadiiMinimumMass, scaleRadiiFailureIsFatal, presetUnphysicalAngularMomenta, presetAngularMomenta, presetAngularMomenta3D, presetOrbits, presetOrbitsSetAll, presetOrbitsAssertAllSet, presetOrbitsBoundOnly, presetNamedReals, presetNamedIntegers, *cosmologyFunctions_, *mergerTreeImporter_, *darkMatterHaloScale_, *darkMatterProfileDMO_, *darkMatterProfileConcentration_, *haloSpinDistribution_, *satelliteMergingTimescales_, *virialOrbit_, *outputTimes_, *darkMatterProfileScaleRadius_, *nodeOperator_, *randomNumberGenerator_, *mergerTreeSeeds_"/>
    !!]

    ! Initialize statuses.
    self%warningNestedHierarchyIssued           =.false.
    self%warningSplitForestNestedHierarchyIssued=.false.
    ! Set initial state indicating if the first tree to process has been found.
    self%foundBeginAt                           =self%beginAt == -1_kind_int8
    ! Initialize split forests counter.
    splitForestUniqueID=mpiCounter()
    ! Get array of output times.
    self%outputTimesCount=self%outputTimes_%count()
    allocate(self%outputTimes(self%outputTimesCount))
    do iOutput=1,self%outputTimesCount
       self%outputTimes(iOutput)=self%outputTimes_%time(iOutput)
    end do
    ! Open the file.
    self%treeNumberOffset=0_c_size_t
    self%fileCurrent     =1
    self%importerOpen    =.true.
    call self%mergerTreeImporter_%open(self%fileNames(self%fileCurrent))
    ! Validate input parameters.
    if (self%presetMergerNodes.and..not.self%presetMergerTimes) then
       message="presetting of merger target nodes requires that merger times also be preset;"//char(10)
       message=message//" try setting [presetMergerTimes]=true."//char(10)
       if (self%mergerTreeImporter_%treesHaveSubhalos() /= booleanTrue) then
          message=message//" Note: presetting merger target nodes and merger times is usually only a good idea when subhalo information is present in the merger trees"
          if (self%mergerTreeImporter_%treesHaveSubhalos() == booleanFalse) then
             message=message//" (which the current trees do not)"
          else
             message=message//" (subhalo presence in the current trees was not specified)"
          end if
       end if
       call Error_Report(message//{introspection:location})
    end if
    ! Warn about lack of branch jumps and subhalo promotions if split forests are being used.
    if (self%forestSizeMaximum > 0_c_size_t .and. .not.(self%allowBranchJumps .and. self%allowSubhaloPromotions)) then
       message=displayMagenta()//'WARNING:'//displayReset()//' large forests may be split for processing but '
       if (                                     .not.self%allowBranchJumps) message=message//' branch jumps'
       if (.not.self%allowSubhaloPromotions.and..not.self%allowBranchJumps) message=message//' and'
       if (.not.self%allowSubhaloPromotions                               ) message=message//' subhalo promotions'
       message=message//' are not allowed - this can result in inconsistent treatment between split and unsplit forest processing'
       call Warn(message)
    end if
    ! Warn if subhalo promotions are allowed, but branch jumps are not.
    if (self%allowSubhaloPromotions.and..not.self%allowBranchJumps) then
       message=displayMagenta()//'WARNING:'//displayReset()//' allowing subhalo promotions while not allowing branch jumps can lead to deadlocking of trees.'//char(10)
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
       message=message//' ---   -----'//{introspection:location}//char(10)
       call Warn(message)
    end if
    ! Warn if subhalo promotions are allowed, but branch jumps are not.
    if (self%allowBranchJumps.and..not.self%presetMergerTimes) then
       message=displayMagenta()//'WARNING:'//displayReset()//' allowing branch jumps while not presetting merger times can lead to tree deadlock if merging occurs prior to the jumped-to node time and before an output time which blocks the jumped-to node''s child.'//char(10)
       message=message//'For example, "a" in the following tree (in which "<==" indicates subhalo host) jumps from "3" to "1", but cannot merge with "1" since "1" still has a child, and that child, "2", cannot reach "1" because it is blocked by the output time, resulting in a deadlock of the tree:'//char(10)//char(10)
       message=message//' ---           ---'//char(10)
       message=message//' |1|<==========|a|'//char(10)
       message=message//' ---           ---'//char(10)
       message=message//'  ^             | '//char(10)
       message=message//'  |             | '//char(10)
       message=message//' ~~~output time~~~'//char(10)
       message=message//'  |             | '//char(10)
       message=message//' ~~~merge  time~~~'//char(10)
       message=message//'  |             | '//char(10)
       message=message//' ---     ---   ---'//char(10)
       message=message//' |2|     |3|<==|a|'//char(10)
       message=message//' ---     ---   ---'//{introspection:location}//char(10)
       call Warn(message)
    end if
    ! Disallow always promoting the most massive progenitor if subhalo promotions are not allowed.
    if (self%alwaysPromoteMostMassive.and..not.self%allowSubhaloPromotions) &
         & call Error_Report('[alwaysPromoteMostMassive]=true requires [allowSubhaloPromotions]=true - in order to always promote the most massive progenitor it must be allowed to promote subhalos'//{introspection:location})
    ! Perform sanity checks if subhalos are not included.
    if (self%mergerTreeImporter_%treesHaveSubhalos() == booleanFalse) then
       if (self%presetMergerTimes   ) call Error_Report('cannot preset merger times as no subhalos are present; try setting [presetMergerTimes]=false'      //{introspection:location})
       if (self%presetMergerNodes   ) call Error_Report('cannot preset merger nodes as no subhalos are present; try setting [presetMergerNodes]=false'      //{introspection:location})
       if (self%presetSubhaloMasses ) call Error_Report('cannot preset subhalo masses as no subhalos are present; try setting [presetSubhaloMasses]=false'  //{introspection:location})
       if (self%presetSubhaloIndices) call Error_Report('cannot preset subhalo indices as no subhalos are present; try setting [presetSubhaloIndices]=false'//{introspection:location})
    end if
    ! Determine if subhalo masses have been included in halo masses.
    if (self%mergerTreeImporter_%treesAreSelfContained() == booleanFalse) call Error_Report('only self-contained trees are supported'//{introspection:location})
    ! Check that position information is present if required.
    if     (                                                                                      &
         &   (                                                                                    &
         &     self%presetPositions                                                               &
         &    .or.                                                                                &
         &     self%presetOrbits                                                                  &
         &   )                                                                                    &
         &  .and.                                                                                 &
         &   .not.self%mergerTreeImporter_%positionsAvailable(positions=.true.,velocities=.true.) &
         & )                                                                                      &
         & call Error_Report(                                                                     &
         &                   "presetting positions requires that both position and                &
         &                    velocity datasets be present in merger tree file"   //              &
         &                   {introspection:location}                                             &
         &                  )
    ! Check that half-mass radius information is present if required.
    if     (                                                                                                                                               &
         &   self%presetScaleRadii                                                                                                                         &
         &  .and.                                                                                                                                          &
         &   .not.self%mergerTreeImporter_%scaleRadiiAvailable()                                                                                           &
         & )  call Error_Report(                                                                                                                           &
         &                      "presetting scale radii requires that at least one of readRadiusHalfMass or scaleRadius datasets be present in merger "//  &
         &                      "tree file; try setting"//char(10)//"  [presetScaleRadii]=false"                                                       //  &
         &                      {introspection:location}                                                                                                   &
         &                     )
    ! Check that angular momentum information is present if required.
    if     (                                                                                                  &
         &   self%presetAngularMomenta                                                                        &
         &  .and.                                                                                             &
         &   .not.                                                                                            &
         &        (                                                                                           &
         &          self%mergerTreeImporter_%          spinAvailable()                                        &
         &         .or.                                                                                       &
         &          self%mergerTreeImporter_%angularMomentaAvailable()                                        &
         &        )                                                                                           &
         & )                                                                                                  &
         & call Error_Report(                                                                                 &
         &                   "presetting angular momenta requires that the spin or angularMomentum datasets   &
         &                    be present in merger tree file; try setting"                                 // &
         &                   char(10)                                                                      // &
         &                   " [presetAngularMomenta]=false"                                               // &
         &                     {introspection:location}                                                       &
         &                  )
    if     (                                                                                                         &
         &   self%presetAngularMomenta3D                                                                             &
         &  .and.                                                                                                    &
         &   .not.                                                                                                   &
         &        (                                                                                                  &
         &          self%mergerTreeImporter_%          spin3DAvailable()                                             &
         &         .or.                                                                                              &
         &          self%mergerTreeImporter_%angularMomenta3DAvailable()                                             &
         &        )                                                                                                  &
         & )                                                                                                         &
         & call Error_Report(                                                                                        &
         &                   "presetting angular momentum vectors requires that the spin or angularMomentum vector   &
         &                    datasets be present in merger tree file; try setting"                               // &
         &                   char(10)                                                                             // &
         &                   " [presetAngularMomenta3D]=false"                                                    // &
         &                     {introspection:location}                                                              &
         &                  )
    ! Create an OpenMP lock that will allow threads to coordinate access to split forest data.
    !$ call OMP_Init_Lock(self%splitForestLock)
    ! Create named datasets if necessary.
    if     (                                    &
         &   size(self%presetNamedReals   ) > 0 &
         &  .or.                                &
         &   size(self%presetNamedIntegers) > 0 &
         &  ) then
       if (size(self%presetNamedReals   ) > 0) then
          allocate(self%indexNamedReals   (size(self%presetNamedReals   )))
          do i=1,size(self%presetNamedReals   )
             !![
	     <addMetaProperty component="basic" name="'preset:'//char(self%presetNamedReals   (i))" type="float"       id="self%indexNamedReals  (i)" isCreator="yes"/>
             !!]
          end do
       end if
       if (size(self%presetNamedIntegers) > 0) then
          allocate(self%indexNamedIntegers(size(self%presetNamedIntegers)))
          do i=1,size(self%presetNamedIntegers)
             !![
	     <addMetaProperty component="basic" name="'preset:'//char(self%presetNamedIntegers(i))" type="longInteger" id="self%indexNamedIntegers(i)" isCreator="yes"/>
             !!]
          end do
       end if
    end if
    return
  end function readConstructorInternal

  subroutine readDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeConstructorRead} merger tree constructor class.
    !!}
    implicit none
    type(mergerTreeConstructorRead), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%mergerTreeImporter_"            />
    <objectDestructor name="self%darkMatterHaloScale_"           />
    <objectDestructor name="self%darkMatterProfileDMO_"          />
    <objectDestructor name="self%darkMatterProfileConcentration_"/>
    <objectDestructor name="self%haloSpinDistribution_"          />
    <objectDestructor name="self%satelliteMergingTimescales_"    />
    <objectDestructor name="self%virialOrbit_"                   />
    <objectDestructor name="self%outputTimes_"                   />
    <objectDestructor name="self%darkMatterProfileScaleRadius_"  />
    <objectDestructor name="self%nodeOperator_"                  />
    <objectDestructor name="self%randomNumberGenerator_"         />
    <objectDestructor name="self%mergerTreeSeeds_"               />
    !!]
    return
  end subroutine readDestructor

  function readConstruct(self,treeNumber,finished) result(tree)
    !!{
    Construct a merger tree by reading its definition from file.
    !!}
    use    :: Array_Utilities           , only : operator(.intersection.)
    use    :: Arrays_Search             , only : searchArrayClosest
    use    :: Functions_Global          , only : State_Retrieve_                  , State_Store_
    use    :: Error                     , only : Component_List                   , Error_Report
    use    :: Galacticus_Nodes          , only : defaultDarkMatterProfileComponent, defaultPositionComponent, defaultSatelliteComponent, defaultSpinComponent, &
          &                                      mergerTree                       , treeNodeList
    use    :: Merger_Tree_Read_Importers, only : nodeData                         , nodeDataMinimal
    use    :: Merger_Tree_State_Store   , only : treeStateStoreSequence
    use    :: Merger_Tree_Walkers       , only : mergerTreeWalkerAllNodes
    use    :: Numerical_Comparison      , only : Values_Agree
    !$ use :: OMP_Lib                   , only : OMP_Unset_Lock
    use    :: Sorting                   , only : sort
    use    :: String_Handling           , only : operator(//)
    use    :: Vectors                   , only : Vector_Magnitude                 , Vector_Product
    implicit none
    type            (mergerTree               ), pointer                              :: tree
    class           (mergerTreeConstructorRead), intent(inout)                        :: self
    integer         (c_size_t                 ), intent(in   )                        :: treeNumber
    logical                                    , intent(  out)                        :: finished
    integer         (c_size_t                 )                                       :: treeNumberInternal
    integer         (kind_int8                ), allocatable, dimension(:  )          :: historyIndex
    double precision                           , allocatable, dimension(:  )          :: historyMass           , historyTime
    double precision                           , allocatable, dimension(:,:)          :: position              , velocity
    class           (nodeDataMinimal          ), allocatable, dimension(:  ), target  :: nodes
    type            (treeNode                 )                             , pointer :: nodeWork
    type            (treeNodeList             ), allocatable, dimension(:  )          :: nodeList
    logical                                    , allocatable, dimension(:  )          :: childIsSubhalo
    double precision                                        , dimension(3  )          :: relativePosition      , relativeVelocity , &
         &                                                                               orbitalAngularMomentum
    integer         (c_size_t                 ), allocatable, dimension(:  )          :: nodeSubset
    integer                                                                           :: isolatedNodeCount
    integer         (c_size_t                 )                                       :: historyCountMaximum   , iNode            , &
         &                                                                               iOutput               , treeNumberMaximum, &
         &                                                                               treeNumberOffset
    logical                                                                           :: returnSplitForest
    type            (mergerTreeWalkerAllNodes )                                       :: treeWalkerAll
    type            (varying_string           )                                       :: message

    ! Snapshot the state of the next tree to read.
    treeStateStoreSequence=treeNumber
    ! Retrieve stored internal state if possible.
    call State_Retrieve_()
    ! Recover the state of the next tree to read.
    treeNumberInternal=treeStateStoreSequence
    ! Determine if we have any split forests to return.
    returnSplitForest=allocated(self%splitForestTreeSize)
    ! Scan trees until we find one to process. This allows us to skip trees until the tree index specified by `beginAt` is found.
    do while (.true.)
       ! Find the maximum tree number in the current file.
       treeNumberMaximum=int(self%mergerTreeImporter_%treeCount(),kind=c_size_t)
       ! Check if we need to move to a new file.
       do while (treeNumber-self%treeNumberOffset > treeNumberMaximum .and. self%fileCurrent < size(self%fileNames))
          self%fileCurrent     =self%fileCurrent     +1
          self%treeNumberOffset=self%treeNumberOffset+treeNumberMaximum
          call self%mergerTreeImporter_%open(self%fileNames(self%fileCurrent))
          treeNumberMaximum=int(self%mergerTreeImporter_%treeCount(),kind=c_size_t)
       end do
       treeNumberOffset=treeNumber-self%treeNumberOffset
       ! If all trees are used up, we're done.
       if (treeNumberOffset > treeNumberMaximum) then
          nullify(tree)
          finished=.true.
          exit
       end if
       ! Test if this is the first forest that we are to process, or if we have already found that forest.
       if (.not.self%foundBeginAt .and. self%mergerTreeImporter_%treeIndex (int(treeNumberOffset)) /= self%beginAt) then
          ! Decrement our internal offset, and try again - this will move us to trying the next forest in the file.
          self%treeNumberOffset=self%treeNumberOffset-1_c_size_t
          cycle
       else
          ! The first forest to process is found - record this.
          self%foundBeginAt=.true.
       end if
       ! Set tree properties.
       allocate(tree)
       ! treeIndex
       tree%index            =self%mergerTreeImporter_%treeIndex (int(treeNumberOffset))
       ! volumeWeight
       self%treeWeightCurrent=self%mergerTreeImporter_%treeWeight(int(treeNumberOffset))
       tree%volumeWeight     =self%treeWeightCurrent
       ! Initialize no events.
       tree%event             => null()
       tree%initializedUntil  =  0.0d0
       tree%isTreeInitialized =  .false.
       call tree%properties%initialize()
       ! Restart the random number sequence for this tree. We use the tree index modulo the largest number representable by
       ! the integer type.
       allocate(tree%randomNumberGenerator_,mold=self%randomNumberGenerator_)
       !$omp critical(mergerTreeConstructReadDeepCopyReset)
       !![
       <deepCopyReset variables="self%randomNumberGenerator_"/>
       <deepCopy source="self%randomNumberGenerator_" destination="tree%randomNumberGenerator_"/>
       <deepCopyFinalize variables="tree%randomNumberGenerator_"/>
       !!]
       !$omp end critical(mergerTreeConstructReadDeepCopyReset)
       call self                 %randomSequenceNonDeterministicWarn(tree)
       call self%mergerTreeSeeds_%set                               (tree)
       ! Store internal state.
       message='Storing state for tree #'
       message=message//treeStateStoreSequence
       call State_Store_(message)
       ! Check if the size of this forest exceeds the maximum allowed.
       if     (                                                                                    &
            &   .not.returnSplitForest                                                             &
            &  .and.                                                                               &
            &   self%mergerTreeImporter_%nodeCount(int(treeNumberOffset)) > self%forestSizeMaximum &
            &  .and.                                                                               &
            &   0                                                         < self%forestSizeMaximum &
            & ) then
#ifdef USEMPI
          call Error_Report('split forest processing is not supported under MPI'                        //{introspection:location})
#else
          call Error_Report('split forest processing is currently broken until MPI support is completed'//{introspection:location})
#endif
          ! Check if the importer supports reading subsets of halos from a forest.
          if (.not.self%mergerTreeImporter_%canReadSubsets()) call Error_Report('forest exceeds maximum allowed size but importer cannot read subsets of halos'//{introspection:location})
          ! Import nodes, and keep only the minimally required data to map the tree structure.
          call self%mergerTreeImporter_%import(int(treeNumberOffset),nodes,structureOnly=.true.)
          ! Find initial root node affinities of all nodes.
          call self%rootNodeAffinitiesInitial(nodes)
          deallocate(nodes)
          returnSplitForest      =.true.
          self%splitForestNextTree    =0
          self%splitForestActiveForest=int(treeNumberOffset)
       end if
       ! Determine subset of nodes to read.
       if (returnSplitForest) then
          ! Move to the next tree.
          self%splitForestNextTree=self%splitForestNextTree+1
          allocate(nodeSubset(self%splitForestTreeSize(self%splitForestNextTree)))
          nodeSubset=self%splitForestMapIndex(self%splitForestTreeStart(self%splitForestNextTree):self%splitForestTreeStart(self%splitForestNextTree)+self%splitForestTreeSize(self%splitForestNextTree)-1)
          call sort(nodeSubset)
       else
          allocate(nodeSubset(1))
          nodeSubset=[-1_c_size_t]
       end if
       ! Read data from the file.
       !![
       <conditionalCall>
       <call>
         call self%mergerTreeImporter_%import(                                                                                                                      &amp;
            &amp;                             int(treeNumberOffset)                                                                                               , &amp;
            &amp;                             nodes                                                                                                               , &amp;
            &amp;                             requireScaleRadii         = self%presetScaleRadii                                                                   , &amp;
            &amp;                             requireAngularMomenta     =(self%presetAngularMomenta     .and.self%mergerTreeImporter_%angularMomentaAvailable  ()), &amp;
            &amp;                             requireAngularMomenta3D   =                                                                                           &amp;
            &amp;                                                        (                                                                                          &amp;
            &amp;                                                           self%presetAngularMomenta3D                                                             &amp;
            &amp;                                                          .or.                                                                                     &amp;
            &amp;                                                           (                                                                                       &amp;
            &amp;                                                             self%presetAngularMomenta                                                             &amp;
            &amp;                                                            .and.                                                                                  &amp;
            &amp;                                                             self%subhaloAngularMomentaMethod == subhaloAngularMomentaMethodSummation              &amp;
            &amp;                                                           )                                                                                       &amp;
            &amp;                                                         )                                                                                         &amp;
            &amp;                                                        .and.                                                                                      &amp;
            &amp;                                                         self%mergerTreeImporter_%angularMomenta3DAvailable()                                    , &amp;
            &amp;                             requireSpin               =(self%presetAngularMomenta     .and.self%mergerTreeImporter_%spinAvailable            ()), &amp;
            &amp;                             requireSpin3D             =(self%presetAngularMomenta3D   .and.self%mergerTreeImporter_%spin3DAvailable          ()), &amp;
            &amp;                             requirePositions          =(self%presetPositions          .or. self%presetOrbits                                   ), &amp;
            &amp;                             nodeSubset                =nodeSubset                                                                                 &amp;
            &amp;                             {conditions}                                                                                                          &amp;
            &amp;                            )
        </call>
        <argument name="requireNamedReals"    value="self%presetNamedReals"    condition="size(self%presetNamedReals   ) > 0"/>
        <argument name="requireNamedIntegers" value="self%presetNamedIntegers" condition="size(self%presetNamedIntegers) > 0"/>
       </conditionalCall>
       !!]
       deallocate(nodeSubset)
       select type (nodes)
       class is (nodeData)
          ! Snap node times to output times if a tolerance has been specified.
          if (self%outputTimeSnapTolerance > 0.0d0) then
             ! Loop over all nodes.
             do iNode=1,size(nodes)
                ! Find closest output time to the node time.
                iOutput=searchArrayClosest(self%outputTimes,nodes(iNode)%nodeTime)
                ! Test if this time is sufficiently close that we should snap the node time to it.
                if (Values_Agree(nodes(iNode)%nodeTime,self%outputTimes(iOutput),relTol=self%outputTimeSnapTolerance)) &
                     & nodes(iNode)%nodeTime=self%outputTimes(iOutput)
             end do
          end if
          ! Sort node indices.
          call self%createNodeIndices      (nodes)
          ! Identify subhalos.
          nodes%isSubhalo=nodes%nodeIndex /= nodes%hostIndex
          ! Build pointers to descendant nodes.
          call self%buildDescendantPointers(nodes)
          ! Find cases where something that was a subhalo stops being a subhalo and prevent them if necessary.
          call self%enforceSubhaloStatus   (nodes)
          ! If necessary, add masses and angular momenta of subhalos to host halos.
          if (.not.self%mergerTreeImporter_%angularMomentaIncludeSubhalos().and.self%subhaloAngularMomentaMethod == subhaloAngularMomentaMethodScale) then
             ! This method requires angular momenta to be available.
             if (.not.self%mergerTreeImporter_%angularMomentaAvailable()) call Error_Report('scaling parent angular momentum for subhalo masses requires angular momenta availability'//{introspection:location})
             do iNode=1,size(nodes)
                if (nodes(iNode)%host%nodeIndex == nodes(iNode)%nodeIndex) then
                   if (self%presetAngularMomenta  )            &
                        & nodes      (iNode)%angularMomentum   &
                        &      =nodes(iNode)%angularMomentum   &
                        &      /nodes(iNode)%nodeMass
                   if (self%presetAngularMomenta3D)            &
                        & nodes      (iNode)%angularMomentum3D &
                        &      =nodes(iNode)%angularMomentum3D &
                        &      /nodes(iNode)%nodeMass
                end if
             end do
          end if
          if (.not.self%mergerTreeImporter_%massesIncludeSubhalos()) then
             do iNode=1,size(nodes)
                if (nodes(iNode)%host%nodeIndex /= nodes(iNode)%nodeIndex) nodes  (iNode)%host%nodeMass &
                     &                                                      =nodes(iNode)%host%nodeMass &
                     &                                                      +nodes(iNode)%nodeMass
             end do
          end if
          if (.not.self%mergerTreeImporter_%angularMomentaIncludeSubhalos().and.self%subhaloAngularMomentaMethod == subhaloAngularMomentaMethodScale) then
             do iNode=1,size(nodes)
                if (nodes(iNode)%host%nodeIndex == nodes(iNode)%nodeIndex) then
                   if (self%presetAngularMomenta  )            &
                        & nodes      (iNode)%angularMomentum   &
                        &      =nodes(iNode)%angularMomentum   &
                        &      *nodes(iNode)%nodeMass
                   if (self%presetAngularMomenta3D)            &
                        & nodes      (iNode)%angularMomentum3D &
                        &      =nodes(iNode)%angularMomentum3D &
                        &      *nodes(iNode)%nodeMass

                end if
             end do
          end if
          if     (                                                                          &
               &  (                                                                         &
               &    self%presetAngularMomenta3D                                             &
               &   .or.                                                                     &
               &    self%presetAngularMomenta                                               &
               &  )                                                                         &
               &  .and.                                                                     &
               &   .not.self%mergerTreeImporter_%angularMomentaIncludeSubhalos()            &
               &  .and.                                                                     &
               &   self%subhaloAngularMomentaMethod == subhaloAngularMomentaMethodSummation &
               & ) then
             ! This method requires 3D angular momenta to be available.
             if (.not.self%mergerTreeImporter_%angularMomenta3DAvailable())                                                                  &
                  & call Error_Report(                                                                                                       &
                  &                   'adding subhalo angular momenta to parent angular momentum requires 3D angular momenta availability'// &
                  &                   {introspection:location}                                                                               &
                  &                  )
             do iNode=1,size(nodes)
                if (nodes(iNode)%host%nodeIndex /= nodes(iNode)%nodeIndex) then
                   ! Find relative position and velocity.
                   relativePosition=nodes(iNode)%position-nodes(iNode)%host%position
                   relativeVelocity=nodes(iNode)%velocity-nodes(iNode)%host%velocity
                   ! Update position/velocity for periodicity and Hubble flow.
                   call self%phaseSpacePositionRealize(nodes(iNode)%nodeTime,relativePosition,relativeVelocity)
                   ! Compute orbital angular momentum of subhalo.
                   orbitalAngularMomentum=+nodes(iNode)%nodeMass                             &
                        &                 *Vector_Product(relativePosition,relativeVelocity)
                   ! Sum orbital and internal angular momenta.
                   nodes        (iNode)%host%angularMomentum3D &
                        & =nodes(iNode)%host%angularMomentum3D &
                        & +nodes(iNode)     %angularMomentum3D &
                        & +orbitalAngularMomentum
                end if
             end do
             ! Update scalar angular momenta.
             do iNode=1,size(nodes)
                if (nodes(iNode)%host%nodeIndex == nodes(iNode)%nodeIndex) &
                     & nodes(iNode)%angularMomentum=Vector_Magnitude(nodes(iNode)%angularMomentum3D)
             end do
          end if
          ! Associate parent pointers with the descendant host.
          call readBuildParentPointers         (     nodes                                          )
          ! Create an array of standard nodes.
          call self%createNodeArray            (tree,nodes,nodeList,isolatedNodeCount,childIsSubhalo)
          ! Assign parent pointers and properties.
          call self%buildIsolatedParentPointers(tree,nodes,nodeList                                 )
          ! Now build child and sibling links.
          call readBuildChildAndSiblingLinks   (     nodes,nodeList                  ,childIsSubhalo)
          ! (Re)assign host tree pointers.
          call readAssignHostTreePointers      (tree                                                )
          ! Assign split forest events.
          call self%assignSplitForestEvents    (     nodes,nodeList                                 )
          ! Check that all required properties exist.
          if (self%presetPositions.or.self%presetOrbits) then
             ! Position and velocity methods are required.
             if     (                                                                                                                                                &
                  &  .not.(                                                                                                                                          &
                  &         defaultPositionComponent%positionIsSettable()                                                                                            &
                  &        .and.                                                                                                                                     &
                  &         defaultPositionComponent%velocityIsSettable()                                                                                            &
                  &       )                                                                                                                                          &
                  & )                                                                                                                                                &
                  & call Error_Report                                                                                                                                &
                  &      (                                                                                                                                           &
                  &       'presetting positions or orbits requires a component that supports position and velocity setting (e.g. set [componentPosition]=preset);'// &
                  &       Component_List(                                                                                                                            &
                  &                      'position'                                                                                                                , &
                  &                       defaultPositionComponent        %             positionAttributeMatch(requireSettable=.true.)                               &
                  &                      .intersection.                                                                                                              &
                  &                       defaultPositionComponent        %             velocityAttributeMatch(requireSettable=.true.)                               &
                  &                     )                                                                                                                         // &
                  &       char(10)                                                                                                                                // &
                  &       'alternatively setting [presetPositions]=false and [presetOrbits]=false will remove the need to store positions and velocities'         // &
                  &       {introspection:location}                                                                                                                   &
                  & )
          end if
          if (self%presetMergerTimes     ) then
             ! Time of merging property is required.
             if (.not.defaultSatelliteComponent%timeOfMergingIsSettable                ())                                                                           &
                  & call Error_Report                                                                                                                                &
                  &      (                                                                                                                                           &
                  &       'presetting merging times requires a component that supports setting of merging times.'                                                 // &
                  &       Component_List(                                                                                                                            &
                  &                      'satellite'                                                                                                               , &
                  &                       defaultSatelliteComponent       %        timeOfMergingAttributeMatch(requireSettable=.true.)                               &
                  &                      )                                                                                                                        // &
                  &       {introspection:location}                                                                                                                   &
                  &      )
          end if
          if (self%presetScaleRadii      ) then
             ! Scale radius property is required.
             if (.not.defaultDarkMatterProfileComponent%scaleIsSettable                ())                                                                           &
                  & call Error_Report                                                                                                                                &
                  &      (                                                                                                                                           &
                  &       'presetting scale radii requires a component that supports setting of scale radii.'                                                     // &
                  &       Component_List(                                                                                                                            &
                  &                      'darkMatterProfile'                                                                                                       , &
                  &                      defaultDarkMatterProfileComponent%                scaleAttributeMatch(requireSettable=.true.)                               &
                  &                     )                                                                                                                         // &
                  &       {introspection:location}                                                                                                                   &
                  &      )
          end if
          if (self%presetAngularMomenta  ) then
             ! Angular momentum property is required.
             if (.not.defaultSpinComponent             %angularMomentumIsSettable      ())                                                                           &
                  & call Error_Report                                                                                                                                &
                  &      (                                                                                                                                           &
                  &       'presetting angular momenta requires a component that supports setting of angular momenta.'                                             // &
                  &       Component_List(                                                                                                                            &
                  &                      'spin'                                                                                                                    , &
                  &                      defaultSpinComponent             %      angularMomentumAttributeMatch(requireSettable=.true.)                               &
                  &                     )                                                                                                                         // &
                  &       {introspection:location}                                                                                                                   &
                  &      )
          end if
          if (self%presetAngularMomenta3D) then
             ! Spin property is required.
             if (.not.defaultSpinComponent             %angularMomentumVectorIsSettable())                                                                           &
                  & call Error_Report                                                                                                                                &
                  &      (                                                                                                                                           &
                  &       'presetting angular momentum vectors requires a component that supports setting of angular momentum vectors.'                           // &
                  &       Component_List(                                                                                                                            &
                  &                      'spinVector'                                                                                                              , &
                  &                      defaultSpinComponent             %angularMomentumVectorAttributeMatch(requireSettable=.true.)                               &
                  &                     )                                                                                                                         // &
                  &       {introspection:location}                                                                                                                   &
                  &      )
          end if
          if (self%presetOrbits     ) then
             ! Orbit property is required.
             if (.not.defaultSatelliteComponent        %virialOrbitIsSettable())                                                                                     &
                  & call Error_Report                                                                                                                                &
                  &      (                                                                                                                                           &
                  &       'presetting orbits requires a component that supports setting of orbits (e.g. [componentSatellite]=preset);'                            // &
                  &       Component_List(                                                                                                                            &
                  &                      'satellite'                                                                                                               , &
                  &                      defaultSatelliteComponent        %          virialOrbitAttributeMatch(requireSettable=.true.)                               &
                  &                     )                                                                                                                         // &
                  &       char(10)                                                                                                                                // &
                  &       'Alternatively, set [presetOrbits]=false to prevent attempts to set orbits)'                                                            // &
                  &       {introspection:location}                                                                                                                   &
                  &      )
          end if
          ! Apply any tree initialization operators.
          treeWalkerAll=mergerTreeWalkerAllNodes(tree,spanForest=.true.)
          do while (treeWalkerAll%next(nodeWork))
             call self%nodeOperator_%nodeTreeInitialize(nodeWork)
          end do
          ! Assign named properties.
          if     (                                    &
               &   size(self%presetNamedReals   ) > 0 &
               &  .or.                                &
               &   size(self%presetNamedIntegers) > 0 &
               & )                                    &
               &                                   call self%assignNamedProperties   (nodes,nodeList                    )
          ! Assign scale radii.
          if     ( self%presetScaleRadii         ) call self%assignScaleRadii        (nodes,nodeList                    )
          ! Assign spin parameters.
          if     (                             &
               &   self%presetAngularMomenta   &
               &  .or.                         &
               &   self%presetAngularMomenta3D &
               & )                                 call self%assignAngularMomenta    (nodes,nodeList                    )
          ! Assign isolated node indices to subhalos.
          call readAssignIsolatedNodeIndices                                         (nodes                             )
          ! Ensure that isolated nodes with progenitors that descend into subhalos have valid primary progenitors.
          call readValidateIsolatedHalos                                             (nodes                             )
          ! Scan subhalos to determine when and how they merge.
          call self%scanForMergers                                                   (nodes,nodeList,historyCountMaximum)
          ! If a split forest was used, but all trees from it have now been processed, remove the split forest data as we no
          ! longer need it at this point.
          if (returnSplitForest) then
             if (self%splitForestNextTree == size(self%splitForestTreeStart)) then
                deallocate(self%splitForestTreeSize )
                deallocate(self%splitForestTreeStart)
                deallocate(self%splitForestPushTo   )
                deallocate(self%splitForestPullFrom )
                deallocate(self%splitForestPushType )
                deallocate(self%splitForestMapIndex )
             end if
          end if
          ! Release the lock on split forest data if necessary as we're finished using it. This allows other threads to begin
          ! using the split forest data.
          !$ if (self%forestSizeMaximum > 0) call OMP_Unset_Lock(self%splitForestLock)
          ! Search for any nodes which were flagged as merging with another node and assign appropriate pointers.
          call self%assignMergers           (nodes,nodeList)
          ! Find cases where something that was a subhalo stops being a subhalo and add events to handle.
          call self%scanForSubhaloPromotions(nodes,nodeList)
          ! Search for subhalos which move between branches/trees.
          call self%scanForBranchJumps      (nodes,nodeList)
          ! Allocate arrays for history building.
          if (allocated(position)) deallocate(position)
          if (allocated(velocity)) deallocate(velocity)
          allocate(historyTime(int(historyCountMaximum)))
          if (self%presetSubhaloIndices                             ) then
             allocate(historyIndex(int(historyCountMaximum)))
          else
             allocate(historyIndex(0 ))
          end if
          if (self%presetSubhaloMasses                              ) then
             allocate(historyMass (int(historyCountMaximum)))
          else
             allocate(historyMass (0 ))
          end if
          if (self%presetPositions    .or.self%presetOrbits) then
             allocate(position    (3,int(historyCountMaximum)))
             allocate(velocity    (3,int(historyCountMaximum)))
          else
             allocate(position    (0,                      0 ))
             allocate(velocity    (0,                      0 ))
          end if
          ! Build subhalo mass histories if required.
          call self%buildSubhaloMassHistories(nodes,nodeList,historyCountMaximum,historyTime,historyIndex,historyMass,position,velocity)
          ! Assign new uniqueIDs to any cloned nodes inserted into the trees.
          call readAssignUniqueIDsToClones(nodeList)
          ! Deallocate history building arrays.
          if (allocated(historyTime)) deallocate(historyTime )
          if (allocated(historyMass)) deallocate(historyIndex)
          if (allocated(historyMass)) deallocate(historyMass )
          if (allocated(position   )) deallocate(position    )
          if (allocated(velocity   )) deallocate(velocity    )
          ! Deallocate the temporary arrays.
          deallocate(nodeList)
          ! Destroy sorted node indices.
          call self%destroyNodeIndices()
       class default
          call Error_Report('nodes arrays is of wrong class'//{introspection:location})
       end select
       ! Deallocate nodes.
       deallocate(nodes)
       ! Indicate that we are not finished.
       exit
    end do
    finished=.not.associated(tree)
    return
  end function readConstruct

  subroutine readCreateNodeIndices(self,nodes)
    !!{
    Create a sorted list of node indices with an index into the original array.
    !!}
    use :: Error                     , only : Warn
    use :: Merger_Tree_Read_Importers, only : nodeDataMinimal
    use :: Sorting                   , only : sortIndex
    use :: String_Handling           , only : operator(//)
    implicit none
    class  (mergerTreeConstructorRead), intent(inout)                :: self
    class  (nodeDataMinimal          ), dimension(:) , intent(inout) :: nodes
    integer                                                          :: iNode
    type   (varying_string           )                               :: message

    ! Build a sorted list of node indices with an index into the original arrays.
    allocate(self%nodeLocations          (size(nodes)))
    allocate(self%nodeIndicesSorted      (size(nodes)))
    allocate(self%descendantLocations    (size(nodes)))
    allocate(self%descendantIndicesSorted(size(nodes)))
    self%nodeLocations      =sortIndex(nodes%nodeIndex      )
    self%descendantLocations=sortIndex(nodes%descendantIndex)
    forall (iNode=1:size(nodes))
       self%nodeIndicesSorted      (iNode)=nodes(self%nodeLocations      (iNode))%nodeIndex
       self%descendantIndicesSorted(iNode)=nodes(self%descendantLocations(iNode))%descendantIndex
    end forall
    do iNode=2,size(nodes)
       if (self%nodeIndicesSorted(iNode) == self%nodeIndicesSorted(iNode-1)) then
          message="WARNING: duplicate node index found in merger tree - this is not allowed ["
          message=message//self%nodeIndicesSorted(iNode)//']'
          call Warn(message)
       end if
    end do
    return
  end subroutine readCreateNodeIndices

  function readNodeLocation(self,nodeIndex)
    !!{
    Return the location in the original array of the given {\normalfont \ttfamily nodeIndex}.
    !!}
    use :: Arrays_Search, only : searchArray
    implicit none
    integer(c_size_t                 )                :: readNodeLocation
    class  (mergerTreeConstructorRead), intent(inout) :: self
    integer(kind_int8                ), intent(in   ) :: nodeIndex
    integer(c_size_t                 )                :: iNode

    iNode=searchArray(self%nodeIndicesSorted,nodeIndex)
    if (iNode < 1 .or. iNode > size(self%nodeIndicesSorted)) then
       readNodeLocation=1
    else
       readNodeLocation=self%nodeLocations(iNode)
    end if
    return
  end function readNodeLocation

  function readDescendantNodeSortIndex(self,descendantIndex)
    !!{
    Return the sort index of the given {\normalfont \ttfamily descendantIndex}.
    !!}
    use :: Arrays_Search, only : searchArray
    implicit none
    integer(c_size_t                 )                :: readDescendantNodeSortIndex
    class  (mergerTreeConstructorRead), intent(inout) :: self
    integer(kind_int8                ), intent(in   ) :: descendantIndex

    readDescendantNodeSortIndex=searchArray(self%descendantIndicesSorted,descendantIndex)
    return
  end function readDescendantNodeSortIndex

  subroutine readDestroyNodeIndices(self)
    !!{
    Destroy the sorted list of node indices.
    !!}
    implicit none
    class(mergerTreeConstructorRead), intent(inout) :: self

    if (allocated(self%nodeLocations          )) deallocate(self%nodeLocations          )
    if (allocated(self%nodeIndicesSorted      )) deallocate(self%nodeIndicesSorted      )
    if (allocated(self%descendantLocations    )) deallocate(self%descendantLocations    )
    if (allocated(self%descendantIndicesSorted)) deallocate(self%descendantIndicesSorted)
    return
  end subroutine readDestroyNodeIndices

  subroutine readBuildDescendantPointers(self,nodes)
    !!{
    Builds pointers from each node to its descendant node.
    !!}
    use :: Display                   , only : displayMessage, verbosityLevelInfo
    use :: Error                     , only : Error_Report
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: String_Handling           , only : operator(//)
    implicit none
    class  (mergerTreeConstructorRead), intent(inout)                        :: self
    class  (nodeData                 ), dimension(:) , intent(inout), target :: nodes
    integer(c_size_t                 )                                       :: iNode  , nodeLocation
    type   (varying_string           )                                       :: message

    do iNode=1,size(nodes)
       ! Does this node have a descendant? And is it staying in this tree?
       if (nodes(iNode)%descendantIndex >= 0.and..not.self%isOnPushList(nodes(iNode))) then
          nodeLocation=self%nodeLocation(nodes(iNode)%descendantIndex)
          if (nodes(nodeLocation)%nodeIndex /= nodes(iNode)%descendantIndex) then
             message='failed to find descendant node: '
             message=message//nodes(iNode)%descendantIndex//' of '//nodes(iNode)%nodeIndex
             call Error_Report(message//{introspection:location})
          end if
          nodes(iNode)%descendant => nodes(nodeLocation)
       else
          nodes(iNode)%descendant => null()
       end if
       if (nodes(iNode)%hostIndex       >= 0) then
          nodeLocation=self%nodeLocation(nodes(iNode)%hostIndex      )
          if (nodes(nodeLocation)%nodeIndex /= nodes(iNode)%hostIndex       ) then
             message='failed to find host node: '
             message=message//nodes(iNode)%hostIndex//' of '//nodes(iNode)%nodeIndex
             if (self%missingHostsAreFatal) then
                call Error_Report(message//{introspection:location})
             else
                message=message//" - resetting this node to be an isolated node"
                call displayMessage(message,verbosity=verbosityLevelInfo)
                nodes(iNode)%hostIndex =  nodes(iNode)%nodeIndex
                nodes(iNode)%host      => nodes(iNode)
                nodes(iNode)%isSubhalo =  .false.
             end if
          else
             nodes(iNode)%host    => nodes(nodeLocation)
          end if
       else
          call Error_Report('negative values are not allowed for hostIndex - if node is self-hosting [i.e. not a subhalo] set hostIndex=nodeIndex'//{introspection:location})
       end if
    end do
    return
  end subroutine readBuildDescendantPointers

  subroutine readEnforceSubhaloStatus(self,nodes)
    !!{
    Ensure that any node which was once a subhalo remains a subhalo.
    !!}
    use :: Error                     , only : Error_Report
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: String_Handling           , only : operator(//)
    implicit none
    class  (mergerTreeConstructorRead), intent(inout)                        :: self
    class  (nodeData                 ), dimension(:) , intent(inout), target :: nodes
    class  (nodeData                 ), pointer                              :: descendantNode, progenitorNode
    integer(c_size_t                 )                                       :: iNode
    logical                                                                  :: failed        , isolatedProgenitorExists
    type   (varying_string           )                                       :: message
    type   (progenitorIterator       )                                       :: progenitors

    ! Return immediately if subhalo promotions are allowed.
    if (self%allowSubhaloPromotions) return
    ! Subhalo promotions are not allowed, so enforce subhalo status.
    do iNode=1,size(nodes)
       if (nodes(iNode)%isSubhalo) then
          descendantNode => nodes(iNode)%descendant
          do while (associated(descendantNode))
             ! Is this node isolated?
             if (.not.descendantNode%isSubhalo) then
                ! Check if there is any isolated node which descends into this node.
                isolatedProgenitorExists=.false.
                call progenitors%descendantSet(self,descendantNode,nodes)
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
             descendantNode => descendantNode%descendant
          end do
       end if
    end do
    ! Check that subhalo enforcement was successful.
    failed=.false.
    do iNode=1,size(nodes)
       ! Find nodes which have no isolated node descending into them.
       isolatedProgenitorExists=.false.
       call progenitors%descendantSet(self,nodes(iNode),nodes)
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
       call Error_Report(message//{introspection:location})
    end if
    return
  end subroutine readEnforceSubhaloStatus

  subroutine readScanForSubhaloPromotions(self,nodes,nodeList)
    !!{
    Scan for cases where a subhalo stops being a subhalo and so must be promoted.
    !!}
    use :: Galacticus_Nodes            , only : nodeEvent                  , nodeEventSubhaloPromotion, treeNode, treeNodeList, &
         &                                      nodeComponentSatellite
    use :: Merger_Tree_Read_Importers  , only : nodeData
    use :: Node_Subhalo_Promotions     , only : nodeSubhaloPromotionPerform
    use :: Satellite_Merging_Timescales, only : satelliteMergeTimeInfinite
    implicit none
    class           (mergerTreeConstructorRead), intent(inout)                              :: self
    class           (nodeData                 ), target       , dimension(:), intent(inout) :: nodes
    type            (treeNodeList             )               , dimension(:), intent(inout) :: nodeList
    class           (nodeData                 ), pointer                                    :: descendantNode          , progenitorNode
    class           (nodeEvent                ), pointer                                    :: newEvent                , pairEvent
    type            (treeNode                 ), pointer                                    :: promotionNode           , node             , &
         &                                                                                     nodeNew
    class           (nodeComponentBasic       ), pointer                                    :: basic
    class           (nodeComponentSatellite   ), pointer                                    :: satellite
    integer         (c_size_t                 )                                             :: iNode
    integer                                                                                 :: i
    logical                                                                                 :: isolatedProgenitorExists, nodeIsMostMassive, &
         &                                                                                     progenitorIsIsolated
    type            (progenitorIterator       )                                             :: progenitors
    double precision                                                                        :: timeSubhaloPromotion

    ! Return immediately if subhalo promotion is not allowed.
    if (.not.self%allowSubhaloPromotions) return
    ! Find subhalos to be promoted.
    do iNode=1,size(nodes)
       if (nodes(iNode)%isSubhalo.and.associated(nodes(iNode)%descendant)) then
          descendantNode => nodes(iNode)%descendant
          ! Is this node isolated?
          if (.not.descendantNode%isSubhalo) then
             ! Check if there is any isolated node which descends into this node, and also whether this is the most massive
             ! subhalo which descends into the descendant.
             isolatedProgenitorExists=.false.
             nodeIsMostMassive       =.true.
             call progenitors%descendantSet(self,descendantNode,nodes)
             do while (progenitors%next(nodes))
                progenitorNode       => progenitors   %current  (nodes)
                progenitorIsIsolated =  progenitorNode%nodeIndex        == progenitorNode%hostIndex
                if (progenitorIsIsolated) isolatedProgenitorExists=.true.
                if (self%alwaysPromoteMostMassive) then
                   ! Find the most massive progenitor. We treat isolated and non-isolated nodes differently, such that if there
                   ! are two nodes of precisely equal mass (which can happen due to the discrete nature of masses in N-body
                   ! simulations for example), one isolated and one not isolated, we treat the isolated one as the primary
                   ! progenitor.
                   if     (                                                                                    &
                        &   (.not.progenitorIsIsolated .and. progenitorNode%nodeMass >  nodes(iNode)%nodeMass) &
                        &  .or.                                                                                &
                        &   (     progenitorIsIsolated .and. progenitorNode%nodeMass >= nodes(iNode)%nodeMass) &
                        & ) nodeIsMostMassive=.false.
                else
                   if     (                                                                                    &
                        &    .not.progenitorIsIsolated .and. progenitorNode%nodeMass >  nodes(iNode)%nodeMass  &
                        & ) nodeIsMostMassive=.false.
                end if
             end do
             ! Determine if we must promote the subhalo.
             if     (                                                                                  &
                  &   .not.progenitors%exist                   () .or.  .not.isolatedProgenitorExists  & ! Only promote a subhalo if no isolated halo exists.
                  &  .or.                                                                              &
                  &   (    self       %alwaysPromoteMostMassive   .and.      nodeIsMostMassive       ) & ! Always promote the most massive halo - isolated or not.
                  & ) then
                if (nodeIsMostMassive) then
                   ! Node is isolated, has no isolated node that descends into it, and our subhalo is the most massive subhalo
                   ! which descends into it, or we are always promoting the most massive progenitor, even if it is a subhalo and
                   ! there are isolated progenitors. Therefore, our subhalo must be promoted to become an isolated halo again.
                   !
                   ! If other isolated progenitors exist, we need to insert a placeholder node as the primary progenitor.
                   if (isolatedProgenitorExists) then
                      allocate(nodeNew)
                      call nodeList(descendantNode%isolatedNodeIndex)%node%copyNodeTo(nodeNew)
                      if (nodeNew%satelliteCount() > 0) then
                         ! Remove any satellite component from the copied node - each branch should have only a single satellite.
                         do i=nodeNew%satelliteCount(),1,-1
                            call nodeNew%satelliteRemove(i)
                         end do
                      end if
                      nodeNew                                         %sibling        => nodeList(descendantNode%isolatedNodeIndex)%node%firstChild
                      nodeNew                                         %parent         => nodeList(descendantNode%isolatedNodeIndex)%node
                      nodeNew                                         %firstSatellite => null()
                      nodeNew                                         %firstChild     => null()
                      nodeNew                                         %mergeTarget    => null()
                      nodeNew                                         %firstMergee    => null()
                      nodeNew                                         %siblingMergee  => null()
                      nodeNew                                         %event          => null()
                      nodeList(descendantNode%isolatedNodeIndex)%node%firstChild      => nodeNew
                      basic                                                           => nodeNew                                        %basic     ()
                      call basic%timeSet(basic%time()*(1.0d0-fractionOffsetTimeClones))
                      promotionNode        =>                                                  nodeNew
                      timeSubhaloPromotion =  basic                                           %    time()
                   else
                      promotionNode        => nodeList      (descendantNode%isolatedNodeIndex)%node
                      timeSubhaloPromotion =  descendantNode                                  %nodeTime
                   end if
                   node                    => nodeList      (nodes(iNode)  %isolatedNodeIndex)%node
                   ! If the node being promoted has a merging time set we unset it now. This was a subhalo which was flagged for
                   ! merging, but we are now promoting it as the primary progenitor of the current node, so it can no longer
                   ! merge.
                   satellite => node%satellite()
                   if (satellite%timeOfMerging() < satelliteMergeTimeInfinite) then
                      call satellite%timeOfMergingSet(satelliteMergeTimeInfinite)
                      if (associated(node%mergeTarget)) then
                         call node%removeFromMergee()
                         nullify(node%mergeTarget)
                      end if
                   end if
                   ! Create the event.
                   allocate(nodeEventSubhaloPromotion ::  newEvent)
                   allocate(nodeEventSubhaloPromotion :: pairEvent)
                   call          node%attachEvent( newEvent)
                   call promotionNode%attachEvent(pairEvent)
                   newEvent %time =  timeSubhaloPromotion
                   newEvent %node => promotionNode
                   newEvent %task => nodeSubhaloPromotionPerform
                   pairEvent%time =  timeSubhaloPromotion
                   pairEvent%node => node
                   pairEvent%task => null()
                   pairEvent%ID   =  newEvent%ID
                else if (self%allowBranchJumps) then
                   ! Node is isolated, has no isolated node that descends into it, and our subhalo is not the most massive subhalo
                   ! which descends into it. Therefore, our subhalo must branch jump if this is allowed.
                   call readCreateBranchJumpEvent(                                                 &
                        &                         nodeList(nodes(iNode)  %isolatedNodeIndex)%node, &
                        &                         nodeList(descendantNode%isolatedNodeIndex)%node, &
                        &                         descendantNode%nodeTime                          &
                        &                        )
                end if
             end if
          end if
       end if
    end do
    return
  end subroutine readScanForSubhaloPromotions

  subroutine readBuildParentPointers(nodes)
    !!{
    Build pointers to node parents.
    !!}
    use :: Error                     , only : Error_Report
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: String_Handling           , only : operator(//)
    implicit none
    class  (nodeData      ), dimension(:), intent(inout), target :: nodes
    class  (nodeData      ), pointer                             :: parentNode
    integer(c_size_t      )                                      :: iNode
    type   (varying_string)                                      :: message

    do iNode=1,size(nodes)
       if (associated(nodes(iNode)%descendant)) then
          if (nodes(iNode)%nodeIndex == nodes(iNode)%host%nodeIndex) then
             ! Find an isolated parent node, by repeatedly jumping from host to host.
             parentNode => nodes(iNode)%descendant%host
             do while (parentNode%isSubhalo)
                if (associated(parentNode,parentNode%host)) then
                   message='node ['
                   message=message//parentNode%nodeIndex//'] flagged as subhalo is self-hosting - exiting to avoid infinite loop'
                   call Error_Report(message//{introspection:location})
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
             call Error_Report(message//{introspection:location})
          end if
       end if
    end do
    return
  end subroutine readBuildParentPointers

  subroutine readCreateNodeArray(self,tree,nodes,nodeList,isolatedNodeCount,childIsSubhalo)
    !!{
    Create an array of standard nodes and associated structures.
    !!}
    use :: Galacticus_Nodes          , only : mergerTree         , treeNode     , treeNodeList
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class  (mergerTreeConstructorRead)                           , intent(inout) :: self
    type   (mergerTree               )                           , intent(inout) :: tree
    class  (nodeData                 )             , dimension(:), intent(inout) :: nodes
    type   (treeNodeList             ), allocatable, dimension(:), intent(inout) :: nodeList
    logical                           , allocatable, dimension(:), intent(inout) :: childIsSubhalo
    integer                                                      , intent(  out) :: isolatedNodeCount
    integer(c_size_t                 )                                           :: iNode
    integer                                                                      :: iIsolatedNode    , initialSatelliteCount
    logical                                                                      :: createNode
    type   (progenitorIterator       )                                           :: progenitors

    ! Determine how many nodes are isolated (i.e. not subhalos).
    isolatedNodeCount=count(.not.nodes%isSubhalo)

    ! Scan here for nodes that are subhalos and have no progenitor. These objects must be
    ! created as satellites within the tree.
    initialSatelliteCount=0
    do iNode=1,size(nodes)
       if (nodes(iNode)%isSubhalo) then
          call progenitors%descendantSet(self,nodes(iNode),nodes)
          if (.not.progenitors%exist()) initialSatelliteCount=initialSatelliteCount+1
      end if
    end do

    ! Allocate nodes.
    allocate(nodeList(isolatedNodeCount+initialSatelliteCount))
    allocate(childIsSubhalo(isolatedNodeCount+initialSatelliteCount))

    ! Create the nodes.
    iIsolatedNode          =0
    nodes%isolatedNodeIndex=nodeReachabilityUnreachable%ID
    do iNode=1,size(nodes)
       createNode=.false.
       if (nodes(iNode)%nodeIndex == nodes(iNode)%host%nodeIndex) then
          createNode=.true.
       else if (nodes(iNode)%isSubhalo) then
          call progenitors%descendantSet(self,nodes(iNode),nodes)
          if (.not.progenitors%exist()) createNode=.true.
       end if
       if (createNode) then
          iIsolatedNode=iIsolatedNode+1
          ! Store a record of where this node goes in the isolated node list.
          nodes(iNode)%isolatedNodeIndex=iIsolatedNode
          nodeList(iIsolatedNode)%node => treeNode(hostTree=tree)
          call nodeList(iIsolatedNode)%node%indexSet(nodes(iNode)%nodeIndex)
          nodes(iNode)%node => nodeList(iIsolatedNode)%node
       end if
    end do
    return
  end subroutine readCreateNodeArray

  subroutine readBuildIsolatedParentPointers(self,tree,nodes,nodeList)
    !!{
    Create parent pointer links between isolated nodes and assign times and masses to those nodes.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : mergerTree  , nodeComponentBasic, treeNodeList
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: String_Handling           , only : operator(//)
    implicit none
    class           (mergerTreeConstructorRead)                       , intent(inout)          :: self
    type            (mergerTree               )                       , intent(inout) , target :: tree
    type            (nodeData                 )         , dimension(:), intent(inout)          :: nodes
    type            (treeNodeList             )         , dimension(:), intent(inout)          :: nodeList
    class           (nodeComponentBasic       ), pointer                                       :: basic
    type            (mergerTree               ), pointer                                       :: treeCurrent
    type            (nodeData                 ), pointer                                       :: parentNode
    integer                                                                                    :: iNode
    integer         (c_size_t                 )                                                :: iIsolatedNode
    type            (varying_string           )                                                :: message
    character       (len=12                   )                                                :: label
    logical                                                                                    :: assignLastIsolatedTime

    do iNode=1,size(nodes)
       ! Only process if this is an isolated node (or an initial satellite).
       if (nodes(iNode)%isolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Check for an isolated node.
          if (nodes(iNode)%nodeIndex == nodes(iNode)%host%nodeIndex) then
             assignLastIsolatedTime=.false.
             if (associated(nodes(iNode)%parent)) then
                nodeList(iIsolatedNode)%node%parent  => nodes(iNode)%parent%node
             else
                nodeList(iIsolatedNode)%node%parent  => null()
                ! Find a tree to attach this base node to. Begin with the original tree passed to us.
                treeCurrent => tree
                ! Check if its nodeBase is already assigned.
                do while (associated(treeCurrent%nodeBase))
                   ! While it is, create the next tree (unless it already exists), then step to it.
                   if (.not.associated(treeCurrent%nextTree)) then
                      allocate(treeCurrent%nextTree)
                      allocate(treeCurrent%nextTree%randomNumberGenerator_,mold=self%randomNumberGenerator_)
                      !$omp critical(mergerTreeConstructReadDeepCopyReset)
                      !![
                      <deepCopyReset variables="self%randomNumberGenerator_"/>
                      <deepCopy source="self%randomNumberGenerator_" destination="treeCurrent%nextTree%randomNumberGenerator_"/>
                      <deepCopyFinalize variables="treeCurrent%nextTree%randomNumberGenerator_"/>
                      !!]
                      !$omp end critical(mergerTreeConstructReadDeepCopyReset)
                      call treeCurrent%nextTree%randomNumberGenerator_%seedSet(seed=treeCurrent%nextTree%index,offset=.true.)
                   end if
                   treeCurrent => treeCurrent%nextTree
                end do
                ! Assign this node as the base node of the current tree.
                treeCurrent   %firstTree         => tree
                treeCurrent   %nodeBase          => nodeList(iIsolatedNode)%node
                if (self%treeIndexToRootNodeIndex) then
                   treeCurrent%index             =  nodes   (iNode        )%nodeIndex
                else
                   treeCurrent%index             =  tree                   %index
                end if
                treeCurrent   %volumeWeight      =  self%treeWeightCurrent
                treeCurrent   %initializedUntil  =  0.0d0
                treeCurrent   %isTreeInitialized =  .false.
                treeCurrent   %event             => null()
                ! Initialize a new random number sequence for this tree, using the sum of the tree index and base node index as the seed increment.
                if (.not.associated(treeCurrent,tree)) call treeCurrent%randomNumberGenerator_%seedSet(seed=tree%index+nodes(iNode)%nodeIndex,offset=.true.)
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
                call Error_Report('initial satellite has no parent defined'//{introspection:location})
             end if
          end if
          ! Assign mass and time. For the case of satellites we also assign the time at
          ! which the satellite was last isolated. Since we do not know this, we simply
          ! set it equal to the current time (which is, obviously, an upper limit).
          if (nodes(iNode)%nodeMass <= 0.0d0) then
             write (label,'(e12.6)') nodes(iNode)%nodeMass
             message='non-positive mass ['//label//'] found for node '
             message=message//nodeList(iIsolatedNode)%node%index()
             call Error_Report(message//{introspection:location})
          end if
          if (nodes(iNode)%nodeTime <= 0.0d0) then
             write (label,'(e12.6)') nodes(iNode)%nodeTime
             message='non-positive time ['//label//'] found for node '
             message=message//nodeList(iIsolatedNode)%node%index()
             call Error_Report(message//{introspection:location})
          end if
          basic => nodeList(iIsolatedNode)%node%basic(autoCreate=.true.)
          call        basic%massSet            (nodes(iNode)%nodeMass)
          call        basic%timeSet            (nodes(iNode)%nodeTime)
          if (assignLastIsolatedTime) &
               & call basic%timeLastIsolatedSet(nodes(iNode)%nodeTime)
       end if
    end do
    return
  end subroutine readBuildIsolatedParentPointers

  subroutine readBuildChildAndSiblingLinks(nodes,nodeList,childIsSubhalo)
    !!{
    Build child and sibling links between nodes.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic, treeNodeList
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class  (nodeData          )             , dimension(:), intent(inout) :: nodes
    type   (treeNodeList      )             , dimension(:), intent(inout) :: nodeList
    logical                    , allocatable, dimension(:), intent(inout) :: childIsSubhalo
    class  (nodeComponentBasic), pointer                                  :: basic            , basicPrimary
    integer                                                               :: iNode
    integer(c_size_t          )                                           :: iIsolatedNode
    logical                                                               :: descendsToSubhalo

    childIsSubhalo=.false.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Check for an isolated node.
          if (nodes(iNode)%nodeIndex == nodes(iNode)%host%nodeIndex) then
             ! Check if the node has a parent.
             if (associated(nodeList(iIsolatedNode)%node%parent)) then
                ! Determine if this node definitely descends to a subhalo - in which case it can never be the primary progenitor.
                descendsToSubhalo=nodes(iNode)%descendantIndex /= nodeList(iIsolatedNode)%node%parent%index()
                ! It does, so set the child pointer of the parent appropriately.
                if (associated(nodeList(iIsolatedNode)%node%parent%firstChild)) then
                   ! A child is already associated. Check if current node does not descend to a subhalo and is more massive.
                   basic        => nodeList(iIsolatedNode)%node                  %basic()
                   basicPrimary => nodeList(iIsolatedNode)%node%parent%firstChild%basic()
                   if (.not.descendsToSubhalo                                          &
                        & .and. (                                                      &
                        &        childIsSubhalo(nodes(iNode)%parent%isolatedNodeIndex) &
                        &         .or.                                                 &
                        &        basic%mass() > basicPrimary%mass()                    &
                        &       )                                                      &
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
    deallocate(childIsSubhalo)
    return
  end subroutine readBuildChildAndSiblingLinks

  subroutine readAssignHostTreePointers(tree)
    !!{
    After tree base nodes have been assigned, walk each tree and set the host tree pointer for each node.
    !!}
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    implicit none
    type(mergerTree              ), intent(inout), target  :: tree
    type(mergerTree              )               , pointer :: treeCurrent
    type(treeNode                )               , pointer :: node
    type(mergerTreeWalkerAllNodes)                         :: treeWalker

    treeCurrent => tree
    do while(associated(treeCurrent))
       treeWalker=mergerTreeWalkerAllNodes(treeCurrent,spanForest=.false.)
       do while (treeWalker%next(node))
          node%hostTree => treeCurrent
       end do
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine readAssignHostTreePointers

  subroutine readAssignScaleRadii(self,nodes,nodeList)
    !!{
    Assign scale radii to nodes.
    !!}
    use :: Display                   , only : displayIndent            , displayMessage                , displayUnindent              , displayVerbosity, &
          &                                   verbosityLevelWarn       , enumerationVerbosityLevelType
    use :: Error                     , only : Error_Report             , errorStatusSuccess
    use :: Galacticus_Nodes          , only : nodeComponentBasic       , nodeComponentDarkMatterProfile, treeNodeList
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: Root_Finder               , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative , rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (mergerTreeConstructorRead        ), target                       , intent(inout) :: self
    class           (nodeData                         )                 , dimension(:), intent(inout) :: nodes
    type            (treeNodeList                     )                 , dimension(:), intent(inout) :: nodeList
    double precision                                   , parameter                                    :: scaleRadiusMaximumAllowed  =100.0d0, toleranceAbsolute  =1.0d-9, &
         &                                                                                               toleranceRelative          =1.0d-9
    logical                                                       , save                              :: excessiveScaleRadiiReported=.false.
    class           (nodeComponentBasic               ), pointer                                      :: basic
    class           (nodeComponentDarkMatterProfile   ), pointer                                      :: darkMatterProfile
    integer                                                                                           :: iNode                              , status
    type             (enumerationVerbosityLevelType   )                                               ::  messageVerbosity
    integer         (c_size_t                         )                                               :: iIsolatedNode
    double precision                                                                                  :: radiusScale
    logical                                                                                           :: excessiveHalfMassRadii             , excessiveScaleRadii       , &
         &                                                                                               useFallbackScaleMethod
    type            (rootFinder                       )           , save                              :: finder
    logical                                                       , save                              :: finderConstructed          =.false.
    !$omp threadprivate(finder,finderConstructed)
    type            (varying_string                   )                                               :: message
    character       (len=16                           )                                               :: label

    ! Initialize our root finder.
    if (.not.finderConstructed) then
       finder           =rootFinder(                                          &
            &                       rootFunction     =readRadiusHalfMassRoot, &
            &                       toleranceAbsolute=toleranceAbsolute     , &
            &                       toleranceRelative=toleranceRelative       &
            &                      )
       finderConstructed=.true.
    end if
    self_ => self
    ! Find the scale radius.
    excessiveScaleRadii   =.false.
    excessiveHalfMassRadii=.false.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Assume that we need to use a fallback method to set halo scale radius.
          useFallbackScaleMethod=.true.
          ! Check if the node is sufficiently massive.
          basic             => nodeList(iIsolatedNode)%node%basic            (                 )
          darkMatterProfile => nodeList(iIsolatedNode)%node%darkMatterProfile(autoCreate=.true.)
          if (basic%mass() >= self%presetScaleRadiiMinimumMass) then
             ! Check if we have scale radii read directly from file.
             if     (                                                                        &
                  &     nodes(iNode)%scaleRadius                                             &
                  &   >                                                                      &
                  &     0.0d0                                                                &
                  &  .and.                                                                   &
                  &     nodes(iNode)%scaleRadius                                             &
                  &   <                                                                      &
                  &     self%darkMatterHaloScale_%radiusVirial(nodeList(iIsolatedNode)%node) &
                  &    /self%presetScaleRadiiConcentrationMinimum                            &
                  &  .and.                                                                   &
                  &     nodes(iNode)%scaleRadius                                             &
                  &   >                                                                      &
                  &     self%darkMatterHaloScale_%radiusVirial(nodeList(iIsolatedNode)%node) &
                  &    /self%presetScaleRadiiConcentrationMaximum                            &
                  & ) then
                ! We do, so simply use them to set the scale radii in tree nodes.
                call darkMatterProfile%scaleSet(nodes(iNode)%scaleRadius)
                useFallbackScaleMethod=.false.
             else if (nodes(iNode)%halfMassRadius > 0.0d0) then
                ! We do not have scale radii read directly. Instead, compute them from half-mass radii.
                ! Set the active node and target half mass radius.
                node_              => nodeList(iIsolatedNode)%node
                darkMatterProfile_ => node_                  %darkMatterProfile()
                basic_             => node_                  %basic            ()
                radiusHalfMass_    =  nodes   (iNode        )%halfMassRadius
                ! Solve for the scale radius.
                call finder%rangeExpand    (                                                                              &
                     &                      rangeExpandDownward          =0.5d0                                         , &
                     &                      rangeExpandUpward            =2.0d0                                         , &
                     &                      rangeDownwardLimit           = self%darkMatterHaloScale_%radiusVirial(node_)  &
                     &                                                    /self%presetScaleRadiiConcentrationMaximum    , &
                     &                      rangeUpwardLimit             = self%darkMatterHaloScale_%radiusVirial(node_)  &
                     &                                                    /self%presetScaleRadiiConcentrationMinimum    , &
                     &                      rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive                 , &
                     &                      rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative                 , &
                     &                      rangeExpandType              =rangeExpandMultiplicative                       &
                     &                     )
                radiusScale=finder%find(rootGuess=radiusHalfMass_,status=status)
                if (status == errorStatusSuccess) then
                   call darkMatterProfile%scaleSet(radiusScale)
                   ! Check for scale radii exceeding the virial radius.
                   if (radiusScale     > self%darkMatterHaloScale_%radiusVirial(node_)) excessiveScaleRadii   =.true.
                   ! Check for half-mass radii exceeding the virial radius.
                   if (radiusHalfMass_ > self%darkMatterHaloScale_%radiusVirial(node_)) excessiveHalfMassRadii=.true.
                   useFallbackScaleMethod=.false.
                else
                   if (self%scaleRadiiFailureIsFatal) then
                      messageVerbosity=displayVerbosity()
                   else
                      messageVerbosity=verbosityLevelWarn
                   end if
                   call displayIndent  ("failed to find scale radius consistent with specified half-mass radius",messageVerbosity)
                   write (label,'(i16)') node_%hostTree%index
                   message="      tree index: "//trim(label)
                   call displayMessage(message,messageVerbosity)
                    write (label,'(i16)') node_%index()
                   message="      node index: "//trim(label)
                   call displayMessage(message,messageVerbosity)
                   write (label,'(e12.6)') self%darkMatterHaloScale_%radiusVirial(node_)
                   message="   virial radius: "//trim(label)
                   call displayMessage(message,messageVerbosity)
                   write (label,'(e12.6)') radiusHalfMass_
                   message="half-mass radius: "//trim(label)
                   call displayMessage(message,messageVerbosity)
                   call displayUnindent("",messageVerbosity)
                   if (self%scaleRadiiFailureIsFatal) then
                      call Error_Report('problem with half-mass radii - see report above - aborting'//{introspection:location})
                   else
                     useFallbackScaleMethod=.true.
                   end if
                end if
             end if
          end if
          if (useFallbackScaleMethod) then
             ! The node mass is below the reliability threshold, or no scale information is available. Set the scale radius using
             ! the fallback method.
             node_ => nodeList(iIsolatedNode)%node
             radiusScale=max(                                                                                                      &
                  &          min(                                                                                                  &
                  &              self%darkMatterProfileScaleRadius_%radius      (node_)                                          , &
                  &              self%darkMatterHaloScale_         %radiusVirial(node_)/self%presetScaleRadiiConcentrationMinimum  &
                  &             )                                                                                                , &
                  &              self%darkMatterHaloScale_         %radiusVirial(node_)/self%presetScaleRadiiConcentrationMaximum  &
                  &         )
             call darkMatterProfile%scaleSet(radiusScale)
          end if
       end if
    end do
    ! Report warning on excessive scale radii if not already done.
    if (excessiveScaleRadii.and..not.excessiveScaleRadiiReported) then
       excessiveScaleRadiiReported=.true.
       call displayMessage('warning - some scale radii exceed the corresponding virial radii - suggests&
            & inconsistent definitions of halo mass/radius'//{introspection:location},verbosityLevelWarn)
    end if

    ! Exit on excessive half mass radii.
    if (excessiveHalfMassRadii) call Error_Report('some half mass radii exceed corresponding virial radii'//{introspection:location})

    return
  end subroutine readAssignScaleRadii

  subroutine readAssignAngularMomenta(self,nodes,nodeList)
    !!{
    Assign angular momenta to nodes.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : nodeComponentBasic, nodeComponentSpin, treeNodeList
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: Mass_Distributions        , only : massDistributionClass
  implicit none
    class           (mergerTreeConstructorRead)                       , intent(inout) :: self
    class           (nodeData                 )         , dimension(:), intent(inout) :: nodes
    type            (treeNodeList             )         , dimension(:), intent(inout) :: nodeList
    class           (nodeComponentBasic       ), pointer                              :: basic
    class           (nodeComponentSpin        ), pointer                              :: spin
    class           (massDistributionClass    ), pointer                              :: massDistribution_
    integer                                                                           :: iNode
    integer         (c_size_t                 )                                       :: iIsolatedNode
    double precision                                                                  :: angularMomentum  , radiusVirial
    double precision                                    , dimension(3)                :: angularMomentum3D

    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Get basic and spin components.
          basic             =>                                         nodeList(iIsolatedNode)%node%basic(                 )
          spin              =>                                         nodeList(iIsolatedNode)%node%spin (autoCreate=.true.)
          radiusVirial      =  self%darkMatterHaloScale_ %radiusVirial(nodeList(iIsolatedNode)%node                         )
          massDistribution_ => self%darkMatterProfileDMO_%get         (nodeList(iIsolatedNode)%node                         )
          if (self%presetAngularMomenta  ) then
             if      (self%mergerTreeImporter_%angularMomentaAvailable()) then
                ! If angular momenta are available directly, use them.
                call spin%angularMomentumSet(nodes(iNode)%angularMomentum)
             else if (self%mergerTreeImporter_%          spinAvailable()) then
                angularMomentum=+nodes(iNode)%spin   &
                     &          *spinNormalization()
                call spin%angularMomentumSet(angularMomentum)
             else
                call Error_Report('no method exists to set angular momenta'//{introspection:location})
             end if
             if (self%presetUnphysicalAngularMomenta.and.spin%angularMomentum() <= 0.0d0) &
                  & call spin%angularMomentumSet(self%haloSpinDistribution_%sample(nodeList(iIsolatedNode)%node)*spinNormalization())
          end if
          if (self%presetAngularMomenta3D) then
             if      (self%mergerTreeImporter_%angularMomenta3DAvailable()) then
                ! If angular momenta are available directly, use them.
                call spin%angularMomentumVectorSet(nodes(iNode)%angularMomentum3D)
             else if (self%mergerTreeImporter_%          spin3DAvailable()) then
                angularMomentum3D=+nodes(iNode)%spin3D   &
                     &            *spinNormalization()
                call spin%angularMomentumVectorSet(angularMomentum3D)
             else
                call Error_Report('no method exists to set vector angular momenta'//{introspection:location})
             end if
          end if
          !![
	  <objectDestructor name="massDistribution_"/>
	  !!]
       end if
    end do
    return

  contains

    double precision function spinNormalization()
      !!{
      Normalization for conversion of spin to angular momentum.
      !!}
      use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
      implicit none
      spinNormalization=+gravitationalConstant_internal                                      &
           &            *basic%mass()**2.5d0                                                 &
           &            /sqrt(abs(massDistribution_%energy(radiusVirial,massDistribution_)))
      return
    end function spinNormalization

  end subroutine readAssignAngularMomenta

  subroutine readAssignNamedProperties(self,nodes,nodeList)
    !!{
    Assign named properties to nodes.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic, treeNodeList
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class  (mergerTreeConstructorRead)                       , intent(inout) :: self
    class  (nodeData                 )         , dimension(:), intent(inout) :: nodes
    type   (treeNodeList             )         , dimension(:), intent(inout) :: nodeList
    class  (nodeComponentBasic       ), pointer                              :: basic
    integer                                                                  :: iNode        , i
    integer(c_size_t                 )                                       :: iIsolatedNode

    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Assign the named properties.
          basic => nodeList(iIsolatedNode)%node%basic()
          if (size(self%presetNamedReals   ) > 0) then
             do i=1,size(self%presetNamedReals   )
                call basic%floatRank0MetaPropertySet      (self%indexNamedReals   (i),nodes(iNode)%reals   (i))
             end do
          end if
          if (size(self%presetNamedIntegers) > 0) then
             do i=1,size(self%presetNamedIntegers)
                call basic%longIntegerRank0MetaPropertySet(self%indexNamedIntegers(i),nodes(iNode)%integers(i))
             end do
          end if
        end if
    end do
    return
  end subroutine readAssignNamedProperties

  double precision function readRadiusHalfMassRoot(radius)
    !!{
    Function used to find scale radius of dark matter halos given their half-mass radius.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Mass_Distributions , only : massDistributionClass
    implicit none
    double precision                       , intent(in   ) :: radius
    class           (massDistributionClass), pointer       :: massDistribution_

    ! Set scale radius to current guess.
    call darkMatterProfile_%scaleSet(radius)
    call Calculations_Reset(node_)
    ! Compute difference between mass fraction enclosed at half mass radius and one half.
    massDistribution_      => self_            %darkMatterProfileDMO_%get(node_          )
    readRadiusHalfMassRoot =  massDistribution_%massEnclosedBySphere     (radiusHalfMass_)/basic_%mass()-0.5d0
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    return
  end function readRadiusHalfMassRoot

  subroutine readAssignIsolatedNodeIndices(nodes)
    !!{
    Assign to each node the number of the corresponding isolated node.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class  (nodeData), dimension(:), intent(inout) :: nodes
    class  (nodeData), pointer                     :: node
    integer(c_size_t)                              :: iIsolatedNode
    integer                                        :: iNode
    logical                                        :: endOfBranch

    ! First make a copy of the currently assigned isolated node indices. These will be used
    ! later to reference the nodes which are the primary node associated with objects in nodeList.
    nodes%primaryIsolatedNodeIndex=nodes%isolatedNodeIndex
    ! Iterate over nodes.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
          iIsolatedNode=nodes(iNode)%isolatedNodeIndex
          ! Find the subset with descendants.
          if (associated(nodes(iNode)%descendant)) then
             ! Select the subset which have a subhalo as a descendant.
             if (nodes(iNode)%descendant%isSubhalo) then
                ! Trace descendants until merging or final time.
                node        => nodes(iNode)%descendant
                endOfBranch =  .false.
                do while (.not.endOfBranch)
                   ! Record that this node was reachable via descendants of an isolated node.
                   if (node%isolatedNodeIndex == nodeReachabilityUnreachable%ID) node%isolatedNodeIndex=nodeReachabilityReachable%ID
                   if (.not.associated(node%descendant)) then
                      ! If there is no descendant then the end of the branch has been reached.
                      endOfBranch=.true.
                   else
                      ! Step to the next descendant.
                      node => node%descendant
                   end if
                end do
             end if
          end if
       end if
    end do
    return
  end subroutine readAssignIsolatedNodeIndices

  subroutine readScanForMergers(self,nodes,nodeList,historyCountMaximum)
    !!{
    Scan for and record mergers between nodes.
    !!}
    use :: Galacticus_Nodes            , only : nodeComponentBasic        , nodeComponentPosition, nodeComponentSatellite, treeNode, &
         &                                      treeNodeList
    use :: Merger_Tree_Read_Importers  , only : nodeData
    use :: Satellite_Merging_Timescales, only : satelliteMergeTimeInfinite
    implicit none
    class           (mergerTreeConstructorRead)                         , intent(inout) :: self
    class           (nodeData                 ), target   , dimension(:), intent(inout) :: nodes
    type            (treeNodeList             )           , dimension(:), intent(inout) :: nodeList
    integer         (c_size_t                 )                         , intent(  out) :: historyCountMaximum
    class           (nodeData                 ), pointer                                :: lastSeenNode                , progenitorNode            , &
         &                                                                                 node
    type            (treeNode                 ), pointer                                :: firstProgenitor             , hostNode                  , &
         &                                                                                 orbitalPartner              , satelliteNode
    integer                                    , parameter                              :: passAssign                =1, passMerge               =2
    class           (nodeComponentBasic       ), pointer                                :: basicChild                  , basic
    class           (nodeComponentSatellite   ), pointer                                :: satellite
    class           (nodeComponentPosition    ), pointer                                :: position                    , positionChild
    integer                                                                             :: iNode                       , pass_
    integer         (c_size_t                 )                                         :: historyCount                , iIsolatedNode
    integer         (kind_int8                )                                         :: progenitorMassMaximumIndex
    logical                                                                             :: branchMerges                , branchTipReached          , &
         &                                                                                 endOfBranch                 , isolatedProgenitorExists  , &
         &                                                                                 nodeWillMerge               , nodeIsMostMassive
    double precision                                                                    :: timeSubhaloMerges           , progenitorMassMaximum
    type            (progenitorIterator       )                                         :: progenitors

    ! Initialize.
    historyCountMaximum  = 0
    nodes%mergesWithIndex=-1
    ! First pass assigns isolated node indices to all descendants, second pass finds mergers.
    do pass_=passAssign,passMerge
       do iNode=1,size(nodes)
          if (nodes(iNode)%primaryIsolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
             iIsolatedNode=nodes(iNode)%primaryIsolatedNodeIndex
             ! Find the subset with descendants.
             if (associated(nodes(iNode)%descendant)) then
                nodeIsMostMassive=.true.
                ! Determine if this is the most massive progenitor.
                if (self%alwaysPromoteMostMassive) then
                   call progenitors%descendantSet(self,nodes(iNode)%descendant,nodes)
                   do while (progenitors%next(nodes))
                      progenitorNode => progenitors%current(nodes)
                      if (progenitorNode%nodeMass > nodes(iNode)%nodeMass) nodeIsMostMassive=.false.
                   end do
                end if
                ! Flag indicating if this is a node for which a merging time should be set.
                nodeWillMerge=.false.
                ! Select the subset which have a subhalo descendant, or which are an initial subhalo.
                if (nodes(iNode)%descendant%isSubhalo.or.nodes(iNode)%isSubhalo) then
                   ! Trace descendants until merging or final time.
                   endOfBranch     =.false.
                   branchTipReached=.false.
                   branchMerges    =.false.
                   historyCount    =0
                   if (nodes(iNode)%isSubhalo) then
                      node => nodes(iNode)
                   else
                      ! Check for an immediate subhalo-subhalo merger.
                      if (self%isSubhaloSubhaloMerger(nodes,nodes(iNode))) then
                         endOfBranch =.true.
                         branchMerges=.true.
                         nodes(iNode)%mergesWithIndex=nodes(iNode)%descendant%nodeIndex
                         historyCount=historyCount+max(0_kind_int8,self%mergerTreeImporter_%subhaloTraceCount(nodes(iNode)))
                      end if
                      lastSeenNode => nodes(iNode)
                      node         => nodes(iNode)%descendant
                   end if
                   do while (.not.endOfBranch)
                      ! Record which isolated node this node belongs to.
                      node%isolatedNodeIndex=iIsolatedNode
                      ! Increment the history count for this branch.
                      historyCount=historyCount+1
                      ! Test the branch.
                      if (.not.associated(node%descendant)) then
                         ! No descendant, indicating tip of branch has been reached
                         branchTipReached            =.true.
                         endOfBranch                 =.true.
                         historyCount                =historyCount+max(0_kind_int8,self%mergerTreeImporter_%subhaloTraceCount(node))
                      else if (.not.node%descendant%isSubhalo) then
                         ! Descendant is not a subhalo, treat as a merging event or a subhalo promotion.
                         endOfBranch                 =.true.
                         historyCount                =historyCount+max(0_kind_int8,self%mergerTreeImporter_%subhaloTraceCount(node))
                         ! Search for any isolated progenitors of the node's descendant.
                         isolatedProgenitorExists=.false.
                         call progenitors%descendantSet(self,node%descendant,nodes)
                         progenitorMassMaximum=-1.0d0
                         progenitorMassMaximumIndex=-1_kind_int8
                         do while (progenitors%next(nodes) .and. .not.isolatedProgenitorExists)
                            progenitorNode => progenitors%current(nodes)
                            isolatedProgenitorExists=(progenitorNode%nodeIndex == progenitorNode%hostIndex)
                            if (progenitorNode%nodeMass > progenitorMassMaximum) then
                               progenitorMassMaximum     =progenitorNode%nodeMass
                               progenitorMassMaximumIndex=progenitorNode%nodeIndex
                            end if
                         end do
                         ! If an isolated progenitor exists, or this is not the most massive subhalo progenitor, this is a merger
                         ! event. If not, it is a subhalo promotion (which will be handled elsewhere).
                         if (isolatedProgenitorExists .or. progenitorMassMaximumIndex /= node%nodeIndex) then
                            branchMerges                =.true.
                            nodes(iNode)%mergesWithIndex=node%descendant%nodeIndex
                            lastSeenNode                => node
                            node                        => node%descendant
                         end if
                      else
                         ! Merges with another subhalo.
                         call progenitors%descendantSet(self,node%descendant,nodes)
                         do while (progenitors%next(nodes))
                            progenitorNode => progenitors%current(nodes)
                            if     (                                                                                       &
                                 &                     progenitorNode%nodeIndex         /= node%nodeIndex                  &
                                 &  .and.              progenitorNode%isolatedNodeIndex /= nodeReachabilityUnreachable%ID  &
                                 &  .and.   associated(progenitorNode%descendant                                         ) &
                                 &  .and.massIsGreater(progenitorNode                   ,  node                          ) &
                                 & ) then
                               ! Another node merges into current node's descendant subhalo and is more massive than current
                               ! node. Therefore, class this as a subhalo-subhalo merger.
                               branchMerges                =  .true.
                               endOfBranch                 =  .true.
                               nodes(iNode)%mergesWithIndex=  progenitorNode%descendant%nodeIndex
                               historyCount                =  historyCount+max(0_kind_int8,self%mergerTreeImporter_%subhaloTraceCount(node))
                               lastSeenNode                => node
                               node                        => node%descendant
                               exit
                            end if
                         end do
                         ! Step to the next descendant.
                         if (.not.endOfBranch) node => node%descendant
                      end if
                   end do
                   ! If on the isolated node index assigning pass, skip to the next halo.
                   if (pass_ == passAssign) cycle
                   ! Only set a merging time if this node is not the primary progenitor of its parent.
                   if (.not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                      ! Record the largest history.
                      historyCountMaximum=max(historyCountMaximum,historyCount)
                      ! Set an appropriate merging time for this subhalo.
                      if      (branchTipReached) then
                         timeSubhaloMerges=satelliteMergeTimeInfinite ! Subhalo never merges, so set merging time to effective infinity.
                      else if (branchMerges    ) then
                         ! Find the time of merging, accounting for any additional (subresolution) time.
                         timeSubhaloMerges=node%nodeTime
                         call self%timeUntilMergingSubresolution(lastSeenNode,nodes,nodeList,iNode,timeSubhaloMerges)
                      else
                         ! Neither the branch tip was reached, not does this branch merge. Therefore, this must be a subhalo which is
                         ! promoted to be an isolated halo. Simply set an infinite merging time as we do not wish this node to merge.
                         timeSubhaloMerges=satelliteMergeTimeInfinite
                      end if
                      ! Flag that this node will merge.
                      nodeWillMerge=.true.
                   end if
                else if (pass_ == passAssign) then
                   ! If on the isolated node index assigning pass, skip to the next halo.
                   cycle
                else if (.not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor() .or. ( self%alwaysPromoteMostMassive .and. .not.nodeIsMostMassive )) then
                   ! Descendant is not a subhalo but this node is not the primary progenitor. Or we are always promoting the most
                   ! massive progenitor, even if it is a subhalo and there are isolated progenitors present. Assume the halo
                   ! merges instantaneously at the time of merging with its parent halo. (That merging event occurs at the time of
                   ! the parent halo.)
                   basic             => nodeList(iIsolatedNode)%node%parent%basic()
                   timeSubhaloMerges =  basic                              %time ()
                   ! Flag that this node will merge.
                   nodeWillMerge=.true.
                   ! Record the node with which the merger occurs.
                   nodes(iNode)%mergesWithIndex=nodes(iNode)%descendant%nodeIndex
                   ! Ensure the history arrays will be large enough to hold data for this node.
                   historyCountMaximum=max(historyCountMaximum,max(0_kind_int8,self%mergerTreeImporter_%subhaloTraceCount(nodes(iNode))))
                   ! Account for any subresolution merging time.
                   call self%timeUntilMergingSubresolution(nodes(iNode),nodes,nodeList,iNode,timeSubhaloMerges)
                end if
                ! Set a merging time and/or orbit if this node will merge.
                if (self%presetMergerTimes) then
                   ! If the node does not merge set an infinite merging time.
                   if (.not.nodeWillMerge) timeSubhaloMerges=satelliteMergeTimeInfinite
                   ! Store the time of merging for this node and all of its primary progenitors.
                   firstProgenitor => nodeList(iIsolatedNode)%node
                   do while (associated(firstProgenitor))
                      satellite => firstProgenitor%satellite(autoCreate=.true.)
                      call satellite%timeOfMergingSet(timeSubhaloMerges)
                      firstProgenitor => firstProgenitor%firstChild
                   end do
                end if
             else
                ! Node has no descendant - it must therefore be an initial subhalo seen only once. A position history must be set
                ! for this, so ensure the history arrays are sufficiently sized.
                historyCountMaximum=max(historyCountMaximum,max(0_kind_int8,self%mergerTreeImporter_%subhaloTraceCount(nodes(iNode))))
             end if
             ! Handle cases where a node jumps to another tree.
             if (pass_ == passMerge .and. self%isOnPushList(nodes(iNode)) .and. self%presetMergerTimes) then
                ! Merger times are to be preset, but this node will be pushed to another tree. We must set its merging time to be
                ! infinite in this case.
                firstProgenitor => nodeList(iIsolatedNode)%node
                do while (associated(firstProgenitor))
                   satellite => firstProgenitor%satellite(autoCreate=.true.)
                   call satellite%timeOfMergingSet(satelliteMergeTimeInfinite)
                   firstProgenitor => firstProgenitor%firstChild
                end do
             end if
             ! Set position and velocity if required.
             if (self%presetPositions) then
                position => nodeList(iIsolatedNode)%node%position(autoCreate=.true.)
                call position%positionSet(nodes(iNode)%position)
                call position%velocitySet(nodes(iNode)%velocity)
                ! Detect if the node parent has no isolated child - in which case one will have been made for it using a
                ! direct copy of itself. Note that this is an ugly solution - once trees can handle nodes with no primary
                ! progenitor (but with secondary progenitors) a cleaner test could be used here.
                if (associated(nodeList(iIsolatedNode)%node%firstChild)) then
                   basic      => nodeList(iIsolatedNode)%node           %basic()
                   basicChild => nodeList(iIsolatedNode)%node%firstChild%basic()
                   if (nodeList(iIsolatedNode)%node%uniqueID() == nodeList(iIsolatedNode)%node%firstChild%uniqueID()) then
                      ! Set the position and velocity of the pseudo-primary progenitor here also.
                      positionChild => nodeList(iIsolatedNode)%node%firstChild%position(autoCreate=.true.)
                      call positionChild%positionSet(nodes(iNode)%position)
                      call positionChild%velocitySet(nodes(iNode)%velocity)
                   end if
                end if
             end if
          end if
       end do
    end do
    ! Set orbits.
    if (self%presetOrbits) then
       iIsolatedNode=0
       do iNode=1,size(nodes)
          if (nodes(iNode)%primaryIsolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
             iIsolatedNode=nodes(iNode)%primaryIsolatedNodeIndex
             ! Find the most massive progenitor.
             nodeIsMostMassive=.true.
             if (self%alwaysPromoteMostMassive) then
                if (self%alwaysPromoteMostMassive.and.associated(nodes(iNode)%descendant)) then
                   call progenitors%descendantSet(self,nodes(iNode)%descendant,nodes)
                   do while (progenitors%next(nodes))
                      progenitorNode => progenitors%current(nodes)
                      if (progenitorNode%nodeMass > nodes(iNode)%nodeMass) nodeIsMostMassive=.false.
                   end do
                end if
             end if
             ! Set the orbit for this halo if required. For this to be required we must first have a parent node. Then, we must be
             ! either not the primary progenitor, or there must be a more massive subhalo which will be promoted.
             satelliteNode => nodeList(iIsolatedNode)%node
             if (associated(satelliteNode%parent) .and. (.not.satelliteNode%isPrimaryProgenitor() .or. (self%alwaysPromoteMostMassive .and. .not.nodeIsMostMassive))) then
                ! Find the orbital partner.
                hostNode => satelliteNode%parent%firstChild
                ! If the parent node has no progenitors, then we are forced to use the parent node itself as the orbital partner.
                if (associated(hostNode)) then
                   orbitalPartner => hostNode     %parent
                else
                   hostNode       => satelliteNode%parent
                   orbitalPartner => hostNode
                end if
                call self%setOrbit(satelliteNode,hostNode,orbitalPartner)
             end if
          end if
       end do
    end if
    return
  end subroutine readScanForMergers

  subroutine readSetOrbit(self,satelliteNode,hostNode,orbitalPartner)
    !!{
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentPosition, nodeComponentSatellite, treeNode
    use :: Kepler_Orbits   , only : keplerOrbit
    use :: Vectors         , only : Vector_Magnitude
    use :: String_Handling , only : operator(//)
    implicit none
    class           (mergerTreeConstructorRead), intent(inout)           :: self
    type            (treeNode                 ), intent(inout)           :: satelliteNode              , hostNode           , &
         &                                                                  orbitalPartner
    class           (nodeComponentBasic       ), pointer                 :: basicOrbitalPartner        , basicSatellite
    class           (nodeComponentPosition    ), pointer                 :: positionSatellite          , positionHost
    class           (nodeComponentSatellite   ), pointer                 :: satelliteSatellite
    type            (keplerOrbit              )                          :: orbit
    logical                                    , parameter               :: acceptUnboundOrbits=.false.
    double precision                                      , dimension(3) :: hostPosition               , relativePosition   , &
         &                                                                  satellitePosition
    double precision                                      , dimension(3) :: hostVelocity               , relativeVelocity   , &
         &                                                                  satelliteVelocity
    double precision                                                     :: radiusApocenter            , radiusPericenter   , &
         &                                                                  radiusVirial
    type            (varying_string           )                          :: message

    ! Get components.
    basicSatellite      => satelliteNode    %basic    (                 )
    positionSatellite   => satelliteNode    %position (                 )
    satelliteSatellite  => satelliteNode    %satellite(autoCreate=.true.)
    positionHost        => hostNode         %position (                 )
    basicOrbitalPartner => orbitalPartner   %basic    (                 )
    ! Get position and velocity.
    satellitePosition   =  positionSatellite%position (                 )
    satelliteVelocity   =  positionSatellite%velocity (                 )
    hostPosition        =       positionHost%position (                 )
    hostVelocity        =       positionHost%velocity (                 )
    ! Find relative position and velocity.
    relativePosition=satellitePosition-hostPosition
    relativeVelocity=satelliteVelocity-hostVelocity
    ! Update position/velocity for periodicity and Hubble flow.
    call self%phaseSpacePositionRealize(basicSatellite%time(),relativePosition,relativeVelocity)
    ! Catch zero separation halos.
    if (Vector_Magnitude(relativePosition) == 0.0d0) then
       if (self%presetOrbitsSetAll) then
          ! The satellite and host have zero separation, so no orbit can be
          ! computed. Since all orbits must be set, choose an orbit at random.
          orbit=self%virialOrbit_%orbit(satelliteNode,hostNode,acceptUnboundOrbits)
          call satelliteSatellite%virialOrbitSet(orbit)
       else
          message='merging halos ['
          message=message//satelliteNode%index()//' and '//hostNode%index()//'] have zero separation'
          call Error_Report(message//{introspection:location})
       end if
    else
       ! Create the orbit.
       orbit=readOrbitConstruct(basicSatellite%mass(),basicOrbitalPartner%mass(),relativePosition,relativeVelocity)
       ! Propagate to the virial radius.
       radiusPericenter=orbit                    %radiusPericenter(              )
       radiusApocenter =orbit                    %radiusApocenter (              )
       radiusVirial    =self%darkMatterHaloScale_%radiusVirial    (orbitalPartner)
       ! Check if the orbit intersects the virial radius.
       if     (                                                              &
            &    radiusVirial >= radiusPericenter                            &
            &  .and.                                                         &
            &   (radiusVirial <= radiusApocenter  .or. .not.orbit%isBound()) &
            &  .and.                                                         &
            &   (.not.self%presetOrbitsBoundOnly  .or.      orbit%isBound()) &
            & ) then
          call orbit%propagate(radiusVirial,infalling=.true.)
          ! Set the orbit.
          call satelliteSatellite%virialOrbitSet(orbit)
          ! If the satellite component supports full phase-space position, set that
          ! also.
          if (satelliteSatellite%positionIsSettable()) call satelliteSatellite%positionSet(relativePosition)
          if (satelliteSatellite%velocityIsSettable()) call satelliteSatellite%velocitySet(relativeVelocity)
       else if (self%presetOrbitsSetAll) then
          ! The given orbit does not cross the virial radius. Since all orbits must be set, choose an orbit at random.
          orbit=self%virialOrbit_%orbit(satelliteNode,hostNode,acceptUnboundOrbits)
          call satelliteSatellite%virialOrbitSet(orbit)
       else if (self%presetOrbitsAssertAllSet) then
          message='virial orbit could not be set for node '
          message=message//satelliteNode%index()//char(10)
          message=message//' -> set [presetOrbitsAssertAllSet]=false to ignore this problem'//char(10)
          message=message//'    (this may lead to other problems)'
          call Error_Report(message//{introspection:location})
       end if
        end if
    return
  end subroutine readSetOrbit

  logical function massIsGreater(node1,node2)
    !!{
    Return true if the mass of {\normalfont \ttfamily node1} is greater than that of {\normalfont \ttfamily node2}. In cases of
    precisely equal masses the tie is broken by considering the node indices (which should never be equal). This ensures that
    there is always a well-defined primary progenitor halo for example.
    !!}
    implicit none
    class(nodeData), intent(in   ) :: node1, node2
    
    massIsGreater=   node1%nodeMass  >  node2%nodeMass  &
         &        .or.                                  &
         &         (                                    &
         &           node1%nodeMass  == node2%nodeMass  &
         &          .and.                               &
         &           node1%nodeIndex >  node2%nodeIndex &
         &         )
    return
  end function massIsGreater

  subroutine readAssignMergers(self,nodes,nodeList)
    !!{
    Assign pointers to merge targets.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : treeNodeList
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: String_Handling           , only : operator(//)
    implicit none
    class  (mergerTreeConstructorRead)              , intent(inout) :: self
    class  (nodeData                 ), dimension(:), intent(inout) :: nodes
    type   (treeNodeList             ), dimension(:), intent(inout) :: nodeList
    type   (treeNode                 ), pointer                     :: nodeMergeTarget
    integer                                                         :: iNode
    integer(c_size_t                 )                              :: jNode
    type   (varying_string           )                              :: message

    if (self%presetMergerNodes) then
       do iNode=1,size(nodes)
          ! Check if this node was flagged as merging with another node.
          if (nodes(iNode)%mergesWithIndex /= -1) then
             ! Search for the node that it merges with.
             jNode=self%nodeLocation(nodes(iNode)%mergesWithIndex)
             if (nodes(jNode)%isolatedNodeIndex <= 0) then
                ! This node does not belong to any isolated halo - this should not happen.
                message='subhalo-subhalo ['
                message=message//nodes(iNode)%nodeIndex//":"//nodes(jNode)%nodeIndex
                message=message//'] merger in which subhalo has no isolated node progenitor - this should not happen'
                call Error_Report(message//{introspection:location})
             else
                ! Check for a cloned copy.
                nodeMergeTarget => nodeList(nodes(jNode)%isolatedNodeIndex)%node
                if     (                                                               &
                     &   associated(nodeMergeTarget%firstChild)                        &
                     &  .and.                                                          &
                     &   nodeMergeTarget%index() == nodeMergeTarget%firstChild%index() &
                     & ) nodeMergeTarget => nodeMergeTarget%firstChild
                ! Set pointer from merging node (a.k.a. the "mergee") to node that will be merged with.
                nodeList(nodes(iNode)%isolatedNodeIndex)%node%mergeTarget => nodeMergeTarget
                ! Make a backward pointer from the merge target to the mergee. Check if the target already has mergees associated with it.
                if (associated(nodeMergeTarget%firstMergee)) then
                   ! It does: unlink them and attached to the "siblingMergee" pointer of the current mergee.
                   nodeList(nodes(iNode)%isolatedNodeIndex)%node%siblingMergee => nodeMergeTarget%firstMergee
                else
                   ! It does not: simply nullify the next mergee pointer of the mergee.
                   nodeList(nodes(iNode)%isolatedNodeIndex)%node%siblingMergee => null()
                end if
                ! Append the mergee as the first mergee on the target node.
                nodeMergeTarget%firstMergee => nodeList(nodes(iNode)%isolatedNodeIndex)%node
             end if
          end if
       end do
    end if
    return
  end subroutine readAssignMergers

  subroutine readScanForBranchJumps(self,nodes,nodeList)
    !!{
    Search for subhalos which move between branches/trees.
    !!}
    use :: Display                   , only : displayMessage, verbosityLevelWarn
    use :: Galacticus_Nodes          , only : treeNodeList
    use :: ISO_Varying_String        , only : assignment(=) , operator(//)
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: String_Handling           , only : operator(//)
    implicit none
    class           (mergerTreeConstructorRead)              , intent(inout)          :: self
    class           (nodeData                 ), dimension(:), intent(inout), target  :: nodes
    type            (treeNodeList             ), dimension(:), intent(inout)          :: nodeList
    class           (nodeData                 )                             , pointer :: currentHost  , descendantNode  , hostDescendant      , jumpToHost, &
         &                                                                               previousNode , isolatedHostNode, isolatedHostHostNode
    integer                                                                           :: iNode
    integer         (c_size_t                 )                                       :: iIsolatedNode
    logical                                                                           :: isMergerEvent, subhaloJumps    , wasMergerEvent
    double precision                                                                  :: timeOfJump
    type            (varying_string           )                                       :: message

    ! If branch jumps are not allowed, simply return.
    if (.not.self%allowBranchJumps) return
    ! Search for subhalos whose descendants live in a different host than that to which their
    ! host descends. These subhalos are jumping between tree branches (or between trees). Add
    ! an event to such nodes to handle the jump.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node (or an initial satellite).
       if (nodes(iNode)%primaryIsolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
          iIsolatedNode=nodes(iNode)%primaryIsolatedNodeIndex
          ! Find those which are a subhalo, or whose descendant is a subhalo.
          descendantNode => null()
          if      (           nodes(iNode)%isSubhalo  ) then
             descendantNode => nodes(iNode)
             previousNode   => nodes(iNode)
          else if (associated(nodes(iNode)%descendant)) then
             if (nodes(iNode)%descendant%isSubhalo) then
                descendantNode => nodes(iNode)%descendant
                previousNode   => nodes(iNode)
             end if
          end if
          ! Check for an immediate subhalo-subhalo merger. If found, nullify the descendant,
          ! so we do not attempt to process this branch.
          if (self%isSubhaloSubhaloMerger(nodes,nodes(iNode))) then
             descendantNode => null()
             currentHost    => readLastHostDescendant(nodes(iNode))
             ! Add a jump if the tree ends before the descendant time.
             if (currentHost%nodeTime <= nodes(iNode)%descendant%nodeTime) then
                timeOfJump     =  currentHost%nodeTime
                jumpToHost     => nodes(iNode)%descendant%host
                do while (jumpToHost%isSubhalo)
                   if (.not.self%warningNestedHierarchyIssued) then
                      message='nested hierarchy detected [node '
                      message=message//nodes(iNode)%descendant%nodeIndex//']'
                      message=message//char(10)//'ignoring as not currently supported'
                      message=message//char(10)//'warning will not be issued again'
                      call displayMessage(message,verbosityLevelWarn)
                      self%warningNestedHierarchyIssued=.true.
                   end if
                   jumpToHost => jumpToHost%host
                end do
                call readCreateBranchJumpEvent(                                                    &
                     &                         nodeList(iIsolatedNode                      )%node, &
                     &                         nodeList(jumpToHost%primaryIsolatedNodeIndex)%node, &
                     &                         timeOfJump                                          &
                     &                        )
             end if
          end if
          ! If a subhalo was found, follow its descent.
          wasMergerEvent=.false.
          do while (associated(descendantNode))
             subhaloJumps=.false.
             timeOfJump  =-1.0d0
             if (descendantNode%isSubhalo.and.associated(descendantNode%descendant)) then
                ! Determine if this is actually a merger event rather than a branch jump.
                ! Assume it is not a merger initially.
                isMergerEvent=.false.
                if (descendantNode%descendant%isSubhalo) then
                   ! Descendant is a subhalo. Check for subhalo-subhalo merger.
                   if (self%isSubhaloSubhaloMerger(nodes,previousNode)) isMergerEvent=.true.
                else
                   ! Descendant is not a subhalo, so this must be a merger event.
                   isMergerEvent=.true.
                end if
                ! If this is a merger event, then check that the current descendant's host has a
                ! descendant that exists beyond the time of the merger. If it does not, then we
                ! still need to allow our node to jump branches (if necessary) as it will not be
                ! able to evolve in the descendantless host.
                wasMergerEvent=isMergerEvent
                if (isMergerEvent) then
                   currentHost => readLastHostDescendant(descendantNode)
                   if (currenthost%nodeTime <= descendantNode%descendant%nodeTime) then
                      isMergerEvent=.false.
                      timeOfJump=currentHost%nodeTime
                   endif
                end if
                ! Proceed only if this is not a merger event.
                if (.not.isMergerEvent) then
                   ! Does this subhalo's descendant live in the host to which the subhalo's host descends.
                   if (.not.associated(descendantNode%host%descendant)) then
                      ! Host has no descendant, so this must be a branch jump.
                      subhaloJumps=.true.
                   else
                      ! In nested hierarchies we must find the isolated node which hosts our node and our node's host.
                      !! First find the ultimate host of our descendant.
                      isolatedHostNode     => descendantNode     %descendant%host
                      do while (associated(isolatedHostNode    %host).and..not.associated(isolatedHostNode    %host,isolatedHostNode    ))
                         isolatedHostNode     => isolatedHostNode    %host
                      end do
                      !! Next find our ultimate host...
                      isolatedHostHostNode => descendantNode%host
                      do while (associated(isolatedHostHostNode%host).and..not.associated(isolatedHostHostNode%host,isolatedHostHostNode))
                         isolatedHostHostNode => isolatedHostHostNode%host
                      end do
                      !! and then find the ultimate host of its descendant.
                      isolatedHostHostNode => isolatedHostHostNode%descendant%host
                      do while (associated(isolatedHostHostNode%host).and..not.associated(isolatedHostHostNode%host,isolatedHostHostNode))
                         isolatedHostHostNode => isolatedHostHostNode%host
                      end do
                      if (isolatedHostNode%nodeIndex /= isolatedHostHostNode%nodeIndex) then
                         ! Host has a descendant, but its host is not the same as our descendant's host.
                         subhaloJumps=.true.
                         ! Check that is not simply a case of the subhalo skipping one or more timesteps before
                         ! reappearing in the expected host.
                         hostDescendant => isolatedHostHostNode
                         if (isolatedHostNode%nodeTime > hostDescendant%nodeTime) then
                            ! Handle cases where the subhalo skipped one or more timesteps.
                            do while (isolatedHostNode%nodeTime > hostDescendant%nodeTime)
                               if (associated(hostDescendant%descendant)) then
                                  hostDescendant => hostDescendant%descendant%host
                                  ! In nested hierarchies, find the isolated host node.
                                  do while (associated(hostDescendant%host).and..not.associated(hostDescendant%host,hostDescendant))
                                     hostDescendant => hostDescendant%host
                                  end do
                               else
                                  exit
                               end if
                            end do
                         else if (isolatedHostNode%nodeTime < hostDescendant%nodeTime) then
                            ! Handle cases where the host skipped one or more timesteps.
                            do while (isolatedHostNode%nodeTime < hostDescendant%nodeTime)
                               if (associated(isolatedHostNode%descendant)) then
                                  isolatedHostNode => isolatedHostNode%descendant%host
                                  ! In nested hierarchies, find the isolated host node.
                                  do while (associated(isolatedHostNode%host).and..not.associated(isolatedHostNode%host,isolatedHostNode))
                                     isolatedHostNode => isolatedHostNode%host
                                  end do
                               else
                                  exit
                               end if
                            end do
                         end if
                         ! Subhalo reappeared in the expected host. This is not a branch jump.
                         if (isolatedHostNode%nodeIndex == hostDescendant%nodeIndex) subhaloJumps=.false.
                      end if
                   end if
                else
                   ! Since this is a merger event, we're finished checking this branch.
                   exit
                end if
             end if
             ! If a jump was detected, create an event.
             if (subhaloJumps) then
                if (timeOfJump < 0.0d0)                   &
                     & timeOfJump=descendantNode%nodeTime
                 jumpToHost => descendantNode%descendant%host
                ! Find an isolated host.
                do while (jumpToHost%isSubhalo)
                   jumpToHost => jumpToHost%host
                end do
                call readCreateBranchJumpEvent(                                                    &
                     &                         nodeList(iIsolatedNode                      )%node, &
                     &                         nodeList(jumpToHost%primaryIsolatedNodeIndex)%node, &
                     &                         timeOfJump                                          &
                     &                        )
             end if
             ! Move to the descendant.
             previousNode   => descendantNode
             descendantNode => descendantNode%descendant
             ! If the descendant is not a subhalo, then we're finished checking this branch.
             if (associated(descendantNode)) then
                if (.not.descendantNode%isSubhalo) exit
             end if
             ! If this was a merger event, then we're finished checking this branch.
             if (wasMergerEvent) exit
          end do
       end if
    end do
    return
  end subroutine readScanForBranchJumps

  function readLastHostDescendant(node) result (currentHost)
    !!{
    Return a pointer to the last descendant that can be reached from {\normalfont \ttfamily node} when descending through hosts.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class(nodeData), pointer               :: currentHost
    class(nodeData), intent(inout), target :: node

    currentHost => node%host
    do while(associated(currentHost%descendant))
       currentHost => currentHost%descendant%host
    end do
    return
  end function readLastHostDescendant

  subroutine readCreateBranchJumpEvent(node,jumpToHost,timeOfJump)
    !!{
    Create a matched-pair of branch jump events in the given nodes.
    !!}
    use :: Galacticus_Nodes , only : nodeEvent       , nodeEventBranchJump, treeNode
    use :: Node_Branch_Jumps, only : Node_Branch_Jump
    implicit none
    type            (treeNode ), intent(inout), target  :: jumpToHost, node
    double precision           , intent(in   )          :: timeOfJump
    class           (nodeEvent)               , pointer :: newEvent  , pairEvent

    allocate(nodeEventBranchJump ::  newEvent)
    allocate(nodeEventBranchJump :: pairEvent)
    call node      %attachEvent( newEvent)
    call jumpToHost%attachEvent(pairEvent)
    newEvent %time =  timeOfJump
    newEvent %node => jumpToHost
    newEvent %task => Node_Branch_Jump
    pairEvent%time =  timeOfJump
    pairEvent%node => node
    pairEvent%task => null()
    pairEvent%ID   =  newEvent%ID
    return
  end subroutine readCreateBranchJumpEvent

  subroutine readBuildSubhaloMassHistories(self,nodes,nodeList,historyCountMaximum,historyTime,historyIndex,historyMass,position,velocity)
    !!{
    Build and attached bound mass histories to subhalos.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : defaultSatelliteComponent, nodeComponentPosition, nodeComponentSatellite, treeNodeList
    use :: Histories                 , only : history                  , longIntegerHistory
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: String_Handling           , only : operator(//)
    implicit none
    class           (mergerTreeConstructorRead)                         , intent(inout) :: self
    class           (nodeData                 ), target , dimension(:  ), intent(inout) :: nodes
    type            (treeNodeList             )         , dimension(:  ), intent(inout) :: nodeList
    integer         (kind_int8                )         , dimension(:  ), intent(inout) :: historyIndex
    integer         (c_size_t                 )                         , intent(in   ) :: historyCountMaximum
    double precision                                    , dimension(:  ), intent(inout) :: historyMass        , historyTime
    double precision                                    , dimension(:,:), intent(inout) :: position           , velocity
    class           (nodeData                 ), pointer                                :: progenitorNode     , node
    class           (nodeComponentSatellite   ), pointer                                :: satellite
    class           (nodeComponentPosition    ), pointer                                :: position_
    integer         (c_size_t                 )                                         :: historyCount       , iIsolatedNode
    integer                                                                             :: iNode
    logical                                                                             :: endOfBranch
    type            (varying_string           )                                         :: message
    type            (history                  )                                         :: subhaloHistory
    type            (longIntegerHistory       )                                         :: subhaloIndexHistory
    type            (progenitorIterator       )                                         :: progenitors

    if (self%presetSubhaloMasses.or.self%presetPositions.or.self%presetSubhaloIndices) then
       ! Check that preset subhalo masses are supported.
       if (self%presetSubhaloMasses .and..not.defaultSatelliteComponent%boundMassHistoryIsSettable()) &
            & call Error_Report('presetting subhalo masses requires a component that supports setting of node bound mass histories'//{introspection:location})
       ! Check that preset subhalo masses are supported.
       if (self%presetSubhaloIndices.and..not.defaultSatelliteComponent%nodeIndexHistoryIsSettable()) &
            & call Error_Report('presetting subhalo indices requires a component that supports setting of node index histories'//{introspection:location})
       ! Get the default cosmology functions object.
       historyBuildNodeLoop: do iNode=1,size(nodes)
          historyBuildIsolatedSelect: if (nodes(iNode)%primaryIsolatedNodeIndex /= nodeReachabilityUnreachable%ID) then
             iIsolatedNode=nodes(iNode)%primaryIsolatedNodeIndex
             ! Set a pointer to the current node - this will be updated if any descendants are traced.
             node => nodes(iNode)
             ! Set initial number of times in the history to zero.
             historyCount=0
             ! Find the subset with descendants.
             historyBuildHasDescendantSelect: if (associated(nodes(iNode)%descendant)) then
                ! Select the subset which have a subhalo as a descendant and are not the primary progenitor or are initial subhalos. Also skip immediate subhalo-subhalo mergers.
                historyBuildSubhaloSelect: if ((nodes(iNode)%descendant%isSubhalo.or.nodes(iNode)%isSubhalo).and..not.self%isSubhaloSubhaloMerger(nodes,nodes(iNode))) then
                   ! Trace descendants until merging or final time.
                   if (nodes(iNode)%isSubhalo) then
                      node => nodes(iNode)
                   else
                      node => nodes(iNode)%descendant
                   end if
                   endOfBranch =.false.
                   historyBuildBranchWalk: do while (.not.endOfBranch)
                      ! Increment the history count for this branch.
                      historyCount=historyCount+1
                      ! Check for history count array size being exceeded.
                      if (historyCount > historyCountMaximum) then
                         message='history array length exceeded for node ['
                         message=message//nodes(iNode)%nodeIndex//'] - this should not happen'
                         call Error_Report(message//{introspection:location})
                      end if
                      ! Store the history.
                      historyTime(historyCount)=node%nodeTime
                      if (self%presetSubhaloIndices) historyIndex(  historyCount)=node%nodeIndex
                      if (self%presetSubhaloMasses ) historyMass (  historyCount)=node%nodeMass
                      if (self%presetPositions     ) position    (:,historyCount)=node%position
                      if (self%presetPositions     ) velocity    (:,historyCount)=node%velocity
                      ! Test the branch.
                      if (.not.associated(node%descendant).or..not.node%descendant%isSubhalo) then
                         ! End of branch reached.                         
                         endOfBranch=.true.
                      else
                         ! Check if merges with another subhalo.
                         call progenitors%descendantSet(self,node%descendant,nodes)
                         do while (progenitors%next(nodes))
                            progenitorNode => progenitors%current(nodes)
                            if     (                                                                                        &
                                 &                      progenitorNode%nodeIndex         /= node%nodeIndex                  &
                                 &  .and.               progenitorNode%isolatedNodeIndex /= nodeReachabilityUnreachable%ID  &
                                 &  .and.    associated(progenitorNode%descendant                                         ) &
                                 &  .and. massIsGreater(progenitorNode                   ,  node                          ) &
                                 & ) then
                               ! Subhalo-subhalo merger.                               
                               endOfBranch =.true.
                               exit
                            end if
                         end do
                         ! Step to the next descendant.
                         if (.not.endOfBranch) node => node%descendant
                      end if
                   end do historyBuildBranchWalk
                   ! Set the mass history for this node.
                   if (self%presetSubhaloMasses) then
                      call subhaloHistory%destroy()
                      call subhaloHistory%create(1,int(historyCount))
                      subhaloHistory%time(:  )=historyTime(1:historyCount)
                      subhaloHistory%data(:,1)=historyMass(1:historyCount)
                      satellite => nodeList(iIsolatedNode)%node%satellite()
                      call satellite%boundMassHistorySet(subhaloHistory)
                   end if
                   ! Set the node index history for this node.
                   if (self%presetSubhaloIndices) then
                      call subhaloIndexHistory%destroy()
                      call subhaloIndexHistory%create(1,int(historyCount))
                      subhaloIndexHistory%time(:  )=historyTime (1:historyCount)
                      subhaloIndexHistory%data(:,1)=historyIndex(1:historyCount)
                      satellite       => nodeList(iIsolatedNode)%node%satellite()
                      call satellite%nodeIndexHistorySet(subhaloIndexHistory)
                   end if
                end if historyBuildSubhaloSelect
             end if historyBuildHasDescendantSelect
             ! Set the position history for this node.
             if (self%presetPositions.and..not.nodeList(iIsolatedNode)%node%isPrimaryProgenitor()) then
                ! Check if particle data is available for this node.
                if (self%mergerTreeImporter_%subhaloTraceCount(node) > 0) then
                   ! Check that arrays are large enough to hold particle data. They should be. If they are not, it's a
                   ! bug.
                   if (historyCount+self%mergerTreeImporter_%subhaloTraceCount(node) > size(historyTime)) then
                      message='history arrays are too small to hold data for node '
                      message=message//nodeList(iIsolatedNode)%node%index()//': ['//historyCount//'+'//self%mergerTreeImporter_%subhaloTraceCount(node)//']='//(historyCount+self%mergerTreeImporter_%subhaloTraceCount(node))//'>'//size(historyTime)
                      call Error_Report(message//{introspection:location})
                   end if
                   ! Read subhalo position trace data.
                   call self%mergerTreeImporter_%subhaloTrace                                                           &
                        & (                                                                                             &
                        &  node                                                                                       , &
                        &  historyTime(  historyCount+1:historyCount+self%mergerTreeImporter_%subhaloTraceCount(node)), &
                        &  position   (:,historyCount+1:historyCount+self%mergerTreeImporter_%subhaloTraceCount(node)), &
                        &  velocity   (:,historyCount+1:historyCount+self%mergerTreeImporter_%subhaloTraceCount(node))  &
                        & )
                   ! Increment the history count for this node.
                   historyCount=historyCount+self%mergerTreeImporter_%subhaloTraceCount(node)
                end if                
                if (historyCount > 0) then
                   call subhaloHistory%destroy()
                   call subhaloHistory%create(6,int(historyCount))
                   subhaloHistory%time(:    )=          historyTime(    1:historyCount)
                   subhaloHistory%data(:,1:3)=transpose(position   (1:3,1:historyCount))
                   subhaloHistory%data(:,4:6)=transpose(velocity   (1:3,1:historyCount))
                   position_ => nodeList(iIsolatedNode)%node%position()
                   call position_%positionHistorySet(subhaloHistory)
                end if
             end if
          end if historyBuildIsolatedSelect
       end do historyBuildNodeLoop
    end if
    return
  end subroutine readBuildSubhaloMassHistories

  subroutine readValidateIsolatedHalos(nodes)
    !!{
    Ensure that nodes have valid primary progenitors.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class  (nodeData          ), dimension(:), intent(inout) :: nodes
    type   (treeNode          ), pointer                     :: nodeNew , nodeSatellite
    class  (nodeComponentBasic), pointer                     :: basicNew
    integer                                                  :: iNode   , i

    ! Search for cases where a node has no progenitors which do not descend into subhalos.
    do iNode=1,size(nodes)
       ! Only process if this is an isolated node.
       if (nodes(iNode)%isolatedNodeIndex /= nodeReachabilityUnreachable%ID .and. .not.nodes(iNode)%isSubhalo) then
          ! Select nodes with parents.
          if (associated(nodes(iNode)%node%parent)) then
             ! Select nodes with subhalo descendants which are also the primary progenitor of their parent.
             if (nodes(iNode)%descendant%isSubhalo.and.associated(nodes(iNode)%node%parent%firstChild,nodes(iNode)%node)) then
                ! Insert a copy of the parent node as its own primary progenitor. This avoids current node being promoted into its
                ! parent even though it is intended to descend into a subhalo. The copy is shifted to a very slightly earlier
                ! time to avoid having two identical halos existing simultaneously (which can be problematic if outputting
                ! quantities which use the node index as a label in dataset names for example).
                allocate(nodeNew)
                call nodes(iNode)%node%parent%copyNodeTo(nodeNew)
                if (nodeNew%satelliteCount() > 0) then
                   ! Remove any satellite component from the copied node - each branch should have only a single satellite.
                   do i=nodeNew%satelliteCount(),1,-1
                      call nodeNew%satelliteRemove(i)
                   end do
                end if
                nodeNew%sibling                     => nodes(iNode)%node
                nodeNew%parent                      => nodes(iNode)%node%parent
                nodeNew%firstChild                  => null()
                nodeNew%mergeTarget                 => null()
                nodeNew%siblingMergee               => null()
                nodes(iNode)%node%parent%firstChild => nodeNew
                basicNew                            => nodeNew%basic()
                call basicNew%timeSet(basicNew%time()*(1.0d0-fractionOffsetTimeClones))
                ! Events remain attached to the original and we do not want to duplicate them.
                nodeNew%event => null()
                ! Any satellites are now attached to the copy.
                nodes(iNode)%node%parent%firstSatellite => null()
                nodeSatellite => nodeNew%firstSatellite
                do while (associated(nodeSatellite))
                   nodeSatellite%parent => nodeNew
                   nodeSatellite        => nodeSatellite%sibling
                end do
             end if
          end if
       end if
    end do
    return
  end subroutine readValidateIsolatedHalos

  subroutine readAssignUniqueIDsToClones(nodeList)
    !!{
    Assign new uniqueID values to any cloned nodes inserted into the trees.
    !!}
    use :: Galacticus_Nodes, only : treeNodeList
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
  end subroutine readAssignUniqueIDsToClones

  logical function readIsSubhaloSubhaloMerger(self,nodes,node)
    !!{
    Returns true if {\normalfont \ttfamily node} undergoes a subhalo-subhalo merger.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class(mergerTreeConstructorRead)              , intent(inout) :: self
    class(nodeData                 ), dimension(:), intent(inout) :: nodes
    class(nodeData                 )              , intent(in   ) :: node
    class(nodeData                 ), pointer                     :: progenitorNode
    type (progenitorIterator       )                              :: progenitors

    readIsSubhaloSubhaloMerger=.false.
    ! Return immediately if there is no descendant. (Since there can be no merger if there is no descendant.)
    if (.not.associated(node%descendant          )) return
    ! Return immediately if descendant is not a subhalo, as this could then not be a subhalo-subhalo merger.
    if (.not.           node%descendant%isSubhalo ) return
    ! Check if node's descendant has any progenitor nodes.
    call progenitors%descendantSet(self,node%descendant,nodes)
    do while (progenitors%next(nodes))
       progenitorNode => progenitors%current(nodes)
       if     (                                                                                       &
            &                     progenitorNode%nodeIndex         /= node%nodeIndex                  &
            &  .and.              progenitorNode%isolatedNodeIndex /= nodeReachabilityUnreachable%ID  &
            &  .and.   associated(progenitorNode%descendant                                         ) &
            &  .and.massIsGreater(progenitorNode                   ,  node                          ) &
            & ) then
          ! It does, so this is a subhalo-subhalo merger.
          readIsSubhaloSubhaloMerger =.true.
          exit
       end if
    end do
    return
  end function readIsSubhaloSubhaloMerger

  !$GLC function attributes unused :: readDumpTree
  subroutine readDumpTree(nodes,highlightNodes,branchRoot)
    !!{
    Dumps the tree structure to a file in a format suitable for processing with \href{http://www.graphviz.org/}{\normalfont \scshape dot}.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class    (nodeData      ), dimension(:), intent(in   ), target   :: nodes
    integer  (kind=kind_int8), dimension(:), intent(in   ), optional :: highlightNodes
    integer  (kind=kind_int8)              , intent(in   ), optional :: branchRoot
    class    (nodeData      ), pointer                               :: node
    integer                                                          :: fileUnit      , iNode
    character(len=20        )                                        :: color         , style
    logical                                                          :: outputNode
    integer  (kind=kind_int8)                                        :: branchRootHost

    ! Open an output file and write the GraphViz opening.
    open(newunit=fileUnit,file='mergerTreeConstructReadTree.gv',status='unknown',form='formatted')
    write (fileUnit,*) 'digraph Tree {'

    ! Identify the host of the branch root.
    if (present(branchRoot)) then
       branchRootHost =  branchRoot
       node       => null()
       do iNode=1,size(nodes)
          if (nodes(iNode)%nodeIndex == branchRootHost) node => nodes(iNode)
       end do
       if (associated(node)) then
          do while (associated(node%host).and..not.associated(node%host,node))
             node => node%host
          end do
          branchRootHost=node%nodeIndex
       end if
    else
       branchRootHost=-1_kind_int8
    end if

    ! Loop over all nodes.
    do iNode=1,size(nodes)
       ! Determine if node is in the branch to be output.
       if (present(branchRoot)) then
          outputNode=.false.
          node => nodes(iNode)
          do while (associated(node%host).and..not.associated(node%host,node))
             node => node%host
          end do
          outputNode=node%nodeIndex == branchRootHost
          do while (.not.outputNode.and.associated(node%descendant))
             node => node%descendant
             do while (associated(node%host).and..not.associated(node%host,node))
                node => node%host
             end do
             outputNode=node%nodeIndex == branchRootHost
          end do
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
             write (fileUnit,'(a,i20.20,a,i20.20,a,f5.2,a,a,a,a,a,f5.2,a)') '"',nodes(iNode)%nodeIndex,'" [shape=box   , label="',nodes(iNode)%nodeIndex,':',nodes(iNode)%nodeTime,'", color=',trim(color),', style=',trim(style),', z=',nodes(iNode)%nodeTime,'];'
             ! If a host node is given, add a link to it as a red line.
             if (associated(nodes(iNode)%host)) write (fileUnit,'(a,i20.20,a,i20.20,a)') '"',nodes(iNode)%nodeIndex,'" -> "',nodes(iNode)%host%nodeIndex,'" [color=red];'
          else
             write (fileUnit,'(a,i20.20,a,i20.20,a,f5.2,a,a,a,a,a,f5.2,a)') '"',nodes(iNode)%nodeIndex,'" [shape=circle, label="',nodes(iNode)%nodeIndex,':',nodes(iNode)%nodeTime,'", color=',trim(color),', style=',trim(style),', z=',nodes(iNode)%nodeTime,'];'
          endif
          ! Make a link to the descendant node using a black line.
          if (associated(nodes(iNode)%descendant)) write (fileUnit,'(a,i20.20,a,i20.20,a)') '"',nodes(iNode)%nodeIndex,'" -> "',nodes(iNode)%descendant%nodeIndex,'" ;'
       end if
    end do

    ! Close the file.
    write (fileUnit,*) '}'
    close(fileUnit)
    return
  end subroutine readDumpTree

  subroutine readTimeUntilMergingSubresolution(self,lastSeenNode,nodes,nodeList,iNode,timeSubhaloMerges)
    !!{
    Compute the additional time until merging after a subhalo is lost from the tree (presumably due to limited resolution).
    !!}
    use :: Display                   , only : displayIndent          , displayMessage       , displayUnindent, verbosityLevelWarn
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : nodeComponentBasic     , nodeComponentPosition, treeNode       , treeNodeList
    use :: Kepler_Orbits             , only : keplerOrbit
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: String_Handling           , only : operator(//)
    use :: Vectors                   , only : Vector_Magnitude
    use :: Calculations_Resets       , only : Calculations_Reset
    implicit none
    class           (mergerTreeConstructorRead)                       , intent(inout) :: self
    class           (nodeData                 )                       , intent(in   ) :: lastSeenNode
    class           (nodeData                 ), target , dimension(:), intent(inout) :: nodes
    type            (treeNodeList             )         , dimension(:), intent(inout) :: nodeList
    integer                                                           , intent(in   ) :: iNode
    double precision                                                  , intent(inout) :: timeSubhaloMerges
    class           (nodeData                 ), pointer                              :: primaryProgenitor    , progenitorNode  , &
         &                                                                               node
    class           (nodeComponentBasic       ), pointer                              :: basic
    class           (nodeComponentPosition    ), pointer                              :: position
    type            (treeNode                 ), pointer                              :: hostNode             , satelliteNode
    double precision                                    , dimension(3)                :: relativePosition     , relativeVelocity
    type            (keplerOrbit              )                                       :: orbit
    double precision                                                                  :: primaryProgenitorMass, timeUntilMerging
    type            (progenitorIterator       )                                       :: progenitors
    character       (len=42                   )                                       :: coordinateLabel
    type            (varying_string           )                                       :: message
    logical                                                                           :: parentIsCloned

    ! Find the nodes that descend into our target node's descendant.
    call progenitors%descendantSet(self,lastSeenNode%descendant,nodes)
    if (progenitors%exist()) then
       ! Determine if the parent node has a clone primary progenitor.
       !! First check if the parent has a child halo.
       if (associated(nodeList(lastSeenNode%isolatedNodeIndex)%node%parent%firstChild)) then
          ! Parent does have a child halo - check if it is a clone.
          parentIsCloned=(nodeList(lastSeenNode%isolatedNodeIndex)%node%parent%uniqueID() == nodeList(lastSeenNode%isolatedNodeIndex)%node%parent%firstChild%uniqueID())
       else
          ! Parent does not have a child halo - this means our current halo must be a subhalo in a childless host. Therefore, the
          ! parent can not have a clone.
          parentIsCloned=.false.
       end if
       ! If parent is cloned, we need to make a temporary progenitor node.
       if (parentIsCloned) then
          allocate(primaryProgenitor)
          basic                          => nodeList(lastSeenNode%isolatedNodeIndex)%node%parent%firstChild%basic   ()
          position                       => nodeList(lastSeenNode%isolatedNodeIndex)%node%parent%firstChild%position()
          primaryProgenitor%nodeIndex    =  nodeList(lastSeenNode%isolatedNodeIndex)%node%parent%firstChild%index   ()
          primaryProgenitor%nodeMass     =  basic                                                          %mass    ()
          if (self%presetPositions) then
             primaryProgenitor%position  =  position                                                       %position()
             primaryProgenitor%velocity  =  position                                                       %velocity()
          end if
       else
          primaryProgenitor     => null()
          primaryProgenitorMass =  0.0d0
          do while (progenitors%next(nodes))
             progenitorNode => progenitors%current(nodes)
             if (progenitorNode%nodeIndex /= lastSeenNode%nodeIndex .and. progenitorNode%nodeMass > primaryProgenitorMass) then
                primaryProgenitorMass =  progenitorNode%nodeMass
                primaryProgenitor     => progenitorNode
             end if
          end do
       end if
       ! Initialize time until merging to zero.
       timeUntilMerging=0.0d0
       ! If position information is available, compute the subresolution orbit.
       if (self%presetPositions) then
          ! Find relative position and velocity.
          relativePosition=lastSeenNode%position-primaryProgenitor%position
          relativeVelocity=lastSeenNode%velocity-primaryProgenitor%velocity
          ! Update position/velocity for periodicity and Hubble flow.
          call self%phaseSpacePositionRealize(lastSeenNode%nodeTime,relativePosition,relativeVelocity)
          ! Catch zero separation halos.
          if (Vector_Magnitude(relativePosition) == 0.0d0) then
             message='merging halos ['
             message=message//lastSeenNode%nodeIndex//' & '//primaryProgenitor%nodeIndex//'] have zero separation'
             call displayIndent  (message                          ,verbosityLevelWarn)
             write (coordinateLabel,'("[",e12.6,",",e12.6,",",e12.6,"]")') primaryProgenitor%position
             message="position [primary  ] = "//trim(coordinateLabel)
             call displayMessage (message                          ,verbosityLevelWarn)
             write (coordinateLabel,'("[",e12.6,",",e12.6,",",e12.6,"]")') lastSeenNode     %position
             message="position [satellite] = "//trim(coordinateLabel)
             call displayMessage (message                          ,verbosityLevelWarn)
             write (coordinateLabel,'("[",e12.6,",",e12.6,",",e12.6,"]")') primaryProgenitor%velocity
             message="velocity [primary  ] = "//trim(coordinateLabel)
             call displayMessage (message                          ,verbosityLevelWarn)
             write (coordinateLabel,'("[",e12.6,",",e12.6,",",e12.6,"]")') lastSeenNode     %velocity
             message="velocity [satellite] = "//trim(coordinateLabel)
             call displayMessage (message                          ,verbosityLevelWarn)
             call displayUnindent('assuming instantaneous merging' ,verbosityLevelWarn)
          else
             ! Create the orbit.
             orbit=readOrbitConstruct(lastSeenNode%nodeMass,primaryProgenitor%nodeMass,relativePosition,relativeVelocity)
             ! Construct temporary nodes.
             satelliteNode                => treeNode()
             hostNode                     => treeNode()
             call    nodeList(lastSeenNode     %isolatedNodeIndex)%node                  %copyNodeTo(satelliteNode,skipEvent=.true.)
             if (parentIsCloned) then
                call nodeList(lastSeenNode     %isolatedNodeIndex)%node%parent%firstChild%copyNodeTo(hostNode     ,skipEvent=.true.)
             else
                call nodeList(primaryProgenitor%isolatedNodeIndex)%node                  %copyNodeTo(hostNode     ,skipEvent=.true.)
             end if
             satelliteNode%parent         => hostNode
             hostNode     %firstSatellite => satelliteNode
             ! Perform a calculation reset as technically these nodes have changed. (Specifically, they may have the same unique
             ! ID as the prior time this function was called, and yet be new copies. Not resetting calculations could result in
             ! the old - now destroyed - copies of these nodes being accessed.)
             call Calculations_Reset(satelliteNode)
             ! Determine the time until merging.
             timeUntilMerging=self%satelliteMergingTimescales_%timeUntilMerging(satelliteNode,orbit)
             ! Clean up.
             call satelliteNode%destroy()
             call hostNode     %destroy()
             deallocate(satelliteNode)
             deallocate(hostNode     )
          end if
       end if
       ! Find the new merging time, and the node with which the merging will occur.
       node              => lastSeenNode%descendant
       timeSubhaloMerges =  timeSubhaloMerges+timeUntilMerging
       do while (associated(node%descendant))
          if (node%descendant%nodeTime > timeSubhaloMerges) then
             nodes(iNode)%mergesWithIndex=node%nodeIndex
             exit
          else
             node => node%descendant
          end if
       end do
       ! Merging time is beyond the end of the tree. Set merging time to infinity.
       if (.not.associated(node%descendant)) nodes(iNode)%mergesWithIndex=node%nodeIndex
       ! Clean up any temporary progenitor.
       if (parentIsCloned) deallocate(primaryProgenitor)
    else
       call Error_Report('no descendants found'//{introspection:location})
    end if
    return
  end subroutine readTimeUntilMergingSubresolution

  function readOrbitConstruct(mass1,mass2,position,velocity) result(orbit)
    !!{
    Construct a Keplerian orbit given body masses, positions, and relative velocities.
    !!}
    use :: Kepler_Orbits           , only : keplerOrbit
    use :: Numerical_Constants_Math, only : Pi
    use :: Vectors                 , only : Vector_Magnitude, Vector_Product
    implicit none
    type            (keplerOrbit)                              :: orbit
    double precision                           , intent(in   ) :: mass1             , mass2
    double precision             , dimension(3), intent(in   ) :: position          , velocity
    double precision             , dimension(3)                :: velocityTangential, velocityRadial   , &
         &                                                        vectorEpsilon1    , vectorEpsilon2   , &
         &                                                        vectorRadial
    double precision                                           :: positionMagnitude , velocityEpsilon1 , &
         &                                                        velocityEpsilon2  , epsilon1Magnitude, &
         &                                                        epsilon

    positionMagnitude  =Vector_Magnitude(position)
    vectorRadial       =+position          &
         &              /positionMagnitude
    velocityRadial     =+Dot_Product(velocity,position) &
         &              *vectorRadial
    velocityTangential =+velocity       &
         &              -velocityRadial
    vectorEpsilon1     =Vector_Product(vectorRadial  ,[0.0d0,0.0d0,1.0d0])
    epsilon1Magnitude=Vector_Magnitude(vectorEpsilon1)
    if (epsilon1Magnitude > 0.0d0) then
       vectorEpsilon2  =Vector_Product(vectorEpsilon1,vectorRadial       )
       vectorEpsilon1  =vectorEpsilon1/Vector_Magnitude(vectorEpsilon1)
       vectorEpsilon2  =vectorEpsilon2/Vector_Magnitude(vectorEpsilon2)
       velocityEpsilon1=Dot_Product(velocityTangential,vectorEpsilon1  )
       velocityEpsilon2=Dot_Product(velocityTangential,vectorEpsilon2  )
       epsilon         =atan2      (velocityEpsilon2  ,velocityEpsilon1)
    else
       epsilon         =Pi/2.0d0
    end if
    call orbit%reset()
    call orbit%massesSet            (       &
         &                           mass1, &
         &                           mass2  &
         &                          )
    call orbit%radiusSet            (                                     positionMagnitude)
    call orbit%velocityRadialSet    (Dot_Product     (velocity,position )/positionMagnitude)
    call orbit%velocityTangentialSet(Vector_Magnitude(velocityTangential)                  )
    call orbit%thetaSet             (acos (position        (3)/positionMagnitude                    ))
    call orbit%phiSet               (atan2(position        (2)                  ,position        (1)))
    call orbit%epsilonSet           (epsilon                                                         )
    return
  end function readOrbitConstruct

  subroutine readPhaseSpacePositionRealize(self,time,position,velocity)
    !!{
    Modify relative positions and velocities to account for both any periodicity of the simulated volume, and for Hubble flow.
    !!}
    use :: Numerical_Constants_Boolean, only : booleanFalse, booleanTrue
    implicit none
    class           (mergerTreeConstructorRead)              , intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    double precision                           , dimension(3), intent(inout) :: position           , velocity
    double precision                                                         :: lengthSimulationBox

    ! Account for periodicity.
    if (self%mergerTreeImporter_%positionsArePeriodic() /= booleanFalse) then
       lengthSimulationBox=self%mergerTreeImporter_%cubeLength(time)
       position=mod(position+0.5d0*lengthSimulationBox,lengthSimulationBox)-0.5d0*lengthSimulationBox
       position=mod(position-0.5d0*lengthSimulationBox,lengthSimulationBox)+0.5d0*lengthSimulationBox
    end if
    ! Account for Hubble flow.
    if (self%mergerTreeImporter_%velocitiesIncludeHubbleFlow() /= booleanTrue) then
       velocity=velocity                                                    &
            &  +position                                                    &
            &  *self%cosmologyFunctions_%hubbleParameterEpochal(time=time)
    end if
    return
  end subroutine readPhaseSpacePositionRealize

  subroutine progenitorIteratorDescendantSet(self,constructor,node,nodes)
    !!{
    Initialize a progenitor iterator object by storing the index of the target {\normalfont \ttfamily node} and finding the location of the first
    progenitor (if any).
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class(progenitorIterator       )              , intent(inout) :: self
    type (mergerTreeConstructorRead), target      , intent(in   ) :: constructor
    class(nodeData                 )              , intent(in   ) :: node
    class(nodeData                 ), dimension(:), intent(in   ) :: nodes

    ! Store a pointer to the tree constructor.
    self%constructor => constructor
    ! Store the index of the target node.
    self%targetIndex     =node%nodeIndex
    ! Assume no progenitors descendants by default.
    self%progenitorsFound=.false.
    ! Find the index of matching nodes in the list sorted by descendant index.
    self%progenitorIndex=self%constructor%descendantNodeSortIndex(node%nodeIndex)
    if (self%progenitorIndex > 0 .and. self%progenitorIndex <= size(nodes)) then
       ! Progenitors may exist, store the location of the first progenitor if found.
       self%progenitorLocation=self%constructor%descendantLocations(self%progenitorIndex)
       if (associated(nodes(self%progenitorLocation)%descendant)) &
            & self%progenitorsFound=(nodes(self%progenitorLocation)%descendant%nodeIndex == node%nodeIndex)
       ! Increment the initial index so that the first call to get the next progenitor can find the first progenitor by
       ! subtracting one from this index.
       self%progenitorIndex=self%progenitorIndex+1
    end if
    return
  end subroutine progenitorIteratorDescendantSet

  logical function progenitorIteratorNext(self,nodes)
    !!{
    Move to the next progenitor using a progenitor iterator object, returning true if the next progenitor exists, false if it
    does not.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class(progenitorIterator)              , intent(inout) :: self
    class(nodeData          ), dimension(:), intent(in   ) :: nodes

    if (self%progenitorsFound) then
       progenitorIteratorNext=.true.
       self%progenitorLocation=-1
       do while (self%progenitorIndex > 0)
          self%progenitorIndex=self%progenitorIndex-1
          if (self%progenitorIndex <= 0) exit
          self%progenitorLocation=self%constructor%descendantLocations(self%progenitorIndex)
          if (associated(nodes(self%progenitorLocation)%descendant)) exit
       end do
       if (.not.associated(nodes(self%progenitorLocation)%descendant)) then
          progenitorIteratorNext=.false.
          return
       end if
       if (      self%progenitorIndex                          == 0               ) progenitorIteratorNext=.false.
       if (nodes(self%progenitorLocation)%descendant%nodeIndex /= self%targetIndex) progenitorIteratorNext=.false.
    else
       progenitorIteratorNext=.false.
    end if
    return
  end function progenitorIteratorNext

  function progenitorIteratorIndex(self,nodes)
    !!{
    Return the node index of the current progenitor in a progenitor iterator object.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    integer(kind=kind_int8    )                              :: progenitorIteratorIndex
    class  (progenitorIterator)              , intent(in   ) :: self
    class  (nodeData          ), dimension(:), intent(in   ) :: nodes

    progenitorIteratorIndex=nodes(self%progenitorLocation)%nodeIndex
    return
  end function progenitorIteratorIndex

  function progenitorIteratorCurrent(self,nodes)
    !!{
    Return a pointer to the current progenitor in a progenitor iterator object.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class(nodeData          ), pointer                             :: progenitorIteratorCurrent
    class(progenitorIterator)              , intent(in   )         :: self
    class(nodeData          ), dimension(:), intent(in   ), target :: nodes

    progenitorIteratorCurrent => nodes(self%progenitorLocation)
    return
  end function progenitorIteratorCurrent

  logical function progenitorIteratorExist(self)
    !!{
    Return true if progenitors exist, false otherwise.
    !!}
    implicit none
    class(progenitorIterator), intent(in   ) :: self

    progenitorIteratorExist=self%progenitorsFound
    return
  end function progenitorIteratorExist

  subroutine readTimingRecord(self,label)
    !!{
    Record timing data.
    !!}
    implicit none
    class    (mergerTreeConstructorRead), intent(inout)               :: self
    character(len=*                    ), intent(in   )               :: label
    real                                , allocatable  , dimension(:) :: timingTimesTmp
    type     (varying_string           ), allocatable  , dimension(:) :: timingLabelsTmp
    real                                                              :: timeNow

    call CPU_Time(timeNow)
    if (allocated(self%timingTimes)) then
       call Move_Alloc(self%timingTimes ,timingTimesTmp )
       call Move_Alloc(self%timingLabels,timingLabelsTmp)
       allocate(self%timingTimes (size(timingTimesTmp )+1))
       allocate(self%timingLabels(size(timingLabelsTmp)+1))
       self%timingTimes (1:size(timingTimesTmp ))=timingTimesTmp
       self%timingLabels(1:size(timingLabelsTmp))=timingLabelsTmp
       deallocate(timingTimesTmp )
       deallocate(timingLabelsTmp)
    else
      allocate(self%timingTimes (1))
      allocate(self%timingLabels(1))
    end if
    self%timingTimes (size(self%timingTimes ))=timeNow
    self%timingLabels(size(self%timingLabels))=trim(label)
    return
  end subroutine readTimingRecord

  subroutine readTimingReport(self)
    !!{
    Report on time taken in various steps of processing merger trees read from file.
    !!}
    use :: Display, only : displayIndent, displayMessage, displayUnindent
    implicit none
    class    (mergerTreeConstructorRead), intent(inout) :: self
    integer                                             :: i        , lengthMaximum
    type     (varying_string           )                :: message
    character(len=12                   )                :: timeTaken

    call displayIndent("Merger tree read processing report:")
    if (allocated(self%timingTimes).and.size(self%timingTimes) > 1) then
       lengthMaximum=maxval(len(self%timingLabels))
       do i=2,size(self%timingTimes)
          write (timeTaken,'(f10.2)') self%timingTimes(i)-self%timingTimes(i-1)
          message=repeat(" ",lengthMaximum-len(self%timingLabels(i)))//self%timingLabels(i)//": "//trim(adjustl(timeTaken))//" s"
          call displayMessage(message)
       end do
    end if
    call displayUnindent("done")
    return
  end subroutine readTimingReport

  subroutine readRootNodeAffinitiesInitial(self,nodes)
    !!{
    Find initial root node affinities for all nodes.
    !!}
    use :: Display                   , only : displayIndent     , displayMessage, displayUnindent, verbosityLevelInfo, &
          &                                   verbosityLevelWarn
    use :: ISO_Varying_String        , only : assignment(=)     , operator(//)  , varying_string
    use :: Merger_Tree_Read_Importers, only : nodeDataMinimal
    use :: Sorting                   , only : sortIndex
    use :: String_Handling           , only : operator(//)
    implicit none
    class  (mergerTreeConstructorRead)                        , intent(inout) :: self
    class  (nodeDataMinimal          ), dimension(         : ), intent(inout) :: nodes
    integer(kind_int8                ), dimension(size(nodes))                :: rootAffinity
    integer(c_size_t                 )                                        :: i                , j                       , &
         &                                                                       treeCount        , pushCount               , &
         &                                                                       k                , progenitorLocation      , &
         &                                                                       forestSizeI      , forestSizeJ
    integer(kind_int8                )                                        :: treeIndexPrevious, treeStartPrevious       , &
         &                                                                       progenitorIndex
    type   (varying_string           )                                        :: message
    logical                                                                   :: nodeIsMostMassive, isolatedProgenitorExists, &
         &                                                                       nodeIsPrimary

    ! Get a unique ID for this split forest.
    self%splitForestUniqueID=splitForestUniqueID%increment()
    ! Build sorted indices into nodes.
    call self%createNodeIndices(nodes)
    ! Initialize root affinities to impossible value.
    rootAffinity=-1_kind_int8
    ! Iterate over nodes.
    do i=1,size(nodes)
       ! Trace through hosts until a self-hosting node is found.
       j=i
       do while (nodes(j)%nodeIndex /= nodes(j)%hostIndex)
          j=self%nodeLocation(nodes(j)%hostIndex)
       end do
       ! Trace descendants until a root node is reached.
       do while (nodes(j)%descendantIndex >= 0)
          ! Jump to descendant.
          if (nodes(j)%descendantIndex >= 0) j=self%nodeLocation(nodes(j)%descendantIndex)
          ! Trace through hosts until a self-hosting node is found.
          do while (nodes(j)%nodeIndex /= nodes(j)%hostIndex)
             j=self%nodeLocation(nodes(j)%hostIndex)
          end do
       end do
       ! Store root affinity.
       rootAffinity(i)=nodes(j)%nodeIndex
    end do
    ! Search for nodes which have a descendant or host with different root affinity and attempt to regroup trees into subforests
    ! of the original forest.
    do i=1,size(nodes)
       ! Process only nodes with a descendant.
       if (nodes(i)%descendantIndex >= 0) then
          ! Check for different root affinity in descendant.
          do k=1,2
             select case (k)
             case (1)
                j=self%nodeLocation(nodes(i)%descendantIndex)
             case (2)
                j=self%nodeLocation(nodes(i)%      hostIndex)
             end select
             if (rootAffinity(i) /= rootAffinity(j)) then
                forestSizeI=count(rootAffinity == rootAffinity(i))
                if (forestSizeI             >= self%forestSizeMaximum) cycle
                forestSizeJ=count(rootAffinity == rootAffinity(j))
                if (forestSizeI+forestSizeJ >  self%forestSizeMaximum) cycle
                where (rootAffinity == rootAffinity(j))
                   rootAffinity=rootAffinity(i)
                end where
              end if
          end do
       end if
    end do
    ! Get a sorted index into the root node affinities.
    allocate(self%splitForestMapIndex(size(nodes)))
    self%splitForestMapIndex=sortIndex(rootaffinity)
    ! Count trees in the forest.
    treeCount        = 0
    treeIndexPrevious=-1
    do i=1,size(nodes)
       if (rootAffinity(self%splitForestMapIndex(i)) /= treeIndexPrevious) then
          treeCount        =treeCount                           +1
          treeIndexPrevious=rootAffinity(self%splitForestMapIndex(i))
       end if
    end do
    ! Identify tree size and start offsets.
    allocate(self%splitForestTreeSize (treeCount))
    allocate(self%splitForestTreeStart(treeCount))
    treeIndexPrevious=-1
    treeStartPrevious= 0
    treeCount        = 0
    do i=1,size(nodes)
       if (rootAffinity(self%splitForestMapIndex(i)) /= treeIndexPrevious) then
          treeCount                      =treeCount                           +1
          treeIndexPrevious              =rootAffinity(self%splitForestMapIndex(i))
          self%splitForestTreeStart(treeCount)=                                 i
       end if
    end do
    if (treeCount > 1) then
       do i=1,treeCount-1
          self%splitForestTreeSize(i)=self%splitForestTreeStart(i+1)-self%splitForestTreeStart(i)
       end do
    end if
    self%splitForestTreeSize(treeCount)=size(nodes)+1-self%splitForestTreeStart(treeCount)
    ! Report.
    call displayIndent('Breaking forest into trees:',verbosityLevelInfo)
    do i=1,treeCount
       message="Tree "
       message=message//i//" of "//treeCount//" contains "//self%splitForestTreeSize(i)//" node"
       if (self%splitForestTreeSize(i) > 1) message=message//"s"
       call displayMessage(message,verbosityLevelInfo)
    end do
    ! Search for nodes which have a descendant with different root affinity.
    pushCount=0
    do i=1,size(nodes)
       ! Process only nodes with a descendant.
       if (nodes(i)%descendantIndex >= 0) then
          ! Check for different root affinity in descendant.
          j=self%nodeLocation(nodes(i)%descendantIndex)
          k=self%nodeLocation(nodes(j)%      hostIndex)
          if     (                                    &
               &   rootAffinity(i) /= rootAffinity(j) &
               &  .or.                                &
               &   rootAffinity(i) /= rootAffinity(k) &
               & ) then
             ! Descendant has different root affinity - this is a cross-tree subhalo promotion event or cross-tree branch jump
             ! event.
             pushCount=pushCount+1
          end if
       end if
    end do
    message="Found "
    message=message//pushCount//" links between trees"
    call displayMessage(message,verbosityLevelInfo)
    ! Build a list of push and pull links.
    allocate(self%splitForestPushTo   (pushCount))
    allocate(self%splitForestPullFrom (pushCount))
    allocate(self%splitForestPushType (pushCount))
    allocate(self%splitForestPushTime (pushCount))
    allocate(self%splitForestIsPrimary(pushCount))
    allocate(self%splitForestPushDone (pushCount))
    allocate(self%splitForestPullDone (pushCount))
    pushCount=0
    do i=1,size(nodes)
       ! Process only nodes with a descendant.
       if (nodes(i)%descendantIndex >= 0) then
          ! Check for different root affinity in descendant.
          j=self%nodeLocation(nodes(i)%descendantIndex)
          ! Test for an inter-tree event. These are identified by a node having a different initial root affinity than its descendant.
          if (rootAffinity(i) /= rootAffinity(j)) then
             ! Determine if our node is the primary progenitor and if an isolated progenitor exists.
             isolatedProgenitorExists=.false.
             nodeIsMostMassive       =.true.
             progenitorIndex         =self%descendantNodeSortIndex(nodes(i)%descendantIndex)
             if (progenitorIndex > 0 .and. progenitorIndex <= size(nodes)) then
                progenitorLocation=self%descendantLocations(progenitorIndex)
                do while (nodes(progenitorLocation)%descendantIndex == nodes(i)%descendantIndex)
                   ! Determine progenitor status.
                   if (nodes(progenitorLocation)%nodeIndex /= nodes(progenitorLocation)%hostIndex) then
                      if     (                                                           &
                           &   nodes(progenitorLocation)%nodeIndex /= nodes(i)%nodeIndex &
                           &  .and.                                                      &
                           &   nodes(progenitorLocation)%nodeMass  >  nodes(i)%nodeMass  &
                           & ) nodeIsMostMassive=.false.
                   else
                      isolatedProgenitorExists=.true.
                   end if
                   ! Move to the next progenitor.
                   progenitorIndex=progenitorIndex-1
                   if (progenitorIndex > 0) then
                      progenitorLocation=self%descendantLocations(progenitorIndex)
                   else
                      exit
                   end if
                end do
             end if
             nodeIsPrimary=nodeIsMostMassive.and..not.isolatedProgenitorExists
             ! Determine the type of event. If the descendant is a subhalo, then this is an inter-tree branch jump. If the
             ! descendant is not a subhalo this is an inter-tree subhalo promotion.
             pushCount=pushCount+1
             if (nodes(j)%nodeIndex == nodes(j)%hostIndex) then
                ! Inter-tree subhalo promotion.
                self%splitForestPushTime (pushCount)=nodes(j)%      nodeTime
                self%splitForestPushTo   (pushCount)=nodes(i)%      nodeIndex
                self%splitForestPullFrom (pushCount)=nodes(i)%descendantIndex
                self%splitForestIsPrimary(pushCount)=nodeIsPrimary
                self%splitForestPushType (pushCount)=pushTypeSubhaloPromotion
                self%splitForestPushDone (pushCount)=.false.
                self%splitForestPullDone (pushCount)=.false.
             else
                ! For non-primary progenitors, detect nested subhalo hierarchy. The descendant node is a subhalo. As nested hierarchies are not currently
                ! handled, we must instead find the isolated host of the descendant and push to that node instead.
                if (.not.nodeIsPrimary) then
                   if (.not.self%warningSplitForestNestedHierarchyIssued) then
                      message='nested hierarchy in split forests detected [node '
                      message=message//nodes(j)%nodeIndex//']'
                      message=message//char(10)//'ignoring as not currently supported'
                      message=message//char(10)//'warning will not be issued again'
                      call displayMessage(message,verbosityLevelWarn)
                      self%warningSplitForestNestedHierarchyIssued=.true.
                   end if
                   do while (nodes(j)%nodeIndex /= nodes(j)%hostIndex)
                      j=self%nodeLocation(nodes(j)%hostIndex)
                   end do
                end if
                ! Inter-tree branch jump.
                k                              =self%nodeLocation(nodes(i)%hostIndex)
                self%splitForestPushTime (pushCount)=              nodes(k)%nodeTime
                self%splitForestPushTo   (pushCount)=              nodes(i)%nodeIndex
                self%splitForestPullFrom (pushCount)=              nodes(j)%nodeIndex
                self%splitForestIsPrimary(pushCount)=nodeIsPrimary
                self%splitForestPushType (pushCount)=pushTypeBranchJump
                self%splitForestPushDone (pushCount)=.false.
                self%splitForestPullDone (pushCount)=.false.
             end if
          end if
       end if
    end do
    call displayUnindent('done',verbosityLevelInfo)
    return
  end subroutine readRootNodeAffinitiesInitial

  logical function readIsOnPushList(self,node)
    !!{
    Return true if the given node is on the current ``push-to'' list of nodes for split forests.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class(mergerTreeConstructorRead), intent(inout) :: self
    class(nodeData                 ), intent(in   ) :: node

    if (allocated(self%splitForestPushTo)) then
       readIsOnPushList=any(self%splitForestPushTo == node%nodeIndex)
    else
       readIsOnPushList=.false.
    end if
    return
  end function readIsOnPushList

  logical function readIsOnPullList(self,node)
    !!{
    Return true if the given node is on the current ``pull-from'' list of nodes for split forests.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class(mergerTreeConstructorRead), intent(inout) :: self
    class(nodeData                 ), intent(in   ) :: node

    if (allocated(self%splitForestPullFrom)) then
       readIsOnPullList=any(self%splitForestPullFrom == node%nodeIndex)
    else
       readIsOnPullList=.false.
    end if
    return
  end function readIsOnPullList

  function readPushListIndex(self,node)
    !!{
    Return the index of the given node in the ``push-to'' list of nodes for split forests.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class  (mergerTreeConstructorRead), intent(inout) :: self
    integer(c_size_t                 )                :: readPushListIndex
    class  (nodeData                 ), intent(in   ) :: node
    integer(c_size_t                 )                :: i

    readPushListIndex=-1_c_size_t
    do i=1,size(self%splitForestPushTo)
       if (self%splitForestPushTo(i) == node%nodeIndex) then
          readPushListIndex=i
          exit
       end if
    end do
    return
  end function readPushListIndex

  function readPullListIndex(self,node,iPull)
    !!{
    Return the index of the given node in the ``pull-from'' list of nodes for split forests.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class  (mergerTreeConstructorRead), intent(inout) :: self
    integer(c_size_t                 )                :: readPullListIndex
    class  (nodeData                 ), intent(in   ) :: node
    integer(c_size_t                 ), intent(in   ) :: iPull
    integer(c_size_t                 )                :: i                , matchesRemaining

    readPullListIndex   =-1_c_size_t
    matchesRemaining=iPull
    do i=1,size(self%splitForestPullFrom)
       if (self%splitForestPullFrom(i) == node%nodeIndex) then
          matchesRemaining=matchesRemaining-1
          if (matchesRemaining == 0) then
             readPullListIndex   =i
             exit
          end if
       end if
    end do
    return
  end function readPullListIndex

  function readPullListCount(self,node)
    !!{
    Return the number of the given node in the ``pull-from'' list of nodes for split forests.
    !!}
    use :: Merger_Tree_Read_Importers, only : nodeData
    implicit none
    class  (mergerTreeConstructorRead), intent(inout) :: self
    integer(c_size_t                 )                :: readPullListCount
    class  (nodeData                 ), intent(in   ) :: node

    if (allocated(self%splitForestPullFrom)) then
       readPullListCount=count(self%splitForestPullFrom == node%nodeIndex)
    else
       readPullListCount=0
    end if
    return
  end function readPullListCount

  subroutine readAssignSplitForestEvents(self,nodes,nodeList)
    !!{
    Assign events to nodes if they jump between trees in a forest.
    !!}
    use :: Display                   , only : displayIndent          , displayMessage     , displayUnindent             , verbosityLevelInfo
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : nodeComponentBasic     , nodeEvent          , nodeEventBranchJumpInterTree, nodeEventSubhaloPromotionInterTree, &
          &                                   treeNode               , treeNodeList
    use :: Merger_Tree_Read_Importers, only : nodeData
    use :: Node_Events_Inter_Tree    , only : Node_Pull_From_Tree    , Node_Push_From_Tree
    use :: String_Handling           , only : operator(//)
    implicit none
    class           (mergerTreeConstructorRead)              , intent(inout), target :: self
    class           (nodeData                 ), dimension(:), intent(inout), target :: nodes
    type            (treeNodeList             ), dimension(:), intent(inout)         :: nodeList
    class           (nodeEvent                ), pointer                             :: newEvent
    class           (nodeData                 ), pointer                             :: node
    type            (treeNode                 ), pointer                             :: nodeNew          , satellite
    class           (nodeComponentBasic       ), pointer                             :: basicNew
    integer                                                                          :: iNode            , i
    integer         (c_size_t                 )                                      :: iIsolatedNode    , progenitorLocation, &
         &                                                                              iPull
    integer         (kind_int8                )                                      :: progenitorIndex
    type            (varying_string           )                                      :: message
    character       (len=12                   )                                      :: label
    logical                                                                          :: nodeIsMostMassive

    !! TODO: Handle cases where branch jumps and/or subhalo promotions are disallowed.

    call displayIndent('Assigning inter-tree events',verbosityLevelInfo)
    do iNode=1,size(nodes)
       ! Process only isolated nodes.
       if (nodes(iNode)%isolatedNodeIndex == nodeReachabilityUnreachable%ID) cycle
       ! Trace through subhalo descendants.
       node => nodes(iNode)
       newEvent => null (     )
       do while (.true.)
          if (self%isOnPushList(node)) then
             iIsolatedNode =  nodes(iNode)%isolatedNodeIndex
             if (.not.self%splitForestPushDone(self%pushListIndex(node))) then
                select case (self%splitForestPushType(self%pushListIndex(node))%ID)
                case (pushTypeSubhaloPromotion%ID)
                   allocate(nodeEventSubhaloPromotionInterTree :: newEvent)
                case (pushTypeBranchJump      %ID)
                   allocate(nodeEventBranchJumpInterTree       :: newEvent)
                end select
                call nodeList(iIsolatedNode)%node%attachEvent(newEvent)
                newEvent%time =  self%splitForestPushTime(self%pushListIndex(node))
                newEvent%node => null()
                newEvent%task => Node_Push_From_Tree
                select type (newEvent)
                type is (nodeEventSubhaloPromotionInterTree)
                   newEvent%splitForestUniqueID =  self%splitForestUniqueID
                   newEvent%pairedNodeID        =  self%splitForestPushTo(self%pushListIndex(node))
                   newEvent%mergeTimeSet        => null()
                type is (nodeEventBranchJumpInterTree      )
                   newEvent%splitForestUniqueID =  self%splitForestUniqueID
                   newEvent%pairedNodeID        =  self%splitForestPushTo(self%pushListIndex(node))
                   newEvent%mergeTimeSet        => null()
                end select
                self%splitForestPushDone(self%pushListIndex(node))=.true.
                write (label,'(f12.8)') self%splitForestPushTime(self%pushListIndex(node))
                message="Attaching push event ["
                message=message//newEvent%ID//"] to node "//node%nodeIndex//" [-->"//self%splitForestPullFrom(self%pushListIndex(node))//"] {ref:"//self%splitForestPushTo(self%pushListIndex(node))//"} at time "//label//" Gyr"
                call displayMessage(message,verbosityLevelInfo)
             end if
          end if
          if (self%isOnPullList(node)) then
             do iPull=1,self%pullListCount(node)
                if (.not.self%splitForestPullDone(self%pullListIndex(node,iPull))) then
                   select case (self%splitForestPushType(self%pullListIndex(node,iPull))%ID)
                   case (pushTypeSubhaloPromotion%ID)
                      allocate(nodeEventSubhaloPromotionInterTree :: newEvent)
                   case (pushTypeBranchJump      %ID)
                      allocate(nodeEventBranchJumpInterTree       :: newEvent)
                   case default
                      call Error_Report('unknown push type'//{introspection:location})
                   end select
                   iIsolatedNode=nodes(iNode)%isolatedNodeIndex
                   if (self%splitForestIsPrimary(self%pullListIndex(node,iPull))) then
                      ! For a subhalo promotion primary progenitor, create a temporary primary progenitor node (unless we have
                      ! previously done so) to which we attach the event. This will later be replaced with our node.
                      if (self%splitForestPushType(self%pullListIndex(node,iPull)) == pushTypeBranchJump) then
                         nodeNew => nodeList(iIsolatedNode)%node
                      else
                         allocate                                    (nodeNew)
                         call nodeList(iIsolatedNode)%node%copyNodeTo(nodeNew)
                         if (nodeNew%satelliteCount() > 0) then
                            ! Remove any satellite component from the copied node - each branch should have only a single satellite.
                            do i=nodeNew%satelliteCount(),1,-1
                               call nodeNew%satelliteRemove(i)
                            end do
                         end if
                         nodeNew%sibling                         => nodeList(iIsolatedNode)%node%firstChild
                         nodeNew%parent                          => nodeList(iIsolatedNode)%node
                         nodeNew%firstChild                      => null()
                         nodeNew%mergeTarget                     => null()
                         nodeNew%siblingMergee                   => null()
                         nodeList(iIsolatedNode)%node%firstChild => nodeNew
                         basicNew                                => nodeNew%basic()
                         call basicNew%timeSet(basicNew%time()*(1.0d0-fractionOffsetTimeClones))
                         ! Events remain attached to the original and we do not want to duplicate them.
                         nodeNew%event => null()
                         ! Any satellites are now attached to the copy.
                         nodeList(iIsolatedNode)%node%firstSatellite => null()
                         satellite => nodeNew%firstSatellite
                         do while (associated(satellite))
                            satellite%parent => nodeNew
                            satellite        => satellite%sibling
                         end do
                      end if
                      select type (newEvent)
                      type is (nodeEventSubhaloPromotionInterTree)
                         newEvent%mergeTimeSet => null()
                      type is (nodeEventBranchJumpInterTree      )
                         newEvent%mergeTimeSet => null()
                         class default
                         call Error_Report('unknown event type'//{introspection:location})
                      end select
                      call nodeNew%attachEvent(newEvent)
                   else
                      ! For a non-primary progenitor, attach the event to the primary progenitor of the node, such that our node
                      ! can later be added as a sibling. If the primary progenitor has no child, create a clone.
                      if (.not.associated(nodeList(iIsolatedNode)%node%firstChild)) then
                         allocate                                    (nodeNew)
                         call nodeList(iIsolatedNode)%node%copyNodeTo(nodeNew)
                         if (nodeNew%satelliteCount() > 0) then
                            ! Remove any satellite component from the copied node - each branch should have only a single satellite.
                            do i=nodeNew%satelliteCount(),1,-1
                               call nodeNew%satelliteRemove(i)
                            end do
                         end if
                         nodeNew%parent                          => nodeList(iIsolatedNode)%node
                         nodeNew%sibling                         => null()
                         nodeNew%firstChild                      => null()
                         nodeNew%mergeTarget                     => null()
                         nodeNew%siblingMergee                   => null()
                         nodeList(iIsolatedNode)%node%firstChild => nodeNew
                         basicNew                                => nodeNew%basic()
                         call basicNew%timeSet(basicNew%time()*(1.0d0-fractionOffsetTimeClones))
                         ! Events remain attached to the original and we do not want to duplicate them.
                         nodeNew%event => null()
                         ! Any satellites are now attached to the copy.
                         nodeList(iIsolatedNode)%node%firstSatellite => null()
                         satellite => nodeNew%firstSatellite
                         do while (associated(satellite))
                            satellite%parent => nodeNew
                            satellite        => satellite%sibling
                         end do
                      end if
                      select type (newEvent)
                      type is (nodeEventSubhaloPromotionInterTree)
                         newEvent%mergeTimeSet => readInterTreeMergeTimeSet
                         newEvent%creator      => self
                      type is (nodeEventBranchJumpInterTree      )
                         newEvent%mergeTimeSet => readInterTreeMergeTimeSet
                         newEvent%creator      => self
                      class default
                         call Error_Report('unknown event type'//{introspection:location})
                      end select
                      call nodeList(iIsolatedNode)%node%firstChild%attachEvent(newEvent)
                   end if
                   newEvent%time =  self%splitForestPushTime(self%pullListIndex(node,iPull))
                   newEvent%node => null()
                   newEvent%task => Node_Pull_From_Tree
                   select type (newEvent)
                   type is (nodeEventSubhaloPromotionInterTree)
                      newEvent%splitForestUniqueID=self%splitForestUniqueID
                      newEvent%pairedNodeID       =self%splitForestPushTo      (self%pullListIndex(node,iPull))
                      newEvent%isPrimary          =self%splitForestIsPrimary   (self%pullListIndex(node,iPull))
                   type is (nodeEventBranchJumpInterTree      )
                      newEvent%splitForestUniqueID=self%splitForestUniqueID
                      newEvent%pairedNodeID       =self%splitForestPushTo      (self%pullListIndex(node,iPull))
                      newEvent%isPrimary          =self%splitForestIsPrimary   (self%pullListIndex(node,iPull))
                   class default
                      call Error_Report('unknown event type'//{introspection:location})
                   end select
                   self%splitForestPullDone(self%pullListIndex(node,iPull))=.true.
                   write (label,'(f12.8)') self%splitForestPushTime(self%pullListIndex(node,iPull))
                   message="Attaching pull event ["
                   message=message//newEvent%ID//"] to node "//node%nodeIndex//" [<--"//self%splitForestPushTo(self%pullListIndex(node,iPull))//"] {ref:"//self%splitForestPushTo(self%pullListIndex(node,iPull))//"} at time "//label//" Gyr"
                   call displayMessage(message,verbosityLevelInfo)
                end if
             end do
          end if
          ! Is this the primary progenitor of its descendant?
          nodeIsMostMassive=.true.
          progenitorIndex  =self%descendantNodeSortIndex(node%descendantIndex)
          if (progenitorIndex > 0 .and. progenitorIndex <= size(nodes)) then
             progenitorLocation=self%descendantLocations(progenitorIndex)
             do while (nodes(progenitorLocation)%descendantIndex == node%descendantIndex)
                ! Determine progenitor status.
                if     (                                                                      &
                     &                 nodes(progenitorLocation)%nodeIndex /= node%nodeIndex  &
                     &  .and.                                                                 &
                     &   massIsGreater(nodes(progenitorLocation)           ,  node          ) &
                     & ) nodeIsMostMassive=.false.
                ! Move to the next progenitor.
                progenitorIndex=progenitorIndex-1
                if (progenitorIndex > 0) then
                   progenitorLocation=self%descendantLocations(progenitorIndex)
                else
                   exit
                end if
             end do
          end if
          ! Move to a subhalo descendant.
          if (associated(node%descendant).and.node%descendant%isSubhalo.and.nodeIsMostMassive) then
             node => node%descendant
             ! If we've reached an isolated node, stop as we will process it in another pass through the loop.
             if (node%isolatedNodeIndex /= nodeReachabilityUnreachable%ID) exit
          else
             exit
          end if
       end do
    end do
    call displayUnindent('done',verbosityLevelInfo)
    return
  end subroutine readAssignSplitForestEvents

  subroutine readInterTreeMergeTimeSet(self,nodeSatellite,nodeHost)
    !!{
    Set the merging time for a node undergoing and inter-tree transfer.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentPosition, nodeComponentSatellite, treeNode
    use :: Kepler_Orbits   , only : keplerOrbit
    use :: String_Handling , only : operator(//)
    implicit none
    class           (*                       ), intent(inout)          :: self
    type            (treeNode                ), intent(inout), target  :: nodeSatellite              , nodeHost
    type            (treeNode                )               , pointer :: nodeTarget
    class           (nodeComponentSatellite  )               , pointer :: satelliteSatellite
    class           (nodeComponentBasic      )               , pointer :: basicSatellite             , basicHost       , &
         &                                                                basicTarget
    class           (nodeComponentPosition   )               , pointer :: positionSatellite           , positionHost
    logical                                   , parameter              :: acceptUnboundOrbits=.false.
    double precision                          , dimension(3)           :: relativePosition           , relativeVelocity
    double precision                                                   :: timeUntilMerging           , radiusPericenter, &
         &                                                                radiusApocenter            , radiusVirial
    type            (keplerOrbit             )                         :: orbit
    type            (varying_string          )                         :: message

    select type (self)
    class is (mergerTreeConstructorRead)
       if (self%presetMergerTimes) then
          ! If merger times are to preset, compute subresolution merging time.
          basicSatellite     => nodeSatellite%basic    (                 )
          basicHost          => nodeHost     %basic    (                 )
          positionSatellite  => nodeSatellite%position (                 )
          positionHost       => nodeHost     %position (                 )
          satelliteSatellite => nodeSatellite%satellite(autoCreate=.true.)
          relativePosition   =  positionSatellite%position()-positionHost%position()
          relativeVelocity   =  positionSatellite%velocity()-positionHost%velocity()
          orbit              =  readOrbitConstruct(basicSatellite%mass(),basicHost%mass(),relativePosition,relativeVelocity)
          timeUntilMerging   =  self%satelliteMergingTimescales_%timeUntilMerging(nodeSatellite,orbit)
          call satelliteSatellite%timeOfMergingSet(timeUntilMerging+basicSatellite%time())
          ! Set target node.
          nodeTarget => nodeHost
          do while (.true.)
             basicTarget => nodeTarget%basic()
             if (basicTarget%time() <= satelliteSatellite%timeOfMerging() .or. .not.associated(nodeTarget%parent)) exit
             nodeTarget => nodeTarget%parent
          end do
          nodeSatellite%mergeTarget   => nodeTarget
          nodeSatellite%siblingMergee => nodeTarget   %firstMergee
          nodeTarget   %firstMergee   => nodeSatellite
          ! Assign virial orbit if necessary.
          if (self%presetOrbits) then
             ! Propagate orbit to the virial radius.
             radiusPericenter     =  orbit               %radiusPericenter (        )
             radiusApocenter      =  orbit               %radiusApocenter  (        )
             radiusVirial         =  self%darkMatterHaloScale_%radiusVirial(nodeHost)
             ! Check if the orbit intersects the virial radius.
             if     (                                                              &
                  &    radiusVirial >= radiusPericenter                            &
                  &  .and.                                                         &
                  &   (radiusVirial <= radiusApocenter  .or. .not.orbit%isBound()) &
                  &  .and.                                                         &
                  &   (.not.self%presetOrbitsBoundOnly  .or.      orbit%isBound()) &
                  & ) then
                call orbit%propagate(radiusVirial,infalling=.true.)
                ! Set the orbit.
                call satelliteSatellite%virialOrbitSet(orbit)
                ! If the satellite component supports full phase-space position, set that also.
                if (satelliteSatellite%positionIsSettable()) call satelliteSatellite%positionSet(relativePosition)
                if (satelliteSatellite%velocityIsSettable()) call satelliteSatellite%velocitySet(relativeVelocity)
             else if (self%presetOrbitsSetAll) then
                ! The given orbit does not cross the virial radius. Since all orbits must be set, choose an orbit at random.
                orbit=self%virialOrbit_%orbit(nodeSatellite,nodeHost,acceptUnboundOrbits)
                call satelliteSatellite%virialOrbitSet(orbit)
             else if (self%presetOrbitsAssertAllSet) then
                message='virial orbit could not be set for node '
                message=message//nodeSatellite%index()//char(10)
                message=message//' -> set [presetOrbitsAssertAllSet]=false to ignore this problem'//char(10)
                message=message//'    (this may lead to other problems)'
                call Error_Report(message//{introspection:location})
             end if
          end if
       else if (self%presetMergerNodes) then
          ! If merger nodes are to be set, set that now.
          nodeSatellite%mergeTarget   => nodeHost
          nodeSatellite%siblingMergee => nodeHost     %firstMergee
          nodeHost     %firstMergee   => nodeSatellite
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine readInterTreeMergeTimeSet
