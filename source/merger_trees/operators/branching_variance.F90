!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a merger tree operator which decomposes the halo-to-halo variance of progenitor and subhalo counts into
  contributions from individual splits in the trees.
  !!}

  use, intrinsic :: ISO_C_Binding             , only : c_size_t
  use            :: Kind_Numbers              , only : kind_int8
  use            :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use            :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !![
  <enumeration docformat="rst">
   <name>branchingVarianceMode</name>
   <description>
   Modes of operation of the branching variance merger tree operator.
   </description>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <visibility>public</visibility>
   <entry label="table" description="Accumulate tables of n₁"                         />
   <entry label="jumps" description="Evaluate the jumps at splits using tables of n₁" />
  </enumeration>
  !!]

  !![
  <mergerTreeOperator name="mergerTreeOperatorBranchingVariance" docformat="rst">
   <description>
   A merger tree operator which measures the split-by-split decomposition of the halo-to-halo variance of counts of
   progenitors and subhalos in merger trees. A tree built by the :cite:t:`cole_hierarchical_2000` algorithm is a Markov
   branching process in the variable :math:`\omega`, the critical overdensity for collapse scaled to the present day (as used
   by that builder). For a count statistic :math:`N`, let :math:`n_1(M,\omega)` be the expected contribution to :math:`N` of
   a single branch of mass :math:`M` at :math:`\omega`. The expectation of :math:`N` conditioned on the tree built so far is
   then a martingale, whose only jumps occur at splits, so that

   .. math::

      \mathrm{Var}(N) = \left\langle \sum_\mathrm{splits} \Delta^2 \right\rangle, \quad \Delta = \sum_i n_1(M_i,\omega) + c -
      n_1\left(\sum_i M_i,\omega\right),

   where the sum over :math:`i` runs over the progenitors created by the split, and :math:`c` is any count added by the split
   itself. Between splits :math:`n_1` also changes deterministically (the "drift", from smooth accretion and the passage of
   :math:`\omega`), such that :math:`N=n_1(\mathrm{root})+\sum\Delta+\sum\mathrm{drift}` holds exactly for every tree.

   Three count statistics are supported, selected by the characters of ``[statistics]``:

   ``A``
     the number of progenitors of mass at least :math:`m` at the redshift ``[redshiftTarget]``, :math:`z_1` (the
     progenitor mass function). A branch is counted where its segment (from a node to its progenitors) crosses
     :math:`\omega_1`, the :math:`\omega` of :math:`z_1`, with its mass there interpolated in :math:`\log M` linearly in
     :math:`\omega` along the segment. Only splits before :math:`z_1` contribute jumps, and they add no count (:math:`c=0`);
   ``B``
     the number of splits anywhere in the tree producing a non-primary progenitor of mass at least :math:`m` (the unevolved
     subhalo mass function summed over all hierarchy depths);
   ``C``
     as ``B``, but counting only splits of the main branch (direct infall onto the main progenitor).

   Each is evaluated for every threshold mass :math:`m` in ``[massThresholds]``. The primary progenitor of a node is its most
   massive progenitor, and the main branch is the chain of primary progenitors from the root.

   The function :math:`n_1` is measured empirically from the trees themselves: by the Markov property every node is a valid
   starting point, so :math:`n_1(M,\omega)` is the mean count of all nodes in a bin of :math:`(\log_{10} M,\omega)`. This
   requires two runs over sets of trees:

   ``[mode]=table``
     the count of every node is accumulated into a table, binned in :math:`\log_{10}M` (bins of width ``[binWidthMass]``
     dex spanning ``[massTableMinimum]`` to ``[massTableMaximum]``) and in :math:`\omega` (bins of width ``[binWidthOmega]``
     spanning ``[redshiftTableMinimum]`` to ``[redshiftTableMaximum]``). Bins are aligned so that the largest mass and the
     smallest :math:`\omega` fall at bin centers. Each bin accumulates the sums needed to fit a plane, :math:`n \approx a +
     b\,\delta\log_{10}M + c\,\delta\omega`, to its nodes by least squares, so that :math:`a` estimates :math:`n_1` at the bin
     center (the mean count of a bin instead estimates it at the centroid of its nodes, which is biased where the density of
     nodes varies across the bin). For statistic ``A`` only nodes before :math:`z_1` (at :math:`\omega&lt;\omega_1`) enter the
     tables, which have their own :math:`\omega` bins, of nearly the same width, whose last edge is :math:`\omega_1`. There
     :math:`n_1` is known exactly, :math:`n_1(M,\omega_1)=\theta(M-m)`, and between the last row of bin centers and
     :math:`\omega_1` it is interpolated linearly in :math:`\omega` to that value (otherwise the bins smear the step in mass at
     :math:`\omega_1`, so that :math:`\mathrm{E}[N|\omega_1]` differs from :math:`N` on average). Trees are divided into three folds by their random number
     seed (modulo 3), and a separate table is accumulated for each fold. The sums are written to the group ``[outputGroupName]`` of the output file, adding to
     any sums already present.
   ``[mode]=jumps``
     the tables are read from the group ``[outputGroupName]`` of the file ``[tableFileName]`` (the output of a ``table``
     run), together with the statistics, thresholds and binning. For each tree, :math:`\Delta` is evaluated at every split
     twice, using the tables of each of the two folds other than the tree's own, and :math:`\Delta^2` is estimated by their
     product. The noise in the two tables is independent, so it does not inflate the estimate, and no tree contributes to a
     table applied to it. Because folds are assigned by seed, a given tree always falls in the same fold, even if the base seed
     or the number of trees differs between the two runs. This requires that the seed of a tree also determines its
     realization. The :cite:t:`cole_hierarchical_2000` builder seeds the realization of each tree with the base seed plus the
     tree index, which is also the seed of the tree under the default ``[mergerTreeSeeds]=incremental``, but not under
     ``[mergerTreeSeeds]=random``: with other seed classes the base seed must be the same in the two runs. The binning is
     reconstructed from the table file, and must reproduce the bins of the ``table`` run (it does not if the cosmology or the
     critical overdensity differ between the runs).

   Nodes outside of the range of the tables are excluded from the tables, and evaluations of :math:`n_1` outside of that range
   take the value at the nearest bin center. The number of such nodes is written to the attribute ``countNodesOutside``, and
   the number of such evaluations (in ``jumps`` mode) to the attribute ``evaluationsOutside`` of each statistic and threshold.

   The tables are shared between threads: the table sums occupy :math:`80\,N_\omega N_M` bytes per fold, statistic and
   threshold (and the tables of a ``jumps`` run :math:`12\,N_\omega N_M` bytes), where :math:`N_\omega` and :math:`N_M` are
   the numbers of bins in :math:`\omega` and mass, however many threads are used. Under MPI each
   process writes its own output file, so the table sums of a ``table`` run are split between those files, and must be summed
   into one file before use in a ``jumps`` run.

   In a ``jumps`` run, per-tree quantities are written to the group ``[outputGroupName]``: the tree index, seed, weight,
   root mass and formation epoch (when the main branch first falls below half of the root mass) in the subgroup ``trees``,
   and, in a subgroup ``statistic&lt;S&gt;/threshold&lt;i&gt;`` for each statistic and threshold, the count :math:`N`,
   :math:`n_1(\mathrm{root})`, :math:`\sum\Delta`, :math:`\sum\mathrm{drift}`, :math:`\sum\Delta^2` (also split into terminal
   splits, whose non-primary progenitors have a total mass of at least :math:`m`, and structural splits; into main-branch and
   side-branch splits; and into splits at or after formation, defined below), and the residual of the identity
   :math:`N=n_1(\mathrm{root})+\sum\Delta+\sum\mathrm{drift}`. Also written are sums for a check of the Markov property (the
   mean of the count minus :math:`n_1` must vanish for roots, for primary and non-primary progenitors created at splits, and
   for nodes created by accretion alone). In both modes the relation between redshift and :math:`\omega` used is written to
   the datasets ``redshift`` and ``omega``.

   The splits at or after formation (``sumDeltaSquaredAfterFormation``) are those at :math:`\omega\le\tau`, where
   :math:`\tau` is the :math:`\omega` of the first main-branch node below half of the root mass. They include the split at
   which the main branch falls below half. Whether a split comes before :math:`\tau` depends only on the tree built up to it,
   so :math:`\tau` is a stopping time, and the mean of this sum is :math:`\mathrm{Var}(\mathrm{E}[N|\mathcal{F}_\tau])`, the
   variance of the count expected given the tree built to :math:`\tau`. The mean of the remainder of :math:`\sum\Delta^2` is
   :math:`\mathrm{E}[\mathrm{Var}(N|\mathcal{F}_\tau)]`. If the main branch never falls below half of the root mass, every
   split is included. The formation epoch written to ``trees`` is instead interpolated in the logarithm of mass between the
   main-branch nodes on either side of half of the root mass, so it lies between :math:`\tau` and the :math:`\omega` of the
   node before.

   A ``jumps`` run also writes, for each statistic and threshold:

   binned sums of :math:`\Delta^2`
     in the subgroup ``binned``, the sum of :math:`\Delta^2` over splits in bins of :math:`\omega` (of width
     ``[binWidthOmegaBinned]``, spanning the tables) and of each of :math:`x=M^\prime/M` (dataset ``massRatio``),
     :math:`M/m` (``massRelativeThreshold``), and :math:`M/M_0` (``massRelativeRoot``), where :math:`M` is the total mass of
     the progenitors created by the split, :math:`M^\prime` that of its non-primary progenitors, and :math:`M_0` the root
     mass. Mass bins are logarithmic, ``[countBinsPerDecadeBinned]`` per decade, spanning the mass range of the tables (for
     :math:`M/m`, over all thresholds); values beyond the bins are placed in the first or last bin. Each dataset is indexed
     by (:math:`\omega` bin, mass bin, class), where the classes are main-branch terminal, main-branch structural,
     side-branch terminal, and side-branch structural splits. Each tree's sum in each bin is added weighted by the volume
     weight of the tree, and the datasets with the suffix ``Squared`` sum the square of that weighted sum, from which errors
     that account for the correlation of splits within a tree follow. The bin edges (in :math:`\log_{10}` for masses) are
     written to datasets of the same names in the subgroup ``binned`` of ``[outputGroupName]``;
   the cumulative variance
     at each redshift in ``[redshiftsCumulative]`` (whose :math:`\omega` are written to the dataset ``omegaCumulative``),
     the expectation of :math:`N` given the tree built to that :math:`\omega`, :math:`\mathrm{E}[N|\omega]` (equal to
     :math:`N` for statistic ``A`` at :math:`\omega\ge\omega_1`), in the dataset
     ``cumulativeExpectation``, and the sum of :math:`\Delta^2` over the splits already made (those at :math:`\omega` no
     greater than that of the grid point, i.e. at redshifts up to it), in ``cumulativeSumDeltaSquared``, each indexed by (redshift, tree). :math:`\mathrm{E}[N|\omega]` is the count added by
     the splits already made plus :math:`n_1` summed over the branches alive at :math:`\omega`, whose masses are
     interpolated in :math:`\log M` linearly in :math:`\omega` along each segment. It is evaluated with the tables of each
     of the two folds other than the tree's own (``cumulativeExpectation`` and ``cumulativeExpectationAlternate``). The
     cumulative variance, :math:`V(\omega)=\mathrm{Var}(\mathrm{E}[N|\omega])`, which should equal the mean of the sum of
     :math:`\Delta^2`, is estimated by the covariance of the two over trees: the variance of either one alone includes the
     noise of its table. The covariance must be taken separately within each fold (trees grouped by ``treeSeed`` modulo 3)
     and the folds then combined: across folds, each table is paired with another fold's, and the mean error of each table
     biases the covariance (negatively, so that it is nonzero even at :math:`z=0`). The covariance removes only the
     second-order table noise. The first-order terms, through the correlation of the tables' errors with the tree's
     quantities and the tables' mean errors, average to zero only over realizations of the tables, so errors estimated by
     resampling trees with the tables fixed are too small; the scatter of the estimates between folds, each made with a
     different pair of tables, measures this noise (as ``mergerTreeVariance.py`` does);
   per-split values
     if ``[treeStrideSplits]`` is positive, for trees whose index is a multiple of it, in the subgroup ``splits``: for each
     split, the tree index, :math:`\Delta`, the estimate of :math:`\Delta^2`, :math:`\omega` and redshift,
     :math:`M`, the masses of the primary and non-primary progenitors, the root mass, :math:`n_1` of the branch before the
     split, the count added by the split, whether it is on the main branch, and its hierarchy depth (the number of
     non-primary links between it and the root).

   The operator must act on unevolved, unpruned trees, so it should precede any pruning operator.
   </description>
   <deepCopy>
     <increment variables="tables%countCopies" atomic="yes"/>
   </deepCopy>
  </mergerTreeOperator>
  !!]

  type :: branchingVarianceTables
     !!{RST
     The :math:`n_1` tables of a branching variance merger tree operator, shared between all copies of the operator.
     !!}
     ! Table sums, indexed by (sum, ω bin, mass bin, fold, case), accumulated in "table" mode.
     double precision, allocatable, dimension(:,:,:,:,:) :: sums
     ! Tables of n₁, indexed by (ω bin, mass bin, fold, case), and flags for empty bins, used in "jumps" mode.
     double precision, allocatable, dimension(:,:,:,:  ) :: nodeMean
     logical         , allocatable, dimension(:,:,:,:  ) :: empty
     ! The numbers of copies of the operator made, and of copies finalized.
     integer                                             :: countCopies=0, countFinalized=0
  end type branchingVarianceTables

  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorBranchingVariance
     !!{RST
     A merger tree operator which decomposes the variance of progenitor and subhalo counts into contributions from splits.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer                           :: cosmologyFunctions_       => null()
     class           (criticalOverdensityClass     ), pointer                           :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass), pointer                           :: cosmologicalMassVariance_ => null()
     type            (enumerationBranchingVarianceModeType)                             :: mode
     integer                                                                            ::                                      countStatistics          , &
          &                                                                                countThresholds                    , countCases               , &
          &                                                                                countBinsMass                      , countBinsOmega           , &
          &                                                                                nodePositionID
     type            (varying_string               )                                    :: statistics                         , outputGroupName          , &
          &                                                                                tableFileName
     double precision                                                                   :: binWidthMass                       , binWidthOmega            , &
          &                                                                                logMassEdgeMinimum                 , omegaEdgeMinimum         , &
          &                                                                                timeNow                            , massTableMinimum         , &
          &                                                                                massTableMaximum                   , redshiftTableMinimum     , &
          &                                                                                redshiftTableMaximum               , omegaSpreadMaximum
     double precision                               , allocatable, dimension(:        ) :: massThresholds                     , logMassEdges             , &
          &                                                                                omegaEdges                         , logMassCenters           , &
          &                                                                                omegaCenters
     ! The target epoch of statistic A, and the ω bins of its tables, which end at that epoch.
     double precision                                                                   :: redshiftTarget                     , omegaTarget              , &
          &                                                                                binWidthOmegaTarget
     integer                                                                            :: countBinsOmegaTarget
     double precision                               , allocatable, dimension(:        ) :: omegaEdgesTarget                   , omegaCentersTarget
     logical                                        , allocatable, dimension(:        ) :: caseIsTarget
     ! The n₁ tables, shared between all copies of the operator.
     type            (branchingVarianceTables      ), pointer                           :: tables                    => null()
     ! Per-tree quantities accumulated in "jumps" mode.
     integer         (c_size_t                     )                                    :: countTrees
     integer         (kind_int8                    ), allocatable, dimension(:        ) :: treeIndex                          , treeSeed
     double precision                               , allocatable, dimension(:        ) :: treeWeight                         , massRoot                 , &
          &                                                                                omegaFormation                     , redshiftFormation
     double precision                               , allocatable, dimension(:,:,:    ) :: treeQuantities
     ! Markov check sums, indexed by (count/sum/sum of squares, node type, case).
     double precision                               , allocatable, dimension(:,:,:    ) :: markovSums
     ! Binned sums of Δ², indexed by (ω bin, mass bin, split class, binned quantity, case), and the bins, used in "jumps" mode.
     integer                                                                            :: countBinsPerDecadeBinned           , countBinsOmegaBinned     , &
          &                                                                                treeStrideSplits                   , countCumulative
     double precision                                                                   :: binWidthOmegaBinned
     integer                                        , allocatable, dimension(:        ) :: countBinsBinned
     double precision                               , allocatable, dimension(:        ) :: logMassBinnedMinimum               , redshiftsCumulative      , &
          &                                                                                omegasCumulative
     double precision                               , allocatable, dimension(:,:,:,:,:) :: binnedSums                         , binnedSumsSquared
     ! Per-tree expectations of N conditioned on the tree built to each ω of the cumulative grid (evaluated with the tables of each
     ! of the two other folds), and sums of Δ² over splits up to that ω, indexed by (ω, case, tree).
     double precision                               , allocatable, dimension(:,:,:    ) :: cumulativeExpectation              , cumulativeSumDeltaSquared, &
          &                                                                                cumulativeExpectationAlternate
     ! Per-split values, indexed by (property, split, case).
     integer         (c_size_t                     ), allocatable, dimension(:        ) :: countSplitsStored
     double precision                               , allocatable, dimension(:,:,:    ) :: splitsReal
     integer         (kind_int8                    ), allocatable, dimension(:,:,:    ) :: splitsInteger
     integer         (kind_int8                    ), allocatable, dimension(:        ) :: evaluationsEmpty                   , evaluationsOutside
     integer         (kind_int8                    )                                    :: countSplits                        , countMultiple            , &
          &                                                                                countPrimaryMismatch               , countNodesOutside
   contains
     !![
     <methods docformat="rst">
       <method description="Return the scaled critical overdensity, :math:`\omega`, at the given time and mass." method="omega" />
       <method description="Construct the bin edges and centers of the :math:`n_1` tables." method="binsConstruct" />
       <method description="Complete the :math:`n_1` tables from their accumulated sums." method="tablesFinalize" />
       <method description="Evaluate :math:`n_1` from a table by bilinear interpolation between bin centers." method="nodeMean" />
       <method description="Find the bin containing a value, or zero if it lies outside of the table." method="binFind" />
       <method description="Return the index of the case for a given statistic and threshold." method="caseIndex" />
     </methods>
     !!]
     final     ::                        branchingVarianceDestructor
     procedure :: operatePreEvolution => branchingVarianceOperatePreEvolution
     procedure :: finalize            => branchingVarianceFinalize
     procedure :: omega               => branchingVarianceOmega
     procedure :: binsConstruct       => branchingVarianceBinsConstruct
     procedure :: tablesFinalize      => branchingVarianceTablesFinalize
     procedure :: nodeMean            => branchingVarianceNodeMean
     procedure :: binFind             => branchingVarianceBinFind
     procedure :: caseIndex           => branchingVarianceCaseIndex
  end type mergerTreeOperatorBranchingVariance

  interface mergerTreeOperatorBranchingVariance
     !!{RST
     Constructors for the branching variance merger tree operator class.
     !!}
     module procedure branchingVarianceConstructorParameters
     module procedure branchingVarianceConstructorInternal
  end interface mergerTreeOperatorBranchingVariance

  !![
  <enumeration docformat="rst">
   <name>tableSum</name>
   <description>
   Sums accumulated in each bin of the :math:`n_1` tables: the number of nodes, the sums of the count, :math:`c`, and its square,
   and the sums of :math:`\delta x`, :math:`\delta y`, :math:`\delta x^2`, :math:`\delta y^2`, :math:`\delta x \delta y`,
   :math:`\delta x c`, and :math:`\delta y c`, where :math:`\delta x` and :math:`\delta y` are the offsets of a node from the
   center of its bin in :math:`\log_{10}` (mass) and :math:`\omega`, in units of the bin widths.
   </description>
   <validator>yes</validator>
   <visibility>private</visibility>
   <indexing>1</indexing>
   <entry label="nodes"/>
   <entry label="c"    />
   <entry label="cc"   />
   <entry label="x"    />
   <entry label="y"    />
   <entry label="xx"   />
   <entry label="yy"   />
   <entry label="xy"   />
   <entry label="xc"   />
   <entry label="yc"   />
  </enumeration>
  !!]

  !![
  <enumeration docformat="rst">
   <name>treeQuantity</name>
   <description>
   Per-tree quantities written in ``jumps`` mode. The labels are the names of the datasets to which they are written, and the
   descriptions are the descriptions of those datasets.
   </description>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <visibility>private</visibility>
   <indexing>1</indexing>
   <entry label="countTotal"                    description="The count, N"                                                                                                     />
   <entry label="n1Root"                        description="The expected count, n₁, at the root"                                                                              />
   <entry label="sumDelta"                      description="The sum of the jumps, Δ, over splits"                                                                             />
   <entry label="sumDrift"                      description="The sum of the drifts over segments"                                                                              />
   <entry label="sumDeltaSquared"               description="The sum of Δ² over splits"                                                                                        />
   <entry label="sumDeltaSquaredTerminal"       description="The sum of Δ² over terminal splits (non-primary progenitor mass at least the threshold)"                          />
   <entry label="sumDeltaSquaredStructural"     description="The sum of Δ² over structural splits (non-primary progenitor mass below the threshold)"                           />
   <entry label="sumDeltaSquaredMain"           description="The sum of Δ² over main-branch splits"                                                                            />
   <entry label="sumDeltaSquaredSide"           description="The sum of Δ² over side-branch splits"                                                                            />
   <entry label="sumDeltaSquaredAfterFormation" description="The sum of Δ² over splits at or after formation (ω ≤ that of the first main-branch node below half the root mass)"/>
   <entry label="telescopeResidual"             description="The residual of the identity N = n₁(root) + ΣΔ + Σ(drift)"                                                        />
  </enumeration>
  !!]

  !![
  <enumeration docformat="rst">
   <name>markovNodeType</name>
   <description>
   Types of node for the check of the Markov property: roots, primary and non-primary progenitors created at splits, and nodes
   created by accretion alone.
   </description>
   <validator>yes</validator>
   <visibility>private</visibility>
   <indexing>1</indexing>
   <entry label="root"     />
   <entry label="primary"  />
   <entry label="secondary"/>
   <entry label="accretion"/>
  </enumeration>
  !!]

  !![
  <enumeration docformat="rst">
   <name>splitClass</name>
   <description>
   Classes of split for the binned sums of :math:`\Delta^2`: on the main branch or a side branch, and terminal (non-primary
   progenitor mass at least the threshold) or structural.
   </description>
   <validator>yes</validator>
   <visibility>private</visibility>
   <indexing>1</indexing>
   <entry label="mainTerminal"  />
   <entry label="mainStructural"/>
   <entry label="sideTerminal"  />
   <entry label="sideStructural"/>
  </enumeration>
  !!]

  !![
  <enumeration docformat="rst">
   <name>binnedQuantity</name>
   <description>
   Mass quantities in whose bins (together with :math:`\omega`) sums of :math:`\Delta^2` are accumulated. The labels are the
   names of the datasets to which they are written, and the descriptions are the descriptions of those datasets.
   </description>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <visibility>private</visibility>
   <indexing>1</indexing>
   <entry label="massRatio"             description="Sums of Δ² in bins of ω and of the mass ratio, x=M′/M, indexed by (ω bin, mass bin, class)"/>
   <entry label="massRelativeThreshold" description="Sums of Δ² in bins of ω and of M/m, indexed by (ω bin, mass bin, class)"                   />
   <entry label="massRelativeRoot"      description="Sums of Δ² in bins of ω and of M/M₀, indexed by (ω bin, mass bin, class)"                  />
  </enumeration>
  !!]

  !![
  <enumeration docformat="rst">
   <name>splitPropertyReal</name>
   <description>
   Real-valued per-split properties. The labels are the names of the datasets to which they are written, and the descriptions
   are the descriptions of those datasets.
   </description>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <visibility>private</visibility>
   <indexing>1</indexing>
   <entry label="delta"         description="The jump, Δ"                                        />
   <entry label="deltaSquared"  description="The estimate of Δ²"                                 />
   <entry label="omega"         description="ω at the split"                                     />
   <entry label="redshift"      description="Redshift of the split"                              />
   <entry label="mass"          description="Total mass of the progenitors, M [M☉]"              />
   <entry label="massPrimary"   description="Mass of the primary progenitor [M☉]"                />
   <entry label="massSecondary" description="Total mass of the non-primary progenitors, M′ [M☉]" />
   <entry label="massRoot"      description="Mass of the root of the tree, M₀ [M☉]"              />
   <entry label="n1Parent"      description="n₁ of the branch before the split"                  />
   <entry label="countAdded"    description="The count added by the split (zero for statistic A)" />
  </enumeration>
  !!]

  !![
  <enumeration docformat="rst">
   <name>splitPropertyInteger</name>
   <description>
   Integer-valued per-split properties. The labels are the names of the datasets to which they are written, and the
   descriptions are the descriptions of those datasets.
   </description>
   <decodeFunction>yes</decodeFunction>
   <validator>yes</validator>
   <visibility>private</visibility>
   <indexing>1</indexing>
   <entry label="treeIndex"      description="Index of the tree"                                                     />
   <entry label="isMainBranch"   description="1 if the split is on the main branch, 0 otherwise"                    />
   <entry label="hierarchyDepth" description="Number of non-primary links between the split and the root of the tree"/>
  </enumeration>
  !!]

  ! Number of folds into which trees are divided.
  integer         , parameter :: countFolds          =3

  ! Number of points in the tabulation of ω as a function of redshift written to the output.
  integer         , parameter :: countOmegaTabulation=10001

  ! Tolerance, in units of the bin width, for agreement between the bins of a "jumps" run and those of its table file.
  double precision, parameter :: toleranceBins       =1.0d-6

contains

  function branchingVarianceConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the branching variance merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeOperatorBranchingVariance)                              :: self
    type            (inputParameters                    ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass            ), pointer                     :: cosmologyFunctions_
    class           (criticalOverdensityClass           ), pointer                     :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), pointer                     :: cosmologicalMassVariance_
    double precision                                     , allocatable  , dimension(:) :: massThresholds           , redshiftsCumulative
    type            (enumerationBranchingVarianceModeType)                             :: mode_
    type            (varying_string                     )                              :: mode                     , statistics          , &
         &                                                                                outputGroupName          , tableFileName
    double precision                                                                   :: massTableMinimum         , massTableMaximum    , &
         &                                                                                redshiftTableMinimum     , redshiftTableMaximum, &
         &                                                                                binWidthMass             , binWidthOmega       , &
         &                                                                                binWidthOmegaBinned      , redshiftTarget
    integer                                                                            :: countBinsPerDecadeBinned , treeStrideSplits

    !![
    <inputParameter docformat="rst">
      <name>mode</name>
      <source>parameters</source>
      <description>
      The mode of operation: ``table`` to accumulate tables of :math:`n_1`, or ``jumps`` to evaluate the jumps at splits using
      tables accumulated previously.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputGroupName</name>
      <source>parameters</source>
      <defaultValue>var_str('branchingVariance')</defaultValue>
      <description>
      The name of the :term:`HDF5` group to which results are written (and, in ``jumps`` mode, from which tables are read).
      </description>
    </inputParameter>
    !!]
    mode_=enumerationBranchingVarianceModeEncode(char(mode),includesPrefix=.false.)
    select case (mode_%ID)
    case (branchingVarianceModeTable%ID)
       !![
       <inputParameter docformat="rst">
         <name>statistics</name>
         <source>parameters</source>
         <defaultValue>var_str('BC')</defaultValue>
         <description>
         The count statistics to measure, any of ``A``, ``B``, and ``C`` (see the description of this operator). Used in
         ``table`` mode only: a ``jumps`` run takes the statistics from its table file.
         </description>
       </inputParameter>
       !!]
       allocate(massThresholds(max(1,parameters%count('massThresholds',zeroIfNotPresent=.true.))))
       !![
       <inputParameter docformat="rst">
         <name>massThresholds</name>
         <source>parameters</source>
         <defaultValue>[1.0d10]</defaultValue>
         <description>
         The threshold masses, :math:`m`, for the count statistics. Used in ``table`` mode only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>massTableMinimum</name>
         <source>parameters</source>
         <description>
         The smallest mass spanned by the :math:`n_1` tables. Nodes of lower mass are excluded from the tables, and evaluations
         of :math:`n_1` below it take the value at the lowest bin center. Node masses can fall somewhat below the mass
         resolution of the trees (to :math:`0.9` times the resolution in the :cite:t:`cole_hierarchical_2000` trees used to
         validate this operator), so a value of about half the resolution is appropriate. Used in ``table`` mode only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>massTableMaximum</name>
         <source>parameters</source>
         <description>
         The largest mass spanned by the :math:`n_1` tables (normally the largest root mass). Used in ``table`` mode only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>redshiftTableMinimum</name>
         <source>parameters</source>
         <defaultValue>0.0d0</defaultValue>
         <description>
         The smallest redshift spanned by the :math:`n_1` tables (normally the redshift of the roots). Used in ``table`` mode
         only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>redshiftTableMaximum</name>
         <source>parameters</source>
         <description>
         The largest redshift spanned by the :math:`n_1` tables (normally the largest redshift reached by the trees). Used in
         ``table`` mode only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>binWidthMass</name>
         <source>parameters</source>
         <defaultValue>0.00625d0</defaultValue>
         <description>
         The width, in dex, of the mass bins of the :math:`n_1` tables. Used in ``table`` mode only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>binWidthOmega</name>
         <source>parameters</source>
         <defaultValue>0.025d0</defaultValue>
         <description>
         The width of the :math:`\omega` bins of the :math:`n_1` tables. Used in ``table`` mode only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>redshiftTarget</name>
         <source>parameters</source>
         <defaultValue>1.0d0</defaultValue>
         <description>
         The redshift, :math:`z_1`, at which progenitors are counted for statistic ``A``. It must lie within the range of the
         tables. Used in ``table`` mode only.
         </description>
       </inputParameter>
       !!]
       tableFileName           =''
       binWidthOmegaBinned     =0.0d0
       countBinsPerDecadeBinned=0
       treeStrideSplits        =0
       allocate(redshiftsCumulative(0))
    case (branchingVarianceModeJumps%ID)
       !![
       <inputParameter docformat="rst">
         <name>tableFileName</name>
         <source>parameters</source>
         <description>
         The name of the file (the output of a ``table`` run) from which tables of :math:`n_1` are read. Used in ``jumps`` mode
         only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>binWidthOmegaBinned</name>
         <source>parameters</source>
         <defaultValue>0.1d0</defaultValue>
         <description>
         The width of the :math:`\omega` bins of the binned sums of :math:`\Delta^2`. Used in ``jumps`` mode only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>countBinsPerDecadeBinned</name>
         <source>parameters</source>
         <defaultValue>10</defaultValue>
         <description>
         The number of bins per decade of the mass bins of the binned sums of :math:`\Delta^2`. Used in ``jumps`` mode only.
         </description>
       </inputParameter>
       !!]
       if (parameters%isPresent('redshiftsCumulative')) then
          allocate(redshiftsCumulative(parameters%count('redshiftsCumulative')))
       else
          allocate(redshiftsCumulative(12                                     ))
       end if
       !![
       <inputParameter docformat="rst">
         <name>redshiftsCumulative</name>
         <source>parameters</source>
         <defaultValue>[0.0d0,0.1d0,0.2d0,0.35d0,0.5d0,0.75d0,1.0d0,1.5d0,2.0d0,3.0d0,4.0d0,6.0d0]</defaultValue>
         <description>
         The redshifts at which the expectation of the count conditioned on the tree built to that epoch is evaluated, for the
         cumulative variance. Used in ``jumps`` mode only.
         </description>
       </inputParameter>
       <inputParameter docformat="rst">
         <name>treeStrideSplits</name>
         <source>parameters</source>
         <defaultValue>0</defaultValue>
         <description>
         If positive, per-split values are written for trees whose index is a multiple of this value (per-split output is large
         for high-resolution trees). If zero, no per-split values are written. Used in ``jumps`` mode only.
         </description>
       </inputParameter>
       !!]
       statistics          =''
       allocate(massThresholds(1))
       massThresholds      =0.0d0
       massTableMinimum    =0.0d0
       massTableMaximum    =0.0d0
       redshiftTableMinimum=0.0d0
       redshiftTableMaximum=0.0d0
       binWidthMass        =0.0d0
       binWidthOmega       =0.0d0
       redshiftTarget      =0.0d0
    end select
    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=mergerTreeOperatorBranchingVariance(mode_,char(statistics),massThresholds,massTableMinimum,massTableMaximum,redshiftTableMinimum,redshiftTableMaximum,redshiftTarget,binWidthMass,binWidthOmega,binWidthOmegaBinned,countBinsPerDecadeBinned,redshiftsCumulative,treeStrideSplits,char(tableFileName),char(outputGroupName),cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function branchingVarianceConstructorParameters

  function branchingVarianceConstructorInternal(mode,statistics,massThresholds,massTableMinimum,massTableMaximum,redshiftTableMinimum,redshiftTableMaximum,redshiftTarget,binWidthMass,binWidthOmega,binWidthOmegaBinned,countBinsPerDecadeBinned,redshiftsCumulative,treeStrideSplits,tableFileName,outputGroupName,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_) result(self)
    !!{RST
    Internal constructor for the branching variance merger tree operator class. In ``jumps`` mode the statistics, thresholds
    and binning of the tables are read from the table file, and the arguments giving them are ignored. In ``table`` mode the
    arguments ``binWidthOmegaBinned``, ``countBinsPerDecadeBinned``, ``redshiftsCumulative``, and ``treeStrideSplits`` are
    ignored.
    !!}
    use :: Error             , only : Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File    , hdf5Group
    use :: ISO_Varying_String, only : var_str     , char     , extract
    implicit none
    type            (mergerTreeOperatorBranchingVariance )                                      :: self
    type            (enumerationBranchingVarianceModeType), intent(in   )                       :: mode
    character       (len=*                               ), intent(in   )                       ::                            statistics              , &
         &                                                                                         tableFileName            , outputGroupName
    double precision                                      , intent(in   ), dimension(:        ) :: massThresholds           , redshiftsCumulative
    double precision                                      , intent(in   )                       :: massTableMinimum         , massTableMaximum        , &
         &                                                                                         redshiftTableMinimum     , redshiftTableMaximum    , &
         &                                                                                         binWidthMass             , binWidthOmega           , &
         &                                                                                         binWidthOmegaBinned      , redshiftTarget
    integer                                               , intent(in   )                       :: countBinsPerDecadeBinned , treeStrideSplits
    class           (cosmologyFunctionsClass             ), intent(in   ), target               :: cosmologyFunctions_
    class           (criticalOverdensityClass            ), intent(in   ), target               :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass       ), intent(in   ), target               :: cosmologicalMassVariance_
    double precision                                      , allocatable  , dimension(:,:,:,:  ) :: sums
    double precision                                      , allocatable  , dimension(:,:,:,:,:) :: sumsAll
    double precision                                      , allocatable  , dimension(:        ) :: redshiftsRead            , omegasRead
    type            (varying_string                      )                                      :: statisticsRead
    integer                                                                                     :: i                        , j                       , &
         &                                                                                         countFoldsRead           , countBinsMassRead       , &
         &                                                                                         countBinsOmegaRead       , countBinsOmegaTargetRead
    double precision                                                                            :: logMassEdgeMinimumRead   , omegaEdgeMinimumRead    , &
         &                                                                                         omegaRecomputed
    !![
    <constructorAssign variables="massTableMinimum, massTableMaximum, redshiftTableMinimum, redshiftTableMaximum, redshiftTarget, binWidthMass, binWidthOmega, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_"/>
    !!]

    self%outputGroupName=var_str(outputGroupName)
    self%tableFileName  =var_str(tableFileName  )
    self%timeNow        =self%cosmologyFunctions_%cosmicTime(1.0d0)
    self%mode           =mode
    ! The tables are shared between all copies of the operator (deep copies copy only the pointer to them). They are freed when
    ! the last copy is finalized.
    allocate(self%tables)
    select case (mode%ID)
    case (branchingVarianceModeTable%ID)
       self%statistics    =var_str(statistics)
       self%massThresholds=massThresholds
       if (any(massThresholds       <= 0.0d0               )                            ) call Error_Report('[massThresholds] must be positive'                        //{introspection:location})
       if (    massTableMinimum     <= 0.0d0                                            ) call Error_Report('[massTableMinimum] must be positive'                      //{introspection:location})
       if (    massTableMaximum     <= massTableMinimum                                 ) call Error_Report('[massTableMaximum] must exceed [massTableMinimum]'        //{introspection:location})
       if (    redshiftTableMinimum <  0.0d0                                            ) call Error_Report('[redshiftTableMinimum] must be non-negative'              //{introspection:location})
       if (    redshiftTableMaximum <= redshiftTableMinimum                             ) call Error_Report('[redshiftTableMaximum] must exceed [redshiftTableMinimum]'//{introspection:location})
       if (    binWidthMass         <= 0.0d0                 .or. binWidthOmega <= 0.0d0) call Error_Report('[binWidthMass] and [binWidthOmega] must be positive'      //{introspection:location})
       call self%binsConstruct()
       ! Binned sums, the cumulative variance, and per-split values are not used in this mode.
       self%binWidthOmegaBinned     =0.0d0
       self%countBinsPerDecadeBinned=0
       self%countBinsOmegaBinned    =0
       self%treeStrideSplits        =0
       self%countCumulative         =0
    case (branchingVarianceModeJumps%ID)
       ! Read the statistics, thresholds, and binning from the table file. HDF5 objects are scoped within blocks so that they
       ! are closed (by finalization) while the HDF5 lock is held.
       !$ call hdf5Access%set()
       block
         type(hdf5File ) :: tableFile
         type(hdf5Group) :: tableGroup
         tableFile=hdf5File(tableFileName,readOnly=.true.)
         if (.not.tableFile%hasGroup(outputGroupName)) call Error_Report('group "'//outputGroupName//'" not found in table file "'//tableFileName//'"'//{introspection:location})
         tableGroup=tableFile%openGroup(outputGroupName)
         call    tableGroup%readAttribute('statistics'          ,statisticsRead           )
         call    tableGroup%readAttribute('massThresholds'      ,self%massThresholds      )
         call    tableGroup%readAttribute('massTableMinimum'    ,self%massTableMinimum    )
         call    tableGroup%readAttribute('massTableMaximum'    ,self%massTableMaximum    )
         call    tableGroup%readAttribute('redshiftTableMinimum',self%redshiftTableMinimum)
         call    tableGroup%readAttribute('redshiftTableMaximum',self%redshiftTableMaximum)
         call    tableGroup%readAttribute('binWidthMass'        ,self%binWidthMass        )
         call    tableGroup%readAttribute('binWidthOmega'       ,self%binWidthOmega       )
         ! Table files written before statistic A was supported lack the target redshift (and so can not include statistic A).
         if (tableGroup%hasAttribute('redshiftTarget')) then
            call tableGroup%readAttribute('redshiftTarget'      ,self%redshiftTarget      )
         else
            self%redshiftTarget=0.0d0
         end if
         if (tableGroup%hasAttribute('countBinsOmegaTarget')) then
            call tableGroup%readAttribute('countBinsOmegaTarget',countBinsOmegaTargetRead )
         else
            countBinsOmegaTargetRead=0
         end if
         call    tableGroup%readAttribute('countFolds'          ,countFoldsRead           )
         call    tableGroup%readAttribute('countBinsMass'       ,countBinsMassRead        )
         call    tableGroup%readAttribute('countBinsOmega'      ,countBinsOmegaRead       )
         call    tableGroup%readAttribute('logMassEdgeMinimum'  ,logMassEdgeMinimumRead   )
         call    tableGroup%readAttribute('omegaEdgeMinimum'    ,omegaEdgeMinimumRead     )
         call    tableGroup%readDataset  ('redshift'            ,redshiftsRead            )
         call    tableGroup%readDataset  ('omega'               ,omegasRead               )
       end block
       !$ call hdf5Access%unset()
       if (countFoldsRead /= countFolds) call Error_Report('table file has the wrong number of folds'//{introspection:location})
       self%statistics=statisticsRead
       ! Construct bins outside of the HDF5 lock (evaluating ω may itself require HDF5 access, e.g. to tabulations).
       call self%binsConstruct()
       ! Check that the bins reproduce those of the table run. ω is recomputed here, so a change in the relation between redshift
       ! and ω (e.g. a different cosmology) would otherwise misalign the tables silently.
       if     (                                                                                        &
            &   self%countBinsMass        /= countBinsMassRead                                         &
            &  .or.                                                                                    &
            &   self%countBinsOmega       /= countBinsOmegaRead                                        &
            &  .or.                                                                                    &
            &   self%countBinsOmegaTarget /= countBinsOmegaTargetRead                                  &
            &  .or.                                                                                    &
            &   abs(self%logMassEdgeMinimum-logMassEdgeMinimumRead) > toleranceBins*self%binWidthMass  &
            &  .or.                                                                                    &
            &   abs(self%omegaEdgeMinimum  -omegaEdgeMinimumRead  ) > toleranceBins*self%binWidthOmega &
            & ) call Error_Report('bins do not match those of the table file - the cosmology or critical overdensity may differ from those of the table run'//{introspection:location})
       do i=1,size(redshiftsRead),max(1,(size(redshiftsRead)-1)/10)
          omegaRecomputed=self%omega(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftsRead(i))),self%massTableMaximum)
          if (abs(omegaRecomputed-omegasRead(i)) > toleranceBins*self%binWidthOmega) &
               & call Error_Report('omega(z) does not match that of the table file - the cosmology or critical overdensity may differ from those of the table run'//{introspection:location})
       end do
       deallocate(redshiftsRead,omegasRead)
       allocate(sumsAll             (tableSumCount,self%countBinsOmega,self%countBinsMass,countFolds,self%countCases))
       allocate(self%tables%nodeMean(              self%countBinsOmega,self%countBinsMass,countFolds,self%countCases))
       allocate(self%tables%empty   (              self%countBinsOmega,self%countBinsMass,countFolds,self%countCases))
       !$ call hdf5Access%set()
       block
         type(hdf5File ) :: tableFile
         type(hdf5Group) :: tableGroup, statisticGroup, caseGroup
         tableFile =hdf5File(tableFileName,readOnly=.true.)
         tableGroup=tableFile%openGroup(outputGroupName)
         do i=1,self%countStatistics
            statisticGroup=tableGroup%openGroup('statistic'//char(extract(self%statistics,i,i)))
            do j=1,self%countThresholds
               caseGroup=statisticGroup%openGroup(char(thresholdGroupName(j)))
               call caseGroup%readDataset('tableSums',sums)
               if     (                                         &
                    &   size(sums,dim=1) /= tableSumCount       &
                    &  .or.                                     &
                    &   size(sums,dim=2) /= self%countBinsOmega &
                    &  .or.                                     &
                    &   size(sums,dim=3) /= self%countBinsMass  &
                    &  .or.                                     &
                    &   size(sums,dim=4) /= countFolds          &
                    & ) call Error_Report('table sums have unexpected shape'//{introspection:location})
               sumsAll(:,:,:,:,self%caseIndex(i,j))=sums
               deallocate(sums)
            end do
         end do
       end block
       !$ call hdf5Access%unset()
       do i=1,self%countCases
          call self%tablesFinalize(sumsAll(:,:,:,:,i),i)
       end do
       deallocate(sumsAll)
       ! Construct the bins of the binned sums of Δ²: in ω, spanning the tables, and in the logarithm of each mass quantity,
       ! spanning the mass range of the tables (for M/m, over all thresholds).
       if (    binWidthOmegaBinned      <= 0.0d0 ) call Error_Report('[binWidthOmegaBinned] must be positive'     //{introspection:location})
       if (    countBinsPerDecadeBinned <= 0     ) call Error_Report('[countBinsPerDecadeBinned] must be positive'//{introspection:location})
       if (    treeStrideSplits         <  0     ) call Error_Report('[treeStrideSplits] must be non-negative'    //{introspection:location})
       if (any(redshiftsCumulative      <  0.0d0)) call Error_Report('[redshiftsCumulative] must be non-negative' //{introspection:location})
       self%binWidthOmegaBinned     =binWidthOmegaBinned
       self%countBinsPerDecadeBinned=countBinsPerDecadeBinned
       self%treeStrideSplits        =treeStrideSplits
       self%countBinsOmegaBinned    =max(ceiling((self%omegaEdges(self%countBinsOmega+1)-self%omegaEdges(1))/binWidthOmegaBinned-toleranceBins),1)
       allocate(self%countBinsBinned     (binnedQuantityCount))
       allocate(self%logMassBinnedMinimum(binnedQuantityCount))
       self%logMassBinnedMinimum(binnedQuantityMassRatio            %ID)=log10(self%massTableMinimum/       self%massTableMaximum )
       self%logMassBinnedMinimum(binnedQuantityMassRelativeThreshold%ID)=log10(self%massTableMinimum/maxval(self%massThresholds  ))
       self%logMassBinnedMinimum(binnedQuantityMassRelativeRoot     %ID)=log10(self%massTableMinimum/       self%massTableMaximum )
       self%countBinsBinned     (binnedQuantityMassRatio            %ID)=max(ceiling((log10(1.0d0                                            )-self%logMassBinnedMinimum(binnedQuantityMassRatio            %ID))*dble(countBinsPerDecadeBinned)-toleranceBins),1)
       self%countBinsBinned     (binnedQuantityMassRelativeThreshold%ID)=max(ceiling((log10(self%massTableMaximum/minval(self%massThresholds))-self%logMassBinnedMinimum(binnedQuantityMassRelativeThreshold%ID))*dble(countBinsPerDecadeBinned)-toleranceBins),1)
       self%countBinsBinned     (binnedQuantityMassRelativeRoot     %ID)=max(ceiling((log10(1.0d0                                            )-self%logMassBinnedMinimum(binnedQuantityMassRelativeRoot     %ID))*dble(countBinsPerDecadeBinned)-toleranceBins),1)
       allocate(self%binnedSums       (self%countBinsOmegaBinned,maxval(self%countBinsBinned),splitClassCount,binnedQuantityCount,self%countCases))
       allocate(self%binnedSumsSquared(self%countBinsOmegaBinned,maxval(self%countBinsBinned),splitClassCount,binnedQuantityCount,self%countCases))
       self%binnedSums       =0.0d0
       self%binnedSumsSquared=0.0d0
       ! Find ω at each redshift of the cumulative grid.
       self%countCumulative    =size(redshiftsCumulative)
       self%redshiftsCumulative=redshiftsCumulative
       allocate(self%omegasCumulative(self%countCumulative))
       do i=1,self%countCumulative
          self%omegasCumulative(i)=self%omega(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftsCumulative(i))),self%massTableMaximum)
       end do
       ! Initialize per-split storage.
       allocate(self%countSplitsStored(self%countCases))
       self%countSplitsStored=0_c_size_t
    end select
    ! Validate statistics.
    if (self%countStatistics == 0) call Error_Report('no statistics given'//{introspection:location})
    do i=1,self%countStatistics
       if (index('ABC',char(extract(self%statistics,i,i))) == 0) call Error_Report('unknown statistic "'//char(extract(self%statistics,i,i))//'" - must be one of "A", "B", or "C"'//{introspection:location})
       if (index(char(extract(self%statistics,1,i-1)),char(extract(self%statistics,i,i))) /= 0) call Error_Report('statistic "'//char(extract(self%statistics,i,i))//'" is repeated'//{introspection:location})
    end do
    ! Allocate accumulators.
    if (self%mode == branchingVarianceModeTable) then
       allocate(self%tables%sums(tableSumCount,self%countBinsOmega,self%countBinsMass,countFolds,self%countCases))
       self%tables%sums=0.0d0
    end if
    allocate(self%markovSums        (3,markovNodeTypeCount,self%countCases))
    allocate(self%evaluationsEmpty  (                      self%countCases))
    allocate(self%evaluationsOutside(                      self%countCases))
    self%markovSums          =0.0d0
    self%evaluationsEmpty    =0_kind_int8
    self%evaluationsOutside  =0_kind_int8
    self%countNodesOutside   =0_kind_int8
    self%countTrees          =0_c_size_t
    self%countSplits         =0_kind_int8
    self%countMultiple       =0_kind_int8
    self%countPrimaryMismatch=0_kind_int8
    self%omegaSpreadMaximum  =0.0d0
    !![
    <addMetaProperty component="basic" name="branchingVariancePosition" type="integer" id="self%nodePositionID" isCreator="yes"/>
    !!]
    return
  end function branchingVarianceConstructorInternal

  subroutine branchingVarianceDestructor(self)
    !!{RST
    Destructor for the branching variance merger tree operator class.
    !!}
    implicit none
    type(mergerTreeOperatorBranchingVariance), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine branchingVarianceDestructor

  function thresholdGroupName(indexThreshold) result(name)
    !!{RST
    Return the name of the output group for the given threshold (numbered from zero).
    !!}
    use :: ISO_Varying_String, only : var_str
    use :: String_Handling   , only : operator(//)
    implicit none
    type   (varying_string)                :: name
    integer                , intent(in   ) :: indexThreshold

    name=var_str('threshold')
    name=name//(indexThreshold-1)
    return
  end function thresholdGroupName

  integer function binIndexClamped(value,minimum,width,countBins)
    !!{RST
    Return the index of the bin of the given ``width``, starting from ``minimum``, that contains ``value``. Values beyond the
    ``countBins`` bins are placed in the first or last bin.
    !!}
    implicit none
    double precision, intent(in   ) :: value    , minimum, &
         &                             width
    integer         , intent(in   ) :: countBins

    binIndexClamped=max(1,min(countBins,int(floor((value-minimum)/width))+1))
    return
  end function binIndexClamped

  integer function branchingVarianceCaseIndex(self,indexStatistic,indexThreshold)
    !!{RST
    Return the index of the case for the given statistic and threshold.
    !!}
    implicit none
    class  (mergerTreeOperatorBranchingVariance), intent(in   ) :: self
    integer                                     , intent(in   ) :: indexStatistic, indexThreshold

    branchingVarianceCaseIndex=(indexStatistic-1)*self%countThresholds+indexThreshold
    return
  end function branchingVarianceCaseIndex

  double precision function branchingVarianceOmega(self,time,mass,node)
    !!{RST
    Return the scaled critical overdensity, :math:`\omega=\delta_\mathrm{c}(t,M)\sigma(M,t_0)/\sigma(M,t)`, as used by the
    :cite:t:`cole_hierarchical_2000` merger tree builder.
    !!}
    implicit none
    class           (mergerTreeOperatorBranchingVariance), intent(inout)           :: self
    double precision                                     , intent(in   )           :: time, mass
    type            (treeNode                           ), intent(inout), optional :: node

    branchingVarianceOmega=+self%criticalOverdensity_     %value       (time=time        ,mass=mass,node=node) &
         &                 *self%cosmologicalMassVariance_%rootVariance(time=self%timeNow,mass=mass          ) &
         &                 /self%cosmologicalMassVariance_%rootVariance(time=time        ,mass=mass          )
    return
  end function branchingVarianceOmega

  subroutine branchingVarianceBinsConstruct(self)
    !!{RST
    Construct the bin edges and centers of the :math:`n_1` tables. Bins are aligned so that the largest mass and the smallest
    :math:`\omega` lie at bin centers, and at least two bins are used in each direction, as needed for interpolation. The
    tables of statistic ``A`` have their own :math:`\omega` bins, of nearly the same width, whose last edge is
    :math:`\omega_1`, the :math:`\omega` of the target epoch.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char        , extract
    implicit none
    class           (mergerTreeOperatorBranchingVariance), intent(inout) :: self
    double precision                                                     :: logMassMinimum, logMassMaximum, &
         &                                                                  omegaMinimum  , omegaMaximum
    integer                                                              :: i             , j

    self%countStatistics   =len(char(self%statistics))
    self%countThresholds   =size(self%massThresholds)
    self%countCases        =self%countStatistics*self%countThresholds
    logMassMinimum         =log10(self%massTableMinimum)
    logMassMaximum         =log10(self%massTableMaximum)
    omegaMinimum           =self%omega(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftTableMinimum)),self%massTableMaximum)
    omegaMaximum           =self%omega(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftTableMaximum)),self%massTableMaximum)
    self%countBinsMass     =max(ceiling((logMassMaximum-logMassMinimum)/self%binWidthMass +0.5d0),2)
    self%countBinsOmega    =max(ceiling((omegaMaximum  -omegaMinimum  )/self%binWidthOmega+0.5d0),2)
    self%logMassEdgeMinimum=(logMassMaximum+0.5d0*self%binWidthMass )-dble(self%countBinsMass)*self%binWidthMass
    self%omegaEdgeMinimum  = omegaMinimum  -0.5d0*self%binWidthOmega
    allocate(self%logMassEdges  (self%countBinsMass +1))
    allocate(self%omegaEdges    (self%countBinsOmega+1))
    allocate(self%logMassCenters(self%countBinsMass   ))
    allocate(self%omegaCenters  (self%countBinsOmega  ))
    do i=1,self%countBinsMass +1
       self%logMassEdges(i)=self%logMassEdgeMinimum+self%binWidthMass *dble(i-1)
    end do
    do i=1,self%countBinsOmega+1
       self%omegaEdges  (i)=self%omegaEdgeMinimum  +self%binWidthOmega*dble(i-1)
    end do
    self%logMassCenters=0.5d0*(self%logMassEdges(2:self%countBinsMass +1)+self%logMassEdges(1:self%countBinsMass ))
    self%omegaCenters  =0.5d0*(self%omegaEdges  (2:self%countBinsOmega+1)+self%omegaEdges  (1:self%countBinsOmega))
    ! Identify the cases of statistic A.
    allocate(self%caseIsTarget(self%countCases))
    do i=1,self%countStatistics
       do j=1,self%countThresholds
          self%caseIsTarget(self%caseIndex(i,j))=extract(self%statistics,i,i) == 'A'
       end do
    end do
    ! Construct the ω bins of statistic A: the smallest ω lies at a bin center, as for the other statistics, and the last edge is
    ! ω₁, where n₁ is known exactly.
    if (any(self%caseIsTarget)) then
       if (self%redshiftTarget <= self%redshiftTableMinimum .or. self%redshiftTarget > self%redshiftTableMaximum) &
            & call Error_Report('[redshiftTarget] must exceed [redshiftTableMinimum] and not exceed [redshiftTableMaximum]'//{introspection:location})
       self%omegaTarget         =self%omega(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftTarget)),self%massTableMaximum)
       self%countBinsOmegaTarget=max(nint((self%omegaTarget-omegaMinimum)/self%binWidthOmega+0.5d0),2)
       self%binWidthOmegaTarget =         (self%omegaTarget-omegaMinimum)/(dble(self%countBinsOmegaTarget)-0.5d0)
       if (self%countBinsOmegaTarget > self%countBinsOmega) call Error_Report('statistic A requires more ω bins than the tables hold'//{introspection:location})
       allocate(self%omegaEdgesTarget  (self%countBinsOmegaTarget+1))
       allocate(self%omegaCentersTarget(self%countBinsOmegaTarget  ))
       do i=1,self%countBinsOmegaTarget
          self%omegaEdgesTarget(i)=omegaMinimum+self%binWidthOmegaTarget*(dble(i)-1.5d0)
       end do
       self%omegaEdgesTarget(self%countBinsOmegaTarget+1)=self%omegaTarget
       self%omegaCentersTarget=0.5d0*(self%omegaEdgesTarget(2:self%countBinsOmegaTarget+1)+self%omegaEdgesTarget(1:self%countBinsOmegaTarget))
    else
       self%omegaTarget         =huge(0.0d0)
       self%countBinsOmegaTarget=0
       self%binWidthOmegaTarget =0.0d0
       allocate(self%omegaEdgesTarget  (0))
       allocate(self%omegaCentersTarget(0))
    end if
    return
  end subroutine branchingVarianceBinsConstruct

  integer function branchingVarianceBinFind(self,edges,value)
    !!{RST
    Find the bin containing ``value`` given bin ``edges``: bins include their lower edge, and the last bin also includes its upper
    edge. Zero is returned for values outside of the table.
    !!}
    implicit none
    class           (mergerTreeOperatorBranchingVariance), intent(in   )               :: self
    double precision                                     , intent(in   ), dimension(:) :: edges
    double precision                                     , intent(in   )               :: value
    integer                                                                            :: lower, upper, middle
    !$GLC attributes unused :: self

    if (value < edges(1) .or. value > edges(size(edges))) then
       branchingVarianceBinFind=0
       return
    end if
    ! Find the number of edges less than or equal to the value by bisection.
    lower=1
    upper=size(edges)
    do while (upper-lower > 1)
       middle=(lower+upper)/2
       if (edges(middle) <= value) then
          lower=middle
       else
          upper=middle
       end if
    end do
    if (edges(upper) <= value) lower=upper
    branchingVarianceBinFind=min(lower,size(edges)-1)
    return
  end function branchingVarianceBinFind

  subroutine branchingVarianceTablesFinalize(self,sums,indexCase)
    !!{RST
    Complete the :math:`n_1` tables of one case from their accumulated sums. In each bin with at least ten nodes, spread across
    the bin in both directions, the value at the bin center is found from a least-squares plane; other populated bins use the
    mean count of their nodes. Empty bins are filled by linear interpolation in mass along their row of :math:`\omega` (extending
    the end values beyond the populated range), and empty rows are filled from the nearest populated row. The tables of statistic
    ``A`` use only their own :math:`\omega` bins, and the rows beyond them are marked empty.
    !!}
    use :: Linear_Algebra, only : determinant3x3
    implicit none
    class           (mergerTreeOperatorBranchingVariance), intent(inout)                     :: self
    double precision                                     , intent(in   ), dimension(:,:,:,:) :: sums
    integer                                              , intent(in   )                     :: indexCase
    double precision                                     , dimension(3,3)                    :: matrix
    double precision                                     , dimension(3  )                    :: vector
    logical                                              , allocatable  , dimension(:      ) :: rowPopulated
    double precision                                                                         :: determinant, countNodes
    integer                                                                                  :: i          , j              , &
         &                                                                                      h          , k              , &
         &                                                                                      jLower     , jUpper         , &
         &                                                                                      iNearest   , countBinsOmega_

    if (self%caseIsTarget(indexCase)) then
       countBinsOmega_=self%countBinsOmegaTarget
    else
       countBinsOmega_=self%countBinsOmega
    end if
    self%tables%empty   (countBinsOmega_+1:self%countBinsOmega,:,:,indexCase)=.true.
    self%tables%nodeMean(countBinsOmega_+1:self%countBinsOmega,:,:,indexCase)=0.0d0
    allocate(rowPopulated(countBinsOmega_))
    do h=1,countFolds
       do j=1,self%countBinsMass
          do i=1,countBinsOmega_
             countNodes                        =sums(tableSumNodes%ID,i,j,h)
             self%tables%empty(i,j,h,indexCase)=countNodes <= 0.0d0
             if (self%tables%empty(i,j,h,indexCase)) then
                self%tables%nodeMean(i,j,h,indexCase)=0.0d0
                cycle
             end if
             self%tables%nodeMean(i,j,h,indexCase)=sums(tableSumC%ID,i,j,h)/countNodes
             matrix(1,:)=[     countNodes         ,sums(tableSumX%ID ,i,j,h),sums(tableSumY%ID ,i,j,h)]
             matrix(2,:)=[sums(tableSumX%ID,i,j,h),sums(tableSumXX%ID,i,j,h),sums(tableSumXY%ID,i,j,h)]
             matrix(3,:)=[sums(tableSumY%ID,i,j,h),sums(tableSumXY%ID,i,j,h),sums(tableSumYY%ID,i,j,h)]
             vector     =[sums(tableSumC%ID,i,j,h),sums(tableSumXC%ID,i,j,h),sums(tableSumYC%ID,i,j,h)]
             determinant=determinant3x3(matrix)
             if (countNodes >= 10.0d0 .and. determinant > 1.0d-3*countNodes**3) then
                ! Solve for the value of the plane at the bin center by Cramer's rule.
                self%tables%nodeMean(i,j,h,indexCase)=determinant3x3(reshape([vector,matrix(:,2),matrix(:,3)],[3,3]))/determinant
             end if
          end do
       end do
       ! Fill empty bins along each row of ω.
       do i=1,countBinsOmega_
          rowPopulated(i)=.not.all(self%tables%empty(i,:,h,indexCase))
          if (.not.rowPopulated(i)) cycle
          jLower=0
          do j=1,self%countBinsMass
             if (.not.self%tables%empty(i,j,h,indexCase)) then
                ! Fill any empty bins between this and the previous populated bin.
                do k=jLower+1,j-1
                   if (jLower == 0) then
                      self%tables%nodeMean(i,k,h,indexCase)=   self%tables%nodeMean(i,j     ,h,indexCase)
                   else
                      self%tables%nodeMean(i,k,h,indexCase)=+  self%tables%nodeMean(i,jLower,h,indexCase)         &
                           &                                +(                                                    &
                           &                                  +self%tables%nodeMean(i,j     ,h,indexCase)         &
                           &                                  -self%tables%nodeMean(i,jLower,h,indexCase)         &
                           &                                 )                                                    &
                           &                                *(self%logMassCenters(k)-self%logMassCenters(jLower)) &
                           &                                /(self%logMassCenters(j)-self%logMassCenters(jLower))
                   end if
                end do
                jLower=j
             end if
          end do
          jUpper=jLower
          do k=jUpper+1,self%countBinsMass
             self%tables%nodeMean(i,k,h,indexCase)=self%tables%nodeMean(i,jUpper,h,indexCase)
          end do
       end do
       ! Fill empty rows from the nearest populated row (the lower one where two are equally near).
       if (any(rowPopulated)) then
          do i=1,countBinsOmega_
             if (rowPopulated(i)) cycle
             iNearest=0
             do k=1,countBinsOmega_
                if (.not.rowPopulated(k)) cycle
                if (iNearest == 0) then
                   iNearest=k
                else if (abs(k-i) < abs(iNearest-i)) then
                   iNearest=k
                end if
             end do
             self%tables%nodeMean(i,:,h,indexCase)=self%tables%nodeMean(iNearest,:,h,indexCase)
          end do
       end if
    end do
    return
  end subroutine branchingVarianceTablesFinalize

  double precision function branchingVarianceNodeMean(self,logMass,omega,indexFold,indexCase,emptyTouched,outside)
    !!{RST
    Evaluate :math:`n_1` from the table of the given fold and case by bilinear interpolation between bin centers, holding it
    constant beyond the outermost centers. ``emptyTouched`` is set to true if any of the bins used was empty (and so filled),
    and ``outside`` is set to true if the point lies outside of the table. For statistic ``A``, :math:`n_1` is known exactly at
    the target epoch, :math:`n_1(M,\omega_1)=\theta(M-m)`, and between the last row of bin centers and :math:`\omega_1` it is
    interpolated linearly in :math:`\omega` to that value.
    !!}
    implicit none
    class           (mergerTreeOperatorBranchingVariance), intent(in   ) :: self
    double precision                                     , intent(in   ) :: logMass        , omega
    integer                                              , intent(in   ) :: indexFold      , indexCase
    logical                                              , intent(  out) :: emptyTouched   , outside
    double precision                                                     :: x              , y        , &
         &                                                                  u              , v        , &
         &                                                                  factor
    integer                                                              :: i              , j        , &
         &                                                                  countBinsOmega_

    if (self%caseIsTarget(indexCase)) then
       outside        =self%binFind(self%logMassEdges,logMass) == 0 .or. self%binFind(self%omegaEdgesTarget,omega) == 0
       y              =(omega-self%omegaCentersTarget(1))/(self%omegaCentersTarget(2)-self%omegaCentersTarget(1))
       countBinsOmega_=self%countBinsOmegaTarget
    else
       outside        =self%binFind(self%logMassEdges,logMass) == 0 .or. self%binFind(self%omegaEdges      ,omega) == 0
       y              =(omega-self%omegaCenters      (1))/(self%omegaCenters      (2)-self%omegaCenters      (1))
       countBinsOmega_=self%countBinsOmega
    end if
    x=(logMass-self%logMassCenters(1))/(self%logMassCenters(2)-self%logMassCenters(1))
    x=max(0.0d0,min(x,dble(self%countBinsMass-1)))
    y=max(0.0d0,min(y,dble(countBinsOmega_   -1)))
    i=min(int(floor(y)),countBinsOmega_   -2)+1
    j=min(int(floor(x)),self%countBinsMass-2)+1
    u=y-dble(i-1)
    v=x-dble(j-1)
    emptyTouched=any(self%tables%empty(i:i+1,j:j+1,indexFold,indexCase))
    branchingVarianceNodeMean=+(1.0d0-u)*(1.0d0-v)*self%tables%nodeMean(i  ,j  ,indexFold,indexCase) &
         &                    +       u *(1.0d0-v)*self%tables%nodeMean(i+1,j  ,indexFold,indexCase) &
         &                    +(1.0d0-u)*       v *self%tables%nodeMean(i  ,j+1,indexFold,indexCase) &
         &                    +       u *       v *self%tables%nodeMean(i+1,j+1,indexFold,indexCase)
    if (self%caseIsTarget(indexCase)) then
       ! Interpolate to the exact value at the target epoch.
       factor=max(0.0d0,min(1.0d0,(omega-self%omegaCentersTarget(countBinsOmega_))/(self%omegaTarget-self%omegaCentersTarget(countBinsOmega_))))
       if (logMass >= log10(self%massThresholds(modulo(indexCase-1,self%countThresholds)+1))) then
          branchingVarianceNodeMean=(1.0d0-factor)*branchingVarianceNodeMean+factor
       else
          branchingVarianceNodeMean=(1.0d0-factor)*branchingVarianceNodeMean
       end if
    end if
    return
  end function branchingVarianceNodeMean

  subroutine branchingVarianceOperatePreEvolution(self,tree)
    !!{RST
    Accumulate tables of :math:`n_1`, or evaluate the jumps at splits, for the given ``tree``.
    !!}
    use            :: Error              , only : Error_Report
    use            :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, treeNode, treeNodeList
    use, intrinsic :: IEEE_Arithmetic    , only : IEEE_Value                   , IEEE_Quiet_NaN
    use            :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    use            :: ISO_Varying_String , only : char                         , extract
    implicit none
    class           (mergerTreeOperatorBranchingVariance), intent(inout), target                     :: self
    type            (mergerTree                         ), intent(inout), target                     :: tree
    type            (mergerTree                         ), pointer                                   :: treeCurrent
    type            (treeNode                           ), pointer                                   :: node                            , nodeChild
    class           (nodeComponentBasic                 ), pointer                                   :: basic
    type            (mergerTreeWalkerIsolatedNodes      )                                            :: treeWalker
    type            (treeNodeList                       ), allocatable  , dimension(:              ) :: nodes
    integer                                              , allocatable  , dimension(:              ) :: parent                          , primaryChild                     , &
         &                                                                                              countChildren                   , hierarchyDepth
    double precision                                     , allocatable  , dimension(:              ) :: mass                            , logMass                          , &
         &                                                                                              omega                           , massChildren                     , &
         &                                                                                              massSecondary                   , omegaEnd                         , &
         &                                                                                              countAdded                      , countNode                        , &
         &                                                                                              delta                           , deltaSquared                     , &
         &                                                                                              drift                           , omegaSplit                       , &
         &                                                                                              countOwn
    logical                                              , allocatable  , dimension(:              ) :: isPrimary                       , isMain                           , &
         &                                                                                              isActive                        , crossesTarget
    double precision                                     , allocatable  , dimension(:,:            ) :: deltaVariant                    , nodeMeanNode                     , &
         &                                                                                              nodeMeanEnd
    double precision                                     , allocatable  , dimension(:,:,:          ) :: treeQuantitiesTmp               , cumulativeExpectationTmp         , &
         &                                                                                              cumulativeSumDeltaSquaredTmp    , cumulativeExpectationAlternateTmp
    double precision                                     , allocatable  , dimension(:,:,:,:        ) :: binnedTree
    double precision                                     , dimension   (binnedQuantityCount         ) :: logMassBinned
    double precision                                     , dimension   (tableSumCount               ) :: contributions
    integer         (kind_int8                          ), allocatable  , dimension(:              ) :: treeIndexTmp                    , treeSeedTmp
    double precision                                     , allocatable  , dimension(:              ) :: treeWeightTmp                   , massRootTmp                      , &
         &                                                                                              omegaFormationTmp               , redshiftFormationTmp
    integer         (c_size_t                           ), parameter                                 :: treeCountIncrement=1000_c_size_t
    integer                                                                                          :: countNodes                      , i                                , &
         &                                                                                              p                               , iStatistic                       , &
         &                                                                                              iThreshold                      , iCase                            , &
         &                                                                                              fold                            , variant                          , &
         &                                                                                              foldTable                       , iBinMass                         , &
         &                                                                                              iBinOmega                       , nodeType                         , &
         &                                                                                              j                               , k                                , &
         &                                                                                              iBin                            , splitClass_
    character       (len=1                              )                                            :: statistic
    double precision                                                                                 :: massThreshold                   , dx                               , &
         &                                                                                              dy                              , omegaFormation_                  , &
         &                                                                                              redshiftFormation_              , factor                           , &
         &                                                                                              massMaximum                     , residual                         , &
         &                                                                                              expectation                     , sumBefore                        , &
         &                                                                                              expectationAlternate            , logMassAlive                     , &
         &                                                                                              omegaStopping
    logical                                                                                          :: emptyTouched                    , outside                          , &
         &                                                                                              storeSplits
    integer         (kind_int8                          )                                            :: seed

    treeCurrent => tree
    do while (associated(treeCurrent))
       if (.not.associated(treeCurrent%nodeBase)) then
          treeCurrent => treeCurrent%nextTree
          cycle
       end if
       ! Assign the tree to a fold by its random number seed.
       if (.not.associated(treeCurrent%randomNumberGenerator_)) call Error_Report('tree has no random number generator, so can not be assigned to a fold'//{introspection:location})
       seed=treeCurrent%randomNumberGenerator_%seed()
       fold=int(modulo(seed,int(countFolds,kind_int8)))+1
       ! Count nodes, recording the position of each in a post-order (depth-first, progenitors before descendants) walk.
       countNodes=0
       treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
       do while (treeWalker%next(node))
          countNodes =  countNodes+1
          basic      => node%basic()
          call basic%integerRank0MetaPropertySet(self%nodePositionID,countNodes)
       end do
       allocate(nodes         (countNodes))
       allocate(parent        (countNodes))
       allocate(primaryChild  (countNodes))
       allocate(countChildren (countNodes))
       allocate(mass          (countNodes))
       allocate(logMass       (countNodes))
       allocate(omega         (countNodes))
       allocate(massChildren  (countNodes))
       allocate(massSecondary (countNodes))
       allocate(omegaSplit    (countNodes))
       allocate(omegaEnd      (countNodes))
       allocate(isPrimary     (countNodes))
       allocate(isMain        (countNodes))
       allocate(isActive      (countNodes))
       allocate(hierarchyDepth(countNodes))
       ! Gather node properties.
       i=0
       treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
       do while (treeWalker%next(node))
          i              =  i+1
          nodes(i)%node  => node
          basic          => node%basic()
          mass   (i)     =  basic%mass()
          logMass(i)     =  log10(mass(i))
          omega  (i)     =  self%omega(basic%time(),mass(i),node)
          if     (                                                 &
               &   self%binFind(self%omegaEdges  ,omega  (i)) == 0 &
               &  .or.                                             &
               &   self%binFind(self%logMassEdges,logMass(i)) == 0 &
               & ) self%countNodesOutside=self%countNodesOutside+1_kind_int8
          if (associated(node%parent)) then
             basic     => node%parent%basic()
             parent(i) =  basic%integerRank0MetaPropertyGet(self%nodePositionID)
          else
             parent(i) =  0
          end if
       end do
       ! Identify the primary (most massive) progenitor of each node, and sum progenitor masses.
       primaryChild =0
       countChildren=0
       massChildren =0.0d0
       isPrimary    =.false.
       do i=1,countNodes
          countChildren(i)=0
          massMaximum     =-huge(0.0d0)
          nodeChild       => nodes(i)%node%firstChild
          do while (associated(nodeChild))
             basic           => nodeChild%basic()
             p               =  basic%integerRank0MetaPropertyGet(self%nodePositionID)
             countChildren(i)=  countChildren(i)+1
             massChildren (i)=  massChildren (i)+mass(p)
             if (mass(p) >= massMaximum) then
                massMaximum    =mass(p)
                primaryChild(i)=     p
             end if
             nodeChild => nodeChild%sibling
          end do
          if (primaryChild(i) > 0) then
             isPrimary(primaryChild(i))=.true.
             if (.not.nodes(primaryChild(i))%node%isPrimaryProgenitor()) self%countPrimaryMismatch=self%countPrimaryMismatch+1_kind_int8
          end if
       end do
       ! Masses of non-primary progenitors, and the ω of each split (and of the end of each segment). Nodes with no progenitors
       ! are given an ω of the split that is larger than any other, so that comparisons with it are well-defined (NaN would raise
       ! a floating point exception).
       massSecondary=0.0d0
       omegaSplit   =huge(0.0d0)
       omegaEnd     =IEEE_Value(0.0d0,IEEE_Quiet_NaN)
       do i=1,countNodes
          if (parent(i) > 0 .and. .not.isPrimary(i)) massSecondary(parent(i))=massSecondary(parent(i))+mass(i)
          if (primaryChild(i) > 0) then
             omegaSplit(i) =  omega(primaryChild(i))
             basic         => nodes(primaryChild(i))%node%basic()
             omegaEnd  (i) =  self%omega(basic%time(),massChildren(i),nodes(i)%node)
          end if
          if (countChildren(i) >= 2) then
             self%countSplits=self%countSplits+1_kind_int8
             if (countChildren(i) > 2) self%countMultiple=self%countMultiple+1_kind_int8
          end if
       end do
       ! Check that the progenitors of each split share a single ω.
       do i=1,countNodes
          if (parent(i) > 0) self%omegaSpreadMaximum=max(self%omegaSpreadMaximum,abs(omega(i)-omegaSplit(parent(i))))
       end do
       ! Main-branch membership: walk from the root (the last node of the walk) toward the leaves.
       ! Hierarchy depth (the number of non-primary links between a node and the root) is found in the same walk.
       isMain                    =.false.
       isMain        (countNodes)=.true.
       hierarchyDepth(countNodes)=0
       do i=countNodes-1,1,-1
          isMain        (i)=isMain        (parent(i)) .and.     isPrimary(i)
          hierarchyDepth(i)=hierarchyDepth(parent(i))+merge(0,1,isPrimary(i))
       end do
       ! Formation epoch: when the main branch first falls below half of the root mass, interpolating ω linearly in the logarithm
       ! of mass. Also found is the stopping time, the ω of the first main-branch node below half of the root mass. Unlike the
       ! interpolated formation epoch, whether a split precedes it depends on nothing beyond it in the tree. If the main branch
       ! never falls below half of the root mass, the stopping time is taken to lie beyond every split.
       omegaStopping     =huge(0.0d0)
       omegaFormation_   =IEEE_Value(0.0d0,IEEE_Quiet_NaN)
       redshiftFormation_=IEEE_Value(0.0d0,IEEE_Quiet_NaN)
       i=countNodes
       do while (primaryChild(i) > 0)
          p=primaryChild(i)
          if (mass(p) < 0.5d0*mass(countNodes)) then
             factor            =+(log(0.5d0*mass(countNodes))-log(mass(i))) &
                  &             /(log(      mass(p         ))-log(mass(i)))
             omegaFormation_   =omega(i)+factor*(omega(p)-omega(i))
             redshiftFormation_=+          redshiftOfNode(nodes(i)%node)  &
                  &             +factor*(+ redshiftOfNode(nodes(p)%node)  &
                  &                      - redshiftOfNode(nodes(i)%node))
             omegaStopping     =omega(p)
             exit
          end if
          i=p
       end do
       ! Store per-tree quantities.
       if (self%mode == branchingVarianceModeJumps) then
          if (.not.allocated(self%treeIndex)) then
             allocate(self%treeIndex                     (                                       treeCountIncrement))
             allocate(self%treeSeed                      (                                       treeCountIncrement))
             allocate(self%treeWeight                    (                                       treeCountIncrement))
             allocate(self%massRoot                      (                                       treeCountIncrement))
             allocate(self%omegaFormation                (                                       treeCountIncrement))
             allocate(self%redshiftFormation             (                                       treeCountIncrement))
             allocate(self%treeQuantities                (     treeQuantityCount,self%countCases,treeCountIncrement))
             allocate(self%cumulativeExpectation         (self%countCumulative  ,self%countCases,treeCountIncrement))
             allocate(self%cumulativeSumDeltaSquared     (self%countCumulative  ,self%countCases,treeCountIncrement))
             allocate(self%cumulativeExpectationAlternate(self%countCumulative  ,self%countCases,treeCountIncrement))
          else if (self%countTrees >= size(self%treeIndex,kind=c_size_t)) then
             call move_alloc(self%treeIndex                ,treeIndexTmp                )
             call move_alloc(self%treeSeed                 ,treeSeedTmp                 )
             call move_alloc(self%treeWeight               ,treeWeightTmp               )
             call move_alloc(self%massRoot                 ,massRootTmp                 )
             call move_alloc(self%omegaFormation           ,omegaFormationTmp           )
             call move_alloc(self%redshiftFormation        ,redshiftFormationTmp        )
             call move_alloc(self%treeQuantities           ,treeQuantitiesTmp           )
             call move_alloc(self%cumulativeExpectation    ,cumulativeExpectationTmp    )
             call move_alloc(self%cumulativeSumDeltaSquared,cumulativeSumDeltaSquaredTmp)
             call move_alloc(self%cumulativeExpectationAlternate,cumulativeExpectationAlternateTmp)
             allocate(self%treeIndex                     (                                     size(treeIndexTmp)+treeCountIncrement))
             allocate(self%treeSeed                      (                                     size(treeIndexTmp)+treeCountIncrement))
             allocate(self%treeWeight                    (                                     size(treeIndexTmp)+treeCountIncrement))
             allocate(self%massRoot                      (                                     size(treeIndexTmp)+treeCountIncrement))
             allocate(self%omegaFormation                (                                     size(treeIndexTmp)+treeCountIncrement))
             allocate(self%redshiftFormation             (                                     size(treeIndexTmp)+treeCountIncrement))
             allocate(self%treeQuantities                (treeQuantityCount   ,self%countCases,size(treeIndexTmp)+treeCountIncrement))
             allocate(self%cumulativeExpectation         (self%countCumulative,self%countCases,size(treeIndexTmp)+treeCountIncrement))
             allocate(self%cumulativeSumDeltaSquared     (self%countCumulative,self%countCases,size(treeIndexTmp)+treeCountIncrement))
             allocate(self%cumulativeExpectationAlternate(self%countCumulative,self%countCases,size(treeIndexTmp)+treeCountIncrement))
             self%treeIndex                     (    1:size(treeIndexTmp))=treeIndexTmp
             self%treeSeed                      (    1:size(treeIndexTmp))=treeSeedTmp
             self%treeWeight                    (    1:size(treeIndexTmp))=treeWeightTmp
             self%massRoot                      (    1:size(treeIndexTmp))=massRootTmp
             self%omegaFormation                (    1:size(treeIndexTmp))=omegaFormationTmp
             self%redshiftFormation             (    1:size(treeIndexTmp))=redshiftFormationTmp
             self%treeQuantities                (:,:,1:size(treeIndexTmp))=treeQuantitiesTmp
             self%cumulativeExpectation         (:,:,1:size(treeIndexTmp))=cumulativeExpectationTmp
             self%cumulativeSumDeltaSquared     (:,:,1:size(treeIndexTmp))=cumulativeSumDeltaSquaredTmp
             self%cumulativeExpectationAlternate(:,:,1:size(treeIndexTmp))=cumulativeExpectationAlternateTmp
             deallocate(treeIndexTmp                     )
             deallocate(treeSeedTmp                      )
             deallocate(treeWeightTmp                    )
             deallocate(massRootTmp                      )
             deallocate(omegaFormationTmp                )
             deallocate(redshiftFormationTmp             )
             deallocate(treeQuantitiesTmp                )
             deallocate(cumulativeExpectationTmp         )
             deallocate(cumulativeSumDeltaSquaredTmp     )
             deallocate(cumulativeExpectationAlternateTmp)
          end if
          self%countTrees                        =self%countTrees+1_c_size_t
          self%treeIndex        (self%countTrees)=treeCurrent%index
          self%treeSeed         (self%countTrees)=seed
          self%treeWeight       (self%countTrees)=treeCurrent%volumeWeight
          self%massRoot         (self%countTrees)=mass(countNodes)
          self%omegaFormation   (self%countTrees)=omegaFormation_
          self%redshiftFormation(self%countTrees)=redshiftFormation_
          allocate(nodeMeanNode(countNodes,2))
          allocate(nodeMeanEnd (countNodes,2))
          allocate(delta       (countNodes  ))
          allocate(deltaSquared(countNodes  ))
          allocate(drift       (countNodes  ))
          allocate(deltaVariant(countNodes,2))
          allocate(binnedTree  (self%countBinsOmegaBinned,maxval(self%countBinsBinned),splitClassCount,binnedQuantityCount))
          storeSplits=self%treeStrideSplits > 0
          if (storeSplits) storeSplits=modulo(treeCurrent%index,int(self%treeStrideSplits,kind_int8)) == 0_kind_int8
       end if
       allocate(countAdded(countNodes))
       allocate(countOwn  (countNodes))
       allocate(countNode (countNodes))
       ! Segments crossing the target epoch of statistic A.
       allocate(crossesTarget(countNodes))
       crossesTarget=omega < self%omegaTarget .and. countChildren > 0 .and. omegaSplit >= self%omegaTarget
       ! Iterate over statistics and thresholds.
       do iStatistic=1,self%countStatistics
          statistic=char(extract(self%statistics,iStatistic,iStatistic))
          do iThreshold=1,self%countThresholds
             iCase        =self%caseIndex(iStatistic,iThreshold)
             massThreshold=self%massThresholds(iThreshold)
             ! Count added by each split (non-primary progenitors above threshold), and the count of each node: summed over its
             ! subtree (statistics A and B) or its own primary-progenitor chain (statistic C). The walk order places progenitors
             ! before their descendants.
             countAdded=0.0d0
             do i=1,countNodes
                if (parent(i) > 0 .and. .not.isPrimary(i) .and. mass(i) >= massThreshold) countAdded(parent(i))=countAdded(parent(i))+1.0d0
             end do
             if (statistic == 'A') then
                ! A segment, from a node to its progenitors, crosses the target epoch if ω of the node is below ω₁ and that of its
                ! progenitors is not. It counts if its mass there, interpolated in log M linearly in ω between the node and the
                ! summed mass of its progenitors, is at least the threshold. Splits add no count.
                countOwn=0.0d0
                do i=1,countNodes
                   if (crossesTarget(i)) then
                      logMassAlive=+logMass(i)                          &
                           &       +(log10(massChildren(i))-logMass(i)) &
                           &       *(self%omegaTarget      -omega  (i)) &
                           &       /(omegaSplit         (i)-omega  (i))
                      if (logMassAlive >= log10(massThreshold)) countOwn(i)=1.0d0
                   end if
                end do
             else
                countOwn=countAdded
             end if
             countNode=countOwn
             do i=1,countNodes-1
                if (statistic /= 'C' .or. isPrimary(i)) countNode(parent(i))=countNode(parent(i))+countNode(i)
             end do
             select case (self%mode%ID)
             case (branchingVarianceModeTable%ID)
                ! Accumulate the count of every node into the table of this tree's fold. For statistic A only nodes before the target
                ! epoch (at ω below ω₁) are included, in the ω bins of that statistic.
                do i=1,countNodes
                   iBinMass =self%binFind(self%logMassEdges,logMass(i))
                   if (statistic == 'A') then
                      if (omega(i) >= self%omegaTarget) cycle
                      iBinOmega=self%binFind(self%omegaEdgesTarget,omega(i))
                   else
                      iBinOmega=self%binFind(self%omegaEdges      ,omega(i))
                   end if
                   if (iBinOmega == 0 .or. iBinMass == 0) cycle
                   dx=(logMass(i)-self%logMassCenters(iBinMass))/self%binWidthMass
                   if (statistic == 'A') then
                      dy=(omega(i)-self%omegaCentersTarget(iBinOmega))/self%binWidthOmegaTarget
                   else
                      dy=(omega(i)-self%omegaCenters      (iBinOmega))/self%binWidthOmega
                   end if
                   contributions=[                             &
                        &         1.0d0                      , &
                        &                     countNode(i)   , &
                        &                     countNode(i)**2, &
                        &         dx                         , &
                        &         dy                         , &
                        &         dx**2                      , &
                        &               dy**2                , &
                        &         dx   *dy                   , &
                        &         dx         *countNode(i)   , &
                        &               dy   *countNode(i)     &
                        &        ]
                   ! The sums are shared between threads, so are updated atomically.
                   do k=1,tableSumCount
                      !$omp atomic
                      self%tables%sums(k,iBinOmega,iBinMass,fold,iCase)=self%tables%sums(k,iBinOmega,iBinMass,fold,iCase)+contributions(k)
                   end do
                end do
             case (branchingVarianceModeJumps%ID)
                ! Splits contributing jumps: all splits for statistic B, those of the main branch for statistic C, and those before
                ! the target epoch for statistic A.
                select case (statistic)
                case ('A')
                   isActive=countChildren >= 2 .and. omegaSplit < self%omegaTarget
                case ('B')
                   isActive=countChildren >= 2
                case ('C')
                   isActive=countChildren >= 2 .and. isMain
                end select
                ! Evaluate the jumps with the tables of each of the two other folds. The first variant also gives n₁ for the
                ! drift, the per-tree sums, the Markov check, and the cumulative variance.
                do variant=1,2
                   foldTable=modulo(fold+variant-1,countFolds)+1
                   do i=1,countNodes
                      ! For statistic A, n₁ is needed only before the target epoch.
                      if (statistic == 'A' .and. omega(i) >= self%omegaTarget) then
                         nodeMeanNode(i,variant)=0.0d0
                         cycle
                      end if
                      nodeMeanNode(i,variant)=self%nodeMean(logMass(i),omega(i),foldTable,iCase,emptyTouched,outside)
                      if (variant == 1 .and. emptyTouched) self%evaluationsEmpty  (iCase)=self%evaluationsEmpty  (iCase)+1_kind_int8
                      if (variant == 1 .and. outside     ) self%evaluationsOutside(iCase)=self%evaluationsOutside(iCase)+1_kind_int8
                   end do
                   nodeMeanEnd(:,variant)=0.0d0
                   do i=1,countNodes
                      if (countChildren(i) == 0) cycle
                      ! For statistic A, n₁ is not needed at the end of a segment which reaches the target epoch.
                      if (statistic == 'A' .and. omegaSplit(i) >= self%omegaTarget) cycle
                      nodeMeanEnd(i,variant)=self%nodeMean(log10(massChildren(i)),omegaEnd(i),foldTable,iCase,emptyTouched,outside)
                      if (variant == 1 .and. emptyTouched) self%evaluationsEmpty  (iCase)=self%evaluationsEmpty  (iCase)+1_kind_int8
                      if (variant == 1 .and. outside     ) self%evaluationsOutside(iCase)=self%evaluationsOutside(iCase)+1_kind_int8
                   end do
                   deltaVariant(:,variant)=0.0d0
                   do i=1,countNodes
                      if (.not.isActive(i)) cycle
                      if (statistic == 'C') then
                         ! Role-dependent n₁: only the primary progenitor carries the main branch forward.
                         deltaVariant(i,variant)=nodeMeanNode(primaryChild(i),variant)
                      else
                         nodeChild => nodes(i)%node%firstChild
                         do while (associated(nodeChild))
                            basic                  => nodeChild%basic()
                            deltaVariant(i,variant)=  deltaVariant(i,variant)+nodeMeanNode(basic%integerRank0MetaPropertyGet(self%nodePositionID),variant)
                            nodeChild              => nodeChild%sibling
                         end do
                      end if
                      if (statistic /= 'A') deltaVariant(i,variant)=deltaVariant(i,variant)+countAdded(i)
                      deltaVariant(i,variant)=deltaVariant(i,variant)-nodeMeanEnd(i,variant)
                   end do
                   if (variant == 1) then
                      ! Deterministic drift of n₁ along each segment, from the node to the end of the segment. For statistic A, a
                      ! segment which crosses the target epoch ends with its own count, and segments after it have no drift.
                      do i=1,countNodes
                         if      (statistic == 'C' .and. .not.isMain       (i)                    ) then
                            drift(i)=0.0d0
                         else if (statistic == 'A' .and.      omega        (i) >= self%omegaTarget) then
                            drift(i)=0.0d0
                         else if (statistic == 'A' .and.      crossesTarget(i)                    ) then
                            drift(i)=countOwn   (i  )-nodeMeanNode(i,1)
                         else if (                            countChildren(i) >  0               ) then
                            drift(i)=nodeMeanEnd(i,1)-nodeMeanNode(i,1)
                         else
                            drift(i)=                -nodeMeanNode(i,1)
                         end if
                      end do
                      ! Markov check: by the Markov property the mean of (count - n₁) must vanish for every type of node.
                      do i=1,countNodes
                         if (nodeMeanNode(i,1) <= 0.0d0) cycle
                         if (parent(i) == 0) then
                            nodeType=markovNodeTypeRoot%ID
                         else if (.not.isPrimary(i)) then
                            nodeType=markovNodeTypeSecondary%ID
                         else if (countChildren(parent(i)) >= 2) then
                            nodeType=markovNodeTypePrimary%ID
                         else
                            nodeType=markovNodeTypeAccretion%ID
                         end if
                         residual=countNode(i)-nodeMeanNode(i,1)
                         self%markovSums(:,nodeType,iCase)=self%markovSums(:,nodeType,iCase)+[1.0d0,residual,residual**2]
                      end do
                      ! Store n₁ at the root.
                      self%treeQuantities(treeQuantityN1Root%ID,iCase,self%countTrees)=nodeMeanNode(countNodes,1)
                   end if
                end do
                delta       =deltaVariant(:,1)
                deltaSquared=deltaVariant(:,1)*deltaVariant(:,2)
                ! Per-tree sums.
                self   %treeQuantities(treeQuantityCountTotal                   %ID,iCase,self%countTrees)=countNode(countNodes)
                self   %treeQuantities(treeQuantitySumDelta                     %ID,iCase,self%countTrees)=sum(delta                                                  )
                self   %treeQuantities(treeQuantitySumDrift                     %ID,iCase,self%countTrees)=sum(drift                                                  )
                self   %treeQuantities(treeQuantitySumDeltaSquared              %ID,iCase,self%countTrees)=sum(deltaSquared                                           )
                self   %treeQuantities(treeQuantitySumDeltaSquaredTerminal      %ID,iCase,self%countTrees)=sum(deltaSquared,mask=      massSecondary >= massThreshold )
                self   %treeQuantities(treeQuantitySumDeltaSquaredStructural    %ID,iCase,self%countTrees)=sum(deltaSquared,mask=.not.(massSecondary >= massThreshold))
                self   %treeQuantities(treeQuantitySumDeltaSquaredMain          %ID,iCase,self%countTrees)=sum(deltaSquared,mask=      isMain                         )
                self   %treeQuantities(treeQuantitySumDeltaSquaredSide          %ID,iCase,self%countTrees)=sum(deltaSquared,mask=.not. isMain                         )
                self   %treeQuantities(treeQuantitySumDeltaSquaredAfterFormation%ID,iCase,self%countTrees)=sum(deltaSquared,mask=      omegaSplit <= omegaStopping    )
                self   %treeQuantities(treeQuantityTelescopeResidual            %ID,iCase,self%countTrees)=+self%treeQuantities(treeQuantityCountTotal%ID,iCase,self%countTrees) &
                     &                                                                                     -self%treeQuantities(treeQuantityN1Root    %ID,iCase,self%countTrees) &
                     &                                                                                     -self%treeQuantities(treeQuantitySumDelta  %ID,iCase,self%countTrees) &
                     &                                                                                     -self%treeQuantities(treeQuantitySumDrift  %ID,iCase,self%countTrees)
                ! Binned sums of Δ². The sums of this tree are added weighted by its volume weight, and their squares are accumulated
                ! for errors (splits within a tree are correlated, so errors are found from the scatter between trees).
                binnedTree=0.0d0
                do i=1,countNodes
                   if (.not.isActive(i)) cycle
                   if (isMain(i)) then
                      if (massSecondary(i) >= massThreshold) then
                         splitClass_=splitClassMainTerminal%ID
                      else
                         splitClass_=splitClassMainStructural%ID
                      end if
                   else
                      if (massSecondary(i) >= massThreshold) then
                         splitClass_=splitClassSideTerminal%ID
                      else
                         splitClass_=splitClassSideStructural%ID
                      end if
                   end if
                   iBinOmega=binIndexClamped(omegaSplit(i),self%omegaEdges(1),self%binWidthOmegaBinned,self%countBinsOmegaBinned)
                   logMassBinned(binnedQuantityMassRatio            %ID)=log10(massSecondary(i)/massChildren(i)  )
                   logMassBinned(binnedQuantityMassRelativeThreshold%ID)=log10(massChildren (i)/massThreshold    )
                   logMassBinned(binnedQuantityMassRelativeRoot     %ID)=log10(massChildren (i)/mass(countNodes))
                   do j=1,binnedQuantityCount
                      iBin                                    =binIndexClamped(logMassBinned(j),self%logMassBinnedMinimum(j),1.0d0/dble(self%countBinsPerDecadeBinned),self%countBinsBinned(j))
                      binnedTree(iBinOmega,iBin,splitClass_,j)=binnedTree(iBinOmega,iBin,splitClass_,j)+deltaSquared(i)
                   end do
                end do
                self%binnedSums       (:,:,:,:,iCase)=+self%binnedSums       (:,:,:,:,iCase) &
                     &                                + treeCurrent%volumeWeight*binnedTree
                self%binnedSumsSquared(:,:,:,:,iCase)=+self%binnedSumsSquared(:,:,:,:,iCase) &
                     &                                +(treeCurrent%volumeWeight*binnedTree)**2
                ! Cumulative variance: at each ω of the grid, the expectation of N conditioned on the tree built to that ω is the
                ! count added by splits already made plus n₁ summed over the branches alive at ω (with masses interpolated in log M
                ! linearly in ω along each segment). It is evaluated with the tables of each of the two other folds: as for Δ², the
                ! covariance of the two over trees estimates the variance without the noise of the tables. Also record the sum of
                ! Δ² over splits already made.
                do k=1,self%countCumulative
                   ! For statistic A the count is fully determined at the target epoch and beyond.
                   if (statistic == 'A' .and. self%omegasCumulative(k) >= self%omegaTarget) then
                      self%cumulativeExpectation         (k,iCase,self%countTrees)=countNode(countNodes)
                      self%cumulativeExpectationAlternate(k,iCase,self%countTrees)=countNode(countNodes)
                      self%cumulativeSumDeltaSquared     (k,iCase,self%countTrees)=sum(deltaSquared)
                      cycle
                   end if
                   expectation         =0.0d0
                   expectationAlternate=0.0d0
                   sumBefore           =0.0d0
                   do i=1,countNodes
                      if (isActive(i) .and. omegaSplit(i) <= self%omegasCumulative(k)) then
                         if (statistic /= 'A') then
                            expectation         =expectation         +countAdded(i)
                            expectationAlternate=expectationAlternate+countAdded(i)
                         end if
                         sumBefore=sumBefore+deltaSquared(i)
                      end if
                      if (countChildren(i) == 0 .or. omega(i) > self%omegasCumulative(k) .or. omegaSplit(i) <= self%omegasCumulative(k)) cycle
                      if (statistic == 'C' .and. .not.isMain(i)) cycle
                      logMassAlive=+logMass(i)                          &
                           &       +(log10(massChildren(i))-logMass(i)) &
                           &       *(self%omegasCumulative(k)-omega(i)) &
                           &       /(omegaSplit           (i)-omega(i))
                      expectation         =+expectation                                                                                                 &
                           &               +self%nodeMean(logMassAlive,self%omegasCumulative(k),modulo(fold  ,countFolds)+1,iCase,emptyTouched,outside)
                      expectationAlternate=+expectationAlternate                                                                                        &
                           &               +self%nodeMean(logMassAlive,self%omegasCumulative(k),modulo(fold+1,countFolds)+1,iCase,emptyTouched,outside)
                   end do
                   self%cumulativeExpectation         (k,iCase,self%countTrees)=expectation
                   self%cumulativeExpectationAlternate(k,iCase,self%countTrees)=expectationAlternate
                   self%cumulativeSumDeltaSquared     (k,iCase,self%countTrees)=sumBefore
                end do
                ! Per-split values.
                if (storeSplits) then
                   call splitsReserve(self%countSplitsStored(iCase)+count(isActive,kind=c_size_t))
                   do i=1,countNodes
                      if (.not.isActive(i)) cycle
                      self%countSplitsStored(iCase)=self%countSplitsStored(iCase)+1_c_size_t
                      self%splitsReal   (splitPropertyRealDelta            %ID,self%countSplitsStored(iCase),iCase)=delta         (                        i         )
                      self%splitsReal   (splitPropertyRealDeltaSquared     %ID,self%countSplitsStored(iCase),iCase)=deltaSquared  (                        i         )
                      self%splitsReal   (splitPropertyRealOmega            %ID,self%countSplitsStored(iCase),iCase)=omegaSplit    (                        i         )
                      self%splitsReal   (splitPropertyRealRedshift         %ID,self%countSplitsStored(iCase),iCase)=redshiftOfNode(nodes     (primaryChild(i))%node  )
                      self%splitsReal   (splitPropertyRealMass             %ID,self%countSplitsStored(iCase),iCase)=massChildren  (                        i         )
                      self%splitsReal   (splitPropertyRealMassPrimary      %ID,self%countSplitsStored(iCase),iCase)=mass          (           primaryChild(i)        )
                      self%splitsReal   (splitPropertyRealMassSecondary    %ID,self%countSplitsStored(iCase),iCase)=massSecondary (                        i         )
                      self%splitsReal   (splitPropertyRealMassRoot         %ID,self%countSplitsStored(iCase),iCase)=mass          (countNodes)
                      self%splitsReal   (splitPropertyRealN1Parent         %ID,self%countSplitsStored(iCase),iCase)=nodeMeanEnd   (                        i       ,1)
                      if (statistic == 'A') then
                         self%splitsReal(splitPropertyRealCountAdded       %ID,self%countSplitsStored(iCase),iCase)=0.0d0
                      else
                         self%splitsReal(splitPropertyRealCountAdded       %ID,self%countSplitsStored(iCase),iCase)=countAdded    (                        i         )
                      end if
                      self%splitsInteger(splitPropertyIntegerTreeIndex     %ID,self%countSplitsStored(iCase),iCase)=treeCurrent%index
                      self%splitsInteger(splitPropertyIntegerIsMainBranch  %ID,self%countSplitsStored(iCase),iCase)=merge(1_kind_int8,0_kind_int8,isMain(i))
                      self%splitsInteger(splitPropertyIntegerHierarchyDepth%ID,self%countSplitsStored(iCase),iCase)=hierarchyDepth(i)
                   end do
                end if
             end select
          end do
       end do
       deallocate(nodes         )
       deallocate(parent        )
       deallocate(primaryChild  )
       deallocate(countChildren )
       deallocate(mass          )
       deallocate(logMass       )
       deallocate(omega         )
       deallocate(massChildren  )
       deallocate(massSecondary )
       deallocate(omegaSplit    )
       deallocate(omegaEnd      )
       deallocate(isPrimary     )
       deallocate(isMain        )
       deallocate(isActive      )
       deallocate(hierarchyDepth)
       deallocate(countAdded    )
       deallocate(countOwn      )
       deallocate(countNode     )
       deallocate(crossesTarget )
       if (self%mode == branchingVarianceModeJumps) then
          deallocate(nodeMeanNode)
          deallocate(nodeMeanEnd )
          deallocate(delta       )
          deallocate(deltaSquared)
          deallocate(drift       )
          deallocate(deltaVariant)
          deallocate(binnedTree  )
       end if
       treeCurrent => treeCurrent%nextTree
    end do
    return

  contains

    double precision function redshiftOfNode(node)
      !!{RST
      Return the redshift of the given node.
      !!}
      implicit none
      type (treeNode          ), intent(inout) :: node
      class(nodeComponentBasic), pointer       :: basic_

      basic_         => node%basic()
      redshiftOfNode =  self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(basic_%time()))
      return
    end function redshiftOfNode

    subroutine splitsReserve(countRequired)
      !!{RST
      Ensure that the per-split storage can hold at least ``countRequired`` splits for every case, growing it geometrically.
      !!}
      implicit none
      integer         (c_size_t ), intent(in   )                 :: countRequired
      double precision           , allocatable, dimension(:,:,:) :: splitsRealTmp
      integer         (kind_int8), allocatable, dimension(:,:,:) :: splitsIntegerTmp
      integer         (c_size_t )                                :: countAllocated

      if (.not.allocated(self%splitsReal)) then
         countAllocated=max(countRequired,treeCountIncrement)
         allocate(self%splitsReal   (splitPropertyRealCount   ,countAllocated,self%countCases))
         allocate(self%splitsInteger(splitPropertyIntegerCount,countAllocated,self%countCases))
      else if (countRequired > size(self%splitsReal,dim=2,kind=c_size_t)) then
         countAllocated=max(countRequired,2_c_size_t*size(self%splitsReal,dim=2,kind=c_size_t))
         call move_alloc(self%splitsReal   ,splitsRealTmp   )
         call move_alloc(self%splitsInteger,splitsIntegerTmp)
         allocate(self%splitsReal   (splitPropertyRealCount   ,countAllocated,self%countCases))
         allocate(self%splitsInteger(splitPropertyIntegerCount,countAllocated,self%countCases))
         self%splitsReal   (:,1:size(splitsRealTmp   ,dim=2),:)=splitsRealTmp
         self%splitsInteger(:,1:size(splitsIntegerTmp,dim=2),:)=splitsIntegerTmp
         deallocate(splitsRealTmp   )
         deallocate(splitsIntegerTmp)
      end if
      return
    end subroutine splitsReserve

  end subroutine branchingVarianceOperatePreEvolution

  subroutine branchingVarianceFinalize(self)
    !!{RST
    Write the accumulated table sums (``table`` mode), or the per-tree quantities and checks (``jumps`` mode), to the output
    file, adding to any already present. The table sums are shared between all copies of the operator, so are written only by
    the last copy to be finalized (which then frees the tables). This assumes that each copy is finalized once, after it has
    processed all of its trees, as by the ``evolveForests`` task.
    !!}
    use    :: Error             , only : Error_Report
    use    :: HDF5              , only : hsize_t
    use    :: HDF5_Access       , only : hdf5Access
    use    :: IO_HDF5           , only : hdf5Group
    use    :: Output_HDF5       , only : outputFile
    use    :: ISO_Varying_String, only : char      , extract
    implicit none
    class           (mergerTreeOperatorBranchingVariance), intent(inout)                       :: self
    integer         (hsize_t                            ), parameter                           :: chunkSize  =1024_hsize_t
    double precision                                     , allocatable  , dimension(:,:,:,:  ) :: sums
    double precision                                     , allocatable  , dimension(:,:,:    ) :: binnedRead
    double precision                                     , allocatable  , dimension(:,:      ) :: markovSums
    double precision                                     , allocatable  , dimension(:        ) :: redshifts               , omegas    , &
         &                                                                                        edges
    integer         (kind_int8                          )                                      :: countRead
    integer                                                                                    :: iStatistic              , iThreshold, &
         &                                                                                        iCase                   , i         , &
         &                                                                                        j
    logical                                                                                    :: preexisting             , lastCopy
    integer                                                                                    :: countFinalized

    ! Determine whether this is the last copy of the operator to be finalized. If the operator was never copied, it is the only
    ! one.
    !$omp critical(branchingVarianceFinalize)
    self%tables%countFinalized=self%tables%countFinalized+1
    lastCopy                  =self%tables%countFinalized >= self%tables%countCopies
    countFinalized            =self%tables%countFinalized
    !$omp end critical(branchingVarianceFinalize)
    ! The tables are freed when the last copy is finalized, so further finalizations (or use) are not possible.
    if (countFinalized > max(self%tables%countCopies,1)) call Error_Report('operator finalized more times than it was copied - its tables have already been freed'//{introspection:location})
    ! Tabulate ω as a function of redshift, outside of the HDF5 lock (evaluating ω may itself require HDF5 access).
    allocate(redshifts(countOmegaTabulation))
    allocate(omegas   (countOmegaTabulation))
    do i=1,countOmegaTabulation
       redshifts(i)=self%redshiftTableMinimum+(self%redshiftTableMaximum-self%redshiftTableMinimum)*dble(i-1)/dble(countOmegaTabulation-1)
       omegas   (i)=self%omega(self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshifts(i))),self%massTableMaximum)
    end do
    ! HDF5 objects are scoped within a block so that they are closed (by finalization) while the HDF5 lock is held.
    !$ call hdf5Access%set()
    block
       type(hdf5Group) :: outputGroup, caseGroup  , treesGroup, statisticGroup, &
            &             binnedGroup, splitsGroup
       outputGroup=outputFile%openGroup(char(self%outputGroupName),'Decomposition of the variance of counts in merger trees into contributions from splits.',objectsOverwritable=.true.,overwriteOverride=.true.)
       preexisting=outputGroup%hasAttribute('statistics')
       if (.not.preexisting) then
          ! Record the statistics, thresholds, and binning, and the relation between redshift and ω.
          call                                              outputGroup%writeAttribute(char(self%statistics                                                        ),'statistics'          )
          call                                              outputGroup%writeAttribute(     self%massThresholds                                                     ,'massThresholds'      )
          call                                              outputGroup%writeAttribute(     self%massTableMinimum                                                   ,'massTableMinimum'    )
          call                                              outputGroup%writeAttribute(     self%massTableMaximum                                                   ,'massTableMaximum'    )
          call                                              outputGroup%writeAttribute(     self%redshiftTableMinimum                                               ,'redshiftTableMinimum')
          call                                              outputGroup%writeAttribute(     self%redshiftTableMaximum                                               ,'redshiftTableMaximum')
          call                                              outputGroup%writeAttribute(     self%binWidthMass                                                       ,'binWidthMass'        )
          call                                              outputGroup%writeAttribute(     self%binWidthOmega                                                      ,'binWidthOmega'       )
          call                                              outputGroup%writeAttribute(     self%redshiftTarget                                                     ,'redshiftTarget'      )
          if (any(self%caseIsTarget)) then
             call                                           outputGroup%writeAttribute(     self%omegaTarget                                                        ,'omegaTarget'         )
             call                                           outputGroup%writeAttribute(     self%countBinsOmegaTarget                                               ,'countBinsOmegaTarget')
          end if
          call                                              outputGroup%writeAttribute(          countFolds                                                         ,'countFolds'          )
          call                                              outputGroup%writeAttribute(     self%logMassEdgeMinimum                                                 ,'logMassEdgeMinimum'  )
          call                                              outputGroup%writeAttribute(     self%omegaEdgeMinimum                                                   ,'omegaEdgeMinimum'    )
          call                                              outputGroup%writeAttribute(     self%countBinsMass                                                      ,'countBinsMass'       )
          call                                              outputGroup%writeAttribute(     self%countBinsOmega                                                     ,'countBinsOmega'      )
          call                                              outputGroup%writeAttribute(char(enumerationBranchingVarianceModeDecode(self%mode,includePrefix=.false.)),'mode'                )
          if (self%mode == branchingVarianceModeJumps) call outputGroup%writeAttribute(char(self%tableFileName                                                     ),'tableFileName'       )
          call outputGroup%writeDataset  (redshifts  ,'redshift'            ,'Redshifts at which ω is tabulated.'                           )
          call outputGroup%writeDataset  (omegas     ,'omega'               ,'ω at the tabulated redshifts (evaluated at massTableMaximum).')
          call outputGroup%writeAttribute(0_kind_int8,'countSplits'                                                                         )
          call outputGroup%writeAttribute(0_kind_int8,'countMultiple'                                                                       )
          call outputGroup%writeAttribute(0_kind_int8,'countPrimaryMismatch'                                                                )
          call outputGroup%writeAttribute(0_kind_int8,'countNodesOutside'                                                                   )
          call outputGroup%writeAttribute(0.0d0      ,'omegaSpreadMaximum'                                                                  )
          if (self%mode == branchingVarianceModeJumps) then
             ! Record the bins of the binned sums of Δ², the cumulative grid, and the tree stride for per-split values.
             binnedGroup=outputGroup%openGroup('binned','Bin edges of the binned sums of Δ².')
             allocate(edges(self%countBinsOmegaBinned+1))
             do i=1,self%countBinsOmegaBinned+1
                edges(i)=self%omegaEdges(1)+self%binWidthOmegaBinned*dble(i-1)
             end do
             call binnedGroup%writeDataset(edges,'omega','Edges of the ω bins.')
             deallocate(edges)
             do j=1,binnedQuantityCount
                allocate(edges(self%countBinsBinned(j)+1))
                do i=1,self%countBinsBinned(j)+1
                   edges(i)=self%logMassBinnedMinimum(j)+dble(i-1)/dble(self%countBinsPerDecadeBinned)
                end do
                call binnedGroup%writeDataset(edges,char(enumerationBinnedQuantityDecode(j,includePrefix=.false.)),'Edges of the bins in the logarithm (base 10) of this mass quantity.')
                deallocate(edges)
             end do
             call outputGroup%writeDataset  (self%redshiftsCumulative,'redshiftCumulative','Redshifts of the grid of the cumulative variance.')
             call outputGroup%writeDataset  (self%omegasCumulative   ,'omegaCumulative'   ,'ω of the grid of the cumulative variance.'        )
             call outputGroup%writeAttribute(self%treeStrideSplits   ,'treeStrideSplits'                                                      )
          end if
       end if
       ! Accumulate counters.
       call outputGroup%readAttribute ('countSplits'         ,countRead)
       call outputGroup%writeAttribute(countRead+self%countSplits         ,'countSplits'         )
       call outputGroup%readAttribute ('countMultiple'       ,countRead)
       call outputGroup%writeAttribute(countRead+self%countMultiple       ,'countMultiple'       )
       call outputGroup%readAttribute ('countPrimaryMismatch',countRead)
       call outputGroup%writeAttribute(countRead+self%countPrimaryMismatch,'countPrimaryMismatch')
       call outputGroup%readAttribute ('countNodesOutside'   ,countRead)
       call outputGroup%writeAttribute(countRead+self%countNodesOutside   ,'countNodesOutside'   )
       block
         double precision :: omegaSpreadMaximum
         call outputGroup%readAttribute ('omegaSpreadMaximum',omegaSpreadMaximum)
         call outputGroup%writeAttribute(max(omegaSpreadMaximum,self%omegaSpreadMaximum),'omegaSpreadMaximum')
       end block
       ! Per-tree quantities common to all statistics.
       if (self%mode == branchingVarianceModeJumps .and. self%countTrees > 0) then
          treesGroup=outputGroup%openGroup('trees','Per-tree quantities.')
          call treesGroup%writeDataset(self%treeIndex        (1:self%countTrees),'treeIndex'        ,'Index of the tree.'                                                     ,appendTo=.true.,chunkSize=chunkSize)
          call treesGroup%writeDataset(self%treeSeed         (1:self%countTrees),'treeSeed'         ,'Random number seed of the tree.'                                        ,appendTo=.true.,chunkSize=chunkSize)
          call treesGroup%writeDataset(self%treeWeight       (1:self%countTrees),'treeWeight'       ,'Volume weight of the tree.'                                             ,appendTo=.true.,chunkSize=chunkSize)
          call treesGroup%writeDataset(self%massRoot         (1:self%countTrees),'massRoot'         ,'Mass of the root node [M☉].'                                            ,appendTo=.true.,chunkSize=chunkSize)
          call treesGroup%writeDataset(self%omegaFormation   (1:self%countTrees),'omegaFormation'   ,'ω at which the main branch first falls below half the root mass.'       ,appendTo=.true.,chunkSize=chunkSize)
          call treesGroup%writeDataset(self%redshiftFormation(1:self%countTrees),'redshiftFormation','Redshift at which the main branch first falls below half the root mass.',appendTo=.true.,chunkSize=chunkSize)
       end if
       ! Per-case quantities.
       do iStatistic=1,self%countStatistics
          do iThreshold=1,self%countThresholds
             iCase    =self%caseIndex(iStatistic,iThreshold)
             statisticGroup=outputGroup   %openGroup('statistic'//char(extract(self%statistics,iStatistic,iStatistic)),'Results for one statistic.'                   ,objectsOverwritable=.true.,overwriteOverride=.true.)
             caseGroup     =statisticGroup%openGroup(char(thresholdGroupName(iThreshold))                             ,'Results for one statistic and threshold mass.',objectsOverwritable=.true.,overwriteOverride=.true.)
             if (.not.caseGroup%hasAttribute('massThreshold')) call caseGroup%writeAttribute(self%massThresholds(iThreshold),'massThreshold')
             select case (self%mode%ID)
             case (branchingVarianceModeTable%ID)
                if (lastCopy) then
                   if (caseGroup%hasDataset('tableSums')) then
                      call caseGroup%readDataset('tableSums',sums)
                      self%tables%sums(:,:,:,:,iCase)=self%tables%sums(:,:,:,:,iCase)+sums
                      deallocate(sums)
                   end if
                   call caseGroup%writeDataset(self%tables%sums(:,:,:,:,iCase),'tableSums','Sums for the n₁ table, indexed by (sum, ω bin, mass bin, fold).')
                end if
             case (branchingVarianceModeJumps%ID)
                if (caseGroup%hasDataset('markovSums')) then
                   call caseGroup%readDataset('markovSums',markovSums)
                   self%markovSums(:,:,iCase)=self%markovSums(:,:,iCase)+markovSums
                   deallocate(markovSums)
                   call caseGroup%readAttribute('evaluationsEmpty'  ,countRead)
                   self%evaluationsEmpty  (iCase)=self%evaluationsEmpty  (iCase)+countRead
                   call caseGroup%readAttribute('evaluationsOutside',countRead)
                   self%evaluationsOutside(iCase)=self%evaluationsOutside(iCase)+countRead
                end if
                call caseGroup%writeDataset  (self%markovSums        (:,:,iCase),'markovSums'        ,'Markov check sums (count, sum and sum of squares of count minus n₁) for roots, primary and non-primary progenitors at splits, and accretion nodes.')
                call caseGroup%writeAttribute(self%evaluationsEmpty  (    iCase),'evaluationsEmpty'                                                                                                                                                       )
                call caseGroup%writeAttribute(self%evaluationsOutside(    iCase),'evaluationsOutside'                                                                                                                                                     )
                if (self%countTrees > 0) then
                   do i=1,treeQuantityCount
                      call caseGroup%writeDataset(self%treeQuantities(i,iCase,1:self%countTrees),char(enumerationTreeQuantityDecode(i,includePrefix=.false.)),char(enumerationTreeQuantityDescription(i)),appendTo=.true.,chunkSize=chunkSize)
                   end do
                   if (self%countCumulative > 0) then
                      call caseGroup%writeDataset(self%cumulativeExpectation         (:,iCase,1:self%countTrees),'cumulativeExpectation'         ,'Expectation of N conditioned on the tree built to each ω of the cumulative grid, indexed by (redshift, tree).'                   ,appendTo=.true.,appendDimension=2)
                      call caseGroup%writeDataset(self%cumulativeExpectationAlternate(:,iCase,1:self%countTrees),'cumulativeExpectationAlternate','As cumulativeExpectation, but evaluated with the n₁ table of the second other fold.'                                             ,appendTo=.true.,appendDimension=2)
                      call caseGroup%writeDataset(self%cumulativeSumDeltaSquared     (:,iCase,1:self%countTrees),'cumulativeSumDeltaSquared'     ,'Sum of Δ² over splits at ω no greater than (i.e. at redshifts up to) each ω of the cumulative grid, indexed by (redshift, tree).',appendTo=.true.,appendDimension=2)
                   end if
                end if
                ! Binned sums of Δ², added to any already present.
                binnedGroup=caseGroup%openGroup('binned','Sums of Δ² in bins.',objectsOverwritable=.true.,overwriteOverride=.true.)
                do j=1,binnedQuantityCount
                   if (binnedGroup%hasDataset(char(enumerationBinnedQuantityDecode(j,includePrefix=.false.)))) then
                      call binnedGroup%readDataset(char(enumerationBinnedQuantityDecode(j,includePrefix=.false.))          ,binnedRead)
                      self%binnedSums       (:,1:self%countBinsBinned(j),:,j,iCase)=self%binnedSums       (:,1:self%countBinsBinned(j),:,j,iCase)+binnedRead
                      deallocate(binnedRead)
                      call binnedGroup%readDataset(char(enumerationBinnedQuantityDecode(j,includePrefix=.false.))//'Squared',binnedRead)
                      self%binnedSumsSquared(:,1:self%countBinsBinned(j),:,j,iCase)=self%binnedSumsSquared(:,1:self%countBinsBinned(j),:,j,iCase)+binnedRead
                      deallocate(binnedRead)
                   end if
                   call binnedGroup%writeDataset(self%binnedSums       (:,1:self%countBinsBinned(j),:,j,iCase),char(enumerationBinnedQuantityDecode(j,includePrefix=.false.))           ,char(enumerationBinnedQuantityDescription(j))                                                             )
                   call binnedGroup%writeDataset(self%binnedSumsSquared(:,1:self%countBinsBinned(j),:,j,iCase),char(enumerationBinnedQuantityDecode(j,includePrefix=.false.))//'Squared','Sums over trees of the square of the weighted sum of each tree in each bin of the corresponding dataset.')
                end do
                ! Per-split values.
                if (self%countSplitsStored(iCase) > 0) then
                   splitsGroup=caseGroup%openGroup('splits','Per-split values.')
                   do i=1,splitPropertyRealCount
                      call splitsGroup%writeDataset(self%splitsReal   (i,1:self%countSplitsStored(iCase),iCase),char(enumerationSplitPropertyRealDecode   (i,includePrefix=.false.)),char(enumerationSplitPropertyRealDescription   (i)),appendTo=.true.,chunkSize=chunkSize)
                   end do
                   do i=1,splitPropertyIntegerCount
                      call splitsGroup%writeDataset(self%splitsInteger(i,1:self%countSplitsStored(iCase),iCase),char(enumerationSplitPropertyIntegerDecode(i,includePrefix=.false.)),char(enumerationSplitPropertyIntegerDescription(i)),appendTo=.true.,chunkSize=chunkSize)
                   end do
                end if
             end select
          end do
       end do
    end block
    call outputFile%flush()
    !$ call hdf5Access%unset()
    ! Free the tables once no copy needs them. The (small) object holding them is not freed, as other copies still point to it.
    if (lastCopy) then
       if (allocated(self%tables%sums    )) deallocate(self%tables%sums    )
       if (allocated(self%tables%nodeMean)) deallocate(self%tables%nodeMean)
       if (allocated(self%tables%empty   )) deallocate(self%tables%empty   )
    end if
    return
  end subroutine branchingVarianceFinalize
