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
  Implements a task which solves the equations for the mean and variance of subhalo counts in merger trees built by the
  :cite:t:`cole_hierarchical_2000` algorithm, without building any trees.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass      , criticalOverdensityClass
  use :: Merger_Tree_Branching     , only : mergerTreeBranchingProbabilityClass

  !![
  <task name="taskMergerTreeBranchingMoments" docformat="rst">
   <description>
   A task which computes the mean, :math:`n_1`, and the variance, :math:`W`, of the count of subhalos in a merger tree built
   by the :cite:t:`cole_hierarchical_2000` algorithm, by solving the equations that they obey, rather than by building
   trees. The trees are those of the ``cole2000`` merger tree builder with the given
   ``[mergerTreeBranchingProbability]``: in each interval :math:`\mathrm{d}\omega` a branch of mass :math:`M` at
   :math:`\omega` (the variable :math:`\omega=\delta_\mathrm{c}(t,M)\sigma(M,t_0)/\sigma(M,t)` of the builder) splits into
   progenitors of mass :math:`M-M^\prime` (the primary progenitor) and :math:`M^\prime\le M/2` with probability
   :math:`f(M^\prime|M,\omega)\mathrm{d}M^\prime\mathrm{d}\omega` for :math:`M^\prime\ge M_\mathrm{res}`, and otherwise loses
   mass to subresolution halos at the rate :math:`\mathrm{d}\ln M/\mathrm{d}\omega=-\phi(M,\omega)`. Here :math:`f` and
   :math:`\phi` are the ``rate`` and ``fractionSubresolution`` of the branching probability class, each multiplied by
   :math:`\sigma(M,t)/\sigma(M,t_0)`, as in the builder. The builder applies accretion in steps of up to its
   ``[accretionLimit]``; this task solves the limit of continuous accretion.

   The statistics, :math:`N`, are:

   ``A``
     the number of progenitors of mass at least :math:`m` at the redshift ``[redshiftTarget]``, :math:`z_1` (the progenitor
     mass function);
   ``B``
     the number of splits anywhere in the tree producing a non-primary progenitor of mass at least :math:`m` (the unevolved
     subhalo mass function summed over all hierarchy depths);
   ``C``
     as ``B``, but counting only splits of the main branch (the chain of primary progenitors from the root).

   These are the statistics of the ``branchingVariance`` merger tree operator. For a node of mass :math:`M` at
   :math:`\omega`, let :math:`n_1(M,\omega)` be the mean of the count from the node's subtree and :math:`W(M,\omega)` its
   variance. A split into :math:`M-M^\prime` and :math:`M^\prime` changes the expected count by the jump
   :math:`\Delta=c+n_1(M-M^\prime)+n_1(M^\prime)-n_1(M)` (for ``C``, :math:`\Delta=c+n_1(M-M^\prime)-n_1(M)`, since only
   the primary progenitor continues the main branch), where, for ``B`` and ``C``, :math:`c=1` if :math:`M^\prime\ge m` and
   :math:`0` otherwise, and, for ``A``, :math:`c=0`.
   The expected count from the subtree is a martingale, so along each branch between splits

   .. math::

      \frac{\mathrm{d}n_1}{\mathrm{d}\omega} = -\int f\Delta\,\mathrm{d}M^\prime,\quad
      \frac{\mathrm{d}W}{\mathrm{d}\omega} = -\int f\left[\Delta^2+W(M-M^\prime)+W(M^\prime)-W(M)\right]\mathrm{d}M^\prime,

   with :math:`\mathrm{d}M/\mathrm{d}\omega=-\phi M`, and the :math:`W(M^\prime)` term absent for ``C``. For ``B`` and
   ``C`` both vanish at ``[redshiftMaximum]`` and for :math:`M\le M_\mathrm{res}`. For ``A``, :math:`n_1=\theta(M-m)` and
   :math:`W=0` at :math:`\omega_1`, the :math:`\omega` of :math:`z_1`. The equations are integrated, from there to the root,
   on a grid uniform in :math:`\ln M` (from :math:`M_\mathrm{res}` to the root mass) and in :math:`\omega` (with
   :math:`\omega_1` on the grid if ``A`` is computed, so that the steps on either side of it differ slightly), along the
   characteristics :math:`\mathrm{d}M/\mathrm{d}\omega=-\phi M` (the paths of branches between splits). Each step, from
   :math:`\omega_{k+1}` to :math:`\omega_k`, follows the characteristic through each point :math:`(M_i,\omega_k)` of the grid
   back to its *foot*: the point at which it crosses the level :math:`\omega_{k+1}` already solved, at the mass
   :math:`M_\mathrm{foot}\le M_i` that a branch of mass :math:`M_i` at :math:`\omega_k` had at :math:`\omega_{k+1}`. The values
   of :math:`n_1`, :math:`W` and their sources at the foot are found by cubic interpolation in :math:`\ln M`, and carried to
   :math:`M_i` by a second-order (Heun) step, which is also used to locate the foot. The integrals over :math:`M^\prime\in[M_\mathrm{res},M/2]` use the trapezoid rule in :math:`\ln
   M^\prime` on the grid, plus the points :math:`M/2` and :math:`m`. The task requires :math:`\omega` to be independent of
   mass, as it is for a mass-independent critical overdensity and scale-independent growth.

   For ``A``, since the progenitors of a branch are less massive than it, :math:`n_1` and :math:`W` vanish below
   :math:`M_\mathrm{c}(\omega)`, the mass on the characteristic that reaches :math:`m` at :math:`\omega_1`, and jump at it to
   :math:`P` and :math:`P(1-P)`, where :math:`P` is the probability of no split along that characteristic. Above it lies a
   layer, of width about :math:`M_\mathrm{res}`, in which no split leaves a progenitor in the support, so that the slope
   of :math:`n_1` changes sharply across it, and :math:`n_1` rises steeply above it where :math:`m\gg M_\mathrm{res}`.
   This structure can be much narrower than the step of the grid, so it is resolved by additional nodes, which follow their
   characteristics from :math:`m` and from offsets above it at :math:`\omega_1` (the first is the boundary itself). They are
   spaced geometrically from :math:`M_\mathrm{res}/4`, by at most half of a step of the grid, up to at least
   :math:`4M_\mathrm{res}` and five steps of the grid above :math:`m`. Within the support, interpolation near the boundary
   uses these nodes and then the points of the grid above them. The integrals over :math:`M^\prime` for ``A`` add to the
   quadrature nodes the images of these nodes, :math:`M^\prime=M_j` and :math:`M-M^\prime=M_j` (with the kernel
   interpolated between the nodes of the grid), and take the steps of the integrand, at :math:`M^\prime=M_\mathrm{c}`
   and :math:`M-M^\prime=M_\mathrm{c}`, as nodes with the values on either side.

   Written to the group ``[outputGroupName]`` of the output file are the grids (datasets ``mass``, ``omega`` and
   ``redshift``), :math:`\sigma(M,t_0)` on the mass grid (``rootVariance``), the rates :math:`\phi`
   (``fractionSubresolution``), :math:`R=\int f\mathrm{d}M^\prime` from the ``probability`` method of the branching
   probability class (``rateSplit``) and from the quadrature used here (``rateSplitQuadrature``), and, in a subgroup
   ``statistic&lt;S&gt;/threshold&lt;i&gt;`` for each statistic and threshold, :math:`n_1` (``countMean``), :math:`W`
   (``countVariance``), and the source of the variance, :math:`\int f\Delta^2\mathrm{d}M^\prime`
   (``sourceDeltaSquared``), on the grid, with the values at the root as attributes (``countMeanRoot``,
   ``countVarianceRoot``). For ``A`` the subgroup also holds :math:`M_\mathrm{c}` (``massBoundary``, zero beyond
   :math:`\omega_1` and where :math:`M_\mathrm{c}` exceeds the root mass), and all of its tables are zero beyond
   :math:`\omega_1`. Two-dimensional datasets are indexed by :math:`\omega` and then mass (in the C order of HDF5).

   If ``[decomposeVariance]`` is true, the task also finds the densities of branches, forward in :math:`\omega` from the root, and
   decomposes the variance into the contributions of splits, as the ``branchingVariance`` merger tree operator measures them in
   trees. The density of branches, :math:`\rho(M,\omega)`, is the expected number of branches per unit :math:`\ln M` alive at
   :math:`\omega`, and the density of the main branch, :math:`\rho_\mathrm{main}(M,\omega)`, is the probability density of its
   mass. They are found by the transpose of the scheme which solves the equations: each step of that scheme is a linear map from
   the solution at :math:`\omega_{k+1}` to that at :math:`\omega_k` (plus the sources), and its transpose maps the numbers of
   branches at the points of the grid at :math:`\omega_k` to those at :math:`\omega_{k+1}` (starting from a single branch at the
   root). For :math:`\rho_\mathrm{main}` only the primary progenitor of each split is followed. The transpose of the interpolation
   to the feet of the characteristics moves each branch along its characteristic and shares it between the four nearest points of
   the grid, so the density at a single point of the grid ripples (typically by much less than a percent, but by more where the
   density varies rapidly, and the narrow peak of the branch from the root which has not yet split rings, with negative values),
   while its integrals against smooth functions are accurate. Because the densities are the exact transpose of the scheme, for
   statistics ``B`` and ``C`` (with :math:`\rho_\mathrm{main}`, as only the main branch counts)

   .. math::

      n_1(M_0,\omega_0) = \langle\textstyle\sum c\rangle_{&lt;\omega} + \int\rho\,n_1\mathrm{d}\ln M,\quad
      W(M_0,\omega_0) = V(\omega) + \int\rho\,W\mathrm{d}\ln M,

   hold to rounding error at every :math:`\omega` of the grid, where :math:`\langle\sum c\rangle_{&lt;\omega}` is the count
   expected from the splits before :math:`\omega`, and :math:`V(\omega)=\mathrm{Var}(\mathrm{E}[N|\omega])` is the
   cumulative variance, the variance contributed by the splits before :math:`\omega`. With the trapezoid rule in
   :math:`\omega` implied by the scheme, :math:`V(\omega)=\int^\omega\mathrm{d}\omega^\prime\int\rho\,s\,\mathrm{d}\ln M`,
   where :math:`s=\int f\Delta^2\mathrm{d}M^\prime` is the source of the variance. For statistic ``A``, which counts all
   branches but is not solved by the same scheme, :math:`\rho` is integrated against its solution (near the boundary of its
   support, by the trapezoid rule over the nodes there), and the identities hold to the accuracy of the solution. The largest
   departures from them over :math:`\omega` are written to the attributes ``identityResidualMean`` and
   ``identityResidualVariance`` of each statistic and threshold.

   The contributions of the splits in each step in :math:`\omega` are then classified, as by the operator, by the mass of the
   parent, :math:`M`, by the mass ratio, :math:`x=M^\prime/M`, and by class: splits of the main branch (weighted by
   :math:`\rho_\mathrm{main}`) or of side branches (weighted by :math:`\rho-\rho_\mathrm{main}`), each either terminal
   (:math:`M^\prime\ge m`) or structural. For statistic ``C`` all splits are of the main branch. Mass ratios are those of a
   lattice, :math:`x=\exp(-d\,\Delta\ln M)` for :math:`d=0,1,\ldots`, where :math:`\Delta\ln M` is the step of the mass
   grid, and the contribution of each quadrature node is assigned to the nearest point of the lattice.

   The densities, :math:`\rho` (``densityBranches``) and :math:`\rho_\mathrm{main}` (``densityMainBranch``), are written on
   the grid as the other tables, and the lattice of mass ratios to ``massRatio``. For each statistic and threshold,
   :math:`V(\omega)` is written to ``varianceCumulative``, and the contributions of splits to ``varianceContributionsMass``
   (by the mass of the parent, on the mass grid) and ``varianceContributionsMassRatio`` (by mass ratio, on the lattice), each
   summed over the steps between consecutive values of :math:`\omega` of the tables, and indexed (in the C order of HDF5) by
   class (main branch terminal, main branch structural, side branch terminal, side branch structural), interval in
   :math:`\omega`, and mass or mass ratio. Each sums to the variance of the count (exactly for statistics ``B`` and ``C``, and,
   for statistic ``A``, to the accuracy of its scheme near the boundary of its support, as do the identities above). Where the
   density rings, next to the root, the contributions include negative values (typically 1 to 2 percent of the variance, balanced by
   larger positive values nearby). For statistic ``A`` with a threshold close to the root mass, whose support then spans only a
   few points of the grid near the root, this ringing makes the decomposition (and :math:`V(\omega)`) inaccurate, which the
   attributes ``identityResidualMean`` and ``identityResidualVariance`` show. The density
   at the lowest point of the grid, :math:`M_\mathrm{res}`, is that of branches which fell below the resolution in the last
   step: they have no count, and are removed in the next step. The densities require the kernel at each
   :math:`\omega` of the grid a second time, so roughly double the cost of the task.
   </description>
  </task>
  !!]
  type, extends(taskClass) :: taskMergerTreeBranchingMoments
     !!{RST
     Implementation of a task which solves the equations for the mean and variance of subhalo counts in merger trees.
     !!}
     private
     class           (cosmologyFunctionsClass            ), pointer                   :: cosmologyFunctions_             => null()
     class           (criticalOverdensityClass           ), pointer                   :: criticalOverdensity_            => null()
     class           (cosmologicalMassVarianceClass      ), pointer                   :: cosmologicalMassVariance_       => null()
     class           (mergerTreeBranchingProbabilityClass), pointer                   :: mergerTreeBranchingProbability_ => null()
     double precision                                                                 :: massRoot                                  , massResolution , &
          &                                                                              redshiftRoot                              , redshiftMaximum, &
          &                                                                              redshiftTarget
     double precision                                     , allocatable, dimension(:) :: massThresholds
     integer                                                                          :: countMassPerDecade                        , countOmega     , &
          &                                                                              outputStrideOmega
     type            (varying_string                     )                            :: statistics                                , outputGroupName
     type            (inputParameters                    )                            :: parameters
     logical                                                                          :: nodeComponentsInitialized       =  .false., decomposeVariance
   contains
     final     ::            mergerTreeBranchingMomentsDestructor
     procedure :: perform => mergerTreeBranchingMomentsPerform
  end type taskMergerTreeBranchingMoments

  interface taskMergerTreeBranchingMoments
     !!{RST
     Constructors for the :galacticus-class:`taskMergerTreeBranchingMoments` task.
     !!}
     module procedure mergerTreeBranchingMomentsConstructorParameters
     module procedure mergerTreeBranchingMomentsConstructorInternal
  end interface taskMergerTreeBranchingMoments

contains

  function mergerTreeBranchingMomentsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`taskMergerTreeBranchingMoments` task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskMergerTreeBranchingMoments     )                              :: self
    type            (inputParameters                    ), intent(inout), target       :: parameters
    class           (cosmologyFunctionsClass            ), pointer                     :: cosmologyFunctions_
    class           (criticalOverdensityClass           ), pointer                     :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), pointer                     :: cosmologicalMassVariance_
    class           (mergerTreeBranchingProbabilityClass), pointer                     :: mergerTreeBranchingProbability_
    type            (inputParameters                    ), pointer                     :: parametersRoot
    double precision                                     , allocatable  , dimension(:) :: massThresholds
    double precision                                                                   :: massRoot                       , massResolution , &
         &                                                                                redshiftRoot                   , redshiftMaximum, &
         &                                                                                redshiftTarget
    integer                                                                            :: countMassPerDecade             , countOmega     , &
         &                                                                                outputStrideOmega
    type            (varying_string                     )                              :: statistics                     , outputGroupName
    logical                                                                            :: decomposeVariance

    ! Ensure the nodes objects are initialized (a node is needed as an argument to the branching probability class).
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       parametersRoot => parameters
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    !![
    <inputParameter docformat="rst">
      <name>massRoot</name>
      <source>parameters</source>
      <defaultValue>1.0d12</defaultValue>
      <description>
      The mass of the root halo, :math:`M_0`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshiftRoot</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>
      The redshift of the root halo.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massResolution</name>
      <source>parameters</source>
      <defaultValue>1.0d9</defaultValue>
      <description>
      The mass resolution of the trees, :math:`M_\mathrm{res}`.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshiftMaximum</name>
      <source>parameters</source>
      <defaultValue>10.0d0</defaultValue>
      <description>
      The redshift at which branches are terminated (the ``[redshiftMaximum]`` of the ``cole2000`` merger tree builder, which
      should be set to the same value for a comparison with trees).
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>statistics</name>
      <source>parameters</source>
      <defaultValue>var_str('BC')</defaultValue>
      <description>
      The count statistics to compute, any of ``A``, ``B`` and ``C`` (see the description of this task).
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshiftTarget</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The redshift, :math:`z_1`, at which progenitors are counted for statistic ``A``.
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
      The threshold masses, :math:`m`, for the count statistics.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>countMassPerDecade</name>
      <source>parameters</source>
      <defaultValue>100</defaultValue>
      <description>
      The number of points per decade of mass in the grid on which the equations are solved.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>countOmega</name>
      <source>parameters</source>
      <defaultValue>2000</defaultValue>
      <description>
      The number of steps in :math:`\omega` from the root to ``[redshiftMaximum]``.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputStrideOmega</name>
      <source>parameters</source>
      <defaultValue>1</defaultValue>
      <description>
      The tables on the grid are written at every ``[outputStrideOmega]``\ :sup:`th` step in :math:`\omega`, which must divide
      ``[countOmega]``.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>decomposeVariance</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, the densities of branches are found, and the variance is decomposed into the contributions of splits (see the
      description of this task). This roughly doubles the cost of the task, and, for each statistic and threshold, adds to the
      memory and output two arrays of :math:`4 N_\omega/N_\mathrm{stride}` by :math:`N_M+1` values, where :math:`N_\omega` is
      ``[countOmega]``, :math:`N_\mathrm{stride}` is ``[outputStrideOmega]``, and :math:`N_M+1` is the number of points of the
      mass grid, so a larger ``[outputStrideOmega]`` reduces it.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>outputGroupName</name>
      <source>parameters</source>
      <defaultValue>var_str('branchingMoments')</defaultValue>
      <description>
      The name of the group in the output file to which results are written.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"             name="cosmologyFunctions_"             source="parameters"/>
    <objectBuilder class="criticalOverdensity"            name="criticalOverdensity_"            source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"       name="cosmologicalMassVariance_"       source="parameters"/>
    <objectBuilder class="mergerTreeBranchingProbability" name="mergerTreeBranchingProbability_" source="parameters"/>
    !!]
    self=taskMergerTreeBranchingMoments(massRoot,redshiftRoot,massResolution,redshiftMaximum,redshiftTarget,statistics,massThresholds,countMassPerDecade,countOmega,outputStrideOmega,decomposeVariance,outputGroupName,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBranchingProbability_,parametersRoot)
    self%nodeComponentsInitialized=.true.
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"            />
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    return
  end function mergerTreeBranchingMomentsConstructorParameters

  function mergerTreeBranchingMomentsConstructorInternal(massRoot,redshiftRoot,massResolution,redshiftMaximum,redshiftTarget,statistics,massThresholds,countMassPerDecade,countOmega,outputStrideOmega,decomposeVariance,outputGroupName,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBranchingProbability_,parameters) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`taskMergerTreeBranchingMoments` task class.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char        , len
    implicit none
    type            (taskMergerTreeBranchingMoments     )                              :: self
    double precision                                     , intent(in   )               :: massRoot                       , massResolution , &
         &                                                                                redshiftRoot                   , redshiftMaximum, &
         &                                                                                redshiftTarget
    double precision                                     , intent(in   ), dimension(:) :: massThresholds
    integer                                              , intent(in   )               :: countMassPerDecade             , countOmega     , &
         &                                                                                outputStrideOmega
    type            (varying_string                     ), intent(in   )               :: statistics                     , outputGroupName
    logical                                              , intent(in   )               :: decomposeVariance
    class           (cosmologyFunctionsClass            ), intent(in   ), target       :: cosmologyFunctions_
    class           (criticalOverdensityClass           ), intent(in   ), target       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), intent(in   ), target       :: cosmologicalMassVariance_
    class           (mergerTreeBranchingProbabilityClass), intent(in   ), target       :: mergerTreeBranchingProbability_
    type            (inputParameters                    ), intent(in   ), target       :: parameters
    integer                                                                            :: i
    character       (len=:                              ), allocatable                 :: statisticsList
    !![
    <constructorAssign variables="massRoot, redshiftRoot, massResolution, redshiftMaximum, redshiftTarget, statistics, massThresholds, countMassPerDecade, countOmega, outputStrideOmega, decomposeVariance, outputGroupName, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *mergerTreeBranchingProbability_"/>
    !!]

    if (len(statistics) < 1) call Error_Report('at least one statistic is required'//{introspection:location})
    statisticsList=char(statistics)
    do i=1,len(statisticsList)
       if (statisticsList(i:i) /= 'A' .and. statisticsList(i:i) /= 'B' .and. statisticsList(i:i) /= 'C') &
            & call Error_Report('statistics must be any of "A", "B" and "C"'//{introspection:location})
       if (index(statisticsList(1:i-1),statisticsList(i:i)) /= 0) &
            & call Error_Report('a statistic is repeated'//{introspection:location})
    end do
    if (massResolution     <= 0.0d0          ) call Error_Report('mass resolution must be positive'               //{introspection:location})
    if (any(massThresholds <= massResolution)) call Error_Report('mass thresholds must exceed the mass resolution'//{introspection:location})
    do i=1,size(massThresholds)-1
       if (any(massThresholds(i+1:) == massThresholds(i))) call Error_Report('a mass threshold is repeated'//{introspection:location})
    end do
    if (massRoot           <= 2.0d0*massResolution) call Error_Report('root mass must exceed twice the mass resolution'//{introspection:location})
    if (redshiftMaximum    <= redshiftRoot        ) call Error_Report('maximum redshift must exceed the root redshift' //{introspection:location})
    if (countMassPerDecade <  1                   ) call Error_Report('at least one mass point per decade is required' //{introspection:location})
    if (countOmega         <  1                   ) call Error_Report('at least one step in ω is required'             //{introspection:location})
    if (outputStrideOmega < 1 .or. mod(countOmega,max(outputStrideOmega,1)) /= 0) call Error_Report('[outputStrideOmega] must be positive and divide [countOmega]'//{introspection:location})
    if (index(statisticsList,'A') /= 0) then
       if (redshiftTarget <= redshiftRoot   ) call Error_Report('target redshift must exceed the root redshift'         //{introspection:location})
       if (redshiftTarget >= redshiftMaximum) call Error_Report('target redshift must be below the maximum redshift'    //{introspection:location})
       if (countOmega     <  2              ) call Error_Report('at least two steps in ω are required for statistic "A"'//{introspection:location})
    end if
    self%parameters=inputParameters(parameters)
    return
  end function mergerTreeBranchingMomentsConstructorInternal

  subroutine mergerTreeBranchingMomentsDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`taskMergerTreeBranchingMoments` task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskMergerTreeBranchingMoments), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"            />
    <objectDestructor name="self%criticalOverdensity_"           />
    <objectDestructor name="self%cosmologicalMassVariance_"      />
    <objectDestructor name="self%mergerTreeBranchingProbability_"/>
    !!]
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine mergerTreeBranchingMomentsDestructor

  subroutine mergerTreeBranchingMomentsPerform(self,status)
    !!{RST
    Solve the equations for the mean and variance of the counts, and write the results to the output file.
    !!}
    use :: Display           , only : displayIndent                    , displayMessage                     , displayUnindent
    use :: Error             , only : Error_Report                     , errorStatusSuccess
    use :: Galacticus_Nodes  , only : nodeComponentBasic               , treeNode
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Group
    use :: ISO_Varying_String, only : char                             , extract                            , len, var_str
    use :: Node_Components   , only : Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize
    use :: Output_HDF5       , only : outputFile
    use :: String_Handling   , only : operator(//)
    implicit none
    class           (taskMergerTreeBranchingMoments     ), intent(inout), target           :: self
    integer                                              , intent(  out), optional         :: status
    double precision                                     , parameter                       :: toleranceMass                   =  1.0d-9, toleranceOmega        =1.0d-6
    ! Placement of the nodes near the boundary of the support of statistic A: the first lies this far above the threshold (in
    ! units of the mass resolution), each spacing is this factor larger than the previous one, up to this fraction of a step of
    ! the grid (so that the nodes are nowhere coarser than the grid), and the last lies at least this many steps of the grid (in
    ! ln M), and this many times the mass resolution, above the threshold (so that the nodes span the layer above the boundary,
    ! whatever the grid). At most this many nodes are placed.
    double precision                                     , parameter                         :: spacingNearFirst                =  0.25d0, growthNear             =1.25d0, &
         &                                                                                      spacingNearGrid                 =  0.50d0, extentNear             =5.00d0, &
         &                                                                                      extentNearResolution            =  4.00d0
    integer                                              , parameter                         :: countNearLimit                  =256
    type            (treeNode                           ), pointer                           :: node
    class           (nodeComponentBasic                 ), pointer                           :: basic
    class           (criticalOverdensityClass           ), pointer    , save                 :: criticalOverdensity_            => null()
    class           (cosmologicalMassVarianceClass      ), pointer    , save                 :: cosmologicalMassVariance_       => null()
    class           (mergerTreeBranchingProbabilityClass), pointer    , save                 :: mergerTreeBranchingProbability_ => null()
    type            (treeNode                           ), pointer    , save                 :: nodeThread                      => null()
    type            (inputParameters                    ), allocatable, save                 :: parameters
    !$omp threadprivate(criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBranchingProbability_,nodeThread,parameters)
    double precision                                     , allocatable  , dimension(:      ) :: logMass                                  , mass                          , &
         &                                                                                      omega                                    , redshift                      , &
         &                                                                                      logMassNode                              , weightNode                    , &
         &                                                                                      kernelCurrent                            , fractionPrevious              , &
         &                                                                                      rateCurrent                              , rateQuadratureCurrent         , &
         &                                                                                      logMassFoot                              , countMeanFoot                 , &
         &                                                                                      countVarianceFoot                        , sourceMeanFoot                , &
         &                                                                                      sourceVarianceFoot                       , countMeanPredicted            , &
         &                                                                                      countVariancePredicted                   , sourceMeanPredicted           , &
         &                                                                                      sourceVariancePredicted                  , sourceSquaredPredicted        , &
         &                                                                                      sourceSquaredCurrent                     , fractionCurrent               , &
         &                                                                                      logMassPrimaryNode                       , countMeanNearPredicted        , &
         &                                                                                      countVarianceNearPredicted               , sourceMeanNearPredicted       , &
         &                                                                                      sourceVarianceNearPredicted              , rootVarianceNow
    double precision                                     , allocatable  , dimension(:,:    ) :: weightsSecondary                         , weightsPrimary                , &
         &                                                                                      countMean                                , countVariance                 , &
         &                                                                                      sourceMean                               , sourceVariance                , &
         &                                                                                      fractionSubresolution                    , rateSplit                     , &
         &                                                                                      rateSplitQuadrature                      , massBoundaryTable             , &
         &                                                                                      offsetNear                               , logMassNear                   , &
         &                                                                                      countMeanNear                            , countVarianceNear             , &
         &                                                                                      sourceMeanNear                           , sourceVarianceNear            , &
         &                                                                                      logMassNearPrevious                      , countMeanNearPrevious         , &
         &                                                                                      countVarianceNearPrevious                , sourceMeanNearPrevious        , &
         &                                                                                      sourceVarianceNearPrevious               , logMassNodeNear               , &
         &                                                                                      kernelNear                               , logMassPrimaryNodeNear
    double precision                                     , allocatable  , dimension(:,:,:  ) :: countMeanTable                           , countVarianceTable            , &
         &                                                                                      sourceSquaredTable
    integer                                              , allocatable  , dimension(:      ) :: nodeStart                                , nodeCount                     , &
         &                                                                                      baseSecondary                            , basePrimary                   , &
         &                                                                                      thresholdNode                            , sideNode                      , &
         &                                                                                      caseThreshold
    integer                                              , allocatable  , dimension(:,:    ) :: countAdded
    double precision                                     , allocatable  , dimension(:      ) :: logMassCandidate
    integer                                              , allocatable  , dimension(:      ) :: thresholdCandidate                       , sideCandidate
    integer                                                                                  :: countCandidates
    logical                                                                                  :: countingOnly                             , hasTarget
    logical                                              , allocatable  , dimension(:      ) :: terminatedPrevious                       , terminatedCurrent
    character       (len=1                              ), allocatable  , dimension(:      ) :: caseStatistic
    double precision                                                                         :: timeNow                                  , timeRoot                      , &
         &                                                                                      timeEarliest                             , omegaRoot                     , &
         &                                                                                      omegaEnd                                 , stepLogMass                   , &
         &                                                                                      stepOmega                                , logMassHalf                   , &
         &                                                                                      rateStepMaximum                          , quadratureErrorMaximum        , &
         &                                                                                      timeTarget                               , omegaTarget                   , &
         &                                                                                      stepOmegaBelow                           , stepOmegaAbove                , &
         &                                                                                      logMassFootTarget                        , offset                        , &
         &                                                                                      spacing
    double precision                                                    , dimension(4      ) :: weightsFoot
    integer                                                                                  :: countMass                                , countCases                    , &
         &                                                                                      countNodes                               , countThresholds               , &
         &                                                                                      i                                        , j                             , &
         &                                                                                      k                                        , iCase                         , &
         &                                                                                      iStatistic                               , iThreshold                    , &
         &                                                                                      kTarget                                  , indexFirstPrevious            , &
         &                                                                                      baseFoot                                 , countNearMaximum              , &
         &                                                                                      countNearActivePrevious
    integer                                              , allocatable  , dimension(:      ) :: countNear                                , countNearActive               , &
         &                                                                                      countNodesNear
    type            (hdf5Group                          )                                    :: outputGroup                              , caseGroup                     , &
         &                                                                                      statisticGroup
    type            (varying_string                     )                                    :: message                                  , groupName
    character       (len=16                             )                                    :: label
    ! Densities of branches, and the decomposition of the variance into the contributions of splits.
    double precision                                     , allocatable  , dimension(:      ) :: kernelNext                               , fractionNext                 , &
         &                                                                                      rateNext                                 , rateQuadratureNext           , &
         &                                                                                      weightsPopulation                        , weightsMain                  , &
         &                                                                                      splitsPopulation                         , splitsMain                   , &
         &                                                                                      rateNode                                 , entryWeight                  , &
         &                                                                                      sourceSquaredTerminal                    , sourceAdded                  , &
         &                                                                                      massRatio
    double precision                                     , allocatable  , dimension(:,:    ) :: densityPopulation                        , densityPopulationPredicted   , &
         &                                                                                      densityMain                              , densityMainPredicted         , &
         &                                                                                      stepVariance                             , stepAdded                    , &
         &                                                                                      aliveMean                                , aliveVariance                , &
         &                                                                                      varianceCumulative                       , countAddedCumulative
    double precision                                     , allocatable  , dimension(:,:,:  ) :: squaredRatio
    double precision                                     , allocatable  , dimension(:,:,:,:) :: contributionsMass                        , contributionsRatio
    logical                                              , allocatable  , dimension(:      ) :: terminatedNext                           , entrySecondary               , &
         &                                                                                      ratioActive
    integer                                              , allocatable  , dimension(:      ) :: entryStart                               , entryNode                    , &
         &                                                                                      entryFill                                , latticeNode
    double precision                                     , allocatable  , dimension(:      ) :: identityResidualMean                     , identityResidualVariance
    integer                                                                                  :: countEntries                             , countBlocks                  , &
         &                                                                                      p                                        , q                            , &
         &                                                                                      l                                        , kEnd

    call displayIndent('Begin task: merger tree branching moments')
    call Node_Components_Thread_Initialize(self%parameters)
    node  => treeNode      (                 )
    basic => node    %basic(autoCreate=.true.)
    ! Cases.
    countThresholds=size(self%massThresholds)
    countCases     =len(self%statistics)*countThresholds
    allocate(caseStatistic(countCases))
    allocate(caseThreshold(countCases))
    iCase=0
    do iStatistic=1,len(self%statistics)
       do iThreshold=1,countThresholds
          iCase               =iCase+1
          caseStatistic(iCase)=extract(self%statistics,iStatistic,iStatistic)
          caseThreshold(iCase)=iThreshold
       end do
    end do
    hasTarget=any(caseStatistic == 'A')
    ! Epochs, and the ω of the builder at the root and at the maximum redshift. The ω grid assumes that ω is independent of
    ! mass, so check that.
    timeNow     =self%cosmologyFunctions_%cosmicTime(1.0d0)
    timeRoot    =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftRoot   ))
    timeEarliest=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftMaximum))
    omegaRoot   =omegaOf(timeRoot    ,self%massRoot)
    omegaEnd    =omegaOf(timeEarliest,self%massRoot)
    if     (                                                                                 &
         &   abs(omegaOf(timeRoot    ,self%massResolution)/omegaRoot-1.0d0) > toleranceOmega &
         &  .or.                                                                             &
         &   abs(omegaOf(timeEarliest,self%massResolution)/omegaEnd -1.0d0) > toleranceOmega &
         & ) call Error_Report('ω depends on mass, which this task does not support'//{introspection:location})
    if (hasTarget) then
       timeTarget =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftTarget))
       omegaTarget=omegaOf(timeTarget,self%massRoot)
       if (abs(omegaOf(timeTarget,self%massResolution)/omegaTarget-1.0d0) > toleranceOmega) &
            & call Error_Report('ω depends on mass, which this task does not support'//{introspection:location})
    else
       timeTarget =timeEarliest
       omegaTarget=omegaEnd
    end if
    ! Grids, uniform in ln M from the resolution to the root mass, and in ω from the root to the maximum redshift.
    countMass  =max(ceiling(log10(self%massRoot/self%massResolution)*dble(self%countMassPerDecade)),3)
    stepLogMass=log(self%massRoot/self%massResolution)/dble(countMass)
    allocate(logMass(0:countMass))
    allocate(mass   (0:countMass))
    do i=0,countMass
       logMass(i)=log(self%massResolution)+dble(i)*stepLogMass
    end do
    logMass(countMass)=log(self%massRoot)
    mass              =exp(logMass)
    ! The ω grid is uniform, except that, if statistic A is computed, ω₁ lies on it, with uniform steps on either side.
    if (hasTarget) then
       kTarget       =max(1,min(self%countOmega-1,nint(dble(self%countOmega)*(omegaTarget-omegaRoot)/(omegaEnd-omegaRoot))))
       stepOmegaBelow=(omegaTarget-omegaRoot  )/dble(                kTarget)
       stepOmegaAbove=(omegaEnd   -omegaTarget)/dble(self%countOmega-kTarget)
    else
       kTarget       =self%countOmega
       stepOmegaBelow=(omegaEnd   -omegaRoot  )/dble(self%countOmega        )
       stepOmegaAbove=stepOmegaBelow
    end if
    allocate(omega   (0:self%countOmega))
    allocate(redshift(0:self%countOmega))
    do k=0,self%countOmega
       if (k <= kTarget) then
          omega(k)=omegaRoot  +dble(k        )*stepOmegaBelow
       else
          omega(k)=omegaTarget+dble(k-kTarget)*stepOmegaAbove
       end if
       if (hasTarget .and. k == kTarget) omega(k)=omegaTarget
       redshift(k)=self%cosmologyFunctions_%redshiftFromExpansionFactor(                                                             &
            &       self%cosmologyFunctions_%expansionFactor            (                                                            &
            &        self%criticalOverdensity_%timeOfCollapse            (criticalOverdensity=omega(k),mass=self%massRoot,node=node) & 
            &                                                           )                                                            &
            &                                                          )
    end do
    ! Quadrature nodes in ln M′ for each parent mass: the grid points up to M/2, the point M/2, and each threshold between the
    ! resolution and M/2 (twice, as the limits of the intervals on either side of the step in the count added by a split).
    ! Grid points within a tolerance of one of the other points are replaced by it.
    allocate(logMassCandidate  (countMass+2*countThresholds+3))
    allocate(thresholdCandidate(countMass+2*countThresholds+3))
    allocate(sideCandidate     (countMass+2*countThresholds+3))
    allocate(nodeStart(0:countMass))
    allocate(nodeCount(0:countMass))
    countNodes=0
    do i=0,countMass
       nodeStart(i)=countNodes+1
       nodeCount(i)=0
       if (mass(i) > 2.0d0*self%massResolution*(1.0d0+toleranceMass)) nodeCount(i)=quadratureNodes(i,count=.true.)
       countNodes  =countNodes+nodeCount(i)
    end do
    allocate(logMassNode     (  countNodes                ))
    allocate(logMassPrimaryNode(countNodes                ))
    allocate(weightNode      (  countNodes                ))
    allocate(thresholdNode   (  countNodes                ))
    allocate(sideNode        (  countNodes                ))
    allocate(baseSecondary   (  countNodes                ))
    allocate(basePrimary     (  countNodes                ))
    allocate(weightsSecondary(4,countNodes                ))
    allocate(weightsPrimary  (4,countNodes                ))
    allocate(countAdded      (  countNodes,countThresholds))
    do i=0,countMass
       if (nodeCount(i) > 0) j=quadratureNodes(i,count=.false.)
    end do
    message=var_str('mass grid: ')//(countMass+1)//' points; ω grid: '
    message=message//(self%countOmega+1)//' points; quadrature nodes: '
    message=message//countNodes
    call displayMessage(message)
    ! For the densities of branches and the decomposition of the variance: the index of each node in a lattice of mass ratios,
    ! x=M′/M=exp(-d Δln M), and, for each point of the grid, the list of the contributions to it of the splits of each parent
    ! (from the stencils of the primary and of the secondary, and from the loss of the parent itself), so that the transpose of
    ! the split operator is a gather in a fixed order (making the results independent of the number of threads).
    if (self%decomposeVariance) then
       allocate(latticeNode(countNodes   ))
       allocate(entryFill  (0:countMass  ))
       allocate(entryStart (0:countMass+1))
       entryFill=0
       do i=0,countMass
          do q=nodeStart(i),nodeStart(i)+nodeCount(i)-1
             latticeNode(q)=max(0,min(countMass,nint((logMass(i)-logMassNode(q))/stepLogMass)))
             entryFill(basePrimary  (q)-1:basePrimary  (q)+2)=entryFill(basePrimary  (q)-1:basePrimary  (q)+2)+1
             entryFill(baseSecondary(q)-1:baseSecondary(q)+2)=entryFill(baseSecondary(q)-1:baseSecondary(q)+2)+1
             entryFill(i                                    )=entryFill(i                                    )+1
          end do
       end do
       entryStart(0)=1
       do p=0,countMass
          entryStart(p+1)=entryStart(p)+entryFill(p)
       end do
       countEntries=entryStart(countMass+1)-1
       allocate(entryNode     (countEntries))
       allocate(entryWeight   (countEntries))
       allocate(entrySecondary(countEntries))
       entryFill=0
       do i=0,countMass
          do q=nodeStart(i),nodeStart(i)+nodeCount(i)-1
             do l=1,4
                call addEntry(basePrimary  (q)-2+l,q,weightsPrimary  (l,q),.false.)
             end do
             do l=1,4
                call addEntry(baseSecondary(q)-2+l,q,weightsSecondary(l,q),.true. )
             end do
             call addEntry(i,q,-1.0d0,.false.)
          end do
       end do
       allocate(massRatio(0:countMass))
       do l=0,countMass
          massRatio(l)=exp(-dble(l)*stepLogMass)
       end do
    end if
    ! Allocate the solution and its tables.
    allocate(countMean                    (0:countMass                                          ,countCases))
    allocate(countVariance                (0:countMass                                          ,countCases))
    allocate(sourceMean                   (0:countMass                                          ,countCases))
    allocate(sourceVariance               (0:countMass                                          ,countCases))
    allocate(countMeanTable               (0:countMass ,0:self%countOmega/self%outputStrideOmega,countCases))
    allocate(countVarianceTable           (0:countMass ,0:self%countOmega/self%outputStrideOmega,countCases))
    allocate(sourceSquaredTable           (0:countMass ,0:self%countOmega/self%outputStrideOmega,countCases))
    allocate(fractionSubresolution        (0:countMass ,0:self%countOmega/self%outputStrideOmega           ))
    allocate(rateSplit                    (0:countMass ,0:self%countOmega/self%outputStrideOmega           ))
    allocate(rateSplitQuadrature          (0:countMass ,0:self%countOmega/self%outputStrideOmega           ))
    allocate(kernelCurrent                (  countNodes                                                    ))
    allocate(fractionPrevious             (0:countMass                                                     ))
    allocate(fractionCurrent              (0:countMass                                                     ))
    allocate(rateCurrent                  (0:countMass                                                     ))
    allocate(rateQuadratureCurrent        (0:countMass                                                     ))
    allocate(terminatedPrevious           (0:countMass                                                     ))
    allocate(terminatedCurrent            (0:countMass                                                     ))
    allocate(logMassFoot                  (0:countMass                                                     ))
    allocate(countMeanFoot                (0:countMass                                                     ))
    allocate(countVarianceFoot            (0:countMass                                                     ))
    allocate(sourceMeanFoot               (0:countMass                                                     ))
    allocate(sourceVarianceFoot           (0:countMass                                                     ))
    allocate(countMeanPredicted           (0:countMass                                                     ))
    allocate(countVariancePredicted       (0:countMass                                                     ))
    allocate(sourceMeanPredicted          (0:countMass                                                     ))
    allocate(sourceVariancePredicted      (0:countMass                                                     ))
    allocate(sourceSquaredPredicted       (0:countMass                                                     ))
    allocate(sourceSquaredCurrent         (0:countMass                                                     ))
    allocate(massBoundaryTable            (             0:self%countOmega/self%outputStrideOmega,countCases))
    if (self%decomposeVariance) then
       countBlocks=self%countOmega/self%outputStrideOmega
       allocate(kernelNext                (  countNodes                                 ))
       allocate(rateNode                  (  countNodes                                 ))
       allocate(fractionNext              (0:countMass                                  ))
       allocate(rateNext                  (0:countMass                                  ))
       allocate(rateQuadratureNext        (0:countMass                                  ))
       allocate(terminatedNext            (0:countMass                                  ))
       allocate(weightsPopulation         (0:countMass                                  ))
       allocate(weightsMain               (0:countMass                                  ))
       allocate(splitsPopulation          (0:countMass                                  ))
       allocate(splitsMain                (0:countMass                                  ))
       allocate(sourceSquaredTerminal     (0:countMass                                  ))
       allocate(sourceAdded               (0:countMass                                  ))
       allocate(densityPopulation         (0:countMass ,0:self%countOmega               ))
       allocate(densityPopulationPredicted(0:countMass ,0:self%countOmega               ))
       allocate(densityMain               (0:countMass ,0:self%countOmega               ))
       allocate(densityMainPredicted      (0:countMass ,0:self%countOmega               ))
       allocate(squaredRatio              (0:countMass ,0:countMass        ,2           ))
       allocate(ratioActive               (0:countMass                                  ))
       allocate(contributionsMass         (0:countMass ,  countBlocks      ,4,countCases))
       allocate(contributionsRatio        (0:countMass ,  countBlocks      ,4,countCases))
       allocate(stepVariance              (             0:self%countOmega-1  ,countCases))
       allocate(stepAdded                 (             0:self%countOmega-1  ,countCases))
       allocate(aliveMean                 (             0:self%countOmega    ,countCases))
       allocate(aliveVariance             (             0:self%countOmega    ,countCases))
       sourceSquaredTerminal     =0.0d0
       sourceAdded               =0.0d0
       densityPopulation         =0.0d0
       densityPopulationPredicted=0.0d0
       densityMain               =0.0d0
       densityMainPredicted      =0.0d0
       squaredRatio              =0.0d0
       ratioActive               =.false.
       contributionsMass         =0.0d0
       contributionsRatio        =0.0d0
       stepVariance              =0.0d0
       stepAdded                 =0.0d0
       aliveMean                 =0.0d0
       aliveVariance             =0.0d0
    end if
    ! Nodes near the boundary of the support of statistic A, which is followed along its characteristic. In a layer of width
    ! about M_res above it no split leaves a progenitor in the support, so that the slope of the solution changes sharply across
    ! the layer. Its width in ln M is M_res/M_c, which can be much less than the step of the grid, so the layer is resolved by
    ! these nodes, which follow their characteristics from offsets above the threshold at the target epoch. The first is the
    ! boundary itself.
    allocate(countNear      (               countCases))
    allocate(countNearActive(               countCases))
    allocate(offsetNear     (countNearLimit,countCases))
    countNear =0
    offsetNear=0.0d0
    do iCase=1,countCases
       if (caseStatistic(iCase) /= 'A') cycle
       countNear(iCase)=1
       offset          =0.0d0
       spacing         =spacingNearFirst
       do while (countNear(iCase) < countNearLimit)
          spacing                          =min(spacing,spacingNearGrid*stepLogMass*(self%massThresholds(caseThreshold(iCase))+offset*self%massResolution)/self%massResolution)
          offset                           =offset          +spacing
          countNear       (iCase)          =countNear(iCase)+1
          offsetNear      (countNear(iCase),iCase)=offset
          if     (                                                                                                                     &
               &   offset                                                                          >= extentNearResolution             &
               &  .and.                                                                                                                &
               &   log(1.0d0+offset*self%massResolution/self%massThresholds(caseThreshold(iCase))) >= extentNear          *stepLogMass &
               & ) exit
          spacing                          =spacing         *growthNear
       end do
    end do
    countNearMaximum=max(1,maxval(countNear))
    allocate(logMassNear                (            countNearMaximum,countCases))
    allocate(countMeanNear              (            countNearMaximum,countCases))
    allocate(countVarianceNear          (            countNearMaximum,countCases))
    allocate(sourceMeanNear             (            countNearMaximum,countCases))
    allocate(sourceVarianceNear         (            countNearMaximum,countCases))
    allocate(logMassNearPrevious        (            countNearMaximum,countCases))
    allocate(countMeanNearPrevious      (            countNearMaximum,countCases))
    allocate(countVarianceNearPrevious  (            countNearMaximum,countCases))
    allocate(sourceMeanNearPrevious     (            countNearMaximum,countCases))
    allocate(sourceVarianceNearPrevious (            countNearMaximum,countCases))
    allocate(countMeanNearPredicted     (            countNearMaximum           ))
    allocate(countVarianceNearPredicted (            countNearMaximum           ))
    allocate(sourceMeanNearPredicted    (            countNearMaximum           ))
    allocate(sourceVarianceNearPredicted(            countNearMaximum           ))
    allocate(countNodesNear             (            countNearMaximum           ))
    allocate(logMassNodeNear            (countMass+2,countNearMaximum           ))
    allocate(kernelNear                 (countMass+2,countNearMaximum           ))
    allocate(logMassPrimaryNodeNear     (countMass+2,countNearMaximum           ))
    countNearActive            =0
    countNodesNear             =0
    logMassNear                =0.0d0
    countMeanNear              =0.0d0
    countVarianceNear          =0.0d0
    sourceMeanNear             =0.0d0
    sourceVarianceNear         =0.0d0
    countMeanNearPredicted     =0.0d0
    countVarianceNearPredicted =0.0d0
    sourceMeanNearPredicted    =0.0d0
    sourceVarianceNearPredicted=0.0d0
    ! The step in ω used in checking the stability at the maximum redshift.
    stepOmega=stepOmegaAbove
    ! The solver runs in a single parallel region, with its loops over the mass grid shared between the threads and the remaining
    ! work done by one thread (a parallel region per loop is prohibitively slow, as each one copies the event hook filters
    ! between threads). Each thread works with its own copies of the objects which compute the rates (they memoize results and
    ! build tables), and its own node. These are threadprivate so that the procedures contained in this one can reach them.
    !$omp parallel private(basic,k,iCase)
    allocate(criticalOverdensity_           ,mold=self%criticalOverdensity_           )
    allocate(cosmologicalMassVariance_      ,mold=self%cosmologicalMassVariance_      )
    allocate(mergerTreeBranchingProbability_,mold=self%mergerTreeBranchingProbability_)
    !$omp critical(taskMergerTreeBranchingMomentsDeepCopy)
    !![
    <deepCopyReset variables="self%criticalOverdensity_ self%cosmologicalMassVariance_ self%mergerTreeBranchingProbability_"/>
    <deepCopy source="self%criticalOverdensity_           " destination="criticalOverdensity_           "/>
    <deepCopy source="self%cosmologicalMassVariance_      " destination="cosmologicalMassVariance_      "/>
    <deepCopy source="self%mergerTreeBranchingProbability_" destination="mergerTreeBranchingProbability_"/>
    <deepCopyFinalize variables="criticalOverdensity_ cosmologicalMassVariance_ mergerTreeBranchingProbability_"/>
    !!]
    !$omp end critical(taskMergerTreeBranchingMomentsDeepCopy)
    !$omp barrier
    allocate(parameters)
    parameters=inputParameters(self%parameters)
    call Node_Components_Thread_Initialize(parameters)
    nodeThread => treeNode        (                 )
    basic      => nodeThread%basic(autoCreate=.true.)
    ! The densities of branches, from the root forward in ω (see the description of this task).
    if (self%decomposeVariance) then
       call levelEvaluate(0,kernelCurrent,fractionCurrent,rateCurrent,rateQuadratureCurrent,terminatedCurrent)
       !$omp single
       densityPopulation(countMass,0)=1.0d0
       densityMain      (countMass,0)=1.0d0
       !$omp end single
       do k=0,self%countOmega-1
          !$omp single
          stepOmega=merge(stepOmegaBelow,stepOmegaAbove,k < kTarget)
          !$omp end single
          call levelEvaluate(k+1,kernelNext,fractionNext,rateNext,rateQuadratureNext,terminatedNext)
          call characteristicFeet(fractionCurrent,terminatedCurrent,fractionNext,terminatedNext)
          ! Splits at ω_k, of the branches which are not terminated there.
          !$omp single
          weightsPopulation=merge(0.0d0,densityPopulation(:,k),terminatedCurrent)
          weightsMain      =merge(0.0d0,densityMain      (:,k),terminatedCurrent)
          !$omp end single
          call splitsTranspose(kernelCurrent,weightsPopulation,.false.,splitsPopulation)
          call splitsTranspose(kernelCurrent,weightsMain      ,.true. ,splitsMain      )
          ! Transport along the characteristics to ω_{k+1}: the predictor (after a full step of splits at ω_k), and the first
          ! part of the corrector (after half of a step).
          !$omp single
          where (terminatedCurrent)
             splitsPopulation=0.0d0
             splitsMain      =0.0d0
          end where
          call transportTranspose(weightsPopulation+      stepOmega*splitsPopulation,densityPopulationPredicted(:,k+1))
          call transportTranspose(weightsMain      +      stepOmega*splitsMain      ,densityMainPredicted      (:,k+1))
          call transportTranspose(weightsPopulation+0.5d0*stepOmega*splitsPopulation,densityPopulation         (:,k+1))
          call transportTranspose(weightsMain      +0.5d0*stepOmega*splitsMain      ,densityMain               (:,k+1))
          !$omp end single
          ! The second part of the corrector: half of a step of splits at ω_{k+1}, of the predicted densities.
          call splitsTranspose(kernelNext,densityPopulationPredicted(:,k+1),.false.,splitsPopulation)
          call splitsTranspose(kernelNext,densityMainPredicted      (:,k+1),.true. ,splitsMain      )
          !$omp single
          densityPopulation    (:,k+1)=densityPopulation(:,k+1)+0.5d0*stepOmega*splitsPopulation
          densityMain          (:,k+1)=densityMain      (:,k+1)+0.5d0*stepOmega*splitsMain
          kernelCurrent               =kernelNext
          fractionCurrent             =fractionNext
          rateCurrent                 =rateNext
          rateQuadratureCurrent       =rateQuadratureNext
          terminatedCurrent           =terminatedNext
          !$omp end single
       end do
    end if
    ! The end condition: no count from branches at the maximum redshift.
    call levelEvaluate(self%countOmega,kernelCurrent,fractionCurrent,rateCurrent,rateQuadratureCurrent,terminatedCurrent)
    !$omp single
    countMean         =0.0d0
    countVariance     =0.0d0
    countMeanTable    =0.0d0
    countVarianceTable=0.0d0
    sourceSquaredTable=0.0d0
    massBoundaryTable =0.0d0
    !$omp end single
    do iCase=1,countCases
       ! Statistic A has no count beyond the target epoch.
       if (caseStatistic(iCase) == 'A') then
          !$omp single
          sourceMean    (:,iCase)=0.0d0
          sourceVariance(:,iCase)=0.0d0
          !$omp end single
       else
          call sources(iCase,kernelCurrent,terminatedCurrent,countMean(:,iCase),countVariance(:,iCase),sourceMean(:,iCase),sourceVariance(:,iCase),sourceSquaredTable(:,self%countOmega/self%outputStrideOmega,iCase))
          if (self%decomposeVariance) call decompositionAccumulate(iCase,self%countOmega,.false.,sourceSquaredTable(:,self%countOmega/self%outputStrideOmega,iCase))
       end if
    end do
    !$omp single
    countMeanTable       (:,self%countOmega/self%outputStrideOmega,:)=countMean
    countVarianceTable   (:,self%countOmega/self%outputStrideOmega,:)=countVariance
    fractionSubresolution(:,self%countOmega/self%outputStrideOmega  )=fractionCurrent
    rateSplit            (:,self%countOmega/self%outputStrideOmega  )=rateCurrent
    rateSplitQuadrature  (:,self%countOmega/self%outputStrideOmega  )=rateQuadratureCurrent
    rateStepMaximum       =maxval(rateQuadratureCurrent)*stepOmega
    if (rateStepMaximum > 1.0d0) call Error_Report('the step in ω is too large for stability (R Δω > 1): increase [countOmega]'//{introspection:location})
    quadratureErrorMaximum=quadratureError(rateCurrent,rateQuadratureCurrent)
    !$omp end single
    ! Step back to the root.
    do k=self%countOmega-1,0,-1
       !$omp single
       fractionPrevious  =fractionCurrent
       terminatedPrevious=terminatedCurrent
       stepOmega         =merge(stepOmegaBelow,stepOmegaAbove,k < kTarget)
       !$omp end single
       call levelEvaluate(k,kernelCurrent,fractionCurrent,rateCurrent,rateQuadratureCurrent,terminatedCurrent)
       !$omp single
       if (mod(k,self%outputStrideOmega) == 0) then
          fractionSubresolution(:,k/self%outputStrideOmega)=fractionCurrent
          rateSplit            (:,k/self%outputStrideOmega)=rateCurrent
          rateSplitQuadrature  (:,k/self%outputStrideOmega)=rateQuadratureCurrent
       end if
       rateStepMaximum           =max(rateStepMaximum       ,maxval(rateQuadratureCurrent)*stepOmega         )
       if (rateStepMaximum > 1.0d0) call Error_Report('the step in ω is too large for stability (R Δω > 1): increase [countOmega]'//{introspection:location})
       quadratureErrorMaximum    =max(quadratureErrorMaximum,quadratureError(rateCurrent,rateQuadratureCurrent))
       !$omp end single
       call characteristicFeet(fractionCurrent,terminatedCurrent,fractionPrevious,terminatedPrevious)
       do iCase=1,countCases
          if (caseStatistic(iCase) == 'A') then
             ! Statistic A is defined only up to the target epoch.
             if (k > kTarget) cycle
             if (k == kTarget) then
                ! The end condition at the target epoch: a count of one for each branch of mass at least the threshold.
                !$omp single
                do j=1,countNear(iCase)
                   logMassNear(j,iCase)=log(self%massThresholds(caseThreshold(iCase))+offsetNear(j,iCase)*self%massResolution)
                end do
                countMeanNear    (:,iCase)=1.0d0
                countVarianceNear(:,iCase)=0.0d0
                countNearActive  (  iCase)=countActive(logMassNear(1:countNear(iCase),iCase))
                where (logMass >= logMassNear(1,iCase) .and. .not.terminatedCurrent)
                   countMean(:,iCase)=1.0d0
                elsewhere
                   countMean(:,iCase)=0.0d0
                end where
                countVariance(:,iCase)=0.0d0
                !$omp end single
                call levelEvaluateNear(k,iCase)
             else
                ! Step the nodes near the boundary back along their characteristics, and predict their values. Nodes above the
                ! largest mass of the grid are no longer needed.
                !$omp single
                countNearActivePrevious            =countNearActive   (  iCase)
                logMassNearPrevious       (:,iCase)=logMassNear       (:,iCase)
                countMeanNearPrevious     (:,iCase)=countMeanNear     (:,iCase)
                countVarianceNearPrevious (:,iCase)=countVarianceNear (:,iCase)
                sourceMeanNearPrevious    (:,iCase)=sourceMeanNear    (:,iCase)
                sourceVarianceNearPrevious(:,iCase)=sourceVarianceNear(:,iCase)
                do j=1,countNearActivePrevious
                   logMassNear(j,iCase)=characteristicStep(logMassNearPrevious(j,iCase))
                end do
                countNearActive(iCase)=countActive(logMassNear(1:countNearActivePrevious,iCase))
                if (countNearActivePrevious > 0) indexFirstPrevious=indexFirst(logMassNearPrevious(countNearActivePrevious,iCase))
                countMeanNearPredicted    =countMeanNearPrevious    (:,iCase)+stepOmega*sourceMeanNearPrevious    (:,iCase)
                countVarianceNearPredicted=countVarianceNearPrevious(:,iCase)+stepOmega*sourceVarianceNearPrevious(:,iCase)
                !$omp end single
                ! Values at the feet of the characteristics. Branches below the boundary have no count. The feet of those on or
                ! above it lie on or above its previous position (up to the accuracy of its step, which is enforced).
                !$omp do private(logMassFootTarget,baseFoot,weightsFoot)
                do i=0,countMass
                   if (outsideSupport(i,iCase)) then
                      countMeanFoot     (i)=0.0d0
                      countVarianceFoot (i)=0.0d0
                      sourceMeanFoot    (i)=0.0d0
                      sourceVarianceFoot(i)=0.0d0
                   else
                      logMassFootTarget=max(logMassFoot(i),logMassNearPrevious(1,iCase))
                      call stencil(logMassFootTarget,baseFoot,weightsFoot)
                      countMeanFoot     (i)=valueSupported(countMean     (:,iCase),logMassNearPrevious(1:countNearActivePrevious,iCase),countMeanNearPrevious     (1:countNearActivePrevious,iCase),indexFirstPrevious,logMassFootTarget,baseFoot,weightsFoot)
                      countVarianceFoot (i)=valueSupported(countVariance (:,iCase),logMassNearPrevious(1:countNearActivePrevious,iCase),countVarianceNearPrevious (1:countNearActivePrevious,iCase),indexFirstPrevious,logMassFootTarget,baseFoot,weightsFoot)
                      sourceMeanFoot    (i)=valueSupported(sourceMean    (:,iCase),logMassNearPrevious(1:countNearActivePrevious,iCase),sourceMeanNearPrevious    (1:countNearActivePrevious,iCase),indexFirstPrevious,logMassFootTarget,baseFoot,weightsFoot)
                      sourceVarianceFoot(i)=valueSupported(sourceVariance(:,iCase),logMassNearPrevious(1:countNearActivePrevious,iCase),sourceVarianceNearPrevious(1:countNearActivePrevious,iCase),indexFirstPrevious,logMassFootTarget,baseFoot,weightsFoot)
                   end if
                end do
                !$omp end do
                ! Heun step, as for the other statistics, for the grid and for the nodes near the boundary.
                !$omp single
                countMeanPredicted    =countMeanFoot    +stepOmega*sourceMeanFoot
                countVariancePredicted=countVarianceFoot+stepOmega*sourceVarianceFoot
                do i=0,countMass
                   if (outsideSupport(i,iCase)) then
                      countMeanPredicted    (i)=0.0d0
                      countVariancePredicted(i)=0.0d0
                   end if
                end do
                !$omp end single
                call levelEvaluateNear(k,iCase)
                call sourcesTarget(iCase,kernelCurrent,terminatedCurrent,countMeanPredicted,countVariancePredicted,countMeanNearPredicted,countVarianceNearPredicted,sourceMeanPredicted,sourceVariancePredicted,sourceSquaredPredicted,sourceMeanNearPredicted,sourceVarianceNearPredicted)
                if (self%decomposeVariance) call decompositionAccumulate(iCase,k,.true.,sourceSquaredPredicted)
                !$omp single
                countMean    (:,iCase)=countMeanFoot    +0.5d0*stepOmega*(sourceMeanFoot    +sourceMeanPredicted    )
                countVariance(:,iCase)=countVarianceFoot+0.5d0*stepOmega*(sourceVarianceFoot+sourceVariancePredicted)
                do i=0,countMass
                   if (outsideSupport(i,iCase)) then
                      countMean    (i,iCase)=0.0d0
                      countVariance(i,iCase)=0.0d0
                   end if
                end do
                countMeanNear    (:,iCase)=countMeanNearPrevious    (:,iCase)+0.5d0*stepOmega*(sourceMeanNearPrevious    (:,iCase)+sourceMeanNearPredicted    )
                countVarianceNear(:,iCase)=countVarianceNearPrevious(:,iCase)+0.5d0*stepOmega*(sourceVarianceNearPrevious(:,iCase)+sourceVarianceNearPredicted)
                !$omp end single
             end if
             call sourcesTarget(iCase,kernelCurrent,terminatedCurrent,countMean(:,iCase),countVariance(:,iCase),countMeanNear(:,iCase),countVarianceNear(:,iCase),sourceMean(:,iCase),sourceVariance(:,iCase),sourceSquaredCurrent,sourceMeanNear(:,iCase),sourceVarianceNear(:,iCase))
             if (self%decomposeVariance) call decompositionAccumulate(iCase,k,.false.,sourceSquaredCurrent)
             !$omp single
             if (mod(k,self%outputStrideOmega) == 0) then
                countMeanTable    (:,k/self%outputStrideOmega,iCase)=countMean    (:,iCase)
                countVarianceTable(:,k/self%outputStrideOmega,iCase)=countVariance(:,iCase)
                sourceSquaredTable(:,k/self%outputStrideOmega,iCase)=sourceSquaredCurrent
                if (countNearActive(iCase) > 0) massBoundaryTable(k/self%outputStrideOmega,iCase)=exp(logMassNear(1,iCase))
             end if
             !$omp end single
             cycle
          end if
          ! Values at the feet of the characteristics. A branch below the resolution at ω_{k+1} contributes nothing.
          !$omp do
          do i=0,countMass
             if (logMassFoot(i) <= logMass(0)) then
                countMeanFoot     (i)=0.0d0
                countVarianceFoot (i)=0.0d0
                sourceMeanFoot    (i)=0.0d0
                sourceVarianceFoot(i)=0.0d0
             else
                countMeanFoot     (i)=interpolateAt(countMean     (:,iCase),logMassFoot(i))
                countVarianceFoot (i)=interpolateAt(countVariance (:,iCase),logMassFoot(i))
                sourceMeanFoot    (i)=interpolateAt(sourceMean    (:,iCase),logMassFoot(i))
                sourceVarianceFoot(i)=interpolateAt(sourceVariance(:,iCase),logMassFoot(i))
             end if
          end do
          !$omp end do
          ! Heun step: predict with the sources at the foot, then correct with the mean of those and the sources at the
          ! predicted values.
          !$omp single
          countMeanPredicted    =countMeanFoot    +stepOmega*sourceMeanFoot
          countVariancePredicted=countVarianceFoot+stepOmega*sourceVarianceFoot
          where (terminatedCurrent)
             countMeanPredicted    =0.0d0
             countVariancePredicted=0.0d0
          end where
          !$omp end single
          call sources(iCase,kernelCurrent,terminatedCurrent,countMeanPredicted,countVariancePredicted,sourceMeanPredicted,sourceVariancePredicted,sourceSquaredPredicted)
          if (self%decomposeVariance) call decompositionAccumulate(iCase,k,.true.,sourceSquaredPredicted)
          !$omp single
          countMean    (:,iCase)=countMeanFoot    +0.5d0*stepOmega*(sourceMeanFoot    +sourceMeanPredicted    )
          countVariance(:,iCase)=countVarianceFoot+0.5d0*stepOmega*(sourceVarianceFoot+sourceVariancePredicted)
          where (terminatedCurrent)
             countMean    (:,iCase)=0.0d0
             countVariance(:,iCase)=0.0d0
          end where
          !$omp end single
          call sources(iCase,kernelCurrent,terminatedCurrent,countMean(:,iCase),countVariance(:,iCase),sourceMean(:,iCase),sourceVariance(:,iCase),sourceSquaredCurrent)
          if (self%decomposeVariance) call decompositionAccumulate(iCase,k,.false.,sourceSquaredCurrent)
          !$omp single
          if (mod(k,self%outputStrideOmega) == 0) then
             countMeanTable    (:,k/self%outputStrideOmega,iCase)=countMean    (:,iCase)
             countVarianceTable(:,k/self%outputStrideOmega,iCase)=countVariance(:,iCase)
             sourceSquaredTable(:,k/self%outputStrideOmega,iCase)=sourceSquaredCurrent
          end if
          !$omp end single
       end do
    end do
    call nodeThread%destroy()
    deallocate(nodeThread)
    !![
    <objectDestructor name="criticalOverdensity_"           />
    <objectDestructor name="cosmologicalMassVariance_"      />
    <objectDestructor name="mergerTreeBranchingProbability_"/>
    !!]
    call Node_Components_Thread_Uninitialize()
    deallocate(parameters)
    !$omp end parallel
    write (label,'(f8.4)') rateStepMaximum
    message='maximum R Δω: '//trim(adjustl(label))
    write (label,'(e10.3)') quadratureErrorMaximum
    message=message//'; maximum fractional difference of quadrature and class split rates: '//trim(adjustl(label))
    call displayMessage(message)
    ! The cumulative variance, V(ω)=Var(E[N|ω]), and the count added by splits up to ω. With the solution of the equations, the
    ! mean and variance of the count are then the sums of these and of the mean and variance expected from the branches alive at
    ! ω (exactly, for statistics B and C, whose densities are found by the transpose of the scheme used for the solution): the
    ! largest departures from these identities are found.
    if (self%decomposeVariance) then
       allocate(varianceCumulative      (0:self%countOmega,countCases))
       allocate(countAddedCumulative    (0:self%countOmega,countCases))
       allocate(identityResidualMean    (                  countCases))
       allocate(identityResidualVariance(                  countCases))
       do iCase=1,countCases
          varianceCumulative  (0,iCase)=0.0d0
          countAddedCumulative(0,iCase)=0.0d0
          do k=1,self%countOmega
             varianceCumulative  (k,iCase)=varianceCumulative  (k-1,iCase)+stepVariance(k-1,iCase)
             countAddedCumulative(k,iCase)=countAddedCumulative(k-1,iCase)+stepAdded   (k-1,iCase)
          end do
          kEnd=self%countOmega
          if (caseStatistic(iCase) == 'A') kEnd=kTarget
          identityResidualMean    (iCase)=maxval(abs(countMeanTable    (countMass,0,iCase)-countAddedCumulative(0:kEnd,iCase)-aliveMean    (0:kEnd,iCase)))
          identityResidualVariance(iCase)=maxval(abs(countVarianceTable(countMass,0,iCase)-varianceCumulative  (0:kEnd,iCase)-aliveVariance(0:kEnd,iCase)))
       end do
    end if
    ! Output.
    allocate(rootVarianceNow(0:countMass))
    do i=0,countMass
       rootVarianceNow(i)=self%cosmologicalMassVariance_%rootVariance(time=timeNow,mass=mass(i))
    end do
    !$ call hdf5Access%set()
    outputGroup=outputFile%openGroup(char(self%outputGroupName),'Mean and variance of subhalo counts in merger trees from the moment equations.')
    call    outputGroup%writeAttribute(self%massRoot                                                                     ,'massRoot'                                                                                            )
    call    outputGroup%writeAttribute(self%massResolution                                                               ,'massResolution'                                                                                      )
    call    outputGroup%writeAttribute(self%redshiftRoot                                                                 ,'redshiftRoot'                                                                                        )
    call    outputGroup%writeAttribute(     omegaRoot                                                                    ,'omegaRoot'                                                                                           )
    call    outputGroup%writeAttribute(     omegaEnd                                                                     ,'omegaEnd'                                                                                            )
    if (hasTarget) then
       call outputGroup%writeAttribute(self%redshiftTarget                                                               ,'redshiftTarget'                                                                                      )
       call outputGroup%writeAttribute(     omegaTarget                                                                  ,'omegaTarget'                                                                                         )
    end if
    call    outputGroup%writeDataset  (     mass                                                                         ,'mass'                 ,'Masses of the grid.'                                                         )
    call    outputGroup%writeDataset  (     rootVarianceNow                                                              ,'rootVariance'         ,'Root variance of the density field at the present day, σ(M,t₀), on the grid.')
    call    outputGroup%writeDataset  (     omega                (0:self%countOmega:self%outputStrideOmega)              ,'omega'                ,'Values of ω of the grid.'                                                    )
    call    outputGroup%writeDataset  (     redshift             (0:self%countOmega:self%outputStrideOmega)              ,'redshift'             ,'Redshifts of the ω grid.'                                                    )
    call    outputGroup%writeDataset  (     fractionSubresolution                                                        ,'fractionSubresolution','Fractional rate of mass loss to subresolution halos, per unit ω.'            )
    call    outputGroup%writeDataset  (     rateSplit                                                                    ,'rateSplit'            ,'Rate of splits, per unit ω, from the branching probability class.'           )
    call    outputGroup%writeDataset  (     rateSplitQuadrature                                                          ,'rateSplitQuadrature'  ,'Rate of splits, per unit ω, from the quadrature used to solve the equations.')
    if (self%decomposeVariance) then
       call outputGroup%writeDataset  (     densityPopulation    (:,0:self%countOmega:self%outputStrideOmega)/stepLogMass,'densityBranches'      ,'Expected number of branches per unit ln M.'                                  )
       call outputGroup%writeDataset  (     densityMain          (:,0:self%countOmega:self%outputStrideOmega)/stepLogMass,'densityMainBranch'    ,'Probability density of the mass of the main branch per unit ln M.'           )
       call outputGroup%writeDataset  (     massRatio                                                                    ,'massRatio'            ,'Mass ratios, M′/M, of the lattice of the decomposition.'                     )
    end if
    do iCase=1,countCases
       statisticGroup=outputGroup%openGroup('statistic'//caseStatistic(iCase),'Results for one statistic.')
       groupName     =var_str('threshold')
       groupName     =groupName//(caseThreshold(iCase)-1)
       caseGroup     =statisticGroup%openGroup(char(groupName),'Results for one statistic and threshold.')
       call caseGroup%writeAttribute(caseStatistic      (                          iCase )                     ,'statistic'                          )
       call caseGroup%writeAttribute(self%massThresholds(            caseThreshold(iCase))                     ,'massThreshold'                      )
       call caseGroup%writeAttribute(countMeanTable     (countMass,0,              iCase )                     ,'countMeanRoot'                      )
       call caseGroup%writeAttribute(countVarianceTable (countMass,0,              iCase )                     ,'countVarianceRoot'                  )
       call caseGroup%writeDataset  (countMeanTable     (:        ,:,              iCase ),'countMean'         ,'Mean count, n₁.'                    )
       call caseGroup%writeDataset  (countVarianceTable (:        ,:,              iCase ),'countVariance'     ,'Variance of the count, W.'          )
       call caseGroup%writeDataset  (sourceSquaredTable (:        ,:,              iCase ),'sourceDeltaSquared','Source of the variance, ∫ f Δ² dM′.')
       if (caseStatistic(iCase) == 'A')                                                                                                      &
            & call caseGroup%writeDataset(massBoundaryTable(:,iCase),'massBoundary','The smallest mass with a count, M_c (zero beyond ω₁).')
       if (self%decomposeVariance) then
          call caseGroup%writeAttribute(identityResidualMean    (iCase)                                         ,'identityResidualMean'                                                                     )
          call caseGroup%writeAttribute(identityResidualVariance(iCase)                                         ,'identityResidualVariance'                                                                 )
          call caseGroup%writeDataset  (varianceCumulative      (0:self%countOmega:self%outputStrideOmega,iCase),'varianceCumulative'            ,'Cumulative variance, V(ω)=Var(E[N|ω]).'                  )
          call caseGroup%writeDataset  (contributionsMass       (:,:,:,iCase)                                   ,'varianceContributionsMass'     ,'Contributions of splits to the variance, by parent mass.')
          call caseGroup%writeDataset  (contributionsRatio      (:,:,:,iCase)                                   ,'varianceContributionsMassRatio','Contributions of splits to the variance, by mass ratio.' )
       end if
       write (label,'(es10.3)') self%massThresholds(caseThreshold(iCase))
       message='statistic '//caseStatistic(iCase)//', m = '//trim(adjustl(label))
       write (label,'(f12.5)') countMeanTable(countMass,0,iCase)
       message=message//': ⟨N⟩ = '//trim(adjustl(label))
       write (label,'(f12.5)') countVarianceTable(countMass,0,iCase)
       message=message//', Var(N) = '//trim(adjustl(label))
       ! The count is exactly zero for a threshold of at least half of the root mass, as no secondary can reach it.
       if (countMeanTable(countMass,0,iCase) > 0.0d0) then
          write (label,'(f8.5)') countVarianceTable(countMass,0,iCase)/countMeanTable(countMass,0,iCase)
       else
          label='n/a'
       end if
       message=message//', Var/⟨N⟩ = '//trim(adjustl(label))
       call displayMessage(message)
    end do
    !$ call hdf5Access%unset()
    call node%destroy()
    deallocate(node)
    call Node_Components_Thread_Uninitialize()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: merger tree branching moments')
    return

  contains

    subroutine addEntry(p_,q_,weight,secondary)
      !!{RST
      Add an entry to the list of contributions to point ``p_`` of the grid from the splits at quadrature node ``q_``.
      !!}
      implicit none
      integer         , intent(in   ) :: p_       , q_
      double precision, intent(in   ) :: weight
      logical         , intent(in   ) :: secondary
      integer                         :: e

      e                =entryStart(p_)+entryFill(p_)
      entryFill  (p_)  =entryFill (p_)+1
      entryNode     (e)=q_
      entryWeight   (e)=weight
      entrySecondary(e)=secondary
      return
    end subroutine addEntry

    double precision function omegaOf(time,mass)
      !!{RST
      Return the ω of the :cite:t:`cole_hierarchical_2000` merger tree builder, :math:`\delta_\mathrm{c}(t,M)\sigma(M,t_0)/\sigma(M,t)`.
      !!}
      implicit none
      double precision, intent(in   ) :: time, mass

      omegaOf=+self%criticalOverdensity_     %value       (time=time   ,mass=mass,node=node) &
           &  *self%cosmologicalMassVariance_%rootVariance(time=timeNow,mass=mass          ) &
           &  /self%cosmologicalMassVariance_%rootVariance(time=time   ,mass=mass          )
      return
    end function omegaOf

    integer function quadratureNodes(i,count)
      !!{RST
      Count, or construct, the quadrature nodes in :math:`\ln M^\prime` for the parent of mass ``mass(i)``: the grid points
      up to :math:`M/2`, the point :math:`M/2`, and each threshold below :math:`M/2`, twice (as the limits of the intervals on
      either side of the step in the count added by a split). Grid points within a tolerance of a threshold or of :math:`M/2`
      are replaced by it.
      !!}
      implicit none
      integer         , intent(in   )                                  :: i
      logical         , intent(in   )                                  :: count
      double precision, dimension(countThresholds                    ) :: logMassThreshold
      integer         , dimension(countMass+2*countThresholds+3      ) :: order
      integer                                                          :: j               , l, q, &
           &                                                              qStart          , m
      logical                                                          :: halfIsThreshold

      logMassHalf     =logMass(i)-log(2.0d0)
      logMassThreshold=log(self%massThresholds)
      countCandidates =0
      countingOnly    =count
      do j=0,countMass
         if (logMass(j) >= logMassHalf-toleranceMass) exit
         if (any(abs(logMass(j)-logMassThreshold) <= toleranceMass)) cycle
         call addCandidate(logMass(j),0,0)
      end do
      halfIsThreshold=.false.
      do l=1,countThresholds
         if      (abs(logMassThreshold(l)-logMassHalf) <= toleranceMass) then
            ! The threshold is at the end point: only the interval below it (where no count is added) is included.
            call addCandidate(logMassHalf        ,l,0)
            call addCandidate(logMassHalf        ,l,1)
            halfIsThreshold=.true.
         else if (logMassThreshold(l) < logMassHalf) then
            call addCandidate(logMassThreshold(l),l,0)
            call addCandidate(logMassThreshold(l),l,1)
         end if
      end do
      if (.not.halfIsThreshold) call addCandidate(logMassHalf,0,0)
      quadratureNodes=countCandidates
      q              =countCandidates
      if (count) return
      ! Sort the candidates by mass (with the two nodes of a threshold in the order below, then above, the step).
      do j=1,q
         order(j)=j
      end do
      do j=2,q
         l=order(j)
         m=j
         do while (m > 1)
            if (precedes(order(m-1),l)) exit
            order(m)=order(m-1)
            m       =m-1
         end do
         order(m)=l
      end do
      qStart=nodeStart(i)
      do j=1,q
         logMassNode  (qStart+j-1)=logMassCandidate  (order(j))
         thresholdNode(qStart+j-1)=thresholdCandidate(order(j))
         sideNode     (qStart+j-1)=sideCandidate     (order(j))
      end do
      do q=qStart,qStart+nodeCount(i)-1
         ! Trapezoid weights in ln M′ (the two nodes at a threshold are separated by an interval of zero width).
         weightNode(q)=0.0d0
         if (q > qStart               ) weightNode(q)=weightNode(q)+0.5d0*(logMassNode(q  )-logMassNode(q-1))
         if (q < qStart+nodeCount(i)-1) weightNode(q)=weightNode(q)+0.5d0*(logMassNode(q+1)-logMassNode(q  ))
         ! Interpolation stencils for the secondary, M′, and the primary, M-M′.
         call stencil(logMassNode       (q),baseSecondary(q),weightsSecondary(:,q))
         logMassPrimaryNode(q)=logMass(i)+log(1.0d0-exp(logMassNode(q)-logMass(i)))
         call stencil(logMassPrimaryNode(q),basePrimary  (q),weightsPrimary  (:,q))
         ! Count added by the split, for each threshold: one if M′ ≥ m (at a threshold node, the side of the step).
         do l=1,countThresholds
            if (thresholdNode(q) == l) then
               countAdded(q,l)=sideNode(q)
            else if (logMassNode(q) > logMassThreshold(l)) then
               countAdded(q,l)=1
            else
               countAdded(q,l)=0
            end if
         end do
      end do
      return
    end function quadratureNodes

    subroutine addCandidate(logMass_,threshold,side)
      !!{RST
      Add a candidate quadrature node.
      !!}
      implicit none
      double precision, intent(in   ) :: logMass_
      integer         , intent(in   ) :: threshold, side

      countCandidates=countCandidates+1
      if (countingOnly) return
      logMassCandidate  (countCandidates)=logMass_
      thresholdCandidate(countCandidates)=threshold
      sideCandidate     (countCandidates)=side
      return
    end subroutine addCandidate

    logical function precedes(a,b)
      !!{RST
      Return true if candidate ``a`` precedes candidate ``b``.
      !!}
      implicit none
      integer, intent(in   ) :: a, b

      if (logMassCandidate(a) /= logMassCandidate(b)) then
         precedes=logMassCandidate(a) <  logMassCandidate(b)
      else
         precedes=sideCandidate   (a) <= sideCandidate   (b)
      end if
      return
    end function precedes

    subroutine stencil(logMass_,base,weights)
      !!{RST
      Find the four-point Lagrange interpolation stencil in :math:`\ln M` at ``logMass_``.
      !!}
      implicit none
      double precision, intent(in   )               :: logMass_
      integer         , intent(  out)               :: base
      double precision, intent(  out), dimension(4) :: weights
      double precision                              :: x

      x   =(logMass_-logMass(0))/stepLogMass
      base=max(1,min(countMass-2,floor(x)))
      x   =x-dble(base)
      weights=[                                      &
           &   -x*(x-1.0d0)*(x-2.0d0)/6.0d0        , &
           &   +(x+1.0d0)*(x-1.0d0)*(x-2.0d0)/2.0d0, &
           &   -(x+1.0d0)*x*(x-2.0d0)/2.0d0        , &
           &   +(x+1.0d0)*x*(x-1.0d0)/6.0d0          &
           &  ]
      return
    end subroutine stencil

    double precision function interpolateAt(values,logMass_)
      !!{RST
      Interpolate ``values`` on the mass grid to ``logMass_``.
      !!}
      implicit none
      double precision, intent(in   ), dimension(0:) :: values
      double precision, intent(in   )                :: logMass_
      double precision               , dimension(4 ) :: weights
      integer                                        :: base

      call stencil(logMass_,base,weights)
      interpolateAt=sum(values(base-1:base+2)*weights)
      return
    end function interpolateAt

    double precision function interpolateFraction(fraction,terminated,logMass_)
      !!{RST
      Interpolate the subresolution mass loss rate to ``logMass_``, replacing values at terminated grid points (where it is
      undefined) by that of the nearest point which is not terminated.
      !!}
      implicit none
      double precision, intent(in   ), dimension(0:        ) :: fraction
      logical         , intent(in   ), dimension(0:        ) :: terminated
      double precision, intent(in   )                        :: logMass_
      double precision               , dimension(0:countMass) :: fractionFilled
      integer                                                :: i

      fractionFilled=fraction
      do i=countMass-1,0,-1
         if (terminated(i)) fractionFilled(i)=fractionFilled(i+1)
      end do
      interpolateFraction=interpolateAt(fractionFilled,logMass_)
      return
    end function interpolateFraction

    subroutine characteristicFeet(fractionCurrent_,terminatedCurrent_,fractionPrevious_,terminatedPrevious_)
      !!{RST
      Find the feet of the characteristics through the points of the grid at :math:`\omega_k`: the mass at :math:`\omega_{k+1}`
      of a branch of mass :math:`M_i` at :math:`\omega_k`, with the subresolution mass loss rate averaged over the step.
      !!}
      implicit none
      double precision, intent(in   ), dimension(0:) :: fractionCurrent_  , fractionPrevious_
      logical         , intent(in   ), dimension(0:) :: terminatedCurrent_, terminatedPrevious_
      double precision                               :: logMassFootTrial_
      integer                                        :: i

      !$omp do
      do i=0,countMass
         if (terminatedCurrent_(i)) then
            logMassFoot(i)=-huge(0.0d0)
         else
            logMassFootTrial_=logMass(i)-fractionCurrent_(i)*stepOmega
            if (logMassFootTrial_ > logMass(1)) then
               logMassFoot(i)=logMass(i)-0.5d0*(fractionCurrent_(i)+interpolateFraction(fractionPrevious_,terminatedPrevious_,logMassFootTrial_))*stepOmega
            else
               logMassFoot(i)=logMassFootTrial_
            end if
         end if
      end do
      !$omp end do
      return
    end subroutine characteristicFeet

    subroutine splitsTranspose(kernel,weights,primaryOnly,result)
      !!{RST
      Apply the transpose of the split operator (the change in the count expected from a branch due to its splits, per unit
      :math:`\omega`, with the kernel ``kernel``) to ``weights``, the numbers of branches at the points of the grid. If
      ``primaryOnly`` is true, only the primary progenitor of each split is followed (as for the main branch).
      !!}
      implicit none
      double precision, intent(in   ), dimension(:  ) :: kernel
      double precision, intent(in   ), dimension(0: ) :: weights
      logical         , intent(in   )                 :: primaryOnly
      double precision, intent(  out), dimension(0: ) :: result
      double precision                                :: sum
      integer                                         :: i          , q, &
           &                                             p          , e

      !$omp do schedule(dynamic)
      do i=0,countMass
         do q=nodeStart(i),nodeStart(i)+nodeCount(i)-1
            rateNode(q)=weights(i)*weightNode(q)*kernel(q)
         end do
      end do
      !$omp end do
      !$omp do
      do p=0,countMass
         sum=0.0d0
         do e=entryStart(p),entryStart(p+1)-1
            if (primaryOnly .and. entrySecondary(e)) cycle
            sum=sum+rateNode(entryNode(e))*entryWeight(e)
         end do
         result(p)=sum
      end do
      !$omp end do
      return
    end subroutine splitsTranspose

    subroutine transportTranspose(values,result)
      !!{RST
      Apply the transpose of the interpolation to the feet of the characteristics to ``values``, the numbers of branches at the
      points of the grid at :math:`\omega_k`, giving their numbers at the points of the grid at :math:`\omega_{k+1}`.
      !!}
      implicit none
      double precision, intent(in   ), dimension(0:) :: values
      double precision, intent(  out), dimension(0:) :: result
      double precision               , dimension(4 ) :: weights
      integer                                        :: i      , base

      result=0.0d0
      do i=0,countMass
         if (logMassFoot(i) <= logMass(0)) cycle
         call stencil(logMassFoot(i),base,weights)
         result(base-1:base+2)=result(base-1:base+2)+values(i)*weights
      end do
      return
    end subroutine transportTranspose

    subroutine decompositionAccumulate(iCase,k,predicted,sourceSquared_)
      !!{RST
      Add the contributions of splits at :math:`\omega_k` to the decomposition of the variance, given the source
      ``sourceSquared_``, :math:`\int f\Delta^2\mathrm{d}M^\prime`, found from the predicted solution at :math:`\omega_k`
      (if ``predicted`` is true; these are weighted by the densities of branches at :math:`\omega_k`, and belong to the step
      from :math:`\omega_k` to :math:`\omega_{k+1}`), or from the solution (weighted by the predicted densities, and belonging
      to the step from :math:`\omega_{k-1}` to :math:`\omega_k`). In the latter case, also find the sums of the mean and
      variance over the branches alive at :math:`\omega_k`. All of the work is done in one shared loop (so with one barrier):
      its first two iterations find the sums over branches alive, and the contributions by parent mass and the totals for the
      step, in a fixed order, and the others find the contributions by mass ratio, one point of the lattice each.
      !!}
      implicit none
      integer         , intent(in   )                :: iCase            , k
      logical         , intent(in   )                :: predicted
      double precision, intent(in   ), dimension(0:) :: sourceSquared_
      double precision                               :: factor           , weightMain     , &
           &                                            weightSide       , sumMainTerminal, &
           &                                            sumMainStructural, sumSideTerminal, &
           &                                            sumSideStructural, sumVariance    , &
           &                                            sumAdded
      integer                                        :: i                , d              , &
           &                                            step             , block
      logical                                        :: mainOnly         , contributes

      ! For statistic C only splits of the main branch count.
      mainOnly=caseStatistic(iCase) == 'C'
      ! The solution at the root has no step before it, so only the sums over branches alive are found.
      contributes=predicted .or. k > 0
      if (predicted) then
         step=k
      else
         step=k-1
      end if
      ! Each source is half of the trapezoid rule over its step.
      if (contributes) then
         factor=0.5d0*merge(stepOmegaBelow,stepOmegaAbove,step < kTarget)
         block =step/self%outputStrideOmega+1
      else
         factor=0.0d0
         block =1
      end if
      !$omp do schedule(static)
      do d=-2,countMass
         if      (d == -2) then
            ! Sums over the branches alive.
            if (.not.predicted) then
               if      (mainOnly                    ) then
                  aliveMean    (k,iCase)=sum(densityMain      (:,k)*countMean    (:,iCase))
                  aliveVariance(k,iCase)=sum(densityMain      (:,k)*countVariance(:,iCase))
               else if (caseStatistic(iCase) == 'A') then
                  call aliveTarget(iCase,k)
               else
                  aliveMean    (k,iCase)=sum(densityPopulation(:,k)*countMean    (:,iCase))
                  aliveVariance(k,iCase)=sum(densityPopulation(:,k)*countVariance(:,iCase))
               end if
            end if
         else if (.not.contributes) then
            cycle
         else if (d == -1) then
            ! Contributions by parent mass, and the totals for the step.
            sumVariance=0.0d0
            sumAdded   =0.0d0
            do i=0,countMass
               call decompositionWeights(i,k,predicted,mainOnly,factor,weightMain,weightSide)
               contributionsMass(i,block,1,iCase)=contributionsMass(i,block,1,iCase)+ weightMain            *                          sourceSquaredTerminal(i)
               contributionsMass(i,block,2,iCase)=contributionsMass(i,block,2,iCase)+ weightMain            *(sourceSquared_       (i)-sourceSquaredTerminal(i))
               contributionsMass(i,block,3,iCase)=contributionsMass(i,block,3,iCase)+            weightSide *                          sourceSquaredTerminal(i)
               contributionsMass(i,block,4,iCase)=contributionsMass(i,block,4,iCase)+            weightSide *(sourceSquared_       (i)-sourceSquaredTerminal(i))
               sumVariance                       =sumVariance                       +(weightMain+weightSide)* sourceSquared_       (i)
               sumAdded                          =sumAdded                          +(weightMain+weightSide)* sourceAdded          (i)
            end do
            stepVariance(step,iCase)=stepVariance(step,iCase)+sumVariance
            stepAdded   (step,iCase)=stepAdded   (step,iCase)+sumAdded
         else
            ! Contributions by mass ratio.
            sumMainTerminal  =0.0d0
            sumMainStructural=0.0d0
            sumSideTerminal  =0.0d0
            sumSideStructural=0.0d0
            do i=d,countMass
               call decompositionWeights(i,k,predicted,mainOnly,factor,weightMain,weightSide)
               sumMainTerminal  =sumMainTerminal  +weightMain*squaredRatio(d,i,1)
               sumMainStructural=sumMainStructural+weightMain*squaredRatio(d,i,2)
               sumSideTerminal  =sumSideTerminal  +weightSide*squaredRatio(d,i,1)
               sumSideStructural=sumSideStructural+weightSide*squaredRatio(d,i,2)
            end do
            contributionsRatio(d,block,1,iCase)=contributionsRatio(d,block,1,iCase)+sumMainTerminal
            contributionsRatio(d,block,2,iCase)=contributionsRatio(d,block,2,iCase)+sumMainStructural
            contributionsRatio(d,block,3,iCase)=contributionsRatio(d,block,3,iCase)+sumSideTerminal
            contributionsRatio(d,block,4,iCase)=contributionsRatio(d,block,4,iCase)+sumSideStructural
         end if
      end do
      !$omp end do
      return
    end subroutine decompositionAccumulate

    subroutine decompositionWeights(i,k,predicted,mainOnly,factor,weightMain,weightSide)
      !!{RST
      Return the weights of the contributions of parents of mass ``mass(i)`` at :math:`\omega_k` from splits of the main branch,
      and of side branches (see ``decompositionAccumulate``). Parents with no contributions in this evaluation are excluded
      (their contributions by mass ratio are stale).
      !!}
      implicit none
      integer         , intent(in   ) :: i               , k
      logical         , intent(in   ) :: predicted       , mainOnly
      double precision, intent(in   ) :: factor
      double precision, intent(  out) :: weightMain      , weightSide
      double precision                :: weightPopulation

      if (.not.ratioActive(i)) then
         weightMain=0.0d0
         weightSide=0.0d0
         return
      end if
      if (predicted) then
         weightMain      =factor*densityMain               (i,k)
         weightPopulation=factor*densityPopulation         (i,k)
      else
         weightMain      =factor*densityMainPredicted      (i,k)
         weightPopulation=factor*densityPopulationPredicted(i,k)
      end if
      if (mainOnly) then
         weightSide=0.0d0
      else
         weightSide=weightPopulation-weightMain
      end if
      return
    end subroutine decompositionWeights

    subroutine aliveTarget(iCase,k)
      !!{RST
      Find the sums of the mean and variance of statistic A over the branches alive at :math:`\omega_k`. These jump at the
      boundary of the support, so, below the cell of the first point of the grid above the nodes near the boundary, the
      density of branches (interpolated from the grid) is integrated against the mean and variance at those nodes by the
      trapezoid rule, from the boundary to that cell.
      !!}
      implicit none
      integer         , intent(in   ) :: iCase           , k
      double precision, dimension(4 ) :: weights
      double precision                :: logMassEnd      , densityPrevious, &
           &                             densityHere     , meanHere       , &
           &                             varianceHere    , meanPrevious   , &
           &                             variancePrevious, logMassPrevious, &
           &                             logMassHere
      integer                         :: countNear_      , indexFirst_    , &
           &                             j               , base
      logical                         :: first

      aliveMean    (k,iCase)=0.0d0
      aliveVariance(k,iCase)=0.0d0
      countNear_           =countNearActive(iCase)
      if (countNear_ == 0) return
      indexFirst_=indexFirst(logMassNear(countNear_,iCase))
      if (indexFirst_ <= countMass) then
         aliveMean    (k,iCase)=sum(densityPopulation(indexFirst_:countMass,k)*countMean    (indexFirst_:countMass,iCase))
         aliveVariance(k,iCase)=sum(densityPopulation(indexFirst_:countMass,k)*countVariance(indexFirst_:countMass,iCase))
         logMassEnd            =logMass(indexFirst_)-0.5d0*stepLogMass
      else
         ! The nodes near the boundary reach into the cell of the root. Its weight is counted with the value on the grid there (as
         ! for the other points of the grid above), and the trapezoid rule ends at the bottom of its cell.
         if (logMass(countMass) >= logMassNear(1,iCase)) then
            aliveMean    (k,iCase)=densityPopulation(countMass,k)*countMean    (countMass,iCase)
            aliveVariance(k,iCase)=densityPopulation(countMass,k)*countVariance(countMass,iCase)
         end if
         logMassEnd=logMass(countMass)-0.5d0*stepLogMass
      end if
      first=.true.
      do j=1,countNear_+1
         if (j <= countNear_) then
            if (logMassNear(j,iCase) >= logMassEnd) cycle
            logMassHere =logMassNear      (j,iCase)
            meanHere    =countMeanNear    (j,iCase)
            varianceHere=countVarianceNear(j,iCase)
         else
            logMassHere =logMassEnd
            call stencil(logMassEnd,base,weights)
            meanHere    =valueSupported(countMean    (:,iCase),logMassNear(1:countNear_,iCase),countMeanNear    (1:countNear_,iCase),indexFirst_,logMassEnd,base,weights)
            varianceHere=valueSupported(countVariance(:,iCase),logMassNear(1:countNear_,iCase),countVarianceNear(1:countNear_,iCase),indexFirst_,logMassEnd,base,weights)
         end if
         densityHere=interpolateAt(densityPopulation(:,k),logMassHere)/stepLogMass
         if (.not.first) then
            aliveMean    (k,iCase)=aliveMean    (k,iCase)+0.5d0*(logMassHere-logMassPrevious)*(densityPrevious*meanPrevious    +densityHere*meanHere    )
            aliveVariance(k,iCase)=aliveVariance(k,iCase)+0.5d0*(logMassHere-logMassPrevious)*(densityPrevious*variancePrevious+densityHere*varianceHere)
         end if
         first           =.false.
         logMassPrevious =logMassHere
         densityPrevious =densityHere
         meanPrevious    =meanHere
         variancePrevious=varianceHere
      end do
      return
    end subroutine aliveTarget

    subroutine levelEvaluate(k,kernel,fraction,rate,rateQuadrature,terminated)
      !!{RST
      Evaluate the kernel of the split integral (:math:`f M^\prime` at each quadrature node), the subresolution mass loss
      rate, and the split rate at :math:`\omega_k`.
      !!}
      implicit none
      integer         , intent(in   )                 :: k
      double precision, intent(  out), dimension(:  ) :: kernel
      double precision, intent(  out), dimension(0: ) :: fraction      , rate        , &
           &                                             rateQuadrature
      logical         , intent(  out), dimension(0: ) :: terminated
      double precision                                :: time          , growthFactor
      integer                                         :: i             , q

      ! The cost of each mass grows with the number of its quadrature nodes, so schedule dynamically.
      !$omp do schedule(dynamic)
      do i=0,countMass
         time        =+criticalOverdensity_     %timeOfCollapse(criticalOverdensity=omega(k),mass=mass(i),node=nodeThread)
         growthFactor=+cosmologicalMassVariance_%rootVariance  (time               =time    ,mass=mass(i)                ) &
              &       /cosmologicalMassVariance_%rootVariance  (time               =timeNow ,mass=mass(i)                )
         if (i == 0) then
            ! A branch at the resolution ends.
            fraction(i)=-1.0d0
         else
            fraction(i)=+mergerTreeBranchingProbability_%fractionSubresolution(mass(i),omega(k),time,self%massResolution,nodeThread) &
                 &      *growthFactor
         end if
         terminated(i)=fraction(i) < 0.0d0
         if (nodeCount(i) > 0 .and. .not.terminated(i)) then
            rate(i)=+mergerTreeBranchingProbability_%probability(mass(i),omega(k),time,self%massResolution,nodeThread) &
                 &  *growthFactor
            do q=nodeStart(i),nodeStart(i)+nodeCount(i)-1
               kernel(q)=+mergerTreeBranchingProbability_%rate(mass(i),omega(k),time,exp(logMassNode(q)),nodeThread) &
                    &    *growthFactor                                                                               &
                    &    *exp(logMassNode(q))
            end do
            rateQuadrature(i)=sum(weightNode(nodeStart(i):nodeStart(i)+nodeCount(i)-1)*kernel(nodeStart(i):nodeStart(i)+nodeCount(i)-1))
         else
            rate          (i)=0.0d0
            rateQuadrature(i)=0.0d0
            if (nodeCount(i) > 0) kernel(nodeStart(i):nodeStart(i)+nodeCount(i)-1)=0.0d0
         end if
      end do
      !$omp end do
      return
    end subroutine levelEvaluate

    double precision function quadratureError(rate,rateQuadrature)
      !!{RST
      Return the maximum fractional difference between the split rates from the quadrature and from the branching probability
      class.
      !!}
      implicit none
      double precision, intent(in   ), dimension(0:) :: rate, rateQuadrature

      quadratureError=maxval(abs(rateQuadrature/rate-1.0d0),mask=rate > 0.0d0)
      return
    end function quadratureError

    subroutine sources(iCase,kernel,terminated,countMean_,countVariance_,sourceMean_,sourceVariance_,sourceSquared_)
      !!{RST
      Evaluate the sources of the mean, :math:`\int f\Delta\mathrm{d}M^\prime`, and of the variance, :math:`\int
      f[\Delta^2+W(M-M^\prime)+W(M^\prime)-W(M)]\mathrm{d}M^\prime`, and :math:`\int f\Delta^2\mathrm{d}M^\prime`, for the
      given mean and variance on the mass grid.
      !!}
      implicit none
      integer         , intent(in   )                 :: iCase
      double precision, intent(in   ), dimension(:  ) :: kernel
      logical         , intent(in   ), dimension(0: ) :: terminated
      double precision, intent(in   ), dimension(0: ) :: countMean_     , countVariance_
      double precision, intent(  out), dimension(0: ) :: sourceMean_    , sourceVariance_ , &
           &                                             sourceSquared_
      double precision                                :: delta          , varianceChildren, &
           &                                             sumMean        , sumVariance     , &
           &                                             sumSquared     , sumTerminal     , &
           &                                             sumAdded       , contribution
      integer                                         :: i              , q               , &
           &                                             latticePrevious

      ! Every thread executes this procedure, so all writes to the (shared) results are made within the shared loop. The sums are
      ! accumulated in local scalars so that threads do not repeatedly write to neighboring array elements. For the decomposition
      ! of the variance, the source ∫ f Δ² dM′ is also split into the contributions of terminal splits (M′ ≥ m) and structural
      ! splits, and by mass ratio, and the count added by splits is found. The nodes of a parent fall on a contiguous range of the
      ! lattice of mass ratios (in decreasing order), so the first contribution to each point of the lattice replaces its value
      ! from the previous evaluation.
      !$omp do schedule(dynamic)
      do i=0,countMass
         sumMean    =0.0d0
         sumVariance=0.0d0
         sumSquared =0.0d0
         sumTerminal=0.0d0
         sumAdded   =0.0d0
         if (self%decomposeVariance) then
            ratioActive    (i)=.not.terminated(i) .and. nodeCount(i) > 0
            latticePrevious   =-1
         end if
         if (terminated(i)) then
            sourceMean_    (i)=0.0d0
            sourceVariance_(i)=0.0d0
            sourceSquared_ (i)=0.0d0
            if (self%decomposeVariance) then
               sourceSquaredTerminal(i)=0.0d0
               sourceAdded          (i)=0.0d0
            end if
            cycle
         end if
         do q=nodeStart(i),nodeStart(i)+nodeCount(i)-1
            if (caseStatistic(iCase) == 'C') then
               delta           =+dble(countAdded(q,caseThreshold(iCase)))                                         &
                    &           +sum(countMean_    (basePrimary  (q)-1:basePrimary  (q)+2)*weightsPrimary  (:,q)) &
                    &           -    countMean_    (i)
               varianceChildren=+sum(countVariance_(basePrimary  (q)-1:basePrimary  (q)+2)*weightsPrimary  (:,q))
            else
               delta           =+dble(countAdded(q,caseThreshold(iCase)))                                         &
                    &           +sum(countMean_    (basePrimary  (q)-1:basePrimary  (q)+2)*weightsPrimary  (:,q)) &
                    &           +sum(countMean_    (baseSecondary(q)-1:baseSecondary(q)+2)*weightsSecondary(:,q)) &
                    &           -    countMean_    (i)
               varianceChildren=+sum(countVariance_(basePrimary  (q)-1:basePrimary  (q)+2)*weightsPrimary  (:,q)) &
                    &           +sum(countVariance_(baseSecondary(q)-1:baseSecondary(q)+2)*weightsSecondary(:,q))
            end if
            sumMean    =sumMean    +weightNode(q)*kernel(q)* delta
            sumVariance=sumVariance+weightNode(q)*kernel(q)*(delta**2+varianceChildren-countVariance_(i))
            sumSquared =sumSquared +weightNode(q)*kernel(q)* delta**2
            if (self%decomposeVariance) then
               contribution=weightNode(q)*kernel(q)*delta**2
               if (latticePrevious < 0) latticePrevious=latticeNode(q)+1
               if (latticeNode(q) < latticePrevious) then
                  squaredRatio(latticeNode(q):latticePrevious-1,i,:)=0.0d0
                  latticePrevious                                   =latticeNode(q)
               end if
               if (countAdded(q,caseThreshold(iCase)) == 1) then
                  sumTerminal                     =sumTerminal                     +contribution
                  squaredRatio(latticeNode(q),i,1)=squaredRatio(latticeNode(q),i,1)+contribution
               else
                  squaredRatio(latticeNode(q),i,2)=squaredRatio(latticeNode(q),i,2)+contribution
               end if
               sumAdded=sumAdded+weightNode(q)*kernel(q)*dble(countAdded(q,caseThreshold(iCase)))
            end if
         end do
         sourceMean_    (i)=sumMean
         sourceVariance_(i)=sumVariance
         sourceSquared_ (i)=sumSquared
         if (self%decomposeVariance) then
            sourceSquaredTerminal(i)=sumTerminal
            sourceAdded          (i)=sumAdded
         end if
      end do
      !$omp end do
      return
    end subroutine sources

    double precision function characteristicStep(logMassPrevious)
      !!{RST
      Return :math:`\ln M` at :math:`\omega_k` of the characteristic whose foot at :math:`\omega_{k+1}`, found as for the points
      of the grid, is ``logMassPrevious``.
      !!}
      implicit none
      double precision, intent(in   ) :: logMassPrevious
      double precision, parameter     :: toleranceStep         =1.0d-13
      integer         , parameter     :: countIterationsMaximum=100
      double precision                :: fractionHere                  , logMassTrial, &
           &                             logMassNew
      integer                         :: iteration

      characteristicStep=logMassPrevious+interpolateFraction(fractionPrevious,terminatedPrevious,min(logMassPrevious,logMass(countMass)))*stepOmega
      do iteration=1,countIterationsMaximum
         fractionHere=interpolateFraction(fractionCurrent,terminatedCurrent,min(characteristicStep,logMass(countMass)))
         logMassTrial=characteristicStep-fractionHere*stepOmega
         if (logMassTrial > logMass(1)) then
            logMassNew=logMassPrevious+0.5d0*(fractionHere+interpolateFraction(fractionPrevious,terminatedPrevious,min(logMassTrial,logMass(countMass))))*stepOmega
         else
            logMassNew=logMassPrevious+fractionHere*stepOmega
         end if
         if (abs(logMassNew-characteristicStep) < toleranceStep) then
            characteristicStep=logMassNew
            return
         end if
         characteristicStep=logMassNew
      end do
      call Error_Report('the step along a characteristic failed to converge'//{introspection:location})
      return
    end function characteristicStep

    integer function countActive(logMassNear_)
      !!{RST
      Return the number of leading nodes near the boundary which lie within the grid.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:) :: logMassNear_

      countActive=0
      do while (countActive < size(logMassNear_))
         if (logMassNear_(countActive+1) > logMass(countMass)) exit
         countActive=countActive+1
      end do
      return
    end function countActive

    logical function outsideSupport(i,iCase)
      !!{RST
      Return true if point ``i`` of the grid lies outside of the support of statistic A for case ``iCase``.
      !!}
      implicit none
      integer, intent(in   ) :: i, iCase

      outsideSupport=terminatedCurrent(i) .or. countNearActive(iCase) == 0
      if (.not.outsideSupport) outsideSupport=logMass(i) < logMassNear(1,iCase)
      return
    end function outsideSupport

    integer function indexFirst(logMassNearLast)
      !!{RST
      Return the index of the first point of the grid used, with the nodes near the boundary (the last of which is at
      ``logMassNearLast``), to interpolate within the support of statistic A: the first at least half of a step of the grid
      above the last node (closer points would make the interpolation ill-conditioned).
      !!}
      implicit none
      double precision, intent(in   ) :: logMassNearLast

      indexFirst=ceiling(min(dble(countMass+1),max(0.0d0,(logMassNearLast+0.5d0*stepLogMass-logMass(0))/stepLogMass)))
      do while (indexFirst > 0)
         if (logMass(indexFirst-1) < logMassNearLast+0.5d0*stepLogMass) exit
         indexFirst=indexFirst-1
      end do
      do while (indexFirst <= countMass)
         if (logMass(indexFirst  ) >= logMassNearLast+0.5d0*stepLogMass) exit
         indexFirst=indexFirst+1
      end do
      return
    end function indexFirst

    double precision function valueSupported(values,logMassNear_,valuesNear_,indexFirst_,logMass_,base,weights)
      !!{RST
      Interpolate ``values`` on the mass grid to ``logMass_``, for statistic A, whose support starts at the first of the nodes
      near its boundary, at ``logMassNear_`` with values ``valuesNear_``: zero below the boundary, from the stencil ``base``,
      ``weights`` where all of its points lie at or above ``indexFirst_``, and otherwise by cubic interpolation in the nodes near
      the boundary followed by the points of the grid from ``indexFirst_``.
      !!}
      implicit none
      double precision, intent(in   ), dimension(0:                  ) :: values
      double precision, intent(in   ), dimension(:                   ) :: logMassNear_, valuesNear_
      double precision, intent(in   )                                  :: logMass_
      integer         , intent(in   )                                  :: indexFirst_ , base
      double precision, intent(in   ), dimension(4                   ) :: weights
      double precision               , dimension(size(logMassNear_)+3) :: nodes       , nodeValues
      double precision                                                 :: weight
      integer                                                          :: countNodes_ , j         , &
           &                                                              l           , lower     , &
           &                                                              upper       , middle

      if (size(logMassNear_) == 0) then
         valueSupported=0.0d0
      else if (logMass_ < logMassNear_(1)) then
         valueSupported=0.0d0
      else if (base-1 >= indexFirst_) then
         valueSupported=sum(values(base-1:base+2)*weights)
      else
         countNodes_                      =size(logMassNear_)
         nodes      (1:countNodes_)=logMassNear_
         nodeValues (1:countNodes_)=valuesNear_
         do j=indexFirst_,min(countMass,indexFirst_+2)
            countNodes_             =countNodes_+1
            nodes      (countNodes_)=logMass(j)
            nodeValues (countNodes_)=values (j)
         end do
         ! Find the interval containing the point, and the four nodes about it.
         lower=1
         upper=countNodes_
         do while (upper-lower > 1)
            middle=(lower+upper)/2
            if (nodes(middle) <= logMass_) then
               lower=middle
            else
               upper=middle
            end if
         end do
         lower=max(1,min(countNodes_-3,lower-1))
         upper=min(countNodes_,lower+3)
         valueSupported=0.0d0
         do j=lower,upper
            weight=1.0d0
            do l=lower,upper
               if (l /= j) weight=weight*(logMass_-nodes(l))/(nodes(j)-nodes(l))
            end do
            valueSupported=valueSupported+weight*nodeValues(j)
         end do
      end if
      return
    end function valueSupported

    subroutine levelEvaluateNear(k,iCase)
      !!{RST
      Construct the quadrature nodes in :math:`\ln M^\prime` for each node near the boundary of the support of statistic A for
      case ``iCase`` (the points of the grid below :math:`M/2`, and :math:`M/2`), and evaluate the kernel there at
      :math:`\omega_k`. The level, ``k``, is an argument, as the variable of the host is private to each thread.
      !!}
      implicit none
      integer         , intent(in   ) :: k    , iCase
      double precision                :: time , growthFactor, &
           &                             mass_
      integer                         :: j    , l           , &
           &                             q

      !$omp do schedule(dynamic)
      do j=1,countNearActive(iCase)
         mass_            =exp(logMassNear(j,iCase))
         countNodesNear(j)=0
         time        =+criticalOverdensity_     %timeOfCollapse(criticalOverdensity=omega(k),mass=mass_,node=nodeThread)
         growthFactor=+cosmologicalMassVariance_%rootVariance  (time               =time    ,mass=mass_                ) &
              &       /cosmologicalMassVariance_%rootVariance  (time               =timeNow ,mass=mass_                )
         ! Branches near the boundary are not allowed to terminate (the branching classes terminate branches only at or below
         ! the resolution, which these nodes never reach).
         if (mergerTreeBranchingProbability_%fractionSubresolution(mass_,omega(k),time,self%massResolution,nodeThread) < 0.0d0) &
              & call Error_Report('a branch near the boundary of the support of statistic A terminates - this is not supported'//{introspection:location})
         if (mass_ <= 2.0d0*self%massResolution*(1.0d0+toleranceMass)) cycle
         do l=0,countMass
            if (logMass(l) >= logMassNear(j,iCase)-log(2.0d0)-toleranceMass) exit
            countNodesNear (  j)=countNodesNear(j)+1
            logMassNodeNear(countNodesNear(j),j)=logMass(l)
         end do
         countNodesNear (  j)=countNodesNear(j)+1
         logMassNodeNear(countNodesNear(j),j)=logMassNear(j,iCase)-log(2.0d0)
         do q=1,countNodesNear(j)
            logMassPrimaryNodeNear(q,j)=logMassNear(j,iCase)+log(1.0d0-exp(logMassNodeNear(q,j)-logMassNear(j,iCase)))
            kernelNear            (q,j)=+mergerTreeBranchingProbability_%rate(mass_,omega(k),time,exp(logMassNodeNear(q,j)),nodeThread) &
                 &                      *growthFactor                                                                                   &
                 &                      *exp(logMassNodeNear(q,j))
         end do
      end do
      !$omp end do
      return
    end subroutine levelEvaluateNear

    function integrandTarget(kernel_,countSecondary,varianceSecondary,countPrimary,variancePrimary,countParent,varianceParent) result(integrand)
      !!{RST
      Return the integrands of the sources of the mean, the variance, and :math:`\int f\Delta^2\mathrm{d}M^\prime`, for statistic
      A, given the kernel and the means and variances of the progenitors and the parent.
      !!}
      implicit none
      double precision, dimension(3)  :: integrand
      double precision, intent(in   ) :: kernel_          , countSecondary, &
           &                             varianceSecondary, countPrimary  , &
           &                             variancePrimary  , countParent   , &
           &                             varianceParent
      double precision                :: delta

      delta    =+countPrimary+countSecondary-countParent
      integrand=+kernel_                                                     &
           &    *[                                                           &
           &      delta                                                    , &
           &      delta**2+variancePrimary+varianceSecondary-varianceParent, &
           &      delta**2                                                   &
           &     ]
      return
    end function integrandTarget

    function quadratureTarget(logMassParent,countParent,varianceParent,logMassNode_,kernel_,logMassPrimaryNode_,countMean_,countVariance_,logMassNear_,countMeanNear_,countVarianceNear_,indexFirst_,baseSecondary_,weightsSecondary_,basePrimary_,weightsPrimary_,logMassThreshold_,squaredRatio_) result(sums)
      !!{RST
      Return the sources of the mean, the variance, and :math:`\int f\Delta^2\mathrm{d}M^\prime`, for statistic A, for a parent
      at ``logMassParent``, by the trapezoid rule in :math:`\ln M^\prime`. The integrand varies rapidly where either progenitor
      is close to the boundary of the support, and steps where it is at the boundary. So to the given quadrature nodes are
      added the images of the nodes near the boundary, :math:`M^\prime=M_j` and :math:`M-M^\prime=M_j`, with the kernel
      interpolated (linearly in its logarithm) between the given nodes, and the steps, at :math:`M^\prime=M_\mathrm{c}` and
      :math:`M-M^\prime=M_\mathrm{c}`, are nodes taken twice, with the values on either side. The interpolation stencils of the
      given nodes are given for the parents on the grid, and found here otherwise. If ``logMassThreshold_`` is given, the fourth
      element of the result is the part of :math:`\int f\Delta^2\mathrm{d}M^\prime` from terminal splits (:math:`M^\prime\ge
      m`), and ``squaredRatio_`` accumulates that integral by mass ratio (on the lattice of the decomposition) for terminal and
      structural splits.
      !!}
      implicit none
      double precision, dimension(4)                                                                  :: sums
      double precision, intent(in   )                                                                 :: logMassParent         , countParent         , &
           &                                                                                             varianceParent
      double precision, intent(in   ), dimension(:  )                                                 :: logMassNode_          , kernel_             , &
           &                                                                                             logMassPrimaryNode_   , logMassNear_        , &
           &                                                                                             countMeanNear_        , countVarianceNear_
      double precision, intent(in   ), dimension(0: )                                                 :: countMean_            , countVariance_
      integer         , intent(in   )                                                                 :: indexFirst_
      integer         , intent(in   ), dimension(:  )                                      , optional :: baseSecondary_        , basePrimary_
      double precision, intent(in   ), dimension(:,:)                                      , optional :: weightsSecondary_     , weightsPrimary_
      double precision, intent(in   )                                                      , optional :: logMassThreshold_
      double precision, intent(inout), dimension(0:,:)                                     , optional :: squaredRatio_
      ! Kinds of node: given, image of a node near the boundary as the secondary, or as the primary, and a step.
      integer         , parameter                                                                     :: kindGiven           =1, kindSecondary     =2, &
           &                                                                                             kindPrimary         =3, kindStep          =4
      double precision               , dimension(size(logMassNode_)+2*size(logMassNear_)+4)           :: logMassMerged         , kernelMerged
      integer                        , dimension(size(logMassNode_)+2*size(logMassNear_)+4)           :: kindMerged            , indexMerged         , &
           &                                                                                             sideMerged
      double precision               , dimension(size(logMassNear_))                                  :: logMassImagePrimary
      integer                        , dimension(size(logMassNear_))                                  :: nearImagePrimary
      double precision               , dimension(2  )                                                 :: logMassStep
      logical                        , dimension(2  )                                                 :: stepSecondary         , stepPrimary
      double precision               , dimension(3  )                                                 :: integrandPrevious     , integrand
      double precision               , dimension(4  )                                                 :: weightsS              , weightsP
      double precision                                                                                :: logMassBoundary_      , logMassPrimary_     , &
           &                                                                                             logMassFirst          , logMassLast         , &
           &                                                                                             logMassCandidate      , countSecondary      , &
           &                                                                                             varianceSecondary     , countPrimary        , &
           &                                                                                             variancePrimary       , kernelLower         , &
           &                                                                                             kernelUpper           , fraction
      integer                                                                                         :: countGiven            , countNear_          , &
           &                                                                                             countSteps            , countMerged         , &
           &                                                                                             qGiven                , jSecondary          , &
           &                                                                                             jPrimary              , iStep               , &
           &                                                                                             kindNext              , l                   , &
           &                                                                                             baseS                 , baseP               , &
           &                                                                                             side                  , countImagesPrimary  , &
           &                                                                                             jStart                , jEnd                , &
           &                                                                                             latticePrevious       , lattice             , &
           &                                                                                             classSplit            , latticeZeroed
      logical                                                                                         :: isSecondaryStep       , isPrimaryStep
      double precision                                                                                :: widthHalf

      sums         =0.0d0
      latticeZeroed=huge(0)
      countGiven   =size(logMassNode_)
      countNear_=size(logMassNear_)
      if (countGiven < 2) return
      logMassFirst    =logMassNode_(1         )
      logMassLast     =logMassNode_(countGiven)
      logMassBoundary_=logMassNear_(1         )
      ! The steps of the integrand: where the secondary, or the primary, is at the boundary.
      countSteps=0
      if (logMassBoundary_ >= logMassFirst-toleranceMass .and. logMassBoundary_ <= logMassLast+toleranceMass) then
         countSteps               =countSteps+1
         logMassStep  (countSteps)=min(max(logMassBoundary_,logMassFirst),logMassLast)
         stepSecondary(countSteps)=.true.
         stepPrimary  (countSteps)=.false.
      end if
      if (logMassParent > logMassBoundary_) then
         logMassPrimary_=logMassParent+log(1.0d0-exp(logMassBoundary_-logMassParent))
         if (logMassPrimary_ >= logMassFirst-toleranceMass .and. logMassPrimary_ <= logMassLast+toleranceMass) then
            logMassPrimary_=min(max(logMassPrimary_,logMassFirst),logMassLast)
            if (countSteps == 1 .and. abs(logMassPrimary_-logMassBoundary_) <= toleranceMass) then
               stepPrimary(1)=.true.
            else
               countSteps               =countSteps+1
               logMassStep  (countSteps)=logMassPrimary_
               stepSecondary(countSteps)=.false.
               stepPrimary  (countSteps)=.true.
            end if
         end if
      end if
      if (countSteps == 2) then
         if (logMassStep(2) < logMassStep(1)) then
            logMassStep  =logMassStep  ([2,1])
            stepSecondary=stepSecondary([2,1])
            stepPrimary  =stepPrimary  ([2,1])
         end if
      end if
      ! Images of the nodes near the boundary (other than the boundary itself, which gives the steps) as the primary, in order of
      ! increasing M′.
      countImagesPrimary=0
      do l=countNear_,2,-1
         if (logMassParent <= logMassNear_(l)) cycle
         logMassCandidate=logMassParent+log(1.0d0-exp(logMassNear_(l)-logMassParent))
         if (logMassCandidate <= logMassFirst .or. logMassCandidate >= logMassLast) cycle
         countImagesPrimary                     =countImagesPrimary+1
         logMassImagePrimary(countImagesPrimary)=logMassCandidate
         nearImagePrimary   (countImagesPrimary)=l
      end do
      ! Merge the given nodes, the images, and the steps (each twice: the side below, then the side above) in order of M′.
      ! Images are the secondary nodes l=2,...,countNear_ within the range, and the primary images found above.
      jStart=2
      jEnd  =countNear_
      countMerged=0
      qGiven     =1
      jSecondary =jStart
      jPrimary   =1
      iStep      =1
      do while (.true.)
         ! Skip images outside of the range of M′.
         do while (jSecondary <= jEnd)
            if (logMassNear_(jSecondary) > logMassFirst .and. logMassNear_(jSecondary) < logMassLast) exit
            jSecondary=jSecondary+1
         end do
         kindNext        =0
         logMassCandidate=huge(0.0d0)
         if (qGiven     <= countGiven         ) then
            if (logMassNode_       (qGiven    ) < logMassCandidate) then
               logMassCandidate=logMassNode_       (qGiven    )
               kindNext        =kindGiven
            end if
         end if
         if (jSecondary <= jEnd               ) then
            if (logMassNear_       (jSecondary) < logMassCandidate) then
               logMassCandidate=logMassNear_       (jSecondary)
               kindNext        =kindSecondary
            end if
         end if
         if (jPrimary   <= countImagesPrimary ) then
            if (logMassImagePrimary(jPrimary  ) < logMassCandidate) then
               logMassCandidate=logMassImagePrimary(jPrimary  )
               kindNext        =kindPrimary
            end if
         end if
         ! A step is placed before any node at the same M′ (its side above follows immediately, so that the nodes at that M′
         ! contribute only intervals of zero width between the two sides).
         if (iStep      <= countSteps         ) then
            if (logMassStep        (iStep     ) <= logMassCandidate) then
               logMassCandidate=logMassStep        (iStep     )
               kindNext        =kindStep
            end if
         end if
         if (kindNext == 0) exit
         select case (kindNext)
         case (kindGiven)
            ! A given node at a step (already placed) is dropped, as the two sides of the step take its place.
            isSecondaryStep=.false.
            if (iStep > 1) isSecondaryStep=abs(logMassCandidate-logMassStep(iStep-1)) <= toleranceMass
            if (.not.isSecondaryStep) then
               countMerged               =countMerged+1
               logMassMerged(countMerged)=logMassCandidate
               kindMerged   (countMerged)=kindGiven
               indexMerged  (countMerged)=qGiven
               sideMerged   (countMerged)=0
               kernelMerged (countMerged)=kernel_(qGiven)
            end if
            qGiven=qGiven+1
         case (kindSecondary,kindPrimary,kindStep)
            ! Interpolate the kernel between the given nodes on either side.
            if      (qGiven == 1         ) then
               kernelMerged(countMerged+1)=kernel_(1         )
            else if (qGiven >  countGiven) then
               kernelMerged(countMerged+1)=kernel_(countGiven)
            else
               kernelLower=kernel_(qGiven-1)
               kernelUpper=kernel_(qGiven  )
               if (logMassNode_(qGiven) > logMassNode_(qGiven-1)) then
                  fraction=(logMassCandidate-logMassNode_(qGiven-1))/(logMassNode_(qGiven)-logMassNode_(qGiven-1))
               else
                  fraction=0.0d0
               end if
               if (kernelLower > 0.0d0 .and. kernelUpper > 0.0d0) then
                  kernelMerged(countMerged+1)=exp(log(kernelLower)+fraction*(log(kernelUpper)-log(kernelLower)))
               else
                  kernelMerged(countMerged+1)=kernelLower+fraction*(kernelUpper-kernelLower)
               end if
            end if
            if (kindNext == kindStep) then
               do side=-1,+1,2
                  countMerged               =countMerged+1
                  logMassMerged(countMerged)=logMassCandidate
                  kindMerged   (countMerged)=kindStep
                  indexMerged  (countMerged)=iStep
                  sideMerged   (countMerged)=side
                  kernelMerged (countMerged)=kernelMerged(countMerged-(side+1)/2)
               end do
               iStep=iStep+1
            else
               countMerged               =countMerged+1
               logMassMerged(countMerged)=logMassCandidate
               kindMerged   (countMerged)=kindNext
               sideMerged   (countMerged)=0
               if (kindNext == kindSecondary) then
                  indexMerged(countMerged)=jSecondary
                  jSecondary              =jSecondary+1
               else
                  indexMerged(countMerged)=jPrimary
                  jPrimary                =jPrimary  +1
               end if
            end if
         end select
      end do
      ! Trapezoid rule over the merged nodes.
      do l=1,countMerged
         isSecondaryStep=.false.
         isPrimaryStep  =.false.
         if (kindMerged(l) == kindStep) then
            isSecondaryStep=stepSecondary(indexMerged(l))
            isPrimaryStep  =stepPrimary  (indexMerged(l))
         end if
         ! The secondary.
         if (isSecondaryStep) then
            if (sideMerged(l) < 0) then
               countSecondary   =0.0d0
               varianceSecondary=0.0d0
            else
               countSecondary   =countMeanNear_    (1)
               varianceSecondary=countVarianceNear_(1)
            end if
         else
            if (kindMerged(l) == kindGiven .and. present(baseSecondary_)) then
               baseS   =baseSecondary_   (  indexMerged(l))
               weightsS=weightsSecondary_(:,indexMerged(l))
            else
               call stencil(logMassMerged(l),baseS,weightsS)
            end if
            countSecondary   =valueSupported(countMean_    ,logMassNear_,countMeanNear_    ,indexFirst_,logMassMerged(l),baseS,weightsS)
            varianceSecondary=valueSupported(countVariance_,logMassNear_,countVarianceNear_,indexFirst_,logMassMerged(l),baseS,weightsS)
         end if
         ! The primary.
         if (isPrimaryStep) then
            if (sideMerged(l) < 0) then
               countPrimary   =countMeanNear_    (1)
               variancePrimary=countVarianceNear_(1)
            else
               countPrimary   =0.0d0
               variancePrimary=0.0d0
            end if
         else
            if      (kindMerged(l) == kindPrimary) then
               ! The primary is a node near the boundary.
               countPrimary   =countMeanNear_    (nearImagePrimary(indexMerged(l)))
               variancePrimary=countVarianceNear_(nearImagePrimary(indexMerged(l)))
            else
               if (kindMerged(l) == kindGiven) then
                  logMassPrimary_=logMassPrimaryNode_(indexMerged(l))
               else
                  logMassPrimary_=logMassParent+log(1.0d0-exp(logMassMerged(l)-logMassParent))
               end if
               if (kindMerged(l) == kindGiven .and. present(basePrimary_)) then
                  baseP   =basePrimary_   (  indexMerged(l))
                  weightsP=weightsPrimary_(:,indexMerged(l))
               else
                  call stencil(logMassPrimary_,baseP,weightsP)
               end if
               countPrimary   =valueSupported(countMean_    ,logMassNear_,countMeanNear_    ,indexFirst_,logMassPrimary_,baseP,weightsP)
               variancePrimary=valueSupported(countVariance_,logMassNear_,countVarianceNear_,indexFirst_,logMassPrimary_,baseP,weightsP)
            end if
         end if
         integrand=integrandTarget(kernelMerged(l),countSecondary,varianceSecondary,countPrimary,variancePrimary,countParent,varianceParent)
         lattice  =max(0,min(countMass,nint((logMassParent-logMassMerged(l))/stepLogMass)))
         if (l > 1) then
            widthHalf=0.5d0*(logMassMerged(l)-logMassMerged(l-1))
            sums(1:3)=sums(1:3)+widthHalf*(integrandPrevious+integrand)
            if (present(logMassThreshold_)) then
               ! Splits are terminal if M′ ≥ m (the nodes include the threshold, so no interval straddles it).
               if (0.5d0*(logMassMerged(l)+logMassMerged(l-1)) >= logMassThreshold_) then
                  classSplit=1
                  sums(4)   =sums(4)+widthHalf*(integrandPrevious(3)+integrand(3))
               else
                  classSplit=2
               end if
               if (present(squaredRatio_)) then
                  ! The lattice points are reached in decreasing order: the first contribution to each replaces its value from the
                  ! previous evaluation.
                  if (latticeZeroed   == huge(0)      ) latticeZeroed=latticePrevious+1
                  if (latticePrevious <  latticeZeroed) then
                     squaredRatio_(latticePrevious:latticeZeroed-1,:)=0.0d0
                     latticeZeroed=latticePrevious
                  end if
                  if (lattice         <  latticeZeroed) then
                     squaredRatio_(lattice        :latticeZeroed-1,:)=0.0d0
                     latticeZeroed=lattice
                  end if
                  squaredRatio_(latticePrevious,classSplit)=squaredRatio_(latticePrevious,classSplit)+widthHalf*integrandPrevious(3)
                  squaredRatio_(lattice        ,classSplit)=squaredRatio_(lattice        ,classSplit)+widthHalf*integrand        (3)
               end if
            end if
         end if
         integrandPrevious=integrand
         latticePrevious  =lattice
      end do
      return
    end function quadratureTarget

    subroutine sourcesTarget(iCase,kernel,terminated,countMean_,countVariance_,countMeanNear_,countVarianceNear_,sourceMean_,sourceVariance_,sourceSquared_,sourceMeanNear_,sourceVarianceNear_)
      !!{RST
      As ``sources``, for statistic A, whose mean and variance vanish below the boundary of the support, and take the values
      ``countMeanNear_`` and ``countVarianceNear_`` at the nodes near it, for the points of the grid and for those nodes.
      !!}
      implicit none
      integer         , intent(in   )                 :: iCase
      double precision, intent(in   ), dimension(:  ) :: kernel
      logical         , intent(in   ), dimension(0: ) :: terminated
      double precision, intent(in   ), dimension(0: ) :: countMean_     , countVariance_
      double precision, intent(in   ), dimension(:  ) :: countMeanNear_ , countVarianceNear_
      double precision, intent(  out), dimension(0: ) :: sourceMean_    , sourceVariance_    , &
           &                                             sourceSquared_
      double precision, intent(inout), dimension(:  ) :: sourceMeanNear_, sourceVarianceNear_
      double precision               , dimension(4  ) :: sums
      integer                                         :: i              , j                  , &
           &                                             qStart         , qEnd               , &
           &                                             countNear_     , indexFirst_

      countNear_ =countNearActive(iCase)
      indexFirst_=countMass+1
      if (countNear_ > 0) indexFirst_=indexFirst(logMassNear(countNear_,iCase))
      ! Every thread executes this procedure, so all writes to the (shared) results are made within the shared loops.
      !$omp do schedule(dynamic)
      do i=0,countMass
         if (self%decomposeVariance) then
            ratioActive          (i)=.false.
            sourceSquaredTerminal(i)=0.0d0
            sourceAdded          (i)=0.0d0
         end if
         if (countNear_ == 0 .or. terminated(i) .or. nodeCount(i) == 0) then
            sourceMean_    (i)=0.0d0
            sourceVariance_(i)=0.0d0
            sourceSquared_ (i)=0.0d0
            cycle
         end if
         if (logMass(i) < logMassNear(1,iCase)) then
            sourceMean_    (i)=0.0d0
            sourceVariance_(i)=0.0d0
            sourceSquared_ (i)=0.0d0
            cycle
         end if
         qStart=nodeStart(i)
         qEnd  =nodeStart(i)+nodeCount(i)-1
         if (self%decomposeVariance) then
            sums=quadratureTarget(                                                                                               &
                 &                logMass(i),countMean_(i),countVariance_(i)                                                   , &
                 &                logMassNode(qStart:qEnd),kernel(qStart:qEnd)                                                 , &
                 &                logMassPrimaryNode(qStart:qEnd)                                                              , &
                 &                countMean_,countVariance_                                                                    , &
                 &                logMassNear(1:countNear_,iCase),countMeanNear_(1:countNear_),countVarianceNear_(1:countNear_), &
                 &                indexFirst_                                                                                  , &
                 &                baseSecondary(qStart:qEnd),weightsSecondary(:,qStart:qEnd)                                   , &
                 &                basePrimary  (qStart:qEnd),weightsPrimary  (:,qStart:qEnd)                                   , &
                 &                log(self%massThresholds(caseThreshold(iCase))),squaredRatio(:,i,:)                             &
                 &               )
            sourceSquaredTerminal(i)=sums(4)
            ratioActive          (i)=.true.
         else
            sums=quadratureTarget(                                                                                               &
                 &                logMass(i),countMean_(i),countVariance_(i)                                                   , &
                 &                logMassNode(qStart:qEnd),kernel(qStart:qEnd)                                                 , &
                 &                logMassPrimaryNode(qStart:qEnd)                                                              , &
                 &                countMean_,countVariance_                                                                    , &
                 &                logMassNear(1:countNear_,iCase),countMeanNear_(1:countNear_),countVarianceNear_(1:countNear_), &
                 &                indexFirst_                                                                                  , &
                 &                baseSecondary(qStart:qEnd),weightsSecondary(:,qStart:qEnd)                                   , &
                 &                basePrimary  (qStart:qEnd),weightsPrimary  (:,qStart:qEnd)                                     &
                 &               )
         end if
         sourceMean_    (i)=sums(1)
         sourceVariance_(i)=sums(2)
         sourceSquared_ (i)=sums(3)
      end do
      !$omp end do
      !$omp do schedule(dynamic)
      do j=1,countNear_
         if (countNodesNear(j) == 0) then
            sourceMeanNear_    (j)=0.0d0
            sourceVarianceNear_(j)=0.0d0
            cycle
         end if
         sums=quadratureTarget(                                                                                               &
              &                logMassNear(j,iCase),countMeanNear_(j),countVarianceNear_(j)                                 , &
              &                logMassNodeNear(1:countNodesNear(j),j),kernelNear(1:countNodesNear(j),j)                     , &
              &                logMassPrimaryNodeNear(1:countNodesNear(j),j)                                                , &
              &                countMean_,countVariance_                                                                    , &
              &                logMassNear(1:countNear_,iCase),countMeanNear_(1:countNear_),countVarianceNear_(1:countNear_), &
              &                indexFirst_                                                                                    &
              &               )
         sourceMeanNear_    (j)=sums(1)
         sourceVarianceNear_(j)=sums(2)
      end do
      !$omp end do
      return
    end subroutine sourcesTarget

  end subroutine mergerTreeBranchingMomentsPerform
