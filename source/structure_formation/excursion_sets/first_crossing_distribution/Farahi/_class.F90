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

!+    Contributions to this file made by: Arya Farahi, Andrew Benson, Christoph Behrens, Xiaolong Du. The pinning of the first
!+    crossing probability and rate tabulations to absolute lattices for issue #1317 was drafted with assistance from Claude, and
!+    reviewed and verified by Andrew Benson.

!!{RST
Implements a excursion set first crossing statistics class using the algorithm of :cite:t:`benson_dark_2012`.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Excursion_Sets_Barriers   , only : excursionSetBarrierClass
  use :: Numerical_Interpolation   , only : interpolator
  use :: Numerical_Ranges          , only : rangeLattice

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahi" docformat="rst">
   <description>
   An excursion set first crossing statistics class using the algorithm of :cite:t:`benson_dark_2012`. For trajectories originating from a point :math:`(S_1,\delta_1)`, the distribution of first crossings of a barrier :math:`B(S)`, :math:`f(S)`, is obtained by finding the solution to the integral equation:

   .. math::
      :label: eq-OldExcursionMethod

       1 =  \int_0^S f(\tilde{S})\mathrm{d}\tilde{S} + \int_{-\infty}^{B(S)} P(\delta,S) \mathrm{d} \delta,

   where :math:`P(\delta,S) \mathrm{d} \delta` is the probability for a trajectory to lie between :math:`\delta` and :math:`\delta + \mathrm{d} \delta` at variance :math:`S`, having originated from the point :math:`(S_1,\delta_1)` having not crossed the barrier at any smaller :math:`\tilde{S} &lt; S`. In the absence of a barrier, :math:`P(\delta,S)` would be equal to :math:`P_0(\delta,S)`. However, since some trajectories will have crossed the barrier at :math:`\tilde{S} &lt; S` we must subtract off their contribution to :math:`P_0(\delta,S)`. Writing the distribution of :math:`\delta` at :math:`S` for trajectories originating at some :math:`(\tilde{\delta},\tilde{S})` as :math:`P^\prime(\delta|S,\tilde{\delta},\tilde{S})`, we can therefore write

   .. math::

      P(\delta,S) = P_0(\delta,S) - \int_{0}^{S} f(\tilde{S}) P^\prime(\delta|S,\tilde{\delta},\tilde{S}) \mathrm{d}\tilde{S},

   For an interesting class of cases, both :math:`P_0(\delta,S)` and :math:`P^\prime(\delta|S,\tilde{\delta},\tilde{S})` are normal distributions, and we can write

   .. math::

      P_0(\delta,S) = \frac{1}{\sqrt{2 \pi S}} \exp\left\{-{\Delta \delta^2 [\delta,\delta_1,S,S_1] \over 2 \Delta S [S,S_1]}\right\},

   and

   .. math::

      P^\prime(\delta|S,\tilde{\delta},\tilde{S}) = \frac{1}{\sqrt{2 \pi \Delta S[S,\tilde{S}]}} \exp\left\{-{\Delta \delta^2 [\delta,\tilde{\delta},S,\tilde{S}] \over 2 \Delta S [S,\tilde{S}]}\right\},

   where we refer to :math:`\Delta \delta[\delta,\tilde{\delta},S,\tilde{S}]` as the "effective offset", and to :math:`\Delta S[S,\tilde{S}] = \mathrm{Var}(S) - \mathrm{Cov}(S,\tilde{S})` as the "residual variance". Note that :math:`\mathrm{Cov}(S,S_1) = 0` (since all trajectories pass through :math:`\delta_1` at :math:`S_1`), and so :math:`\Delta S[S,S_1] = \mathrm{Var}(S)`.

   For a standard Weiner process (such as applies to the standard case considered in excursion set theory, namely uncorrelated and unconstrained steps), we have trivially that

   .. math::

      \Delta \delta[\delta,\tilde{\delta},S,\tilde{S}] = \delta - \tilde{\delta},

   and

   .. math::

      \Delta S[S,\tilde{S}] = S - \tilde{S},

   since the Weiner process is invariant under translations of the starting point.

   Using the above results, we can rewrite eqn. (:eq:`eq-OldExcursionMethod`):

   .. math::

      1 = \int_0^S f(\tilde{S})\mathrm{d}\tilde{S} + \int_{-\infty}^{B(S)} \left[ P_0(\delta,S) - \int_{0}^{S} f(\tilde{S})
      P^\prime(\delta|S,\tilde{\delta},\tilde{S}) \mathrm{d} \delta \right] ,

   in general and, for the case of a Gaussian distribution:

   .. math::

      1 = \int_0^S f(\tilde{S})\mathrm{d}\tilde{S} + \int_{-\infty}^{B(S)} \left[ \frac{1}{\sqrt{2 \pi \Delta S[S,S_1]}}
      \exp\left(-\frac{\Delta \delta^2[\delta,\delta_1,S,S_1]}{2 \Delta S[S,S_1]}\right) - \int_{0}^{S} f(\tilde{S}) \frac{1}{\sqrt{2 \pi \Delta S [S,\tilde{S}]}}
      \exp\left(-\frac{\Delta \delta^2[\delta,B(\tilde{S}),S,\tilde{S}]}{2 \Delta S [S,\tilde{S}]}\right)\mathrm{d}\tilde{S} \right] \mathrm{d} \delta .

   The integral over :math:`\mathrm{d}\delta` can be carried out analytically to give:

   .. math::
      :label: eq-NewExcursionMethod

       1 = \int_0^S f(\tilde{S})\mathrm{d}\tilde{S}+ \hbox{erf}\left\{\frac{\Delta \delta [B(S),\delta_1,S,S_1]}{\sqrt{2\Delta S[S,S_1]}}\right\} - \int_{0}^{S} f(\tilde{S})
       \hbox{erf}\left\{\frac{\Delta \delta [B(S),B(\tilde{S}),S,\tilde{S}]}{\sqrt{2 \Delta S [S,\tilde{S}]}}\right\} \mathrm{d}S^{\prime\prime}.

   We now discretize eqn. (:eq:`eq-NewExcursionMethod`). Specifically, we divide the :math:`S` space into :math:`N` intervals defined by the points:

   .. math::

      S_i = \left\{ \begin{array}{ll}
                     0 &amp; \hbox{if } i=0 \\
                     \sum_{j=0}^{i-1} \Delta S_j &amp; \hbox{if } i &gt; 1.
                    \end{array}
            \right.

   Note that :math:`f(0)=0` by definition, so :math:`f(S_0)=0` always. We choose :math:`\Delta S_i = S_\mathrm{max}/N` (i.e. uniform spacing in :math:`S`) when computing first crossing distributions, and :math:`\Delta S_i \propto S_i` (i.e. uniform spacing in :math:`\log(S)`) when computing first crossing rates.

   Discretizing the integrals in eqn. (:eq:`eq-NewExcursionMethod`) gives:

   .. math::
      :label: eq-Des1

      \int_0^{S_i} f(\tilde{S})\d \tilde{S} = \sum_{j=0}^{i-1} \frac{f(S_j) + f(S_{j+1})}{2} \Delta S_j

   and:

   .. math::
      :label: eq-Des2

      \int_{0}^{S_i} f(\tilde{S}) \hbox{erf}\left\{\frac{\Delta \delta [B(S),B(\tilde{S}),S,\tilde{S}]}{\sqrt{2 \Delta S[S,\tilde{S}]}}\right\} \d \tilde{S} =
      \sum_{j=0}^{i-1} \frac{1}{2} \left(f(S_j) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_j), S_i, S_j]}{\sqrt{2 \Delta S[S_i,S_j]}}\right\} + f(S_{j+1})
      \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_{j+1}),S_i,S_{j+1}]}{\sqrt{2 \Delta S[S_i,S_{j+1}]}}\right\} \right) \Delta S_j.

   We can now rewrite eqn. (:eq:`eq-NewExcursionMethod`) in discretized form:

   .. math::
      :label: eq-DesFinal1

      1 = \sum_{j=0}^{i-1} \frac{f(S_j) + f(S_{j+1})}{2} \Delta S_j + \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2 \Delta S[S_i,S_1]}}\right\} -
      \frac{1}{2} \sum_{j=0}^{i-1} \left( f(S_j) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_j),S_i,S_j]}{\sqrt{2 \Delta S[S_i,S_j]}}\right\} + f(S_{j+1})
      \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_{j+1}),S_i,S_{j+1}]}{\sqrt{2 \Delta S[S_i,S_{j+1}]}}\right\} \right) \Delta S_j.

   Solving eqn. (:eq:`eq-DesFinal1`) for :math:`f(S_i)`:

   .. math::
      :label: eq-DesFinal11

      \left( \frac{1}{2} - \frac{1}{2} \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_i), S_i, S_i]}{\sqrt{2 \Delta S[S_i,S_i]}}\right\} \right) \Delta S_{i-1}
       f(S_i) &amp; = 1 - \sum_{j=0}^{i-2} \frac{f(S_j) + f(S_{j+1})}{2} \Delta S_j - \frac{f(S_{i-1})}{2} \Delta S_{i-1} -
       \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2 \Delta S[S_i,S_1]}}\right\} \nonumber \\
      &amp;  + \frac{1}{2} \sum_{j=0}^{i-2} \left( f(S_j) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_j),S_i,S_j]}{\sqrt{2 \Delta S [S_i,S_j]}}\right\} +
      f(S_{j+1}) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_{j+1}),S_i,S_{j+1}]}{\sqrt{2 \Delta S[S_i,S_{j+1}]}}\right\} \right)\Delta S_j \nonumber \\
      &amp;  + \frac{1}{2} f(S_{i-1}) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_{i-1}),S_i,S_{i-1}]}{\sqrt{2 \Delta S [S_i,S_{i-1}]}}\right\} \Delta S_{i-1}.

   For all barriers that we consider:

   .. math::

      \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_i),S_i,S_i]}{\sqrt{2 \Delta S[S_i,S_i]}}\right\} = 0.

   We can then simplify eqn. (:eq:`eq-DesFinal11`):

   .. math::
      :label: eq-DesFinal2

      f(S_i) &amp; = {2 \over \Delta S_{i-1}}\left[1 - \sum_{j=0}^{i-2} \frac{f(S_j) + f(S_{j+1})}{2} \Delta S_j -
         \frac{f(S_{i-1})}{2} \Delta S_{i-1} - \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2 \Delta S [S_i,S_1] }}\right\} \right.  \nonumber \\
      &amp;  + \frac{1}{2} \sum_{j=0}^{i-2} \left( f(S_j) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_j),S_i,S_j]}{\sqrt{2 \Delta S [S_i,S_j]}}\right\} +
      f(S_{j+1}) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_{j+1}),S_i,S_{j+1}]}{\sqrt{2 \Delta S [S_i,S_{j+1}]}}\right\} \right)\Delta S_j \nonumber \\
      &amp;  \left. + \frac{1}{2} f(S_{i-1}) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_{i-1}),S_i,S_{i-1}]}{\sqrt{2 \Delta S [S_i,S_{i-1}]}}\right\} \Delta
       S_{i-1}\right].

   Consolidating terms in the summations:

   .. math::
      :label: eq-DesFinal2a

        f(S_i) = {2 \over \Delta S_{i-1}}\left[1 - \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2\Delta S[S_i,S_1]}}\right\} - \sum_{j=0}^{i-1}
        \left( 1-\hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_j),S_i,S_j]}{\sqrt{2 \Delta S [S_i,S_j]}}\right\} \right) f(S_j) {\Delta S_{j-1} + \Delta S_j
        \over 2} \right].

   In the case of constant :math:`\Delta S_j(=\Delta S)` this can be simplified further:

   .. math::
      :label: eq-DesFinal3

        f(S_i) = {2 \over \Delta S}\left[1 - \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2\Delta S [S_i,S_1]}}\right\}\right] - 2 \sum_{j=0}^{i-1}
        \left(1- \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_j),S_i,S_j]}{\sqrt{2 \Delta S[S_i,S_j]}}\right\} \right) f(S_j).

   In either case (i.e. eqns. :eq:`eq-DesFinal2a` and :eq:`eq-DesFinal3`) solution proceeds recursively: :math:`f(S_0)=0` by definition, :math:`f(S_1)` depends only on the known barrier and :math:`f(S_0)`, :math:`f(S_i)` depends only on the known barrier and :math:`f(S_{&lt;i})`.

   The first crossing rate is computed using the same method but with an effective barrier which is offset by the position of the progenitor in the :math:`(\delta,S)` plane, plus a small shift in time. The non-crossing rate, :math:`g(S_\mathrm{max})`, defined as the rate at which trajectories reach the maximum possible variance, :math:`S_\mathrm{max}`, without ever crossing the barrier---is computed directly by integrating over the first crossing rate distribution, i.e. :math:`g(S_\mathrm{max}) = 1 -\int_0^{S_\mathrm{max}} f(\tilde{S}) \mathrm{d}\tilde{S}`. Note that since the numerical integration occurs only up to a finite maximum :math:`S_\mathrm{max}`, a non-zero non-crossing rate will be computed for CDM-like barriers even though in reality they should have zero non-crossing rate. As such, use of this method for such barriers is not recommended.
   </description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingClass) :: excursionSetFirstCrossingFarahi
     !!{RST
     An excursion set first crossing statistics class using the algorithm of :cite:t:`benson_dark_2012`.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer                       :: cosmologyFunctions_              => null()
     class           (excursionSetBarrierClass     ), pointer                       :: excursionSetBarrier_             => null()
     class           (cosmologicalMassVarianceClass), pointer                       :: cosmologicalMassVariance_        => null()
     ! Variables used in tabulation the first crossing function.
     double precision                                                               :: timeMaximum                                , timeMinimum                               , &
          &                                                                            varianceMaximum
     integer                                                                        :: countTime                                  , countVariance
     ! Lattices to which the probability tabulation is pinned. These are the source of truth for its extent: the limits and
     ! counts above are kept in step with them, since they are read in many places.
     type            (rangeLattice                 )                                :: latticeVariance                            , latticeTime
     double precision                               , allocatable, dimension(:,:)   :: firstCrossingProbability
     double precision                               , allocatable, dimension(:  )   :: time                                       , variance
     double precision                                                               :: varianceStep
     logical                                                                        :: tableInitialized                 =  .false., fileNameInitialized                       , &
          &                                                                            varianceIsUnlimited
     type            (interpolator                 ), allocatable                   :: interpolatorTime                           , interpolatorVariance
     ! Variables used in tabulation the first crossing rate function.
     double precision                                                               :: timeMaximumRate                            , timeMinimumRate                           , &
          &                                                                            varianceMaximumRate                        , massMinimumRateNonCrossing
     integer                                                                        :: countTimeRate                              , countVarianceProgenitorRate               , &
          &                                                                            countVarianceCurrentRate                   , countVarianceCurrentRateNonCrossing
     ! Lattices to which the rate tabulation is pinned. As for the probability tabulation these are the source of truth for its
     ! extent. Only two of its four axes lie on a lattice: `varianceCurrentRate`, which is linear in variance, and `timeRate`,
     ! which is logarithmic in time. The remaining two axes are built by `varianceRange` with a spacing which varies smoothly
     ! from logarithmic to linear, so that every interior point moves when either bound moves and no lattice can represent
     ! them. They are made deterministic instead by pinning the bounds which generate them: the upper bound is shared with
     ! `latticeVarianceCurrentRate`, while `latticeVarianceMinimumRate` exists solely to pin the lower bound.
     type            (rangeLattice                 )                                :: latticeVarianceCurrentRate                 , latticeTimeRate                           , &
          &                                                                            latticeVarianceMinimumRate
     double precision                               , allocatable, dimension(:,:,:) :: firstCrossingRate
     double precision                               , allocatable, dimension(:,:  ) :: nonCrossingRate
     double precision                               , allocatable, dimension(:    ) :: timeRate                                   , varianceProgenitorRate                    , &
          &                                                                            varianceCurrentRate                        , varianceCurrentRateNonCrossing
     logical                                                                        :: tableInitializedRate             =  .false., retabulateRateNonCrossing
     type            (interpolator                 ), allocatable                   :: interpolatorTimeRate                       , interpolatorVarianceRate                  , &
          &                                                                            interpolatorVarianceCurrentRate            , interpolatorVarianceCurrentRateNonCrossing
     ! File name used to store tabulations.
     type            (varying_string               )                                :: fileName
     logical                                                                        :: useFile
     ! Tabulation resolutions.
     integer                                                                        :: varianceNumberPerUnitProbability           , varianceNumberPerUnit                     , &
          &                                                                            timeNumberPerDecade                        , varianceNumberPerDecade                   , &
          &                                                                            varianceNumberPerDecadeNonCrossing
     ! The fractional step in time used to compute barrier crossing rates.
     double precision                                                               :: fractionalTimeStep
     ! Record of variance and time in previous call to rate functions.
     double precision                                                               :: timePreviousRate                           , variancePreviousRate
     double precision                                            , dimension(0:1)   :: hTimeRate                                  , hVarianceRate
     integer         (c_size_t                     )                                :: iTimeRate                                  , iVarianceRate
   contains
     !![
     <methods docformat="rst">
       <method description="Interpolate in the tabulated excursion set barrier crossing rates."                                            method="rateInterpolate"           />
       <method description="Interpolate in the tabulated excursion set barrier non-crossing rates."                                        method="rateNonCrossingInterpolate"/>
       <method description="Tabulate excursion set barrier crossing rates ensuring that they span the given progenitor variance and time." method="rateTabulate"              />
       <method description="Build a range of variances at which to tabulate the excursion set solutions."                                  method="varianceRange"             />
       <method description="Return the maximum variance to which to tabulate."                                                             method="varianceLimit"             />
       <method description="Return a hard upper limit on the variance to which to tabulate, or `huge(0)` if there is none."                method="varianceLimitHard"         />
       <method description="Return the lattice on which the variance of the current halo is tabulated for barrier crossing rates."         method="rateVarianceLattice"       />
       <method description="Compute the residual variance between two points."                                                             method="varianceResidual"          />
       <method description="Compute the effective offset between two points."                                                              method="offsetEffective"           />
       <method description="Read excursion set solutions from file."                                                                       method="fileRead"                  />
       <method description="Write excursion set solutions to file."                                                                        method="fileWrite"                 />
       <method description="Initialize the file name for storing excursion set data."                                                      method="fileNameInitialize"        />
     </methods>
     !!]
     final     ::                               farahiDestructor
     procedure :: probability                => farahiProbability
     procedure :: rate                       => farahiRate
     procedure :: rateNonCrossing            => farahiRateNonCrossing
     procedure :: rateInterpolate            => farahiRateInterpolate
     procedure :: rateNonCrossingInterpolate => farahiRateNonCrossingInterpolate
     procedure :: rateTabulate               => farahiRateTabulate
     procedure :: fileRead                   => farahiFileRead
     procedure :: fileWrite                  => farahiFileWrite
     procedure :: fileNameInitialize         => farahiFileNameInitialize
     procedure :: varianceRange              => farahiVarianceRange
     procedure :: varianceLimit              => farahiVarianceLimit
     procedure :: varianceLimitHard          => farahiVarianceLimitHard
     procedure :: rateVarianceLattice        => farahiRateVarianceLattice
     procedure :: varianceResidual           => farahiVarianceResidual
     procedure :: offsetEffective            => farahiOffsetEffective
  end type excursionSetFirstCrossingFarahi

  interface excursionSetFirstCrossingFarahi
     !!{RST
     Constructors for the Farahi excursion set barrier class.
     !!}
     module procedure farahiConstructorParameters
     module procedure farahiConstructorInternal
  end interface excursionSetFirstCrossingFarahi

  ! Parameters controlling tabulation range
  double precision, parameter :: redshiftMaximum       =30.0d0, redshiftMinimum=0.0d0
  ! The number of anchor points per unit of variance to which the pinned variance axes are aligned. The solution for the first
  ! crossing distribution is a recursion over variance, so its cost grows as the square of the number of points tabulated, and
  ! the variance requested is often of order unity - anchoring to whole units could therefore double the work for no gain in
  ! accuracy. Tenths of a unit bound the overshoot to ten percent of a unit while leaving the set of possible ranges coarse.
  integer         , parameter :: varianceAnchorsPerUnit=10
  ! The lattice on which the lower bound of the transition-spaced variance axes of the rate tabulation is pinned. This lattice
  ! carries no tabulated points - it exists only to make that one bound reproducible - so its density is simply its anchor
  ! interval. Tenths of a decade are chosen because the bound is not a request, but a quantity derived from the barrier and from
  ! σ(M): the tabulation of the latter is not itself pinned, so the bound can move in its final bits from run to run, and
  ! anchoring is what keeps such a movement from shifting every point of both axes. Whole decades would do that too, but the
  ! progenitor variance axis is tabulated at `varianceNumberPerDecade` points per decade - four hundred by default - and its
  ! solution is a recursion whose cost grows as the square of its length, so lowering the bound by a whole decade is expensive.
  integer         , parameter :: varianceMinimumAnchorsPerDecade=10

contains

  function farahiConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the Farahi excursion set class first crossing class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (excursionSetFirstCrossingFarahi)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class           (excursionSetBarrierClass       ), pointer       :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass  ), pointer       :: cosmologicalMassVariance_
    double precision                                                 :: fractionalTimeStep
    integer                                                          :: varianceNumberPerUnitProbability  , varianceNumberPerUnit  , &
         &                                                              timeNumberPerDecade               , varianceNumberPerDecade, &
         &                                                              varianceNumberPerDecadeNonCrossing
    type            (varying_string                 )                :: fileName
    logical                                                          :: varianceIsUnlimited

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <defaultValue>var_str('none')</defaultValue>
      <source>parameters</source>
      <description>
      The name of the file to/from which tabulations of barrier first crossing probabilities should be written/read. If set to "``none``" tables will not be stored.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>fractionalTimeStep</name>
      <defaultValue>0.01d0</defaultValue>
      <source>parameters</source>
      <description>
      The fractional time step used when computing barrier crossing rates (i.e. the step used in finite difference calculations).
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>varianceNumberPerUnitProbability</name>
      <defaultValue>1000</defaultValue>
      <source>parameters</source>
      <description>
      The number of points to tabulate per unit variance for first crossing probabilities.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>varianceNumberPerUnit</name>
      <defaultValue>40</defaultValue>
      <source>parameters</source>
      <description>
      The number of tabulation points per unit of :math:`\sigma^2` used when building the rate look-up table for the Farahi excursion-set first-crossing distribution; higher values improve interpolation accuracy at the cost of memory and initialization time.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>varianceNumberPerDecade</name>
      <defaultValue>400</defaultValue>
      <source>parameters</source>
      <description>
      The number of points to tabulate per decade of progenitor variance for first crossing rates.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>varianceNumberPerDecadeNonCrossing</name>
      <defaultValue>40</defaultValue>
      <source>parameters</source>
      <description>
      The number of points to tabulate per decade of progenitor variance for non-crossing rates.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>timeNumberPerDecade</name>
      <defaultValue>10</defaultValue>
      <source>parameters</source>
      <description>
      The number of tabulation points per decade of cosmic time used when building the first-crossing rate look-up table as a function of time; higher values improve temporal interpolation accuracy for rapidly evolving cosmologies.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>varianceIsUnlimited</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>
      If true, the variance is assumed to have no upper limit (e.g. as in the case of :term:`CDM`). This allows the tabulated solutions to be extended arbitrarily. Otherwise, tables are extended to encompass just the range of variance requested.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="excursionSetBarrier"      name="excursionSetBarrier_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=excursionSetFirstCrossingFarahi(fractionalTimeStep,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,varianceNumberPerDecadeNonCrossing,timeNumberPerDecade,varianceIsUnlimited,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="excursionSetBarrier_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function farahiConstructorParameters

  function farahiConstructorInternal(fractionalTimeStep,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,varianceNumberPerDecadeNonCrossing,timeNumberPerDecade,varianceIsUnlimited,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_) result(self)
    !!{RST
    Internal constructor for the Farahi excursion set class first crossing class.
    !!}
    implicit none
    type            (excursionSetFirstCrossingFarahi)                        :: self
    double precision                                 , intent(in   )         :: fractionalTimeStep
    integer                                          , intent(in   )         :: varianceNumberPerUnitProbability  , varianceNumberPerUnit  , &
         &                                                                      timeNumberPerDecade               , varianceNumberPerDecade, &
         &                                                                      varianceNumberPerDecadeNonCrossing
    logical                                          , intent(in   )         :: varianceIsUnlimited
    type            (varying_string                 ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetBarrierClass       ), intent(in   ), target :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="fractionalTimeStep, fileName, varianceNumberPerUnitProbability, varianceNumberPerUnit, varianceNumberPerDecade, varianceNumberPerDecadeNonCrossing, timeNumberPerDecade, varianceIsUnlimited, *cosmologyFunctions_, *excursionSetBarrier_, *cosmologicalMassVariance_"/>
    !!]

    self%tableInitialized          =.false.
    self%tableInitializedRate      =.false.
    self%retabulateRateNonCrossing =.false.
    self%timeMaximum               =-huge(0.0d0)
    self%timeMinimum               =+huge(0.0d0)
    self%varianceMaximum           =      0.0d0
    self%timeMaximumRate           =-huge(0.0d0)
    self%timeMinimumRate           =+huge(0.0d0)
    self%varianceMaximumRate       =-huge(0.0d0)
    self%massMinimumRateNonCrossing=      0.0d0
    self%timePreviousRate          =-huge(0.0d0)
    self%variancePreviousRate      =-huge(0.0d0)
    self%useFile                   =(self%fileName /= 'none')
    self%fileNameInitialized       =.not.self%useFile
    return
  end function farahiConstructorInternal

  subroutine farahiFileNameInitialize(self)
    use :: File_Utilities, only : Directory_Make, File_Path          , File_Name_Expand
    use :: Input_Paths   , only : inputPath     , pathTypeDataDynamic
    implicit none
    class(excursionSetFirstCrossingFarahi), intent(inout) :: self

    if (self%fileNameInitialized) return
    ! Build an automatic file name based on the descriptor for this object.
    if (self%fileName == "auto") &
         & self%fileName=inputPath(pathTypeDataDynamic)//'largeScaleStructure/excursionSets/'//self%objectType()//'_'//self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)//'.hdf5'
    ! Expand the file name.
    self%fileName=File_Name_Expand(self%fileName)
    ! Ensure directory exists.
    call Directory_Make(File_Path(self%fileName))
    self%fileNameInitialized=.true.
    return
  end subroutine farahiFileNameInitialize

  subroutine farahiDestructor(self)
    !!{RST
    Destructor for the Farahi excursion set first crossing class.
    !!}
    implicit none
    type(excursionSetFirstCrossingFarahi), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%excursionSetBarrier_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    return
  end subroutine farahiDestructor

  double precision function farahiProbability(self,variance,time,node)
    !!{RST
    Return the excursion set barrier at the given variance and time.
    !!}
    use :: Display         , only : displayCounter              , displayCounterClear  , displayIndent       , displayMessage, &
          &                         displayUnindent             , verbosityLevelWorking
    use :: Error_Functions , only : Error_Function_Complementary
    use :: File_Utilities  , only : File_Lock                   , File_Unlock          , lockDescriptor      , File_Exists
    use :: Kind_Numbers    , only : kind_dble                   , kind_quad
    use :: MPI_Utilities   , only : mpiBarrier                  , mpiSelf
    use :: Numerical_Ranges, only : Range_Pinned                , Range_Lattice_Offset , gridSchemePerUnit   , gridSchemePerDecade
    use :: Table_Labels    , only : extrapolationTypeFix
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)                 :: self
    double precision                                 , intent(in   )                 :: variance                               , time
    type            (treeNode                       ), intent(inout)                 :: node
    double precision                                 ,                dimension(0:1) :: hTime                                  , hVariance
    double precision                                 , parameter                     :: toleranceRelativeVariance       =1.0d-6
    type            (rangeLattice                   )                                :: latticeVariance                        , latticeTime
    double precision                                                                 :: timeSeed
    double precision                                 , allocatable  , dimension(:,:) :: firstCrossingProbabilityPrevious
    integer                                                                          :: offsetVariance                         , offsetTime       , &
         &                                                                              countVariancePrevious                  , countTimePrevious
    logical                                                                          :: reuseSolutions
    class           (excursionSetBarrierClass       ), pointer                       :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass  ), pointer                       :: cosmologicalMassVariance_
    double precision                                 , allocatable  , dimension( : ) :: barrier
    double precision                                                                 :: barrierTest
    logical                                                                          :: makeTable
    integer         (c_size_t                       )                                :: iTime                                  , iVariance       , &
         &                                                                              loopCount                              , loopCountTotal  , &
         &                                                                              i                                      , j               , &
         &                                                                              jTime                                  , jVariance       , &
         &                                                                              iVarianceStart
    double precision                                                                 :: probabilityCrossingPrior
    real            (kind_quad                      )                                :: offsetEffective                        , varianceResidual
    character       (len=6                          )                                :: label
    type            (varying_string                 )                                :: message
    type            (lockDescriptor                 )                                :: fileLock

    ! Note that this solver follows the convention used through Galacticus that σ(M) grows following linear theory. That is:
    !
    !  • the root-variance of the density field smoothed on a mass scale M is a function of time, σ(M,t) = σ(M,t₀) D(t)/D(t₀),
    !    where D(t) is the linear growth factor (which may also be scale-dependent);
    !  • the critical overdensity for collapse does not include a factor of the linear growth factor, i.e. δ_c ≅ 1.686 at all
    !    epochs (varying only due to the weak dependence on the epoch-dependent cosmological parameters).
    !
    ! This differs from standard treatments of the excursion set problem in which typically the root-variance, σ(M), is evaluated
    ! at z=0, and the critical overdensity for collapse is replaced with δ_c(t)/D(t). Mathematically these two approaches are
    ! equivalent, but it can be important to keep these distinctions in mind.

    ! In this function the following translations between internal variable names and math symbols are used:
    !
    !  Sᵢ                 = self%variance                (i      )
    !  B(Sᵢ)              =      barrier                 (i      )
    !  f(Sᵢ,t)            = self%firstCrossingProbability(i,iTime)
    !  Δδ[t,S₁,S₂,δ₁,δ₂] = self%offsetEffective         (self%time(iTime),0,S1,S2,0,barrier1,barrier2)
    !  ΔS[t,S₁,S₂]       = self%varianceResidual        (self%time(iTime),0,S1,S2                    )
    
    ! Read tables from file if possible.
    call self%fileNameInitialize()
    if (self%useFile .and. .not.self%tableInitialized .and. File_Exists(self%fileName)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(self%fileName,fileLock,lockIsShared=.true.)
       call self%fileRead()
       call File_Unlock(fileLock)
    end if
    ! Construct the table if necessary.
    makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+toleranceRelativeVariance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
#ifdef USEMPI
    if (self%coordinatedMPI_) call mpiBarrier()
#endif
    if (makeTable) then
       !$omp critical(farahiProbabilityTabulate)
       ! Attempt to read the file again now that we are within the critical section. If another thread made the file while we were waiting we may be able to skip building the table.
       if (self%useFile .and. File_Exists(self%fileName)) then
          call File_Lock(self%fileName,fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+toleranceRelativeVariance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
       if (makeTable) then
          ! Construct the lattices on which we will solve for the first crossing distribution, pinning both axes to absolute
          ! lattices so that the tabulation - and hence every value interpolated from it - is independent of the variance and
          ! time at which it happened to be first requested, and so that it can be extended without recomputing solutions
          ! already found.
          !
          ! The variance axis is linear and must always begin at zero, since the solution is a recursion which starts from
          ! f(S=0)=0. Zero is therefore included among the target values, so that the range is built to cover it rather than
          ! merely as a window about the requested variance, and is imposed as a hard lower limit as well, so that the safety
          ! margin cannot carry the range below it. No seed is needed: every range begins at lattice index zero, so any two of
          ! them overlap and one always contains the other. It is anchored to tenths of a unit rather than to whole units: the
          ! solution is a recursion over variance, so its cost grows as the square of the number of points, and the variance
          ! requested can be of order unity, where anchoring to a whole unit could double the work. The margin is one anchor
          ! interval, which also guarantees that the requested variance lies strictly below the final tabulated point - a point
          ! which is forced to zero as a boundary condition (see below) and so is not a solution.
          latticeVariance=Range_Pinned(                                                                                   &
               &                                     [0.0d0,variance]                                                   , &
               &                                     self%varianceNumberPerUnitProbability                              , &
               &                                     gridSchemePerUnit                                                  , &
               &                      marginOffset  =1.0d0/dble(varianceAnchorsPerUnit)                                 , &
               &                      limitMinimum  =0.0d0                                                              , &
               &                      anchorEvery   =max(1,self%varianceNumberPerUnitProbability/varianceAnchorsPerUnit), &
               &                      latticeCurrent=self%latticeVariance                                                 &
               &                     )
          ! The time axis is logarithmic, with a safety margin of a factor of two at each end - the range this class used before
          ! being pinned. It is anchored to half decades: a whole decade of cosmic time spans most of the history of the
          ! universe.
          !
          ! The axis is seeded at its upper end only, with twice the present age of the universe, which is exactly the bound the
          ! unpinned code imposed. Seeding only one end is enough to keep any two tabulations mergeable: both then reach at least
          ! that time, so both contain a common upper region and must overlap. The seed is supplied through `rangeCurrent`, which
          ! is applied after the safety margin, so it is not itself inflated by that factor. Note that this is a degenerate range
          ! - the same value at both ends - precisely so that it can raise the upper bound without ever lowering the lower one.
          ! Seeding the lower end as well, over the whole redshift range used by the rate tabulation, was measured to take this
          ! tabulation from seven epochs to thirty-one and to roughly double the cost of `tests.excursion_sets`, for no benefit
          ! in a model, which requests times spanning that range in any case.
          timeSeed    =+2.0d0*self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
          latticeTime=Range_Pinned(                                                  &
               &                                  [time]                           , &
               &                                  self%timeNumberPerDecade         , &
               &                                  gridSchemePerDecade              , &
               &                   marginFactor  =2.0d0                            , &
               &                   rangeCurrent  =[timeSeed,timeSeed]              , &
               &                   anchorEvery   =max(1,self%timeNumberPerDecade/2), &
               &                   latticeCurrent=self%latticeTime                   &
               &                  )
          ! Determine whether the solutions already found can be carried over. Each time slice is solved independently, and
          ! within a slice the solution is a recursion in which f(Sᵢ) depends only on the barrier and on f(S₍<ᵢ₎), so a solution
          ! remains valid wherever the axes are extended around it. The one exception is the final point of the variance axis,
          ! which is forced to zero as a boundary condition rather than solved; when the variance axis grows it must be solved
          ! for like any other point, so it is deliberately excluded from what is carried over.
          reuseSolutions=           self%tableInitialized                      &
               &         .and.                                                 &
               &                    self%latticeVariance         %isDefined()  &
               &         .and.                                                 &
               &                    self%latticeTime             %isDefined()  &
               &         .and.                                                 &
               &          allocated(self%firstCrossingProbability            )
          countVariancePrevious=0
          countTimePrevious    =0
          offsetVariance       =0
          offsetTime           =0
          if (reuseSolutions) then
             countVariancePrevious=self%latticeVariance%count-1
             countTimePrevious    =self%latticeTime    %count
             offsetVariance       =Range_Lattice_Offset(self%latticeVariance,latticeVariance)
             offsetTime           =Range_Lattice_Offset(self%latticeTime    ,latticeTime    )
             call move_alloc(self%firstCrossingProbability,firstCrossingProbabilityPrevious)
          end if
          if (allocated(self%variance                )) deallocate(self%variance                )
          if (allocated(self%time                    )) deallocate(self%time                    )
          if (allocated(self%firstCrossingProbability)) deallocate(self%firstCrossingProbability)
          self%latticeVariance=latticeVariance
          self%latticeTime    =latticeTime
          self%countVariance  =latticeVariance%count-1
          self%countTime      =latticeTime    %count
          self%varianceMaximum=latticeVariance%maximum()
          self%timeMinimum    =latticeTime    %minimum()
          self%timeMaximum    =latticeTime    %maximum()
          allocate(self%variance                (0:self%countVariance               ))
          allocate(self%time                    (                     self%countTime))
          allocate(self%firstCrossingProbability(0:self%countVariance,self%countTime))
          ! Take the abscissae from the lattices rather than by subdividing the ranges, so that they are bit-identical to those
          ! of any other tabulation built on the same lattices. Note that the step in variance likewise comes from the lattice:
          ! taking it as the difference of two neighbouring points is not invariant under a shift of the range.
          self%time        =latticeTime    %values()
          self%variance    =latticeVariance%values()
          self%varianceStep=latticeVariance%step  ()
          ! Loop through the table and solve for the first crossing distribution.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayIndent("solving for excursion set barrier crossing probabilities",verbosityLevelWorking)
             message="    time: "
             write (label,'(f6.3)') self%timeMinimum
             message=message//label//" to "
             write (label,'(f6.3)') self%timeMaximum
             message=message//label
             call displayMessage(message,verbosityLevelWorking)
             message="variance: "
             write (label,'(f6.3)') self%varianceMaximum
             message=message//label
             call displayMessage(message,verbosityLevelWorking)
#ifdef USEMPI
          end if
#endif
#ifdef USEMPI
          if (mpiSelf%isMaster() .and. self%coordinatedMPI_) then
             loopCountTotal=(int(self%countTime,kind=c_size_t)/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t)*(int(self%countVariance-1,kind=c_size_t)*int(self%countVariance,kind=c_size_t))/2_c_size_t
          else
#endif
             loopCountTotal= int(self%countTime,kind=c_size_t)                                               *(int(self%countVariance-1,kind=c_size_t)*int(self%countVariance,kind=c_size_t))/2_c_size_t
#ifdef USEMPI
          end if
#endif
          loopCount=0
#ifdef USEMPI
          if (self%coordinatedMPI_) self%firstCrossingProbability=0.0d0
#endif
          ! Carry over the solutions already found. The variance axis always begins at lattice index zero, so it is extended
          ! only upwards and the surviving solutions are offset along the time axis alone. The final point of the previous
          ! variance axis is excluded: it was forced to zero as a boundary condition rather than solved, and must now be solved
          ! for like any other point.
          if (reuseSolutions) then
             if (countVariancePrevious >= 1)                                                                                    &
                  & self%firstCrossingProbability        (0:countVariancePrevious-1,offsetTime+1:offsetTime+countTimePrevious)= &
                  &      firstCrossingProbabilityPrevious(0:countVariancePrevious-1,           1:           countTimePrevious)
             deallocate(firstCrossingProbabilityPrevious)
          end if
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range in which we are interested.
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMinimum,node,rateCompute=.false.)
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMaximum,node,rateCompute=.false.)
          ! Enter an OpenMP parallel region. Each thread will solve for the first crossing distribution at a different epoch.
          !$omp parallel private(iTime,i,j,iVarianceStart,probabilityCrossingPrior,excursionSetBarrier_,cosmologicalMassVariance_,barrier,offsetEffective,varianceResidual) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
          ! Create threadprivate copies of the barrier and mas variance objects.
          allocate(excursionSetBarrier_     ,mold=self%excursionSetBarrier_     )
          allocate(cosmologicalMassVariance_,mold=self%cosmologicalMassVariance_)
          !$omp critical(excursionSetsSolverFarahiDeepCopy)
          !![
          <deepCopyReset variables="self%excursionSetBarrier_ self%cosmologicalMassVariance_"/>
          <deepCopy source="self%excursionSetBarrier_"      destination="excursionSetBarrier_"     />
          <deepCopy source="self%cosmologicalMassVariance_" destination="cosmologicalMassVariance_"/>
          <deepCopyFinalize variables="excursionSetBarrier_ cosmologicalMassVariance_"/>
          !!]
          !$omp end critical(excursionSetsSolverFarahiDeepCopy)
          !$omp barrier
          allocate(barrier(0:self%countVariance))
          !$omp do schedule(dynamic)
          do iTime=1,self%countTime
#ifdef USEMPI
             if (self%coordinatedMPI_ .and. mod(iTime-1,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
             ! Construct a table of barrier values as a function of variance. The whole table is needed even where solutions are
             ! carried over, since the recursion below sums over every Sⱼ < Sᵢ.
             do i=0,self%countVariance
                barrier(i)=excursionSetBarrier_%barrier(self%variance(i),self%time(iTime),node,rateCompute=.false.)
             end do
             ! Determine where the solution for this epoch must begin. Where the whole epoch was carried over, only those points
             ! above the previous variance axis - together with its final point, which was never solved for - remain to be
             ! found; the first two points are boundary conditions which were carried over with the rest.
             if (reuseSolutions .and. iTime > int(offsetTime,kind=c_size_t) .and. iTime <= int(offsetTime+countTimePrevious,kind=c_size_t)) then
                iVarianceStart=max(2_c_size_t,int(countVariancePrevious,kind=c_size_t))
             else
                iVarianceStart=2_c_size_t
                ! Compute the first-crossing rate at the first entry in the table of variances.
                offsetEffective                       =+self%offsetEffective (self%time(iTime),0.0_kind_quad,real(self%variance(1),kind_quad),0.0_kind_quad,0.0_kind_quad,real(barrier(1),kind_quad),0.0_kind_quad,cosmologicalMassVariance_)
                varianceResidual                      =+self%varianceResidual(self%time(iTime),0.0_kind_quad,real(self%variance(1),kind_quad),0.0_kind_quad                                                       ,cosmologicalMassVariance_)
                self%firstCrossingProbability(0,iTime)=+0.0d0
                self%firstCrossingProbability(1,iTime)=+2.0d0                                                                    &
                     &                                 *real(                                                                    &
                     &                                       Error_Function_Complementary(                                       &
                     &                                                                    +offsetEffective                       &
                     &                                                                    /sqrt(2.0_kind_quad*varianceResidual)  &
                     &                                                                   )                                     , &
                     &                                        kind=kind_dble                                                     &
                     &                                       )                                                                   &
                     &                                 /self%varianceStep
             end if
             ! Iterate over variance, computing the first crossing distribution at each value.
             do i=iVarianceStart,self%countVariance
                ! Coordinate MPI processes.
#ifdef USEMPI
                if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                   call displayCounter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityLevelWorking)
#ifdef USEMPI
                end if
#endif
                !$omp atomic
                loopCount=loopCount+(i-1)
                ! Here we sum over the solution for all Sⱼ < Sᵢ to find the fraction of trajectories which crossed the barrier at
                ! S < Sᵢ at are below the barrier at Sᵢ.
                probabilityCrossingPrior=0.0d0
                do j=1,i-1
                   offsetEffective         =+self%offsetEffective (self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),real(self%variance(j),kind_quad),0.0_kind_quad,real(barrier(i),kind_quad),real(barrier(j),kind_quad),cosmologicalMassVariance_)
                   varianceResidual        =+self%varianceResidual(self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),real(self%variance(j),kind_quad)                                                                    ,cosmologicalMassVariance_)
                   probabilityCrossingPrior=+probabilityCrossingPrior                                                 &
                        &                   +self%firstCrossingProbability(j,iTime)                                   &
                        &                   *real(                                                                    &
                        &                         Error_Function_Complementary(                                       &
                        &                                                      +offsetEffective                       &
                        &                                                      /sqrt(2.0_kind_quad*varianceResidual)  &
                        &                                                     )                                     , &
                        &                         kind=kind_dble                                                      &
                        &                        )
                end do
                ! Evaluate the first crossing probability at Sᵢ.
                offsetEffective                       =self%offsetEffective (self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),0.0_kind_quad,0.0_kind_quad,real(barrier(i),kind_quad),0.0_kind_quad,cosmologicalMassVariance_)
                varianceResidual                      =self%varianceResidual(self%time(iTime),0.0_kind_quad,real(self%variance(i),kind_quad),0.0_kind_quad                                                       ,cosmologicalMassVariance_)
                self%firstCrossingProbability(i,iTime)=max(                                                                          &
                     &                                     +0.0d0,                                                                   &
                     &                                     +2.0d0                                                                    &
                     &                                     *real(                                                                    &
                     &                                           Error_Function_Complementary(                                       &
                     &                                                                        +offsetEffective                       &
                     &                                                                        /sqrt(2.0_kind_quad*varianceResidual)  &
                     &                                                                       )                                     , &
                     &                                           kind=kind_dble                                                      &
                     &                                          )                                                                    &
                     &                                     /self%varianceStep                                                        &
                     &                                     -2.0d0                                                                    &
                     &                                     *probabilityCrossingPrior                                                 &
                     &                                    )
             end do
             ! Force the probability at maximum variance to zero.
             self%firstCrossingProbability(self%countVariance,iTime)=0.0d0
          end do
          !$omp end do
          !![
          <objectDestructor name="excursionSetBarrier_"     />
          <objectDestructor name="cosmologicalMassVariance_"/>
          !!]
          deallocate(barrier)
          !$omp end parallel
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayCounterClear(verbosityLevelWorking)
             call displayUnindent("done",verbosityLevelWorking)
#ifdef USEMPI
          end if
          if (self%coordinatedMPI_) then
             call mpiBarrier()
             self%firstCrossingProbability=mpiSelf%sum(self%firstCrossingProbability)
          end if
#endif
          ! Build the interpolators.
          if (allocated(self%interpolatorVariance)) deallocate(self%interpolatorVariance)
          if (allocated(self%interpolatorTime    )) deallocate(self%interpolatorTime    )
          allocate(self%interpolatorVariance)
          allocate(self%interpolatorTime    )
          self%interpolatorVariance=interpolator(self%variance,extrapolationType=extrapolationTypeFix)
          self%interpolatorTime    =interpolator(self%time    ,extrapolationType=extrapolationTypeFix)
          ! Record that the table is now built.
          self%tableInitialized=.true.
          ! Write the table to file if possible.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             if (self%useFile) then
                call File_Lock(self%fileName,fileLock,lockIsShared=.false.)
                call self%fileWrite()
                call File_Unlock(fileLock,sync=.true.)
             end if
#ifdef USEMPI
          end if
#endif
       end if
       !$omp end critical(farahiProbabilityTabulate)
    end if
    ! Get interpolating factors.
    call self%interpolatorTime    %linearFactors(time    ,iTime    ,hTime    )
    call self%interpolatorVariance%linearFactors(variance,iVariance,hVariance)
    ! Compute first crossing probability by interpolating in the tabulated solutions.
    farahiProbability=0.0d0
    do jTime=0,1
       do jVariance=0,1
          farahiProbability=+farahiProbability                                                &
               &            +hTime                        (                            jTime) &
               &            *hVariance                    (            jVariance            ) &
               &            *self%firstCrossingProbability(iVariance-1+jVariance,iTime+jTime)
       end do
    end do
    return
  end function farahiProbability

  double precision function farahiRate(self,variance,varianceProgenitor,time,node)
    !!{RST
    Return the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)  :: self
    double precision                                 , intent(in   )  :: variance, varianceProgenitor, &
         &                                                               time
    type            (treeNode                       ), intent(inout)  :: node

    ! Note that this solver follows the convention used through Galacticus that σ(M) grows following linear theory. That is:
    !
    !  • the root-variance of the density field smoothed on a mass scale M is a function of time, σ(M,t) = σ(M,t₀) D(t)/D(t₀),
    !    where D(t) is the linear growth factor (which may also be scale-dependent);
    !  • the critical overdensity for collapse does not include a factor of the linear growth factor, i.e. δ_c ≅ 1.686 at all
    !    epochs (varying only due to the weak dependence on the epoch-dependent cosmological parameters).
    !
    ! This differs from standard treatments of the excursion set problem in which typically the root-variance, σ(M), is evaluated
    ! at z=0, and the critical overdensity for collapse is replaced with δ_c(t)/D(t). Mathematically these two approaches are
    ! equivalent, but it can be important to keep these distinctions in mind.

    if (varianceProgenitor <= variance) then
       ! For progenitor variances less than or equal to the original variance, return zero.
       farahiRate=0.0d0
    else
       ! Otherwise, interpolate in the tabulated solutions.
       farahiRate=self%rateInterpolate(variance,varianceProgenitor,time,node)
    end if
    return
  end function farahiRate

  double precision function farahiRateInterpolate(self,variance,varianceProgenitor,time,node)
    !!{RST
    Interpolate in the tabulated excursion set barrier crossing rate.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)  :: self
    double precision                                 , intent(in   )  :: variance           , varianceProgenitor, &
         &                                                               time
    type            (treeNode                       ), intent(inout)  :: node
    double precision                                 , dimension(0:1) :: hVarianceProgenitor
    integer                                                           :: jVarianceProgenitor, jTime             , &
         &                                                               jVariance
    integer         (c_size_t                       )                 :: iVarianceProgenitor

    ! Ensure that the rate is tabulated.
    call self%rateTabulate(varianceProgenitor,time,node)
    ! For progenitor variances greater than the maximum allowed variance, return zero.
    if (varianceProgenitor > self%varianceMaximumRate) then
       farahiRateInterpolate=0.0d0
       return
    end if
    ! Get interpolation in time.
    if (time /= self%timePreviousRate) then
       self%timePreviousRate    =time
       call self%interpolatorTimeRate           %linearFactors(time    ,self%iTimeRate    ,self%hTimeRate    )
    end if
    ! Get interpolation in variance.
    if (variance /= self%variancePreviousRate) then
       self%variancePreviousRate=variance
       call self%interpolatorVarianceCurrentRate%linearFactors(variance,self%iVarianceRate,self%hVarianceRate)
    end if
    ! Get interpolation in progenitor variance.
    iVarianceProgenitor=self%interpolatorVarianceRate%locate(varianceProgenitor-variance)
    ! Catch cases where the maximum variance is approached.
    if (self%varianceProgenitorRate(iVarianceProgenitor)+variance > self%varianceMaximumRate) then
       ! Force the rate to drop to zero at the maximum variance. (Necessary because we will not have a tabulated point precisely
       ! at the maximum variance.)
       hVarianceProgenitor=[                                                                                           &
            &               +1.0d0                                                                                     &
            &               -((     varianceProgenitor -variance)-self%varianceProgenitorRate(iVarianceProgenitor-1))  &
            &               /((self%varianceMaximumRate-variance)-self%varianceProgenitorRate(iVarianceProgenitor-1)), &
            &               +0.0d0                                                                                     &
            &              ]
    else
       call self%interpolatorVarianceRate%linearWeights(varianceProgenitor-variance,iVarianceProgenitor,hVarianceProgenitor)
    end if
    ! Compute first crossing probability by interpolating.
    farahiRateInterpolate=0.0d0
    do jTime=0,1
       do jVariance=0,1
          do jVarianceProgenitor=0,1
             farahiRateInterpolate=+farahiRateInterpolate                                                                                                   &
                  &                +self%hTimeRate          (                                                                                        jTime) &
                  &                *self%hVarianceRate      (                                                               jVariance                     ) &
                  &                *     hVarianceProgenitor(                      jVarianceProgenitor                                                    ) &
                  &                *self%firstCrossingRate  (iVarianceProgenitor-1+jVarianceProgenitor,self%iVarianceRate-1+jVariance,self%iTimeRate+jTime)
          end do
       end do
    end do
    return
  end function farahiRateInterpolate

  double precision function farahiRateNonCrossing(self,variance,massMinimum,time,node)
    !!{RST
    Return the rate for excursion set non-crossing.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout) :: self
    double precision                                 , intent(in   ) :: time       , variance, &
         &                                                              massMinimum
    type            (treeNode                       ), intent(inout) :: node

    ! Note that this solver follows the convention used through Galacticus that σ(M) grows following linear theory. That is:
    !
    !  • the root-variance of the density field smoothed on a mass scale M is a function of time, σ(M,t) = σ(M,t₀) D(t)/D(t₀),
    !    where D(t) is the linear growth factor (which may also be scale-dependent);
    !  • the critical overdensity for collapse does not include a factor of the linear growth factor, i.e. δ_c ≅ 1.686 at all
    !    epochs (varying only due to the weak dependence on the epoch-dependent cosmological parameters).
    !
    ! This differs from standard treatments of the excursion set problem in which typically the root-variance, σ(M), is evaluated
    ! at z=0, and the critical overdensity for collapse is replaced with δ_c(t)/D(t). Mathematically these two approaches are
    ! equivalent, but it can be important to keep these distinctions in mind.
    farahiRateNonCrossing=self%rateNonCrossingInterpolate(variance,massMinimum,time,node)
    return
  end function farahiRateNonCrossing

  double precision function farahiRateNonCrossingInterpolate(self,variance,massMinimum,time,node)
    !!{RST
    Interpolate the rate for excursion set non-crossing.
    !!}
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)  :: self
    double precision                                 , intent(in   )  :: time                           , variance , &
         &                                                               massMinimum
    type            (treeNode                       ), intent(inout)  :: node
    double precision                                 , parameter      :: toleranceRelativeMass   =2.5d-3
    double precision                                 , dimension(0:1) :: hVarianceRateNonCrossing
    double precision                                                  :: varianceMaximum
    integer         (c_size_t                       )                 :: iVarianceRateNonCrossing
    integer                                                           :: jTime                          , jVariance

    ! If the minimum mass used in computing non-crossing rates changes, retabulate non-crossing rates. This should only happen if
    ! the growth rate of σ(M) is mass dependent. We must also retabulate if non-crossing rates have no yet been tabulated.
    if (Values_Differ(self%massMinimumRateNonCrossing,massMinimum,relTol=toleranceRelativeMass)) then
       self%retabulateRateNonCrossing =.true.
       self%massMinimumRateNonCrossing=massMinimum
    end if
    ! Ensure that the rate is tabulated.
    varianceMaximum=self%cosmologicalMassVariance_%rootVariance(massMinimum,time)**2
    call self%rateTabulate(varianceMaximum,time,node)
    ! Get interpolation in time.
    if (time /= self%timePreviousRate) then
       self%timePreviousRate=time
       call self%interpolatorTimeRate%linearFactors(time,self%iTimeRate,self%hTimeRate)
    end if
    ! Get interpolation in variance.
    call self%interpolatorVarianceCurrentRateNonCrossing%linearFactors(                              &
         &                                                             max(                          &
         &                                                                 varianceMaximum-variance, &
         &                                                                 0.0d0                     &
         &                                                                )                        , &
         &                                                             iVarianceRateNonCrossing    , &
         &                                                             hVarianceRateNonCrossing      &
         &                                                            )
    ! Compute non-crossing probability by interpolating.
    farahiRateNonCrossingInterpolate=0.0d0
    do jTime=0,1
       do jVariance=0,1
          farahiRateNonCrossingInterpolate=+farahiRateNonCrossingInterpolate                                                     &
               &                           +self%hTimeRate           (                                                    jTime) &
               &                           *hVarianceRateNonCrossing (                           jVariance                     ) &
               &                           *self%nonCrossingRate     (iVarianceRateNonCrossing-1+jVariance,self%iTimeRate+jTime)
       end do
    end do
    return
  end function farahiRateNonCrossingInterpolate

  subroutine farahiRateTabulate(self,varianceProgenitor,time,node)
    !!{RST
    Tabulate the excursion set crossing rate.
    !!}
    use :: Display         , only : displayCounter              , displayCounterClear  , displayIndent       , displayMessage, &
          &                         displayUnindent             , verbosityLevelWorking
    use :: Error_Functions , only : Error_Function_Complementary
    use :: File_Utilities  , only : File_Lock                   , File_Unlock          , lockDescriptor      , File_Exists
    use :: Kind_Numbers    , only : kind_dble                   , kind_quad
    use :: MPI_Utilities   , only : mpiBarrier                  , mpiSelf
    use :: Numerical_Ranges, only : Range_Pinned                , Range_Lattice_Offset , gridSchemePerDecade
    use :: Table_Labels    , only : extrapolationTypeFix
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)                   :: self
    double precision                                 , intent(in   )                   :: time                             , varianceProgenitor
    type            (treeNode                       ), intent(inout)                   :: node
    double precision                                 , parameter                       :: varianceMinimumDefault    =1.0d-2
    double precision                                 , parameter                       :: varianceTolerance         =1.0d-6
    double precision                                 , parameter                       :: massLarge                 =1.0d16
    real            (kind=kind_quad                 ), allocatable  , dimension(:    ) :: firstCrossingRateQuad            , varianceCurrentRateQuad   , &
         &                                                                                varianceProgenitorRateQuad       , barrierRateQuad
    type            (rangeLattice                   )                                  :: latticeVarianceCurrentRate       , latticeTimeRate           , &
         &                                                                                latticeVarianceMinimumRate
    double precision                                 , allocatable  , dimension(:,:,:) :: firstCrossingRatePrevious
    double precision                                 , allocatable  , dimension(:,:  ) :: nonCrossingRatePrevious
    integer                                                                            :: offsetTimeRate                   , countTimeRatePrevious     , &
         &                                                                                countTimeRateSolved              , countTimeRateSolvedNonCrossing
    logical                                                                            :: reuseSolutions                   , reuseSolutionsNonCrossing
    double precision                                                                   :: barrierRateTest
    class           (excursionSetBarrierClass       ), pointer                         :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass  ), pointer                         :: cosmologicalMassVariance_
#ifdef USEMPI
    integer                                                                            :: taskCount
#endif
    logical                                                                            :: makeTable
    integer         (c_size_t                       )                                  :: loopCount                        , loopCountTotal
    integer                                                                            :: i                                , iTime                     , &
         &                                                                                iVariance                        , j                         , &
         &                                                                                iCompute                         , countVarianceCurrentRate
    double precision                                                                   :: timeProgenitor                   , varianceMinimumRate       , &
         &                                                                                massProgenitor                   , varianceMaximumRateLimit
    double precision                                                , dimension(2    ) :: timeSeedRate
    character       (len=64                         )                                  :: label
    type            (varying_string                 )                                  :: message
    type            (lockDescriptor                 )                                  :: fileLock
    real            (kind=kind_quad                 )                                  :: crossingFraction                 , barrierProgenitorEffective, &
         &                                                                                probabilityCrossingPrior         , varianceStepRate          , &
         &                                                                                barrier                          , growthFactorEffective     , &
         &                                                                                varianceResidual                 , offsetEffective

    ! Note that this solver follows the convention used through Galacticus that σ(M) grows following linear theory. That is:
    !
    !  • the root-variance of the density field smoothed on a mass scale M is a function of time, σ(M,t) = σ(M,t₀) D(t)/D(t₀),
    !    where D(t) is the linear growth factor (which may also be scale-dependent);
    !  • the critical overdensity for collapse does not include a factor of the linear growth factor, i.e. δ_c ≅ 1.686 at all
    !    epochs (varying only due to the weak dependence on the epoch-dependent cosmological parameters).
    !
    ! This differs from standard treatments of the excursion set problem in which typically the root-variance, σ(M), is evaluated
    ! at z=0, and the critical overdensity for collapse is replaced with δ_c(t)/D(t). Mathematically these two approaches are
    ! equivalent, but it can be important to keep these distinctions in mind.

    ! In this function the following translations between internal variable names and math symbols are used:
    !
    !   S₁                = varianceCurrent
    !   S̃                = varianceProgenitor  +varianceCurrent
    !   S                 = varianceIntermediate+varianceCurrent
    !   B(Sᵢ)              = barrier(i)
    !   f(Sᵢ,t)            = self%firstCrossingProbability(i,iTime)
    !   Δδ[t,S₁,S₂,δ₁,δ₂] = self%offsetEffective         (self%time(iTime),0,S1,S2,0,barrier1,barrier2)
    !   ΔS[t,S₁,S₂]       = self%varianceResidual        (self%time(iTime),0,S1,S2                    )

    ! Note that the variables "varianceIntermediate" and "varianceProgenitor" are defined to be the variances in excess of S₁ - which is why they
    ! appear with "varianceCurrent" added to them in the above.
    !
    ! This function is used in the calculation of the distribution of δ at some S for trajectories originating from (S₁,δ₁) and
    ! which did not cross the barrier at any intermediate variance. As such suffixes in variable names have the following
    ! meanings:
    !
    !   "Current"      - refers to the current halo being considered for branching, i.e. the halo existing at point (S₁,δ₁);
    !   "Progenitor"   - refers to the potential progenitor halo being considered, i.e. the halo corresponding to some variance S > S₁;
    !   "Intermediate" - refers to the intermediate variance, S̃ (with S₁ < S̃ < S);
    !   "Quad"         - refers to a quantity computed in quad-precision.

    ! Determine if we need to make the table.
    ! Read tables from file if possible.
    call self%fileNameInitialize()
    if (self%useFile .and. .not.self%tableInitializedRate .and. File_Exists(self%fileName)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(self%fileName,fileLock,lockIsShared=.true.)
       call self%fileRead()
       call File_Unlock(fileLock)
    end if
    makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
#ifdef USEMPI
    if (self%coordinatedMPI_) call mpiBarrier()
#endif
    if (makeTable.or.self%retabulateRateNonCrossing) then
       !$omp critical(farahiRateTabulate)
       ! Attempt to read the file again now that we are within the critical section. If another thread made the file while we were waiting we may be able to skip building the table.
       if (self%useFile .and. File_Exists(self%fileName)) then
          call File_Lock(self%fileName,fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
       if (makeTable.or.self%retabulateRateNonCrossing) then
          reuseSolutions           =.false.
          reuseSolutionsNonCrossing=.false.
          countTimeRatePrevious    =0
          offsetTimeRate           =0
          if (makeTable) then
             ! Construct the lattices on which the rate tabulation is built, pinning both of the axes which can be pinned so
             ! that the tabulation - and hence every rate interpolated from it - is independent of the variance and time at
             ! which it happened to be first requested. The variance axis of the current halo is built by `rateVarianceLattice`,
             ! which is shared with the subclasses which solve for rates.
             latticeVarianceCurrentRate=self%rateVarianceLattice(varianceProgenitor)
             self%varianceMaximumRate  =latticeVarianceCurrentRate%maximum()
             ! The time axis is logarithmic, with the factor of two safety margin at each end which this class used before being
             ! pinned, and is anchored to half decades - a whole decade of cosmic time spanning most of the history of the
             ! universe. Unlike the probability tabulation it is seeded at both ends, over the full range of redshifts, since
             ! the unpinned code seeded it that way already and so there is no additional cost to weigh.
             timeSeedRate(1)           =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
             timeSeedRate(2)           =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
             latticeTimeRate           =Range_Pinned(                                                    &
                  &                                                [time]                             ,  &
                  &                                                self%timeNumberPerDecade           ,  &
                  &                                                gridSchemePerDecade                ,  &
                  &                                 marginFactor  =2.0d0                              ,  &
                  &                                 rangeCurrent  =timeSeedRate                       ,  &
                  &                                 anchorEvery   =max(1,self%timeNumberPerDecade/2)  ,  &
                  &                                 latticeCurrent=self%latticeTimeRate                  &
                  &                                )
             self%timeMinimumRate=latticeTimeRate%minimum()
             self%timeMaximumRate=latticeTimeRate%maximum()
             ! Set the default minimum variance.
             varianceMinimumRate=varianceMinimumDefault
             ! Next reduce the variance if necessary such that the typical amplitude of fluctuations is less (by a factor of 10) than
             ! the effective barrier height at zero variance for the minimum and maximum times that we must consider. We use some
             ! suitably large mass to estimate the growth of fluctuations on large scales (since we can't assume infinitely large
             ! scales).
             allocate(excursionSetBarrier_     ,mold=self%excursionSetBarrier_     )
             allocate(cosmologicalMassVariance_,mold=self%cosmologicalMassVariance_)
             !$omp critical(excursionSetsSolverFarahiDeepCopy)
             !![
             <deepCopyReset variables="self%excursionSetBarrier_ self%cosmologicalMassVariance_"/>
             <deepCopy source="self%excursionSetBarrier_"      destination="excursionSetBarrier_"     />
             <deepCopy source="self%cosmologicalMassVariance_" destination="cosmologicalMassVariance_"/>
             <deepCopyFinalize variables="excursionSetBarrier_ cosmologicalMassVariance_"/>
             !!]
             !$omp end critical(excursionSetsSolverFarahiDeepCopy)
             growthFactorEffective          =+cosmologicalMassVariance_%rootVariance(massLarge,self%timeMaximumRate                                ) &
                  &                          /cosmologicalMassVariance_%rootVariance(massLarge,self%timeMaximumRate*(1.0d0-self%fractionalTimeStep))
             varianceMinimumRate            =min(                                                                                                                      &
                  &                              +varianceMinimumRate                                                                                                , &
                  &                              +1.0d-2                                                                                                               &
                  &                              *(                                                                                                                    &
                  &                                +excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate*(1.0d0-self%fractionalTimeStep),node,rateCompute=.true.)  &
                  &                                *dble(growthFactorEffective)                                                                                        &
                  &                                -excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate                                ,node,rateCompute=.true.)  &
                  &                               )**2                                                                                                                 &
                  &                             )
             !![
             <objectDestructor name="excursionSetBarrier_"     />
             <objectDestructor name="cosmologicalMassVariance_"/>
             !!]
             ! Pin the minimum variance. It is not a request but a quantity derived from the barrier and from σ(M) at the
             ! maximum time, so - the maximum time now being pinned - it is already deterministic in exact arithmetic. It is
             ! pinned nonetheless because the tabulation of σ(M) is not itself pinned, so the value above can move in its final
             ! bits between runs; both of the axes generated from it are spaced by a rule under which every interior point moves
             ! when either bound does, so an unpinned bound would leave those axes reproducible only to within that movement.
             ! The maximum variance is passed alongside it purely to guarantee a range of at least two points, and the safety
             ! margin is suppressed since the bound is not approached from below by anything.
             latticeVarianceMinimumRate=Range_Pinned(                                                                      &
                  &                                                [varianceMinimumRate,self%varianceMaximumRate]        , &
                  &                                                varianceMinimumAnchorsPerDecade                       , &
                  &                                                gridSchemePerDecade                                   , &
                  &                                 marginFactor  =1.0d0                                                 , &
                  &                                 anchorEvery   =1                                                     , &
                  &                                 latticeCurrent=self%latticeVarianceMinimumRate                         &
                  &                                )
             varianceMinimumRate=latticeVarianceMinimumRate%minimum()
             ! Determine whether the solutions already found can be carried over. Reuse is possible only where the tabulation is
             ! extended in time alone: the two transition-spaced variance axes are generated from the pinned bounds above, and
             ! every one of their interior points moves if either bound moves, so nothing survives a change in either. Where they
             ! do not move, each epoch is solved independently of every other, and so whole time slices carry over.
             reuseSolutions=      self%tableInitializedRate                                                                        &
                  &         .and. self%latticeTimeRate           %isDefined  ()                                                    &
                  &         .and. self%latticeVarianceCurrentRate%isDefined  ()                                                    &
                  &         .and. self%latticeVarianceMinimumRate%isDefined  ()                                                    &
                  &         .and. self%latticeVarianceCurrentRate%indexMinimum ==      latticeVarianceCurrentRate%indexMinimum      &
                  &         .and. self%latticeVarianceCurrentRate%count        ==      latticeVarianceCurrentRate%count             &
                  &         .and. self%latticeVarianceMinimumRate%indexMinimum ==      latticeVarianceMinimumRate%indexMinimum
             ! The non-crossing rates are tabulated on a different variance axis, but on the same time axis, and so carry over
             ! under the same conditions - unless the minimum mass at which they are evaluated has changed, which invalidates
             ! every one of them.
             reuseSolutionsNonCrossing=reuseSolutions .and. .not.self%retabulateRateNonCrossing
             if (reuseSolutions) then
                countTimeRatePrevious=self%latticeTimeRate%count
                offsetTimeRate       =Range_Lattice_Offset(self%latticeTimeRate,latticeTimeRate)
                                               call move_alloc(self%firstCrossingRate,firstCrossingRatePrevious)
                if (reuseSolutionsNonCrossing) call move_alloc(self%nonCrossingRate  ,nonCrossingRatePrevious  )
             end if
             if (allocated(self%varianceProgenitorRate        )) deallocate(self%varianceProgenitorRate        )
             if (allocated(self%varianceCurrentRate           )) deallocate(self%varianceCurrentRate           )
             if (allocated(self%varianceCurrentRateNonCrossing)) deallocate(self%varianceCurrentRateNonCrossing)
             if (allocated(self%timeRate                      )) deallocate(self%timeRate                      )
             if (allocated(self%firstCrossingRate             )) deallocate(self%firstCrossingRate             )
             self%latticeVarianceCurrentRate         =latticeVarianceCurrentRate
             self%latticeTimeRate                    =latticeTimeRate
             self%latticeVarianceMinimumRate         =latticeVarianceMinimumRate
             self%countTimeRate                      =latticeTimeRate           %count
             self%countVarianceCurrentRate           =latticeVarianceCurrentRate%count-1
             self%countVarianceProgenitorRate        =int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecade           ))+1
             self%countVarianceCurrentRateNonCrossing=int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecadeNonCrossing))+1
             allocate(self%varianceProgenitorRate        (0:self%countVarianceProgenitorRate                                                              ))
             allocate(self%varianceCurrentRate           (                                   0:self%countVarianceCurrentRate                              ))
             allocate(self%varianceCurrentRateNonCrossing(                                   0:self%countVarianceCurrentRateNonCrossing                   ))
             allocate(self%timeRate                      (                                                                              self%countTimeRate))
             allocate(self%firstCrossingRate             (0:self%countVarianceProgenitorRate,0:self%countVarianceCurrentRate           ,self%countTimeRate))
             ! For the variance table, the zeroth point is always zero, higher points are distributed uniformly in variance. The
             ! abscissae of the two axes which lie on lattices are taken from those lattices rather than by subdividing their
             ! ranges, so that they are bit-identical to those of any other tabulation built on the same lattices.
             self%varianceProgenitorRate        (0                                         )=0.0d0
             self%varianceProgenitorRate        (1:self%countVarianceProgenitorRate        )=self%varianceRange(varianceMinimumRate,self%varianceMaximumRate,self%countVarianceProgenitorRate        ,exponent=1.0d0)
             self%varianceCurrentRate           (0:self%countVarianceCurrentRate           )=latticeVarianceCurrentRate%values()
             self%varianceCurrentRateNonCrossing(0                                         )=0.0d0
             self%varianceCurrentRateNonCrossing(1:self%countVarianceCurrentRateNonCrossing)=self%varianceRange(varianceMinimumRate,self%varianceMaximumRate,self%countVarianceCurrentRateNonCrossing,exponent=1.0d0)
             ! The time table is logarithmically distributed in time.
             self%timeRate                                                                  =latticeTimeRate           %values()
          end if
          if (allocated(self%nonCrossingRate)) deallocate(self%nonCrossingRate)
          allocate(self%nonCrossingRate(0:self%countVarianceCurrentRateNonCrossing,self%countTimeRate))
          ! Allocate temporary arrays used in quad-precision solver for barrier crossing rates.
          allocate(varianceProgenitorRateQuad(0:self%countVarianceProgenitorRate))
          varianceProgenitorRateQuad=self%varianceProgenitorRate
          ! Loop through the table and solve for the first crossing distribution.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayIndent("solving for excursion set barrier crossing rates",verbosityLevelWorking)
             message="    time: "
             write (label,'(f6.3)') self%timeMinimumRate
             message=message//trim(label)//" to "
             write (label,'(f6.3)') self%timeMaximumRate
             message=message//trim(label)
             call displayMessage(message,verbosityLevelWorking)
             message="variance: "
             write (label,'(f9.3)') self%varianceMaximumRate
             message=message//trim(label)
             call displayMessage(message,verbosityLevelWorking)
#ifdef USEMPI
          end if
#endif
          ! Count only those epochs which must actually be solved for - those carried over from the previous tabulation are
          ! skipped below, and counting them would leave the progress report stuck short of completion.
          countTimeRateSolved           =self%countTimeRate
          countTimeRateSolvedNonCrossing=self%countTimeRate
          if (reuseSolutions           ) countTimeRateSolved           =countTimeRateSolved           -countTimeRatePrevious
          if (reuseSolutionsNonCrossing) countTimeRateSolvedNonCrossing=countTimeRateSolvedNonCrossing-countTimeRatePrevious
          loopCountTotal   =+int(countTimeRateSolvedNonCrossing,kind=c_size_t)*int(self%countVarianceCurrentRateNonCrossing+1,kind=c_size_t)
          if (makeTable) then
             loopCountTotal=+loopCountTotal                                                                                                 &
                  &         +int(countTimeRateSolved           ,kind=c_size_t)*int(self%countVarianceCurrentRate           +1,kind=c_size_t)
          end if
#ifdef USEMPI
          if (mpiSelf%isMaster() .and. self%coordinatedMPI_) then
             loopCountTotal= loopCountTotal/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t
          end if
#endif
          loopCount=0
#ifdef USEMPI
          if (self%coordinatedMPI_) then
             if (makeTable) then
                self%firstCrossingRate=0.0d0
             end if
             self%nonCrossingRate=0.0d0
          end if
          taskCount=-1
#endif
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range in which we are interested.
          barrierRateTest=self%excursionSetBarrier_%barrier(self%varianceMaximumRate,self%timeMinimumRate*(1.0d0-self%fractionalTimeStep),node,rateCompute=.true.)
          barrierRateTest=self%excursionSetBarrier_%barrier(self%varianceMaximumRate,self%timeMaximumRate                                ,node,rateCompute=.true.)
          ! Begin an OpenMP parallel region. Each parallel thread will compute first crossing rates for a different epoch.
          !$omp parallel private(iTime,timeProgenitor,iVariance,varianceStepRate,i,j,iCompute,probabilityCrossingPrior,crossingFraction,barrier,barrierProgenitorEffective,firstCrossingRateQuad,excursionSetBarrier_,cosmologicalMassVariance_,barrierRateQuad,varianceCurrentRateQuad,massProgenitor,growthFactorEffective,offsetEffective,varianceResidual,varianceMaximumRateLimit) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
          ! Create threadprivate copies of the barrier and mass variance objects.
          allocate(excursionSetBarrier_     ,mold=self%excursionSetBarrier_     )
          allocate(cosmologicalMassVariance_,mold=self%cosmologicalMassVariance_)
          !$omp critical(excursionSetsSolverFarahiDeepCopy)
          !![
          <deepCopyReset variables="self%excursionSetBarrier_ self%cosmologicalMassVariance_"/>
          <deepCopy source="self%excursionSetBarrier_"      destination="excursionSetBarrier_"     />
          <deepCopy source="self%cosmologicalMassVariance_" destination="cosmologicalMassVariance_"/>
          <deepCopyFinalize variables="excursionSetBarrier_ cosmologicalMassVariance_"/>
          !!]
          !$omp end critical(excursionSetsSolverFarahiDeepCopy)
          !$omp barrier
          allocate(barrierRateQuad(self%countVarianceProgenitorRate))
          ! In the first run, the first crossing rates are computed. In the second run, the non-crossing rates are computed on
          ! slightly different grid points.
          do iCompute=1,2
             if (iCompute == 1) then
                if (.not.makeTable) cycle
                countVarianceCurrentRate=self%countVarianceCurrentRate
                allocate(varianceCurrentRateQuad(0:self%countVarianceCurrentRate           ))
                varianceCurrentRateQuad =self%varianceCurrentRate
             else
                countVarianceCurrentRate=self%countVarianceCurrentRateNonCrossing
                if (allocated(varianceCurrentRateQuad)) deallocate(varianceCurrentRateQuad)
                allocate(varianceCurrentRateQuad(0:self%countVarianceCurrentRateNonCrossing))
             end if
             do iTime=1,self%countTimeRate
                ! Skip those epochs whose solutions are carried over from the previous tabulation. The test involves no
                ! thread-private state, so every thread of the team takes the same branch and the `!$omp do` below is reached by
                ! all of them or by none.
                if     (                                                                                     &
                     &        iTime >  offsetTimeRate                                                        &
                     &  .and. iTime <= offsetTimeRate+countTimeRatePrevious                                  &
                     &  .and. merge(reuseSolutions,reuseSolutionsNonCrossing,iCompute == 1)                   &
                     & ) cycle
                if (.not.allocated(firstCrossingRateQuad)) allocate(firstCrossingRateQuad(0:self%countVarianceProgenitorRate))
                ! Compute a suitable progenitor time.
                timeProgenitor=self%timeRate(iTime)*(1.0d0-self%fractionalTimeStep)
                if (iCompute == 1) then
                   varianceMaximumRateLimit=self%varianceMaximumRate
                else
                   if (self%massMinimumRateNonCrossing > 0.0d0) then
                      varianceMaximumRateLimit=cosmologicalMassVariance_%rootVariance(self%massMinimumRateNonCrossing,self%timeRate(iTime))**2
                   else
                      varianceMaximumRateLimit=self%varianceMaximumRate
                   end if
                   ! For computing non-crossing rates, the results are tabulated with respect to S_max-S so that interpolation
                   ! is more accurate when S approaches S_max.
                   do iVariance=0,countVarianceCurrentRate
                      varianceCurrentRateQuad(iVariance)=max(varianceMaximumRateLimit-self%varianceCurrentRateNonCrossing(iVariance),0.0d0)
                   end do
                end if
                ! Iterate over the variance of the current halo.
                !$omp do schedule(dynamic)
                do iVariance=0,countVarianceCurrentRate
#ifdef USEMPI
                   taskCount=taskCount+1
                   if (self%coordinatedMPI_ .and. mod(taskCount,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
#ifdef USEMPI
                   if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                      call displayCounter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityLevelWorking)
#ifdef USEMPI
                   end if
#endif
                   !$omp atomic
                   loopCount=loopCount+1_c_size_t
                   ! Construct the barrier table for progenitor halos.
                   do i=1,self%countVarianceProgenitorRate
                      massProgenitor          =+cosmologicalMassVariance_%mass        (real(sqrt(+varianceProgenitorRateQuad(i)+varianceCurrentRateQuad(iVariance)),kind=8),self%timeRate      (iTime)                        )
                      growthFactorEffective   =+cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,self%timeRate      (iTime)                        ) &
                           &                   /cosmologicalMassVariance_%rootVariance(           massProgenitor                                                           ,     timeProgenitor                               )
                      barrierRateQuad      (i)=+excursionSetBarrier_     %barrier     (real(     +varianceProgenitorRateQuad(i)+varianceCurrentRateQuad(iVariance) ,kind=8),     timeProgenitor       ,node,rateCompute=.true.) &
                           &                   *growthFactorEffective
                   end do
                   ! For zero variance, the rate is initialized to zero.
                   firstCrossingRateQuad(0)=0.0_kind_quad
                   ! Compute the step in variance across this first grid point.
                   varianceStepRate=+varianceProgenitorRateQuad(1) &
                        &           -varianceProgenitorRateQuad(0)
                   ! Compute the barrier for the descendant.
                   barrier=real(excursionSetBarrier_%barrier(real(varianceCurrentRateQuad(iVariance),kind=8),self%timeRate(iTime),node,rateCompute=.true.),kind=kind_quad)
                   ! Compute the first crossing distribution at the first grid point.
                   if (varianceProgenitorRateQuad(1)+varianceCurrentRateQuad(iVariance) >= varianceMaximumRateLimit) then
                      firstCrossingRateQuad(1)= 0.0_kind_quad
                   else
                      offsetEffective         =self%offsetEffective (self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(1),0.0_kind_quad,barrier,barrierRateQuad(1)-barrier,0.0_kind_quad,cosmologicalMassVariance_)
                      varianceResidual        =self%varianceResidual(self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(1),0.0_kind_quad                                                 ,cosmologicalMassVariance_)
                      firstCrossingRateQuad(1)=+2.0_kind_quad                                                      &
                           &                   *Error_Function_Complementary(                                      &
                           &                                                 +offsetEffective                      &
                           &                                                 /sqrt(2.0_kind_quad*varianceResidual) &
                           &                                                )                                      &
                           &                   /varianceStepRate
                   end if
                   ! Iterate over remaining progenitor variances, solving for the first crossing rate at each.
                   do i=2,self%countVarianceProgenitorRate
                      if (varianceProgenitorRateQuad(i)+varianceCurrentRateQuad(iVariance) >= varianceMaximumRateLimit) then
                         firstCrossingRateQuad(i)=0.0_kind_quad
                      else
                         barrierProgenitorEffective=+barrierRateQuad(i) &
                              &                     -barrier
                         probabilityCrossingPrior  =+0.0_kind_quad
                         ! Sum the contributions from trajectories which crossed the barrier at some smaller progenitor variance.
                         do j=1,i-1
                            varianceStepRate        =(varianceProgenitorRateQuad(j+1)-varianceProgenitorRateQuad(j-1))/2.0_kind_quad
                            offsetEffective         =self%offsetEffective (self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),varianceProgenitorRateQuad(j),barrier,barrierProgenitorEffective,barrierRateQuad(j)-barrier,cosmologicalMassVariance_)
                            varianceResidual        =self%varianceResidual(self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),varianceProgenitorRateQuad(j)                                                              ,cosmologicalMassVariance_)
                            probabilityCrossingPrior=+probabilityCrossingPrior                                           &
                                 &                   +firstCrossingRateQuad(j)                                           &
                                 &                   *varianceStepRate                                                   &
                                 &                   *Error_Function_Complementary(                                      &
                                 &                                                 +offsetEffective                      &
                                 &                                                 /sqrt(2.0_kind_quad*varianceResidual) &
                                 &                                                )
                         end do
                         varianceStepRate        =varianceProgenitorRateQuad(i)-varianceProgenitorRateQuad(i-1)
                         offsetEffective         =self%offsetEffective (self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),0.0_kind_quad,barrier,barrierProgenitorEffective,0.0_kind_quad,cosmologicalMassVariance_)
                         varianceResidual        =self%varianceResidual(self%timeRate(iTime),varianceCurrentRateQuad(iVariance),varianceProgenitorRateQuad(i),0.0_kind_quad                                                 ,cosmologicalMassVariance_)
                         firstCrossingRateQuad(i)=max(                                                                      &
                              &                       +0.0_kind_quad,                                                       &
                              &                       +(                                                                    &
                              &                         +2.0_kind_quad                                                      &
                              &                         *Error_Function_Complementary(                                      &
                              &                                                       +offsetEffective                      &
                              &                                                       /sqrt(2.0_kind_quad*varianceResidual) &
                              &                                                      )                                      &
                              &                         -2.0_kind_quad*probabilityCrossingPrior                             &
                              &                        )                                                                    &
                              &                       /varianceStepRate                                                     &
                              &                      )
                      end if
                   end do
                   if (iCompute == 1) then
                      ! Store the compute crossing rate in our table.
                      self%firstCrossingRate(:,iVariance,iTime)=real(firstCrossingRateQuad,kind=kind_dble)
                      ! Divide through by the time step to get the rate of barrier crossing.
                      self%firstCrossingRate(:,iVariance,iTime)=+self%firstCrossingRate (:,iVariance,iTime) &
                           &                                    /self%timeRate          (            iTime) &
                           &                                    /self%fractionalTimeStep
                   else
                      ! Compute the fraction of trajectories which never cross the barrier.
                      crossingFraction=0.0_kind_quad
                      do j=0,self%countVarianceProgenitorRate-1
                         if (varianceCurrentRateQuad(iVariance)+varianceProgenitorRateQuad(j+1) <= varianceMaximumRateLimit) then
                            varianceStepRate=+varianceProgenitorRateQuad(j+1) &
                                 &           -varianceProgenitorRateQuad(j)
                            crossingFraction=+crossingFraction             &
                                 &           +0.5_kind_quad                &
                                 &           *(                            &
                                 &              firstCrossingRateQuad(j  ) &
                                 &             +firstCrossingRateQuad(j+1) &
                                 &            )                            &
                                 &           *varianceStepRate
                         end if
                      end do
                      ! Compute the rate for trajectories which never cross the barrier.
                      self%nonCrossingRate(iVariance,iTime)=real(                                     &
                           &                                     +max(                                &
                           &                                          1.0_kind_quad-crossingFraction, &
                           &                                          0.0_kind_quad                   &
                           &                                         )                                &
                           &                                     /self%timeRate(iTime)                &
                           &                                     /self%fractionalTimeStep           , &
                           &                                     kind=kind_dble                       &
                           &                                    )
                   end if
                end do
                !$omp end do
             end do
          end do
          !![
          <objectDestructor name="excursionSetBarrier_"     />
          <objectDestructor name="cosmologicalMassVariance_"/>
          !!]
          ! Deallocate work arrays
          deallocate(barrierRateQuad        )
          deallocate(varianceCurrentRateQuad)
          if (allocated(firstCrossingRateQuad)) deallocate(firstCrossingRateQuad)
          !$omp end parallel
          deallocate(varianceProgenitorRateQuad)
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayCounterClear(       verbosityLevelWorking)
             call displayUnindent    ("done",verbosityLevelWorking)
#ifdef USEMPI
          end if
          if (self%coordinatedMPI_) then
             call mpiBarrier()
             if (makeTable) then
                self%firstCrossingRate=mpiSelf%sum(self%firstCrossingRate)
             end if
             self%nonCrossingRate=mpiSelf%sum(self%nonCrossingRate)
          end if
#endif
          ! Carry over the solutions already found, into the epochs skipped above. This is done after the reduction over MPI
          ! processes rather than before the solution: under coordinated MPI each process holds only the epochs it solved for
          ! and the reduction is a sum over processes, which would otherwise count a carried-over epoch once per process.
          if (reuseSolutions           ) then
             self%firstCrossingRate(:,:,offsetTimeRate+1:offsetTimeRate+countTimeRatePrevious)=firstCrossingRatePrevious
             deallocate(firstCrossingRatePrevious)
          end if
          if (reuseSolutionsNonCrossing) then
             self%nonCrossingRate  (  :,offsetTimeRate+1:offsetTimeRate+countTimeRatePrevious)=nonCrossingRatePrevious
             deallocate(nonCrossingRatePrevious  )
          end if
          ! Build the interpolators.
          if (allocated(self%interpolatorVarianceRate                  )) deallocate(self%interpolatorVarianceRate                  )
          if (allocated(self%interpolatorVarianceCurrentRate           )) deallocate(self%interpolatorVarianceCurrentRate           )
          if (allocated(self%interpolatorVarianceCurrentRateNonCrossing)) deallocate(self%interpolatorVarianceCurrentRateNonCrossing)
          if (allocated(self%interpolatorTimeRate                      )) deallocate(self%interpolatorTimeRate                      )
          allocate(self%interpolatorVarianceRate                  )
          allocate(self%interpolatorVarianceCurrentRate           )
          allocate(self%interpolatorVarianceCurrentRateNonCrossing)
          allocate(self%interpolatorTimeRate                      )
          self%interpolatorVarianceRate                  =interpolator(self%varianceProgenitorRate        ,extrapolationType=extrapolationTypeFix)
          self%interpolatorVarianceCurrentRate           =interpolator(self%varianceCurrentRate           ,extrapolationType=extrapolationTypeFix)
          self%interpolatorVarianceCurrentRateNonCrossing=interpolator(self%varianceCurrentRateNonCrossing,extrapolationType=extrapolationTypeFix)
          self%interpolatorTimeRate                      =interpolator(self%timeRate                      ,extrapolationType=extrapolationTypeFix)
          ! Set previous variance and time to unphysical values to force recompute of interpolation factors on next call.
          self%variancePreviousRate=-1.0d0
          self%timePreviousRate    =-1.0d0
          ! Record that the table is now built.
          self%tableInitializedRate     =.true.
          self%retabulateRateNonCrossing=.false.
          ! Write the table to file if possible.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             if (self%useFile) then
                call File_Lock(self%fileName,fileLock,lockIsShared=.false.)
                call self%fileWrite()
                call File_Unlock(fileLock,sync=.true.)
             end if
#ifdef USEMPI
          end if
#endif
       end if
       !$omp end critical(farahiRateTabulate)
    end if
    return
  end subroutine farahiRateTabulate

  subroutine farahiFileRead(self)
    !!{RST
    Read tabulated data on excursion set first crossing probabilities from file.
    !!}
    use :: Display           , only : displayIndent       , displayMessage     , displayUnindent, verbosityLevelWorking
    use :: File_Utilities    , only : File_Exists         , File_Name_Expand
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File            , hdf5Group
    use :: ISO_Varying_String, only : operator(//)        , var_str            , varying_string
    use :: Numerical_Ranges  , only : gridSchemePerUnit   , gridSchemePerDecade
    use :: String_Handling   , only : operator(//)
    use :: Table_Labels      , only : extrapolationTypeFix
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)                   :: self
    double precision                                 , allocatable  , dimension(:    ) :: varianceCurrentTmp           , varianceTmp    , &
         &                                                                                varianceCurrentNonCrossingTmp
    double precision                                 , allocatable  , dimension(:,:  ) :: firstCrossingProbabilityTmp  , nonCrossingRate
    double precision                                 , allocatable  , dimension(:,:,:) :: firstCrossingRateTmp
    double precision                                                                   :: massMinimumRateNonCrossing   , varianceMinimumRate
    type            (varying_string                 )                                  :: message
    character       (len=32                         )                                  :: label

    ! Initialize file name.
    call self%fileNameInitialize()
    ! Return immediately if the file does not exist.
    if (.not.File_Exists(self%fileName)) return
    ! Open the data file.
    !$ call hdf5Access%set()
    hdf5FileScope: block
      type(hdf5File ) :: dataFile
      type(hdf5Group) :: dataGroup
      dataFile=hdf5File(self%fileName,readOnly=.true.)
      ! Check if the standard table is populated.
      if (dataFile%hasGroup('probability')) then
         ! Deallocate arrays if necessary.
         if (allocated(self%variance                )) deallocate(self%variance                )
         if (allocated(self%time                    )) deallocate(self%time                    )
         if (allocated(self%firstCrossingProbability)) deallocate(self%firstCrossingProbability)
         ! Read the datasets.
         dataGroup=dataFile%openGroup("probability")
         call dataGroup%readDataset('variance'                ,varianceTmp                )
         call dataGroup%readDataset('time'                    ,self%time                  )
         call dataGroup%readDataset('firstCrossingProbability',firstCrossingProbabilityTmp)
         ! Set table sizes and limits.
         self%countVariance=size(varianceTmp)-1
         self%countTime    =size(self%time  )
         ! Transfer to tables.
         allocate(self%variance                (0:self%countVariance               ))
         allocate(self%firstCrossingProbability(0:self%countVariance,self%countTime))
         self%variance                         (0:self%countVariance  )=varianceTmp                (1:self%countVariance+1  )
         self%firstCrossingProbability         (0:self%countVariance,:)=firstCrossingProbabilityTmp(1:self%countVariance+1,:)
         deallocate(varianceTmp                )
         deallocate(firstCrossingProbabilityTmp)
         ! Place the restored tabulation on its lattices, and take the abscissae from them rather than from the file so that
         ! they are bit-identical to those of any other tabulation built on the same lattices. A file which does not record
         ! lattices, or whose lattices do not match the datasets stored alongside them or the grid densities this object would
         ! use, is not usable: everything read from it is discarded and the tabulation is rebuilt, rather than leaving a
         ! partially restored tabulation behind which could neither be extended nor trusted.
         call farahiLatticeRead(dataGroup,'variance',gridSchemePerUnit  ,self%varianceNumberPerUnitProbability,self%countVariance+1,self%latticeVariance)
         call farahiLatticeRead(dataGroup,'time'    ,gridSchemePerDecade,self%timeNumberPerDecade             ,self%countTime      ,self%latticeTime    )
         if (self%latticeVariance%isDefined() .and. self%latticeTime%isDefined()) then
            self%variance        =self%latticeVariance%values ()
            self%time            =self%latticeTime    %values ()
            self%varianceStep    =self%latticeVariance%step   ()
            self%timeMinimum     =self%latticeTime    %minimum()
            self%timeMaximum     =self%latticeTime    %maximum()
            self%varianceMaximum =self%latticeVariance%maximum()
            self%tableInitialized=.true.
         else
            deallocate(self%variance                )
            deallocate(self%time                    )
            deallocate(self%firstCrossingProbability)
            self%tableInitialized=.false.
         end if
      end if
      if (self%tableInitialized) then
         ! Build the interpolators.
         if (allocated(self%interpolatorVariance)) deallocate(self%interpolatorVariance)
         if (allocated(self%interpolatorTime    )) deallocate(self%interpolatorTime    )
         allocate(self%interpolatorVariance)
         allocate(self%interpolatorTime    )
         self%interpolatorVariance=interpolator(self%variance,extrapolationType=extrapolationTypeFix)
         self%interpolatorTime    =interpolator(self%time    ,extrapolationType=extrapolationTypeFix)
         ! Report.
         message=var_str('read excursion set first crossing probability from: ')//self%fileName
         call displayIndent  (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%timeMinimum
         message=var_str('    time minimum: ')//label//' Gyr'
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%timeMaximum
         message=var_str('    time maximum: ')//label//' Gyr'
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%varianceMaximum
         message=var_str('variance maximum: ')//label
         call displayMessage (message,verbosityLevelWorking)
         message=var_str('      table size: ')//size(self%time)//' ⨉ '//size(self%variance)
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(f7.3)') dble(sizeof(self%time)+sizeof(self%variance)+sizeof(self%firstCrossingProbability))/1024.0d0**3
         message=var_str('     memory size: ')//label//' Gib'
         call displayMessage (message,verbosityLevelWorking)
         call displayUnindent(''     ,verbosityLevelWorking)
      end if
      ! Check if the rate table is populated.
      if (dataFile%hasGroup('rate')) then
         ! Deallocate arrays if necessary.
         if (allocated(self%varianceProgenitorRate        )) deallocate(self%varianceProgenitorRate        )
         if (allocated(self%varianceCurrentRate           )) deallocate(self%varianceCurrentRate           )
         if (allocated(self%varianceCurrentRateNonCrossing)) deallocate(self%varianceCurrentRateNonCrossing)
         if (allocated(self%timeRate                      )) deallocate(self%timeRate                      )
         if (allocated(self%firstCrossingRate             )) deallocate(self%firstCrossingRate             )
         if (allocated(self%nonCrossingRate               )) deallocate(self%nonCrossingRate               )
         ! Read the datasets.
         dataGroup=dataFile%openGroup("rate")
         call dataGroup%readDataset  ('varianceProgenitor'        ,varianceTmp                  )
         call dataGroup%readDataset  ('varianceCurrent'           ,varianceCurrentTmp           )
         call dataGroup%readDataset  ('varianceCurrentNonCrossing',varianceCurrentNonCrossingTmp)
         call dataGroup%readDataset  ('time'                      ,self%timeRate                )
         call dataGroup%readDataset  ('firstCrossingRate'         ,firstCrossingRateTmp         )
         call dataGroup%readDataset  ('nonCrossingRate'           ,nonCrossingRate              )
         call dataGroup%readAttribute('massMinimumRateNonCrossing',massMinimumRateNonCrossing   )
         if (self%massMinimumRateNonCrossing == massMinimumRateNonCrossing) then
            self%retabulateRateNonCrossing=.false.
         end if
         ! Set table sizes and limits.
         self%countVarianceProgenitorRate        =size(varianceTmp                  )-1
         self%countVarianceCurrentRate           =size(varianceCurrentTmp           )-1
         self%countVarianceCurrentRateNonCrossing=size(varianceCurrentNonCrossingTmp)-1
         self%countTimeRate                      =size(self%timeRate                )
         ! Transfer to tables.
         allocate(self%varianceProgenitorRate        (0:self%countVarianceProgenitorRate                                                              ))
         allocate(self%varianceCurrentRate           (                                   0:self%countVarianceCurrentRate                              ))
         allocate(self%varianceCurrentRateNonCrossing(                                   0:self%countVarianceCurrentRateNonCrossing                   ))
         allocate(self%firstCrossingRate             (0:self%countVarianceProgenitorRate,0:self%countVarianceCurrentRate           ,self%countTimeRate))
         allocate(self%nonCrossingRate               (                                   0:self%countVarianceCurrentRateNonCrossing,self%countTimeRate))
         self%varianceProgenitorRate        (0:self%countVarianceProgenitorRate                                             )=varianceTmp                  (1:self%countVarianceProgenitorRate+1                                               )
         self%varianceCurrentRate           (                                   0:self%countVarianceCurrentRate             )=varianceCurrentTmp           (                                     1:self%countVarianceCurrentRate           +1  )
         self%varianceCurrentRateNonCrossing(                                   0:self%countVarianceCurrentRateNonCrossing  )=varianceCurrentNonCrossingTmp(                                     1:self%countVarianceCurrentRateNonCrossing+1  )
         self%firstCrossingRate             (0:self%countVarianceProgenitorRate,0:self%countVarianceCurrentRate           ,:)=firstCrossingRateTmp         (1:self%countVarianceProgenitorRate+1,1:self%countVarianceCurrentRate           +1,:)
         self%nonCrossingRate               (                                   0:self%countVarianceCurrentRateNonCrossing,:)=nonCrossingRate              (                                     1:self%countVarianceCurrentRateNonCrossing+1,:)
         deallocate(varianceTmp                  )
         deallocate(varianceCurrentTmp           )
         deallocate(varianceCurrentNonCrossingTmp)
         ! Place the restored tabulation on its lattices, taking the abscissae of the two axes which lie on lattices from those
         ! lattices rather than from the file, so that they are bit-identical to those of any other tabulation built on the same
         ! lattices. As for the probability tabulation, a file which does not record the lattices, or whose lattices do not match
         ! the datasets stored alongside them or the grid densities this object would use, is not usable and everything read from
         ! it is discarded. The two transition-spaced axes carry no lattice of their own, so they are checked instead against the
         ! lengths which the pinned bounds imply: were they to disagree, the tabulation could not be extended without moving
         ! points which the stored solutions were computed at.
         call farahiLatticeRead(dataGroup,'varianceCurrent',gridSchemePerUnit  ,self%varianceNumberPerUnit     ,self%countVarianceCurrentRate+1,self%latticeVarianceCurrentRate)
         call farahiLatticeRead(dataGroup,'time'           ,gridSchemePerDecade,self%timeNumberPerDecade       ,self%countTimeRate             ,self%latticeTimeRate           )
         call farahiLatticeRead(dataGroup,'varianceMinimum',gridSchemePerDecade,varianceMinimumAnchorsPerDecade,lattice=self%latticeVarianceMinimumRate                        )
         if     (                                                                                                                                                        &
              &        self%latticeVarianceCurrentRate%isDefined()                                                                                                       &
              &  .and. self%latticeTimeRate           %isDefined()                                                                                                       &
              &  .and. self%latticeVarianceMinimumRate%isDefined()                                                                                                       &
              & ) then
            self%varianceMaximumRate=self%latticeVarianceCurrentRate%maximum()
            varianceMinimumRate     =self%latticeVarianceMinimumRate%minimum()
            self%tableInitializedRate=       self%countVarianceProgenitorRate         == int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecade           ))+1 &
                 &                    .and.  self%countVarianceCurrentRateNonCrossing == int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecadeNonCrossing))+1
         else
            self%tableInitializedRate=.false.
         end if
         if (self%tableInitializedRate) then
            self%varianceCurrentRate=self%latticeVarianceCurrentRate%values ()
            self%timeRate           =self%latticeTimeRate           %values ()
            self%timeMinimumRate    =self%latticeTimeRate           %minimum()
            self%timeMaximumRate    =self%latticeTimeRate           %maximum()
         else
            deallocate(self%varianceProgenitorRate        )
            deallocate(self%varianceCurrentRate           )
            deallocate(self%varianceCurrentRateNonCrossing)
            deallocate(self%timeRate                      )
            deallocate(self%firstCrossingRate             )
            deallocate(self%nonCrossingRate               )
         end if
      end if
      if (self%tableInitializedRate) then
         ! Build the interpolators.
         if (allocated(self%interpolatorVarianceRate                  )) deallocate(self%interpolatorVarianceRate                  )
         if (allocated(self%interpolatorVarianceCurrentRate           )) deallocate(self%interpolatorVarianceCurrentRate           )
         if (allocated(self%interpolatorVarianceCurrentRateNonCrossing)) deallocate(self%interpolatorVarianceCurrentRateNonCrossing)
         if (allocated(self%interpolatorTimeRate                      )) deallocate(self%interpolatorTimeRate                      )
         allocate(self%interpolatorVarianceRate                  )
         allocate(self%interpolatorVarianceCurrentRate           )
         allocate(self%interpolatorVarianceCurrentRateNonCrossing)
         allocate(self%interpolatorTimeRate                      )
         self%interpolatorVarianceRate                  =interpolator(self%varianceProgenitorRate        ,extrapolationType=extrapolationTypeFix)
         self%interpolatorVarianceCurrentRate           =interpolator(self%varianceCurrentRate           ,extrapolationType=extrapolationTypeFix)
         self%interpolatorVarianceCurrentRateNonCrossing=interpolator(self%varianceCurrentRateNonCrossing,extrapolationType=extrapolationTypeFix)
         self%interpolatorTimeRate                      =interpolator(self%timeRate                      ,extrapolationType=extrapolationTypeFix)
         ! Report.
         message=var_str('read excursion set first crossing rates from: ')//self%fileName
         call displayIndent  (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%timeMinimumRate
         message=var_str('    time minimum: ')//label//' Gyr'
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%timeMaximumRate
         message=var_str('    time maximum: ')//label//' Gyr'
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%varianceMaximumRate
         message=var_str('variance maximum: ')//label
         call displayMessage (message,verbosityLevelWorking)
         message=var_str('      table size: ')//size(self%timeRate)//' ⨉ '//size(self%varianceProgenitorRate)//' ⨉ '//size(self%varianceCurrentRate)
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(f7.3)') dble(sizeof(self%timeRate)+sizeof(self%varianceProgenitorRate)+sizeof(self%varianceCurrentRate)+sizeof(self%firstCrossingRate)+sizeof(self%nonCrossingRate))/1024.0d0**3
         message=var_str('     memory size: ')//label//' Gib'
         call displayMessage (message,verbosityLevelWorking)
         call displayUnindent(''     ,verbosityLevelWorking)
      end if
    end block hdf5FileScope
    !$ call hdf5Access%unset()
    return
  end subroutine farahiFileRead

  subroutine farahiFileWrite(self)
    !!{RST
    Write tabulated data on excursion set first crossing probabilities to file.
    !!}
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File     , hdf5Group
    use :: Display           , only : displayIndent, displayMessage, displayUnindent, verbosityLevelWorking
    use :: HDF5              , only : hsize_t
    use :: ISO_Varying_String, only : operator(//) , var_str       , varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    class    (excursionSetFirstCrossingFarahi), intent(inout) :: self
    type     (varying_string                 )                :: message
    character(len=32                         )                :: label

    ! Don't write anything if neither table is initialized.
    if (.not.(self%tableInitialized.or.self%tableInitializedRate)) return
    ! Initialize file name.
    call self%fileNameInitialize()
    ! Open the data file.
    !$ call hdf5Access%set()
    hdf5FileScope: block
      type(hdf5File ) :: dataFile
      type(hdf5Group) :: dataGroup
      dataFile=hdf5File(self%fileName,overWrite=.true.,chunkSize=100_hsize_t,compressionLevel=9)
      ! Check if the standard table is populated.
      if (self%tableInitialized) then
         dataGroup=dataFile%openGroup("probability")
         call dataGroup%writeDataset(self%variance                ,'variance'                ,'The variance at which results are tabulated.'                         )
         call dataGroup%writeDataset(self%time                    ,'time'                    ,'The cosmic times at which results are tabulated.'                     )
         call dataGroup%writeDataset(self%firstCrossingProbability,'firstCrossingProbability','The probability of first crossing as a function of variance and time.')
         ! Record the lattices alongside the tabulation. The bounds follow from the lattice, but the converse does not, and it
         ! is the lattice which determines whether what is read back can be extended to cover a new range without recomputation.
         call farahiLatticeWrite(dataGroup,'variance',self%latticeVariance)
         call farahiLatticeWrite(dataGroup,'time'    ,self%latticeTime    )
         ! Report.
         message=var_str('write excursion set first crossing probability to: ')//self%fileName
         call displayIndent  (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%timeMinimum
         message=var_str('    time minimum: ')//label//' Gyr'
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%timeMaximum
         message=var_str('    time maximum: ')//label//' Gyr'
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%varianceMaximum
         message=var_str('variance maximum: ')//label
         call displayMessage (message,verbosityLevelWorking)
         message=var_str('      table size: ')//size(self%time)//' ⨉ '//size(self%variance)
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(f7.3)') dble(sizeof(self%time)+sizeof(self%variance)+sizeof(self%firstCrossingProbability))/1024.0d0**3
         message=var_str('     memory size: ')//label//' Gib'
         call displayUnindent(''     ,verbosityLevelWorking)
      end if
      ! Check if the rate table is populated.
      if (self%tableInitializedRate) then
         dataGroup=dataFile%openGroup("rate")
         call dataGroup%writeDataset  (self%varianceProgenitorRate        ,'varianceProgenitor'        ,'The variance at which results are tabulated.'                               )
         call dataGroup%writeDataset  (self%varianceCurrentRate           ,'varianceCurrent'           ,'The variance of the base halo at which first crossing rates are tabulated.' )
         call dataGroup%writeDataset  (self%varianceCurrentRateNonCrossing,'varianceCurrentNonCrossing','The variance of the base halo at which non-crossing rates are tabulated.'   )
         call dataGroup%writeDataset  (self%timeRate                      ,'time'                      ,'The cosmic times at which results are tabulated.'                           )
         call dataGroup%writeDataset  (self%firstCrossingRate             ,'firstCrossingRate'         ,'The probability rate of first crossing as a function of variances and time.')
         call dataGroup%writeDataset  (self%nonCrossingRate               ,'nonCrossingRate'           ,'The probability rate of non crossing as a function of variance and time.'   )
         call dataGroup%writeAttribute(self%massMinimumRateNonCrossing    ,'massMinimumRateNonCrossing'                                                                              )
         ! Record the lattices alongside the tabulation. Only two of the four axes lie on a lattice; the third lattice carries no
         ! tabulated points at all, and records the pinned lower bound from which the remaining two axes are generated. All three
         ! are needed to establish whether what is read back can be extended without recomputation.
         call farahiLatticeWrite(dataGroup,'varianceCurrent',self%latticeVarianceCurrentRate)
         call farahiLatticeWrite(dataGroup,'time'           ,self%latticeTimeRate           )
         call farahiLatticeWrite(dataGroup,'varianceMinimum',self%latticeVarianceMinimumRate)
         ! Report.
         message=var_str('wrote excursion set first crossing rates to: ')//self%fileName
         call displayIndent  (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%timeMinimumRate
         message=var_str('    time minimum: ')//label//' Gyr'
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%timeMaximumRate
         message=var_str('    time maximum: ')//label//' Gyr'
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(e22.16)') self%varianceMaximumRate
         message=var_str('variance maximum: ')//label
         call displayMessage (message,verbosityLevelWorking)
         message=var_str('      table size: ')//size(self%timeRate)//' ⨉ '//size(self%varianceProgenitorRate)//' ⨉ '//size(self%varianceCurrentRate)
         call displayMessage (message,verbosityLevelWorking)
         write (label,'(f7.3)') dble(sizeof(self%timeRate)+sizeof(self%varianceProgenitorRate)+sizeof(self%varianceCurrentRate)+sizeof(self%firstCrossingRate)+sizeof(self%nonCrossingRate))/1024.0d0**3
         message=var_str('     memory size: ')//label//' Gib'
         call displayMessage (message,verbosityLevelWorking)
         call displayUnindent(''     ,verbosityLevelWorking)
      end if
    end block hdf5FileScope
    !$ call hdf5Access%unset()
    return
  end subroutine farahiFileWrite

  subroutine farahiLatticeWrite(dataGroup,axisName,lattice)
    !!{RST
    Record the ``rangeLattice`` on which an axis of a stored tabulation is built, as attributes named for that axis.
    !!}
    use :: IO_HDF5, only : hdf5Group
    implicit none
    type     (hdf5Group   ), intent(inout) :: dataGroup
    character(len=*       ), intent(in   ) :: axisName
    type     (rangeLattice), intent(in   ) :: lattice

    call dataGroup%writeAttribute(lattice%scheme%ID   ,axisName//'GridScheme'  )
    call dataGroup%writeAttribute(lattice%pointsPer   ,axisName//'PointsPer'   )
    call dataGroup%writeAttribute(lattice%indexMinimum,axisName//'IndexMinimum')
    call dataGroup%writeAttribute(lattice%count       ,axisName//'Count'       )
    return
  end subroutine farahiLatticeWrite

  subroutine farahiLatticeRead(dataGroup,axisName,scheme,pointsPer,countExpected,lattice)
    !!{RST
    Restore the ``rangeLattice`` on which an axis of a stored tabulation was built. The lattice is returned undefined---which
    the caller must treat as the tabulation being unusable---unless the file records one, and that lattice is self-consistent,
    uses the gridding scheme and density of points which this object would use, and has the same number of points as the
    datasets stored alongside it. Older files, written before the lattices were recorded, therefore report an undefined lattice
    rather than being misread.

    ``countExpected`` is omitted for a lattice which carries no tabulated points---one recorded solely to pin a bound---since
    there is then no dataset whose length it should match.
    !!}
    use :: IO_HDF5         , only : hdf5Group
    use :: Numerical_Ranges, only : enumerationGridSchemeType
    implicit none
    type     (hdf5Group                ), intent(inout)           :: dataGroup
    character(len=*                    ), intent(in   )           :: axisName
    type     (enumerationGridSchemeType), intent(in   )           :: scheme
    integer                             , intent(in   )           :: pointsPer
    integer                             , intent(in   ), optional :: countExpected
    type     (rangeLattice             ), intent(  out)           :: lattice
    integer                                                       :: schemeStored, pointsPerStored, &
         &                                                           indexMinimum, count_

    lattice=rangeLattice()
    if     (                                                       &
         &   .not.dataGroup%hasAttribute(axisName//'GridScheme'  ) &
         &  .or.                                                   &
         &   .not.dataGroup%hasAttribute(axisName//'PointsPer'   ) &
         &  .or.                                                   &
         &   .not.dataGroup%hasAttribute(axisName//'IndexMinimum') &
         &  .or.                                                   &
         &   .not.dataGroup%hasAttribute(axisName//'Count'       ) &
         & ) return
    call dataGroup%readAttribute(axisName//'GridScheme'  ,schemeStored   )
    call dataGroup%readAttribute(axisName//'PointsPer'   ,pointsPerStored)
    call dataGroup%readAttribute(axisName//'IndexMinimum',indexMinimum   )
    call dataGroup%readAttribute(axisName//'Count'       ,count_         )
    ! Comparing the stored scheme against the one expected is stronger than merely checking that it is a valid member of the
    ! enumeration, so no separate validity test is needed.
    if (enumerationGridSchemeType(schemeStored) /= scheme       ) return
    if (pointsPerStored                         /= pointsPer    ) return
    if (present(countExpected)) then
       if (count_                               /= countExpected) return
    end if
    lattice=rangeLattice(enumerationGridSchemeType(schemeStored),pointsPerStored,indexMinimum,count_)
    if (.not.lattice%isDefined()) lattice=rangeLattice()
    return
  end subroutine farahiLatticeRead

  function farahiVarianceRange(self,rangeMinimum,rangeMaximum,rangeNumber,exponent) result (rangeValues)
    !!{RST
    Builds a numerical range between ``rangeMinimum`` and ``rangeMaximum`` using ``rangeNumber`` points with spacing that varies from logarithmic to linear spacing with the transition point controlled by ``exponent``. Specifically, suppose we have :math:`N=`\ ``rangeNumber`` points in the range, from :math:`S_\mathrm{min}=`\ ``rangeMinimum`` to :math:`S_\mathrm{max}=`\ ``rangeMaximum``. We define :math:`f_i=(i-1)/(N-1)` where :math:`i` runs from :math:`1` to :math:`N`. We then define:

    .. math::

       f_i = { \int_{S_\mathrm{min}}^{S_i} x^{n_i} \mathrm{d} x \over \int_{S_\mathrm{min}}^{S_\mathrm{max}} x^{n_i} \mathrm{d} x},

    and solve for :math:`S_i` to find the :math:`i^\mathrm{th}` range value. If :math:`n_i=0` then this will give :math:`S_i` linearly spaced between :math:`S_\mathrm{min}` and :math:`S_\mathrm{max}`, while if :math:`n_i=-1` this will give :math:`S_i` logarithmically spaced between :math:`S_\mathrm{min}` and :math:`S_\mathrm{max}`. Therefore, if we make :math:`n_i` vary from :math:`-1` to :math:`0` at :math:`i` ranges from :math:`1` to :math:`N` we will get a smooth transition from logarithmic to linear spacing. We choose to use :math:`n_i=-1+f_i^\alpha` where :math:`\alpha=`\ ``exponent`` is a supplied argument.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)          :: self
    double precision                                 , intent(in   )          :: rangeMaximum , rangeMinimum    , &
         &                                                                       exponent
    integer                                          , intent(in   )          :: rangeNumber
    double precision                                 , dimension(rangeNumber) :: rangeValues
    integer                                                                   :: iRange
    double precision                                                          :: fractionRange, integrandExponent
    !$GLC attributes unused :: self

    do iRange=1,rangeNumber
       fractionRange    =dble(iRange-1)/dble(rangeNumber-1)
       integrandExponent=-1.0d0+fractionRange**exponent
       if (integrandExponent == -1.0d0) then
          rangeValues(iRange)=exp(log(rangeMinimum)                          +log(rangeMaximum                           /rangeMinimum                           )*fractionRange)
       else
          rangeValues(iRange)=   (    rangeMinimum**(1.0d0+integrandExponent)+   (rangeMaximum**(1.0d0+integrandExponent)-rangeMinimum**(1.0d0+integrandExponent))*fractionRange)**(1.00/(1.0d0+integrandExponent))
       end if
    end do
    return
  end function farahiVarianceRange

  double precision function farahiVarianceLimit(self,varianceProgenitor)
    !!{RST
    Return the maximum variance which the rate tabulation is required to reach for the given progenitor variance.

    Note that this is the variance *requested*, and takes no account of the variance already tabulated: the union with the
    range in use is taken on the lattice itself, in exact integer arithmetic. Folding it in here instead would apply the
    safety margin to a bound which had already been given one, so that the tabulated range crept upwards on every retabulation
    and no solution already found could ever be reused.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout) :: self
    double precision                                 , intent(in   ) :: varianceProgenitor

    if (self%varianceIsUnlimited) then
       farahiVarianceLimit=2.0d0*varianceProgenitor
    else
       farahiVarianceLimit=      varianceProgenitor
    end if
    return
  end function farahiVarianceLimit

  double precision function farahiVarianceLimitHard(self)
    !!{RST
    Return a hard upper limit on the variance to which the rate tabulation may extend, beyond which the solution is not merely
    unnecessary but undefined. There is no such limit for this class, so ``huge(0)`` is returned to indicate its absence---the
    tabulated range is then bounded only by what is requested, plus the safety margin.
    !!}
    implicit none
    class(excursionSetFirstCrossingFarahi), intent(inout) :: self
    !$GLC attributes unused :: self

    farahiVarianceLimitHard=huge(0.0d0)
    return
  end function farahiVarianceLimitHard

  function farahiRateVarianceLattice(self,varianceProgenitor) result(lattice)
    !!{RST
    Return the lattice on which the variance of the current halo is tabulated when solving for barrier crossing rates. This
    axis is linear in variance and begins at zero, so it lies on a ``perUnit`` lattice; zero is included among the target
    values, and imposed as a hard lower limit as well so that the safety margin cannot carry the range below it. No seed is
    needed, since every range begins at lattice index zero and so any two of them overlap. The axis is anchored to tenths of a
    unit for the same reason as the probability tabulation: the progenitor solution at each of its points is a recursion, so
    the cost is linear in the length of this axis and quadratic in that of the progenitor axis, and whole-unit anchoring on a
    variance of order unity would be extravagant.

    The maximum variance of this axis sets all three of the variance axes of the rate tabulation, the other two being
    generated between it and the pinned minimum variance. It is therefore built here, once, rather than separately by each
    class which solves for rates.
    !!}
    use :: Numerical_Ranges, only : Range_Pinned, gridSchemePerUnit
    implicit none
    type            (rangeLattice                   )                :: lattice
    class           (excursionSetFirstCrossingFarahi), intent(inout) :: self
    double precision                                 , intent(in   ) :: varianceProgenitor
    double precision                                                 :: varianceMaximumHard

    varianceMaximumHard=self%varianceLimitHard()
    if (varianceMaximumHard < huge(0.0d0)) then
       lattice=Range_Pinned(                                                                         &
            &                              [0.0d0,self%varianceLimit(varianceProgenitor)]          , &
            &                              self%varianceNumberPerUnit                              , &
            &                              gridSchemePerUnit                                       , &
            &               marginOffset  =1.0d0/dble(varianceAnchorsPerUnit)                      , &
            &               limitMinimum  =0.0d0                                                   , &
            &               limitMaximum  =varianceMaximumHard                                     , &
            &               anchorEvery   =max(1,self%varianceNumberPerUnit/varianceAnchorsPerUnit), &
            &               latticeCurrent=self%latticeVarianceCurrentRate                           &
            &              )
    else
       lattice=Range_Pinned(                                                                         &
            &                              [0.0d0,self%varianceLimit(varianceProgenitor)]          , &
            &                              self%varianceNumberPerUnit                              , &
            &                              gridSchemePerUnit                                       , &
            &               marginOffset  =1.0d0/dble(varianceAnchorsPerUnit)                      , &
            &               limitMinimum  =0.0d0                                                   , &
            &               anchorEvery   =max(1,self%varianceNumberPerUnit/varianceAnchorsPerUnit), &
            &               latticeCurrent=self%latticeVarianceCurrentRate                           &
            &              )
    end if
    return
  end function farahiRateVarianceLattice

  function farahiVarianceResidual(self,time,varianceCurrent,varianceProgenitor,varianceIntermediate,cosmologicalMassVariance_) result(varianceResidual)
    !!{RST
    Return the residual variance between two points for a standard Weiner process.
    !!}
    use :: Kind_Numbers, only : kind_quad
    implicit none
    real            (kind_quad                      )                :: varianceResidual
    class           (excursionSetFirstCrossingFarahi), intent(inout) :: self
    real            (kind_quad                      ), intent(in   ) :: varianceCurrent          , varianceIntermediate, &
         &                                                              varianceProgenitor
    double precision                                 , intent(in   ) :: time
    class           (cosmologicalMassVarianceClass  ), intent(inout) :: cosmologicalMassVariance_
    !$GLC attributes unused :: self, varianceCurrent, time, cosmologicalMassVariance_
    
    ! Note that this solver follows the convention used through Galacticus that σ(M) grows following linear theory. That is:
    !
    !  • the root-variance of the density field smoothed on a mass scale M is a function of time, σ(M,t) = σ(M,t₀) D(t)/D(t₀),
    !    where D(t) is the linear growth factor (which may also be scale-dependent);
    !  • the critical overdensity for collapse does not include a factor of the linear growth factor, i.e. δ_c ≅ 1.686 at all
    !    epochs (varying only due to the weak dependence on the epoch-dependent cosmological parameters).
    !
    ! This differs from standard treatments of the excursion set problem in which typically the root-variance, σ(M), is evaluated
    ! at z=0, and the critical overdensity for collapse is replaced with δ_c(t)/D(t). Mathematically these two approaches are
    ! equivalent, but it can be important to keep these distinctions in mind.

    ! In this function the following translations between internal variable names and math symbols are used:
    !
    !   S₁ = varianceCurrent
    !   S̃ = varianceProgenitor  +varianceCurrent
    !   S  = varianceIntermediate+varianceCurrent
    !   δ₁ = deltaCurrent
    !   δ̃  = deltaProgenitor     +deltaCurrent
    !   δ  = deltaIntermediate   +deltaCurrent
    !    
    ! Note that the variables "varianceIntermediate" and "varianceProgenitor" are defined to be the variances in excess of S₁ -
    ! which is why they appear with "varianceCurrent" added to them in the above.
    !
    ! This function is used in the calculation of the distribution of δ at some S for trajectories originating from (S₁,δ₁) and
    ! which did not cross the barrier at any intermediate variance. As such suffixes in variable names have the following
    ! meanings:
    !
    !   "Current"      - refers to the current halo being considered for branching, i.e. the halo existing at point (S₁,δ₁);
    !   "Progenitor"   - refers to the potential progenitor halo being considered, i.e. the halo corresponding to some variance S > S₁;
    !   "Intermediate" - refers to the intermediate variance, S̃ (with S₁ < S̃ < S).
    varianceResidual=+varianceProgenitor   &
         &           -varianceIntermediate
    return
  end function farahiVarianceResidual

  function farahiOffsetEffective(self,time,varianceCurrent,varianceProgenitor,varianceIntermediate,deltaCurrent,deltaProgenitor,deltaIntermediate,cosmologicalMassVariance_) result(offsetEffective)
    !!{RST
    Return the residual variance between two points for a standard Weiner process.
    !!}
    use :: Kind_Numbers, only : kind_quad
    implicit none
    real            (kind_quad                      )                :: offsetEffective
    class           (excursionSetFirstCrossingFarahi), intent(inout) :: self
    real            (kind_quad                      ), intent(in   ) :: deltaCurrent             , deltaIntermediate , &
         &                                                              deltaProgenitor          , varianceCurrent   , &
         &                                                              varianceIntermediate     , varianceProgenitor
    double precision                                 , intent(in   ) :: time
    class           (cosmologicalMassVarianceClass  ), intent(inout) :: cosmologicalMassVariance_
    !$GLC attributes unused :: self, deltaCurrent, varianceCurrent, varianceIntermediate, varianceProgenitor, time, cosmologicalMassVariance_
    
    ! Note that this solver follows the convention used through Galacticus that σ(M) grows following linear theory. That is:
    !
    !  • the root-variance of the density field smoothed on a mass scale M is a function of time, σ(M,t) = σ(M,t₀) D(t)/D(t₀),
    !    where D(t) is the linear growth factor (which may also be scale-dependent);
    !  • the critical overdensity for collapse does not include a factor of the linear growth factor, i.e. δ_c ≅ 1.686 at all
    !    epochs (varying only due to the weak dependence on the epoch-dependent cosmological parameters).
    !
    ! This differs from standard treatments of the excursion set problem in which typically the root-variance, σ(M), is evaluated
    ! at z=0, and the critical overdensity for collapse is replaced with δ_c(t)/D(t). Mathematically these two approaches are
    ! equivalent, but it can be important to keep these distinctions in mind.

    ! In this function the following translations between internal variable names and math symbols are used:
    !
    !   S₁ = varianceCurrent
    !   S̃ = varianceProgenitor  +varianceCurrent
    !   S  = varianceIntermediate+varianceCurrent
    !    
    ! Note that the variables "varianceIntermediate" and "varianceProgenitor" are defined to be the variances in excess of S₁ -
    ! which is why they appear with "varianceCurrent" added to them in the above.
    !
    ! This function is used in the calculation of the distribution of δ at some S for trajectories originating from (S₁,δ₁) and
    ! which did not cross the barrier at any intermediate variance. As such suffixes in variable names have the following
    ! meanings:
    !
    !   "Current"      - refers to the current halo being considered for branching, i.e. the halo existing at point (S₁,δ₁);
    !   "Progenitor"   - refers to the potential progenitor halo being considered, i.e. the halo corresponding to some variance S > S₁;
    !   "Intermediate" - refers to the intermediate variance, S̃ (with S₁ < S̃ < S).
    offsetEffective=+deltaProgenitor   &
         &          -deltaIntermediate
    return
  end function farahiOffsetEffective
