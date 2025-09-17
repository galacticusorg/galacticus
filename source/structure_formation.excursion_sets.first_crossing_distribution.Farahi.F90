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

!+    Contributions to this file made by: Arya Farahi, Andrew Benson, Christoph Behrens, Xiaolong Du.

!!{
Implements a excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Excursion_Sets_Barriers   , only : excursionSetBarrierClass
  use :: Numerical_Interpolation   , only : interpolator

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahi">
   <description>
    An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}. For trajectories originating
    from a point $(S_1,\delta_1)$, the distribution of first crossings of a barrier $B(S)$, $f(S)$, is obtained by finding the
    solution to the integral equation:    
    \begin{equation}
      1 =  \int_0^S f(\tilde{S})\mathrm{d}\tilde{S} + \int_{-\infty}^{B(S)} P(\delta,S) \mathrm{d} \delta,
     \label{eq:OldExcursionMethod}
    \end{equation}
    where $P(\delta,S) \mathrm{d} \delta$ is the probability for a trajectory to lie between $\delta$ and $\delta + \mathrm{d}
    \delta$ at variance $S$, having originated from the point $(S_1,\delta_1)$ having not crossed the barrier at any smaller
    $\tilde{S} &lt; S$. In the absence of a barrier, $P(\delta,S)$ would be equal to $P_0(\delta,S)$. However, since some
    trajectories will have crossed the barrier at $\tilde{S} &lt; S$ we must subtract off their contribution to
    $P_0(\delta,S)$. Writing the distribution of $\delta$ at $S$ for trajectories originating at some $(\tilde{\delta},\tilde{S})$
    as $P^\prime(\delta|S,\tilde{\delta},\tilde{S})$, we can therefore write
    \begin{equation}
      P(\delta,S) = P_0(\delta,S) - \int_{0}^{S} f(\tilde{S}) P^\prime(\delta|S,\tilde{\delta},\tilde{S}) \mathrm{d}\tilde{S},
    \end{equation}

    For an interesting class of cases, both $P_0(\delta,S)$ and $P^\prime(\delta|S,\tilde{\delta},\tilde{S})$ are normal
    distributions, and we can write
    \begin{equation}
     P_0(\delta,S) = \frac{1}{\sqrt{2 \pi S}} \exp\left\{-{\Delta \delta^2 [\delta,\delta_1,S,S_1] \over 2 \Delta S [S,S_1]}\right\},
    \end{equation}
    and
    \begin{equation}
     P^\prime(\delta|S,\tilde{\delta},\tilde{S}) = \frac{1}{\sqrt{2 \pi \Delta S[S,\tilde{S}]}} \exp\left\{-{\Delta \delta^2 [\delta,\tilde{\delta},S,\tilde{S}] \over 2 \Delta S [S,\tilde{S}]}\right\},
    \end{equation}
    where we refer to $\Delta \delta[\delta,\tilde{\delta},S,\tilde{S}]$ as the ``effective offset'', and to $\Delta
    S[S,\tilde{S}] = \mathrm{Var}(S) - \mathrm{Cov}(S,\tilde{S})$ as the ``residual variance''. Note that $\mathrm{Cov}(S,S_1) =
    0$ (since all trajectories pass through $\delta_1$ at $S_1$), and so $\Delta S[S,S_1] = \mathrm{Var}(S)$.

    For a standard Weiner process (such as applies to the standard case considered in excursion set theory, namely uncorrelated
    and unconstrained steps), we have trivially that
    \begin{equation}
     \Delta \delta[\delta,\tilde{\delta},S,\tilde{S}] = \delta - \tilde{\delta},
    \end{equation}
    and
    \begin{equation}
     \Delta S[S,\tilde{S}] = S - \tilde{S},
    \end{equation}
    since the Weiner process is invariant under translations of the starting point.

    Using the above results, we can rewrite eqn.~(\ref{eq:OldExcursionMethod}):
    \begin{equation}
      1 = \int_0^S f(\tilde{S})\mathrm{d}\tilde{S} + \int_{-\infty}^{B(S)} \left[ P_0(\delta,S) - \int_{0}^{S} f(\tilde{S})
      P^\prime(\delta|S,\tilde{\delta},\tilde{S}) \mathrm{d} \delta \right] ,
    \end{equation}
    in general and, for the case of a Gaussian distribution:
    \begin{equation}
      1 = \int_0^S f(\tilde{S})\mathrm{d}\tilde{S} + \int_{-\infty}^{B(S)} \left[ \frac{1}{\sqrt{2 \pi \Delta S[S,S_1]}}
      \exp\left(-\frac{\Delta \delta^2[\delta,\delta_1,S,S_1]}{2 \Delta S[S,S_1]}\right) - \int_{0}^{S} f(\tilde{S}) \frac{1}{\sqrt{2 \pi \Delta S [S,\tilde{S}]}}
      \exp\left(-\frac{\Delta \delta^2[\delta,B(\tilde{S}),S,\tilde{S}]}{2 \Delta S [S,\tilde{S}]}\right)\mathrm{d}\tilde{S} \right] \mathrm{d} \delta .
    \end{equation}
    The integral over $\mathrm{d}\delta$ can be carried out analytically to give:
    \begin{equation}
     1 = \int_0^S f(\tilde{S})\mathrm{d}\tilde{S}+ \hbox{erf}\left\{\frac{\Delta \delta [B(S),\delta_1,S,S_1]}{\sqrt{2\Delta S[S,S_1]}}\right\} - \int_{0}^{S} f(\tilde{S})
     \hbox{erf}\left\{\frac{\Delta \delta [B(S),B(\tilde{S}),S,\tilde{S}]}{\sqrt{2 \Delta S [S,\tilde{S}]}}\right\} \mathrm{d}S^{\prime\prime}.
    \label{eq:NewExcursionMethod}
    \end{equation}
    We now discretize eqn.~(\ref{eq:NewExcursionMethod}). Specifically, we divide the $S$ space into $N$ intervals defined by
    the points:
    \begin{equation}
     S_i = \left\{ \begin{array}{ll}
                    0 &amp; \hbox{if } i=0 \\
                    \sum_{j=0}^{i-1} \Delta S_j &amp; \hbox{if } i &gt; 1.
                   \end{array}
           \right.
    \end{equation}
    Note that $f(0)=0$ by definition, so $f(S_0)=0$ always. We choose $\Delta S_i = S_\mathrm{max}/N$ (i.e. uniform spacing in
    $S$) when computing first crossing distributions, and $\Delta S_i \propto S_i$ (i.e. uniform spacing in $\log(S)$) when
    computing first crossing rates.
    
    Discretizing the integrals in eqn.~(\ref{eq:NewExcursionMethod}) gives:
    \begin{equation} \label{eq:Des1}
     \int_0^{S_i} f(\tilde{S})\d \tilde{S} = \sum_{j=0}^{i-1} \frac{f(S_j) + f(S_{j+1})}{2} \Delta S_j
    \end{equation}
    and:
    \begin{equation} \label{eq:Des2}
     \int_{0}^{S_i} f(\tilde{S}) \hbox{erf}\left\{\frac{\Delta \delta [B(S),B(\tilde{S}),S,\tilde{S}]}{\sqrt{2 \Delta S[S,\tilde{S}]}}\right\} \d \tilde{S} =
     \sum_{j=0}^{i-1} \frac{1}{2} \left(f(S_j) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_j), S_i, S_j]}{\sqrt{2 \Delta S[S_i,S_j]}}\right\} + f(S_{j+1})
     \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_{j+1}),S_i,S_{j+1}]}{\sqrt{2 \Delta S[S_i,S_{j+1}]}}\right\} \right) \Delta S_j.
    \end{equation}
    We can now rewrite eqn.~(\ref{eq:NewExcursionMethod}) in discretized form:
    \begin{equation} \label{eq:DesFinal1}
     1 = \sum_{j=0}^{i-1} \frac{f(S_j) + f(S_{j+1})}{2} \Delta S_j + \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2 \Delta S[S_i,S_1]}}\right\} -
     \frac{1}{2} \sum_{j=0}^{i-1} \left( f(S_j) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_j),S_i,S_j]}{\sqrt{2 \Delta S[S_i,S_j]}}\right\} + f(S_{j+1})
     \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_{j+1}),S_i,S_{j+1}]}{\sqrt{2 \Delta S[S_i,S_{j+1}]}}\right\} \right) \Delta S_j.
    \end{equation}
    Solving eqn.~(\ref{eq:DesFinal1}) for $f(S_i)$:
    \begin{eqnarray} \label{eq:DesFinal11}
     \left( \frac{1}{2} - \frac{1}{2} \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_i), S_i, S_i]}{\sqrt{2 \Delta S[S_i,S_i]}}\right\} \right) \Delta S_{i-1}
     f(S_i) &amp;=&amp; 1 - \sum_{j=0}^{i-2} \frac{f(S_j) + f(S_{j+1})}{2} \Delta S_j - \frac{f(S_{i-1})}{2} \Delta S_{i-1} -
     \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2 \Delta S[S_i,S_1]}}\right\} \nonumber\\
    &amp; &amp; + \frac{1}{2} \sum_{j=0}^{i-2} \left( f(S_j) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_j),S_i,S_j]}{\sqrt{2 \Delta S [S_i,S_j]}}\right\} +
    f(S_{j+1}) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_{j+1}),S_i,S_{j+1}]}{\sqrt{2 \Delta S[S_i,S_{j+1}]}}\right\} \right)\Delta S_j \nonumber \\
     &amp; &amp; + \frac{1}{2} f(S_{i-1}) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_{i-1}),S_i,S_{i-1}]}{\sqrt{2 \Delta S [S_i,S_{i-1}]}}\right\} \Delta S_{i-1}.
    \end{eqnarray}
    For all barriers that we consider:
    \begin{equation} 
    \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_i),S_i,S_i]}{\sqrt{2 \Delta S[S_i,S_i]}}\right\} = 0.
    \end{equation}
    We can then simplify eqn.~(\ref{eq:DesFinal11}):
    \begin{eqnarray} \label{eq:DesFinal2}
       f(S_i) &amp;=&amp; {2 \over \Delta S_{i-1}}\left[1 - \sum_{j=0}^{i-2} \frac{f(S_j) + f(S_{j+1})}{2} \Delta S_j -
       \frac{f(S_{i-1})}{2} \Delta S_{i-1} - \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2 \Delta S [S_i,S_1] }}\right\} \right.  \nonumber\\
    &amp; &amp; + \frac{1}{2} \sum_{j=0}^{i-2} \left( f(S_j) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_j),S_i,S_j]}{\sqrt{2 \Delta S [S_i,S_j]}}\right\} +
    f(S_{j+1}) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_{j+1}),S_i,S_{j+1}]}{\sqrt{2 \Delta S [S_i,S_{j+1}]}}\right\} \right)\Delta S_j \nonumber \\
     &amp; &amp; \left. + \frac{1}{2} f(S_{i-1}) \hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_{i-1}),S_i,S_{i-1}]}{\sqrt{2 \Delta S [S_i,S_{i-1}]}}\right\} \Delta
     S_{i-1}\right].
    \end{eqnarray}
    Consolidating terms in the summations:
    \begin{equation} \label{eq:DesFinal2a}
       f(S_i) = {2 \over \Delta S_{i-1}}\left[1 - \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2\Delta S[S_i,S_1]}}\right\} - \sum_{j=0}^{i-1}
       \left( 1-\hbox{erf}\left\{\frac{\Delta \delta [B(S_i) , B(S_j),S_i,S_j]}{\sqrt{2 \Delta S [S_i,S_j]}}\right\} \right) f(S_j) {\Delta S_{j-1} + \Delta S_j
       \over 2} \right].
    \end{equation}
    In the case of constant $\Delta S_j(=\Delta S)$ this can be simplified further:
    \begin{equation} \label{eq:DesFinal3}
       f(S_i) = {2 \over \Delta S}\left[1 - \hbox{erf}\left\{\frac{\Delta \delta [B(S_i),\delta_1,S_i,S_1]}{\sqrt{2\Delta S [S_i,S_1]}}\right\}\right] - 2 \sum_{j=0}^{i-1}
       \left(1- \hbox{erf}\left\{\frac{\Delta \delta [B(S_i), B(S_j),S_i,S_j]}{\sqrt{2 \Delta S[S_i,S_j]}}\right\} \right) f(S_j).
    \end{equation}
    
    In either case (i.e. eqns.~\ref{eq:DesFinal2a} and \ref{eq:DesFinal3}) solution proceeds recursively: $f(S_0)=0$ by
    definition, $f(S_1)$ depends only on the known barrier and $f(S_0)$, $f(S_i)$ depends only on the known barrier and
    $f(S_{&lt;i})$.
    
    The first crossing rate is computed using the same method but with an effective barrier which is offset by the position of the
    progenitor in the $(\delta,S)$ plane, plus a small shift in time. The non-crossing rate, $g(S_\mathrm{max})$, defined as the
    rate at which trajectories reach the maximum possible variance, $S_\mathrm{max}$, without ever crossing the barrier---is
    computed directly by integrating over the first crossing rate distribution, i.e. $g(S_\mathrm{max}) = 1
    -\int_0^{S_\mathrm{max}} f(\tilde{S}) \mathrm{d}\tilde{S}$. Note that since the numerical integration occurs only up to a
    finite maximum $S_\mathrm{max}$, a non-zero non-crossing rate will be computed for CDM-like barriers even though in reality
    they should have zero non-crossing rate. As such, use of this method for such barriers is not recommended.
   </description>
  </excursionSetFirstCrossing>
  !!]
  type, extends(excursionSetFirstCrossingClass) :: excursionSetFirstCrossingFarahi
     !!{
     An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer                       :: cosmologyFunctions_              => null()
     class           (excursionSetBarrierClass     ), pointer                       :: excursionSetBarrier_             => null()
     class           (cosmologicalMassVarianceClass), pointer                       :: cosmologicalMassVariance_        => null()
     ! Variables used in tabulation the first crossing function.
     double precision                                                               :: timeMaximum                                , timeMinimum                               , &
          &                                                                            varianceMaximum
     integer                                                                        :: countTime                                  , countVariance
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
     <methods>
       <method description="Interpolate in the tabulated excursion set barrier crossing rates."                                            method="rateInterpolate"           />
       <method description="Interpolate in the tabulated excursion set barrier non-crossing rates."                                        method="rateNonCrossingInterpolate"/>
       <method description="Tabulate excursion set barrier crossing rates ensuring that they span the given progenitor variance and time." method="rateTabulate"              />
       <method description="Build a range of variances at which to tabulate the excursion set solutions."                                  method="varianceRange"             />
       <method description="Return the maximum variance to which to tabulate."                                                             method="varianceLimit"             />
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
     procedure :: varianceResidual           => farahiVarianceResidual
     procedure :: offsetEffective            => farahiOffsetEffective
  end type excursionSetFirstCrossingFarahi

  interface excursionSetFirstCrossingFarahi
     !!{
     Constructors for the Farahi excursion set barrier class.
     !!}
     module procedure farahiConstructorParameters
     module procedure farahiConstructorInternal
  end interface excursionSetFirstCrossingFarahi

  ! Parameters controlling tabulation range
  double precision, parameter :: redshiftMaximum=30.0d0 , redshiftMinimum=0.0d0

contains

  function farahiConstructorParameters(parameters) result(self)
    !!{
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
    <inputParameter>
      <name>fileName</name>
      <defaultValue>var_str('none')</defaultValue>
      <source>parameters</source>
      <description>The name of the file to/from which tabulations of barrier first crossing probabilities should be written/read. If set to ``{\normalfont \ttfamily none}'' tables will not be stored.</description>
    </inputParameter>
    <inputParameter>
      <name>fractionalTimeStep</name>
      <defaultValue>0.01d0</defaultValue>
      <source>parameters</source>
      <description>The fractional time step used when computing barrier crossing rates (i.e. the step used in finite difference calculations).</description>
    </inputParameter>
    <inputParameter>
      <name>varianceNumberPerUnitProbability</name>
      <defaultValue>1000</defaultValue>
      <source>parameters</source>
      <description>The number of points to tabulate per unit variance for first crossing probabilities.</description>
    </inputParameter>
    <inputParameter>
      <name>varianceNumberPerUnit</name>
      <defaultValue>40</defaultValue>
      <source>parameters</source>
      <description>The number of points to tabulate per unit variance for first crossing rates.</description>
    </inputParameter>
    <inputParameter>
      <name>varianceNumberPerDecade</name>
      <defaultValue>400</defaultValue>
      <source>parameters</source>
      <description>The number of points to tabulate per decade of progenitor variance for first crossing rates.</description>
    </inputParameter>
    <inputParameter>
      <name>varianceNumberPerDecadeNonCrossing</name>
      <defaultValue>40</defaultValue>
      <source>parameters</source>
      <description>The number of points to tabulate per decade of progenitor variance for non-crossing rates.</description>
    </inputParameter>
    <inputParameter>
      <name>timeNumberPerDecade</name>
      <defaultValue>10</defaultValue>
      <source>parameters</source>
      <description>The number of points to tabulate per decade of time.</description>
    </inputParameter>
    <inputParameter>
      <name>varianceIsUnlimited</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, the variance is assumed to have no upper limit (e.g. as in the case of \gls{cdm}). This allows the tabulated solutions to be extended arbitrarily. Otherwise, tables are extended to encompass just the range of variance requested.</description>
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
    !!{
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
    use :: File_Utilities, only : Directory_Make, File_Path
    use :: Input_Paths   , only : inputPath     , pathTypeDataDynamic
    implicit none
    class(excursionSetFirstCrossingFarahi), intent(inout) :: self

    if (self%fileNameInitialized) return
    ! Build an automatic file name based on the descriptor for this object.
    if (self%fileName == "auto") &
         & self%fileName=inputPath(pathTypeDataDynamic)//'largeScaleStructure/excursionSets/'//self%objectType()//'_'//self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)//'.hdf5'
    ! Ensure directory exists.
    call Directory_Make(File_Path(self%fileName))
    self%fileNameInitialized=.true.
    return
  end subroutine farahiFileNameInitialize

  subroutine farahiDestructor(self)
    !!{
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
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    use :: Display         , only : displayCounter              , displayCounterClear  , displayIndent       , displayMessage, &
          &                         displayUnindent             , verbosityLevelWorking
    use :: Error_Functions , only : Error_Function_Complementary
    use :: File_Utilities  , only : File_Lock                   , File_Unlock          , lockDescriptor
    use :: Kind_Numbers    , only : kind_dble                   , kind_quad
    use :: MPI_Utilities   , only : mpiBarrier                  , mpiSelf
    use :: Numerical_Ranges, only : Make_Range                  , rangeTypeLinear      , rangeTypeLogarithmic
    use :: Table_Labels    , only : extrapolationTypeFix
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)                 :: self
    double precision                                 , intent(in   )                 :: variance                        , time
    type            (treeNode                       ), intent(inout)                 :: node
    double precision                                 ,                dimension(0:1) :: hTime                           , hVariance
    double precision                                 , parameter                     :: toleranceRelativeVariance=1.0d-6
    class           (excursionSetBarrierClass       ), pointer                       :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass  ), pointer                       :: cosmologicalMassVariance_
    double precision                                 , allocatable  , dimension( : ) :: barrier
    double precision                                                                 :: barrierTest
    logical                                                                          :: makeTable
    integer         (c_size_t                       )                                :: iTime                           , iVariance       , &
         &                                                                              loopCount                       , loopCountTotal  , &
         &                                                                              i                               , j               , &
         &                                                                              jTime                           , jVariance
    double precision                                                                 :: probabilityCrossingPrior
    real            (kind_quad                      )                                :: offsetEffective                 , varianceResidual
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
    if (self%useFile.and..not.self%tableInitialized) then
       call self%fileNameInitialize()
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
       if (self%useFile) then
          call self%fileNameInitialize()
          call File_Lock(self%fileName,fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+toleranceRelativeVariance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
       if (makeTable) then
          ! Construct the table of variance on which we will solve for the first crossing distribution.
          if (allocated(self%variance                )) deallocate(self%variance                )
          if (allocated(self%time                    )) deallocate(self%time                    )
          if (allocated(self%firstCrossingProbability)) deallocate(self%firstCrossingProbability)
          self%varianceMaximum   =max(self%varianceMaximum,variance)
          self%countVariance=int(self%varianceMaximum*dble(self%varianceNumberPerUnitProbability))
          if (self%tableInitialized) then
             self%timeMinimum=min(      self%timeMinimum                                          ,0.5d0*time)
             self%timeMaximum=max(      self%timeMaximum                                          ,2.0d0*time)
          else
             self%timeMinimum=                                                                     0.5d0*time
             self%timeMaximum=max(2.0d0*self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0),2.0d0*time)
          end if
          self%countTime=max(2,int(log10(self%timeMaximum/self%timeMinimum)*dble(self%timeNumberPerDecade))+1)
          allocate(self%variance                (0:self%countVariance               ))
          allocate(self%time                    (                     self%countTime))
          allocate(self%firstCrossingProbability(0:self%countVariance,self%countTime))
          self%time        =Make_Range(self%timeMinimum,self%timeMaximum    ,self%countTime      ,rangeType=rangeTypeLogarithmic)
          self%variance    =Make_Range(0.0d0           ,self%varianceMaximum,self%countVariance+1,rangeType=rangeTypeLinear     )
          self%varianceStep=+self%variance(1) &
               &            -self%variance(0)
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
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range in which we are interested.
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMinimum,node,rateCompute=.false.)
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMaximum,node,rateCompute=.false.)
          ! Enter an OpenMP parallel region. Each thread will solve for the first crossing distribution at a different epoch.
          !$omp parallel private(iTime,i,j,probabilityCrossingPrior,excursionSetBarrier_,cosmologicalMassVariance_,barrier,offsetEffective,varianceResidual) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
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
          allocate(barrier(0:self%countVariance))
          !$omp do schedule(dynamic)
          do iTime=1,self%countTime
#ifdef USEMPI
             if (self%coordinatedMPI_ .and. mod(iTime-1,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
             ! Construct a table of barrier values as a function of variance.
             do i=0,self%countVariance
                barrier(i)=excursionSetBarrier_%barrier(self%variance(i),self%time(iTime),node,rateCompute=.false.)
             end do
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
             ! Iterate over variance, computing the first crossing distribution at each value.
             do i=2,self%countVariance
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
                call File_Unlock(fileLock)
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
    !!{
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
    !!{
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
    !!{
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
    !!{
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
    !!{
    Tabulate the excursion set crossing rate.
    !!}
    use :: Display         , only : displayCounter              , displayCounterClear  , displayIndent       , displayMessage, &
          &                         displayUnindent             , verbosityLevelWorking
    use :: Error_Functions , only : Error_Function_Complementary
    use :: File_Utilities  , only : File_Lock                   , File_Unlock          , lockDescriptor
    use :: Kind_Numbers    , only : kind_dble                   , kind_quad
    use :: MPI_Utilities   , only : mpiBarrier                  , mpiSelf
    use :: Numerical_Ranges, only : Make_Range                  , rangeTypeLinear      , rangeTypeLogarithmic
    use :: Table_Labels    , only : extrapolationTypeFix
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)               :: self
    double precision                                 , intent(in   )               :: time                             , varianceProgenitor
    type            (treeNode                       ), intent(inout)               :: node
    double precision                                 , parameter                   :: varianceMinimumDefault    =1.0d-2
    double precision                                 , parameter                   :: varianceTolerance         =1.0d-6
    double precision                                 , parameter                   :: massLarge                 =1.0d16
    real            (kind=kind_quad                 ), allocatable  , dimension(:) :: firstCrossingRateQuad            , varianceCurrentRateQuad   , &
         &                                                                            varianceProgenitorRateQuad       , barrierRateQuad
    double precision                                                               :: barrierRateTest
    class           (excursionSetBarrierClass       ), pointer                     :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass  ), pointer                     :: cosmologicalMassVariance_
#ifdef USEMPI
    integer                                                                        :: taskCount
#endif
    logical                                                                        :: makeTable
    integer         (c_size_t                       )                              :: loopCount                        , loopCountTotal
    integer                                                                        :: i                                , iTime                     , &
         &                                                                            iVariance                        , j                         , &
         &                                                                            iCompute                         , countVarianceCurrentRate
    double precision                                                               :: timeProgenitor                   , varianceMinimumRate       , &
         &                                                                            massProgenitor                   , varianceMaximumRateLimit
    character       (len=64                         )                              :: label
    type            (varying_string                 )                              :: message
    type            (lockDescriptor                 )                              :: fileLock
    real            (kind=kind_quad                 )                              :: crossingFraction                 , barrierProgenitorEffective, &
         &                                                                            probabilityCrossingPrior         , varianceStepRate          , &
         &                                                                            barrier                          , growthFactorEffective     , &
         &                                                                            varianceResidual                 , offsetEffective

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
    if (self%useFile.and..not.self%tableInitializedRate) then
       call self%fileNameInitialize()
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
       if (self%useFile) then
          call File_Lock(self%fileName,fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
       if (makeTable.or.self%retabulateRateNonCrossing) then
          if (makeTable) then
             if (allocated(self%varianceProgenitorRate        )) deallocate(self%varianceProgenitorRate        )
             if (allocated(self%varianceCurrentRate           )) deallocate(self%varianceCurrentRate           )
             if (allocated(self%varianceCurrentRateNonCrossing)) deallocate(self%varianceCurrentRateNonCrossing)
             if (allocated(self%timeRate                      )) deallocate(self%timeRate                      )
             if (allocated(self%firstCrossingRate             )) deallocate(self%firstCrossingRate             )
             if (self%tableInitializedRate) then
                self%timeMinimumRate=min(self%timeMinimumRate,0.5d0*time)
                self%timeMaximumRate=max(self%timeMaximumRate,2.0d0*time)
                self%countTimeRate  =int(log10(self%timeMaximumRate/self%timeMinimumRate)*dble(self%timeNumberPerDecade))+1
             else
                self%timeMinimumRate=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
                self%timeMaximumRate=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
                self%timeMinimumRate=min(self%timeMinimumRate,0.5d0*time)
                self%timeMaximumRate=max(self%timeMaximumRate,2.0d0*time)
                self%countTimeRate  =max(int(log10(self%timeMaximumRate/self%timeMinimumRate)*dble(self%timeNumberPerDecade))+1,2)
             end if
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
             self%varianceMaximumRate                =self%varianceLimit(varianceProgenitor)
             self%countVarianceProgenitorRate        =int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecade           ))+1
             self%countVarianceCurrentRate           =int(self%varianceMaximumRate*dble(self%varianceNumberPerUnit))
             self%countVarianceCurrentRateNonCrossing=int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecadeNonCrossing))+1
             allocate(self%varianceProgenitorRate        (0:self%countVarianceProgenitorRate                                                              ))
             allocate(self%varianceCurrentRate           (                                   0:self%countVarianceCurrentRate                              ))
             allocate(self%varianceCurrentRateNonCrossing(                                   0:self%countVarianceCurrentRateNonCrossing                   ))
             allocate(self%timeRate                      (                                                                              self%countTimeRate))
             allocate(self%firstCrossingRate             (0:self%countVarianceProgenitorRate,0:self%countVarianceCurrentRate           ,self%countTimeRate))
             ! For the variance table, the zeroth point is always zero, higher points are distributed uniformly in variance.
             self%varianceProgenitorRate        (0                                         )=0.0d0
             self%varianceProgenitorRate        (1:self%countVarianceProgenitorRate        )=self%varianceRange(varianceMinimumRate ,self%varianceMaximumRate,self%countVarianceProgenitorRate          ,exponent =1.0d0              )
             self%varianceCurrentRate           (0:self%countVarianceCurrentRate           )=Make_Range        (0.0d0               ,self%varianceMaximumRate,self%countVarianceCurrentRate           +1,rangeType=rangeTypeLinear    )
             self%varianceCurrentRateNonCrossing(0                                         )=0.0d0
             self%varianceCurrentRateNonCrossing(1:self%countVarianceCurrentRateNonCrossing)=self%varianceRange(varianceMinimumRate ,self%varianceMaximumRate,self%countVarianceCurrentRateNonCrossing  ,exponent =1.0d0              )
             ! The time table is logarithmically distributed in time.
             self%timeRate                                                                  =Make_Range        (self%timeMinimumRate,self%timeMaximumRate   ,self%countTimeRate                        ,rangeType=rangeTypeLogarithmic)
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
          loopCountTotal   =+int(self%countTimeRate,kind=c_size_t)*int(self%countVarianceCurrentRateNonCrossing+1,kind=c_size_t)
          if (makeTable) then
             loopCountTotal=+loopCountTotal                                                                                      &
                  &         +int(self%countTimeRate,kind=c_size_t)*int(self%countVarianceCurrentRate           +1,kind=c_size_t)
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
                   ! For computing non-crossing rates, the results are tabulated with respect to $S_{\rm max}-S$ so that interpolation
                   ! is more accurate when $S$ approaches $S_{\rm max}$.
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
                call File_Unlock(fileLock)
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
    !!{
    Read tabulated data on excursion set first crossing probabilities from file.
    !!}
    use :: Display           , only : displayIndent       , displayMessage  , displayUnindent, verbosityLevelWorking
    use :: File_Utilities    , only : File_Exists         , File_Name_Expand
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : operator(//)        , var_str         , varying_string
    use :: String_Handling   , only : operator(//)
    use :: Table_Labels      , only : extrapolationTypeFix
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)                   :: self
    type            (hdf5Object                     )                                  :: dataFile                     , dataGroup
    double precision                                 , allocatable  , dimension(:    ) :: varianceCurrentTmp           , varianceTmp    , &
         &                                                                                varianceCurrentNonCrossingTmp
    double precision                                 , allocatable  , dimension(:,:  ) :: firstCrossingProbabilityTmp  , nonCrossingRate
    double precision                                 , allocatable  , dimension(:,:,:) :: firstCrossingRateTmp
    double precision                                                                   :: massMinimumRateNonCrossing
    type            (varying_string                 )                                  :: message
    character       (len=32                         )                                  :: label

    ! Initialize file name.
    call self%fileNameInitialize()
    ! Return immediately if the file does not exist.
    if (.not.File_Exists(self%fileName)) return
    ! Open the data file.
    !$ call hdf5Access%set()
    call dataFile%openFile(self%fileName)
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
       call dataGroup%close()
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
       ! Set table limits.
       self%timeMinimum     =+self%time    (                 1)
       self%timeMaximum     =+self%time    (self%    countTime)
       self%varianceMaximum =+self%variance(self%countVariance)
       self%varianceStep    =+self%variance(                 1) &
            &                -self%variance(                 0)
       self%tableInitialized=.true.
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
       call dataGroup%close()
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
       ! Set table limits.
       self%varianceMaximumRate =self%varianceProgenitorRate(self%countVarianceProgenitorRate)
       self%timeMinimumRate     =self%timeRate              (                               1)
       self%timeMaximumRate     =self%timeRate              (self%countTimeRate              )
       self%tableInitializedRate=.true.
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
    ! Close the data file.
    call dataFile%close()
    !$ call hdf5Access%unset()
    return
  end subroutine farahiFileRead

  subroutine farahiFileWrite(self)
    !!{
    Write tabulated data on excursion set first crossing probabilities to file.
    !!}
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: Display           , only : displayIndent, displayMessage, displayUnindent, verbosityLevelWorking
    use :: HDF5              , only : hsize_t
    use :: ISO_Varying_String, only : operator(//) , var_str       , varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    class    (excursionSetFirstCrossingFarahi), intent(inout) :: self
    type     (hdf5Object                     )                :: dataFile, dataGroup
    type     (varying_string                 )                :: message
    character(len=32                         )                :: label

    ! Don't write anything if neither table is initialized.
    if (.not.(self%tableInitialized.or.self%tableInitializedRate)) return
    ! Initialize file name.
    call self%fileNameInitialize()
    ! Open the data file.
    !$ call hdf5Access%set()
    call dataFile%openFile(self%fileName,overWrite=.true.,chunkSize=100_hsize_t,compressionLevel=9)
    ! Check if the standard table is populated.
    if (self%tableInitialized) then
       dataGroup=dataFile%openGroup("probability")
       call dataGroup%writeDataset(self%variance                ,'variance'                ,'The variance at which results are tabulated.'                         )
       call dataGroup%writeDataset(self%time                    ,'time'                    ,'The cosmic times at which results are tabulated.'                     )
       call dataGroup%writeDataset(self%firstCrossingProbability,'firstCrossingProbability','The probability of first crossing as a function of variance and time.')
       call dataGroup%close()
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
       call dataGroup%close()
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
    ! Close the data file.
    call dataFile%close()
    !$ call hdf5Access%unset()
    return
  end subroutine farahiFileWrite

  function farahiVarianceRange(self,rangeMinimum,rangeMaximum,rangeNumber,exponent) result (rangeValues)
    !!{
    Builds a numerical range between {\normalfont \ttfamily rangeMinimum} and {\normalfont \ttfamily rangeMaximum} using
    {\normalfont \ttfamily rangeNumber} points with spacing that varies from logarithmic to linear spacing with the transition
    point controlled by {\normalfont \ttfamily exponent}. Specifically, suppose we have $N=${\normalfont \ttfamily rangeNumber}
    points in the range, from $S_\mathrm{min}=${\normalfont \ttfamily rangeMinimum} to $S_\mathrm{max}=${\normalfont \ttfamily
    rangeMaximum}. We define $f_i=(i-1)/(N-1)$ where $i$ runs from $1$ to $N$. We then define:
    \begin{equation}
     f_i = { \int_{S_\mathrm{min}}^{S_i} x^{n_i} \mathrm{d} x \over \int_{S_\mathrm{min}}^{S_\mathrm{max}} x^{n_i} \mathrm{d} x},
    \end{equation}
    and solve for $S_i$ to find the $i^\mathrm{th}$ range value. If $n_i=0$ then this will give $S_i$ linearly spaced between
    $S_\mathrm{min}$ and $S_\mathrm{max}$, while if $n_i=-1$ this will give $S_i$ logarithmically spaced between
    $S_\mathrm{min}$ and $S_\mathrm{max}$. Therefore, if we make $n_i$ vary from $-1$ to $0$ at $i$ ranges from $1$ to $N$ we
    will get a smooth transition from logarithmic to linear spacing. We choose to use $n_i=-1+f_i^\alpha$ where
    $\alpha=${\normalfont \ttfamily exponent} is a supplied argument.
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
    !!{
    Return the maximum variance to which to tabulate.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout) :: self
    double precision                                 , intent(in   ) :: varianceProgenitor

    if (self%varianceIsUnlimited) then
       farahiVarianceLimit=max(self%varianceMaximumRate,2.0d0*varianceProgenitor)
    else
       farahiVarianceLimit=max(self%varianceMaximumRate,      varianceProgenitor)
    end if
    return
  end function farahiVarianceLimit

  function farahiVarianceResidual(self,time,varianceCurrent,varianceProgenitor,varianceIntermediate,cosmologicalMassVariance_) result(varianceResidual)
    !!{
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
    !!{
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
