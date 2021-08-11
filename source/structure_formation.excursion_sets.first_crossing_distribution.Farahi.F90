!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}.
!!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Excursion_Sets_Barriers   , only : excursionSetBarrierClass
  use :: Numerical_Interpolation   , only : interpolator

  !![
  <excursionSetFirstCrossing name="excursionSetFirstCrossingFarahi">
   <description>
    An excursion set first crossing statistics class using the algorithm of \cite{benson_dark_2012}, which proceeds by finding
    the solution to the integral equation:
    \begin{equation}
      1 =  \int_0^S f(S^\prime)\mathrm{d}S^\prime + \int_{-\infty}^{B(S)} P(\delta,S) \mathrm{d} \delta,
     \label{eq:OldExcursionMethod}
    \end{equation}
    where $P(\delta,S) \mathrm{d} \delta$ is the probability for a trajectory to lie between $\delta$ and $\delta + \mathrm{d}
    \delta$ at variance $S$. In the absence of a barrier, $P(\delta,S)$ would be equal to $P_0(\delta,S)$ which is simply a
    Gaussian distribution with variance $S$:
    \begin{equation}
      P_0(\delta,S) = \frac{1}{\sqrt{2 \pi S}} \exp\left(-{\delta^2 \over 2 S}\right).
      \label{eq:Gaussian}
    \end{equation}
    Since the barrier absorbs any random walks which cross is at smaller $S$, the actual $P(\delta,S)$ must therefore be given
    by:
    \begin{equation}
       P(\delta,S) = P_0(\delta,S) - \int_{0}^{S} f(S^\prime) P_0[\delta - B(S^\prime),S - S^\prime]\mathrm{d}S^\prime .
     \label{eq:Displaced}
    \end{equation}
    In the second term on the right hand side of eqn.~(\ref{eq:Displaced}) represents the $P_0[\delta - B(S^\prime),S -
    S^\prime]$ term represents the distribution of random trajectories orginating from the point $(S,B(S))$. The integral
    therefore gives the fraction of trajectories which crossed the barrier at $S&lt;S^\prime$ and which can now be found at
    $(S,\delta)$.
    
    Using this result, we can rewrite eqn.~(\ref{eq:OldExcursionMethod}):
    \begin{equation}
      1 = \int_0^S f(S^\prime)\mathrm{d}S^\prime + \int_{-\infty}^{B(S)} \left[ P_0(\delta,S) - \int_{0}^{S} f(S^\prime)
      P_0(\delta - B(S^\prime),S - S^\prime)\mathrm{d}S^\prime )\right] \mathrm{d} \delta ,
    \end{equation}
    in general and, for the Gaussian distribution of eqn.~(\ref{eq:Gaussian}):
    \begin{equation}
      1 = \int_0^S f(S^\prime)\mathrm{d}S^\prime + \int_{-\infty}^{B(S)} \left[ \frac{1}{\sqrt{2 \pi S}}
      \exp\left(-\frac{\delta^2}{2 S}\right) - \int_{0}^{S} f(S^\prime) \frac{1}{\sqrt{2 \pi (S-S^\prime)}}
      \exp\left(-\frac{[\delta - B(S^\prime)]^2}{2 (S-S^\prime)}\right)\mathrm{d}S^\prime \right] \mathrm{d} \delta .
    \end{equation}
    The integral over $\mathrm{d}\delta$ can be carried out analytically to give:
    \begin{equation}
     1 = \int_0^S f(S^\prime)\mathrm{d}S^\prime+ \hbox{erf}\left[\frac{B(S)}{\sqrt{2S}}\right] - \int_{0}^{S} f(S^\prime)
     \hbox{erf}\left[\frac{B(S) - B(S^\prime)}{\sqrt{2 (S-S^\prime)}}\right] \mathrm{d}S^{\prime\prime}.
    \label{eq:NewExcursionMethod}
    \end{equation}
    We now discretize eqn.~(\ref{eq:NewExcursionMethod}). Specifically, we divide the $S$ space into $N$ intervals defined by
    the points:
    \begin{equation}
      S_i = \left\{ \begin{array}{ll}
                     0 &amp; \hbox{if } i=0 \\
                     \sum_0^{i-1} \Delta S_i &amp; \hbox{if } i &gt; 1.
                    \end{array}
            \right.
    \end{equation}
    Note that $f(0)=0$ by definition, so $f(S_0)=0$ always. We choose $\Delta S_i = S_\mathrm{max}/N$ (i.e. uniform spacing in
    $S$) when computing first crossing distributions, and $\Delta S_i \propto S_i$ (i.e. uniform spacing in $\log(S)$) when
    computing first crossing rates.
    
    Discretizing the integrals in eqn.~(\ref{eq:NewExcursionMethod}) gives:
    \begin{equation} \label{eq:Des1}
     \int_0^{S_j} f(S^\prime)\d S^\prime = \sum_{i=0}^{j-1} \frac{f(S_i) + f(S_{i+1})}{2} \Delta S_i
    \end{equation}
    and:
    \begin{equation} \label{eq:Des2}
     \int_{0}^{S_j} f(S^\prime) \hbox{erf}\left[\frac{B(S) - B(S^\prime)}{\sqrt{2 (S-S^\prime)}}\right] \d S^\prime =
     \sum_{i=0}^{j-1} \frac{1}{2} \left(f(S_i) \hbox{erf}\left[\frac{B(S_j) - B(S_i)}{\sqrt{2 (S_j-S_i)}}\right] + f(S_{i+1})
     \hbox{erf}\left[\frac{B(S_j) - B(S_{i+1})}{\sqrt{2 (S_j-S_{i+1})}}\right] \right) \Delta S_i.
    \end{equation}
    We can now rewrite eqn.~(\ref{eq:NewExcursionMethod}) in discretized form:
    \begin{equation} \label{eq:DesFinal1}
     1 = \sum_{i=0}^{j-1} \frac{f(S_i) + f(S_{i+1})}{2} \Delta S_i + \hbox{erf}\left[\frac{B(S_j)}{\sqrt{2S_j}}\right] -
     \frac{1}{2} \sum_{i=0}^{j-1} \left( f(S_i) \hbox{erf}\left[\frac{B(S_j) - B(S_i)}{\sqrt{2 (S_j-S_i)}}\right] + f(S_{i+1})
     \hbox{erf}\left[\frac{B(S_j) - B(S_{i+1})}{\sqrt{2 (S_j-S_{i+1})}}\right] \right) \Delta S_i.
    \end{equation}
    Solving eqn.~(\ref{eq:DesFinal1}) for $f(S_j)$:
    \begin{eqnarray} \label{eq:DesFinal11}
     \left( \frac{1}{2} - \frac{1}{2} \hbox{erf}\left[\frac{B(S_j) - B(S_j)}{\sqrt{2 (S_j-S_j)}}\right] \right) \Delta S_{j-1}
     f(S_j) &amp;=&amp; 1 - \sum_{i=0}^{j-2} \frac{f(S_i) + f(S_{i+1})}{2} \Delta S_i - \frac{f(S_{j-1})}{2} \Delta S_{j-1} -
     \hbox{erf}\left\{\frac{B(S_j)}{\sqrt{2S_j}}\right\} \nonumber\\
    &amp; &amp; + \frac{1}{2} \sum_{i=0}^{j-2} \left( f(S_i) \hbox{erf}\left\{\frac{[B(S_j) - B(S_i)]}{\sqrt{2 (S_j-S_i)}}\right\} +
    f(S_{i+1}) \hbox{erf}\left\{\frac{[B(S_j) - B(S_{i+1})]}{\sqrt{2 (S_j-S_{i+1})}}\right\} \right)\Delta S_i \nonumber \\
     &amp; &amp; + \frac{1}{2} f(S_{j-1}) \hbox{erf}\left\{\frac{[B(S_j) - B(S_{j-1})]}{\sqrt{2 (S_j-S_{j-1})}}\right\} \Delta S_{j-1}.
    \end{eqnarray}
    For all barriers that we consider:
    \begin{equation} 
    \hbox{erf}\left[\frac{B(S_j) - B(S_j)}{\sqrt{2 (S_j-S_j)}}\right] = 0.
    \end{equation}
    We can then simplify eqn.~(\ref{eq:DesFinal11}):
    \begin{eqnarray} \label{eq:DesFinal2}
       f(S_j) &amp;=&amp; {2 \over \Delta S_{j-1}}\left[1 - \sum_{i=0}^{j-2} \frac{f(S_i) + f(S_{i+1})}{2} \Delta S_i -
       \frac{f(S_{j-1})}{2} \Delta S_{j-1} - \hbox{erf}\left\{\frac{B(S_j)}{\sqrt{2S_j}}\right\} \right.  \nonumber\\
    &amp; &amp; + \frac{1}{2} \sum_{i=0}^{j-2} \left( f(S_i) \hbox{erf}\left\{\frac{[B(S_j) - B(S_i)]}{\sqrt{2 (S_j-S_i)}}\right\} +
    f(S_{i+1}) \hbox{erf}\left\{\frac{[B(S_j) - B(S_{i+1})]}{\sqrt{2 (S_j-S_{i+1})}}\right\} \right)\Delta S_i \nonumber \\
     &amp; &amp; \left. + \frac{1}{2} f(S_{j-1}) \hbox{erf}\left\{\frac{[B(S_j) - B(S_{j-1})]}{\sqrt{2 (S_j-S_{j-1})}}\right\} \Delta
     S_{j-1}\right].
    \end{eqnarray}
    Consolidating terms in the summations:
    \begin{equation} \label{eq:DesFinal2a}
       f(S_j) = {2 \over \Delta S_{j-1}}\left[1 - \hbox{erf}\left\{\frac{B(S_j)}{\sqrt{2S_j}}\right\} - \sum_{i=0}^{j-1}
       \left\{ 1-\hbox{erf}\left[\frac{B(S_j) - B(S_i)}{\sqrt{2 (S_j-S_i)}}\right] \right\} f(S_i) {\Delta S_{i-1} + \Delta S_i
       \over 2} \right].
    \end{equation}
    In the case of constant $\Delta S_i(=\Delta S)$ this can be simplified further:
    \begin{equation} \label{eq:DesFinal3}
       f(S_j) = {2 \over \Delta S}\left[1 - \hbox{erf}\left\{\frac{B(S_j)}{\sqrt{2S_j}}\right\}\right] - 2 \sum_{i=0}^{j-1}
       \left\{1- \hbox{erf}\left[\frac{B(S_j) - B(S_i)}{\sqrt{2 (S_j-S_i)}}\right] \right\} f(S_i).
    \end{equation}
    
    In either case (i.e. eqns.~\ref{eq:DesFinal2a} and \ref{eq:DesFinal3}) solution proceeds recursively: $f(S_0)=0$ by
    definition, $f(S_1)$ depends only on the known barrier and $f(S_0)$, $f(S_j)$ depends only on the known barrier and
    $f(S_{&lt;j})$.
    
    The first crossing rate is computed using the same method but with an effective barrier which is offset by the position of
    the progenitor in the $(\delta,S)$ plane, plus a small shift in time. The non-crossing rate is computed directly by
    integrating over the first crossing rate distribution. Note that since the numerical integration occurs only up to a finite
    maximum $S$, a non-zero non-crossing rate will be computed for CDM-like barriers even though in reality they should have
    zero non-crossing rate. As such, use of this method for such barriers is not recommended.
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
     double precision                                                               :: timeMaximum                               , timeMinimum             , &
          &                                                                            varianceMaximum
     integer                                                                        :: timeTableCount                            , varianceTableCount
     double precision                               , allocatable, dimension(:,:)   :: firstCrossingProbabilityTable
     double precision                               , allocatable, dimension(:  )   :: timeTable                                 , varianceTable
     double precision                                                               :: varianceTableStep
     logical                                                                        :: tableInitialized                          , fileNameInitialized
     type            (interpolator                 ), allocatable                   :: interpolatorTime                          , interpolatorVariance
     ! Variables used in tabulation the first crossing rate function.
     double precision                                                               :: timeMaximumRate                           , timeMinimumRate         , &
          &                                                                            varianceMaximumRate
     integer                                                                        :: timeTableCountRate                        , varianceTableCountRate  , &
          &                                                                            varianceTableCountRateBase
     double precision                               , allocatable, dimension(:,:,:) :: firstCrossingTableRate
     double precision                               , allocatable, dimension(:,:  ) :: nonCrossingTableRate
     double precision                               , allocatable, dimension(:    ) :: timeTableRate                             , varianceTableRate       , &
          &                                                                            varianceTableRateBase
     logical                                                                        :: tableInitializedRate
     type            (interpolator                 ), allocatable                   :: interpolatorTimeRate                      , interpolatorVarianceRate, &
          &                                                                            interpolatorVarianceRateBase
     ! File name used to store tabulations.
     type            (varying_string               )                                :: fileName
     logical                                                                        :: useFile
     ! Tabulation resolutions.
     integer                                                                        :: varianceNumberPerUnitProbability          , varianceNumberPerUnit   , &
          &                                                                            timeNumberPerDecade                       , varianceNumberPerDecade
     ! The fractional step in time used to compute barrier crossing rates.
     double precision                                                               :: timeStepFractional
     ! Record of variance and time in previous call to rate functions.
     double precision                                                               :: timeRatePrevious                          , varianceRatePrevious
     double precision                                            , dimension(0:1)   :: hTimeRate                                 , hVarianceRate
     integer         (c_size_t                     )                                :: iTimeRate                                 , iVarianceRate
   contains
     !![
     <methods>
       <method description="Tabulate excursion set barrier crossing rates ensuring that they span the given progenitor variance and time." method="rateTabulate" />
       <method description="Build a range of variances at which to tabulate the excursion set solutions." method="varianceRange" />
       <method description="Read excursion set solutions from file." method="fileRead" />
       <method description="Write excursion set solutions to file." method="fileWrite" />
       <method description="Initialize the file name for storing excursion set data." method="fileNameInitialize" />
     </methods>
     !!]
     final     ::                       farahiDestructor
     procedure :: probability        => farahiProbability
     procedure :: rate               => farahiRate
     procedure :: rateNonCrossing    => farahiRateNonCrossing
     procedure :: rateTabulate       => farahiRateTabulate
     procedure :: fileRead           => farahiFileRead
     procedure :: fileWrite          => farahiFileWrite
     procedure :: fileNameInitialize => farahiFileNameInitialize
     procedure :: varianceRange      => farahiVarianceRange
  end type excursionSetFirstCrossingFarahi

  interface excursionSetFirstCrossingFarahi
     !!{
     Constructors for the Farahi excursion set barrier class.
     !!}
     module procedure farahiConstructorParameters
     module procedure farahiConstructorInternal
  end interface excursionSetFirstCrossingFarahi

  ! Parameters controlling tabulation range
  double precision                , parameter :: farahiRateRedshiftMaximum=30.0d0 , farahiRateRedshiftMinimum=0.0d0

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
    double precision                                                 :: timeStepFractional
    integer                                                          :: varianceNumberPerUnitProbability, varianceNumberPerUnit  , &
         &                                                              timeNumberPerDecade             , varianceNumberPerDecade
    type            (varying_string                 )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <defaultValue>var_str('none')</defaultValue>
      <source>parameters</source>
      <description>The name of the file to/from which tabulations of barrier first crossing probabilities should be written/read. If set to ``{\normalfont \ttfamily none}'' tables will not be stored.</description>
    </inputParameter>
    <inputParameter>
      <name>timeStepFractional</name>
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
      <name>timeNumberPerDecade</name>
      <defaultValue>10</defaultValue>
      <source>parameters</source>
      <description>The number of points to tabulate per decade of time.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="excursionSetBarrier"      name="excursionSetBarrier_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=excursionSetFirstCrossingFarahi(timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="excursionSetBarrier_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function farahiConstructorParameters

  function farahiConstructorInternal(timeStepFractional,fileName,varianceNumberPerUnitProbability,varianceNumberPerUnit,varianceNumberPerDecade,timeNumberPerDecade,cosmologyFunctions_,excursionSetBarrier_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the Farahi excursion set class first crossing class.
    !!}
    implicit none
    type            (excursionSetFirstCrossingFarahi)                        :: self
    double precision                                 , intent(in   )         :: timeStepFractional
    integer                                          , intent(in   )         :: varianceNumberPerUnitProbability, varianceNumberPerUnit  , &
         &                                                                      timeNumberPerDecade             , varianceNumberPerDecade
    type            (varying_string                 ), intent(in   )         :: fileName
    class           (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class           (excursionSetBarrierClass       ), intent(in   ), target :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass  ), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="timeStepFractional, fileName, varianceNumberPerUnitProbability, varianceNumberPerUnit, varianceNumberPerDecade, timeNumberPerDecade, *cosmologyFunctions_, *excursionSetBarrier_, *cosmologicalMassVariance_"/>
    !!]

    self%tableInitialized                  =.false.
    self%tableInitializedRate              =.false.
    self%timeMaximum                       =-huge(0.0d0)
    self%timeMinimum                       =+huge(0.0d0)
    self%varianceMaximum                   =      0.0d0
    self%timeMaximumRate                   =-huge(0.0d0)
    self%timeMinimumRate                   =+huge(0.0d0)
    self%varianceMaximumRate               =-huge(0.0d0)
    self%timeRatePrevious                  =-huge(0.0d0)
    self%varianceRatePrevious              =-huge(0.0d0)
    self%useFile                           =(self%fileName /= 'none')
    self%fileNameInitialized               =.not.self%useFile
    return
  end function farahiConstructorInternal

  subroutine farahiFileNameInitialize(self)
    use :: File_Utilities  , only : Directory_Make, File_Name_Expand   , File_Path
    use :: Galacticus_Paths, only : galacticusPath, pathTypeDataDynamic
    implicit none
    class(excursionSetFirstCrossingFarahi), intent(inout) :: self

    if (self%fileNameInitialized) return
    ! Build an automatic file name based on the descriptor for this object.
    if (self%fileName == "auto") &
         & self%fileName=galacticusPath(pathTypeDataDynamic)//'largeScaleStructure/excursionSets/'//self%objectType()//'_'//self%hashedDescriptor(includeSourceDigest=.true.)//'.hdf5'
    ! Expand file name.
    self%fileName=File_Name_Expand(char(self%fileName))
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
    use :: Display          , only : displayCounter , displayCounterClear  , displayIndent       , displayMessage, &
          &                          displayUnindent, verbosityLevelWorking
    use :: Error_Functions  , only : erfApproximate
    use :: File_Utilities   , only : File_Lock      , File_Unlock          , lockDescriptor
    use :: Kind_Numbers     , only : kind_dble      , kind_quad
    use :: MPI_Utilities    , only : mpiBarrier     , mpiSelf
    use :: Memory_Management, only : allocateArray  , deallocateArray
    use :: Numerical_Ranges , only : Make_Range     , rangeTypeLinear      , rangeTypeLogarithmic
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)                 :: self
    double precision                                 , intent(in   )                 :: variance                     , time
    type            (treeNode                       ), intent(inout)                 :: node
    double precision                                 ,                dimension(0:1) :: hTime                        , hVariance
    double precision                                 , parameter                     :: varianceTableTolerance=1.0d-6
    class           (excursionSetBarrierClass       ), pointer                       :: excursionSetBarrier_
    double precision                                 , allocatable  , dimension( : ) :: barrierTable
    double precision                                                                 :: barrierTest
    logical                                                                          :: makeTable
    integer         (c_size_t                       )                                :: iTime                        , iVariance     , &
         &                                                                              loopCount                    , loopCountTotal, &
         &                                                                              i                            , j             , &
         &                                                                              jTime                        , jVariance
    double precision                                                                 :: sigma1f
    character       (len=6                          )                                :: label
    type            (varying_string                 )                                :: message
    type            (lockDescriptor                 )                                :: fileLock

    ! Read tables from file if possible.
    if (self%useFile.and..not.self%tableInitialized) then
       call self%fileNameInitialize()
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
       call self%fileRead()
       call File_Unlock(fileLock)
    end if
    ! Construct the table if necessary.
    makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+varianceTableTolerance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
#ifdef USEMPI
    if (self%coordinatedMPI_) call mpiBarrier()
#endif
    if (makeTable) then
       !$omp critical(farahiProbabilityTabulate)
       ! Attempt to read the file again now that we are within the critical section. If another thread made the file while we were waiting we may be able to skip building the table.
       if (self%useFile) then
          call self%fileNameInitialize()
          call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitialized.or.(variance > self%varianceMaximum*(1.0d0+varianceTableTolerance)).or.(time < self%timeMinimum).or.(time > self%timeMaximum)
       if (makeTable) then
          ! Construct the table of variance on which we will solve for the first crossing distribution.
          if (allocated(self%varianceTable                )) call deallocateArray(self%varianceTable                )
          if (allocated(self%timeTable                    )) call deallocateArray(self%timeTable                    )
          if (allocated(self%firstCrossingProbabilityTable)) call deallocateArray(self%firstCrossingProbabilityTable)
          self%varianceMaximum   =max(self%varianceMaximum,variance)
          self%varianceTableCount=int(self%varianceMaximum*dble(self%varianceNumberPerUnitProbability))
          if (self%tableInitialized) then
             self%timeMinimum=min(      self%timeMinimum                                          ,0.5d0*time)
             self%timeMaximum=max(      self%timeMaximum                                          ,2.0d0*time)
          else
             self%timeMinimum=                                                                     0.5d0*time
             self%timeMaximum=max(2.0d0*self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0),2.0d0*time)
          end if
          self%timeTableCount=max(2,int(log10(self%timeMaximum/self%timeMinimum)*dble(self%timeNumberPerDecade))+1)
          call allocateArray(self%varianceTable                ,[1+self%varianceTableCount                    ],lowerBounds=[0  ])
          call allocateArray(self%timeTable                    ,[                          self%timeTableCount]                  )
          call allocateArray(self%firstCrossingProbabilityTable,[1+self%varianceTableCount,self%timeTableCount],lowerBounds=[0,1])
          self%timeTable        =Make_Range(self%timeMinimum,self%timeMaximum    ,self%timeTableCount      ,rangeType=rangeTypeLogarithmic)
          self%varianceTable    =Make_Range(0.0d0           ,self%varianceMaximum,self%varianceTableCount+1,rangeType=rangeTypeLinear     )
          self%varianceTableStep=self%varianceTable(1)-self%varianceTable(0)
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
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
             loopCountTotal=(int(self%timeTableCount,kind=c_size_t)/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t)*(int(self%varianceTableCount-1,kind=c_size_t)*int(self%varianceTableCount,kind=c_size_t))/2_c_size_t
          else
#endif
             loopCountTotal= int(self%timeTableCount,kind=c_size_t)                                               *(int(self%varianceTableCount-1,kind=c_size_t)*int(self%varianceTableCount,kind=c_size_t))/2_c_size_t
#ifdef USEMPI
          end if
#endif
          loopCount=0
#ifdef USEMPI
          if (self%coordinatedMPI_) self%firstCrossingProbabilityTable=0.0d0
#endif
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range we are intereseted in.
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMinimum,node,rateCompute=.false.)
          barrierTest=self%excursionSetBarrier_%barrier(self%varianceMaximum,self%timeMaximum,node,rateCompute=.false.)
          !$omp parallel private(iTime,i,j,sigma1f,excursionSetBarrier_,barrierTable) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
          allocate(excursionSetBarrier_,mold=self%excursionSetBarrier_)
          !$omp critical(excursionSetsSolverFarahiDeepCopy)
          !![
          <deepCopyReset variables="self%excursionSetBarrier_"/>
          <deepCopy source="self%excursionSetBarrier_" destination="excursionSetBarrier_"/>
          <deepCopyFinalize variables="excursionSetBarrier_"/>
          !!]
          !$omp end critical(excursionSetsSolverFarahiDeepCopy)
          call allocateArray(barrierTable,[1+self%varianceTableCount],lowerBounds=[0])
          !$omp do schedule(dynamic)
          do iTime=1,self%timeTableCount
#ifdef USEMPI
             if (self%coordinatedMPI_ .and. mod(iTime-1,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
             ! Construct the barrier table.
             do i=0,self%varianceTableCount
                barrierTable(i)=excursionSetBarrier_%barrier(self%varianceTable(i),self%timeTable(iTime),node,rateCompute=.false.)
             end do
             self%firstCrossingProbabilityTable(0,iTime)=0.0d0
             self%firstCrossingProbabilityTable(1,iTime)=                                                              &
                  &                                 real(                                                              &
                  &                                        +2.0_kind_quad                                              &
                  &                                      *(                                                            &
                  &                                        +1.0_kind_quad                                              &
                  &                                        -erfApproximate(                                            &
                  &                                                        +barrierTable(1)                            &
                  &                                                        /sqrt(2.0_kind_quad*self%varianceTable(1))  &
                  &                                                       )                                            &
                  &                                       )                                                            &
                  &                                      /self%varianceTableStep                                     , &
                  &                                      kind=kind_dble                                                &
                  &                                     )
             do i=2,self%varianceTableCount
#ifdef USEMPI
                if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
                   call displayCounter(int(100.0d0*dble(loopCount)/dble(loopCountTotal)),loopCount==0,verbosityLevelWorking)
#ifdef USEMPI
                end if
#endif
                !$omp atomic
                loopCount=loopCount+(i-1)
                sigma1f  =0.0d0
                do j=1,i-1
                   sigma1f=+sigma1f                                                                                  &
                        &  +self%firstCrossingProbabilityTable(j,iTime)                                              &
                        &  *real(                                                                                    &
                        &        +1.0_kind_quad                                                                      &
                        &        -erfApproximate(                                                                    &
                        &                        +(                                                                  &
                        &                           +barrierTable(i)                                                 &
                        &                           -barrierTable(j)                                                 &
                        &                         )                                                                  &
                        &                        /sqrt(2.0_kind_quad*(self%varianceTable(i)-self%varianceTable(j)))  &
                        &                       )                                                                  , &
                        &        kind=kind_dble                                                                      &
                        &       )
                end do
                self%firstCrossingProbabilityTable(i,iTime)=+real(                                                                   &
                     &                                            +max(                                                              &
                     &                                                 +0.0_kind_quad,                                               &
                     &                                                 +2.0_kind_quad                                                &
                     &                                                 *(                                                            &
                     &                                                   +1.0_kind_quad                                              &
                     &                                                   -erfApproximate(                                            &
                     &                                                                   +barrierTable(i)                            &
                     &                                                                   /sqrt(2.0_kind_quad*self%varianceTable(i))  &
                     &                                                                  )                                            &
                     &                                                  )                                                            &
                     &                                                 /self%varianceTableStep                                       &
                     &                                                 -2.0_kind_quad                                                &
                     &                                                 *sigma1f                                                      &
                     &                                                )                                                            , &
                     &                                            kind=kind_dble                                                     &
                     &                                           )
             end do
             ! Force the probability at maximum variance to zero.
             self%firstCrossingProbabilityTable(self%varianceTableCount,iTime)=0.0d0
          end do
          !$omp end do
          !![
          <objectDestructor name="excursionSetBarrier_"/>
          !!]
          call deallocateArray(barrierTable)
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
             self%firstCrossingProbabilityTable=mpiSelf%sum(self%firstCrossingProbabilityTable)
          end if
#endif
          ! Build the interpolators.
          if (allocated(self%interpolatorVariance)) deallocate(self%interpolatorVariance)
          if (allocated(self%interpolatorTime    )) deallocate(self%interpolatorTime    )
          allocate(self%interpolatorVariance)
          allocate(self%interpolatorTime    )
          self%interpolatorVariance=interpolator(self%varianceTable)
          self%interpolatorTime    =interpolator(self%timeTable    )
          ! Record that the table is now built.
          self%tableInitialized=.true.
          ! Write the table to file if possible.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             if (self%useFile) then
                call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
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
    call self%interpolatorTime%linearFactors    (time    ,iTime    ,hTime    )
    call self%interpolatorVariance%linearFactors(variance,iVariance,hVariance)
    ! Compute first crossing probability by interpolating.
    farahiProbability=0.0d0
    do jTime=0,1
       do jVariance=0,1
          farahiProbability=+farahiProbability                                                     &
               &            +hTime                             (                            jTime) &
               &            *hVariance                         (            jVariance            ) &
               &            *self%firstCrossingProbabilityTable(iVariance-1+jVariance,iTime+jTime)
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
    double precision                                 , intent(in   )  :: variance           , varianceProgenitor, &
         &                                                               time
    type            (treeNode                       ), intent(inout)  :: node
    double precision                                 , dimension(0:1) :: hVarianceProgenitor
    integer                                                           :: jVarianceProgenitor, jTime             , &
         &                                                               jVariance
    integer         (c_size_t                       )                 :: iVarianceProgenitor

    ! For progenitor variances less than or equal to the original variance, return zero.
    if (varianceProgenitor <= variance) then
       farahiRate=0.0d0
       return
    end if
    ! Ensure that the rate is tabulated.
    call self%rateTabulate(varianceProgenitor,time,node)
    ! For progenitor variances greater than the maximum allowed variance, return zero.
    if (varianceProgenitor > self%varianceMaximumRate) then
       farahiRate=0.0d0
       return
    end if
    ! Get interpolation in time.
    if (time /= self%timeRatePrevious) then
       self%timeRatePrevious    =time
       call self%interpolatorTimeRate        %linearFactors(time    ,self%iTimeRate    ,self%hTimeRate    )
    end if
    ! Get interpolation in variance.
    if (variance /= self%varianceRatePrevious) then
       self%varianceRatePrevious=variance
       call self%interpolatorVarianceRateBase%linearFactors(variance,self%iVarianceRate,self%hVarianceRate)
    end if
    ! Get interpolation in progenitor variance.
    iVarianceProgenitor=self%interpolatorVarianceRate%locate(varianceProgenitor-variance)
    ! Catch cases where the maximum variance is approached.
    if (self%varianceTableRate(iVarianceProgenitor)+variance > self%varianceMaximumRate) then
       ! Force the rate to drop to zero at the maximum variance. (Necessary because we will not have a tabulated point precisely
       ! at the maximum variance.)
       hVarianceProgenitor=[                                                                                      &
            &               +1.0d0                                                                                &
            &               -((     varianceProgenitor -variance)-self%varianceTableRate(iVarianceProgenitor-1))  &
            &               /((self%varianceMaximumRate-variance)-self%varianceTableRate(iVarianceProgenitor-1)), &
            &               +0.0d0                                                                                &
            &              ]
    else
       call self%interpolatorVarianceRate%linearWeights(varianceProgenitor-variance,iVarianceProgenitor,hVarianceProgenitor)
    end if
    ! Compute first crossing probability by interpolating.
    farahiRate=0.0d0
    do jTime=0,1
       do jVariance=0,1
          do jVarianceProgenitor=0,1
             farahiRate=+farahiRate                                                                                                                 &
                  &     +self%hTimeRate             (                                                                                        jTime) &
                  &     *self%hVarianceRate         (                                                               jVariance                     ) &
                  &     *hVarianceProgenitor        (                      jVarianceProgenitor                                                    ) &
                  &     *self%firstCrossingTableRate(iVarianceProgenitor-1+jVarianceProgenitor,self%iVarianceRate-1+jVariance,self%iTimeRate+jTime)
          end do
       end do
    end do
    return
  end function farahiRate

  double precision function farahiRateNonCrossing(self,variance,time,node)
    !!{
    Return the rate for excursion set non-crossing.
    !!}
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout) :: self
    double precision                                 , intent(in   ) :: time , variance
    type            (treeNode                       ), intent(inout) :: node
    integer                                                          :: jTime, jVariance

    ! Ensure that the rate is tabulated.
    call self%rateTabulate(variance,time,node)
    ! Get interpolation in time.
    if (time /= self%timeRatePrevious) then
       self%timeRatePrevious    =time
       call self%interpolatorTimeRate        %linearFactors(time    ,self%iTimeRate    ,self%hTimeRate    )
    end if
    ! Get interpolation in variance.
    if (variance /= self%varianceRatePrevious) then
       self%varianceRatePrevious=variance
       call self%interpolatorVarianceRateBase%linearFactors(variance,self%iVarianceRate,self%hVarianceRate)
    end if
    ! Compute non-crossing probability by interpolating.
    farahiRateNonCrossing=0.0d0
    do jTime=0,1
       do jVariance=0,1
          farahiRateNonCrossing=+farahiRateNonCrossing                                                          &
               &                +self%hTimeRate           (                                              jTime) &
               &                *self%hVarianceRate       (                     jVariance                     ) &
               &                *self%nonCrossingTableRate(self%iVarianceRate-1+jVariance,self%iTimeRate+jTime)
       end do
    end do
    return
  end function farahiRateNonCrossing

  subroutine farahiRateTabulate(self,varianceProgenitor,time,node)
    !!{
    Tabulate the excursion set crossing rate.
    !!}
    use :: Display          , only : displayCounter , displayCounterClear  , displayIndent       , displayMessage, &
          &                          displayUnindent, verbosityLevelWorking
    use :: Error_Functions  , only : erfApproximate
    use :: File_Utilities   , only : File_Lock      , File_Unlock          , lockDescriptor
    use :: Kind_Numbers     , only : kind_dble      , kind_quad
    use :: MPI_Utilities    , only : mpiBarrier     , mpiSelf
    use :: Memory_Management, only : allocateArray  , deallocateArray
    use :: Numerical_Ranges , only : Make_Range     , rangeTypeLinear      , rangeTypeLogarithmic
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)               :: self
    double precision                                 , intent(in   )               :: time                             , varianceProgenitor
    type            (treeNode                       ), intent(inout)               :: node
    double precision                                 , parameter                   :: varianceMinimumDefault    =1.0d-2
    double precision                                 , parameter                   :: varianceTolerance         =1.0d-6
    double precision                                 , parameter                   :: massLarge                 =1.0d16
    real            (kind=kind_quad                 ), allocatable  , dimension(:) :: firstCrossingTableRateQuad       , varianceTableRateBaseQuad, &
         &                                                                            varianceTableRateQuad            , barrierTableRateQuad
    double precision                                                               :: barrierRateTest
    class           (excursionSetBarrierClass       ), pointer                     :: excursionSetBarrier_
    class           (cosmologicalMassVarianceClass  ), pointer                     :: cosmologicalMassVariance_
#ifdef USEMPI
    integer                                                                        :: taskCount
#endif
    logical                                                                        :: makeTable
    integer         (c_size_t                       )                              :: loopCount                        , loopCountTotal
    integer                                                                        :: i                                , iTime                    , &
         &                                                                            iVariance                        , j
    double precision                                                               :: timeProgenitor                   , varianceMinimumRate      , &
         &                                                                            massProgenitor
    character       (len=6                          )                              :: label
    type            (varying_string                 )                              :: message
    type            (lockDescriptor                 )                              :: fileLock
    real            (kind=kind_quad                 )                              :: crossingFraction                 , effectiveBarrierInitial  , &
         &                                                                            sigma1f                          , varianceTableStepRate    , &
         &                                                                            barrier                          , growthFactorEffective

    ! Determine if we need to make the table.
    ! Read tables from file if possible.
    if (self%useFile.and..not.self%tableInitializedRate) then
       call self%fileNameInitialize()
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
       call self%fileRead()
       call File_Unlock(fileLock)
    end if
    makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
#ifdef USEMPI
    if (self%coordinatedMPI_) call mpiBarrier()
#endif
    if (makeTable) then
       !$omp critical(farahiRateTabulate)
       ! Attempt to read the file again now that we are within the critical section. If another thread made the file while we were waiting we may be able to skip building the table.
       if (self%useFile) then
          call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
          call self%fileRead()
          call File_Unlock(fileLock)
       end if
       makeTable=.not.self%tableInitializedRate.or.(varianceProgenitor > self%varianceMaximumRate*(1.0d0+varianceTolerance)).or.(time < self%timeMinimumRate).or.(time > self%timeMaximumRate)
       if (makeTable) then
          if (allocated(self%varianceTableRate     )) call deallocateArray(self%varianceTableRate     )
          if (allocated(self%varianceTableRateBase )) call deallocateArray(self%varianceTableRateBase )
          if (allocated(self%timeTableRate         )) call deallocateArray(self%timeTableRate         )
          if (allocated(self%firstCrossingTableRate)) call deallocateArray(self%firstCrossingTableRate)
          if (allocated(self%nonCrossingTableRate  )) call deallocateArray(self%nonCrossingTableRate  )
          if (self%tableInitializedRate) then
             self%timeMinimumRate   =min(self%timeMinimumRate,0.5d0*time)
             self%timeMaximumRate   =max(self%timeMaximumRate,2.0d0*time)
             self%timeTableCountRate=int(log10(self%timeMaximumRate/self%timeMinimumRate)*dble(self%timeNumberPerDecade))+1
          else
             self%timeMinimumRate   =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(farahiRateRedshiftMaximum))
             self%timeMaximumRate   =self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(farahiRateRedshiftMinimum))
             self%timeMinimumRate   =min(self%timeMinimumRate,0.5d0*time)
             self%timeMaximumRate   =max(self%timeMaximumRate,2.0d0*time)
             self%timeTableCountRate=max(int(log10(self%timeMaximumRate/self%timeMinimumRate)*dble(self%timeNumberPerDecade))+1,2)
          end if
          ! Set the default minimum variance.
          varianceMinimumRate       =varianceMinimumDefault
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
               &                          /cosmologicalMassVariance_%rootVariance(massLarge,self%timeMaximumRate*(1.0d0-self%timeStepFractional))
          varianceMinimumRate            =min(                                                                                                                      &
               &                              +varianceMinimumRate                                                                                                , &
               &                              +1.0d-2                                                                                                               &
               &                              *(                                                                                                                    &
               &                                +excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate*(1.0d0-self%timeStepFractional),node,rateCompute=.true.)  &
               &                                *dble(growthFactorEffective)                                                                                        &
               &                                -excursionSetBarrier_%barrier(+0.0d0,self%timeMaximumRate                                ,node,rateCompute=.true.)  &
               &                               )**2                                                                                                                 &
               &                             )
          !![
          <objectDestructor name="excursionSetBarrier_"     />
          <objectDestructor name="cosmologicalMassVariance_"/>
          !!]
          self%varianceMaximumRate       =max(self%varianceMaximumRate,varianceProgenitor)
          self%varianceTableCountRate    =int(log10(self%varianceMaximumRate/varianceMinimumRate)*dble(self%varianceNumberPerDecade))+1
          self%varianceTableCountRateBase=int(self%varianceMaximumRate*dble(self%varianceNumberPerUnit))
          call allocateArray(self%varianceTableRate     ,[1+self%varianceTableCountRate                                                          ],lowerBounds=[0    ])
          call allocateArray(self%varianceTableRateBase ,[                              1+self%varianceTableCountRateBase                        ],lowerBounds=[0    ])
          call allocateArray(self%timeTableRate         ,[                                                                self%timeTableCountRate]                    )
          call allocateArray(self%firstCrossingTableRate,[1+self%varianceTableCountRate,1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[0,0,1])
          call allocateArray(self%nonCrossingTableRate  ,[                              1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[  0,1])
          ! For the variance table, the zeroth point is always zero, higher points are distributed uniformly in variance.
          self%varianceTableRate    (0                                )=0.0d0
          self%varianceTableRate    (1:self%varianceTableCountRate    )=self%varianceRange(varianceMinimumRate,self%varianceMaximumRate,self%varianceTableCountRate      ,exponent =1.0d0          )
          self%varianceTableRateBase(0:self%varianceTableCountRateBase)=Make_Range        (0.0d0              ,self%varianceMaximumRate,self%varianceTableCountRateBase+1,rangeType=rangeTypeLinear)
          ! Allocate temporary arrays used in quad-precision solver for barrier crossing rates.
          allocate(varianceTableRateQuad     (0:self%varianceTableCountRate    ))
          varianceTableRateQuad    =self%varianceTableRate
          allocate(varianceTableRateBaseQuad (0:self%varianceTableCountRateBase))
          varianceTableRateBaseQuad=self%varianceTableRateBase
          ! The time table is logarithmically distributed in time.
          self%timeTableRate=Make_Range(self%timeMinimumRate,self%timeMaximumRate,self%timeTableCountRate,rangeType=rangeTypeLogarithmic)
          ! Loop through the table and solve for the first crossing distribution.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayIndent("solving for excursion set barrier crossing rates",verbosityLevelWorking)
             message="    time: "
             write (label,'(f6.3)') self%timeMinimumRate
             message=message//label//" to "
             write (label,'(f6.3)') self%timeMaximumRate
             message=message//label
             call displayMessage(message,verbosityLevelWorking)
             message="variance: "
             write (label,'(f6.3)') self%varianceMaximumRate
             message=message//label
             call displayMessage(message,verbosityLevelWorking)
#ifdef USEMPI
          end if
#endif
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
             loopCountTotal=int(self%timeTableCountRate,kind=c_size_t)*int(self%varianceTableCountRateBase+1,kind=c_size_t)/int(mpiSelf%count(),kind=c_size_t)+1_c_size_t
          else
#endif
             loopCountTotal=int(self%timeTableCountRate,kind=c_size_t)*int(self%varianceTableCountRateBase+1,kind=c_size_t)
#ifdef USEMPI
          end if
#endif
          loopCount=0
#ifdef USEMPI
          if (self%coordinatedMPI_) then
             self%firstCrossingTableRate=0.0d0
             self%nonCrossingTableRate  =0.0d0
          end if
          taskCount=-1
#endif
          ! Make a call to the barrier function at maximum variance for the minimum and maximum times so that the barrier function
          ! is initialized and covers the whole range we are intereseted in.
          barrierRateTest=self%excursionSetBarrier_%barrier(self%varianceMaximumRate,self%timeMinimumRate*(1.0d0-self%timeStepFractional),node,rateCompute=.true.)
          barrierRateTest=self%excursionSetBarrier_%barrier(self%varianceMaximumRate,self%timeMaximumRate                                ,node,rateCompute=.true.)
          !$omp parallel private(iTime,timeProgenitor,iVariance,varianceTableStepRate,i,j,sigma1f,crossingFraction,barrier,effectiveBarrierInitial,firstCrossingTableRateQuad,excursionSetBarrier_,cosmologicalMassVariance_,barrierTableRateQuad,massProgenitor,growthFactorEffective) if (.not.mpiSelf%isActive() .or. .not.self%coordinatedMPI_)
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
          call allocateArray(barrierTableRateQuad,[self%varianceTableCountRate])
          !$omp do schedule(dynamic)
          do iTime=1,self%timeTableCountRate
             if (.not.allocated(firstCrossingTableRateQuad)) allocate(firstCrossingTableRateQuad(0:self%varianceTableCountRate))
             ! Compute a suitable progenitor time.
             timeProgenitor=self%timeTableRate(iTime)*(1.0d0-self%timeStepFractional)
             ! Loop through the starting variances.
             do iVariance=0,self%varianceTableCountRateBase
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
                ! Construct the barrier table.
                do i=1,self%varianceTableCountRate
                   massProgenitor         =+cosmologicalMassVariance_%mass        (real(+varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance),kind=8),self%timeTableRate(iTime)                        )
                   growthFactorEffective  =+cosmologicalMassVariance_%rootVariance(      massProgenitor                                                       ,self%timeTableRate(iTime)                        ) &
                        &                  /cosmologicalMassVariance_%rootVariance(      massProgenitor                                                       ,     timeProgenitor                              )
                   barrierTableRateQuad(i)=+excursionSetBarrier_     %barrier     (real(+varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance),kind=8),     timeProgenitor      ,node,rateCompute=.true.) &
                        &                  *growthFactorEffective
                end do
                !$omp atomic
                loopCount=loopCount+1_c_size_t
                ! For zero variance, the rate is initialized to zero.
                firstCrossingTableRateQuad(0)=0.0_kind_quad
                ! Compute the step in variance across this first grid cell.
                varianceTableStepRate=varianceTableRateQuad(1)-varianceTableRateQuad(0)
                ! Compute the barrier for the descendent.
                barrier=real(excursionSetBarrier_%barrier(real(varianceTableRateBaseQuad(iVariance),kind=8),self%timeTableRate(iTime),node,rateCompute=.true.),kind=kind_quad)
                ! Compute the first crossing distribution at the first grid cell.
                if (varianceTableRateQuad(1)+varianceTableRateBaseQuad(iVariance) > self%varianceMaximumRate) then
                   firstCrossingTableRateQuad(1)= 0.0_kind_quad
                else
                   firstCrossingTableRateQuad(1)=+2.0_kind_quad                                                  &
                        &                        *(                                                              &
                        &                          +1.0_kind_quad                                                &
                        &                          -erfApproximate(                                              &
                        &                                          +(                                            &
                        &                                            +barrierTableRateQuad(1)                    &
                        &                                            -barrier                                    &
                        &                                           )                                            &
                        &                                          /sqrt(2.0_kind_quad*varianceTableRateQuad(1)) &
                        &                                         )                                              &
                        &                         )                                                              &
                        &                        /varianceTableStepRate
                end if
                do i=2,self%varianceTableCountRate
                   if (varianceTableRateQuad(i)+varianceTableRateBaseQuad(iVariance) > self%varianceMaximumRate) then
                      firstCrossingTableRateQuad(i)=0.0_kind_quad
                   else
                      effectiveBarrierInitial=+barrierTableRateQuad(i)-barrier
                      sigma1f                =+0.0_kind_quad
                      do j=1,i-1
                         varianceTableStepRate=(varianceTableRateQuad(j+1)-varianceTableRateQuad(j-1))/2.0_kind_quad
                         sigma1f=+sigma1f                                                                                   &
                              &  +firstCrossingTableRateQuad(j)                                                             &
                              &  *varianceTableStepRate                                                                     &
                              &  *(                                                                                         &
                              &    +1.0_kind_quad                                                                           &
                              &    -erfApproximate(                                                                         &
                              &                    +(                                                                       &
                              &                      +effectiveBarrierInitial                                               &
                              &                      -barrierTableRateQuad(j)                                               &
                              &                      +barrier                                                               &
                              &                     )                                                                       &
                              &                    /sqrt(2.0_kind_quad*(varianceTableRateQuad(i)-varianceTableRateQuad(j))) &
                              &                   )                                                                         &
                              &   )
                      end do
                      varianceTableStepRate=varianceTableRateQuad(i)-varianceTableRateQuad(i-1)
                      firstCrossingTableRateQuad(i)=max(                                                                  &
                           &                            +0.0_kind_quad,                                                   &
                           &                            +(                                                                &
                           &                              +2.0_kind_quad                                                  &
                           &                              *(                                                              &
                           &                                +1.0_kind_quad                                                &
                           &                                -erfApproximate(                                              &
                           &                                                +effectiveBarrierInitial                      &
                           &                                                /sqrt(2.0_kind_quad*varianceTableRateQuad(i)) &
                           &                                               )                                              &
                           &                               )                                                              &
                           &                              -2.0_kind_quad*sigma1f                                          &
                           &                             )                                                                &
                           &                            /varianceTableStepRate                                            &
                           &                           )
                   end if
                end do
                ! Compute the fraction of trajectories which never cross the barrier.
                crossingFraction=0.0_kind_quad
                do j=0,self%varianceTableCountRate-1
                   varianceTableStepRate=varianceTableRateQuad(j+1)-varianceTableRateQuad(j)
                   crossingFraction=+crossingFraction                  &
                        &           +0.5_kind_quad                     &
                        &           *(                                 &
                        &              firstCrossingTableRateQuad(j  ) &
                        &             +firstCrossingTableRateQuad(j+1) &
                        &            )                                 &
                        &           *varianceTableStepRate
                end do
                ! Compute the rate for trajectories which never cross the barrier.
                self%nonCrossingTableRate(iVariance,iTime)=real(                                   &
                     &                                          +(1.0_kind_quad-crossingFraction)  &
                     &                                          /self%timeTableRate(iTime)         &
                     &                                          /self%timeStepFractional         , &
                     &                                          kind=kind_dble                     &
                     &                                         )
                ! Store the compute crossing rate in our table.
                self%firstCrossingTableRate(:,iVariance,iTime)=real(firstCrossingTableRateQuad,kind=kind_dble)

             end do
             ! Divide through by the time step to get the rate of barrier crossing.
             self%firstCrossingTableRate(:,:,iTime)=+self%firstCrossingTableRate(:,:,iTime) &
                  &                                 /self%timeTableRate         (    iTime) &
                  &                                 /self%timeStepFractional
          end do
          !$omp end do
          !![
          <objectDestructor name="excursionSetBarrier_"     />
          <objectDestructor name="cosmologicalMassVariance_"/>
          !!]
          call deallocateArray(barrierTableRateQuad)
          !$omp end parallel
          ! Deallocate work arrays.
          deallocate(varianceTableRateBaseQuad )
          deallocate(varianceTableRateQuad     )
          if (allocated(firstCrossingTableRateQuad)) deallocate(firstCrossingTableRateQuad)
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             call displayCounterClear(       verbosityLevelWorking)
             call displayUnindent     ("done",verbosityLevelWorking)
#ifdef USEMPI
          end if
          if (self%coordinatedMPI_) then
             call mpiBarrier()
             self%firstCrossingTableRate=mpiSelf%sum(self%firstCrossingTableRate)
             self%  nonCrossingTableRate=mpiSelf%sum(self%  nonCrossingTableRate)
          end if
#endif
          ! Build the interpolators.
          if (allocated(self%interpolatorVarianceRate    )) deallocate(self%interpolatorVarianceRate)
          if (allocated(self%interpolatorVarianceRateBase)) deallocate(self%interpolatorVarianceRateBase)
          if (allocated(self%interpolatorTimeRate        )) deallocate(self%interpolatorTimeRate        )
          allocate(self%interpolatorVarianceRate    )
          allocate(self%interpolatorVarianceRateBase)
          allocate(self%interpolatorTimeRate        )
          self%interpolatorVarianceRate    =interpolator(self%varianceTableRate    )
          self%interpolatorVarianceRateBase=interpolator(self%varianceTableRateBase)
          self%interpolatorTimeRate        =interpolator(self%timeTableRate        )
          ! Set previous variance and time to unphysical values to force recompute of interpolation factors on next call.
          self%varianceRatePrevious=-1.0d0
          self%timeRatePrevious    =-1.0d0
          ! Record that the table is now built.
          self%tableInitializedRate=.true.
          ! Write the table to file if possible.
#ifdef USEMPI
          if (mpiSelf%isMaster() .or. .not.self%coordinatedMPI_) then
#endif
             if (self%useFile) then
                call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
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
    use :: Display           , only : displayIndent, displayMessage  , displayUnindent, verbosityLevelWorking
    use :: File_Utilities    , only : File_Exists  , File_Name_Expand
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : operator(//) , var_str         , varying_string
    use :: Memory_Management , only : allocateArray, deallocateArray
    implicit none
    class           (excursionSetFirstCrossingFarahi), intent(inout)                   :: self
    type            (hdf5Object                     )                                  :: dataFile                   , dataGroup
    double precision                                 , allocatable  , dimension(:    ) :: varianceTableBaseTmp       , varianceTableTmp
    double precision                                 , allocatable  , dimension(:,:  ) :: firstCrossingProbabilityTmp, nonCrossingTableRate
    double precision                                 , allocatable  , dimension(:,:,:) :: firstCrossingRateTmp
    type            (varying_string                 )                                  :: message
    character       (len=32                         )                                  :: label

    ! Initialize file name.
    call self%fileNameInitialize()
    ! Return immediately if the file does not exist.
    if (.not.File_Exists(char(self%fileName))) return
    ! Open the data file.
    !$ call hdf5Access%set()
    call dataFile%openFile(char(self%fileName))
    ! Check if the standard table is populated.
    if (dataFile%hasGroup('probability')) then
       ! Deallocate arrays if necessary.
       if (allocated(self%varianceTable                )) call deallocateArray(self%varianceTable                )
       if (allocated(self%timeTable                    )) call deallocateArray(self%timeTable                    )
       if (allocated(self%firstCrossingProbabilityTable)) call deallocateArray(self%firstCrossingProbabilityTable)
       ! Read the datasets.
       dataGroup=dataFile%openGroup("probability")
       call dataGroup%readDataset('variance'                ,varianceTableTmp           )
       call dataGroup%readDataset('time'                    ,self%timeTable             )
       call dataGroup%readDataset('firstCrossingProbability',firstCrossingProbabilityTmp)
       call dataGroup%close()
       ! Set table sizes and limits.
       self%varianceTableCount=size(varianceTableTmp)-1
       self%timeTableCount    =size(self%timeTable  )
       ! Transfer to tables.
       call allocateArray(self%varianceTable                ,[1+self%varianceTableCount                    ],lowerBounds=[0  ])
       call allocateArray(self%firstCrossingProbabilityTable,[1+self%varianceTableCount,self%timeTableCount],lowerBounds=[0,1])
       self%varianceTable                (0:self%varianceTableCount  )=varianceTableTmp           (1:self%varianceTableCount+1  )
       self%firstCrossingProbabilityTable(0:self%varianceTableCount,:)=firstCrossingProbabilityTmp(1:self%varianceTableCount+1,:)
       call deallocateArray(varianceTableTmp           )
       call deallocateArray(firstCrossingProbabilityTmp)
       ! Set table limits.
       self%timeMinimum      =self%timeTable    (                      1)
       self%timeMaximum      =self%timeTable    (self%    timeTableCount)
       self%varianceMaximum  =self%varianceTable(self%varianceTableCount)
       self%varianceTableStep=self%varianceTable(1)-self%varianceTable(0)
       self%tableInitialized =.true.
       ! Build the interpolators.
       if (allocated(self%interpolatorVariance)) deallocate(self%interpolatorVariance)
       if (allocated(self%interpolatorTime    )) deallocate(self%interpolatorTime    )
       allocate(self%interpolatorVariance)
       allocate(self%interpolatorTime    )
       self%interpolatorVariance=interpolator(self%varianceTable)
       self%interpolatorTime    =interpolator(self%timeTable    )
       ! Report.
       message=var_str('read excursion set first crossing probability from: ')//char(self%fileName)
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
       call displayUnindent(''     ,verbosityLevelWorking)
    end if
    ! Check if the rate table is populated.
    if (dataFile%hasGroup('rate')) then
       ! Deallocate arrays if necessary.
       if (allocated(self%varianceTableRate     )) call deallocateArray(self%varianceTableRate     )
       if (allocated(self%varianceTableRateBase )) call deallocateArray(self%varianceTableRateBase )
       if (allocated(self%timeTableRate         )) call deallocateArray(self%timeTableRate         )
       if (allocated(self%firstCrossingTableRate)) call deallocateArray(self%firstCrossingTableRate)
       if (allocated(self%nonCrossingTableRate  )) call deallocateArray(self%nonCrossingTableRate  )
       ! Read the datasets.
       dataGroup=dataFile%openGroup("rate")
       call dataGroup%readDataset('variance'         ,varianceTableTmp    )
       call dataGroup%readDataset('varianceBase'     ,varianceTableBaseTmp)
       call dataGroup%readDataset('time'             ,self%timeTableRate  )
       call dataGroup%readDataset('firstCrossingRate',firstCrossingRateTmp)
       call dataGroup%readDataset('nonCrossingRate'  ,nonCrossingTableRate)
       call dataGroup%close()
       ! Set table sizes and limits.
       self%varianceTableCountRate    =size(varianceTableTmp    )-1
       self%varianceTableCountRateBase=size(varianceTableBaseTmp)-1
       self%timeTableCountRate        =size(self%timeTableRate  )
       ! Transfer to tablse.
       call allocateArray(self%varianceTableRate     ,[1+self%varianceTableCountRate                                                          ],lowerBounds=[0    ])
       call allocateArray(self%varianceTableRateBase ,[                              1+self%varianceTableCountRateBase                        ],lowerBounds=[  0  ])
       call allocateArray(self%firstCrossingTableRate,[1+self%varianceTableCountRate,1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[0,0,1])
       call allocateArray(self%nonCrossingTableRate  ,[                              1+self%varianceTableCountRateBase,self%timeTableCountRate],lowerBounds=[  0,1])
       self%varianceTableRate     (0:self%varianceTableCountRate                                    )=varianceTableTmp    (1:self%varianceTableCountRate+1                                      )
       self%varianceTableRateBase (                              0:self%varianceTableCountRateBase  )=varianceTableBaseTmp(                                1:self%varianceTableCountRateBase+1  )
       self%firstCrossingTableRate(0:self%varianceTableCountRate,0:self%varianceTableCountRateBase,:)=firstCrossingRateTmp(1:self%varianceTableCountRate+1,1:self%varianceTableCountRateBase+1,:)
       self%nonCrossingTableRate  (                              0:self%varianceTableCountRateBase,:)=nonCrossingTableRate(                                1:self%varianceTableCountRateBase+1,:)
       call deallocateArray(varianceTableTmp    )
       call deallocateArray(varianceTableBaseTmp)
       ! Set table limits.
       self%varianceMaximumRate =self%varianceTableRate(self%varianceTableCountRate)
       self%timeMinimumRate     =self%timeTableRate    (                          1)
       self%timeMaximumRate     =self%timeTableRate    (    self%timeTableCountRate)
       self%tableInitializedRate=.true.
       ! Build the interpolators.
       if (allocated(self%interpolatorVarianceRate    )) deallocate(self%interpolatorVarianceRate)
       if (allocated(self%interpolatorVarianceRateBase)) deallocate(self%interpolatorVarianceRateBase)
       if (allocated(self%interpolatorTimeRate        )) deallocate(self%interpolatorTimeRate        )
       allocate(self%interpolatorVarianceRate    )
       allocate(self%interpolatorVarianceRateBase)
       allocate(self%interpolatorTimeRate        )
       self%interpolatorVarianceRate    =interpolator(self%varianceTableRate    )
       self%interpolatorVarianceRateBase=interpolator(self%varianceTableRateBase)
       self%interpolatorTimeRate        =interpolator(self%timeTableRate        )
       ! Report.
       message=var_str('read excursion set first crossing rates from: ')//char(self%fileName)
       call displayIndent  (message,verbosityLevelWorking)
       write (label,'(e22.16)') self%timeMinimumRate
       message=var_str('    time minimum: ')//label//' Gyr'
       call displayMessage (message,verbosityLevelWorking)
       write (label,'(e22.16)') self%timeMaximumRate
       message=var_str('    time maximum: ')//label//' Gyr'
       write (label,'(e22.16)') self%varianceMaximumRate
       message=var_str('variance minimum: ')//label
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
    call dataFile%openFile(char(self%fileName),overWrite=.true.,chunkSize=100_hsize_t,compressionLevel=9)
    ! Check if the standard table is populated.
    if (self%tableInitialized) then
       dataGroup=dataFile%openGroup("probability")
       call dataGroup%writeDataset(self%varianceTable                ,'variance'                ,'The variance at which results are tabulated.'                         )
       call dataGroup%writeDataset(self%timeTable                    ,'time'                    ,'The cosmic times at which results are tabulated.'                     )
       call dataGroup%writeDataset(self%firstCrossingProbabilityTable,'firstCrossingProbability','The probability of first crossing as a function of variance and time.')
       call dataGroup%close()
       ! Report.
       message=var_str('write excursion set first crossing probability to: ')//char(self%fileName)
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
       call displayUnindent(''     ,verbosityLevelWorking)
    end if
    ! Check if the rate table is populated.
    if (self%tableInitializedRate) then
       dataGroup=dataFile%openGroup("rate")
       call dataGroup%writeDataset(self%varianceTableRate     ,'variance'         ,'The variance at which results are tabulated.'                               )
       call dataGroup%writeDataset(self%varianceTableRateBase ,'varianceBase'     ,'The variance of the base halo at which results are tabulated.'              )
       call dataGroup%writeDataset(self%timeTableRate         ,'time'             ,'The cosmic times at which results are tabulated.'                           )
       call dataGroup%writeDataset(self%firstCrossingTableRate,'firstCrossingRate','The probability rate of first crossing as a function of variances and time.')
       call dataGroup%writeDataset(self%nonCrossingTableRate  ,'nonCrossingRate'  ,'The probability rate of non crossing as a function of variance and time.')
       call dataGroup%close()
       ! Report.
       message=var_str('wrote excursion set first crossing rates to: ')//char(self%fileName)
       call displayIndent  (message,verbosityLevelWorking)
       write (label,'(e22.16)') self%timeMinimumRate
       message=var_str('    time minimum: ')//label//' Gyr'
       call displayMessage (message,verbosityLevelWorking)
       write (label,'(e22.16)') self%timeMaximumRate
       message=var_str('    time maximum: ')//label//' Gyr'
       write (label,'(e22.16)') self%varianceMaximumRate
       message=var_str('variance minimum: ')//label
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

