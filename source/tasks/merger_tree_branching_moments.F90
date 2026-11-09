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

   ``B``
     the number of splits anywhere in the tree producing a non-primary progenitor of mass at least :math:`m` (the unevolved
     subhalo mass function summed over all hierarchy depths);
   ``C``
     as ``B``, but counting only splits of the main branch (the chain of primary progenitors from the root).

   These are the statistics ``B`` and ``C`` of the ``branchingVariance`` merger tree operator. For a node of mass :math:`M`
   at :math:`\omega`, let :math:`n_1(M,\omega)` be the mean of the count from the node's subtree and :math:`W(M,\omega)` its
   variance. A split into :math:`M-M^\prime` and :math:`M^\prime` changes the expected count by the jump
   :math:`\Delta=c+n_1(M-M^\prime)+n_1(M^\prime)-n_1(M)` (for ``C``, :math:`\Delta=c+n_1(M-M^\prime)-n_1(M)`, since only
   the primary progenitor continues the main branch), where :math:`c=1` if :math:`M^\prime\ge m` and :math:`0` otherwise.
   The expected count from the subtree is a martingale, so along each branch between splits

   .. math::

      \frac{\mathrm{d}n_1}{\mathrm{d}\omega} = -\int f\Delta\,\mathrm{d}M^\prime,\quad
      \frac{\mathrm{d}W}{\mathrm{d}\omega} = -\int f\left[\Delta^2+W(M-M^\prime)+W(M^\prime)-W(M)\right]\mathrm{d}M^\prime,

   with :math:`\mathrm{d}M/\mathrm{d}\omega=-\phi M`, and the :math:`W(M^\prime)` term absent for ``C``. Both vanish at
   ``[redshiftMaximum]`` and for :math:`M\le M_\mathrm{res}`. The equations are integrated, from there to the root, on a
   grid uniform in :math:`\ln M` (from :math:`M_\mathrm{res}` to the root mass) and in :math:`\omega`, along the
   characteristics :math:`\mathrm{d}M/\mathrm{d}\omega=-\phi M` with a second-order (Heun) step and cubic interpolation in
   :math:`\ln M`. The integrals over :math:`M^\prime\in[M_\mathrm{res},M/2]` use the trapezoid rule in :math:`\ln
   M^\prime` on the grid, plus the points :math:`M/2` and :math:`m`. The task requires :math:`\omega` to be independent of
   mass, as it is for a mass-independent critical overdensity and scale-independent growth.

   Written to the group ``[outputGroupName]`` of the output file are the grids (datasets ``mass``, ``omega`` and
   ``redshift``), the rates :math:`\phi` (``fractionSubresolution``), :math:`R=\int f\mathrm{d}M^\prime` from the
   ``probability`` method of the branching probability class (``rateSplit``) and from the quadrature used here
   (``rateSplitQuadrature``), and, in a subgroup ``statistic&lt;S&gt;/threshold&lt;i&gt;`` for each statistic and threshold,
   :math:`n_1` (``countMean``), :math:`W` (``countVariance``), and the source of the variance,
   :math:`\int f\Delta^2\mathrm{d}M^\prime` (``sourceDeltaSquared``), on the grid, with the values at the root as
   attributes (``countMeanRoot``, ``countVarianceRoot``). Two-dimensional datasets are indexed by :math:`\omega` and then
   mass (in the C order of HDF5).
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
     double precision                                                                 :: massRoot                                  , massResolution        , &
          &                                                                              redshiftRoot                              , redshiftMaximum
     double precision                                     , allocatable, dimension(:) :: massThresholds
     integer                                                                          :: countMassPerDecade                        , countOmega            , &
          &                                                                              outputStrideOmega
     type            (varying_string                     )                            :: statistics                                , outputGroupName
     type            (inputParameters                    )                            :: parameters
     logical                                                                          :: nodeComponentsInitialized       =  .false.
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
         &                                                                                redshiftRoot                   , redshiftMaximum
    integer                                                                            :: countMassPerDecade             , countOmega     , &
         &                                                                                outputStrideOmega
    type            (varying_string                     )                              :: statistics                     , outputGroupName

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
      The count statistics to compute, any of ``B`` and ``C`` (see the description of this task).
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
    self=taskMergerTreeBranchingMoments(massRoot,redshiftRoot,massResolution,redshiftMaximum,statistics,massThresholds,countMassPerDecade,countOmega,outputStrideOmega,outputGroupName,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBranchingProbability_,parametersRoot)
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

  function mergerTreeBranchingMomentsConstructorInternal(massRoot,redshiftRoot,massResolution,redshiftMaximum,statistics,massThresholds,countMassPerDecade,countOmega,outputStrideOmega,outputGroupName,cosmologyFunctions_,criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBranchingProbability_,parameters) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`taskMergerTreeBranchingMoments` task class.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char        , len
    implicit none
    type            (taskMergerTreeBranchingMoments     )                              :: self
    double precision                                     , intent(in   )               :: massRoot                       , massResolution , &
         &                                                                                redshiftRoot                   , redshiftMaximum
    double precision                                     , intent(in   ), dimension(:) :: massThresholds
    integer                                              , intent(in   )               :: countMassPerDecade             , countOmega     , &
         &                                                                                outputStrideOmega
    type            (varying_string                     ), intent(in   )               :: statistics                     , outputGroupName
    class           (cosmologyFunctionsClass            ), intent(in   ), target       :: cosmologyFunctions_
    class           (criticalOverdensityClass           ), intent(in   ), target       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass      ), intent(in   ), target       :: cosmologicalMassVariance_
    class           (mergerTreeBranchingProbabilityClass), intent(in   ), target       :: mergerTreeBranchingProbability_
    type            (inputParameters                    ), intent(in   ), target       :: parameters
    integer                                                                            :: i
    character       (len=:                              ), allocatable                 :: statisticsList
    !![
    <constructorAssign variables="massRoot, redshiftRoot, massResolution, redshiftMaximum, statistics, massThresholds, countMassPerDecade, countOmega, outputStrideOmega, outputGroupName, *cosmologyFunctions_, *criticalOverdensity_, *cosmologicalMassVariance_, *mergerTreeBranchingProbability_"/>
    !!]

    if (len(statistics) < 1) call Error_Report('at least one statistic is required'//{introspection:location})
    statisticsList=char(statistics)
    do i=1,len(statisticsList)
       if (statisticsList(i:i) /= 'B' .and. statisticsList(i:i) /= 'C') &
            & call Error_Report('statistics must be any of "B" and "C"'//{introspection:location})
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
    double precision                                     , parameter                       :: toleranceMass                   =  1.0d-9, toleranceOmega     =1.0d-6
    type            (treeNode                           ), pointer                         :: node
    class           (nodeComponentBasic                 ), pointer                         :: basic
    class           (criticalOverdensityClass           ), pointer    , save               :: criticalOverdensity_            => null()
    class           (cosmologicalMassVarianceClass      ), pointer    , save               :: cosmologicalMassVariance_       => null()
    class           (mergerTreeBranchingProbabilityClass), pointer    , save               :: mergerTreeBranchingProbability_ => null()
    type            (treeNode                           ), pointer    , save               :: nodeThread                      => null()
    type            (inputParameters                    ), allocatable, save               :: parameters
    !$omp threadprivate(criticalOverdensity_,cosmologicalMassVariance_,mergerTreeBranchingProbability_,nodeThread,parameters)
    double precision                                     , allocatable  , dimension(:    ) :: logMass                                  , mass                  , &
         &                                                                                    omega                                    , redshift              , &
         &                                                                                    logMassNode                              , weightNode            , &
         &                                                                                    kernelCurrent                            , fractionPrevious      , &
         &                                                                                    rateCurrent                              , rateQuadratureCurrent , &
         &                                                                                    logMassFoot                              , countMeanFoot         , &
         &                                                                                    countVarianceFoot                        , sourceMeanFoot        , &
         &                                                                                    sourceVarianceFoot                       , countMeanPredicted    , &
         &                                                                                    countVariancePredicted                   , sourceMeanPredicted   , &
         &                                                                                    sourceVariancePredicted                  , sourceSquaredPredicted, &
         &                                                                                    sourceSquaredCurrent                     , fractionCurrent
    double precision                                     , allocatable  , dimension(:,:  ) :: weightsSecondary                         , weightsPrimary        , &
         &                                                                                    countMean                                , countVariance         , &
         &                                                                                    sourceMean                               , sourceVariance        , &
         &                                                                                    fractionSubresolution                    , rateSplit             , &
         &                                                                                    rateSplitQuadrature
    double precision                                     , allocatable  , dimension(:,:,:) :: countMeanTable                           , countVarianceTable    , &
         &                                                                                    sourceSquaredTable
    integer                                              , allocatable  , dimension(:    ) :: nodeStart                                , nodeCount             , &
         &                                                                                    baseSecondary                            , basePrimary           , &
         &                                                                                    thresholdNode                            , sideNode              , &
         &                                                                                    caseThreshold
    integer                                              , allocatable  , dimension(:,:  ) :: countAdded
    double precision                                     , allocatable  , dimension(:    ) :: logMassCandidate
    integer                                              , allocatable  , dimension(:    ) :: thresholdCandidate                       , sideCandidate
    integer                                                                                :: countCandidates
    logical                                                                                :: countingOnly
    logical                                              , allocatable  , dimension(:    ) :: terminatedPrevious                       , terminatedCurrent
    character       (len=1                              ), allocatable  , dimension(:    ) :: caseStatistic
    double precision                                                                       :: timeNow                                  , timeRoot              , &
         &                                                                                    timeEarliest                             , omegaRoot             , &
         &                                                                                    omegaEnd                                 , stepLogMass           , &
         &                                                                                    stepOmega                                , logMassHalf           , &
         &                                                                                    logMassFootTrial                         , rateStepMaximum       , &
         &                                                                                    quadratureErrorMaximum
    integer                                                                                :: countMass                                , countCases            , &
         &                                                                                    countNodes                               , countThresholds       , &
         &                                                                                    i                                        , j                     , &
         &                                                                                    k                                        , iCase                 , &
         &                                                                                    iStatistic                               , iThreshold
    type            (hdf5Group                          )                                  :: outputGroup                              , caseGroup             , &
         &                                                                                    statisticGroup
    type            (varying_string                     )                                  :: message                                  , groupName
    character       (len=16                             )                                  :: label

    call displayIndent('Begin task: merger tree branching moments')
    call Node_Components_Thread_Initialize(self%parameters)
    node  => treeNode      (                 )
    basic => node    %basic(autoCreate=.true.)
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
    stepOmega         =(omegaEnd-omegaRoot)/dble(self%countOmega)
    allocate(omega   (0:self%countOmega))
    allocate(redshift(0:self%countOmega))
    do k=0,self%countOmega
       omega   (k)=omegaRoot+dble(k)*stepOmega
       redshift(k)=self%cosmologyFunctions_%redshiftFromExpansionFactor(                                                             &
            &       self%cosmologyFunctions_%expansionFactor            (                                                            &
            &        self%criticalOverdensity_%timeOfCollapse            (criticalOverdensity=omega(k),mass=self%massRoot,node=node) & 
            &                                                           )                                                            &
            &                                                          )
    end do
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
    ! Allocate the solution and its tables.
    allocate(countMean              (0:countMass                                          ,countCases))
    allocate(countVariance          (0:countMass                                          ,countCases))
    allocate(sourceMean             (0:countMass                                          ,countCases))
    allocate(sourceVariance         (0:countMass                                          ,countCases))
    allocate(countMeanTable         (0:countMass ,0:self%countOmega/self%outputStrideOmega,countCases))
    allocate(countVarianceTable     (0:countMass ,0:self%countOmega/self%outputStrideOmega,countCases))
    allocate(sourceSquaredTable     (0:countMass ,0:self%countOmega/self%outputStrideOmega,countCases))
    allocate(fractionSubresolution  (0:countMass ,0:self%countOmega/self%outputStrideOmega           ))
    allocate(rateSplit              (0:countMass ,0:self%countOmega/self%outputStrideOmega           ))
    allocate(rateSplitQuadrature    (0:countMass ,0:self%countOmega/self%outputStrideOmega           ))
    allocate(kernelCurrent          (  countNodes                                                    ))
    allocate(fractionPrevious       (0:countMass                                                     ))
    allocate(fractionCurrent        (0:countMass                                                     ))
    allocate(rateCurrent            (0:countMass                                                     ))
    allocate(rateQuadratureCurrent  (0:countMass                                                     ))
    allocate(terminatedPrevious     (0:countMass                                                     ))
    allocate(terminatedCurrent      (0:countMass                                                     ))
    allocate(logMassFoot            (0:countMass                                                     ))
    allocate(countMeanFoot          (0:countMass                                                     ))
    allocate(countVarianceFoot      (0:countMass                                                     ))
    allocate(sourceMeanFoot         (0:countMass                                                     ))
    allocate(sourceVarianceFoot     (0:countMass                                                     ))
    allocate(countMeanPredicted     (0:countMass                                                     ))
    allocate(countVariancePredicted (0:countMass                                                     ))
    allocate(sourceMeanPredicted    (0:countMass                                                     ))
    allocate(sourceVariancePredicted(0:countMass                                                     ))
    allocate(sourceSquaredPredicted (0:countMass                                                     ))
    allocate(sourceSquaredCurrent   (0:countMass                                                     ))
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
    ! The end condition: no count from branches at the maximum redshift.
    call levelEvaluate(self%countOmega,kernelCurrent,fractionCurrent,rateCurrent,rateQuadratureCurrent,terminatedCurrent)
    !$omp single
    countMean    =0.0d0
    countVariance=0.0d0
    !$omp end single
    do iCase=1,countCases
       call sources(iCase,kernelCurrent,terminatedCurrent,countMean(:,iCase),countVariance(:,iCase),sourceMean(:,iCase),sourceVariance(:,iCase),sourceSquaredTable(:,self%countOmega/self%outputStrideOmega,iCase))
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
       ! Feet of the characteristics: the mass at ω_{k+1} of a branch of mass M_i at ω_k, with the subresolution mass loss
       ! rate averaged over the step.
       !$omp do private(logMassFootTrial)
       do i=0,countMass
          if (terminatedCurrent(i)) then
             logMassFoot(i)=-huge(0.0d0)
          else
             logMassFootTrial=logMass(i)-fractionCurrent(i)*stepOmega
             if (logMassFootTrial > logMass(1)) then
                logMassFoot(i)=logMass(i)-0.5d0*(fractionCurrent(i)+interpolateFraction(fractionPrevious,terminatedPrevious,logMassFootTrial))*stepOmega
             else
                logMassFoot(i)=logMassFootTrial
             end if
          end if
       end do
       !$omp end do
       do iCase=1,countCases
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
          !$omp single
          countMean    (:,iCase)=countMeanFoot    +0.5d0*stepOmega*(sourceMeanFoot    +sourceMeanPredicted    )
          countVariance(:,iCase)=countVarianceFoot+0.5d0*stepOmega*(sourceVarianceFoot+sourceVariancePredicted)
          where (terminatedCurrent)
             countMean    (:,iCase)=0.0d0
             countVariance(:,iCase)=0.0d0
          end where
          !$omp end single
          call sources(iCase,kernelCurrent,terminatedCurrent,countMean(:,iCase),countVariance(:,iCase),sourceMean(:,iCase),sourceVariance(:,iCase),sourceSquaredCurrent)
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
    ! Output.
    !$ call hdf5Access%set()
    outputGroup=outputFile%openGroup(char(self%outputGroupName),'Mean and variance of subhalo counts in merger trees from the moment equations.')
    call outputGroup%writeAttribute(self%massRoot                                                       ,'massRoot'                                                                                            )
    call outputGroup%writeAttribute(self%massResolution                                                 ,'massResolution'                                                                                      )
    call outputGroup%writeAttribute(self%redshiftRoot                                                   ,'redshiftRoot'                                                                                        )
    call outputGroup%writeAttribute(     omegaRoot                                                      ,'omegaRoot'                                                                                           )
    call outputGroup%writeAttribute(     omegaEnd                                                       ,'omegaEnd'                                                                                            )
    call outputGroup%writeDataset  (     mass                                                           ,'mass'                 ,'Masses of the grid.'                                                         )
    call outputGroup%writeDataset  (     omega                (0:self%countOmega:self%outputStrideOmega),'omega'                ,'Values of ω of the grid.'                                                    )
    call outputGroup%writeDataset  (     redshift             (0:self%countOmega:self%outputStrideOmega),'redshift'             ,'Redshifts of the ω grid.'                                                    )
    call outputGroup%writeDataset  (     fractionSubresolution                                          ,'fractionSubresolution','Fractional rate of mass loss to subresolution halos, per unit ω.'            )
    call outputGroup%writeDataset  (     rateSplit                                                      ,'rateSplit'            ,'Rate of splits, per unit ω, from the branching probability class.'           )
    call outputGroup%writeDataset  (     rateSplitQuadrature                                            ,'rateSplitQuadrature'  ,'Rate of splits, per unit ω, from the quadrature used to solve the equations.')
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
         call stencil(logMassNode(q)                                       ,baseSecondary(q),weightsSecondary(:,q))
         call stencil(logMass(i)+log(1.0d0-exp(logMassNode(q)-logMass(i))),basePrimary  (q),weightsPrimary  (:,q))
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
      double precision, intent(in   ), dimension(0: ) :: countMean_    , countVariance_
      double precision, intent(  out), dimension(0: ) :: sourceMean_   , sourceVariance_ , &
           &                                             sourceSquared_
      double precision                                :: delta         , varianceChildren, &
           &                                             sumMean       , sumVariance     , &
           &                                             sumSquared
      integer                                         :: i             , q

      ! Every thread executes this procedure, so all writes to the (shared) results are made within the shared loop. The sums are
      ! accumulated in local scalars so that threads do not repeatedly write to neighboring array elements.
      !$omp do schedule(dynamic)
      do i=0,countMass
         sumMean    =0.0d0
         sumVariance=0.0d0
         sumSquared =0.0d0
         if (terminated(i)) then
            sourceMean_    (i)=0.0d0
            sourceVariance_(i)=0.0d0
            sourceSquared_ (i)=0.0d0
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
         end do
         sourceMean_    (i)=sumMean
         sourceVariance_(i)=sumVariance
         sourceSquared_ (i)=sumSquared
      end do
      !$omp end do
      return
    end subroutine sources

  end subroutine mergerTreeBranchingMomentsPerform
