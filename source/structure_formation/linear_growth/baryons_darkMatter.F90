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

  !+    Contributions to this file made by: Andrew Benson. The pinning of the linear growth tabulation to absolute lattices for
  !+    issue #1317 was drafted with assistance from Claude, and reviewed and verified by Andrew Benson.

  !!{RST
  An implementation of linear growth of cosmological structure in models containing baryons and dark matter. Assumes no growth of radiation perturbations.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctions      , cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParameters     , cosmologyParametersClass     , hubbleUnitsTime
  use :: File_Utilities            , only : lockDescriptor
  use :: Intergalactic_Medium_State, only : intergalacticMediumState, intergalacticMediumStateClass
  use :: Numerical_Ranges          , only : rangeLattice
  use :: Tables                    , only : table2DLogLogLin

  !![
  <linearGrowth name="linearGrowthBaryonsDarkMatter" docformat="rst">
   <description>
   Linear growth of cosmological density perturbations in models containing both baryons and collisionless dark matter, computed by numerically integrating the coupled growth equations. Radiation perturbation growth is neglected. The integration is initialized at the redshift ``[redshiftInitial]`` and can use CAMB to set transfer function wavenumber sampling.
   </description>
   <deepCopy>
    <functionClass variables="linearGrowthCollisionlessMatter_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="linearGrowthCollisionlessMatter_"/>
   </stateStorable>
  </linearGrowth>
  !!]
  type, extends(linearGrowthClass) :: linearGrowthBaryonsDarkMatter
     !!{RST
     A linear growth of cosmological structure contrast class in models containing baryons and dark matter. Assumes no growth of radiation perturbations.
     !!}
     private
     logical                                                                      :: tableInitialized                 =  .false., darkMatterOnlyInitialConditions
     double precision                                                             :: tableTimeMinimum                           , tableTimeMaximum                               , &
          &                                                                          tableWavenumberMinimum                     , tableWavenumberMaximum                         , &
          &                                                                          fractionDarkMatter                         , fractionBaryons                                , &
          &                                                                          normalizationMatterDominated               , redshiftInitial                                , &
          &                                                                          redshiftInitialDelta
     integer                                                                      :: cambCountPerDecade
     type            (table2DLogLogLin               )                            :: growthFactor
     ! Lattices to which the two axes of the tabulation are pinned. These are the source of truth for its extent; the limits
     ! above are derived from them. The time axis begins at the epoch from which the growth equations are integrated, so its
     ! lower bound never moves and the axis is only ever extended upward.
     type            (rangeLattice                   )                            :: latticeTime                                , latticeWavenumber
     ! The state of the growth equations at the final tabulated epoch, one entry per wavenumber, together with the factor by
     ! which each wavenumber was normalized. This is all that is needed to resume the integration when the time axis is
     ! extended, and it costs five numbers per wavenumber rather than the four further tables which retaining the state at
     ! every epoch would require.
     !
     ! The state is held *unnormalized* - as the integration itself produced it - and deliberately not divided by
     ! `normalizationFactor`. The equations are linear, so a normalized solution is also a solution and could be integrated
     ! directly; but their solver applies an absolute tolerance as well as a relative one, and an absolute tolerance is not
     ! invariant under a change of scale. Resuming from normalized values would therefore integrate the Jeans-suppressed baryon
     ! modes - which are smaller than the dark matter modes by some seven orders of magnitude - against an effective tolerance
     ! coarser by the normalization factor, and the result would differ from that of an uninterrupted integration by of order
     ! `10⁻⁴` of their amplitude. Resuming from exactly the numbers an uninterrupted integration would have held at that epoch
     ! instead reproduces it bit for bit, which is what makes the tabulation a function of the range alone rather than of the
     ! sequence of requests which produced it.
     double precision                                 , allocatable, dimension(:) :: valueDarkMatterFinal                      , derivativeDarkMatterFinal                       , &
               &                                                                     valueBaryonsFinal                         , derivativeBaryonsFinal                          , &
               &                                                                     normalizationFactor
     type            (varying_string                 )                            :: fileName
     class           (cosmologyParametersClass       ), pointer                   :: cosmologyParameters_             => null() , cosmologyParametersInitialConditions_ => null()
     class           (cosmologyFunctionsClass        ), pointer                   :: cosmologyFunctions_              => null()
     class           (intergalacticMediumStateClass  ), pointer                   :: intergalacticMediumState_        => null()
     type            (linearGrowthCollisionlessMatter), pointer                   :: linearGrowthCollisionlessMatter_ => null()
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulate linear growth factor."              method="retabulate" />
       <method description="Read the tabulated mass variance from file." method="fileWrite"  />
       <method description="Read the tabulated mass variance from file." method="fileRead"   />
       <method description="Return true if the table must be remade."    method="remakeTable"/>
     </methods>
     !!]
     final     ::                                         baryonsDarkMatterDestructor
     procedure :: value                                => baryonsDarkMatterValue
     procedure :: logarithmicDerivativeExpansionFactor => baryonsDarkMatterLogarithmicDerivativeExpansionFactor
     procedure :: logarithmicDerivativeWavenumber      => baryonsDarkMatterLogarithmicDerivativeWavenumber
     procedure :: retabulate                           => baryonsDarkMatterRetabulate
     procedure :: isWavenumberDependent                => baryonsDarkMatterIsWavenumberDependent
     procedure :: remakeTable                          => baryonsDarkMatterRemakeTable
     procedure :: fileWrite                            => baryonsDarkMatterFileWrite
     procedure :: fileRead                             => baryonsDarkMatterFileRead
  end type linearGrowthBaryonsDarkMatter

  interface linearGrowthBaryonsDarkMatter
     !!{RST
     Constructors for the :galacticus-class:`linearGrowthBaryonsDarkMatter` linear growth class.
     !!}
     module procedure baryonsDarkMatterConstructorParameters
     module procedure baryonsDarkMatterConstructorInternal
  end interface linearGrowthBaryonsDarkMatter

  ! Tolerance parameter used to ensure times do not exceed that at the Big Crunch.
  double precision                               , parameter               :: timeToleranceRelative               =1.0d-4

  ! Reference wavenumber used when no wavenumber is specified. Small enough that it should be into the regime where baryon
  ! suppression is negligible, but large enough to avoid the BAO region.
  double precision                               , parameter               :: wavenumberReference                 =1.0d+0

  ! Indices of tables for baryons and dark matter.
  integer                                        , parameter               :: indexDarkMatter                     =1
  integer                                        , parameter               :: indexBaryons                        =2

  ! Seed ranges for the tabulation - the ranges spanned when nothing requires a wider one. Because they are fixed, and derived
  ! from no request, every tabulation contains them, so any two tabulations overlap and can always be merged. The wavenumber
  ! seed already spans whole decades.
  double precision                               , parameter               :: timeMaximumSeed                     =2.0d+1
  double precision                               , parameter               :: wavenumberMinimumSeed               =1.0d-3                             , wavenumberMaximumSeed=1.0d+4

  ! Densities of tabulation points, and the intervals - in lattice steps - to which the bounds of the two axes are pinned. The
  ! time axis is pinned to the lattice points themselves: it carries a thousand points per decade, so anchoring to whole
  ! decades could add a thousand epochs, each of which costs an integration step at every wavenumber. The wavenumber axis is
  ! pinned to whole decades, which is the granularity of its seed range and of the factor-of-two margin applied to a request.
  integer                                        , parameter               :: growthTablePointsPerDecadeTime      =1000
  integer                                        , parameter               :: growthTablePointsPerDecadeWavenumber= 100
  integer                                        , parameter               :: anchorEveryTime                     =   1
  integer                                        , parameter               :: anchorEveryWavenumber               =growthTablePointsPerDecadeWavenumber

  ! Lock used for file access.
  type            (lockDescriptor               )                          :: fileLock

  ! Variables used in ODE solving to allow for parallelism.
  double precision                                                         :: wavenumber_
  class           (cosmologyFunctionsClass      ), pointer                 :: cosmologyFunctions_
  class           (intergalacticMediumStateClass), pointer                 :: intergalacticMediumState_
  !$omp threadprivate(wavenumber_, cosmologyFunctions_,intergalacticMediumState_)

contains

  function baryonsDarkMatterConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`linearGrowthBaryonsDarkMatter` linear growth class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (linearGrowthBaryonsDarkMatter)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_           , cosmologyParametersInitialConditions_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (intergalacticMediumStateClass), pointer       :: intergalacticMediumState_
    double precision                                               :: redshiftInitial                , redshiftInitialDelta
    integer                                                        :: cambCountPerDecade
    logical                                                        :: darkMatterOnlyInitialConditions

    !![
    <inputParameter docformat="rst">
      <name>redshiftInitial</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <description>
      The initial redshift from which integration of linear growth should be begin.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>redshiftInitialDelta</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The initial step in redshift used to estimate growth rates of perturbations using finite differencing.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>cambCountPerDecade</name>
      <source>parameters</source>
      <defaultValue>0</defaultValue>
      <description>
      The number of points per decade of wavenumber to compute in the CAMB transfer function. A value of 0 allows CAMB to choose what it considers to be optimal spacing of wavenumbers.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>darkMatterOnlyInitialConditions</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, set the initial conditions for baryonic modes using the dark matter mode initial conditions.
      </description>
    </inputParameter>
    <objectBuilder    class="cosmologyParameters"      name="cosmologyParameters_"                  source="parameters"                                                     />
    <objectBuilder    class="cosmologyFunctions"       name="cosmologyFunctions_"                   source="parameters"                                                     />
    <objectBuilder    class="intergalacticMediumState" name="intergalacticMediumState_"             source="parameters"                                                     />
    !!]
    if (parameters%isPresent('cosmologyParametersInitialConditions')) then
       !![
       <objectBuilder class="cosmologyParameters"      name="cosmologyParametersInitialConditions_" source="parameters" parameterName="cosmologyParametersInitialConditions"/>
       !!]
    else
       !![
       <objectBuilder class="cosmologyParameters"      name="cosmologyParametersInitialConditions_" source="parameters"                                                     />
       !!]
    end if
    self=baryonsDarkMatterConstructorInternal(redshiftInitial,redshiftInitialDelta,cambCountPerDecade,darkMatterOnlyInitialConditions,cosmologyParameters_,cosmologyParametersInitialConditions_,cosmologyFunctions_,intergalacticMediumState_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"                 />
    <objectDestructor name="cosmologyParametersInitialConditions_"/>
    <objectDestructor name="cosmologyFunctions_"                  />
    <objectDestructor name="intergalacticMediumState_"            />
    !!]
    return
  end function baryonsDarkMatterConstructorParameters

  function baryonsDarkMatterConstructorInternal(redshiftInitial,redshiftInitialDelta,cambCountPerDecade,darkMatterOnlyInitialConditions,cosmologyParameters_,cosmologyParametersInitialConditions_,cosmologyFunctions_,intergalacticMediumState_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`linearGrowthBaryonsDarkMatter` linear growth class.
    !!}
    use :: File_Utilities, only : Directory_Make, File_Path
    use :: Error         , only : Error_Report
    use :: Input_Paths   , only : inputPath     , pathTypeDataDynamic
    implicit none
    type            (linearGrowthBaryonsDarkMatter)                           :: self
    double precision                                          , intent(in   ) :: redshiftInitial                , redshiftInitialDelta
    integer                                                   , intent(in   ) :: cambCountPerDecade
    logical                                                   , intent(in   ) :: darkMatterOnlyInitialConditions
    class           (cosmologyParametersClass     ), target   , intent(in   ) :: cosmologyParameters_           , cosmologyParametersInitialConditions_
    class           (cosmologyFunctionsClass      ), target   , intent(in   ) :: cosmologyFunctions_
    class           (intergalacticMediumStateClass), target   , intent(in   ) :: intergalacticMediumState_
    double precision                                                          :: timeNow
    !![
    <constructorAssign variables="redshiftInitial, redshiftInitialDelta, cambCountPerDecade, darkMatterOnlyInitialConditions, *cosmologyParameters_, *cosmologyParametersInitialConditions_, *cosmologyFunctions_, *intergalacticMediumState_"/>
    !!]

    ! The extent of the tabulation is not set here: it is determined entirely by the lattices to which the two axes are pinned,
    ! which `baryonsDarkMatterRetabulate` builds from the seed ranges above and from whatever is requested. Nothing reads the
    ! limits before then, since `baryonsDarkMatterRemakeTable` consults them only once the table is initialized.
    self%tableInitialized      =.false.
    ! Validate initial redshifts.
    if (redshiftInitialDelta > redshiftInitial) call Error_Report('[redshiftInitialDelta] ≤ [redshiftInitial] is required'//{introspection:location})
    ! Compute dark matter and baryon fractions.
    self%fractionDarkMatter=(+self%cosmologyParameters_%OmegaMatter()-self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter()
    self%fractionBaryons   =(                                        +self%cosmologyParameters_%OmegaBaryon())/self%cosmologyParameters_%OmegaMatter()
    ! Build a linear growth object of the "collisionlessMatter" class which we will use to derive the matter-dominated phase normalization
    ! factor. This is used to figure out the amplitude of a growing mode which grows as the expansion factor at early times, as is
    ! needed for calculations of critical overdensity. The initial conditions we use from CAMB are not pure growing modes, so we
    ! can't compute the normalization from them directly.
    allocate(self%linearGrowthCollisionlessMatter_)
    !![
    <referenceConstruct isResult="yes" owner="self" object="linearGrowthCollisionlessMatter_" constructor="linearGrowthCollisionlessMatter(cosmologyParameters_,cosmologyFunctions_)"/>
    !!]
    timeNow=self%cosmologyFunctions_%cosmicTime(1.0d0)
    self%normalizationMatterDominated=+self%linearGrowthCollisionlessMatter_%value(timeNow,normalize=normalizeMatterDominated) &
         &                            /self%linearGrowthCollisionlessMatter_%value(timeNow                                   )
    self%fileName                    =inputPath(pathTypeDataDynamic)                                                       // &
         &                            'largeScaleStructure/'                                                               // &
         &                            self%objectType      (                                                              )// &
         &                            '_'                                                                                  // &
         &                            self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &                            '.hdf5'
    call Directory_Make(File_Path(self%fileName))
    return
  end function baryonsDarkMatterConstructorInternal

  subroutine baryonsDarkMatterDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`linearGrowthBaryonsDarkMatter` linear growth class.
    !!}
    implicit none
    type (linearGrowthBaryonsDarkMatter), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"                 />
    <objectDestructor name="self%cosmologyParametersInitialConditions_"/>
    <objectDestructor name="self%cosmologyFunctions_"                  />
    <objectDestructor name="self%intergalacticMediumState_"            />
    <objectDestructor name="self%linearGrowthCollisionlessMatter_"     />
    !!]
    call self%growthFactor%destroy()
    return
  end subroutine baryonsDarkMatterDestructor

  subroutine baryonsDarkMatterRetabulate(self,time,wavenumber)
    !!{RST
    Returns the linear growth factor :math:`D(a)` for expansion factor ``aExpansion``, normalized such that :math:`D(1)=1` for a baryons plus dark matter plus cosmological constant cosmology.
    !!}
    use    :: File_Utilities       , only : File_Lock                       , File_Unlock
    use    :: Error                , only : Error_Report
    use    :: Interface_GSL        , only : GSL_Success
    use    :: Interfaces_CAMB      , only : Interface_CAMB_Transfer_Function
    use    :: Numerical_ODE_Solvers, only : odeSolver
    use    :: Numerical_Ranges     , only : Range_Pinned                    , Range_Lattice_Offset, gridSchemePerDecade
    !$ use :: OMP_Lib              , only : omp_lock_kind
    use    :: Table_Labels         , only : extrapolationTypeAbort          , extrapolationTypeFix
    use    :: Tables               , only : table1DGeneric
    implicit none
    class           (linearGrowthBaryonsDarkMatter), intent(inout)              :: self
    double precision                               , intent(in   )              :: time
    double precision                               , intent(in   ), optional    :: wavenumber
    double precision                               , parameter                  :: odeToleranceAbsolute          =1.0d-10, odeToleranceRelative            =1.0d-10
    double precision                               , dimension(4)               :: growthFactorODEVariables
    double precision                               , dimension(2)               :: redshiftsInitial                      , timesInitial
    double precision                               , dimension(4)               :: timeTarget
    double precision                               , dimension(2)               :: wavenumberSeed
    type            (rangeLattice                 )                             :: latticeTime                           , latticeWavenumber
    logical                                        , dimension(:) , allocatable :: wavenumberIsCarried
    logical                                        , dimension(:,:), allocatable :: isComputed
    double precision                               , dimension(:) , allocatable :: normalizationFactorPrevious           , valueDarkMatterPrevious                 , &
         &                                                                         derivativeDarkMatterPrevious          , valueBaryonsPrevious                    , &
         &                                                                         derivativeBaryonsPrevious
    type            (odeSolver                    ), save         , allocatable :: solver
    !$omp threadprivate(solver)
    integer                                                                     :: i                                     , j                                       , &
         &                                                                         iStart                                , jPrevious                               , &
         &                                                                         iPresent                              , offsetWavenumber                        , &
         &                                                                         offsetTime                            , countTimePrevious                       , &
         &                                                                         countWavenumberPrevious
    double precision                                                            :: growthFactorDerivativeBaryons         , growthFactorDerivativeDarkMatter        , &
         &                                                                         timeNow                               , wavenumberLogarithmic                   , &
         &                                                                         timePresent                           , timeBigCrunch                           , &
         &                                                                         timeInitialNominal                    , hPresent                                , &
         &                                                                         linearGrowthFactorPresent             , wavenumberTarget
    logical                                                                     :: carryOver
    integer                                                                     :: growthTableNumberPoints
    type            (table1DGeneric               )                             :: transferFunctionDarkMatter            , transferFunctionBaryons
    integer                                                                     :: countWavenumbers
    !$ integer      (omp_lock_kind                )                             :: lockBaryons                           , lockDarkMatter

    ! Check if we need to recompute our table.
    if (self%remakeTable(time)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
       call self%fileRead()
       call File_Unlock(fileLock,sync=.false.)
    end if
    if (self%remakeTable(time)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.false.)
       ! Find the present-day epoch.
       timePresent=self%cosmologyFunctions_%cosmicTime(1.0d0,collapsingPhase=self%cosmologyParameters_%HubbleConstant() < 0.0d0)
       ! Find the epoch at which the growth equations are to be initialized, as implied by the chosen initial redshift. This is
       ! only the *nominal* epoch: the tabulation is pinned to an absolute lattice, and the integration necessarily begins at
       ! the first point of that lattice, so the actual initial epoch is the lattice point at or below this one and the initial
       ! redshift is derived from it below. The two differ by at most one lattice step - a part in `10³` of a decade in time -
       ! and the nominal epoch depends only on the parameters of this object, so the lattice, and hence the epoch at which CAMB
       ! is asked for initial conditions, is the same for every tabulation this object ever builds.
       timeInitialNominal=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(self%redshiftInitial))
       if (time        < timeInitialNominal) call Error_Report('requested epoch is before the chosen initial epoch'//{introspection:location})
       if (timePresent < timeInitialNominal) call Error_Report('present epoch is before the chosen initial epoch'  //{introspection:location})
       timeBigCrunch     =self%cosmologyFunctions_%timeBigCrunch()
       if (timeBigCrunch > 0.0d0 .and. timeInitialNominal > timeBigCrunch) &
            & call Error_Report('Big Crunch occurs before the chosen initial epoch'//{introspection:location})
       ! Pin the time axis. The request itself is the target - never the range already tabulated, which `latticeCurrent` unions
       ! in afterwards - since folding the current range into the target would apply the margin to an already-margined bound
       ! and ratchet the range upward on every retabulation. No further margin is applied here (`marginFactor=1`) because the
       ! factor of two on the requested epoch already is one.
       timeTarget=[timeInitialNominal,timePresent,2.0d0*time,timeMaximumSeed]
       if (timeBigCrunch > 0.0d0) then
          ! A Big Crunch exists - avoid attempting to tabulate times beyond this epoch.
          latticeTime      =Range_Pinned(                                                            &
               &                                        timeTarget                                 , &
               &                                        growthTablePointsPerDecadeTime             , &
               &                                        gridSchemePerDecade                        , &
               &                         marginFactor  =1.0d0                                      , &
               &                         anchorEvery   =anchorEveryTime                            , &
               &                         limitMaximum  =(1.0d0-timeToleranceRelative)*timeBigCrunch, &
               &                         latticeCurrent=self%latticeTime                             &
               &                        )
       else
          latticeTime      =Range_Pinned(                                                            &
               &                                        timeTarget                                 , &
               &                                        growthTablePointsPerDecadeTime             , &
               &                                        gridSchemePerDecade                        , &
               &                         marginFactor  =1.0d0                                      , &
               &                         anchorEvery   =anchorEveryTime                            , &
               &                         latticeCurrent=self%latticeTime                             &
               &                        )
       end if
       ! Pin the wavenumber axis. The seed range is supplied as `rangeCurrent` so that it is never inflated by the margin, and
       ! so that every tabulation contains it; the requested wavenumber carries the factor-of-two margin which this class has
       ! always applied.
       wavenumberSeed      =[wavenumberMinimumSeed,wavenumberMaximumSeed]
       if (present(wavenumber)) then
          wavenumberTarget =wavenumber
       else
          wavenumberTarget =wavenumberReference
       end if
       latticeWavenumber   =Range_Pinned(                                                      &
            &                                           [wavenumberTarget]                   , &
            &                                            growthTablePointsPerDecadeWavenumber, &
            &                                            gridSchemePerDecade                 , &
            &                            marginFactor  =2.0d0                                , &
            &                            anchorEvery   =anchorEveryWavenumber                , &
            &                            rangeCurrent  =wavenumberSeed                       , &
            &                            latticeCurrent=self%latticeWavenumber                 &
            &                           )
       ! The integration begins at the first point of the pinned time axis; derive the initial redshifts from it. The first
       ! epoch is taken from the lattice directly, rather than by converting the redshift back to a time, so that it is exactly
       ! the abscissa at which the initial conditions are stored.
       timesInitial    (1)=latticeTime%minimum()
       redshiftsInitial(1)=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timesInitial(1)))
       redshiftsInitial(2)=redshiftsInitial(1)-self%redshiftInitialDelta
       timesInitial    (2)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftsInitial(2)))
       ! Get the initial conditions from CAMB.
       call transferFunctionDarkMatter%destroy()
       call transferFunctionBaryons   %destroy()
       call Interface_CAMB_Transfer_Function(self%cosmologyParametersInitialConditions_,redshiftsInitial,latticeWavenumber%maximum(),latticeWavenumber%maximum(),countPerDecade=self%cambCountPerDecade,transferFunctionDarkMatter=transferFunctionDarkMatter,transferFunctionBaryons=transferFunctionBaryons)
       ! Record the extent of the tabulation already in hand, so that the columns carried over into the extended one can be
       ! identified, and then extend onto the new lattices. Carrying over is possible only if the time axis still begins where
       ! it did, since a column is resumed from the top of its carried block and must therefore reach back to the initial
       ! epoch; by construction it always does, as the lower bound of the time axis is a function of the parameters alone.
       countTimePrevious      =0
       countWavenumberPrevious=0
       offsetWavenumber       =0
       offsetTime             =0
       carryOver              =       self%latticeTime      %isDefined()        &
            &                  .and.  self%latticeWavenumber%isDefined()        &
            &                  .and.  allocated(self%normalizationFactor      ) &
            &                  .and.  allocated(self%valueDarkMatterFinal     ) &
            &                  .and.  allocated(self%derivativeDarkMatterFinal) &
            &                  .and.  allocated(self%valueBaryonsFinal        ) &
            &                  .and.  allocated(self%derivativeBaryonsFinal   )
       if (carryOver) then
          countTimePrevious      =self%latticeTime      %count
          countWavenumberPrevious=self%latticeWavenumber%count
          offsetTime             =Range_Lattice_Offset(self%latticeTime      ,latticeTime      )
          offsetWavenumber       =Range_Lattice_Offset(self%latticeWavenumber,latticeWavenumber)
          if (offsetTime /= 0) call Error_Report('the time axis no longer begins at the initial epoch'//{introspection:location})
          call Move_Alloc(self%normalizationFactor      ,normalizationFactorPrevious      )
          call Move_Alloc(self%valueDarkMatterFinal     ,valueDarkMatterPrevious          )
          call Move_Alloc(self%derivativeDarkMatterFinal,derivativeDarkMatterPrevious     )
          call Move_Alloc(self%valueBaryonsFinal        ,valueBaryonsPrevious             )
          call Move_Alloc(self%derivativeBaryonsFinal   ,derivativeBaryonsPrevious        )
       end if
       call self%growthFactor%extend(latticeTime,latticeWavenumber,isComputed,tableCount=2,extrapolationTypeX=extrapolationTypeAbort,extrapolationTypeY=extrapolationTypeFix)
       self%latticeTime           =latticeTime
       self%latticeWavenumber     =latticeWavenumber
       growthTableNumberPoints    =latticeTime      %count
       countWavenumbers           =latticeWavenumber%count
       self%tableTimeMinimum      =timesInitial(1)
       self%tableTimeMaximum      =latticeTime      %maximum()
       self%tableWavenumberMinimum=latticeWavenumber%minimum()
       self%tableWavenumberMaximum=latticeWavenumber%maximum()
       if (allocated(self%normalizationFactor      )) deallocate(self%normalizationFactor      )
       if (allocated(self%valueDarkMatterFinal     )) deallocate(self%valueDarkMatterFinal     )
       if (allocated(self%derivativeDarkMatterFinal)) deallocate(self%derivativeDarkMatterFinal)
       if (allocated(self%valueBaryonsFinal        )) deallocate(self%valueBaryonsFinal        )
       if (allocated(self%derivativeBaryonsFinal   )) deallocate(self%derivativeBaryonsFinal   )
       allocate(self%normalizationFactor      (countWavenumbers))
       allocate(self%valueDarkMatterFinal     (countWavenumbers))
       allocate(self%derivativeDarkMatterFinal(countWavenumbers))
       allocate(self%valueBaryonsFinal        (countWavenumbers))
       allocate(self%derivativeBaryonsFinal   (countWavenumbers))
       allocate(wavenumberIsCarried           (countWavenumbers))
       do j=1,countWavenumbers
          jPrevious             =j-offsetWavenumber
          wavenumberIsCarried(j)=carryOver .and. jPrevious >= 1 .and. jPrevious <= countWavenumberPrevious
       end do
       ! Interpolating factors for the present epoch, at which the growth factor is normalized to unity. They are found once,
       ! since the epoch is the same for every wavenumber, and are applied *within* a single column rather than through the
       ! table's two-dimensional interpolation. Confining them to the column matters: a column computed afresh alongside
       ! carried-over columns which are already normalized would otherwise pick up a contribution from a neighbour on a
       ! different normalization, at the level of the rounding error in the interpolating weight - which would make its own
       ! normalization depend on which columns happened to have been carried over.
       if      (log(timePresent) <  self%growthFactor%xv(1                      )) then
          iPresent=1
       else if (log(timePresent) >= self%growthFactor%xv(growthTableNumberPoints)) then
          iPresent=growthTableNumberPoints-1
       else
          iPresent=int((log(timePresent)-self%growthFactor%xv(1))/latticeTime%stepLogarithmic())+1
       end if
       hPresent=(log(timePresent)-self%growthFactor%xv(iPresent))/latticeTime%stepLogarithmic()
       ! Iterate over wavenumber.
       !$ call OMP_Init_Lock(lockBaryons   )
       !$ call OMP_Init_Lock(lockDarkMatter)
       ! Note that `linearGrowthFactorPresent`, `iStart` and `jPrevious` must be private: each is written once per wavenumber,
       ! and the wavenumbers are shared out between the threads. (`linearGrowthFactorPresent` was formerly an array indexed by
       ! wavenumber, and so was safely shared; it is a scalar here because each column is now normalized as it is finished.)
       ! `iPresent` and `hPresent`, by contrast, are evaluated once before this region and only read within it.
       !$omp parallel private(i,j,wavenumberLogarithmic,growthFactorDerivativeDarkMatter,growthFactorDerivativeBaryons,timeNow,growthFactorODEVariables,iStart,jPrevious,linearGrowthFactorPresent)
       allocate(cosmologyFunctions_      ,mold=self%cosmologyFunctions_      )
       allocate(intergalacticMediumState_,mold=self%intergalacticMediumState_)
       allocate(solver                                                       )
       !$omp critical(linearGrowthBaryonsDrkMttrDeepCopy)
       !![
       <deepCopyReset variables="self%cosmologyFunctions_ self%intergalacticMediumState_"/>
       <deepCopy source="self%cosmologyFunctions_"       destination="cosmologyFunctions_"      />
       <deepCopy source="self%intergalacticMediumState_" destination="intergalacticMediumState_"/>
       <deepCopyFinalize variables="cosmologyFunctions_ intergalacticMediumState_"/>
       !!]
       !$omp end critical(linearGrowthBaryonsDrkMttrDeepCopy)
       !$omp barrier
       solver=odeSolver(4_c_size_t,growthFactorODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)
       !$omp do
       do j=1,countWavenumbers
          wavenumber_          =self%growthFactor%y(j)
          wavenumberLogarithmic=log(wavenumber_)
          if (wavenumberIsCarried(j)) then
             ! This wavenumber was tabulated on an earlier pass. Resume its integration from the top of the block carried over,
             ! from exactly the state which that integration left - unnormalized, so that these steps are taken over precisely
             ! the numbers an uninterrupted integration would have held here.
             iStart                     =countTimePrevious
             jPrevious                  =j-offsetWavenumber
             growthFactorODEVariables(1)=valueDarkMatterPrevious     (jPrevious)
             growthFactorODEVariables(2)=derivativeDarkMatterPrevious(jPrevious)
             growthFactorODEVariables(3)=valueBaryonsPrevious        (jPrevious)
             growthFactorODEVariables(4)=derivativeBaryonsPrevious   (jPrevious)
             self%normalizationFactor(j)=normalizationFactorPrevious (jPrevious)
          else
             ! Solve ODE to get corresponding expansion factors. Initialize with solution from CAMB.
             iStart                          =1
             !$ call OMP_Set_Lock  (lockDarkMatter)
             call    self%growthFactor%populate(exp(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1)),1,j,table=indexDarkMatter)
             growthFactorDerivativeDarkMatter=(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=2)-transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1))*exp(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1))/(timesInitial(2)-timesInitial(1))
             !$ call OMP_Unset_Lock(lockDarkMatter)
             if (self%darkMatterOnlyInitialConditions) then
                ! Initial conditions are to be taken as those of pure dark matter, so the baryon perturbation begins equal to
                ! the dark matter perturbation, and grows at the same initial rate.
                !$ call OMP_Set_Lock  (lockDarkMatter)
                call self%growthFactor%populate(exp(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1)),1,j,table=indexBaryons   )
                growthFactorDerivativeBaryons=(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=2)-transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1))*exp(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1))/(timesInitial(2)-timesInitial(1))
                !$ call OMP_Unset_Lock(lockDarkMatter)
             else
                !$ call OMP_Set_Lock  (lockBaryons   )
                call self%growthFactor%populate(exp(transferFunctionBaryons   %interpolate(wavenumberLogarithmic,table=1)),1,j,table=indexBaryons   )
                growthFactorDerivativeBaryons=(transferFunctionBaryons   %interpolate(wavenumberLogarithmic,table=2)-transferFunctionBaryons   %interpolate(wavenumberLogarithmic,table=1))*exp(transferFunctionBaryons   %interpolate(wavenumberLogarithmic,table=1))/(timesInitial(2)-timesInitial(1))
                !$ call OMP_Unset_Lock(lockBaryons   )
             end if
             ! Take the initial state from the table, into which both perturbations have just been written.
             growthFactorODEVariables(1)=self%growthFactor%zv(1,j,indexDarkMatter)
             growthFactorODEVariables(2)=growthFactorDerivativeDarkMatter
             growthFactorODEVariables(3)=self%growthFactor%zv(1,j,indexBaryons   )
             growthFactorODEVariables(4)=growthFactorDerivativeBaryons
          end if
          ! Carry the state of the equations forward in the solution vector itself, rather than reading the values back from the
          ! table between steps: where the integration is being resumed the table holds normalized values, which are not the
          ! state the equations were left in.
          do i=iStart+1,growthTableNumberPoints
             timeNow=self%growthFactor%x(i-1)
             call solver             %solve   (timeNow,self%growthFactor%x(i),growthFactorODEVariables                             )
             call self  %growthFactor%populate(                               growthFactorODEVariables(1),i,j,table=indexDarkMatter)
             call self  %growthFactor%populate(                               growthFactorODEVariables(3),i,j,table=indexBaryons   )
          end do
          ! Retain the state of the equations at the final epoch, so that the integration can be resumed if the time axis is
          ! later extended.
          self%valueDarkMatterFinal     (j)=growthFactorODEVariables(1)
          self%derivativeDarkMatterFinal(j)=growthFactorODEVariables(2)
          self%valueBaryonsFinal        (j)=growthFactorODEVariables(3)
          self%derivativeBaryonsFinal   (j)=growthFactorODEVariables(4)
       end do
       !$omp end do
       !![
       <objectDestructor name="cosmologyFunctions_"      />
       <objectDestructor name="intergalacticMediumState_"/>
       !!]
       deallocate(solver)
       !$omp do
       do j=1,countWavenumbers
          ! Normalize to growth factor of unity at present day. A column carried over from an earlier pass keeps the factor it
          ! was normalized by, and only the epochs newly integrated for it are divided by that factor: the factor is fixed at
          ! the present day, which lies within the block carried over, so it cannot have changed, and leaving the carried epochs
          ! untouched keeps them exactly as they were.
          if (wavenumberIsCarried(j)) then
             iStart                     =countTimePrevious
             linearGrowthFactorPresent  =self%normalizationFactor(j)
          else
             iStart                     =0
             linearGrowthFactorPresent  =+self%growthFactor%zv(iPresent  ,j,indexDarkMatter)*(1.0d0-hPresent) &
                  &                      +self%growthFactor%zv(iPresent+1,j,indexDarkMatter)*       hPresent
             self%normalizationFactor(j)=linearGrowthFactorPresent
          end if
          do i=iStart+1,growthTableNumberPoints
             call self%growthFactor%populate(self%growthFactor%zv(i,j,indexDarkMatter)/linearGrowthFactorPresent,i,j,table=indexDarkMatter)
             call self%growthFactor%populate(self%growthFactor%zv(i,j,indexBaryons   )/linearGrowthFactorPresent,i,j,table=indexBaryons   )
          end do
       end do
       !$omp end do
       !$omp end parallel
       !$ call OMP_Destroy_Lock(lockBaryons   )
       !$ call OMP_Destroy_Lock(lockDarkMatter)
       self%tableInitialized=.true.
       ! Store file.
       call self%fileWrite()
       call File_Unlock(fileLock)
    end if
    return

  contains

    integer function growthFactorODEs(time,values,derivatives)
      !!{RST
      System of differential equations to solve for the growth factor.
      !!}
      use :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial, hydrogenByMassPrimordial
      use :: Numerical_Constants_Atomic      , only : electronMass          , massHeliumAtom          , massHydrogenAtom
      use :: Numerical_Constants_Physical    , only : boltzmannsConstant
      use :: Numerical_Constants_Prefixes    , only : kilo
      implicit none
      double precision              , intent(in   ) :: time
      double precision, dimension(:), intent(in   ) :: values
      double precision, dimension(:), intent(  out) :: derivatives
      double precision                              :: perturbationDarkMatter             , perturbationBaryons             , &
           &                                           perturbationDarkMatterDerivative1st, perturbationBaryonsDerivative1st, &
           &                                           perturbationDarkMatterDerivative2nd, perturbationBaryonsDerivative2nd, &
           &                                           expansionFactor                    , wavenumberJeans                 , &
           &                                           massParticleMean                   , speedSound

      ! Get expansion factor.
      expansionFactor                    =+                                 cosmologyFunctions_      %expansionFactor (time)
      ! Compute the instantaneous Jeans wavenumber. Note that we want the comoving wavenumber here, so we multiply by the expansion factor.
      massParticleMean                   =+(hydrogenByMassPrimordial*(1.0d0+intergalacticMediumState_%electronFraction(time)*electronMass/massHydrogenAtom)                 +heliumByMassPrimordial               ) &
           &                              /(hydrogenByMassPrimordial*(1.0d0+intergalacticMediumState_%electronFraction(time)                              )/massHydrogenAtom+heliumByMassPrimordial/massHeliumAtom)
      speedSound                         =+sqrt(                                                    &
           &                                    +boltzmannsConstant                                 &
           &                                    *intergalacticMediumState_%temperature       (time) &
           &                                    /massParticleMean                                   &
           &                                   )                                                    &
           &                              /kilo
      if (speedSound > 0.0d0) then
         wavenumberJeans                 =+sqrt(                                                    &
              &                                 +1.5d0                                              &
              &                                 *(                                                  &
              &                                   +cosmologyFunctions_%hubbleParameterEpochal(time) &
              &                                   /speedSound                                       &
              &                                  )**2                                               &
              &                                )                                                    &
              &                           *expansionFactor
      else
         wavenumberJeans=huge(0.0d0)
      end if
      ! Compute perturbations and their derivatives.
      perturbationDarkMatter             =+    values(1)
      perturbationDarkMatterDerivative1st=+    values(2)
      perturbationBaryons                =+    values(3)
      perturbationBaryonsDerivative1st   =+    values(4)
      perturbationDarkMatterDerivative2nd=+    1.5d0                                                                      &
           &                              *    cosmologyFunctions_%expansionRate     (                expansionFactor)**2 &
           &                              *    cosmologyFunctions_%omegaMatterEpochal(expansionFactor=expansionFactor)    &
           &                              *    (                                                                          &
           &                                    +self%    fractionDarkMatter                                              &
           &                                    *     perturbationDarkMatter                                              &
           &                                    +self%    fractionBaryons                                                 &
           &                                    *     perturbationBaryons                                                 &
           &                                   )                                                                          &
           &                              -    2.0d0                                                                      &
           &                              *abs(cosmologyFunctions_%expansionRate     (                expansionFactor))   &
           &                              *    perturbationDarkMatterDerivative1st
      perturbationBaryonsDerivative2nd   =+    1.5d0                                                                      &
           &                              *    cosmologyFunctions_%expansionRate     (                expansionFactor)**2 &
           &                              *    cosmologyFunctions_%omegaMatterEpochal(expansionFactor=expansionFactor)    &
           &                              *    (                                                                          &
           &                                    +self%    fractionDarkMatter                                              &
           &                                    *     perturbationDarkMatter                                              &
           &                                    +self%    fractionBaryons                                                 &
           &                                    *     perturbationBaryons                                                 &
           &                                    -(                                                                        &
           &                                      +wavenumber_                                                            &
           &                                      /wavenumberJeans                                                        &
           &                                     )**2                                                                     &
           &                                    *     perturbationBaryons                                                 &
           &                                   )                                                                          &
           &                              -    2.0d0                                                                      &
           &                              *abs(cosmologyFunctions_%expansionRate     (                expansionFactor))   &
           &                              *    perturbationBaryonsDerivative1st
      ! Set the derivatives in the ODE arrays.
      derivatives    (1)=perturbationDarkMatterDerivative1st
      derivatives    (2)=perturbationDarkMatterDerivative2nd
      derivatives    (3)=perturbationBaryonsDerivative1st
      derivatives    (4)=perturbationBaryonsDerivative2nd
      growthFactorODEs  =GSL_Success
    end function growthFactorODEs

  end subroutine baryonsDarkMatterRetabulate

  double precision function baryonsDarkMatterValue(self,time,expansionFactor,collapsing,normalize,component,wavenumber)
    !!{RST
    Return the linear growth factor at the given epoch.
    !!}
    implicit none
    class           (linearGrowthBaryonsDarkMatter), intent(inout)           :: self
    double precision                               , intent(in   ), optional :: time      , expansionFactor
    logical                                        , intent(in   ), optional :: collapsing
    type            (enumerationNormalizeType     ), intent(in   ), optional :: normalize
    type            (enumerationComponentType     ), intent(in   ), optional :: component
    double precision                               , intent(in   ), optional :: wavenumber
    double precision                                                         :: time_     , wavenumber_
    !![
    <optionalArgument name="normalize" defaultsTo="normalizePresentDay" />
    !!]
    !$GLC attributes unused :: component

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_,wavenumber)
    ! Get wavenumber.
    if (present(wavenumber)) then
       wavenumber_=wavenumber
    else
       wavenumber_=wavenumberReference
    end if
    ! Interpolate to get the linear growth factor.
    baryonsDarkMatterValue=self%growthFactor%interpolate(time_,wavenumber_,table=indexDarkMatter)
    ! Normalize.
    select case (normalize_%ID)
    case (normalizeMatterDominated%ID)
       ! Normalize such that the radii of long-wavelength, growing perturbations behave as the scale factor at early times.
       baryonsDarkMatterValue=+baryonsDarkMatterValue            &
            &                 *self%normalizationMatterDominated
    end select
    return
  end function baryonsDarkMatterValue

  double precision function baryonsDarkMatterLogarithmicDerivativeExpansionFactor(self,time,expansionFactor,collapsing,component,wavenumber)
    !!{RST
    Return the logarithmic gradient of linear growth factor with respect to expansion factor at the given epoch.
    !!}
    implicit none
    class           (linearGrowthBaryonsDarkMatter), intent(inout)           :: self
    double precision                               , intent(in   ), optional :: time       , expansionFactor
    logical                                        , intent(in   ), optional :: collapsing
    type            (enumerationComponentType     ), intent(in   ), optional :: component
    double precision                               , intent(in   ), optional :: wavenumber
    double precision                                                         :: time_      , expansionFactor_, &
         &                                                                      wavenumber_
    !$GLC attributes unused :: component

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_,expansionFactorOut=expansionFactor_)
    ! Remake the table if necessary.
    call self%retabulate(time_,wavenumber)
    ! Get wavenumber.
    if (present(wavenumber)) then
       wavenumber_=wavenumber
    else
       wavenumber_=wavenumberReference
    end if
    ! Interpolate to get the expansion factor.
    baryonsDarkMatterLogarithmicDerivativeExpansionFactor=+self%growthFactor       %interpolateGradient(time_           ,wavenumber_,table=indexDarkMatter,dim=1) &
         &                                                /self%growthFactor       %interpolate        (time_           ,wavenumber_,table=indexDarkMatter      ) &
         &                                                /self%cosmologyFunctions_%expansionRate      (expansionFactor_                                        )
    return
  end function baryonsDarkMatterLogarithmicDerivativeExpansionFactor

  double precision function baryonsDarkMatterLogarithmicDerivativeWavenumber(self,time,expansionFactor,collapsing,component,wavenumber)
    !!{RST
    Return the logarithmic gradient of linear growth factor with respect to expansion factor at the given epoch.
    !!}
    implicit none
    class           (linearGrowthBaryonsDarkMatter), intent(inout)           :: self
    double precision                               , intent(in   ), optional :: time       , expansionFactor
    logical                                        , intent(in   ), optional :: collapsing
    type            (enumerationComponentType     ), intent(in   ), optional :: component
    double precision                               , intent(in   ), optional :: wavenumber
    double precision                                                         :: time_      , expansionFactor_, &
         &                                                                      wavenumber_
    !$GLC attributes unused :: component

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_,expansionFactorOut=expansionFactor_)
    ! Remake the table if necessary.
    call self%retabulate(time_,wavenumber)
    ! Get wavenumber.
    if (present(wavenumber)) then
       wavenumber_=wavenumber
    else
       wavenumber_=wavenumberReference
    end if
    ! Interpolate to get the expansion factor.
    baryonsDarkMatterLogarithmicDerivativeWavenumber=+self%growthFactor%interpolateGradient(time_,wavenumber_,table=indexDarkMatter,dim=2) &
         &                                           /self%growthFactor%interpolate        (time_,wavenumber_,table=indexDarkMatter      ) &
         &                                           *                                            wavenumber_
    return
  end function baryonsDarkMatterLogarithmicDerivativeWavenumber

  logical function baryonsDarkMatterIsWavenumberDependent(self,component)
    !!{RST
    Return true indicating that the growth function is wavenumber-dependent.
    !!}
    implicit none
    class(linearGrowthBaryonsDarkMatter), intent(inout)           :: self
    type (enumerationComponentType     ), intent(in   ), optional :: component
    !$GLC attributes unused :: self, component

    baryonsDarkMatterIsWavenumberDependent=.true.
    return
  end function baryonsDarkMatterIsWavenumberDependent

  logical function baryonsDarkMatterRemakeTable(self,time)
    !!{RST
    Determine if the table should be remade.
    !!}
    implicit none
    class           (linearGrowthBaryonsDarkMatter), intent(inout) :: self
    double precision                               , intent(in   ) :: time

    if (self%tableInitialized) then
       baryonsDarkMatterRemakeTable=time > self%tableTimeMaximum
    else
       baryonsDarkMatterRemakeTable=.true.
    end if
    return
  end function baryonsDarkMatterRemakeTable

  subroutine baryonsDarkMatterFileRead(self)
    !!{RST
    Read tabulated data on linear growth factor from file.
    !!}
    use :: Display       , only : displayMessage        , verbosityLevelWorking
    use :: File_Utilities, only : File_Exists
    use :: HDF5_Access   , only : hdf5Access
    use :: IO_HDF5       , only : hdf5File
    use :: Table_Labels  , only : extrapolationTypeAbort, extrapolationTypeFix
    implicit none
    class           (linearGrowthBaryonsDarkMatter), intent(inout)               :: self
    double precision                               , dimension(:,:), allocatable :: growthFactorDarkMatter, growthFactorBaryons
    double precision                               , dimension(:  ), allocatable :: normalizationFactor   , valueDarkMatter    , &
         &                                                                          derivativeDarkMatter  , valueBaryons       , &
         &                                                                          derivativeBaryons
    type            (rangeLattice                 )                              :: latticeTime           , latticeWavenumber
    logical                                        , dimension(:,:), allocatable :: isComputed
    type            (hdf5File                     )                              :: dataFile

    ! Return immediately if the file does not exist.
    if (.not.File_Exists(self%fileName)) return
    call displayMessage('reading D(k,t) data from: '//self%fileName,verbosityLevelWorking)
    !$ call hdf5Access%set()
    ! The file is opened read-only: a read-write open takes an exclusive lock, which aborts any concurrent reader.
    dataFile=hdf5File(self%fileName,overWrite=.false.,readOnly=.true.)
    ! Recover the lattices on which the stored tabulation was built. A file which records none, or which records lattices this
    ! object would not use, is ignored - as is one written before the derivatives needed to resume the integration were stored.
    call baryonsDarkMatterLatticeRead(dataFile,'time'      ,growthTablePointsPerDecadeTime      ,latticeTime      )
    call baryonsDarkMatterLatticeRead(dataFile,'wavenumber',growthTablePointsPerDecadeWavenumber,latticeWavenumber)
    if (latticeTime%isDefined() .and. latticeWavenumber%isDefined()) then
       call dataFile%readDataset('growthFactorDarkMatter',growthFactorDarkMatter)
       call dataFile%readDataset('growthFactorBaryons'   ,growthFactorBaryons   )
       call dataFile%readDataset('normalizationFactor'   ,normalizationFactor   )
       call dataFile%readDataset('valueDarkMatterFinal'  ,valueDarkMatter       )
       call dataFile%readDataset('derivativeDarkMatter'  ,derivativeDarkMatter  )
       call dataFile%readDataset('valueBaryonsFinal'     ,valueBaryons          )
       call dataFile%readDataset('derivativeBaryons'     ,derivativeBaryons     )
    end if
    !$ call hdf5Access%unset()
    if (.not.latticeTime%isDefined() .or. .not.latticeWavenumber%isDefined()) return
    ! Reject a stored tabulation whose datasets do not match the lattices recorded alongside them.
    if     (                                                               &
         &   size(growthFactorDarkMatter,dim=1) /= latticeTime      %count &
         &  .or.                                                           &
         &   size(growthFactorDarkMatter,dim=2) /= latticeWavenumber%count &
         &  .or.                                                           &
         &   size(growthFactorBaryons   ,dim=1) /= latticeTime      %count &
         &  .or.                                                           &
         &   size(growthFactorBaryons   ,dim=2) /= latticeWavenumber%count &
         &  .or.                                                           &
         &   size(normalizationFactor         ) /= latticeWavenumber%count &
         &  .or.                                                           &
         &   size(valueDarkMatter             ) /= latticeWavenumber%count &
         &  .or.                                                           &
         &   size(derivativeDarkMatter        ) /= latticeWavenumber%count &
         &  .or.                                                           &
         &   size(valueBaryons                ) /= latticeWavenumber%count &
         &  .or.                                                           &
         &   size(derivativeBaryons           ) /= latticeWavenumber%count &
         & ) return
    if (self%tableInitialized) call self%growthFactor%destroy()
    ! Build the table on the lattices recorded in the file, rather than by subdividing the ranges they span, so that its
    ! abscissae are bit-identical to those of any other tabulation built on the same lattices.
    call self%growthFactor%extend  (                                           &
         &                                              latticeTime          , &
         &                                              latticeWavenumber    , &
         &                                              isComputed           , &
         &                          tableCount        =2                     , &
         &                          extrapolationTypeX=extrapolationTypeAbort, &
         &                          extrapolationTypeY=extrapolationTypeFix    &
         &                         )
    call self%growthFactor%populate(growthFactorDarkMatter,table=indexDarkMatter)
    call self%growthFactor%populate(growthFactorBaryons   ,table=indexBaryons   )
    self%latticeTime           =latticeTime
    self%latticeWavenumber     =latticeWavenumber
    self%tableTimeMinimum      =latticeTime      %minimum()
    self%tableTimeMaximum      =latticeTime      %maximum()
    self%tableWavenumberMinimum=latticeWavenumber%minimum()
    self%tableWavenumberMaximum=latticeWavenumber%maximum()
    if (allocated(self%normalizationFactor      )) deallocate(self%normalizationFactor      )
    if (allocated(self%valueDarkMatterFinal     )) deallocate(self%valueDarkMatterFinal     )
    if (allocated(self%derivativeDarkMatterFinal)) deallocate(self%derivativeDarkMatterFinal)
    if (allocated(self%valueBaryonsFinal        )) deallocate(self%valueBaryonsFinal        )
    if (allocated(self%derivativeBaryonsFinal   )) deallocate(self%derivativeBaryonsFinal   )
    call Move_Alloc(normalizationFactor ,self%normalizationFactor      )
    call Move_Alloc(valueDarkMatter     ,self%valueDarkMatterFinal     )
    call Move_Alloc(derivativeDarkMatter,self%derivativeDarkMatterFinal)
    call Move_Alloc(valueBaryons        ,self%valueBaryonsFinal        )
    call Move_Alloc(derivativeBaryons   ,self%derivativeBaryonsFinal   )
    deallocate(growthFactorDarkMatter)
    deallocate(growthFactorBaryons   )
    self%tableInitialized=.true.
    return
  end subroutine baryonsDarkMatterFileRead

  subroutine baryonsDarkMatterLatticeWrite(dataFile,axisName,lattice)
    !!{RST
    Record the ``rangeLattice`` on which an axis of the stored tabulation is built, as attributes named for that axis.
    !!}
    use :: IO_HDF5, only : hdf5File
    implicit none
    type     (hdf5File    ), intent(inout) :: dataFile
    character(len=*       ), intent(in   ) :: axisName
    type     (rangeLattice), intent(in   ) :: lattice

    call dataFile%writeAttribute(lattice%scheme%ID   ,axisName//'GridScheme'  )
    call dataFile%writeAttribute(lattice%pointsPer   ,axisName//'PointsPer'   )
    call dataFile%writeAttribute(lattice%indexMinimum,axisName//'IndexMinimum')
    call dataFile%writeAttribute(lattice%count       ,axisName//'Count'       )
    return
  end subroutine baryonsDarkMatterLatticeWrite

  subroutine baryonsDarkMatterLatticeRead(dataFile,axisName,pointsPer,lattice)
    !!{RST
    Restore the ``rangeLattice`` on which an axis of the stored tabulation was built. The lattice is returned undefined---which
    the caller must treat as the tabulation being unusable---unless the file records one which is self-consistent and which
    uses the density of points that this object would use, so that a file written before the lattices were recorded reports an
    undefined lattice rather than being misread.
    !!}
    use :: IO_HDF5         , only : hdf5File
    use :: Numerical_Ranges, only : enumerationGridSchemeType, gridSchemePerDecade
    implicit none
    type     (hdf5File    ), intent(inout) :: dataFile
    character(len=*       ), intent(in   ) :: axisName
    integer                , intent(in   ) :: pointsPer
    type     (rangeLattice), intent(  out) :: lattice
    integer                                :: schemeStored, pointsPerStored, &
         &                                    indexMinimum, count_

    lattice=rangeLattice()
    if     (                                                      &
         &   .not.dataFile%hasAttribute(axisName//'GridScheme'  ) &
         &  .or.                                                  &
         &   .not.dataFile%hasAttribute(axisName//'PointsPer'   ) &
         &  .or.                                                  &
         &   .not.dataFile%hasAttribute(axisName//'IndexMinimum') &
         &  .or.                                                  &
         &   .not.dataFile%hasAttribute(axisName//'Count'       ) &
         & ) return
    call dataFile%readAttribute(axisName//'GridScheme'  ,schemeStored   )
    call dataFile%readAttribute(axisName//'PointsPer'   ,pointsPerStored)
    call dataFile%readAttribute(axisName//'IndexMinimum',indexMinimum   )
    call dataFile%readAttribute(axisName//'Count'       ,count_         )
    ! Comparing the stored scheme against the one expected is stronger than merely checking that it is a valid member of the
    ! enumeration, so no separate validity test is needed.
    if (enumerationGridSchemeType(schemeStored) /= gridSchemePerDecade) return
    if (pointsPerStored                         /= pointsPer          ) return
    lattice=rangeLattice(enumerationGridSchemeType(schemeStored),pointsPerStored,indexMinimum,count_)
    if (.not.lattice%isDefined()) lattice=rangeLattice()
    return
  end subroutine baryonsDarkMatterLatticeRead

  subroutine baryonsDarkMatterFileWrite(self)
    !!{RST
    Write tabulated data on linear growth factor to file.
    !!}
    use :: Display    , only : displayMessage, verbosityLevelWorking
    use :: HDF5       , only : hsize_t
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5File
    implicit none
    class(linearGrowthBaryonsDarkMatter), intent(inout) :: self
    type (hdf5File                     )                :: dataFile

    ! Open the data file.
    call displayMessage('writing D(k,t) data to: '//self%fileName,verbosityLevelWorking)
    !$ call hdf5Access%set()
    dataFile=hdf5File(self%fileName,overWrite=.true.,chunkSize=100_hsize_t,compressionLevel=9)
    call dataFile%writeDataset  (reshape(self%growthFactor             %zs(table=indexDarkMatter),[self%growthFactor%size(dim=1),self%growthFactor%size(dim=2)]),'growthFactorDarkMatter')
    call dataFile%writeDataset  (reshape(self%growthFactor             %zs(table=indexBaryons   ),[self%growthFactor%size(dim=1),self%growthFactor%size(dim=2)]),'growthFactorBaryons'   )
    ! Store the state of the growth equations at the final tabulated epoch - unnormalized, as the integration left it - together
    ! with the factor by which each wavenumber was normalized, so that the integration can be resumed exactly if the time axis
    ! is later extended.
    call dataFile%writeDataset  (        self%normalizationFactor                                                                                               ,'normalizationFactor'   )
    call dataFile%writeDataset  (        self%valueDarkMatterFinal                                                                                              ,'valueDarkMatterFinal'  )
    call dataFile%writeDataset  (        self%derivativeDarkMatterFinal                                                                                         ,'derivativeDarkMatter'  )
    call dataFile%writeDataset  (        self%valueBaryonsFinal                                                                                                 ,'valueBaryonsFinal'     )
    call dataFile%writeDataset  (        self%derivativeBaryonsFinal                                                                                            ,'derivativeBaryons'     )
    ! Record the lattices on which the two axes are built. The bounds formerly stored alongside them are not: each is a function
    ! of the lattices, and is recovered from them when the file is read, so that a restored tabulation cannot come to be
    ! described differently from a freshly built one.
    call baryonsDarkMatterLatticeWrite(dataFile,'time'      ,self%latticeTime      )
    call baryonsDarkMatterLatticeWrite(dataFile,'wavenumber',self%latticeWavenumber)
    !$ call hdf5Access%unset()
    return
  end subroutine baryonsDarkMatterFileWrite
