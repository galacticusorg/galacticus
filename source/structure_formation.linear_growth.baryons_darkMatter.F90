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
  An implementation of linear growth of cosmological structure in models containing baryons and dark matter. Assumes no growth
  of radiation perturbations.
  !!}

  use :: Cosmology_Functions       , only : cosmologyFunctions      , cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParameters     , cosmologyParametersClass     , hubbleUnitsTime
  use :: File_Utilities            , only : lockDescriptor
  use :: Intergalactic_Medium_State, only : intergalacticMediumState, intergalacticMediumStateClass
  use :: Tables                    , only : table2DLogLogLin

  !![
  <linearGrowth name="linearGrowthBaryonsDarkMatter">
   <description>Linear growth of cosmological structure in models containing baryons and dark matter. Assumes no growth of radiation perturbations.</description>
   <deepCopy>
    <functionClass variables="linearGrowthCollisionlessMatter_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="linearGrowthCollisionlessMatter_"/>
   </stateStorable>
  </linearGrowth>
  !!]
  type, extends(linearGrowthClass) :: linearGrowthBaryonsDarkMatter
     !!{
     A linear growth of cosmological structure contrast class in models containing baryons and dark matter. Assumes no growth
     of radiation perturbations.
     !!}
     private
     logical                                                    :: tableInitialized                 =  .false., darkMatterOnlyInitialConditions
     double precision                                           :: tableTimeMinimum                           , tableTimeMaximum                               , &
          &                                                        tableWavenumberMinimum                     , tableWavenumberMaximum                         , &
          &                                                        fractionDarkMatter                         , fractionBaryons                                , &
          &                                                        normalizationMatterDominated               , redshiftInitial                                , &
          &                                                        redshiftInitialDelta
     integer                                                    :: cambCountPerDecade
     type            (table2DLogLogLin               )          :: growthFactor
     type            (varying_string                 )          :: fileName
     class           (cosmologyParametersClass       ), pointer :: cosmologyParameters_             => null() , cosmologyParametersInitialConditions_ => null()
     class           (cosmologyFunctionsClass        ), pointer :: cosmologyFunctions_              => null()
     class           (intergalacticMediumStateClass  ), pointer :: intergalacticMediumState_        => null()
     type            (linearGrowthCollisionlessMatter), pointer :: linearGrowthCollisionlessMatter_ => null()
   contains
     !![
     <methods>
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
     !!{
     Constructors for the \refClass{linearGrowthBaryonsDarkMatter} linear growth class.
     !!}
     module procedure baryonsDarkMatterConstructorParameters
     module procedure baryonsDarkMatterConstructorInternal
  end interface linearGrowthBaryonsDarkMatter

  ! Tolerance parameter used to ensure times do not exceed that at the Big Crunch.
  double precision                               , parameter               :: timeToleranceRelative    =1.0d-4

  ! Reference wavenumber used when no wavenumber is specified. Small enough that it should be into the regime where baryon
  ! suppression is negligible, but large enough to avoid the BAO region.
  double precision                               , parameter               :: wavenumberReference      =1.0d+0

  ! Indices of tables for baryons and dark matter.
  integer                                        , parameter               :: indexDarkMatter               =1
  integer                                        , parameter               :: indexBaryons                  =2

  ! Lock used for file access.
  type            (lockDescriptor               )                          :: fileLock

  ! Variables used in ODE solving to allow for parallelism.
  double precision                                                         :: wavenumber_
  class           (cosmologyFunctionsClass      ), pointer                 :: cosmologyFunctions_
  class           (intergalacticMediumStateClass), pointer                 :: intergalacticMediumState_
  !$omp threadprivate(wavenumber_, cosmologyFunctions_,intergalacticMediumState_)

contains

  function baryonsDarkMatterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{linearGrowthBaryonsDarkMatter} linear growth class which takes a parameter set as input.
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
    <inputParameter>
      <name>redshiftInitial</name>
      <source>parameters</source>
      <defaultValue>100.0d0</defaultValue>
      <description>The initial redshift from which integration of linear growth should be begin.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftInitialDelta</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>The initial step in redshift used to estimate growth rates of perturbations using finite differencing.</description>
    </inputParameter>
    <inputParameter>
      <name>cambCountPerDecade</name>
      <source>parameters</source>
      <defaultValue>0</defaultValue>
      <description>The number of points per decade of wavenumber to compute in the CAMB transfer function. A value of 0 allows CAMB to choose what it considers to be optimal spacing of wavenumbers.</description>
    </inputParameter>
    <inputParameter>
      <name>darkMatterOnlyInitialConditions</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, set the initial conditions for baryonic modes using the dark matter mode initial conditions.</description>
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
    !!{
    Internal constructor for the \refClass{linearGrowthBaryonsDarkMatter} linear growth class.
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
    double precision                                                          :: timeBigCrunch                  , timeNow
    !![
    <constructorAssign variables="redshiftInitial, redshiftInitialDelta, cambCountPerDecade, darkMatterOnlyInitialConditions, *cosmologyParameters_, *cosmologyParametersInitialConditions_, *cosmologyFunctions_, *intergalacticMediumState_"/>
    !!]

    self%tableInitialized      =.false.
    self%tableTimeMinimum      =1.0d+0
    self%tableTimeMaximum      =2.0d+1
    self%tableWavenumberMinimum=1.0d-3
    self%tableWavenumberMaximum=1.0d+4
    timeBigCrunch              =self%cosmologyFunctions_%timeBigCrunch()
    if (timeBigCrunch > 0.0d0) then
       ! A Big Crunch exists - avoid attempting to tabulate times beyond this epoch.
       if (self%tableTimeMinimum > timeBigCrunch) self%tableTimeMinimum= 0.5d0                       *timeBigCrunch
       if (self%tableTimeMaximum > timeBigCrunch) self%tableTimeMaximum=(1.0d0-timeToleranceRelative)*timeBigCrunch
    end if
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
    self%fileName              =inputPath(pathTypeDataDynamic)                                                       // &
         &                      'largeScaleStructure/'                                                               // &
         &                      self%objectType      (                                                              )// &
         &                      '_'                                                                                  // &
         &                      self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &                      '.hdf5'
    call Directory_Make(File_Path(self%fileName))
    return
  end function baryonsDarkMatterConstructorInternal

  subroutine baryonsDarkMatterDestructor(self)
    !!{
    Destructor for the \refClass{linearGrowthBaryonsDarkMatter} linear growth class.
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
    !!{
    Returns the linear growth factor $D(a)$ for expansion factor {\normalfont \ttfamily aExpansion}, normalized such that
    $D(1)=1$ for a baryons plus dark matter plus cosmological constant cosmology.
    !!}
    use    :: File_Utilities       , only : File_Lock                       , File_Unlock
    use    :: Error                , only : Error_Report
    use    :: Interface_GSL        , only : GSL_Success
    use    :: Interfaces_CAMB      , only : Interface_CAMB_Transfer_Function
    use    :: Numerical_ODE_Solvers, only : odeSolver
    !$ use :: OMP_Lib              , only : omp_lock_kind
    use    :: Table_Labels         , only : extrapolationTypeAbort          , extrapolationTypeFix
    use    :: Tables               , only : table1DGeneric
    implicit none
    class           (linearGrowthBaryonsDarkMatter), intent(inout)              :: self
    double precision                               , intent(in   )              :: time
    double precision                               , intent(in   ), optional    :: wavenumber
    double precision                               , parameter                  :: odeToleranceAbsolute          =   1.0d-10, odeToleranceRelative                = 1.0d-10
    integer                                        , parameter                  :: growthTablePointsPerDecadeTime=1000      , growthTablePointsPerDecadeWavenumber=100
    double precision                               , dimension(4)               :: growthFactorODEVariables
    double precision                               , dimension(2)               :: redshiftsInitial                         , timesInitial
    double precision                               , dimension(:) , allocatable :: linearGrowthFactorPresent
    type            (odeSolver                    )               , allocatable :: solver
    integer                                                                     :: i                                        , j
    double precision                                                            :: growthFactorDerivativeBaryons            , growthFactorDerivativeDarkMatter             , &
         &                                                                         timeNow                                  , wavenumberLogarithmic                        , &
         &                                                                         timePresent                              , timeBigCrunch
    integer                                                                     :: growthTableNumberPoints
    type            (table1DGeneric               )                             :: transferFunctionDarkMatter               , transferFunctionBaryons
    integer                                                                     :: countWavenumbers
    !$ integer      (omp_lock_kind                )                             :: lockBaryons                              , lockDarkMatter

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
       ! Find the time corresponding to our CAMB starting redshift.
       redshiftsInitial(1)=+self%redshiftInitial
       redshiftsInitial(2)=+self%redshiftInitial      &
            &              -self%redshiftInitialDelta
       timesInitial    (1)=+self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftsInitial(1)))
       timesInitial    (2)=+self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftsInitial(2)))
       if (time        < timesInitial(1)) call Error_Report('requested epoch is before the chosen initial epoch'//{introspection:location})
       if (timePresent < timesInitial(1)) call Error_Report('present epoch is before the chosen initial epoch'  //{introspection:location})
       ! Find minimum and maximum times to tabulate.
       self%tableTimeMinimum=timesInitial(1)
       self%tableTimeMaximum=max(self%tableTimeMaximum,max(timePresent,2.0d0*time))
       timeBigCrunch        =self%cosmologyFunctions_%timeBigCrunch()
       if (timeBigCrunch > 0.0d0) then
          ! A Big Crunch exists - avoid attempting to tabulate times beyond this epoch.
          if (self%tableTimeMinimum > timeBigCrunch) call Error_Report('Big Crunch occurs before the chosen initial epoch'//{introspection:location})
          if (self%tableTimeMaximum > timeBigCrunch) self%tableTimeMaximum=(1.0d0-timeToleranceRelative)*timeBigCrunch
       end if
       ! Find minimum and maximum wavenumbers to tabulate.
       if (present(wavenumber)) then
          self%tableWavenumberMinimum=min(self%tableWavenumberMinimum,0.5d0*wavenumber)
          self%tableWavenumberMaximum=max(self%tableWavenumberMaximum,2.0d0*wavenumber)
       end if
       ! Get the initial conditions from CAMB.
       call transferFunctionDarkMatter%destroy()
       call transferFunctionBaryons   %destroy()
       call Interface_CAMB_Transfer_Function(self%cosmologyParametersInitialConditions_,redshiftsInitial,self%tableWavenumberMaximum,self%tableWavenumberMaximum,countPerDecade=self%cambCountPerDecade,transferFunctionDarkMatter=transferFunctionDarkMatter,transferFunctionBaryons=transferFunctionBaryons)
       ! Determine number of points to tabulate.
       growthTableNumberPoints=int(log10(self%tableTimeMaximum/self%tableTimeMinimum)*dble(growthTablePointsPerDecadeTime))
       ! Destroy current table.
       call self%growthFactor%destroy()
       ! Create table.
       countWavenumbers=int(dble(growthTablePointsPerDecadeWavenumber)*log10(self%tableWavenumberMaximum/self%tableWavenumberMinimum))+1
       call self%growthFactor%create(self%tableTimeMinimum,self%tableTimeMaximum,growthTableNumberPoints,self%tableWavenumberMinimum,self%tableWavenumberMaximum,countWavenumbers,tableCount=2,extrapolationTypeX=extrapolationTypeAbort,extrapolationTypeY=extrapolationTypeFix)
       allocate(linearGrowthFactorPresent(countWavenumbers))
       ! Iterate over wavenumber.
       !$ call OMP_Init_Lock(lockBaryons   )
       !$ call OMP_Init_Lock(lockDarkMatter)
       !$omp parallel private(i,j,wavenumberLogarithmic,growthFactorDerivativeDarkMatter,growthFactorDerivativeBaryons,timeNow,growthFactorODEVariables,solver)
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
       !$omp do
       do j=1,countWavenumbers
          wavenumber_          =self%growthFactor%y(j)
          wavenumberLogarithmic=log(wavenumber_)
          ! Solve ODE to get corresponding expansion factors. Initialize with solution from CAMB.
          !$ call OMP_Set_Lock  (lockDarkMatter)
          call    self%growthFactor%populate(exp(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1)),1,j,table=indexDarkMatter)
          growthFactorDerivativeDarkMatter=(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=2)-transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1))*exp(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1))/(timesInitial(2)-timesInitial(1))
          !$ call OMP_Unset_Lock(lockDarkMatter)
          if (self%darkMatterOnlyInitialConditions) then
             !$ call OMP_Set_Lock  (lockDarkMatter)
             call self%growthFactor%populate(exp(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1)),1,j,table=indexDarkMatter)
             growthFactorDerivativeBaryons=(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=2)-transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1))*exp(transferFunctionDarkMatter%interpolate(wavenumberLogarithmic,table=1))/(timesInitial(2)-timesInitial(1))
             !$ call OMP_Unset_Lock(lockDarkMatter)
          else
             !$ call OMP_Set_Lock  (lockBaryons   )
             call self%growthFactor%populate(exp(transferFunctionBaryons   %interpolate(wavenumberLogarithmic,table=1)),1,j,table=indexBaryons   )
             growthFactorDerivativeBaryons=(transferFunctionBaryons   %interpolate(wavenumberLogarithmic,table=2)-transferFunctionBaryons   %interpolate(wavenumberLogarithmic,table=1))*exp(transferFunctionBaryons   %interpolate(wavenumberLogarithmic,table=1))/(timesInitial(2)-timesInitial(1))
             !$ call OMP_Unset_Lock(lockBaryons   )
          end if
          solver=odeSolver(4_c_size_t,growthFactorODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
          do i=2,growthTableNumberPoints
             timeNow                    =self%growthFactor                    %x(i-1                        )
             growthFactorODEVariables(1)=self%growthFactor                    %z(i-1,j,table=indexDarkMatter)
             growthFactorODEVariables(2)=growthFactorDerivativeDarkMatter
             growthFactorODEVariables(3)=self%growthFactor                    %z(i-1,j,table=indexBaryons   )
             growthFactorODEVariables(4)=growthFactorDerivativeBaryons
             call solver             %solve   (timeNow,self%growthFactor%x(i),growthFactorODEVariables                             )
             call self  %growthFactor%populate(                               growthFactorODEVariables(1),i,j,table=indexDarkMatter)
             call self  %growthFactor%populate(                               growthFactorODEVariables(3),i,j,table=indexBaryons   )
             growthFactorDerivativeDarkMatter=growthFactorODEVariables(2)
             growthFactorDerivativeBaryons   =growthFactorODEVariables(4)
          end do
       end do
       !$omp end do
       !![
       <objectDestructor name="cosmologyFunctions_"      />
       <objectDestructor name="intergalacticMediumState_"/>
       !!]
       !$omp barrier
       !$omp single
       ! Get present day growth factor at every wavenumber.
       do j=1,countWavenumbers
          linearGrowthFactorPresent(j)=self%growthFactor%interpolate(timePresent,self%growthFactor%y(j))
       end do
       !$omp end single
       !$omp barrier
       !$omp do
       do j=1,countWavenumbers
          ! Normalize to growth factor of unity at present day.
          do i=1,growthTableNumberPoints
             call self%growthFactor%populate(self%growthFactor%z(i,j,table=indexDarkMatter)/linearGrowthFactorPresent(j),i,j,table=indexDarkMatter)
             call self%growthFactor%populate(self%growthFactor%z(i,j,table=indexBaryons   )/linearGrowthFactorPresent(j),i,j,table=indexBaryons   )
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
      !!{
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
    !!{
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
    !!{
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
    !!{
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
    !!{
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
    !!{
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
    !!{
    Read tabulated data on linear growth factor from file.
    !!}
    use :: Display       , only : displayMessage        , verbosityLevelWorking
    use :: File_Utilities, only : File_Exists
    use :: HDF5_Access   , only : hdf5Access
    use :: IO_HDF5       , only : hdf5Object
    use :: Table_Labels  , only : extrapolationTypeAbort, extrapolationTypeFix
    implicit none
    class           (linearGrowthBaryonsDarkMatter), intent(inout)               :: self
    double precision                               , dimension(:,:), allocatable :: growthFactorDarkMatter, growthFactorBaryons
    type            (hdf5Object                   )                              :: dataFile

    ! Return immediately if the file does not exist.
    if (.not.File_Exists(char(self%fileName))) return
    call displayMessage('reading D(k,t) data from: '//self%fileName,verbosityLevelWorking)
    if (self%tableInitialized) call self%growthFactor%destroy()
    !$ call hdf5Access%set()
    dataFile=hdf5Object(char(self%fileName),overWrite=.false.)
    call dataFile%readDataset  ('growthFactorDarkMatter',                growthFactorDarkMatter)
    call dataFile%readDataset  ('growthFactorBaryons'   ,                growthFactorBaryons   )
    call dataFile%readAttribute('wavenumberMinimum'     ,          self%tableWavenumberMinimum )
    call dataFile%readAttribute('wavenumberMaximum'     ,          self%tableWavenumberMaximum )
    call dataFile%readAttribute('timeMinimum'           ,          self%tableTimeMinimum       )
    call dataFile%readAttribute('timeMaximum'           ,          self%tableTimeMaximum       )
    !$ call hdf5Access%unset()
    call self%growthFactor%create  (                                                                                                               &
         &                                             self%tableTimeMinimum      ,self%tableTimeMaximum      ,size(growthFactorDarkMatter,dim=1), &
         &                                             self%tableWavenumberMinimum,self%tableWavenumberMaximum,size(growthFactorDarkmatter,dim=2), &
         &                          tableCount        =2                                                                                         , &
         &                          extrapolationTypeX=extrapolationTypeAbort                                                                    , &
         &                          extrapolationTypeY=extrapolationTypeFix                                                                        &
         &                         )
    call self%growthFactor%populate(growthFactorDarkMatter,table=indexDarkMatter)
    call self%growthFactor%populate(growthFactorBaryons   ,table=indexBaryons   )
    deallocate(growthFactorDarkMatter)
    deallocate(growthFactorBaryons   )
    self%tableInitialized=.true.
    return
  end subroutine baryonsDarkMatterFileRead

  subroutine baryonsDarkMatterFileWrite(self)
    !!{
    Write tabulated data on linear growth factor to file.
    !!}
    use :: Display    , only : displayMessage, verbosityLevelWorking
    use :: HDF5       , only : hsize_t
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Object
    implicit none
    class(linearGrowthBaryonsDarkMatter), intent(inout) :: self
    type (hdf5Object                   )                :: dataFile

    ! Open the data file.
    call displayMessage('writing D(k,t) data to: '//self%fileName,verbosityLevelWorking)
    !$ call hdf5Access%set()
    dataFile=hdf5Object(char(self%fileName),overWrite=.true.,chunkSize=100_hsize_t,compressionLevel=9)
    call dataFile%writeDataset  (reshape(self%growthFactor          %zs(table=indexDarkMatter),[self%growthFactor%size(dim=1),self%growthFactor%size(dim=2)]),          'growthFactorDarkMatter'                       )
    call dataFile%writeDataset  (reshape(self%growthFactor          %zs(table=indexBaryons   ),[self%growthFactor%size(dim=1),self%growthFactor%size(dim=2)]),          'growthFactorBaryons'                          )
    call dataFile%writeAttribute(        self%tableWavenumberMinimum                                                                                         ,          'wavenumberMinimum'                            )
    call dataFile%writeAttribute(        self%tableWavenumberMaximum                                                                                         ,          'wavenumberMaximum'                            )
    call dataFile%writeAttribute(        self%tableTimeMinimum                                                                                               ,          'timeMinimum'                                  )
    call dataFile%writeAttribute(        self%tableTimeMaximum                                                                                               ,          'timeMaximum'                                  )
    !$ call hdf5Access%unset()
    return
  end subroutine baryonsDarkMatterFileWrite
