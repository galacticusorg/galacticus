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
  A spherical collapse solver class for universes consisting of collisionless matter and a cosmological constant.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
  use :: Linear_Growth      , only : linearGrowth      , linearGrowthClass

  !![
  <sphericalCollapseSolver name="sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt">
   <description>A spherical collapse solver for universes consisting of collisionless matter and a cosmological constant.</description>
  </sphericalCollapseSolver>
  !!]
  type, extends(sphericalCollapseSolverClass) :: sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt
     !!{
     A spherical collapse solver for universes consisting of collisionless matter and a cosmological constant.
     !!}
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_         => null()
     class(linearGrowthClass      ), pointer :: linearGrowth_               => null()
     type (varying_string         )          :: fileNameCriticalOverdensity          , fileNameVirialDensityContrast, &
          &                                     fileNameRadiusTurnaround
   contains
     !![
     <methods>
       <method description="Restore a tabulated solution from file." method="restoreTable"/>
       <method description="Store a tabulated solution to file."     method="storeTable"  />
       <method description="Construct a tabulated solution."         method="tabulate"    />
     </methods>
     !!]
     final     ::                          cllsnlssMttCsmlgclCnstntDestructor
     procedure :: criticalOverdensity   => cllsnlssMttCsmlgclCnstntCriticalOverdensity
     procedure :: virialDensityContrast => cllsnlssMttCsmlgclCnstntVirialDensityContrast
     procedure :: radiusTurnaround      => cllsnlssMttCsmlgclCnstntRadiusTurnaround
     procedure :: linearNonlinearMap    => cllsnlssMttCsmlgclCnstntLinearNonlinearMap
     procedure :: restoreTable          => cllsnlssMttCsmlgclCnstntRestoreTable
     procedure :: storeTable            => cllsnlssMttCsmlgclCnstntStoreTable
     procedure :: tabulate              => cllsnlssMttCsmlgclCnstntTabulate
  end type sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt

  interface sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt
     !!{
     Constructors for the \refClass{sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt} spherical collapse solver class.
     !!}
     module procedure cllsnlssMttCsmlgclCnstntConstructorParameters
     module procedure cllsnlssMttCsmlgclCnstntConstructorInternal
  end interface sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt

  ! Enumeration of types of calculation to perform.
  !![
  <enumeration>
   <name>cllsnlssMttCsmlgclCnstntClcltn</name>
   <description>Enumeration of calculation types to be performed by the spherical collapse solver.</description>
   <entry label="criticalOverdensity"  />
   <entry label="virialDensityContrast"/>
   <entry label="radiusTurnaround"     />
  </enumeration>
  !!]

  ! Resolution of tabulated solutions.
  integer         , parameter :: tablePointsPerDecade=1000

  ! Variables used in root finding.
  double precision            :: OmegaDarkEnergyEpochal, OmegaMatterEpochal, &
       &                         amplitudePerturbation , hubbleTimeEpochal , &
       &                         time_                 , timeTarget        , &
       &                         radiusMaximum
  !$omp threadprivate(OmegaDarkEnergyEpochal,OmegaMatterEpochal,amplitudePerturbation,hubbleTimeEpochal,time_,timeTarget,radiusMaximum)

contains

  function cllsnlssMttCsmlgclCnstntConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt} spherical collapse solver class that takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)                :: self
    type (inputParameters                                 ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                         ), pointer       :: cosmologyFunctions_
    class(linearGrowthClass                               ), pointer       :: linearGrowth_

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="linearGrowth"       name="linearGrowth_"       source="parameters"/>
    !!]
    self=sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt(cosmologyFunctions_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="linearGrowth_"      />
    !!]
    return
  end function cllsnlssMttCsmlgclCnstntConstructorParameters

  function cllsnlssMttCsmlgclCnstntConstructorInternal(cosmologyFunctions_,linearGrowth_) result(self)
    !!{
    Internal constructor for the \refClass{sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt} spherical collapse solver class.
    !!}
    use :: Input_Paths       , only : inputPath   , pathTypeDataDynamic
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    type (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)                                  :: self
    class(cosmologyFunctionsClass                         ), intent(in   ), target           :: cosmologyFunctions_
    class(linearGrowthClass                               ), intent(in   ), target, optional :: linearGrowth_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *linearGrowth_"/>
    !!]
    
    self%fileNameCriticalOverdensity  =inputPath(pathTypeDataDynamic)                                                       // &
         &                             'largeScaleStructure/'                                                               // &
         &                             self%objectType      (                                                              )// &
         &                             'CriticalOverdensity_'                                                               // &
         &                             self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &                             '.hdf5'
    self%fileNameVirialDensityContrast=inputPath(pathTypeDataDynamic)                                                       // &
         &                             'largeScaleStructure/'                                                               // &
         &                             self%objectType      (                                                              )// &
         &                             'VirialDensityContrast_'                                                             // &
         &                             self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &                             '.hdf5'
    self%fileNameRadiusTurnaround     =inputPath(pathTypeDataDynamic)                                                       // &
         &                             'largeScaleStructure/'                                                               // &
         &                             self%objectType      (                                                              )// &
         &                             'TurnaroundRadius_'                                                                  // &
         &                             self%hashedDescriptor(includeSourceDigest=.true.,includeFileModificationTimes=.true.)// &
         &                             '.hdf5'
    return
  end function cllsnlssMttCsmlgclCnstntConstructorInternal

  subroutine cllsnlssMttCsmlgclCnstntDestructor(self)
    !!{
    Destructor for the \refClass{sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt} spherical collapse solver class.
    !!}
    implicit none
    type(sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    <objectDestructor name="self%linearGrowth_"      />
    !!]
    return
  end subroutine cllsnlssMttCsmlgclCnstntDestructor

  subroutine cllsnlssMttCsmlgclCnstntCriticalOverdensity(self,time,tableStore,criticalOverdensity_)
    !!{
    Compute the critical overdensity for collapse for the spherical collapse model.
    !!}
    use :: Error         , only : errorStatusSuccess
    use :: File_Utilities, only : File_Lock         , File_Unlock, lockDescriptor, Directory_Make, &
         &                        File_Path
    use :: Tables        , only : table1D
    implicit none
    class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)             , intent(inout) :: self
    double precision                                                               , intent(in   ) :: time
    logical                                                                        , intent(in   ) :: tableStore
    class           (table1D                                         ), allocatable, intent(inout) :: criticalOverdensity_
    integer                                                                                        :: status
    type            (lockDescriptor                                  )                             :: fileLock

    if (tableStore) then
       call Directory_Make(char(File_Path(char(self%fileNameCriticalOverdensity)))                             )
       call File_Lock     (               char(self%fileNameCriticalOverdensity)  ,fileLock,lockIsShared=.false.)
    end if
    call    self%restoreTable(time,criticalOverdensity_,self%fileNameCriticalOverdensity                 ,tableStore,status)
    if (status /= errorStatusSuccess) then
       call self%tabulate    (time,criticalOverdensity_,cllsnlssMttCsmlgclCnstntClcltnCriticalOverdensity                  )
       call self%storeTable  (     criticalOverdensity_,self%fileNameCriticalOverdensity                 ,tableStore       )
    end if
    if (tableStore) call File_Unlock(fileLock)
    return
  end subroutine cllsnlssMttCsmlgclCnstntCriticalOverdensity

  subroutine cllsnlssMttCsmlgclCnstntVirialDensityContrast(self,time,tableStore,virialDensityContrast_)
    !!{
    Tabulate the virial density contrast for the spherical collapse model.
    !!}
    use :: Error         , only : errorStatusSuccess
    use :: File_Utilities, only : File_Lock         , File_Unlock, lockDescriptor, Directory_Make, &
         &                        File_Path
    use :: Tables        , only : table1D
    implicit none
    class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)             , intent(inout) :: self
    double precision                                                               , intent(in   ) :: time
    logical                                                                        , intent(in   ) :: tableStore
    class           (table1D                                         ), allocatable, intent(inout) :: virialDensityContrast_
    integer                                                                                        :: status
    type            (lockDescriptor                                  )                             :: fileLock

    if (tableStore) then
       call Directory_Make(char(File_Path(char(self%fileNameVirialDensityContrast)))                              )
       call File_Lock     (               char(self%fileNameVirialDensityContrast)  ,fileLock,lockIsShared=.false.)
    end if
    call    self%restoreTable(time,virialDensityContrast_,self%fileNameVirialDensityContrast                 ,tableStore,status)
    if (status /= errorStatusSuccess) then
       call self%tabulate    (time,virialDensityContrast_,cllsnlssMttCsmlgclCnstntClcltnVirialDensityContrast                  )
       call self%storeTable  (     virialDensityContrast_,self%fileNameVirialDensityContrast                 ,tableStore       )
    end if
    if (tableStore) call File_Unlock(fileLock)
    return
  end subroutine cllsnlssMttCsmlgclCnstntVirialDensityContrast

  subroutine cllsnlssMttCsmlgclCnstntRadiusTurnaround(self,time,tableStore,radiusTurnaround_)
    !!{
    Tabulate the ratio of turnaround to virial radii for the spherical collapse model.
    !!}
    use :: Error         , only : errorStatusSuccess
    use :: File_Utilities, only : File_Lock         , File_Unlock, lockDescriptor, Directory_Make, &
         &                        File_Path
    use :: Tables        , only : table1D
    implicit none
    class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)             , intent(inout) :: self
    double precision                                                               , intent(in   ) :: time
    logical                                                                        , intent(in   ) :: tableStore
    class           (table1D                                         ), allocatable, intent(inout) :: radiusTurnaround_
    integer                                                                                        :: status
    type            (lockDescriptor                                  )                             :: fileLock

    if (tableStore) then
       call Directory_Make(char(File_Path(char(self%fileNameRadiusTurnaround)))                             )
       call File_Lock     (               char(self%fileNameRadiusTurnaround)  ,fileLock,lockIsShared=.false.)
    end if
    call    self%restoreTable(time,radiusTurnaround_,self%fileNameRadiusTurnaround                 ,tableStore,status)
    if (status /= errorStatusSuccess) then
       call self%tabulate    (time,radiusTurnaround_,cllsnlssMttCsmlgclCnstntClcltnRadiusTurnaround                  )
       call self%storeTable  (     radiusTurnaround_,self%fileNameRadiusTurnaround                 ,tableStore       )
    end if
    if (tableStore) call File_Unlock(fileLock)
    return
  end subroutine cllsnlssMttCsmlgclCnstntRadiusTurnaround

  subroutine cllsnlssMttCsmlgclCnstntTabulate(self,time,sphericalCollapse_,calculationType)
    !!{
    Tabulate spherical collapse solutions for $\delta_\mathrm{crit}$, $\Delta_\mathrm{vir}$, or $R_\mathrm{ta}/R_\mathrm{vir}$ vs. time.
    !!}
    use :: Cosmology_Functions, only : timeToleranceRelativeBigCrunch
    use :: Error              , only : Error_Report
    use :: Kind_Numbers       , only : kind_dble
    use :: Linear_Growth      , only : normalizeMatterDominated
    use :: Root_Finder        , only : rangeExpandMultiplicative     , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: Tables             , only : table1DLogarithmicLinear
    implicit none
    class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)             , intent(inout) :: self
    double precision                                                               , intent(in   ) :: time
    type            (enumerationCllsnlssMttCsmlgclCnstntClcltnType   )             , intent(in   ) :: calculationType
    class           (table1D                                         ), allocatable, intent(inout) :: sphericalCollapse_
    double precision                                                  , parameter                  :: toleranceAbsolute         =0.0d+0 , toleranceRelative         =1.0d-9
    type            (rootFinder                                      ), save                       :: finder
    logical                                                           , save                       :: finderConstructed         =.false.
    !$omp threadprivate(finder,finderConstructed)
    integer                                                                                        :: countTimes                        , i
    double precision                                                                               :: expansionFactor                   , epsilonPerturbation              , &
         &                                                                                            epsilonPerturbationMaximum        , epsilonPerturbationMinimum       , &
         &                                                                                            eta                               , normalization                    , &
         &                                                                                            radiiRatio                        , radiusMaximum                    , &
         &                                                                                            timeBigCrunch                     , timeMinimum                      , &
         &                                                                                            timeMaximum
    double complex                                                                                 :: a                                 , b                                , &
         &                                                                                            c                                 , d                                , &
         &                                                                                            Delta

    ! Find minimum and maximum times to tabulate.
    if (allocated(sphericalCollapse_)) then
       ! Use currently tabulated range as the starting point.
       timeMinimum=sphericalCollapse_%x(+1)
       timeMaximum=sphericalCollapse_%x(-1)
    else
       ! Specify an initial default range.
       timeMinimum= 1.0d0
       timeMaximum=20.0d0
    end if
    ! Expand the range to ensure the requested time is included.
    timeMinimum=min(timeMinimum,time/2.0d0)
    timeMaximum=max(timeMaximum,time*2.0d0)
    timeBigCrunch=self%cosmologyFunctions_%timeBigCrunch()
    ! If a Big Crunch exists - avoid attempting to tabulate times beyond this epoch.
    if (timeBigCrunch > 0.0d0) then
       if (timeMinimum > timeBigCrunch) timeMinimum= 0.5d0                                *timeBigCrunch
       if (timeMaximum > timeBigCrunch) timeMaximum=(1.0d0-timeToleranceRelativeBigCrunch)*timeBigCrunch
    end if
    ! Determine number of points to tabulate.
    countTimes=int(log10(timeMaximum/timeMinimum)*dble(tablePointsPerDecade))
    ! Deallocate table if currently allocated.
    if (allocated(sphericalCollapse_)) then
       call sphericalCollapse_%destroy()
       deallocate(sphericalCollapse_)
    end if
    allocate(table1DLogarithmicLinear :: sphericalCollapse_)
    select type (sphericalCollapse_)
    type is (table1DLogarithmicLinear)
       call sphericalCollapse_%create(timeMinimum,timeMaximum,countTimes)
       ! Solve ODE to get corresponding overdensities.
       do i=1,countTimes
          time_=sphericalCollapse_%x(i)
          ! Get the current expansion factor.
          expansionFactor=self%cosmologyFunctions_%expansionFactor(time_)
          ! Determine the largest (i.e. least negative) value of ε for which a perturbation can collapse.
          if (self%cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor) > 0.0d0) then
             epsilonPerturbationMaximum=-(                                                                                     &
                  &                       +27.0d0                                                                              &
                  &                       / 4.0d0                                                                              &
                  &                       *self%cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor)    &
                  &                       *self%cosmologyFunctions_%omegaMatterEpochal    (expansionFactor=expansionFactor)**2 &
                  &                      )**(1.0d0/3.0d0)
          else
             epsilonPerturbationMaximum=-1.0d-6
          end if
          ! Estimate a suitably negative minimum value for ε.
          epsilonPerturbationMinimum=-10.0d0
          ! Compute epochal cosmological parameters.
          OmegaMatterEpochal    =    self%cosmologyFunctions_%omegaMatterEpochal    (expansionFactor=expansionFactor)
          OmegaDarkEnergyEpochal=    self%cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor)
          hubbleTimeEpochal     =abs(self%cosmologyFunctions_%expansionRate         (                expansionFactor))
          ! Find the value of ε for which the perturbation just collapses at this time.
          if (.not.finderConstructed) then
             finder=rootFinder(                                                                                &
                  &            rootFunction                 =cllsnlssMttCsmlgclCnstntPerturbationCollapseRoot, &
                  &            toleranceAbsolute            =toleranceAbsolute                               , &
                  &            toleranceRelative            =toleranceRelative                               , &
                  &            rangeExpandUpward            =0.5d0                                           , &
                  &            rangeExpandDownward          =2.0d0                                           , &
                  &            rangeExpandType              =rangeExpandMultiplicative                       , &
                  &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive                   , &
                  &            rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative                     &
                  &           )
             finderConstructed=.true.
          end if
          epsilonPerturbation=finder%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
          ! Compute the required quantity for this perturbation.
          select case (calculationType%ID)
          case (cllsnlssMttCsmlgclCnstntClcltnCriticalOverdensity%ID)
             ! Critical linear overdensity.
             if (.not.associated(self%linearGrowth_)) call Error_Report('no linearGrowth object was supplied'//{introspection:location})
             normalization=self%linearGrowth_%value(time_,normalize=normalizeMatterDominated)/expansionFactor
             call sphericalCollapse_%populate(                           &
                  &                           +normalization             &
                  &                           *0.6d0                     &
                  &                           *(                         &
                  &                             +1.0d0                   &
                  &                             -OmegaMatterEpochal      &
                  &                             -OmegaDarkEnergyEpochal  &
                  &                             -epsilonPerturbation     &
                  &                            )                         &
                  &                           /OmegaMatterEpochal      , &
                  &                           i                          &
                  &                           )
          case (cllsnlssMttCsmlgclCnstntClcltnVirialDensityContrast%ID,cllsnlssMttCsmlgclCnstntClcltnRadiusTurnaround%ID)
             ! Compute the maximum radius of the perturbation.
             radiusMaximum=+cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum(epsilonPerturbation)
             ! Find the η-factor (see Lahav et al. 1991) which measures the dark energy contribution to the energy of the
             ! perturbation.
             eta          =+2.0d0                    &
                  &        *(                        &
                  &          +OmegaDarkEnergyEpochal &
                  &          /OmegaMatterEpochal     &
                  &         )*radiusMaximum**3
             ! Handle the open universe case separately.
             if (OmegaDarkEnergyEpochal == 0.0d0) then
                radiiRatio=0.5d0
             else
                ! Coefficients of the cubic energy equation.
                a=cmplx(  2.0d0*eta ,0.0d0,kind=kind_dble)
                b=cmplx(  0.0d0     ,0.0d0,kind=kind_dble)
                c=cmplx(-(2.0d0+eta),0.0d0,kind=kind_dble)
                d=cmplx(  1.0d0     ,0.0d0,kind=kind_dble)
                ! Solve the cubic energy equation.
                Delta     =+(                                        &
                     &       +sqrt( 3.0d0                          ) &
                     &       *sqrt(27.0d0*a**4*d**2+4.0d0*a**3*c**3) &
                     &       -      9.0d0*a**2*d                     &
                     &      )**(1.0d0/3.0d0)
                radiiRatio=+real(                                                                                                        &
                     &           +cmplx(1.0d0,-sqrt(3.0d0),kind=kind_dble)*c            /2.0d0**(2.0d0/3.0d0)/3.0d0**(1.0d0/3.0d0)/Delta &
                     &           -cmplx(1.0d0,+sqrt(3.0d0),kind=kind_dble)*Delta/2.0d0/a/2.0d0**(1.0d0/3.0d0)/3.0d0**(2.0d0/3.0d0)       &
                     &          )
             end if
             select case (calculationType%ID)
             case (cllsnlssMttCsmlgclCnstntClcltnVirialDensityContrast%ID)
                call sphericalCollapse_%populate(                                                                                              &
                     &                           1.0d0/(radiiRatio*cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum(epsilonPerturbation))**3, &
                     &                           i                                                                                             &
                     &                          )
             case (cllsnlssMttCsmlgclCnstntClcltnRadiusTurnaround%ID)
                call sphericalCollapse_%populate(                                                                                              &
                     &                           1.0d0/ radiiRatio                                                                           , &
                     &                           i                                                                                             &
                     &                          )
             end select
          end select
       end do
    end select
    return
  end subroutine cllsnlssMttCsmlgclCnstntTabulate

  double precision function cllsnlssMttCsmlgclCnstntPerturbationCollapseRoot(epsilonPerturbation)
    !!{
    Root function used to determine when the collapse time for a perturbation of amplitude {\normalfont \ttfamily
    epsilonPerturbation} collapses at the current time.
    !!}
    implicit none
    double precision, intent(in   ) :: epsilonPerturbation

    ! Evaluate the root function.
    cllsnlssMttCsmlgclCnstntPerturbationCollapseRoot=time_Collapse(epsilonPerturbation)-time_
    return
  end function cllsnlssMttCsmlgclCnstntPerturbationCollapseRoot

  double precision function cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum(epsilonPerturbation)
    !!{
    Find the maximum radius of a perturbation with initial curvature {\normalfont \ttfamily epsilonPerturbation}.
    !!}
    use :: Root_Finder, only : rootFinder
    implicit none
    double precision            , intent(in   ) :: epsilonPerturbation
    double precision            , parameter     :: toleranceAbsolute               =0.0d0  , toleranceRelative               =1.0d-9
    type            (rootFinder), save          :: finder
    logical                     , save          :: finderConstructed               =.false.
    !$omp threadprivate(finder,finderConstructed)
    double precision                            :: expansionFactorTurnaroundMaximum        , expansionFactorTurnaroundMinimum

    if (OmegaDarkEnergyEpochal == 0.0d0) then
       ! No cosmological constant - use the simple analytic solution.
       cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum=-OmegaMatterEpochal &
            &                                            /epsilonPerturbation
    else
       ! Cosmological constant - use root finder.
       amplitudePerturbation=+epsilonPerturbation
       expansionFactorTurnaroundMinimum =-OmegaMatterEpochal       &
            &                            /amplitudePerturbation
       expansionFactorTurnaroundMaximum =+(                        &
            &                              +OmegaMatterEpochal     &
            &                              /OmegaDarkEnergyEpochal &
            &                              /2.0d0                  &
            &                             )**(1.0d0/3.0d0)
       if      (cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximumRoot(expansionFactorTurnaroundMaximum) > 0.0d0) then
          ! If the root function is not negative at the upper limit for the expansion factor at turnaround it is due to rounding
          ! errors in the calculation of the upper limit for the expansion factor at turnaround which implies that the upper limit
          ! for the expansion factor at turnaround is very close to the actual root.
          cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum=expansionFactorTurnaroundMaximum
       else if (cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximumRoot(expansionFactorTurnaroundMinimum) < 0.0d0)  then
          ! If the root function is not positive at the lower limit for the expansion factor at turnaround it is due to rounding
          ! errors in the calculation of the lower limit for the expansion factor at turnaround which implies that the lower limit
          ! for the expansion factor at turnaround is very close to the actual root.
          cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum=expansionFactorTurnaroundMinimum
       else
          if (.not.finderConstructed) then
             finder=rootFinder(                                                                         &
                  &            rootFunction     =cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximumRoot, &
                  &            toleranceAbsolute=toleranceAbsolute                                    , &
                  &            toleranceRelative=toleranceRelative                                      &
                  &           )
             finderConstructed=.true.
          end if
          cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum=finder%find(rootRange=[expansionFactorTurnaroundMinimum,expansionFactorTurnaroundMaximum])
       end if
    end if
    return
  end function cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum

  double precision function cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximumRoot(radiusMaximum)
    !!{
    Function used in root finding to determine the maximum expansion radius of the perturbation. Evaluates the expansion speed
    of the perturbation which must be zero at the maximum radius.
    !!}
    double precision, intent(in   ) :: radiusMaximum

    cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximumRoot=+OmegaMatterEpochal     &
         &                                                /radiusMaximum                                  &
         &                                                +amplitudePerturbation  &
         &                                                +OmegaDarkEnergyEpochal &
         &                                                *radiusMaximum**2
    return
  end function cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximumRoot

  double precision function time_Collapse(epsilonPerturbation)
    use :: Numerical_Integration, only : integrator
    implicit none
    real   (c_double  ), intent(in   ) :: epsilonPerturbation
    type   (integrator)                :: integrator_
    real   (c_double  ), parameter     :: radiusMinimum        =0.0d0
    real   (c_double  ), parameter     :: numericalLimitEpsilon=1.0d-4
    real   (c_double  )                :: radiusUpperNumerical        , radiusTurnaround, &
         &                                timeTurnaround

    ! Find the maximum radius of the perturbation.
    radiusTurnaround       =+cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum(epsilonPerturbation)
    ! Compute maximum value of a for numerical integration.
    radiusUpperNumerical=+(                       &
         &                 +1.0d0                 &
         &                 -numericalLimitEpsilon &
         &                )                       &
         &               *radiusTurnaround
    ! Share the ε parameter.
    amplitudePerturbation=epsilonPerturbation
    ! Integrate the perturbation equation from size zero to maximum size to get the time to maximum expansion.
    integrator_   = integrator           (cllsnlssMttCsmlgclCnstntPerturbationIntegrand,toleranceRelative   =1.0d-6,hasSingularities=.true.)
    timeTurnaround=+integrator_%integrate(radiusMinimum                                ,radiusUpperNumerical                               ) &
         &         /                      hubbleTimeEpochal
    ! Add on analytic correction for region close to the turnaround radius.
    timeTurnaround=+timeTurnaround                  &
         &         -2.0d0                           &
         &         *sqrt(                           &
         &               +OmegaMatterEpochal        &
         &               /radiusUpperNumerical      &
         &               +epsilonPerturbation       &
         &               +OmegaDarkEnergyEpochal    &
         &               *radiusUpperNumerical  **2 &
         &              )                           &
         &         /    (                           &
         &               +2.0d0                     &
         &               *OmegaDarkEnergyEpochal    &
         &               *radiusUpperNumerical      &
         &               -OmegaMatterEpochal        &
         &               /radiusUpperNumerical  **2 &
         &              )                           &
         &         /      hubbleTimeEpochal
    ! Time to collapse is twice the time to maximum expansion.
    time_Collapse=+2.0d0          &
         &        *timeTurnaround
    return
  end function time_Collapse

  double precision function cllsnlssMttCsmlgclCnstntPerturbationIntegrand(radius)
    !!{
    Integrand function giving the $\mathrm{d}t/\mathrm{d}r$ for the perturbation.
    !!}
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: sqrtArgument

    ! Compute the integrand.
    sqrtArgument=+OmegaMatterEpochal               &
         &       +amplitudePerturbation *radius    &
         &       +OmegaDarkEnergyEpochal*radius**3
    if (sqrtArgument > 0.0d0) then
       cllsnlssMttCsmlgclCnstntPerturbationIntegrand=+sqrt(              &
            &                                              +radius       &
            &                                              /sqrtArgument &
            &                                             )
    else
       cllsnlssMttCsmlgclCnstntPerturbationIntegrand=+0.0d0
    end if
    return
  end function cllsnlssMttCsmlgclCnstntPerturbationIntegrand

  subroutine cllsnlssMttCsmlgclCnstntLinearNonlinearMap(self,time,linearNonlinearMap_)
    !!{
    Tabulate the mapping between linear and nonlinear overdensity for the spherical collapse model.
    !!}
    use :: Array_Utilities      , only : Array_Reverse
    use :: Arrays_Search        , only : searchArray
    use :: Error                , only : Error_Report
    use :: Linear_Growth        , only : normalizeMatterDominated
    use :: Numerical_Integration, only : integrator
    use :: Numerical_Ranges     , only : Make_Range               , rangeTypeLinear              , rangeTypeLogarithmic
    use :: Root_Finder          , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt), intent(inout)              :: self
    double precision                                                  , intent(in   )              :: time
    class           (table2DLinLinLin                                ), intent(inout)              :: linearNonlinearMap_
    integer                                                           , parameter                  :: tableIncrement               =100
    integer                                                           , parameter                  :: timesPerDecade               = 10
    integer                                                           , parameter                  :: overdensityLinearCount       =500
    double precision                                                  , parameter                  :: numericalLimitEpsilon        =  1.0d-4
    double precision                                                  , parameter                  :: toleranceAbsolute            =  0.0d+0, toleranceRelative         =1.0d-9
    double precision                                                  , parameter                  :: expansionFactorMinimum       =  1.0d-2, expansionFactorMaximum    =1.0d+0
    double precision                                                  , parameter                  :: overdensityLinearMinimum     = -5.0d+0, overdensityLinearMaximum  =2.0d+0
    double precision                                                  , allocatable, dimension(:)  :: overdensityLinear                     , overdensityNonlinear              , &
         &                                                                                            overdensityLinearTmp                  , overdensityNonLinearTmp           , &
         &                                                                                            times                                 , overdensitiesLinear
    double precision                                                                               :: expansionFactor                       , epsilonPerturbationMaximum        , &
         &                                                                                            epsilonPerturbationCollapsed          , radiusNow                         , &
         &                                                                                            epsilonPerturbation                   , epsilonPerturbationMinimum        , &
         &                                                                                            timeMaximum                           , radiusUpperLimit                  , &
         &                                                                                            normalization                         , overdensityNonlinear_             , &
         &                                                                                            timesMinimum                          , timesMaximum
    logical                                                                                        :: finderPerturbationConstructed         , finderRadiusConstructed
    type            (integrator                                      )                             :: integrator_
    type            (rootFinder                                      )                             :: finderPerturbation                    , finderRadius
    integer                                                                                        :: i                                     , timeCount                         , &
         &                                                                                            iOverdensityLinear                    , iOverdensity                      , &
         &                                                                                            iTime

    ! Check that we have a linear growth object.
    if (.not.associated(self%linearGrowth_)) call Error_Report('no linearGrowth object was supplied'//{introspection:location})
    ! Set initial state of root finder objects.
    finderPerturbationConstructed=.false. 
    finderRadiusConstructed      =.false.
    ! Find a suitable range of times to tabulate, and generate an array of times.
    timesMinimum      =min(0.5d0*time,self%cosmologyFunctions_%cosmicTime(expansionFactorMinimum))
    timesMaximum      =max(2.0d0*time,self%cosmologyFunctions_%cosmicTime(expansionFactorMaximum))
    timeCount         =int(log10(timesMaximum/timesMinimum)*dble(timesPerDecade))+1
    times             =Make_Range(timesMinimum,timesMaximum,timeCount,rangeTypeLogarithmic)
    ! Generate a range of linear overdensities.
    overdensitiesLinear=Make_Range(overdensityLinearMinimum,overdensityLinearMaximum,overdensityLinearCount,rangeTypeLinear)
    ! Create the table.
    call linearNonlinearMap_%create(overdensitiesLinear,times)
    ! Construct an integrator.
    integrator_=integrator(cllsnlssMttCsmlgclCnstntPerturbationIntegrand,toleranceRelative=1.0d-6,hasSingularities=.true.)
    ! Iterate over times.
    do iTime=1,timeCount
       ! Get the current expansion factor.
       time_           =times(iTime)
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time_)
       ! Determine the largest (i.e. least negative) value of ε for which a perturbation can collapse.
       if (self%cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor) > 0.0d0) then
          epsilonPerturbationMaximum=-(                                                                                     &
               &                       +27.0d0                                                                              &
               &                       / 4.0d0                                                                              &
               &                       *self%cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor)    &
               &                       *self%cosmologyFunctions_%omegaMatterEpochal    (expansionFactor=expansionFactor)**2 &
               &                      )**(1.0d0/3.0d0)
       else
          epsilonPerturbationMaximum=-1.0d-6
       end if
       ! Estimate a suitably negative minimum value for ε.
       epsilonPerturbationMinimum=-10.0d0
       ! Compute cosmological parameters at this epoch.
       OmegaMatterEpochal    =self%cosmologyFunctions_%omegaMatterEpochal    (expansionFactor=expansionFactor)
       OmegaDarkEnergyEpochal=self%cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor)
       hubbleTimeEpochal     =self%cosmologyFunctions_%expansionRate         (                expansionFactor)
       ! Find the value of ε for which the perturbation just collapses at this time.
       if (.not.finderPerturbationConstructed) then
          finderPerturbation=rootFinder(                                                                                &
               &                        rootFunction                 =cllsnlssMttCsmlgclCnstntPerturbationCollapseRoot, &
               &                        toleranceAbsolute            =toleranceAbsolute                               , &
               &                        toleranceRelative            =toleranceRelative                               , &
               &                        rangeExpandUpward            =0.5d0                                           , &
               &                        rangeExpandDownward          =2.0d0                                           , &
               &                        rangeExpandType              =rangeExpandMultiplicative                       , &
               &                        rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive                   , &
               &                        rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative                     &
               &                       )
          finderPerturbationConstructed=.true.
       end if
       epsilonPerturbationCollapsed=finderPerturbation%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
       ! For non-collapsed regions, ε will be greater then that for a collapsed perturbation. Step through values until
       ! sufficiently low non-linear overdensity is reached.
       epsilonPerturbation=epsilonPerturbationCollapsed
       i=0
       do while (.true.)
          i                  =i                  +1
          epsilonPerturbation=epsilonPerturbation+1.0d-2*abs(epsilonPerturbationCollapsed)
          ! Share the ε parameter.
          amplitudePerturbation=epsilonPerturbation
          ! For collapsing perturbations, find the time of maximum radius.
          if (epsilonPerturbation > epsilonPerturbationMaximum) then
             ! This perturbation will not collapse. Maximum radius is reached at infinite time.
             radiusMaximum=huge(1.0d0)
             timeTarget   =     time_
          else
             ! This perturbation will collapse. Find the maximum radius.
             radiusMaximum=cllsnlssMttCsmlgclCnstntRadiusPerturbationMaximum(epsilonPerturbation)
             ! Compute maximum value of a for numerical integration.
             radiusUpperLimit=(1.0d0-numericalLimitEpsilon)*radiusMaximum
             ! Integrate the perturbation equation from size zero to maximum size to get the time to maximum expansion, adding on the
             ! analytic correction for the region close to maximum expansion.
             timeMaximum     =+integrator_%integrate(0.0d0,radiusUpperLimit) &
                  &           /hubbleTimeEpochal                             &
                  &           -2.0d0                                         &
                  &           *sqrt(                                         &
                  &                 +OmegaMatterEpochal                      &
                  &                 /radiusUpperLimit                        &
                  &                 +epsilonPerturbation                     &
                  &                 +OmegaDarkEnergyEpochal                  &
                  &                 *radiusUpperLimit      **2               &
                  &                )                                         &
                  &           /(                                             &
                  &             +2.0d0                                       &
                  &             *OmegaDarkEnergyEpochal                      &
                  &             *radiusUpperLimit                            &
                  &             -OmegaMatterEpochal                          &
                  &             /radiusUpperLimit          **2               &
                  &            )                                             &
                  &           /hubbleTimeEpochal
             ! Set the target time
             if (timeMaximum > time_) then
                ! Expanding phase.
                timeTarget=+      time_
             else
                ! Collapsing phase.
                timeTarget=+2.0d0*timeMaximum      &
                     &                 -      time_
             end if
          end if
          ! Solve for the radius at the present time.
          if (.not.finderRadiusConstructed) then
             finderRadius=rootFinder(                                                                  &
                  &                  rootFunction                 =cllsnlssMttCsmlgclCnstntRadiusRoot, &
                  &                  toleranceAbsolute            =toleranceAbsolute                 , &
                  &                  toleranceRelative            =toleranceRelative                 , &
                  &                  rangeExpandDownward          =0.5d0                             , &
                  &                  rangeExpandUpward            =2.0d0                             , &
                  &                  rangeExpandType              =rangeExpandMultiplicative         , &
                  &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative     , &
                  &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive       &
                  &                 )
             finderRadiusConstructed=.true.
          end if
          if (epsilonPerturbation <= epsilonPerturbationMaximum .and. cllsnlssMttCsmlgclCnstntRadiusRoot(radiusMaximum) < 0.0d0) then
             ! Perturbation is close to maximum radius. Adopt this as the solution.
             radiusNow=radiusMaximum
          else
             ! Find the current radius.
             radiusNow=finderRadius%find(rootGuess=1.0d0)
          end if
          normalization=+self%linearGrowth_%value(time_,normalize=normalizeMatterDominated) &
               &        /                         expansionFactor
          if (.not.allocated(overdensityLinear)) then
             allocate(overdensityLinear   (tableIncrement))
             allocate(overdensityNonLinear(tableIncrement))
          else if (i > size(overdensityLinear)) then
             call move_alloc(overdensityLinear   ,overdensityLinearTmp   )
             call move_alloc(overdensityNonLinear,overdensityNonLinearTmp)
             allocate(overdensityLinear   (size(overdensityLinearTmp   )+tableIncrement))
             allocate(overdensityNonLinear(size(overdensityNonLinearTmp)+tableIncrement))
             overdensityLinear   (1:size(overdensityLinearTmp   ))=overdensityLinearTmp
             overdensityNonLinear(1:size(overdensityNonLinearTmp))=overdensityNonLinearTmp
          end if
          overdensityLinear   (i)=+normalization            &
               &                  *0.6d0                    &
               &                  *(                        &
               &                    +1.0d0                  &
               &                    -OmegaMatterEpochal     &
               &                    -OmegaDarkEnergyEpochal &
               &                    -epsilonPerturbation    &
               &                   )                        &
               &                  /OmegaMatterEpochal
          overdensityNonLinear(i)=+1.0d0                    &
               &                  /radiusNow**3             &
               &                  -1.0d0
          if (overdensityNonLinear(i) <= -0.99d0) exit
       end do
       ! Reverse the arrays such that we have overdensity increasing.
       overdensityLinearTmp   =Array_Reverse(overdensityLinear   (1:i))
       overdensityNonLinearTmp=Array_Reverse(overdensityNonLinear(1:i))
       deallocate(overdensityLinear   )
       deallocate(overdensityNonLinear)
       call move_alloc(overdensityLinearTmp   ,overdensityLinear   )
       call move_alloc(overdensityNonLinearTmp,overdensityNonLinear)
       ! Populate the table.
       do iOverdensity=1,overdensityLinearCount
          ! Test for out of range overdensity.
          if      (overdensitiesLinear(iOverdensity) < overdensityLinear(1)) then
             ! Tabulated overdensity is lower than any we've computed. Use the lowest nonlinear overdensity.
             overdensityNonLinear_=overdensityNonLinear(1)
          else if (overdensitiesLinear(iOverdensity) > overdensityLinear(i)) then
             ! Tabulated overdensity exceeds any we've computed, so this overdensity is already collapsed. Use highest nonlinear
             ! overdensity.
             overdensityNonLinear_=overdensityNonLinear(i)
          else
             ! Find the tabulated in those computed and interpolate.
             iOverdensityLinear   =int(searchArray(overdensityLinear,overdensitiesLinear(iOverdensity)))
             overdensityNonLinear_=+  overdensityNonLinear(iOverdensityLinear  ) &
                  &                +(                                            &
                  &                  +overdensityNonLinear(iOverdensityLinear+1) &
                  &                  -overdensityNonLinear(iOverdensityLinear  ) &
                  &                 )                                            &
                  &                *(                                            &
                  &                  +overdensitiesLinear (iOverdensity        ) &
                  &                  -overdensityLinear   (iOverdensityLinear  ) &
                  &                 )                                            &
                  &                /(                                            &
                  &                  +overdensityLinear   (iOverdensityLinear+1) &
                  &                  -overdensityLinear   (iOverdensityLinear  ) &
                  &                 )
          end if
          ! Populate this point in the table.
          call linearNonlinearMap_%populate(overdensityNonLinear_,iOverdensity,iTime)
       end do
    end do
    return
  end subroutine cllsnlssMttCsmlgclCnstntLinearNonlinearMap

  double precision function cllsnlssMttCsmlgclCnstntRadiusRoot(radiusNow)
    !!{
    Root function used in solving for the radius of a perturbation.
    !!}
    use :: Numerical_Integration, only : integrator
    implicit none
    double precision            , intent(in   ) :: radiusNow
    double precision            , parameter     :: numericalLimitEpsilon=1.0d-4
    type            (integrator)                :: integrator_
    double precision                            :: radiusUpperLimit

    radiusUpperLimit                  =min(                                              &
         &                                 +(1.0d0-numericalLimitEpsilon)*radiusMaximum, &
         &                                 +                              radiusNow      &
         &                                )
    integrator_                       = integrator(cllsnlssMttCsmlgclCnstntPerturbationIntegrand,toleranceRelative=1.0d-6,hasSingularities =.true.)
    cllsnlssMttCsmlgclCnstntRadiusRoot=+integrator_%integrate(0.0d0,radiusUpperLimit) &
         &                             /hubbleTimeEpochal                             &
         &                             -timeTarget
    if (radiusUpperLimit < radiusNow) then
       cllsnlssMttCsmlgclCnstntRadiusRoot       =+cllsnlssMttCsmlgclCnstntRadiusRoot                                                                             &
            &                        -2.0d0*sqrt(+OmegaMatterEpochal/radiusUpperLimit   +      OmegaDarkEnergyEpochal*radiusUpperLimit**2+amplitudePerturbation) &
            &                        /          (-OmegaMatterEpochal/radiusUpperLimit**2+2.0d0*OmegaDarkEnergyEpochal*radiusUpperLimit                         ) &
            &                        /hubbleTimeEpochal
       if (radiusNow < radiusMaximum)                                                                                                                            &
            & cllsnlssMttCsmlgclCnstntRadiusRoot=+cllsnlssMttCsmlgclCnstntRadiusRoot                                                                             &
            &                        + 2.0d0*sqrt(+OmegaMatterEpochal/radiusNow         +      OmegaDarkEnergyEpochal*radiusNow       **2+amplitudePerturbation) &
            &                        /(           -OmegaMatterEpochal/radiusNow      **2+2.0d0*OmegaDarkEnergyEpochal*radiusNow                                ) &
            &                        /hubbleTimeEpochal
    end if
    return
  end function cllsnlssMttCsmlgclCnstntRadiusRoot

  subroutine cllsnlssMttCsmlgclCnstntRestoreTable(self,time,restoredTable,fileName,tableStore,status)
    !!{
    Attempt to restore a table from file.
    !!}
    use :: Error             , only : errorStatusFail, errorStatusSuccess
    use :: File_Utilities    , only : File_Exists
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char           , varying_string
    use :: Tables            , only : table1D        , table1DLogarithmicLinear
    implicit none
    class           (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt)             , intent(inout) :: self
    double precision                                                               , intent(in   ) :: time
    class           (table1D                                         ), allocatable, intent(inout) :: restoredTable
    type            (varying_string                                  )             , intent(in   ) :: fileName
    logical                                                                        , intent(in   ) :: tableStore
    integer                                                                        , intent(  out) :: status
    type            (hdf5Object                                      )                             :: file
    double precision                                                  , allocatable, dimension(:)  :: timeTable    , valueTable
    !$GLC attributes unused :: self

    status=errorStatusFail
    if (.not.tableStore) return
    if (File_Exists(fileName)) then
       !$ call hdf5Access%set()
       file=hdf5Object(char(fileName))
       call file%readDataset('time',timeTable)
       if     (                                    &
            &   timeTable(1              ) <= time &
            &  .and.                               &
            &   timeTable(size(timeTable)) >= time &
            & ) then
          call file%readDataset('value',valueTable)
          ! Deallocate table if currently allocated.
          if (allocated(restoredTable)) then
             call restoredTable%destroy()
             deallocate(restoredTable)
          end if
          allocate(table1DLogarithmicLinear :: restoredTable)
          select type (restoredTable)
          type is (table1DLogarithmicLinear)
             call restoredTable%create  (timeTable (1),timeTable(size(timeTable)),size(timeTable))
             call restoredTable%populate(valueTable                                              )
          end select
          status=errorStatusSuccess
       end if
       !$ call hdf5Access%unset()
    end if
    return
  end subroutine cllsnlssMttCsmlgclCnstntRestoreTable

  subroutine cllsnlssMttCsmlgclCnstntStoreTable(self,storeTable,fileName,tableStore)
    !!{
    Store a table to file.
    !!}
    use :: File_Utilities    , only : Directory_Make, File_Path
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char          , varying_string
    use :: Tables            , only : table1D
    implicit none
    class  (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt), intent(inout) :: self
    class  (table1D                                         ), intent(in   ) :: storeTable
    type   (varying_string                                  ), intent(in   ) :: fileName
    logical                                                  , intent(in   ) :: tableStore
    type   (hdf5Object                                      )                :: file
    !$GLC attributes unused :: self

    if (.not.tableStore) return
    call Directory_Make(char(File_Path(char(fileName))))
    !$ call hdf5Access%set()
    file=hdf5Object(char(fileName),overWrite=.true.,readOnly=.false.)
    call file%writeDataset(        storeTable%xs()                     ,'time' )
    call file%writeDataset(reshape(storeTable%ys(),[storeTable%size()]),'value')
    !$ call hdf5Access%unset()
    return
  end subroutine cllsnlssMttCsmlgclCnstntStoreTable
