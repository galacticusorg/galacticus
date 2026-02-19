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
  A spherical collapse solver class for universes consisting of collisionless matter and dark energy.
  !!}

  ! Enumeration of radii at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.
  !![
  <enumeration>
   <name>cllsnlssMttrDarkEnergyFixedAt</name>
   <description>Enumeration of radii at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="undefined"    />
   <entry label="turnaround"   />
   <entry label="virialization"/>
  </enumeration>
  !!]

  !![
  <sphericalCollapseSolver name="sphericalCollapseSolverCllsnlssMttrDarkEnergy">
   <description>A spherical collapse solver for universes consisting of collisionless matter and dark energy.</description>
  </sphericalCollapseSolver>
  !!]
  type, extends(sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt) :: sphericalCollapseSolverCllsnlssMttrDarkEnergy
     !!{
     A spherical collapse solver for universes consisting of collisionless matter and dark energy.
     !!}
     private
     type(enumerationCllsnlssMttrDarkEnergyFixedAtType) :: energyFixedAt
   contains
     procedure :: linearNonlinearMap => cllsnlssMttrDarkEnergyLinearNonlinearMap
     procedure :: tabulate           => cllsnlssMttrDarkEnergyTabulate
  end type sphericalCollapseSolverCllsnlssMttrDarkEnergy

  interface sphericalCollapseSolverCllsnlssMttrDarkEnergy
     !!{
     Constructors for the \refClass{sphericalCollapseSolverCllsnlssMttrDarkEnergy} spherical collapse solver class.
     !!}
     module procedure cllsnlssMttrDarkEnergyConstructorParameters
     module procedure cllsnlssMttrDarkEnergyConstructorInternal
  end interface sphericalCollapseSolverCllsnlssMttrDarkEnergy

  ! Pointer to the default cosmology functions object.
  class           (cosmologyFunctionsClass), pointer   :: cosmologyFunctions_            => null()
  !$omp threadprivate(cosmologyFunctions_)

  ! Fraction of current expansion factor to use as initial time in perturbation dynamics solver.
  double precision                         , parameter :: expansionFactorInitialFraction =  1.0d-6

  ! Variables used in root finding.
  double precision                                     :: perturbationRadiusInitial
  !$omp threadprivate(perturbationRadiusInitial)

contains

  function cllsnlssMttrDarkEnergyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{sphericalCollapseSolverCllsnlssMttrDarkEnergy} spherical collapse solver class that takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (sphericalCollapseSolverCllsnlssMttrDarkEnergy)                :: self
    type   (inputParameters                              ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                      ), pointer       :: cosmologyFunctions_
    class  (linearGrowthClass                            ), pointer       :: linearGrowth_
    type   (varying_string                               )                :: energyFixedAt

    !![
    <inputParameter>
      <name>energyFixedAt</name>
      <source>parameters</source>
      <defaultValue>var_str('turnaround')</defaultValue>
      <description>The radius at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="linearGrowth"       name="linearGrowth_"       source="parameters"/>
    !!]
    self=sphericalCollapseSolverCllsnlssMttrDarkEnergy(enumerationCllsnlssMttrDarkEnergyFixedAtEncode(char(energyFixedAt),includesPrefix=.false.),cosmologyFunctions_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    <objectDestructor name="linearGrowth_"      />
    !!]
    return
  end function cllsnlssMttrDarkEnergyConstructorParameters

  function cllsnlssMttrDarkEnergyConstructorInternal(energyFixedAt,cosmologyFunctions_,linearGrowth_) result(self)
    !!{
    Internal constructor for the \refClass{sphericalCollapseSolverCllsnlssMttrDarkEnergy} spherical collapse solver class.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath   , pathTypeDataDynamic
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    type (sphericalCollapseSolverCllsnlssMttrDarkEnergy)                                  :: self
    type (enumerationCllsnlssMttrDarkEnergyFixedAtType ), intent(in   )                   :: energyFixedAt
    class(cosmologyFunctionsClass                      ), intent(in   ), target           :: cosmologyFunctions_
    class(linearGrowthClass                            ), intent(in   ), target, optional :: linearGrowth_
    !![
    <constructorAssign variables="energyFixedAt, *cosmologyFunctions_, *linearGrowth_"/>
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
    if (.not.enumerationCllsnlssMttrDarkEnergyFixedAtIsValid(energyFixedAt)) call Error_Report('invalid energyFixedAt'//{introspection:location})
    return
  end function cllsnlssMttrDarkEnergyConstructorInternal

  subroutine cllsnlssMttrDarkEnergyTabulate(self,time,sphericalCollapse_,calculationType)
    !!{
    Tabulate spherical collapse solutions for $\delta_\mathrm{crit}$, $\Delta_\mathrm{vir}$, or $R_\mathrm{ta}/R_\mathrm{vir}$ vs. time.
    !!}
    use :: Display      , only : displayCounter           , displayCounterClear          , displayIndent                , displayUnindent, &
          &                      verbosityLevelWorking
    use :: Error        , only : Error_Report
    use :: Linear_Growth, only : normalizeMatterDominated
    use :: Root_Finder  , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: Tables       , only : table1DLogarithmicLinear
    implicit none
    class           (sphericalCollapseSolverCllsnlssMttrDarkEnergy)             , intent(inout) :: self
    double precision                                                            , intent(in   ) :: time
    type            (enumerationCllsnlssMttCsmlgclCnstntClcltnType)             , intent(in   ) :: calculationType
    class           (table1D                                      ), allocatable, intent(inout) :: sphericalCollapse_
    class           (linearGrowthClass                            ), pointer                    :: linearGrowth_                  => null()
    double precision                                               , parameter                  :: toleranceAbsolute              =  0.0d0  , toleranceRelative              =1.0d-9
    double precision                                               , dimension(2)               :: timeRange
    type            (rootFinder                                   ), save                       :: finderAmplitudePerturbation              , finderExpansionMaximum
    logical                                                                                     :: finderAmplitudeConstructed     =  .false., finderExpansionConstructed     =.false.
    !$omp threadprivate(finderAmplitudePerturbation,finderExpansionMaximum,finderAmplitudeConstructed,finderExpansionConstructed)
    integer                                                                                     :: countTimes                               , iTime                                  , &
         &                                                                                         iCount
    double precision                                                                            :: expansionFactor                          , epsilonPerturbation                    , &
         &                                                                                         epsilonPerturbationMaximum               , epsilonPerturbationMinimum             , &
         &                                                                                         densityContrastExpansionMaximum          , expansionFactorExpansionMaximum        , &
         &                                                                                         radiusExpansionMaximum                   , maximumExpansionTime                   , &
         &                                                                                         normalization                            , q                                      , &
         &                                                                                         timeEnergyFixed                          , timeInitial                            , &
         &                                                                                         y                                        , timeMinimum                            , &
         &                                                                                         timeMaximum
    double complex                                                                              :: a                                        , b                                      , &
         &                                                                                         x
    type            (varying_string                               )                             :: message
    character       (len=13                                       )                             :: label

    ! Validate.
    if (calculationType == cllsnlssMttCsmlgclCnstntClcltnCriticalOverdensity .and. .not.associated(self%linearGrowth_)) call Error_Report('linearGrowth object was not provided'//{introspection:location})
    ! Find minimum and maximum times to tabulate.
    if (allocated(sphericalCollapse_)) then
       ! Use currently tabulated range as the starting point.
       timeMinimum=sphericalCollapse_%x(+1)
       timeMaximum=sphericalCollapse_%x(-1)
    else
       ! Specify an initial default range.
       timeMinimum= 0.1d0
       timeMaximum=20.0d0
    end if
    ! Expand the range to ensure the requested time is included.
    timeMinimum=min(timeMinimum,time/2.0d0)
    timeMaximum=max(timeMaximum,time*2.0d0)
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
       ! Create the table.
       call sphericalCollapse_%create(timeMinimum,timeMaximum,countTimes)
       ! Solve ODE to get corresponding overdensities.
       message="Solving spherical collapse model for dark energy universe for "
       write (label,'(e12.6)') timeMinimum
       message=message//trim(adjustl(label))//" ≤ t/Gyr ≤ "
       write (label,'(e12.6)') timeMaximum
       message=message//trim(adjustl(label))
       call displayIndent(message,verbosity=verbosityLevelWorking)
       iCount=0
       call displayCounter(                                 &
            &                        iCount               , &
            &              isNew    =.true.               , &
            &              verbosity=verbosityLevelWorking  &
            &             )
       !$omp parallel private(expansionFactor,epsilonPerturbationMaximum,epsilonPerturbationMinimum,epsilonPerturbation,timeInitial,timeRange,maximumExpansionTime,expansionFactorExpansionMaximum,q,y,timeEnergyFixed,a,b,x,linearGrowth_)       
       !$omp critical(sphericalCollapseSolveCllnlssMttrDrkEnrgyDeepCopy)
       allocate(cosmologyFunctions_,mold=self%cosmologyFunctions_)
       !![
       <deepCopyReset variables="self%cosmologyFunctions_"/>
       <deepCopy source="self%cosmologyFunctions_" destination="cosmologyFunctions_"/>
       <deepCopyFinalize variables="cosmologyFunctions_"/>
       !!]
       if (calculationType == cllsnlssMttCsmlgclCnstntClcltnCriticalOverdensity) then
          allocate(linearGrowth_,mold=self%linearGrowth_)
          !![
          <deepCopyReset variables="self%linearGrowth_"/>
          <deepCopy source="self%linearGrowth_" destination="linearGrowth_"/>
          <deepCopyFinalize variables="linearGrowth_"/>
          !!]
       else
          linearGrowth_ => null()
       end if
       !$omp end critical(sphericalCollapseSolveCllnlssMttrDrkEnrgyDeepCopy)
       !$omp do schedule(dynamic)
       do iTime=1,countTimes
          call displayCounter(                                                        &
               &                        int(100.0d0*dble(iCount-1)/dble(countTimes)), &
               &              isNew    =.false.                                     , &
               &              verbosity=verbosityLevelWorking                         &
               &             )
          ! Get the current expansion factor.
          expansionFactor=cosmologyFunctions_%expansionFactor(sphericalCollapse_%x(iTime))
          ! In the case of dark energy we cannot (easily) determine the largest (i.e. least negative) value of ε for which a
          ! perturbation can collapse. So, use no perturbation.
          epsilonPerturbationMaximum=  0.0d0
          ! Estimate a suitably negative minimum value for ε.
          epsilonPerturbationMinimum=-10.0d0
          ! Evaluate cosmological parameters at the present time.
          OmegaMatterEpochal    =cosmologyFunctions_%omegaMatterEpochal    (expansionFactor=expansionFactor)
          OmegaDarkEnergyEpochal=cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor)
          hubbleTimeEpochal     =cosmologyFunctions_%expansionRate         (                expansionFactor)
          time_                  =sphericalCollapse_%x                     (                iTime          )
          ! Check dark energy equation of state is within acceptable range.
          if (cosmologyFunctions_%equationOfStateDarkEnergy(time=time_) >= -1.0d0/3.0d0) &
               & call Error_Report('ω<-⅓ required'//{introspection:location})
          ! Find the value of ε for which the perturbation just collapses at this time.
          if (.not.finderAmplitudeConstructed) then
             finderAmplitudePerturbation=rootFinder(                                                       &
                  &                                 rootFunction=cllsnlssMttrDarkEnergyRadiusPerturbation, &
                  &                                 toleranceAbsolute=toleranceAbsolute                  , &
                  &                                 toleranceRelative=toleranceRelative                    &
                  &                                )
             finderAmplitudeConstructed=.true.
          end if
          epsilonPerturbation=finderAmplitudePerturbation%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
          select case (calculationType%ID)
          case (cllsnlssMttCsmlgclCnstntClcltnCriticalOverdensity%ID)
             ! Critical linear overdensity.
             normalization=+linearGrowth_%value(time_,normalize=normalizeMatterDominated) &
                  &        /                    expansionFactor
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
                  &                           iTime                      &
                  &                          )
          case (cllsnlssMttCsmlgclCnstntClcltnVirialDensityContrast%ID,cllsnlssMttCsmlgclCnstntClcltnRadiusTurnaround%ID)
             ! Find the epoch of maximum expansion for the perturbation.
             if (.not.finderExpansionConstructed) then
                finderExpansionMaximum=rootFinder(                                                                   &
                     &                            rootFunction     =cllsnlssMttrDarkEnergyExpansionRatePerturbation, &
                     &                            toleranceAbsolute=toleranceAbsolute                              , &
                     &                            toleranceRelative=toleranceRelative                                &
                     &                           )
                finderExpansionConstructed=.true.
             end if
             call finderExpansionMaximum%rangeExpand (                                                             &
                  &                                   rangeExpandDownward          =1.0d0-1.0d-2                 , &
                  &                                   rangeExpandUpward            =1.0d0+1.0d-2                 , &
                  &                                   rangeExpandType              =rangeExpandMultiplicative    , &
                  &                                   rangeUpwardLimit             =time_                        , &
                  &                                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
                  &                                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative  &
                  &                                  )
             amplitudePerturbation=epsilonPerturbation
             ! Compute the corresponding time of maximum expansion.
             timeInitial                    =cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactor(time_)*expansionFactorInitialFraction)
             ! Guess that the time of maximum expansion occurred at close to half of the current time.
             timeRange                      =[0.45d0,0.55d0]*time_
             maximumExpansionTime           =finderExpansionMaximum%find(rootRange=timeRange)
             expansionFactorExpansionMaximum=cosmologyFunctions_%expansionFactor(maximumExpansionTime)
             ! Solve the dynamics of the perturbation to find the radius at the point of maximum expansion.
             call cllsnlssMttrDarkEnergyPerturbationDynamicsSolver(epsilonPerturbation,maximumExpansionTime,radiusExpansionMaximum)
             ! Compute the density contrast of the perturbation at maximum expansion.
             densityContrastExpansionMaximum=(expansionFactorExpansionMaximum/expansionFactor/radiusExpansionMaximum)**3
             ! Solve the cubic equation (Percival, 2005, A&A, 443, 819, eqn. 38) to give the ratio of virial to turnaround radii,
             ! x.
             q=      cosmologyFunctions_%omegaDarkEnergyEpochal(time=maximumExpansionTime) &
                  & /cosmologyFunctions_%omegaMatterEpochal    (time=maximumExpansionTime) &
                  & /densityContrastExpansionMaximum
             y=      expansionFactorExpansionMaximum**cosmologyFunctions_%exponentDarkEnergy(time=maximumExpansionTime) &
                  & /expansionFactor                **cosmologyFunctions_%exponentDarkEnergy(time=time_               )
             select case (self%energyFixedAt%ID)
             case (cllsnlssMttrDarkEnergyFixedAtTurnaround   %ID)
                timeEnergyFixed=maximumExpansionTime
             case (cllsnlssMttrDarkEnergyFixedAtVirialization%ID)
                timeEnergyFixed=time_
             case default
                call Error_Report('unrecognized epoch'//{introspection:location})
             end select
             a=1.0d0-(1.0d0+3.0d0*cosmologyFunctions_%equationOfStateDarkEnergy(time=timeEnergyFixed))*q/2.0d0
             b=      (1.0d0+3.0d0*cosmologyFunctions_%equationOfStateDarkEnergy(time=time_          ))*q/y
             ! Check for special cases.
             if (q == 0.0d0) then
                ! No dark energy, the ratio of radii is always ½.
                x=0.5d0
             else
                ! General case, use the appropriate root of the cubic equation.
                x=                                                                                                        &
                     & +(0.0d0,0.5d0)*sqrt(3.0d0)                                                                         &
                     & *(                                                                                                 &
                     &   1.0d0/b*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(+1.0d0/3.0d0)/ 6.0d0  &
                     &  +2.0d0*a*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(-1.0d0/3.0d0)         &
                     & )                                                                                                  &
                     & -(1.0d0/b*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(+1.0d0/3.0d0)/12.0d0) &
                     & +(      a*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(-1.0d0/3.0d0)       )
             end if
             select case (calculationType%ID)
             case (cllsnlssMttCsmlgclCnstntClcltnVirialDensityContrast%ID)
                call sphericalCollapse_%populate(                                           &
                     &                           1.0d0/(dble(x)*radiusExpansionMaximum)**3, &
                     &                           iTime                                      &
                     &                          )
             case (cllsnlssMttCsmlgclCnstntClcltnRadiusTurnaround     %ID)
                call sphericalCollapse_%populate(                                           &
                     &                           1.0d0/ dble(x)                           , &
                     &                           iTime                                      &
                     &                          )
             end select
          end select
          !$omp atomic
          iCount=iCount+1
       end do
       !$omp end do
       !![
       <objectDestructor name="cosmologyFunctions_"/>
       !!]
       if (calculationType == cllsnlssMttCsmlgclCnstntClcltnCriticalOverdensity) then
          !![
	  <objectDestructor name="linearGrowth_"/>
          !!]
       end if
       !$omp end parallel
       call displayCounterClear(       verbosity=verbosityLevelWorking)
       call displayUnindent    ('done',verbosity=verbosityLevelWorking)
    end select
    return
  end subroutine cllsnlssMttrDarkEnergyTabulate

  double precision function cllsnlssMttrDarkEnergyRadiusPerturbation(epsilonPerturbation)
    !!{
    Return the radius of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    amplitude {\normalfont \ttfamily epsilonPerturbation}.
    !!}
    implicit none
    double precision, intent(in   ) :: epsilonPerturbation

    call cllsnlssMttrDarkEnergyPerturbationDynamicsSolver(epsilonPerturbation,time_,cllsnlssMttrDarkEnergyRadiusPerturbation)
    return
  end function cllsnlssMttrDarkEnergyRadiusPerturbation

  double precision function cllsnlssMttrDarkEnergyExpansionRatePerturbation(time)
    !!{
    Return the expansion rate of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    amplitude {\normalfont \ttfamily epsilonPerturbation}.
    !!}
    implicit none
    double precision, intent(in   ) :: time

    call cllsnlssMttrDarkEnergyPerturbationDynamicsSolver(amplitudePerturbation,time,expansionRatePerturbation=cllsnlssMttrDarkEnergyExpansionRatePerturbation)
    return
  end function cllsnlssMttrDarkEnergyExpansionRatePerturbation

  subroutine cllsnlssMttrDarkEnergyPerturbationDynamicsSolver(epsilonPerturbation,time,radiusPerturbation,expansionRatePerturbation)
    !!{
    Integrate the dynamics of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    amplitude {\normalfont \ttfamily epsilonPerturbation}.
    !!}
    use :: Interface_GSL        , only : GSL_Success
    use :: Numerical_ODE_Solvers, only : odeSolver
    implicit none
    double precision                                       , intent(in   )           :: epsilonPerturbation                 , time
    double precision                                       , intent(  out), optional :: expansionRatePerturbation           , radiusPerturbation
    integer         (c_size_t ), parameter                                           :: countProperties          =2_c_size_t
    double precision           , dimension(countProperties)                          :: propertyValues
    double precision           , parameter                                           :: odeToleranceAbsolute     =0.0d0     , odeToleranceRelative            =1.0d-12
    type            (odeSolver)                                                      :: solver
    double precision                                                                 :: expansionFactorInitial               , expansionRatePerturbationInitial        , &
         &                                                                              overdensityInitial                   , timeInitial
    integer                                                                          :: odeStatus

    ! Specify a sufficiently early time.
    expansionFactorInitial=expansionFactorInitialFraction
    ! Find the corresponding cosmic time.
    timeInitial=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactor(time)*expansionFactorInitial)
    ! Find the overdensity of the perturbation at early time (Percival, 2005, A&A, 443, 819, eqn. 25).
    overdensityInitial=+0.6d0                    &
         &             *(                        &
         &               +1.0d0                  &
         &               -OmegaMatterEpochal     &
         &               -OmegaDarkEnergyEpochal &
         &               -epsilonPerturbation    &
         &              )                        &
         &             *(                        &
         &               +3.0d0                  &
         &               /2.0d0                  &
         &               /OmegaMatterEpochal     &
         &               *hubbleTimeEpochal      &
         &               *timeInitial            &
         &              )**(2.0d0/3.0d0)
    ! Find the perturbation radius at this early time. This is, by construction, just the initial expansion factor.
    perturbationRadiusInitial=+expansionFactorInitial
    ! Find the perturbation expansion rate at early time (Percival, 2005, A&A, 443, 819, eqn. 22).
    expansionRatePerturbationInitial=+hubbleTimeEpochal            &
         &                           *sqrt(                        &
         &                                 +OmegaMatterEpochal     &
         &                                 /expansionFactorInitial &
         &                                 +epsilonPerturbation    &
         &                                )
    ! Set initial conditions.
    propertyValues=[                                  &
         &          perturbationRadiusInitial       , &
         &          expansionRatePerturbationInitial  &
         &         ]
    ! Evolve if the requested time is after the initial time.
    if (time > timeInitial) then
       ! Solve the ODE to find the perturbation radius at the present day.
       solver=odeSolver(countProperties,cllsnlssMttrDarkEnergyPerturbationODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
       call solver%solve(timeInitial,time,propertyValues,status=odeStatus)
       ! If the ODE solver did not succeed, it is because the perturbation collapsed to zero radius (causing a divergence). This
       ! means it collapsed prior to the current time. We extrapolate to negative radius (using the velocity at the final step) to
       ! permit our root finder to locate the point at which collapse occurs at the current time.
       if (odeStatus /= GSL_Success) propertyValues(1)=+propertyValues(1) &
            &                                          +propertyValues(2) &
            &                                          *(                 &
            &                                            +time            &
            &                                            -timeInitial     &
            &                                           )
    end if
    ! Return the computed quantities.
    if (present(radiusPerturbation       )) radiusPerturbation       =propertyValues(1)
    if (present(expansionRatePerturbation)) expansionRatePerturbation=propertyValues(2)
    return
  end subroutine cllsnlssMttrDarkEnergyPerturbationDynamicsSolver

  integer function cllsnlssMttrDarkEnergyPerturbationODEs(time,y,dydt)
    !!{
    Differential equations describing the evolution of spherical perturbations in a universe containing collisionless dark matter and dark energy.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision, intent(in   )               :: time
    double precision, intent(in   ), dimension(:) :: y
    double precision, intent(  out), dimension(:) :: dydt
    double precision                              :: expansionFactor

    if (y(1) <= 0.0d0) then
       dydt(1:2)=0.0d0
    else
       expansionFactor=+cosmologyFunctions_%expansionFactor(time) &
            &          /cosmologyFunctions_%expansionFactor(time_)
       dydt(1)=+y(2)
       dydt(2)=-0.5d0                                                                                                                                                              &
            &  *y(1)                                                                                                                                                               &
            &  *hubbleTimeEpochal**2                                                                                                                                               &
            &  *(                                                                                                                                                                  &
            &    +                                                                       OmegaMatterEpochal    /y(1)           **3                                                 &
            &    +(3.0d0*cosmologyFunctions_%equationOfStateDarkEnergy(time=time)+1.0d0)*OmegaDarkEnergyEpochal*expansionFactor**cosmologyFunctions_%exponentDarkEnergy(time=time) &
            &   )
    end if
    ! Return success.
    cllsnlssMttrDarkEnergyPerturbationODEs=GSL_Success
    return
  end function cllsnlssMttrDarkEnergyPerturbationODEs

  subroutine cllsnlssMttrDarkEnergyLinearNonlinearMap(self,time,linearNonlinearMap_)
    !!{
    Tabulate the mapping between linear and nonlinear overdensity for the spherical collapse model.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (sphericalCollapseSolverCllsnlssMttrDarkEnergy), intent(inout) :: self
    double precision                                               , intent(in   ) :: time
    class           (table2DLinLinLin                             ), intent(inout) :: linearNonlinearMap_
    !$GLC attributes unused :: self, time, linearNonlinearMap_

    call Error_Report('linear-nonlinear mapping is not supported by this class'//{introspection:location})
    return
  end subroutine cllsnlssMttrDarkEnergyLinearNonlinearMap
