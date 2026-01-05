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

  !!{
  A spherical collapse solver class for universes consisting of baryons, collisionless matter, and dark energy.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass

  !![
  <sphericalCollapseSolver name="sphericalCollapseSolverBaryonsDarkMatterDarkEnergy">
   <description>A spherical collapse solver for universes consisting of baryons, collisionless matter, and dark energy.</description>
  </sphericalCollapseSolver>
  !!]
  type, extends(sphericalCollapseSolverCllsnlssMttrDarkEnergy) :: sphericalCollapseSolverBaryonsDarkMatterDarkEnergy
     !!{
     A spherical collapse solver for universes consisting of baryons, collisionless matter, and dark energy.
     !!}
     private
     class  (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     logical                                    :: baryonsCluster
     integer                                    :: tablePointsPerOctave
   contains
     final     ::             baryonsDarkMatterDarkEnergyDestructor
     procedure :: tabulate => baryonsDarkMatterDarkEnergyTabulate
  end type sphericalCollapseSolverBaryonsDarkMatterDarkEnergy

  interface sphericalCollapseSolverBaryonsDarkMatterDarkEnergy
     !!{
     Constructors for the \refClass{sphericalCollapseSolverBaryonsDarkMatterDarkEnergy} spherical collapse solver class.
     !!}
     module procedure baryonsDarkMatterDarkEnergyConstructorParameters
     module procedure baryonsDarkMatterDarkEnergyConstructorInternal
  end interface sphericalCollapseSolverBaryonsDarkMatterDarkEnergy

  ! Variables used in root finding.
  double precision :: OmegaBaryonEpochal
  logical          :: baryonsCluster
  !$omp threadprivate(OmegaBaryonEpochal,baryonsCluster)

contains

  function baryonsDarkMatterDarkEnergyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{sphericalCollapseSolverBaryonsDarkMatterDarkEnergy} spherical collapse solver class that takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (sphericalCollapseSolverBaryonsDarkMatterDarkEnergy)                :: self
    type   (inputParameters                                   ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                           ), pointer       :: cosmologyFunctions_
    class  (cosmologyParametersClass                          ), pointer       :: cosmologyParameters_
    type   (varying_string                                    )                :: energyFixedAt
    logical                                                                    :: baryonsCluster
    integer                                                                    :: tablePointsPerOctave

    !![
    <inputParameter>
      <name>baryonsCluster</name>
      <source>parameters</source>
      <description>If true baryons are assumed to cluster in the same way as collisionless matter. If false, baryons are assumed not to cluster at all.</description>
    </inputParameter>
    <inputParameter>
      <name>energyFixedAt</name>
      <source>parameters</source>
      <defaultValue>var_str('turnaround')</defaultValue>
      <description>The radius at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.</description>
    </inputParameter>
    <inputParameter>
      <name>tablePointsPerOctave</name>
      <source>parameters</source>
      <defaultValue>300</defaultValue>
      <description>The number of points per octave of time at which to tabulate solutions.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=sphericalCollapseSolverBaryonsDarkMatterDarkEnergy(baryonsCluster,tablePointsPerOctave,enumerationCllsnlssMttrDarkEnergyFixedAtEncode(char(energyFixedAt),includesPrefix=.false.),cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function baryonsDarkMatterDarkEnergyConstructorParameters

  function baryonsDarkMatterDarkEnergyConstructorInternal(baryonsCluster,tablePointsPerOctave,energyFixedAt,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{sphericalCollapseSolverBaryonsDarkMatterDarkEnergy} spherical collapse solver class.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath                      , pathTypeDataDynamic
    use :: ISO_Varying_String, only : operator(//)
    use :: Linear_Growth     , only : linearGrowthCollisionlessMatter, linearGrowthNonClusteringBaryonsDarkMatter
    implicit none
    type   (sphericalCollapseSolverBaryonsDarkMatterDarkEnergy)                        :: self
    logical                                                    , intent(in   )         :: baryonsCluster
    integer                                                    , intent(in   )         :: tablePointsPerOctave
    type   (enumerationCllsnlssMttrDarkEnergyFixedAtType      ), intent(in   )         :: energyFixedAt
    class  (cosmologyFunctionsClass                           ), intent(in   ), target :: cosmologyFunctions_
    class  (cosmologyParametersClass                          ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="baryonsCluster, tablePointsPerOctave, energyFixedAt, *cosmologyFunctions_, *cosmologyParameters_"/>
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
    if (baryonsCluster) then
       allocate(linearGrowthCollisionlessMatter            :: self%linearGrowth_)
    else
       allocate(linearGrowthNonClusteringBaryonsDarkMatter :: self%linearGrowth_)
    end if
    select type (linearGrowth_ => self%linearGrowth_)
    type is (linearGrowthCollisionlessMatter           )
       !![
       <referenceConstruct object="linearGrowth_" constructor="linearGrowthCollisionlessMatter           (self%cosmologyParameters_,self%cosmologyFunctions_)"/>
       !!]
    type is (linearGrowthNonClusteringBaryonsDarkMatter)
       !![
       <referenceConstruct object="linearGrowth_" constructor="linearGrowthNonClusteringBaryonsDarkMatter(self%cosmologyParameters_,self%cosmologyFunctions_)"/>
       !!]
    end select
    return
  end function baryonsDarkMatterDarkEnergyConstructorInternal

  subroutine baryonsDarkMatterDarkEnergyDestructor(self)
    !!{
    Destructor for the \refClass{sphericalCollapseSolverBaryonsDarkMatterDarkEnergy} spherical collapse solver class.
    !!}
    implicit none
    type(sphericalCollapseSolverBaryonsDarkMatterDarkEnergy), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine baryonsDarkMatterDarkEnergyDestructor

  subroutine baryonsDarkMatterDarkEnergyTabulate(self,time,sphericalCollapse_,calculationType)
    !!{
    Tabulate spherical collapse solutions for $\delta_\mathrm{crit}$, $\Delta_\mathrm{vir}$, or $R_\mathrm{ta}/R_\mathrm{vir}$ vs. time.
    !!}
    use :: Display    , only : displayCounter           , displayCounterClear          , displayIndent                , displayUnindent, &
          &                    verbosityLevelWorking
    use :: Error      , only : Error_Report
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    use :: Tables     , only : table1DLogarithmicLinear
    implicit none
    class           (sphericalCollapseSolverBaryonsDarkMatterDarkEnergy)             , intent(inout)  :: self
    double precision                                                                 , intent(in   )  :: time
    type            (enumerationCllsnlssMttCsmlgclCnstntClcltnType     )             , intent(in   )  :: calculationType
    class           (table1D                                           ), allocatable, intent(inout)  :: sphericalCollapse_
    class           (linearGrowthClass                                 ), pointer                     :: linearGrowth_                  => null()
    double precision                                                    , parameter                   :: toleranceAbsolute              =  0.0d0  , toleranceRelative              =1.0d-12
    double precision                                                                 , dimension(2  ) :: timeRange
    double precision                                                    , allocatable, dimension(:  ) :: timesPrevious
    double precision                                                    , allocatable, dimension(:,:) :: valuesPrevious
    type            (rootFinder                                        ), save                        :: finderPerturbationInitial                , finderExpansionMaximum
    logical                                                             , save                        :: finderPerturbationConstructed  =  .false., finderExpansionConstructed     =.false.
    !$omp threadprivate(finderPerturbationInitial,finderExpansionMaximum,finderPerturbationConstructed,finderExpansionConstructed)
    integer                                                                                           :: countTimes                               , iTime                                  , &
         &                                                                                               iCount                                   , iTimeMinimum                           , &
         &                                                                                               iTimeMaximum                             , countTimesEffective
    double precision                                                                                  :: expansionFactor                          , epsilonPerturbation                    , &
         &                                                                                               epsilonPerturbationMaximum               , epsilonPerturbationMinimum             , &
         &                                                                                               densityContrastExpansionMaximum          , expansionFactorExpansionMaximum        , &
         &                                                                                               radiusExpansionMaximum                   , timeExpansionMaximum                   , &
         &                                                                                               normalization                            , q                                      , &
         &                                                                                               timeEnergyFixed                          , timeInitial                            , &
         &                                                                                               y                                        , timeMinimum                            , &
         &                                                                                               timeMaximum                              , r                                      , &
         &                                                                                               z                                        , fractionDarkMatter
    double complex                                                                                    :: a                                        , b                                      , &
         &                                                                                               x
    type            (varying_string                                    )                              :: message
    character       (len=13                                            )                              :: label

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
    ! Round to the nearest factor of 2.
    timeMinimum=2.0d0**floor  (log(timeMinimum)/log(2.0d0))
    timeMaximum=2.0d0**ceiling(log(timeMaximum)/log(2.0d0))
    ! Determine number of points to tabulate.
    countTimes=nint(log(timeMaximum/timeMinimum)/log(2.0d0)*dble(self%tablePointsPerOctave))
    ! Copy baryon clustering option to module-scope.
    baryonsCluster=self%baryonsCluster
    ! Deallocate table if currently allocated.
    if (allocated(sphericalCollapse_)) then
       ! Store the current solution.
       timesPrevious      =sphericalCollapse_%xs()
       valuesPrevious     =sphericalCollapse_%ys()
       iTimeMinimum       =nint(log(timesPrevious(                 1 )/timeMinimum)/log(2.0d0)*self%tablePointsPerOctave)+1
       iTimeMaximum       =nint(log(timesPrevious(size(timesPrevious))/timeMinimum)/log(2.0d0)*self%tablePointsPerOctave)+1
       countTimesEffective=countTimes-(iTimeMaximum-iTimeMinimum+1)
       ! Destroy the table.
       call sphericalCollapse_%destroy()
       deallocate(sphericalCollapse_)
    else
       iTimeMinimum       =+huge(0)
       iTimeMaximum       =-huge(0)
       countTimesEffective=countTimes
    end if
    allocate(table1DLogarithmicLinear :: sphericalCollapse_)
    select type (sphericalCollapse_)
    type is (table1DLogarithmicLinear)
       ! Create the table.
       call sphericalCollapse_%create(timeMinimum,timeMaximum,countTimes)
       ! Solve ODE to get corresponding overdensities.
       message="Solving spherical collapse model for baryons + dark matter + dark energy universe for "
       select case (calculationType%ID)
       case (cllsnlssMttCsmlgclCnstntClcltnCriticalOverdensity  %ID)
          message=message//"critical overdensity"
       case (cllsnlssMttCsmlgclCnstntClcltnVirialDensityContrast%ID)
          message=message//"virial density contrast"
       case (cllsnlssMttCsmlgclCnstntClcltnRadiusTurnaround     %ID)
          message=message//"turnaround radius"
       end select
       if (baryonsCluster) then
          message=message//" (baryons do cluster)"
       else
          message=message//" (baryons do not cluster)"
       end if
       message=message//" for "
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
       !$omp parallel private(expansionFactor,epsilonPerturbationMaximum,epsilonPerturbationMinimum,epsilonPerturbation,timeInitial,timeRange,timeExpansionMaximum,expansionFactorExpansionMaximum,q,y,r,z,timeEnergyFixed,a,b,x,linearGrowth_)
       !$omp critical(sphericalCollapseSolverBrynsDrkMttrDrkEnrgyDeepCopy)
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
       !$omp end critical(sphericalCollapseSolverBrynsDrkMttrDrkEnrgyDeepCopy)
       !$omp do schedule(dynamic)
       do iTime=1,countTimes
          if (iTime >= iTimeMinimum .and. iTime <= iTimeMaximum) then
             call sphericalCollapse_%populate(                                        &
                  &                           valuesPrevious(iTime-iTimeMinimum+1,1), &
                  &                                          iTime                    &
                  &                          )
          else
             call displayCounter(                                                       &
                  &              int(100.0d0*dble(iCount-1)/dble(countTimesEffective)), &
                  &              isNew=.false.                                        , &
                  &              verbosity=verbosityLevelWorking                        &
                  &             )
             ! Get the current expansion factor.
             expansionFactor=cosmologyFunctions_%expansionFactor(sphericalCollapse_%x(iTime))
             ! Initial guess for the range of the initial perturbation amplitude. Since we expect a collapsing perturbation to have
             ! linear theory amplitude of order unity at the time of collapse, and since linear perturbations grow proportional to
             ! the expansion factor in an Einstein-de Sitter universe with no baryons, we use an initial guess for the lower and
             ! upper limits which are a multiple of our starting expansion factor.
             epsilonPerturbationMinimum=1.0d-1*expansionFactorInitialFraction
             epsilonPerturbationMaximum=1.0d+1*expansionFactorInitialFraction
             ! Evaluate cosmological parameters at the present time.
             OmegaMatterEpochal    =+     cosmologyFunctions_%omegaMatterEpochal     (expansionFactor=expansionFactor)
             if (self%baryonsCluster) then
                OmegaBaryonEpochal =+0.0d0
             else
                OmegaBaryonEpochal =+                          OmegaMatterEpochal                                      &
                     &              *self%cosmologyParameters_%OmegaBaryon           (                               ) &
                     &              /self%cosmologyParameters_%OmegaMatter           (                               )
             end if
             OmegaDarkEnergyEpochal=+     cosmologyFunctions_ %omegaDarkEnergyEpochal(expansionFactor=expansionFactor)
             hubbleTimeEpochal     =+     cosmologyFunctions_ %expansionRate         (                expansionFactor)
             time_                 =+     sphericalCollapse_  %x                     (                iTime          )
             ! Check dark energy equation of state is within acceptable range.
             if (cosmologyFunctions_%equationOfStateDarkEnergy(time=time_) >= -1.0d0/3.0d0) &
                  & call Error_Report('ω<-⅓ required'//{introspection:location})
             ! Find the value of ε for which the perturbation just collapses at this time.
             if (.not.finderPerturbationConstructed) then
                finderPerturbationInitial=rootFinder(                                                                 &
                     &                               rootFunction     =baryonsDarkMatterDarkEnergyRadiusPerturbation, &
                     &                               toleranceAbsolute=toleranceAbsolute                            , &
                     &                               toleranceRelative=toleranceRelative                            , &
                     &                               rangeExpandUpward=2.0d0                                        , &
                     &                               rangeExpandType  =rangeExpandMultiplicative                      &
                     &                              )
                finderPerturbationConstructed=.true.
             end if
             epsilonPerturbation=finderPerturbationInitial%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
             select case (calculationType%ID)
             case (cllsnlssMttCsmlgclCnstntClcltnCriticalOverdensity%ID)
                ! Critical linear overdensity.
                normalization=+linearGrowth_%value                                                                                (time_) &
                     &        /linearGrowth_%value(                                                                                       &
                     &                             cosmologyFunctions_%cosmicTime(                                                        &
                     &                                                            +expansionFactorInitialFraction                         &
                     &                                                            *cosmologyFunctions_            %expansionFactor(time_) &
                     &                                                           )                                                        &
                     &                            )
                call sphericalCollapse_%populate(                                   &
                     &                           normalization*epsilonPerturbation, &
                     &                           iTime                              &
                     &                          )
             case (cllsnlssMttCsmlgclCnstntClcltnVirialDensityContrast%ID,cllsnlssMttCsmlgclCnstntClcltnRadiusTurnaround%ID)
                ! Find the epoch of maximum expansion for the perturbation.
                if (.not.finderExpansionConstructed) then
                   finderExpansionMaximum=rootFinder(                                                                   &
                        &                            rootFunction=baryonsDarkMatterDarkEnergyExpansionRatePerturbation, &
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
                timeExpansionMaximum           =finderExpansionMaximum%find(rootRange=timeRange)
                expansionFactorExpansionMaximum=cosmologyFunctions_%expansionFactor(timeExpansionMaximum)
                ! Solve the dynamics of the perturbation to find the radius at the point of maximum expansion.
                call baryonsDarkMatterDarkEnergyPerturbationDynamicsSolver(epsilonPerturbation,timeExpansionMaximum,radiusExpansionMaximum)
                ! Compute the density contrast of the perturbation at maximum expansion.
                densityContrastExpansionMaximum=(expansionFactorExpansionMaximum/expansionFactor/radiusExpansionMaximum)**3
                ! Solve the cubic equation (Percival, 2005, A&A, 443, 819, eqn. 38; but modified to include the effects of baryons)
                ! to give the ratio of virial to turnaround radii, x.
                select case (self%energyFixedAt%ID)
                case (cllsnlssMttrDarkEnergyFixedAtTurnaround   %ID)
                   timeEnergyFixed=timeExpansionMaximum
                case (cllsnlssMttrDarkEnergyFixedAtVirialization%ID)
                   timeEnergyFixed=time_
                case default
                   call Error_Report('unrecognized epoch'//{introspection:location})
                end select
                if (self%baryonsCluster) then
                   q                 =     +cosmologyFunctions_%omegaDarkEnergyEpochal(time=timeExpansionMaximum) &
                        &                  /cosmologyFunctions_%omegaMatterEpochal    (time=timeExpansionMaximum) &
                        &                  /densityContrastExpansionMaximum
                   y                 =     +expansionFactorExpansionMaximum**cosmologyFunctions_%exponentDarkEnergy(time=timeExpansionMaximum) &
                        &                  /expansionFactor                **cosmologyFunctions_%exponentDarkEnergy(time=time_               )
                   a                 =+1.0d0-(1.0d0+3.0d0*cosmologyFunctions_%equationOfStateDarkEnergy(time=timeEnergyFixed))*q/2.0d0
                   b                 =      +(1.0d0+3.0d0*cosmologyFunctions_%equationOfStateDarkEnergy(time=time_          ))*q/y
                else
                   fractionDarkMatter=+(                                         &
                        &               +self%cosmologyParameters_%OmegaMatter() &
                        &               -self%cosmologyParameters_%OmegaBaryon() &
                        &              )                                         &
                        &             /  self%cosmologyParameters_%OmegaMatter()
                   q                 =+cosmologyFunctions_%omegaDarkEnergyEpochal(time=timeExpansionMaximum) &
                        &             /cosmologyFunctions_%omegaMatterEpochal    (time=timeExpansionMaximum) &
                        &             /fractionDarkMatter                                                    &
                        &             /densityContrastExpansionMaximum
                   y                 = expansionFactorExpansionMaximum**cosmologyFunctions_%exponentDarkEnergy(time=timeExpansionMaximum) &
                        &             /expansionFactor                **cosmologyFunctions_%exponentDarkEnergy(time=time_               )
                   r                 =+  self%cosmologyParameters_%OmegaBaryon() &
                        &             /(                                         &
                        &               +self%cosmologyParameters_%OmegaMatter() &
                        &               -self%cosmologyParameters_%OmegaBaryon() &
                        &              )                                         &
                        &             /densityContrastExpansionMaximum
                   z                 =+(                                    &
                        &               +expansionFactorExpansionMaximum    &
                        &               /expansionFactor                    &
                        &              )**3
                   a                 =+1.0d0+r        -(1.0d0+3.0d0*cosmologyFunctions_%equationOfStateDarkEnergy(time=timeEnergyFixed))*q/2.0d0
                   b                 =      -r/z/2.0d0+(1.0d0+3.0d0*cosmologyFunctions_%equationOfStateDarkEnergy(time=time_          ))*q/y
                end if
                x      =+(0.0d0,0.5d0)*sqrt(3.0d0)                                                                          &
                     &  *(                                                                                                  &
                     &    +1.0d0/b*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(+1.0d0/3.0d0)/ 6.0d0  &
                     &    +2.0d0*a*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(-1.0d0/3.0d0)         &
                     &   )                                                                                                  &
                     &  - (1.0d0/b*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(+1.0d0/3.0d0)/12.0d0) &
                     &  + (      a*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(-1.0d0/3.0d0)       )
                select case (calculationType%ID)
                case (cllsnlssMttCsmlgclCnstntClcltnVirialDensityContrast%ID)
                   ! The density contrast calculated as Δ=1/(x Rmax)³ is Δ=ρvir/⟨ρDM⟩ - i.e. the density of the virialized dark
                   ! matter perturbation relative to the mean density of dark matter. However, what we want (for the definition used
                   ! by Galacticus) is the density of the perturbation relative to the total mean density. So we perform that
                   ! conversion here.
                   call sphericalCollapse_%populate(                                            &
                        &                           +(                                          &
                        &                             +self%cosmologyParameters_%OmegaMatter()  &
                        &                             -self%cosmologyParameters_%OmegaBaryon()  &
                        &                            )                                          &
                        &                           /  self%cosmologyParameters_%OmegaMatter()  &
                        &                           /(dble(x)*radiusExpansionMaximum)**3      , &
                        &                           iTime                                       &
                        &                          )
                case (cllsnlssMttCsmlgclCnstntClcltnRadiusTurnaround    %ID)
                   call sphericalCollapse_%populate(                                            &
                        &                           1.0d0/ dble(x)                            , &
                        &                           iTime                                       &
                        &                          )
                end select
             end select
             !$omp atomic
             iCount=iCount+1
          end if
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
       call displayUnindent     ('done',verbosity=verbosityLevelWorking)
    end select
    return
  end subroutine baryonsDarkMatterDarkEnergyTabulate

  double precision function baryonsDarkMatterDarkEnergyRadiusPerturbation(epsilonPerturbation)
    !!{
    Return the radius of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    amplitude {\normalfont \ttfamily epsilonPerturbation}.
    !!}
    implicit none
    double precision, intent(in   ) :: epsilonPerturbation

    call baryonsDarkMatterDarkEnergyPerturbationDynamicsSolver(epsilonPerturbation,time_,baryonsDarkMatterDarkEnergyRadiusPerturbation)
    return
  end function baryonsDarkMatterDarkEnergyRadiusPerturbation

  double precision function baryonsDarkMatterDarkEnergyExpansionRatePerturbation(time)
    !!{
    Return the expansion rate of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    amplitude {\normalfont \ttfamily epsilonPerturbation}.
    !!}
    implicit none
    double precision, intent(in   ) :: time

    call baryonsDarkMatterDarkEnergyPerturbationDynamicsSolver(amplitudePerturbation,time,expansionRatePerturbation=baryonsDarkMatterDarkEnergyExpansionRatePerturbation)
    return
  end function baryonsDarkMatterDarkEnergyExpansionRatePerturbation

  subroutine baryonsDarkMatterDarkEnergyPerturbationDynamicsSolver(perturbationOverdensityInitial,time,radiusPerturbation,expansionRatePerturbation)
    !!{
    Integrate the dynamics of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    amplitude {\normalfont \ttfamily epsilonPerturbation}.
    !!}
    use :: Error                , only : Error_Report
    use :: Interface_GSL        , only : GSL_Success
    use :: Numerical_ODE_Solvers, only : odeSolver
    implicit none
    double precision           , intent(in   )                        :: perturbationOverdensityInitial           , time
    double precision           , intent(  out)             , optional :: expansionRatePerturbation                , radiusPerturbation
    integer         (c_size_t ), parameter                            :: countProperties               =2_c_size_t
    double precision           , dimension(countProperties)           :: propertyValues
    double precision           , parameter                            :: odeToleranceAbsolute          =0.0d0     , odeToleranceRelative             =1.0d-12
    type            (odeSolver)                                       :: solver
    double precision                                                  :: expansionFactorInitial                   , radiusDifferenceGrowthRateInitial        , &
         &                                                               timeInitial                              , exponent                                 , &
         &                                                               radiusDifferenceInitial
    integer                                                           :: odeStatus

    ! Validate the perturbation overdensity.
    if (perturbationOverdensityInitial < 0.0d+0) call Error_Report('initial overdensity of perturbation should be non-negative'//{introspection:location})
    if (perturbationOverdensityInitial > 1.0d-3) call Error_Report('initial overdensity of perturbation should be small'       //{introspection:location})
    ! Specify a sufficiently early time.
    expansionFactorInitial=expansionFactorInitialFraction
    ! Find the corresponding cosmic time.
    timeInitial=cosmologyFunctions_%cosmicTime(expansionFactorInitial*cosmologyFunctions_%expansionFactor(time_))
    ! Determine the initial growth rate of the overdensity assuming the matter-dominated phase growth factor is proportional to the expansion factor.
    if (baryonsCluster) then
       ! Baryons are assumed to cluster, so the usual matter dominated solution of D(a)=a applies.
       exponent=1.0d0
    else
       ! Baryons do not cluster, the growth factor is D(a)=a^{3p/2} with p=[sqrt(1+24fₓ)-1]/6
       exponent=0.25d0*(sqrt(1.0d0+24.0d0*(OmegaMatterEpochal-OmegaBaryonEpochal)/OmegaMatterEpochal)-1.0d0)
    end if
    ! We solve the equations for the evolution of y(t) = a_p(t) - a(t) - where a_p(t) is the radius of our perturbation, and a(t)
    ! is the expansion factor. Using y(t) as our variable allows us to maintain precision at early times when y(t) is very small -
    ! if we instead solve directly for a_p(t) numerical issues prevent a numerically robust solution from being obtained.
    !! Find the initial radius parameter, y(t) of the perturbation.
    radiusDifferenceInitial          =+      expansionFactorInitial                                                &
         &                            *(                                                                           &
         &                              -( 1.0d0/  3.0d0)                       *perturbationOverdensityInitial    &
         &                              +( 2.0d0/  9.0d0)                       *perturbationOverdensityInitial**2 &
         &                              -(14.0d0/ 81.0d0)                       *perturbationOverdensityInitial**3 &
         &                              +(35.0d0/243.0d0)                       *perturbationOverdensityInitial**4 &
         &                              -(91.0d0/729.0d0)                       *perturbationOverdensityInitial**5 &
         &                             )
    ! Find the initial growth rate of the perturbation radius.
    radiusDifferenceGrowthRateInitial=+      hubbleTimeEpochal                                                     &
         &                            *sqrt(                                                                       &
         &                                  +OmegaMatterEpochal                                                    &
         &                                  /expansionFactorInitial**3                                             &
         &                                 )                                                                       &
         &                            *      expansionFactorInitial                                                &
         &                            *(                                                                           &
         &                              -( 1.0d0/  3.0d0)*(1.0d0+      exponent)*perturbationOverdensityInitial    &
         &                              +( 2.0d0/  9.0d0)*(1.0d0+2.0d0*exponent)*perturbationOverdensityInitial**2 &
         &                              -(14.0d0/ 81.0d0)*(1.0d0+3.0d0*exponent)*perturbationOverdensityInitial**3 &
         &                              +(35.0d0/243.0d0)*(1.0d0+4.0d0*exponent)*perturbationOverdensityInitial**4 &
         &                              -(91.0d0/729.0d0)*(1.0d0+5.0d0*exponent)*perturbationOverdensityInitial**5 &
         &                             )
    ! Set initial conditions.
    propertyValues=[radiusDifferenceInitial,radiusDifferenceGrowthRateInitial]
    ! Evolve if the requested time is after the initial time.
    if (time > timeInitial) then
       ! Solve the ODE to find the perturbation radius at the present day.
       solver=odeSolver(countProperties,baryonsDarkMatterDarkEnergyPerturbationODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
       call solver%solve(timeInitial,time,propertyValues,status=odeStatus)
       ! If the ODE solver did not succeed, it is because the perturbation collapsed to zero radius (causing a divergence). This
       ! means it collapsed prior to the current time. We extrapolate to negative radius (using the velocity at the final step) to
       ! permit our root finder to locate the point at which collapse occurs at the current time.
       if (odeStatus /= GSL_Success) propertyValues(1)=+propertyValues(1) &
            &                                          +propertyValues(2) &
            &                                          *(                 &
            &                                            +time            &
            &                                            -timeInitial     &
            &                                          )
    end if
    ! Return the radius and/or expansion rate of the perturbation. Note that here we add back the expansion factor (or its growth
    ! rate) since we've solved the the radius and expansion rate of the perturbation relative to the expansion factor.
    if (present(radiusPerturbation       )) radiusPerturbation       =+propertyValues                                                       (   1)  &
         &                                                            +                                  cosmologyFunctions_%expansionFactor(time)  &
         &                                                            /                                  cosmologyFunctions_%expansionFactor(time_)
    if (present(expansionRatePerturbation)) expansionRatePerturbation=+propertyValues                                                       (   2)  &
         &                                                            +cosmologyFunctions_%expansionRate(cosmologyFunctions_%expansionFactor(time)) &
         &                                                            *                                  cosmologyFunctions_%expansionFactor(time)  &
         &                                                            /                                  cosmologyFunctions_%expansionFactor(time_)
    return
  end subroutine baryonsDarkMatterDarkEnergyPerturbationDynamicsSolver

  integer function baryonsDarkMatterDarkEnergyPerturbationODEs(time,y,dydt)
    !!{
    Differential equations describing the evolution of spherical perturbations in a universe containing baryons, collisionless dark matter and dark energy.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision, intent(in   )               :: time
    double precision, intent(in   ), dimension(:) :: y
    double precision, intent(  out), dimension(:) :: dydt
    double precision                              :: expansionFactor

    expansionFactor  =+cosmologyFunctions_%expansionFactor(time) &
         &            /cosmologyFunctions_%expansionFactor(time_)
    if (y(1) <= -expansionFactor) then
       dydt(1:2)=0.0d0
    else
       dydt(1  )=+  y(2)
       dydt(2  )=-0.5d0                                                                                                                                                                                                                    &
            &    *hubbleTimeEpochal**2                                                                                                                                                                                                     &
            &    *(                                                                                                                                                                                                                        &
            &      +y(1)                                                                                                                                                                                                                   &
            &      *(                                                                                                                                                                                                                      &
            &         +(OmegaMatterEpochal    -OmegaBaryonEpochal)                                                                       *(y(1)+expansionFactor)**(-3)                                                                     &
            &         +                        OmegaBaryonEpochal                                                                        *                              expansionFactor**(-3)                                              &
            &         + OmegaDarkEnergyEpochal                    *(3.0d0*cosmologyFunctions_%equationOfStateDarkEnergy(time=time)+1.0d0)*                              expansionFactor**cosmologyFunctions_%exponentDarkEnergy(time=time) &
            &        )                                                                                                                                                                                                                     &
            &      +expansionFactor                                                                                                                                                                                                        &
            &      *(                                                                                                                                                                                                                      &
            &         +(OmegaMatterEpochal    -OmegaBaryonEpochal)                                                                       *((y(1)+expansionFactor)**(-3)-expansionFactor**(-3))                                             &
            &        )                                                                                                                                                                                                                     &
            &     )
    end if
    ! Return success.
    baryonsDarkMatterDarkEnergyPerturbationODEs=GSL_Success
    return
  end function baryonsDarkMatterDarkEnergyPerturbationODEs
