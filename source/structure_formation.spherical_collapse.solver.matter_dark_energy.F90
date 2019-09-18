!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% A spherical collapse solver class for universes consisting of collisionless matter and dark energy.

  !# <sphericalCollapseSolver name="sphericalCollapseSolverMatterDarkEnergy">
  !#  <description>A spherical collapse solver for universes consisting of collisionless matter and dark energy.</description>
  !# </sphericalCollapseSolver>
  type, extends(sphericalCollapseSolverMatterLambda) :: sphericalCollapseSolverMatterDarkEnergy
     !% A spherical collapse solver for universes consisting of collisionless matter and dark energy.
     private
     integer :: energyFixedAt
   contains
     procedure :: linearNonlinearMap => matterDarkEnergyLinearNonlinearMap
     procedure :: tabulate           => matterDarkEnergyTabulate
  end type sphericalCollapseSolverMatterDarkEnergy

  interface sphericalCollapseSolverMatterDarkEnergy
     !% Constructors for the {\normalfont \ttfamily matterDarkEnergy} spherical collapse solver class.
     module procedure matterDarkEnergyConstructorParameters
     module procedure matterDarkEnergyConstructorInternal
  end interface sphericalCollapseSolverMatterDarkEnergy

  ! Enumeration of radii at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.
  !# <enumeration>
  !#  <name>matterDarkEnergyFixedAt</name>
  !#  <description>Enumeration of radii at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <validator>yes</validator>
  !#  <visibility>public</visibility>
  !#  <entry label="undefined"    />
  !#  <entry label="turnaround"   />
  !#  <entry label="virialization"/>
  !# </enumeration>

  ! Pointer to the default cosmology functions object.
  class           (cosmologyFunctionsClass), pointer   :: matterDarkEnergyCosmologyFunctions_
  !$omp threadprivate(matterDarkEnergyCosmologyFunctions_)

  ! Fraction of current expansion factor to use as initial time in perturbation dynamics solver.
  double precision                         , parameter :: matterDarkEnergyExpansionFactorInitialFraction=1.0d-6

  ! Variables used in root finding.
  double precision                                     :: matterDarkEnergyPerturbationRadiusInitial
  !$omp threadprivate(matterDarkEnergyPerturbationRadiusInitial)
  
contains

  function matterDarkEnergyConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily matterDarkEnergy} spherical collapse solver class that takes a parameter set as
    !% input.
    use Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (sphericalCollapseSolverMatterDarkEnergy)                :: self
    type   (inputParameters                        ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class  (linearGrowthClass                      ), pointer       :: linearGrowth_
    type   (varying_string                         )                :: energyFixedAt

    !# <inputParameter>
    !#   <name>energyFixedAt</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('turnaround')</defaultValue>
    !#   <description>The radius at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.</description>
    !#   <type>integer</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !# <objectBuilder class="linearGrowth"       name="linearGrowth_"       source="parameters"/>
    self=sphericalCollapseSolverMatterDarkEnergy(enumerationMatterDarkEnergyFixedAtEncode(char(energyFixedAt),includesPrefix=.false.),cosmologyFunctions_,linearGrowth_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"/>
    !# <objectDestructor name="linearGrowth_"      />
    return
  end function matterDarkEnergyConstructorParameters

  function matterDarkEnergyConstructorInternal(energyFixedAt,cosmologyFunctions_,linearGrowth_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily matterDarkEnergy} spherical collapse solver class.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type   (sphericalCollapseSolverMatterDarkEnergy)                                  :: self
    integer                                         , intent(in   )                   :: energyFixedAt
    class  (cosmologyFunctionsClass                ), intent(in   ), target           :: cosmologyFunctions_
    class  (linearGrowthClass                      ), intent(in   ), target, optional :: linearGrowth_
    !# <constructorAssign variables="energyFixedAt, *cosmologyFunctions_, *linearGrowth_"/>

    if (.not.enumerationMatterDarkEnergyFixedAtIsValid(energyFixedAt)) call Galacticus_Error_Report('invalid energyFixedAt'//{introspection:location})
    return
  end function matterDarkEnergyConstructorInternal

  subroutine matterDarkEnergyTabulate(self,time,sphericalCollapse_,calculationType)
    !% Tabulate spherical collapse solutions for $\delta_\mathrm{crit}$, $\Delta_\mathrm{vir}$, or $R_\mathrm{ta}/R_\mathrm{vir}$ vs. time.
    use Root_Finder       , only : rootFinder               , rangeExpandMultiplicative  , rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    use Galacticus_Error  , only : Galacticus_Error_Report
    use Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Counter   , Galacticus_Display_Counter_Clear, &
         &                         verbosityWorking
    use Tables            , only : table1DLogarithmicLinear
    use Linear_Growth     , only : normalizeMatterDominated
    implicit none
    class           (sphericalCollapseSolverMatterDarkEnergy)             , intent(inout) :: self
    double precision                                                      , intent(in   ) :: time
    integer                                                               , intent(in   ) :: calculationType
    class           (table1D                                ), allocatable, intent(inout) :: sphericalCollapse_
    class           (linearGrowthClass                      ), pointer                    :: linearGrowth_
    double precision                                         , parameter                  :: toleranceAbsolute              =0.0d0, toleranceRelative              =1.0d-9
    double precision                                         , dimension(2)               :: timeRange
    type            (rootFinder                             ), save                       :: finderAmplitudePerturbation          , finderExpansionMaximum
    !$omp threadprivate(finderAmplitudePerturbation,finderExpansionMaximum)
    integer                                                                               :: countTimes                           , iTime                                 , &
         &                                                                                   iCount
    double precision                                                                      :: expansionFactor                      , epsilonPerturbation                   , &
         &                                                                                   epsilonPerturbationMaximum           , epsilonPerturbationMinimum            , &
         &                                                                                   densityContrastExpansionMaximum      , expansionFactorExpansionMaximum       , &
         &                                                                                   radiusExpansionMaximum               , maximumExpansionTime                  , &
         &                                                                                   normalization                        , q                                     , &
         &                                                                                   timeEnergyFixed                      , timeInitial                           , &
         &                                                                                   y                                    , timeMinimum                           , &
         &                                                                                   timeMaximum
    double complex                                                                        :: a                                    , b                                     , &
         &                                                                                   x
    type            (varying_string                         )                             :: message
    character       (len=13                                 )                             :: label

    ! Validate.
    if (calculationType == matterLambdaCalculationCriticalOverdensity .and. .not.associated(self%linearGrowth_)) call Galacticus_Error_Report('linearGrowth object was not provided'//{introspection:location})
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
    countTimes=int(log10(timeMaximum/timeMinimum)*dble(matterLambdaTablePointsPerDecade))
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
       call Galacticus_Display_Indent(message,verbosity=verbosityWorking)
       iCount=0
       call Galacticus_Display_Counter(                            &
            &                                    iCount          , &
            &                          isNew    =.true.          , &
            &                          verbosity=verbosityWorking  &
            &                         )
       !$omp parallel private(expansionFactor,epsilonPerturbationMaximum,epsilonPerturbationMinimum,epsilonPerturbation,timeInitial,timeRange,maximumExpansionTime,expansionFactorExpansionMaximum,q,y,timeEnergyFixed,a,b,x,linearGrowth_)
       allocate(matterDarkEnergyCosmologyFunctions_,mold=self%cosmologyFunctions_)
       !# <deepCopy source="self%cosmologyFunctions_" destination="matterDarkEnergyCosmologyFunctions_"/>
       if (calculationType == matterLambdaCalculationCriticalOverdensity) then
          allocate(linearGrowth_,mold=self%linearGrowth_)
          !# <deepCopy source="self%linearGrowth_" destination="linearGrowth_"/>
       end if
       !$omp do schedule(dynamic)
       do iTime=1,countTimes
          call Galacticus_Display_Counter(                                              &
               &                          int(100.0d0*dble(iCount-1)/dble(countTimes)), &
               &                          isNew=.false.                               , &
               &                          verbosity=verbosityWorking                    &
               &                         )
          ! Get the current expansion factor.
          expansionFactor=matterDarkEnergyCosmologyFunctions_%expansionFactor(sphericalCollapse_%x(iTime))
          ! In the case of dark energy we cannot (easily) determine the largest (i.e. least negative) value of ε for which a
          ! perturbation can collapse. So, use no perturbation.
          epsilonPerturbationMaximum=  0.0d0
          ! Estimate a suitably negative minimum value for ε.
          epsilonPerturbationMinimum=-10.0d0
          ! Evaluate cosmological parameters at the present time.
          matterLambdaOmegaMatterEpochal    =matterDarkEnergyCosmologyFunctions_%omegaMatterEpochal    (expansionFactor=expansionFactor)
          matterLambdaOmegaDarkEnergyEpochal=matterDarkEnergyCosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor)
          matterLambdaHubbleTimeEpochal     =matterDarkEnergyCosmologyFunctions_%expansionRate         (                expansionFactor)
          matterLambdaTime                  =sphericalCollapse_                 %x                     (                iTime          )
          ! Check dark energy equation of state is within acceptable range.
          if (matterDarkEnergyCosmologyFunctions_%equationOfStateDarkEnergy(time=matterLambdaTime) >= -1.0d0/3.0d0) &
               & call Galacticus_Error_Report('ω<-⅓ required'//{introspection:location})
          ! Find the value of ε for which the perturbation just collapses at this time.
          if (.not.finderAmplitudePerturbation%isInitialized()) then
             call finderAmplitudePerturbation%rootFunction(matterDarkEnergyRadiusPerturbation                 )
             call finderAmplitudePerturbation%tolerance   (toleranceAbsolute,toleranceRelative)
          end if
          epsilonPerturbation=finderAmplitudePerturbation%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
          select case (calculationType)
          case (matterLambdaCalculationCriticalOverdensity)
             ! Critical linear overdensity.
             normalization=+linearGrowth_%value(matterLambdaTime,normalize=normalizeMatterDominated) &
                  &        /                    expansionFactor
             call sphericalCollapse_%populate(                                       &
                  &                           +normalization                         &
                  &                           *0.6d0                                 &
                  &                           *(                                     &
                  &                             +1.0d0                               &
                  &                             -matterLambdaOmegaMatterEpochal      &
                  &                             -matterLambdaOmegaDarkEnergyEpochal  &
                  &                             -epsilonPerturbation                 &
                  &                            )                                     &
                  &                           /matterLambdaOmegaMatterEpochal      , &
                  &                           iTime                                  &
                  &                          )
          case (matterLambdaCalculationVirialDensityContrast,matterLambdaCalculationRadiusTurnaround)
             ! Find the epoch of maximum expansion for the perturbation.
             if (.not.finderExpansionMaximum%isInitialized()) then
                call finderExpansionMaximum%rootFunction(matterDarkEnergyExpansionRatePerturbation)
                call finderExpansionMaximum%tolerance   (toleranceAbsolute,toleranceRelative      )
             end if
             call finderExpansionMaximum%rangeExpand (                                                             &
                  &                                   rangeExpandDownward          =1.0d0-1.0d-2                 , &
                  &                                   rangeExpandUpward            =1.0d0+1.0d-2                 , &
                  &                                   rangeExpandType              =rangeExpandMultiplicative    , &
                  &                                   rangeUpwardLimit             =matterLambdaTime             , &
                  &                                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
                  &                                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative  &
                  &                                  )
             matterLambdaAmplitudePerturbation=epsilonPerturbation
             ! Compute the corresponding time of maximum expansion.
             timeInitial                    =matterDarkEnergyCosmologyFunctions_%cosmicTime(matterDarkEnergyCosmologyFunctions_%expansionFactor(matterLambdaTime)*matterDarkEnergyExpansionFactorInitialFraction)
             ! Guess that the time of maximum expansion occurred at close to half of the current time.
             timeRange                      =[0.45d0,0.55d0]*matterLambdaTime
             maximumExpansionTime           =finderExpansionMaximum%find(rootRange=timeRange)
             expansionFactorExpansionMaximum=matterDarkEnergyCosmologyFunctions_%expansionFactor(maximumExpansionTime)
             ! Solve the dynamics of the perturbation to find the radius at the point of maximum expansion.
             call matterDarkEnergyPerturbationDynamicsSolver(epsilonPerturbation,maximumExpansionTime,radiusExpansionMaximum)
             ! Compute the density contrast of the perturbation at maximum expansion.
             densityContrastExpansionMaximum=(expansionFactorExpansionMaximum/expansionFactor/radiusExpansionMaximum)**3
             ! Solve the cubic equation (Percival, 2005, A&A, 443, 819, eqn. 38) to give the ratio of virial to turnaround radii,
             ! x.
             q=      matterDarkEnergyCosmologyFunctions_%omegaDarkEnergyEpochal(time=maximumExpansionTime) &
                  & /matterDarkEnergyCosmologyFunctions_%omegaMatterEpochal    (time=maximumExpansionTime) &
                  & /densityContrastExpansionMaximum
             y=      expansionFactorExpansionMaximum**matterDarkEnergyCosmologyFunctions_%exponentDarkEnergy(time=maximumExpansionTime) &
                  & /expansionFactor                **matterDarkEnergyCosmologyFunctions_%exponentDarkEnergy(time=matterLambdaTime    )
             select case (self%energyFixedAt)
             case (matterDarkEnergyFixedAtTurnaround   )
                timeEnergyFixed=maximumExpansionTime
             case (matterDarkEnergyFixedAtVirialization)
                timeEnergyFixed=matterLambdaTime
             case default
                call Galacticus_Error_Report('unrecognized epoch'//{introspection:location})
             end select
             a=1.0d0-(1.0d0+3.0d0*matterDarkEnergyCosmologyFunctions_%equationOfStateDarkEnergy(time=timeEnergyFixed ))*q/2.0d0
             b=      (1.0d0+3.0d0*matterDarkEnergyCosmologyFunctions_%equationOfStateDarkEnergy(time=matterLambdaTime))*q/y
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
             select case (calculationType)
             case (matterLambdaCalculationVirialDensityContrast)
                call sphericalCollapse_%populate(                                           &
                     &                           1.0d0/(dble(x)*radiusExpansionMaximum)**3, &
                     &                           iTime                                      &
                     &                          )
             case (matterLambdaCalculationRadiusTurnaround)
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
       !# <objectDestructor name="matterDarkEnergyCosmologyFunctions_"/>
       !$omp end parallel
       call Galacticus_Display_Counter_Clear(       verbosity=verbosityWorking)
       call Galacticus_Display_Unindent     ('done',verbosity=verbosityWorking)
    end select
    return
  end subroutine matterDarkEnergyTabulate

  double precision function matterDarkEnergyRadiusPerturbation(epsilonPerturbation)
    !% Return the radius of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\normalfont \ttfamily epsilonPerturbation}.
    implicit none
    double precision, intent(in   ) :: epsilonPerturbation

    call matterDarkEnergyPerturbationDynamicsSolver(epsilonPerturbation,matterLambdaTime,matterDarkEnergyRadiusPerturbation)
    return
  end function matterDarkEnergyRadiusPerturbation

  double precision function matterDarkEnergyExpansionRatePerturbation(time)
    !% Return the expansion rate of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\normalfont \ttfamily epsilonPerturbation}.
    implicit none
    double precision, intent(in   ) :: time

    call matterDarkEnergyPerturbationDynamicsSolver(matterLambdaAmplitudePerturbation,time,expansionRatePerturbation=matterDarkEnergyExpansionRatePerturbation)
    return
  end function matterDarkEnergyExpansionRatePerturbation

  subroutine matterDarkEnergyPerturbationDynamicsSolver(epsilonPerturbation,time,radiusPerturbation,expansionRatePerturbation)
    !% Integrate the dynamics of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\normalfont \ttfamily epsilonPerturbation}.
    use FGSL         , only : FGSL_Success
    use ODEIV2_Solver, only : ODEIV2_Solve  , ODEIV2_Solver_Free
    use FODEIV2      , only : fodeiv2_system, fodeiv2_driver
    implicit none
    double precision                                                , intent(in   )           :: epsilonPerturbation            , time
    double precision                                                , intent(  out), optional :: expansionRatePerturbation      , radiusPerturbation
    integer                             , parameter                                           :: countProperties          =2
    double precision                    , dimension(countProperties)                          :: propertyValues
    double precision                    , parameter                                           :: odeToleranceAbsolute     =0.0d0, odeToleranceRelative            =1.0d-12
    type            (fodeiv2_system    )                                                      :: ode2System
    type            (fodeiv2_driver    )                                                      :: ode2Driver
    logical                                                                                   :: odeReset
    double precision                                                                          :: expansionFactorInitial         , expansionRatePerturbationInitial        , &
         &                                                                                       overdensityInitial             , timeInitial
    integer                                                                                   :: odeStatus

    ! Specify a sufficiently early time.
    expansionFactorInitial=matterDarkEnergyExpansionFactorInitialFraction
    ! Find the corresponding cosmic time.
    timeInitial=matterDarkEnergyCosmologyFunctions_%cosmicTime(matterDarkEnergyCosmologyFunctions_%expansionFactor(time)*expansionFactorInitial)
    ! Find the overdensity of the perturbation at early time (Percival, 2005, A&A, 443, 819, eqn. 25).
    overdensityInitial=+0.6d0                                &
         &             *(                                    &
         &               +1.0d0                              &
         &               -matterLambdaOmegaMatterEpochal     &
         &               -matterLambdaOmegaDarkEnergyEpochal &
         &               -epsilonPerturbation                &
         &              )                                    &
         &             *(                                    &
         &               +3.0d0                              &
         &               /2.0d0                              &
         &               /matterLambdaOmegaMatterEpochal     &
         &               *matterLambdaHubbleTimeEpochal      &
         &               *timeInitial                        &
         &              )**(2.0d0/3.0d0)
    ! Find the perturbation radius at this early time. This is, by construction, just the initial expansion factor.
    matterDarkEnergyPerturbationRadiusInitial=+expansionFactorInitial
    ! Find the perturbation expansion rate at early time (Percival, 2005, A&A, 443, 819, eqn. 22).
    expansionRatePerturbationInitial         =+matterLambdaHubbleTimeEpochal        &
         &                                    *sqrt(                                &
         &                                          +matterLambdaOmegaMatterEpochal &
         &                                          /expansionFactorInitial         &
         &                                          +epsilonPerturbation            &
         &                                         )
    ! Set initial conditions.
    propertyValues=[                                           &
         &          matterDarkEnergyPerturbationRadiusInitial, &
         &          expansionRatePerturbationInitial           &
         &         ]
    ! Evolve if the requested time is after the initial time.
    if (time > timeInitial) then
       ! Solve the ODE to find the perturbation radius at the present day.
       odeReset=.true.
       call ODEIV2_Solve(                                                  &
            &                   ode2Driver,ode2System                    , &
            &                   timeInitial,time                         , &
            &                   countProperties                          , &
            &                   propertyValues                           , &
            &                   matterDarkEnergyPerturbationODEs         , &
            &                   odeToleranceAbsolute,odeToleranceRelative, &
            &                   reset=odeReset                           , &
            &                   odeStatus=odeStatus                        &
            &                  )
       call ODEIV2_Solver_Free(ode2Driver,ode2System)
       ! If the ODE solver did not succeed, it is because the perturbation collapsed to zero radius (causing a divergence). This
       ! means it collapsed prior to the current time. We extrapolate to negative radius (using the velocity at the final step) to
       ! permit our root finder to locate the point at which collapse occurs at the current time.
       if (odeStatus /= FGSL_Success) propertyValues(1)=+propertyValues(1) &
            &                                           +propertyValues(2) &
            &                                           *(                 &
            &                                             +time            &
            &                                             -timeInitial     &
            &                                            )
    end if
    ! Return the computed quantities.
    if (present(radiusPerturbation       )) radiusPerturbation       =propertyValues(1)
    if (present(expansionRatePerturbation)) expansionRatePerturbation=propertyValues(2)
    return
  end subroutine matterDarkEnergyPerturbationDynamicsSolver

  integer function matterDarkEnergyPerturbationODEs(time,y,dydt)
    !% Differential equations describing the evolution of spherical perturbations in a universe containing collisionless dark matter and dark energy.
    use FGSL, only : FGSL_Success
    implicit none
    double precision, intent(in   )               :: time
    double precision, intent(in   ), dimension(:) :: y
    double precision, intent(  out), dimension(:) :: dydt
    double precision                              :: expansionFactor

    if (y(1) <= 0.0d0) then
       dydt(1:2)=0.0d0
    else
       expansionFactor=+matterDarkEnergyCosmologyFunctions_%expansionFactor(            time) &
            &          /matterDarkEnergyCosmologyFunctions_%expansionFactor(matterLambdaTime)
       dydt(1)=+y(2)
       dydt(2)=-0.5d0                                                                                                                                                                                                          &
            &  *y(1)                                                                                                                                                                                                           &
            &  *matterLambdaHubbleTimeEpochal**2                                                                                                                                                                               &
            &  *(                                                                                                                                                                                                              &
            &    +                                                                                       matterLambdaOmegaMatterEpochal    /y(1)           **3                                                                 &
            &    +(3.0d0*matterDarkEnergyCosmologyFunctions_%equationOfStateDarkEnergy(time=time)+1.0d0)*matterLambdaOmegaDarkEnergyEpochal*expansionFactor**matterDarkEnergyCosmologyFunctions_%exponentDarkEnergy(time=time) &
            &   )
    end if
    ! Return success.
    matterDarkEnergyPerturbationODEs=FGSL_Success
    return
  end function matterDarkEnergyPerturbationODEs

  subroutine matterDarkEnergyLinearNonlinearMap(self,time,linearNonlinearMap_)
    !% Tabulate the mapping between linea rna dnonlinear overdensity for the spherical collapse model. 
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (sphericalCollapseSolverMatterDarkEnergy), intent(inout) :: self
    double precision                                         , intent(in   ) :: time
    class           (table2DLinLinLin                       ), intent(inout) :: linearNonlinearMap_
    !GCC$ attributes unused :: self, time, linearNonlinearMap_

    call Galacticus_Error_Report('linear-nonlinear mapping is not supported by this class'//{introspection:location})
    return
  end subroutine matterDarkEnergyLinearNonlinearMap
