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

!% Contains a module which implements calculations of spherical top hat collapse in cosmologies containing matter and dark energy.

module Spherical_Collapse_Matter_Dark_Energy
  !% Implements calculations of spherical top hat collapse in cosmologies containing matter and dark energy.
  use FGSL               , only : FGSL_Success
  use ISO_Varying_String
  use Cosmology_Functions
  implicit none
  private
  public :: Spherical_Collapse_Dark_Energy_Critical_Overdensity_Tabulate, Spherical_Collapse_Dark_Energy_Virial_Density_Contrast_Tabulate, &
       &    Spherical_Collapse_Dark_Energy_Turnaround_Radius_Tabulate

  ! Variables to hold the tabulated critical overdensity data.
  integer                                  , parameter     :: deltaTableNPointsPerDecade    =100

  ! Variables used in root finding.
  double precision                                         :: OmegaDE                              , OmegaM                      , &
       &                                                      epsilonPerturbationShared            , hubbleParameterInvGyr       , &
       &                                                      perturbationRadiusInitial            , tNow
  !$omp threadprivate(OmegaDE,OmegaM,hubbleParameterInvGyr,tNow,epsilonPerturbationShared,perturbationRadiusInitial)
  
  ! Fraction of current expansion factor to use as initial time in perturbation dynamics solver.
  double precision                         , parameter     :: expansionFactorInitialFraction=1.0d-6

  ! Calculation types.
  integer                                  , parameter     :: calculationDeltaCrit          =0     , calculationDeltaVirial=1, &
       &                                                      calculationTurnaround         =2

  ! Pointer to the default cosmology functions object.
  class           (cosmologyFunctionsClass), allocatable   :: cosmologyFunctions__
  !$omp threadprivate(cosmologyFunctions__)

  ! Enumeration of radii at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.
  !# <enumeration>
  !#  <name>darkEnergySphericalCollapseEnergyFixedAt</name>
  !#  <description>Enumeration of radii at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <validator>yes</validator>
  !#  <visibility>public</visibility>
  !#  <entry label="turnaround"   />
  !#  <entry label="virialization"/>
  !# </enumeration>

contains

  subroutine Spherical_Collapse_Dark_Energy_Critical_Overdensity_Tabulate(time,deltaCritTable,cosmologyFunctions_,linearGrowth_)
    !% Tabulate the critical overdensity for collapse for the spherical collapse model.
    use Tables
    use Linear_Growth
    implicit none
    double precision                                      , intent(in   ) :: time
    class           (table1D                ), allocatable, intent(inout) :: deltaCritTable
    class           (cosmologyFunctionsClass)             , intent(inout) :: cosmologyFunctions_    
    class           (linearGrowthClass      )             , intent(inout) :: linearGrowth_    

    call Make_Table(time,deltaCritTable,calculationDeltaCrit,cosmologyFunctions_,linearGrowth_=linearGrowth_)
    return
  end subroutine Spherical_Collapse_Dark_Energy_Critical_Overdensity_Tabulate

  subroutine Spherical_Collapse_Dark_Energy_Virial_Density_Contrast_Tabulate(time,energyFixedAt,deltaVirialTable,cosmologyFunctions_)
    !% Tabulate the virial density contrast for the spherical collapse model.
    use Tables
    implicit none
    double precision                                      , intent(in   ) :: time
    integer                                               , intent(in   ) :: energyFixedAt
    class           (table1D                ), allocatable, intent(inout) :: deltaVirialTable
    class           (cosmologyFunctionsClass)             , intent(inout) :: cosmologyFunctions_    

    call Make_Table(time,deltaVirialTable,calculationDeltaVirial,cosmologyFunctions_,energyFixedAt=energyFixedAt)
    return
  end subroutine Spherical_Collapse_Dark_Energy_Virial_Density_Contrast_Tabulate

  subroutine Spherical_Collapse_Dark_Energy_Turnaround_Radius_Tabulate(time,energyFixedAt,turnaroundTable,cosmologyFunctions_)
    !% Tabulate the ratio of turnaround to virial radiii for the spherical collapse model.
    use Tables
    implicit none
    double precision                                      , intent(in   ) :: time
    integer                                               , intent(in   ) :: energyFixedAt
    class           (table1D                ), allocatable, intent(inout) :: turnaroundTable
    class           (cosmologyFunctionsClass)             , intent(inout) :: cosmologyFunctions_    

    call Make_Table(time,turnaroundTable,calculationTurnaround,cosmologyFunctions_,energyFixedAt=energyFixedAt)
    return
  end subroutine Spherical_Collapse_Dark_Energy_Turnaround_Radius_Tabulate

  subroutine Make_Table(time,deltaTable,calculationType,cosmologyFunctions_,linearGrowth_,energyFixedAt)
    !% Tabulate $\delta_\mathrm{crit}$ or $\Delta_\mathrm{vir}$ vs. time.
    use Linear_Growth
    use Root_Finder
    use Tables
    use Galacticus_Error
    use Galacticus_Display
    use Input_Parameters
    implicit none
    double precision                                       , intent(in   ) :: time
    integer                                                , intent(in   ) :: calculationType
    class           (table1D                ), allocatable , intent(inout) :: deltaTable
    class           (cosmologyFunctionsClass)              , intent(inout) :: cosmologyFunctions_    
    class           (linearGrowthClass      ), optional    , intent(inout) :: linearGrowth_    
    integer                                  , optional    , intent(in   ) :: energyFixedAt
    class           (linearGrowthClass      ), allocatable                 :: linearGrowth__
    double precision                         , parameter                   :: toleranceAbsolute              =0.0d0, toleranceRelative              =1.0d-9
    double precision                         , dimension(2)                :: timeRange
    type            (rootFinder             ), save                        :: finder                               , maximumExpansionFinder
    !$omp threadprivate(finder,maximumExpansionFinder)
    integer                                                                :: deltaTableNumberPoints               , iTime                                 , &
         &                                                                    iCount
    double precision                                                       :: aExpansionNow                        , epsilonPerturbation                   , &
         &                                                                    epsilonPerturbationMaximum           , epsilonPerturbationMinimum            , &
         &                                                                    maximumExpansionDensityContrast      , maximumExpansionExpansionFactor       , &
         &                                                                    maximumExpansionRadius               , maximumExpansionTime                  , &
         &                                                                    normalization                        , q                                     , &
         &                                                                    timeEnergyFixed                      , timeInitial                           , &
         &                                                                    y                                    , deltaTableTimeMinimum                 , &
         &                                                                    deltaTableTimeMaximum
    double complex                                                         :: a                                    , b                                     , &
         &                                                                    x
    type            (varying_string         )                              :: message
    character       (len=13                 )                              :: label

    ! Validate input.
    select case (calculationType)
    case (calculationDeltaCrit)
       if (.not.present(linearGrowth_)) call Galacticus_Error_Report('linearGrowth_ object must be provided for calcualtion of critical overdensity'//{introspection:location})
    case (calculationDeltaVirial,calculationTurnaround)
       if (.not.present(energyFixedAt)) call Galacticus_Error_Report('energyFixedAt must be provided for calcualtion of virial properties'          //{introspection:location})
    end select
    ! Find minimum and maximum times to tabulate.
    if (allocated(deltaTable)) then
       ! Use currently tabulated range as the starting point.
       deltaTableTimeMinimum=deltaTable%x(+1)
       deltaTableTimeMaximum=deltaTable%x(-1)
    else
       ! Specify an initial default range.
       deltaTableTimeMinimum= 0.1d0
       deltaTableTimeMaximum=20.0d0
    end if
    ! Expand the range to ensure the requested time is included.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)
    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(log10(deltaTableTimeMaximum/deltaTableTimeMinimum)*dble(deltaTableNPointsPerDecade))
    ! Deallocate table if currently allocated.
    if (allocated(deltaTable)) then
       call deltaTable%destroy()
       deallocate(deltaTable)
    end if
    allocate(table1DLogarithmicLinear :: deltaTable)
    select type (deltaTable)
    type is (table1DLogarithmicLinear)
       ! Create the table.
       call deltaTable%create(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints)
       ! Solve ODE to get corresponding overdensities.
       message="Solving spherical collapse model for dark energy universe for "
       write (label,'(e12.6)') deltaTableTimeMinimum
       message=message//trim(adjustl(label))//" ≤ t/Gyr ≤ "
       write (label,'(e12.6)') deltaTableTimeMaximum
       message=message//trim(adjustl(label))
       call Galacticus_Display_Indent(message,verbosity=verbosityWorking)
       iCount=0
       call Galacticus_Display_Counter(                            &
            &                                    iCount          , &
            &                          isNew    =.true.          , &
            &                          verbosity=verbosityWorking  &
            &                         )
       !$omp parallel private(aExpansionNow,epsilonPerturbationMaximum,epsilonPerturbationMinimum,epsilonPerturbation,timeInitial,timeRange,maximumExpansionTime,maximumExpansionExpansionFactor,q,y,timeEnergyFixed,a,b,x,linearGrowth__)       
       allocate(cosmologyFunctions__,mold=cosmologyFunctions_)
       !# <deepCopy source="cosmologyFunctions_" destination="cosmologyFunctions__"/>
       if (calculationType == calculationDeltaCrit) then
          allocate(linearGrowth__,mold=linearGrowth_)
          !# <deepCopy source="linearGrowth_" destination="linearGrowth__"/>
       end if
       !$omp do schedule(dynamic)
       do iTime=1,deltaTableNumberPoints
          call Galacticus_Display_Counter(                                                         &
               &                          int(100.0d0*dble(iCount-1)/dble(deltaTableNumberPoints)), &
               &                          isNew=.false.                                          , &
               &                          verbosity=verbosityWorking                               &
               &                         )
          ! Get the current expansion factor.
          aExpansionNow=cosmologyFunctions__%expansionFactor(deltaTable%x(iTime))
          ! In the case of dark energy we cannot (easily) determine the largest (i.e. least negative) value of epsilonPerturbation
          ! for which a perturbation can collapse. So, use no perturbation.
          epsilonPerturbationMaximum=  0.0d0
          ! Estimate a suitably negative minimum value for epsilon.
          epsilonPerturbationMinimum=-10.0d0
          ! Evaluate cosmological parameters at the present time.
          OmegaM               =cosmologyFunctions__%omegaMatterEpochal    (expansionFactor=aExpansionNow)
          OmegaDE              =cosmologyFunctions__%omegaDarkEnergyEpochal(expansionFactor=aExpansionNow)
          hubbleParameterInvGyr=cosmologyFunctions__%expansionRate         (                aExpansionNow)
          tNow                 =deltaTable%x(iTime)
          ! Check dark energy equation of state is within acceptable range.
          if (cosmologyFunctions__%equationOfStateDarkEnergy(time=tNow) >= -1.0d0/3.0d0) &
               & call Galacticus_Error_Report('ω<-1/3 required'//{introspection:location})
          ! Find the value of epsilon for which the perturbation just collapses at this time.
          if (.not.finder%isInitialized()) then
             call finder%rootFunction(radiusPerturbation                 )
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          end if
          epsilonPerturbation=finder%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
          select case (calculationType)
          case (calculationDeltaCrit)
             ! Critical linear overdensity.
             normalization=linearGrowth__%value(tNow,normalize=normalizeMatterDominated)/linearGrowth__%value(tNow)/aExpansionNow
             call deltaTable%populate(                                                                       &
                  &                   normalization*0.6d0*(1.0d0-OmegaM-OmegaDE-epsilonPerturbation)/OmegaM, &
                  &                   iTime                                                                  &
                  &                  )
          case (calculationDeltaVirial,calculationTurnaround)
             ! Find the epoch of maximum expansion for the perturbation.
             if (.not.maximumExpansionFinder%isInitialized()) then
                call maximumExpansionFinder%rootFunction(expansionRatePerturbation          )
                call maximumExpansionFinder%tolerance   (toleranceAbsolute,toleranceRelative)
             end if
             call maximumExpansionFinder%rangeExpand (                                                             &
                  &                                   rangeExpandDownward          =1.0d0-1.0d-2                 , &
                  &                                   rangeExpandUpward            =1.0d0+1.0d-2                 , &
                  &                                   rangeExpandType              =rangeExpandMultiplicative    , &
                  &                                   rangeUpwardLimit             =tNow                         , &
                  &                                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
                  &                                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative  &
                  &                                  )
             epsilonPerturbationShared=epsilonPerturbation
             ! Compute the corresponding time of maximum expansion.
             timeInitial                    =cosmologyFunctions__%cosmicTime(cosmologyFunctions__%expansionFactor(tNow)*expansionFactorInitialFraction)
             ! Guess that the time of maximum expansion occurred at close to half of the current time.
             timeRange                      =[0.45d0,0.55d0]*tNow
             maximumExpansionTime           =maximumExpansionFinder%find(rootRange=timeRange)
             maximumExpansionExpansionFactor=cosmologyFunctions__%expansionFactor(maximumExpansionTime)
             ! Solve the dynamics of the perturbation to find the radius at the point of maximum expansion.
             call Perturbation_Dynamics_Solver(epsilonPerturbation,maximumExpansionTime,maximumExpansionRadius)
             ! Compute the density contrast of the perturbation at maximum expansion.
             maximumExpansionDensityContrast=(maximumExpansionExpansionFactor/aExpansionNow/maximumExpansionRadius)**3
             ! Solve the cubic equation (Percival, 2005, A&A, 443, 819, eqn. 38) to give the ratio of virial to turnaround radii,
             ! x.
             q=      cosmologyFunctions__%omegaDarkEnergyEpochal(time=maximumExpansionTime) &
                  & /cosmologyFunctions__%omegaMatterEpochal    (time=maximumExpansionTime) &
                  & /maximumExpansionDensityContrast
             y=      maximumExpansionExpansionFactor**cosmologyFunctions__%exponentDarkEnergy(time=maximumExpansionTime) &
                  & /aExpansionNow                  **cosmologyFunctions__%exponentDarkEnergy(time=tNow                )
             select case (energyFixedAt)
             case (darkEnergySphericalCollapseEnergyFixedAtTurnaround   )
                timeEnergyFixed=maximumExpansionTime
             case (darkEnergySphericalCollapseEnergyFixedAtVirialization)
                timeEnergyFixed=tNow
             case default
                call Galacticus_Error_Report('unrecognized epoch'//{introspection:location})
             end select
             a=1.0d0-(1.0d0+3.0d0*cosmologyFunctions__%equationOfStateDarkEnergy(time=timeEnergyFixed))*q/2.0d0
             b=      (1.0d0+3.0d0*cosmologyFunctions__%equationOfStateDarkEnergy(time=tNow           ))*q/y
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
             case (calculationDeltaVirial)
                call deltaTable%populate(                                           &
                     &                   1.0d0/(dble(x)*maximumExpansionRadius)**3, &
                     &                   iTime                                      &
                     &                  )
             case (calculationTurnaround)
                call deltaTable%populate(                                           &
                     &                   1.0d0/ dble(x)                           , &
                     &                   iTime                                      &
                     &                  )
             end select
          end select
          !$omp atomic
          iCount=iCount+1
       end do
       !$omp end do
       deallocate(cosmologyFunctions__)
       !$omp end parallel
       call Galacticus_Display_Counter_Clear(verbosity=verbosityWorking)
       call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
    end select
    return
  end subroutine Make_Table

  double precision function radiusPerturbation(epsilonPerturbation)
    !% Return the radius of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\normalfont \ttfamily epsilonPerturbation}.
    implicit none
    double precision, intent(in   ) :: epsilonPerturbation

    call Perturbation_Dynamics_Solver(epsilonPerturbation,tNow,radiusPerturbation)
    return
  end function radiusPerturbation

  double precision function expansionRatePerturbation(time)
    !% Return the expansion rate of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\normalfont \ttfamily epsilonPerturbation}.
    implicit none
    double precision, intent(in   ) :: time

    call Perturbation_Dynamics_Solver(epsilonPerturbationShared,time,perturbationExpansionRate=expansionRatePerturbation)
    return
  end function expansionRatePerturbation

  subroutine Perturbation_Dynamics_Solver(epsilonPerturbation,time,perturbationRadius,perturbationExpansionRate)
    !% Integrate the dynamics of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\normalfont \ttfamily epsilonPerturbation}.
    use ODEIV2_Solver
    use FODEIV2
    implicit none
    double precision                                            , intent(in   )           :: epsilonPerturbation                 , time
    double precision                                            , intent(  out), optional :: perturbationExpansionRate           , perturbationRadius
    integer                             , parameter                                       :: nProperties                   =2
    double precision                    , dimension(nProperties)                          :: propertyValues
    double precision                    , parameter                                       :: odeToleranceAbsolute          =0.0d0, odeToleranceRelative            =1.0d-12
    type            (fodeiv2_system    )                                                  :: ode2System
    type            (fodeiv2_driver    )                                                  :: ode2Driver
    logical                                                                               :: odeReset
    double precision                                                                      :: expansionFactorInitial              , perturbationExpansionRateInitial        , &
         &                                                                                   perturbationOverdensityInitial      , timeInitial
    integer                                                                               :: odeStatus

    ! Specify a sufficiently early time.
    expansionFactorInitial=expansionFactorInitialFraction
    ! Find the corresponding cosmic time.
    timeInitial=cosmologyFunctions__%cosmicTime(cosmologyFunctions__%expansionFactor(time)*expansionFactorInitial)
    ! Find the overdensity of the perturbation at early time (Percival, 2005, A&A, 443, 819, eqn. 25).
    perturbationOverdensityInitial=0.6d0*(1.0d0-OmegaM-OmegaDE-epsilonPerturbation)*((3.0d0/2.0d0/OmegaM)*(hubbleParameterInvGyr&
         &*timeInitial))**(2.0d0/3.0d0)
    ! Find the perturbation radius at this early time. This is, by construction, just the initial expansion factor.
    perturbationRadiusInitial=expansionFactorInitial
    ! Find the perturbation expansion rate at early time (Percival, 2005, A&A, 443, 819, eqn. 22).
    perturbationExpansionRateInitial=hubbleParameterInvGyr*sqrt(OmegaM/expansionFactorInitial+epsilonPerturbation)
    ! Set initial conditions.
    propertyValues=[perturbationRadiusInitial,perturbationExpansionRateInitial]
    ! Evolve if the requested time is after the initial time.
    if (time > timeInitial) then
       ! Solve the ODE to find the perturbation radius at the present day.
       odeReset=.true.
       call ODEIV2_Solve(                                                  &
            &                   ode2Driver,ode2System                    , &
            &                   timeInitial,time                         , &
            &                   nProperties                              , &
            &                   propertyValues                           , &
            &                   perturbationODEs                         , &
            &                   odeToleranceAbsolute,odeToleranceRelative, &
            &                   reset=odeReset                           , &
            &                   odeStatus=odeStatus                        &
            &                  )
       call ODEIV2_Solver_Free(ode2Driver,ode2System)
       ! If the ODE solver did not succeed, it is because the perturbation collapsed to zero radius (causing a divergence). This
       ! means it collapsed prior to the current time. We extrapolate to negative radius (using the velocity at the final step) to
       ! permit our root finder to locate the point at which collapse occurs at the current time.
       if (odeStatus /= FGSL_Success) propertyValues(1)=propertyValues(1)+propertyValues(2)*(time-timeInitial)
    end if
    ! Return the computed quantities.
    if (present(perturbationRadius       )) perturbationRadius       =propertyValues(1)
    if (present(perturbationExpansionRate)) perturbationExpansionRate=propertyValues(2)
    return
  end subroutine Perturbation_Dynamics_Solver

  integer function perturbationODEs(time,y,dydt)
    implicit none
    double precision, intent(in   )               :: time
    double precision, intent(in   ), dimension(:) :: y
    double precision, intent(  out), dimension(:) :: dydt
    double precision                              :: expansionFactor

    if (y(1) <= 0.0d0) then
       dydt(1:2)=0.0d0
    else
       expansionFactor=cosmologyFunctions__%expansionFactor(time)/cosmologyFunctions__%expansionFactor(tNow)
       dydt(1)=y(2)
       dydt(2)=-0.5d0*y(1)*(hubbleParameterInvGyr**2)*(OmegaM/y(1)**3+(3.0d0*cosmologyFunctions__%equationOfStateDarkEnergy(time=time)&
            &+1.0d0)*OmegaDE*expansionFactor**cosmologyFunctions__%exponentDarkEnergy(time=time))
    end if
    ! Return success.
    perturbationODEs=FGSL_Success
    return
  end function perturbationODEs

end module Spherical_Collapse_Matter_Dark_Energy
