!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use FGSL
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  implicit none
  private
  public :: Spherical_Collapse_Dark_Energy_Delta_Critical_Initialize, Spherical_Collapse_Dark_Energy_Critical_Overdensity,&
       & Spherical_Collapse_Dark_Energy_Delta_Virial_Initialize, Spherical_Collapse_Dark_Energy_Virial_Density_Contrast,&
       & Spherical_Collapse_Matter_Dark_Energy_State_Store, Spherical_Collapse_Matter_Dark_Energy_State_Retrieve

  ! Variables to hold the tabulated critical overdensity data.
  double precision                            :: deltaTableTimeMaximum                                    =20.0d0, deltaTableTimeMinimum =1.0d0
  integer                         , parameter :: deltaTableNPointsPerDecade                               =100

  ! Variables used in root finding.
  double precision                            :: OmegaDE                                                         , OmegaM                       , &
       &                                         epsilonPerturbationShared                                       , hubbleParameterInvGyr        , &
       &                                         perturbationRadiusInitial                                       , tNow

  ! Fraction of current expansion factor to use as initial time in perturbation dynamics solver.
  double precision                , parameter :: expansionFactorInitialFraction                           =1.0d-6

  ! Calculation types.
  integer                         , parameter :: calculationDeltaCrit                                     =0     , calculationDeltaVirial=1

  ! Parameter controlling the epoch at which the energy of a perturbation is fixed when computing virial overdensities.
  type            (varying_string)            :: virialDensityContrastSphericalTopHatDarkEnergyFixEnergyAt

contains

  !# <criticalOverdensityMethod>
  !#  <unitName>Spherical_Collapse_Dark_Energy_Delta_Critical_Initialize</unitName>
  !# </criticalOverdensityMethod>
  subroutine Spherical_Collapse_Dark_Energy_Delta_Critical_Initialize(criticalOverdensityMethod,Critical_Overdensity_Tabulate)
    !% Initializes the $\delta_{\rm crit}$ calculation for the spherical collapse module.
    implicit none
    type     (varying_string                                     ), intent(in   )          :: criticalOverdensityMethod
    procedure(Spherical_Collapse_Dark_Energy_Critical_Overdensity), intent(inout), pointer :: Critical_Overdensity_Tabulate

    if (criticalOverdensityMethod == 'sphericalTopHatDarkEnergy') Critical_Overdensity_Tabulate => Spherical_Collapse_Dark_Energy_Critical_Overdensity
    return
  end subroutine Spherical_Collapse_Dark_Energy_Delta_Critical_Initialize

  !# <virialDensityContrastMethod>
  !#  <unitName>Spherical_Collapse_Dark_Energy_Delta_Virial_Initialize</unitName>
  !# </virialDensityContrastMethod>
  subroutine Spherical_Collapse_Dark_Energy_Delta_Virial_Initialize(virialDensityContrastMethod,Virial_Density_Contrast_Tabulate)
    !% Initializes the $\Delta_{\rm vir}$ calculation for the spherical collapse module.
    use Input_Parameters
    type     (varying_string                                        ), intent(in   )          :: virialDensityContrastMethod
    procedure(Spherical_Collapse_Dark_Energy_Virial_Density_Contrast), intent(inout), pointer :: Virial_Density_Contrast_Tabulate

    if (virialDensityContrastMethod == 'sphericalTopHatDarkEnergy') then
       Virial_Density_Contrast_Tabulate => Spherical_Collapse_Dark_Energy_Virial_Density_Contrast
       ! Read parameters controlling the calculation.
       !@ <inputParameter>
       !@   <name>virialDensityContrastSphericalTopHatDarkEnergyFixEnergyAt</name>
       !@   <defaultValue>turnaround</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Selects the epoch at which the energy of a spherical top hat perturbation in a dark energy cosmology should be
       !@    ``fixed'' for the purposes of computing virial density contrasts. (See the discussion in
       !@    \citealt{percival_cosmological_2005}; \S8.)
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('virialDensityContrastSphericalTopHatDarkEnergyFixEnergyAt',virialDensityContrastSphericalTopHatDarkEnergyFixEnergyAt,defaultValue='turnaround')
    end if
    return
  end subroutine Spherical_Collapse_Dark_Energy_Delta_Virial_Initialize

  subroutine Spherical_Collapse_Dark_Energy_Critical_Overdensity(time,deltaCritTable)
    !% Tabulate the critical overdensity for collapse for the spherical collapse model.
    use Tables
    implicit none
    double precision                      , intent(in   ) :: time
    class           (table1D), allocatable, intent(inout) :: deltaCritTable

    !$omp critical(Spherical_Collapse_Make_Table)
    call Make_Table(time,deltaCritTable,calculationDeltaCrit)
    !$omp end critical(Spherical_Collapse_Make_Table)

    return
  end subroutine Spherical_Collapse_Dark_Energy_Critical_Overdensity

  subroutine Spherical_Collapse_Dark_Energy_Virial_Density_Contrast(time,deltaVirialTable)
    !% Tabulate the virial density contrast for the spherical collapse model.
    use Tables
    implicit none
    double precision                      , intent(in   ) :: time
    class           (table1D), allocatable, intent(inout) :: deltaVirialTable

    !$omp critical(Spherical_Collapse_Make_Table)
    call Make_Table(time,deltaVirialTable,calculationDeltaVirial)
    !$omp end critical(Spherical_Collapse_Make_Table)
    return
  end subroutine Spherical_Collapse_Dark_Energy_Virial_Density_Contrast

  subroutine Make_Table(time,deltaTable,calculationType)
    !% Tabulate $\delta_{\rm crit}$ or $\Delta_{\rm vir}$ vs. time.
    use Linear_Growth
    use Cosmology_Functions
    use Root_Finder
    use Tables
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    double precision                             , intent(in   ) :: time
    integer                                      , intent(in   ) :: calculationType
    class           (table1D       ), allocatable, intent(inout) :: deltaTable
    double precision                , parameter                  :: toleranceAbsolute              =0.0d0, toleranceRelative              =1.0d-9
    type            (rootFinder    ), save                       :: finder                               , maximumExpansionFinder
    !$omp threadprivate(finder,maximumExpansionFinder)
    integer                                                      :: deltaTableNumberPoints               , iTime
    double precision                                             :: aExpansionNow                        , epsilonPerturbation                    , &
         &                                                          epsilonPerturbationMaximum           , epsilonPerturbationMinimum             , &
         &                                                          maximumExpansionDensityContrast      , maximumExpansionExpansionFactor        , &
         &                                                          maximumExpansionRadius               , maximumExpansionTime                   , &
         &                                                          y                                    , normalization                          , &
         &                                                          q                                    , timeEnergyFixed                        , &
         &                                                          timeInitial
    double complex :: a,b,x
    type     (varying_string) :: message
    character(len=7         ) :: label

    ! Find minimum and maximum times to tabulate.
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
       write (label,'(f7.3)') deltaTableTimeMinimum
       message=message//trim(adjustl(label))//" ≤ t/Gyr ≤ "
       write (label,'(f7.3)') deltaTableTimeMaximum
       message=message//trim(adjustl(label))
       call Galacticus_Display_Indent(message,verbosity=verbosityWorking)
       do iTime=1,deltaTableNumberPoints
          call Galacticus_Display_Counter(                                                         &
               &                          int(100.0d0*dble(iTime-1)/dble(deltaTableNumberPoints)), &
               &                          isNew=(iTime==1)                                       , &
               &                          verbosity=verbosityWorking                               &
               &                         )
          ! Get the current expansion factor.
          aExpansionNow=Expansion_Factor(deltaTable%x(iTime))
          ! In the case of dark energy we cannot (easily) determine the largest (i.e. least negative) value of epsilonPerturbation
          ! for which a perturbation can collapse. So, use no perturbation.
          epsilonPerturbationMaximum=  0.0d0
          ! Estimate a suitably negative minimum value for epsilon.
          epsilonPerturbationMinimum=-10.0d0
          ! Evaluate cosmological parameters at the present time.
          OmegaM               =Omega_Matter_Total(aExpansion=aExpansionNow)
          OmegaDE              =Omega_Dark_Energy (aExpansion=aExpansionNow)
          hubbleParameterInvGyr=Expansion_Rate    (           aExpansionNow)
          tNow                 =deltaTable%x(iTime)
          ! Check dark energy equation of state is within acceptable range.
          if (Cosmology_Dark_Energy_Equation_Of_State(time=tNow) >= -1.0d0/3.0d0) &
               & call Galacticus_Error_Report('Make_Table','ω<-1/3 required')
          ! Find the value of epsilon for which the perturbation just collapses at this time.
          if (.not.finder%isInitialized()) then
             call finder%rootFunction(radiusPerturbation                 )
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          end if
          epsilonPerturbation=finder%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
          ! Compute the corresponding critical overdensity.
          normalization=Linear_Growth_Factor(tNow,normalize=normalizeMatterDominated)/Linear_Growth_Factor(tNow)/aExpansionNow
          select case (calculationType)
          case (calculationDeltaCrit)
             ! Critical linear overdensity.
             call deltaTable%populate(                                                                       &
                  &                   normalization*0.6d0*(1.0d0-OmegaM-OmegaDE-epsilonPerturbation)/OmegaM, &
                  &                   iTime                                                                  &
                  &                  )
          case (calculationDeltaVirial)
             ! Find the epoch of maximum expansion for the perturbation.
             if (.not.maximumExpansionFinder%isInitialized()) then
                call maximumExpansionFinder%rootFunction(expansionRatePerturbation          )
                call maximumExpansionFinder%tolerance   (toleranceAbsolute,toleranceRelative)
             end if
             epsilonPerturbationShared=epsilonPerturbation
             ! Compute the corresponding time of maximum expansion.
             timeInitial                    =Cosmology_Age(Expansion_Factor(tNow)*expansionFactorInitialFraction)
             maximumExpansionTime           =maximumExpansionFinder%find(rootRange=[timeInitial,tNow])
             maximumExpansionExpansionFactor=Expansion_Factor(maximumExpansionTime)
             ! Solve the dynamics of the perturbation to find the radius at the point of maximum expansion.
             call Perturbation_Dynamics_Solver(epsilonPerturbation,maximumExpansionTime,maximumExpansionRadius)
             ! Compute the density contrast of the perturbation at maximum expansion.
             maximumExpansionDensityContrast=(maximumExpansionExpansionFactor/aExpansionNow/maximumExpansionRadius)**3
             ! Solve the cubic equation (Percival, 2005, A&A, 443, 819, eqn. 38) to give the ratio of virial to turnaround radii,
             ! x.
             q=      Omega_Dark_Energy (tCosmological=maximumExpansionTime) &
                  & /Omega_Matter_Total(tCosmological=maximumExpansionTime) &
                  & /maximumExpansionDensityContrast
             y=      maximumExpansionExpansionFactor**Cosmology_Dark_Energy_Exponent(time=maximumExpansionTime) &
                  & /aExpansionNow                  **Cosmology_Dark_Energy_Exponent(time=tNow                )
             select case (char(virialDensityContrastSphericalTopHatDarkEnergyFixEnergyAt))
             case ('turnaround'   )
                timeEnergyFixed=maximumExpansionTime
             case ('virialization')
                timeEnergyFixed=tNow
             case default
                call Galacticus_Error_Report('Make_Table','unrecognized epoch')
             end select
             a=1.0d0-(1.0d0+3.0d0*Cosmology_Dark_Energy_Equation_Of_State(time=timeEnergyFixed))*q/2.0d0
             b=(1.0d0+3.0d0*Cosmology_Dark_Energy_Equation_Of_State(time=tNow))*q/y
             ! Check for special cases.
             if (q == 0.0d0) then
                ! No dark energy, the ratio of radii is always 1/2.
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
             call deltaTable%populate(                                           &
                  &                   1.0d0/(dble(x)*maximumExpansionRadius)**3, &
                  &                   iTime                                      &
                  &                  )
          end select
       end do
       call Galacticus_Display_Counter_Clear(verbosity=verbosityWorking)
       call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
    end select
    return
  end subroutine Make_Table

  double precision function radiusPerturbation(epsilonPerturbation)
    !% Return the radius of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\tt epsilonPerturbation}.
    implicit none
    double precision, intent(in   ) :: epsilonPerturbation

    call Perturbation_Dynamics_Solver(epsilonPerturbation,tNow,radiusPerturbation)
    return
  end function radiusPerturbation

  double precision function expansionRatePerturbation(time)
    !% Return the expansion rate of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\tt epsilonPerturbation}.
    implicit none
    double precision, intent(in   ) :: time

    call Perturbation_Dynamics_Solver(epsilonPerturbationShared,time,perturbationExpansionRate=expansionRatePerturbation)
    return
  end function expansionRatePerturbation

  subroutine Perturbation_Dynamics_Solver(epsilonPerturbation,time,perturbationRadius,perturbationExpansionRate)
    !% Integrate the dynamics of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\tt epsilonPerturbation}.
    use ODEIV2_Solver
    use FODEIV2
    use Cosmology_Functions
    implicit none
    double precision                                            , intent(in   )           :: epsilonPerturbation                 , time
    double precision                                            , intent(  out), optional :: perturbationExpansionRate           , perturbationRadius
    integer                             , parameter                                       :: nProperties                   =2
    double precision                    , dimension(nProperties)                          :: propertyValues
    double precision                    , parameter                                       :: odeToleranceAbsolute          =0.0d0, odeToleranceRelative            =1.0d-12
    type            (fodeiv2_system    )                                                  :: ode2System
    type            (fodeiv2_driver    )                                                  :: ode2Driver
    type            (c_ptr             )                                                  :: parameterPointer
    logical                                                                               :: odeReset
    double precision                                                                      :: expansionFactorInitial              , perturbationExpansionRateInitial         , &
         &                                                                                   perturbationOverdensityInitial      , timeInitial
    integer                                                                               :: odeStatus

    ! Specify a sufficiently early time.
    expansionFactorInitial=expansionFactorInitialFraction
    ! Find the corresponding cosmic time.
    timeInitial=Cosmology_Age(Expansion_Factor(time)*expansionFactorInitial)
    ! Find the overdensity of the perturbation at early time (Percival, 2005, A&A, 443, 819, eqn. 25).
    perturbationOverdensityInitial=0.6d0*(1.0d0-OmegaM-OmegaDE-epsilonPerturbation)*((3.0d0/2.0d0/OmegaM)*(hubbleParameterInvGyr&
         &*timeInitial))**(2.0d0/3.0d0)
    ! Find the perturbation radius at this early time. This is, by construction, just the initial expansion factor.
    perturbationRadiusInitial=expansionFactorInitial
    ! Find the perturbation expansion rate at early time (Percival, 2005, A&A, 443, 819, eqn. 22).
    perturbationExpansionRateInitial=hubbleParameterInvGyr*sqrt(OmegaM/expansionFactorInitial+epsilonPerturbation)
    ! Set initial conditions.
    propertyValues=[perturbationRadiusInitial,perturbationExpansionRateInitial]
    ! Evovle if the requested time is after the initial time.
    if (time > timeInitial) then
       ! Solve the ODE to find the perturbation radius at the present day.
       odeReset=.true.
       call ODEIV2_Solve(                                                  &
            &                   ode2Driver,ode2System                    , &
            &                   timeInitial,time                         , &
            &                   nProperties                              , &
            &                   propertyValues                           , &
            &                   perturbationODEs                         , &
            &                   parameterPointer                         , &
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

  function perturbationODEs(time,y,dydt,parameterPointer) bind(c)
    use Cosmology_Functions
    implicit none
    integer         (kind=c_int   )                       :: perturbationODEs
    real            (kind=c_double)               , value :: time
    real            (kind=c_double), intent(in   )        :: y               (2)
    real            (kind=c_double), intent(  out)        :: dydt            (2)
    type            (c_ptr        )               , value :: parameterPointer
    double precision                                      :: expansionFactor

    if (y(1) <= 0.0d0) then
       dydt=0.0d0
    else
       expansionFactor=Expansion_Factor(time)/Expansion_Factor(tNow)
       dydt(1)=y(2)
       dydt(2)=-0.5d0*y(1)*(hubbleParameterInvGyr**2)*(OmegaM/y(1)**3+(3.0d0*Cosmology_Dark_Energy_Equation_Of_State(time=time)&
            &+1.0d0)*OmegaDE*expansionFactor**Cosmology_Dark_Energy_Exponent(time=time))
    end if
    ! Return success.
    perturbationODEs=FGSL_Success
    return
  end function perturbationODEs

  !# <galacticusStateStoreTask>
  !#  <unitName>Spherical_Collapse_Matter_Dark_Energy_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Spherical_Collapse_Matter_Dark_Energy_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    write (stateFile) deltaTableTimeMinimum,deltaTableTimeMaximum
    return
  end subroutine Spherical_Collapse_Matter_Dark_Energy_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Spherical_Collapse_Matter_Dark_Energy_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Spherical_Collapse_Matter_Dark_Energy_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    read (stateFile) deltaTableTimeMinimum,deltaTableTimeMaximum
    return
  end subroutine Spherical_Collapse_Matter_Dark_Energy_State_Retrieve

end module Spherical_Collapse_Matter_Dark_Energy
