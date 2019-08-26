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

!% Contains a module which implements calculations of spherical top hat collapse in cosmologies containing baryons and dark matter.

module Spherical_Collapse_BDM
  !% Implements calculations of spherical top hat collapse in cosmologies containing matter and dark energy.
  use FGSL               , only : FGSL_Success
  use Cosmology_Functions, only : cosmologyFunctionsClass
  implicit none
  private
  public :: Spherical_Collapse_BDM_Critical_Overdensity_Tabulate, Spherical_Collapse_BDM_Virial_Density_Contrast_Tabulate, &
       &    Spherical_Collapse_BDM_Turnaround_Radius_Tabulate

  ! Variables to hold the tabulated critical overdensity data.
  integer                                  , parameter     :: deltaTableNPointsPerDecade    =100

  ! Variables used in root finding.
  double precision                                         :: OmegaDE                              , OmegaM                  , &
       &                                                      epsilonPerturbationShared            , hubbleParameterInvGyr   , &
       &                                                      perturbationRadiusInitial            , tNow                    , &
       &                                                      OmegaB
  logical                                                  :: baryonsCluster_
  !$omp threadprivate(OmegaDE,OmegaM,OmegaB,hubbleParameterInvGyr,tNow,epsilonPerturbationShared,perturbationRadiusInitial,baryonsCluster_)
  
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
  !#  <name>sphericalCollapseEnergyFixedAt</name>
  !#  <description>Enumeration of radii at which the energy of a spherical top-hat perturbation in a dark energy cosmology can be considered to be fixed.</description>
  !#  <encodeFunction>yes</encodeFunction>
  !#  <validator>yes</validator>
  !#  <visibility>public</visibility>
  !#  <entry label="turnaround"   />
  !#  <entry label="virialization"/>
  !# </enumeration>

contains
  
  subroutine Spherical_Collapse_BDM_Critical_Overdensity_Tabulate(time,baryonsCluster,deltaCritTable,cosmologyParameters_,cosmologyFunctions_)
    !% Tabulate the critical overdensity for collapse for the spherical collapse model.
    use Tables              , only : table1D
    use Cosmology_Parameters, only : cosmologyParametersClass
    implicit none
    double precision                                       , intent(in   ) :: time
    logical                                                , intent(in   ) :: baryonsCluster
    class           (table1D                 ), allocatable, intent(inout) :: deltaCritTable
    class           (cosmologyParametersClass)            , intent(inout) :: cosmologyParameters_    
    class           (cosmologyFunctionsClass )             , intent(inout) :: cosmologyFunctions_    

    call Make_Table(time,baryonsCluster,deltaCritTable,calculationDeltaCrit,cosmologyParameters_,cosmologyFunctions_)
    return
  end subroutine Spherical_Collapse_BDM_Critical_Overdensity_Tabulate

  subroutine Spherical_Collapse_BDM_Virial_Density_Contrast_Tabulate(time,baryonsCluster,energyFixedAt,deltaVirialTable,cosmologyParameters_,cosmologyFunctions_)
    !% Tabulate the virial density contrast for the spherical collapse model.
    use Tables              , only : table1D
    use Cosmology_Parameters, only : cosmologyParametersClass
    implicit none
    double precision                                       , intent(in   ) :: time
    logical                                                , intent(in   ) :: baryonsCluster
    integer                                                , intent(in   ) :: energyFixedAt
    class           (table1D                 ), allocatable, intent(inout) :: deltaVirialTable
    class           (cosmologyParametersClass)             , intent(inout) :: cosmologyParameters_    
    class           (cosmologyFunctionsClass )             , intent(inout) :: cosmologyFunctions_    

    call Make_Table(time,baryonsCluster,deltaVirialTable,calculationDeltaVirial,cosmologyParameters_,cosmologyFunctions_,energyFixedAt=energyFixedAt)
    return
  end subroutine Spherical_Collapse_BDM_Virial_Density_Contrast_Tabulate

  subroutine Spherical_Collapse_BDM_Turnaround_Radius_Tabulate(time,baryonsCluster,energyFixedAt,turnaroundTable,cosmologyParameters_,cosmologyFunctions_)
    !% Tabulate the ratio of turnaround to virial radiii for the spherical collapse model.
    use Tables              , only : table1D
    use Cosmology_Parameters, only : cosmologyParametersClass
    implicit none
    double precision                                       , intent(in   ) :: time
    logical                                                , intent(in   ) :: baryonsCluster
    integer                                                , intent(in   ) :: energyFixedAt
    class           (table1D                 ), allocatable, intent(inout) :: turnaroundTable
    class           (cosmologyParametersClass)             , intent(inout) :: cosmologyParameters_    
    class           (cosmologyFunctionsClass )             , intent(inout) :: cosmologyFunctions_    

    call Make_Table(time,baryonsCluster,turnaroundTable,calculationTurnaround,cosmologyParameters_,cosmologyFunctions_,energyFixedAt=energyFixedAt)
    return
  end subroutine Spherical_Collapse_BDM_Turnaround_Radius_Tabulate

  subroutine Make_Table(time,baryonsCluster,deltaTable,calculationType,cosmologyParameters_,cosmologyFunctions_,energyFixedAt)
    !% Tabulate $\delta_\mathrm{crit}$ or $\Delta_\mathrm{vir}$ vs. time.
    use Root_Finder         , only : rootFinder               , rangeExpandMultiplicative  , rangeExpandSignExpectNegative   , rangeExpandSignExpectPositive
    use Tables              , only : table1D                  , table1DLogarithmicLinear
    use Galacticus_Error    , only : Galacticus_Error_Report
    use Galacticus_Display  , only : Galacticus_Display_Indent, Galacticus_Display_Unindent, Galacticus_Display_Counter      , Galacticus_Display_Counter_Clear, &
         &                           verbosityWorking
    use Cosmology_Parameters, only : cosmologyParametersClass
    use Linear_Growth       , only : linearGrowthClass        , linearGrowthSimple         , linearGrowthNonClusteringBaryons
    use ISO_Varying_String  , only : varying_string           , assignment(=)              , operator(//)
    implicit none
    double precision                                        , intent(in   ) :: time
    logical                                                 , intent(in   ) :: baryonsCluster
    integer                                                 , intent(in   ) :: calculationType
    class           (table1D                 ), allocatable , intent(inout) :: deltaTable
    class           (cosmologyParametersClass)              , intent(inout) :: cosmologyParameters_    
    class           (cosmologyFunctionsClass )              , intent(inout) :: cosmologyFunctions_    
    integer                                   , optional    , intent(in   ) :: energyFixedAt
    class           (linearGrowthClass       ), pointer                     :: linearGrowth_
    double precision                          , parameter                   :: toleranceAbsolute              =0.0d0, toleranceRelative              =1.0d-9
    double precision                          , dimension(2)                :: timeRange
    type            (rootFinder              ), save                        :: finder                               , maximumExpansionFinder
    !$omp threadprivate(finder,maximumExpansionFinder)
    integer                                                                 :: deltaTableNumberPoints               , iTime                                 , &
         &                                                                     iCount
    double precision                                                        :: aExpansionNow                        , epsilonPerturbation                   , &
         &                                                                     epsilonPerturbationMaximum           , epsilonPerturbationMinimum            , &
         &                                                                     maximumExpansionDensityContrast      , maximumExpansionExpansionFactor       , &
         &                                                                     maximumExpansionRadius               , maximumExpansionTime                  , &
         &                                                                     normalization                        , q                                     , &
         &                                                                     timeEnergyFixed                      , timeInitial                           , &
         &                                                                     y                                    , deltaTableTimeMinimum                 , &
         &                                                                     deltaTableTimeMaximum                , r                                     , &
         &                                                                     z                                    , darkMatterFraction
    double complex                                                          :: a                                    , b                                     , &
         &                                                                     x
    type            (varying_string          )                              :: message
    character       (len=13                  )                              :: label

    ! Validate input.
    select case (calculationType)
     case (calculationDeltaVirial,calculationTurnaround)
       if (.not.present(energyFixedAt)) call Galacticus_Error_Report('energyFixedAt must be provided for calcualtion of virial properties'//{introspection:location})
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
    ! Copy baryon clustering option to module-scope.
    baryonsCluster_=baryonsCluster
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
       message="Solving spherical collapse model for baryons + dark matter + dark energy universe for "
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
       !$omp parallel private(aExpansionNow,epsilonPerturbationMaximum,epsilonPerturbationMinimum,epsilonPerturbation,timeInitial,timeRange,maximumExpansionTime,maximumExpansionExpansionFactor,q,y,timeEnergyFixed,a,b,x,linearGrowth_)       
       allocate(cosmologyFunctions__,mold=cosmologyFunctions_)
       !# <deepCopy source="cosmologyFunctions_" destination="cosmologyFunctions__"/>
       if (calculationType == calculationDeltaCrit) then
          if (baryonsCluster) then
             allocate(linearGrowthSimple :: linearGrowth_)
          else
             allocate(linearGrowthNonClusteringBaryons :: linearGrowth_)
          end if
          select type (linearGrowth_)
          type is (linearGrowthSimple              )
             linearGrowth_=linearGrowthSimple              (cosmologyParameters_,cosmologyFunctions__)
          type is (linearGrowthNonClusteringBaryons)
             linearGrowth_=linearGrowthNonClusteringBaryons(cosmologyParameters_,cosmologyFunctions__)
          end select
       end if
       !$omp do schedule(dynamic)
       do iTime=1,deltaTableNumberPoints
          call Galacticus_Display_Counter(                                                          &
               &                          int(100.0d0*dble(iCount-1)/dble(deltaTableNumberPoints)), &
               &                          isNew=.false.                                           , &
               &                          verbosity=verbosityWorking                                &
               &                         )
          ! Get the current expansion factor.
          aExpansionNow=cosmologyFunctions__%expansionFactor(deltaTable%x(iTime))
          ! Initial guess for the range of the initial perturbation amplitude. This must be positive in order for the perturbation
          ! to collapse - this sets the lower limit. Since we expect a collapsing perturbation to have linear theory amplitude of
          ! order unity at the time of collapse, and since linear perturbations grow proportional to the expansion factor in an
          ! Einstein-de Sitter universe with no baryons, we use an initial guess for the upper limit which is a multiple of our
          ! starting expansion factor.
          epsilonPerturbationMinimum=0.0d+0
          epsilonPerturbationMaximum=1.0d+1*expansionFactorInitialFraction
          ! Evaluate cosmological parameters at the present time.
          OmegaM               =+cosmologyFunctions__%omegaMatterEpochal    (expansionFactor=aExpansionNow)
          if (baryonsCluster) then
             OmegaB            =+0.0d0
          else
             OmegaB            =+OmegaM                                                                     &
                  &             *cosmologyParameters_%OmegaBaryon           (                             ) &
                  &             /cosmologyParameters_%OmegaMatter           (                             )
          end if
          OmegaDE              =+cosmologyFunctions__%omegaDarkEnergyEpochal(expansionFactor=aExpansionNow)
          hubbleParameterInvGyr=+cosmologyFunctions__%expansionRate         (                aExpansionNow)
          tNow                 =+deltaTable%x(iTime)
          ! Check dark energy equation of state is within acceptable range.
          if (cosmologyFunctions__%equationOfStateDarkEnergy(time=tNow) >= -1.0d0/3.0d0) &
               & call Galacticus_Error_Report('ω<-⅓ required'//{introspection:location})
          ! Find the value of epsilon for which the perturbation just collapses at this time.
          if (.not.finder%isInitialized()) then
             call finder%rootFunction(radiusPerturbation                 )
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
             call finder%rangeExpand(rangeExpandUpward=2.0d0,rangeExpandType=rangeExpandMultiplicative)
          end if
          epsilonPerturbation=finder%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
          select case (calculationType)
          case (calculationDeltaCrit)
             ! Critical linear overdensity.
             normalization=linearGrowth_%value(tNow)/linearGrowth_%value(cosmologyFunctions__%cosmicTime(expansionFactorInitialFraction*cosmologyFunctions__%expansionFactor(tNow)))
             call deltaTable%populate(                                   &
                  &                   normalization*epsilonPerturbation, &
                  &                   iTime                              &
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
             call Perturbation_Dynamics_Solver(epsilonPerturbation,maximumExpansionTime,tNow,maximumExpansionRadius)
             ! Compute the density contrast of the perturbation at maximum expansion.
             maximumExpansionDensityContrast=(maximumExpansionExpansionFactor/aExpansionNow/maximumExpansionRadius)**3
             ! Solve the cubic equation (Percival, 2005, A&A, 443, 819, eqn. 38) to give the ratio of virial to turnaround radii,
             ! x.
             select case (energyFixedAt)
             case (sphericalCollapseEnergyFixedAtTurnaround   )
                timeEnergyFixed=maximumExpansionTime
             case (sphericalCollapseEnergyFixedAtVirialization)
                timeEnergyFixed=tNow
             case default
                call Galacticus_Error_Report('unrecognized epoch'//{introspection:location})
             end select
             if (baryonsCluster) then
                q=      cosmologyFunctions__%omegaDarkEnergyEpochal(time=maximumExpansionTime) &
                     & /cosmologyFunctions__%omegaMatterEpochal    (time=maximumExpansionTime) &
                     & /maximumExpansionDensityContrast
                y=      maximumExpansionExpansionFactor**cosmologyFunctions__%exponentDarkEnergy(time=maximumExpansionTime) &
                     & /aExpansionNow                  **cosmologyFunctions__%exponentDarkEnergy(time=tNow                )
                a=1.0d0-(1.0d0+3.0d0*cosmologyFunctions__%equationOfStateDarkEnergy(time=timeEnergyFixed))*q/2.0d0
                b=      (1.0d0+3.0d0*cosmologyFunctions__%equationOfStateDarkEnergy(time=tNow           ))*q/y
             else
                darkMatterFraction=(cosmologyParameters_%OmegaMatter()-cosmologyParameters_%OmegaBaryon())/cosmologyParameters_%OmegaMatter()
                q=      cosmologyFunctions__%omegaDarkEnergyEpochal(time=maximumExpansionTime) &
                     & /cosmologyFunctions__%omegaMatterEpochal    (time=maximumExpansionTime) &
                     & /darkMatterFraction                                                     &
                     & /maximumExpansionDensityContrast
                y=      maximumExpansionExpansionFactor**cosmologyFunctions__%exponentDarkEnergy(time=maximumExpansionTime) &
                     & /aExpansionNow                  **cosmologyFunctions__%exponentDarkEnergy(time=tNow                )
                
                r=     +  cosmologyParameters_%OmegaBaryon() &
                     & /(                                    &
                     &   +cosmologyParameters_%OmegaMatter() &
                     &   -cosmologyParameters_%OmegaBaryon() & 
                     &  )                                    &
                     & /maximumExpansionDensityContrast
                z=     +(                                    &
                     &   +maximumExpansionExpansionFactor    &
                     &   /aExpansionNow                      &
                     &  )**3
                a=     +1.0d0+r        -(1.0d0+3.0d0*cosmologyFunctions__%equationOfStateDarkEnergy(time=timeEnergyFixed))*q/2.0d0
                b=           -r/z/2.0d0+(1.0d0+3.0d0*cosmologyFunctions__%equationOfStateDarkEnergy(time=tNow           ))*q/y
             end if
             x=                                                                                                        &
                  & +(0.0d0,0.5d0)*sqrt(3.0d0)                                                                         &
                  & *(                                                                                                 &
                  &   1.0d0/b*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(+1.0d0/3.0d0)/ 6.0d0  &
                  &  +2.0d0*a*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(-1.0d0/3.0d0)         &
                  & )                                                                                                  &
                  & -(1.0d0/b*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(+1.0d0/3.0d0)/12.0d0) &
                  & +(      a*((54.0d0+6.0d0*sqrt(3.0d0)*sqrt((16.0d0*a**3+27.0d0*b)/b))*b**2)**(-1.0d0/3.0d0)       )
             select case (calculationType)
             case (calculationDeltaVirial)
                ! The density contrast calculated as Δ=1/(x Rmax)³ is Δ=ρvir/⟨ρDM⟩ - i.e. the density of the virialized dark
                ! matter perturbation relative to the mean density of dark matter. However, what we want (for the definition used
                ! by Galacticus) is the density of the perturbation relative to the total mean density. So we perform that
                ! conversion here.
                call deltaTable%populate(                                      &
                     &                   +(                                    &
                     &                     +cosmologyParameters_%OmegaMatter() &
                     &                     -cosmologyParameters_%OmegaBaryon() &
                     &                   )                                     &
                     &                   /  cosmologyParameters_%OmegaMatter() &
                     &                   /(dble(x)*maximumExpansionRadius)**3, &
                     &                   iTime                                 &
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

    call Perturbation_Dynamics_Solver(epsilonPerturbation,tNow,tNow,radiusPerturbation)
    return
  end function radiusPerturbation

  double precision function expansionRatePerturbation(time)
    !% Return the expansion rate of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\normalfont \ttfamily epsilonPerturbation}.
    implicit none
    double precision, intent(in   ) :: time

    call Perturbation_Dynamics_Solver(epsilonPerturbationShared,time,tNow,perturbationExpansionRate=expansionRatePerturbation)
    return
  end function expansionRatePerturbation

  subroutine Perturbation_Dynamics_Solver(perturbationOverdensityInitial,time,tNow,perturbationRadius,perturbationExpansionRate)
    !% Integrate the dynamics of a spherical top-hat perturbation in a dark energy universe given an initial perturbation
    !% amplitude {\normalfont \ttfamily epsilonPerturbation}.
    use ODEIV2_Solver
    use FODEIV2
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    double precision                                            , intent(in   ) :: perturbationOverdensityInitial      , time                                    , &
         &                                                                         tNow
    double precision                    , optional              , intent(  out) :: perturbationExpansionRate           , perturbationRadius
    integer                             , parameter                             :: nProperties                   =2
    double precision                    , dimension(nProperties)                :: propertyValues
    double precision                    , parameter                             :: odeToleranceAbsolute          =0.0d0, odeToleranceRelative            =1.0d-12
    type            (fodeiv2_system    )                                        :: ode2System
    type            (fodeiv2_driver    )                                        :: ode2Driver
    logical                                                                     :: odeReset
    double precision                                                            :: expansionFactorInitial              , perturbationExpansionRateInitial        , &
         &                                                                         timeInitial                         , perturbationGrowthRate                  , &
         &                                                                         exponent
    integer                                                                     :: odeStatus

    ! Validate the perturbation overdensity.
    if (perturbationOverdensityInitial < 0.0d+0) call Galacticus_Error_Report('initial overdensity of perturbation should be non-negative'//{introspection:location})
    if (perturbationOverdensityInitial > 1.0d-3) call Galacticus_Error_Report('initial overdensity of perturbation should be small'       //{introspection:location})
    ! Specify a sufficiently early time.
    expansionFactorInitial=expansionFactorInitialFraction
    ! Find the corresponding cosmic time.
    timeInitial=cosmologyFunctions__%cosmicTime(expansionFactorInitial*cosmologyFunctions__%expansionFactor(tNow))
    ! Determine the initial growth rate of the overdensity assuming the matter-dominated phase growth factor is proportional to the expansion factor.
    if (baryonsCluster_) then
       ! Baryons are assumed to cluster, so the usual matter dominated solution of D(a)=a applies.
       perturbationGrowthRate=perturbationOverdensityInitial*hubbleParameterInvGyr*sqrt(OmegaM/expansionFactorInitial)
    else
       ! Baryons do not cluster, the growth factor is D(a)=a^{3p/2} with p=[sqrt(1+24fₓ)-1]/6
       exponent=0.25d0*(sqrt(1.0d0+24.0d0*(OmegaM-OmegaB)/OmegaM)-1.0d0)
       perturbationGrowthRate=perturbationOverdensityInitial*exponent*hubbleParameterInvGyr*sqrt(OmegaM/expansionFactorInitial)
    end if
    ! Find the initial radius of the perturbation.
    perturbationRadiusInitial=expansionFactorInitial*(1.0d0-perturbationOverdensityInitial/3.0d0)
    ! Find the initial growth rate of the perturbation radius.
    perturbationExpansionRateInitial=perturbationRadiusInitial*hubbleParameterInvGyr*sqrt(OmegaM/expansionFactorInitial**3)*(1.0d0-perturbationOverdensityInitial/3.0d0/(1.0d0+perturbationOverdensityInitial))
    ! Set initial conditions.
    propertyValues=[perturbationRadiusInitial,perturbationExpansionRateInitial]
    ! Evolve if the requested time is after the initial time.
    if (time > timeInitial) then
       ! Solve the ODE to find the perturbation radius at the present day.
       odeReset=.true.
       call ODEIV2_Solve(                                &
            &                      ode2Driver          , &
            &                      ode2System          , &
            &                      timeInitial         , &
            &                      time                , &
            &                      nProperties         , &
            &                      propertyValues      , &
            &                      perturbationODEs    , &
            &                      odeToleranceAbsolute, &
            &                      odeToleranceRelative, &
            &            reset    =odeReset            , &
            &            odeStatus=odeStatus             &
            &           )
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
    !% Differential equations describing the collapse of a spherical perturbation in a universe with non-clustering baryons.
    implicit none
    double precision, intent(in   )               :: time
    double precision, intent(in   ), dimension(:) :: y
    double precision, intent(  out), dimension(:) :: dydt
    double precision                              :: expansionFactor

    if (y(1) <= 0.0d0) then
       dydt(1:2)=0.0d0
    else
       expansionFactor  =+cosmologyFunctions__%expansionFactor(time) &
            &            /cosmologyFunctions__%expansionFactor(tNow)
       dydt          (1)=+y(2)
       dydt          (2)=-0.5d0*y(1)*(hubbleParameterInvGyr**2)                                                                                                                             &
            &            *(                                                                                                                                                                 &
            &               +(OmegaM-OmegaB)                                                                        /y              (1)**3                                                  &
            &               +        OmegaB                                                                         /expansionFactor   **3                                                  &
            &               + OmegaDE       *(3.0d0*cosmologyFunctions__%equationOfStateDarkEnergy(time=time)+1.0d0)*expansionFactor   **cosmologyFunctions__%exponentDarkEnergy(time=time) &
            &              )
    end if
    ! Return success.
    perturbationODEs=FGSL_Success
    return
  end function perturbationODEs

end module Spherical_Collapse_BDM
