!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements calculations of spherical top hat collapse in cosmologies containing matter and a
!% cosmological constant.

module Spherical_Collapse_Matter_Lambda
  use            :: FGSL         , only : fgsl_function, fgsl_integration_workspace
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Spherical_Collapse_Matter_Lambda_Critical_Overdensity_Tabulate, Spherical_Collape_Matter_Lambda_Delta_Virial_Tabulate, &
       &    Spherical_Collapse_Matter_Lambda_Nonlinear_Mapping

  ! Variables to hold the tabulated critical overdensity data.
  integer         , parameter :: deltaTableNPointsPerDecade=1000

  ! Variables used in root finding.
  double precision            :: OmegaDE                          , OmegaM                      , &
       &                         epsilonPerturbationShared        , hubbleParameterInvGyr       , &
       &                         tNow                             , timeTarget, radiusMaximum
  !$omp threadprivate(OmegaDE,OmegaM,epsilonPerturbationShared,hubbleParameterInvGyr,tNow,timeTarget,radiusMaximum)
  
  ! Calculation types.
  integer         , parameter :: calculationDeltaCrit          =0     , calculationDeltaVirial=1

contains

  subroutine Spherical_Collapse_Matter_Lambda_Critical_Overdensity_Tabulate(time,deltaCritTable,cosmologyFunctions_,linearGrowth_)
    !% Tabulate the critical overdensity for collapse for the spherical collapse model.
    use Tables
    use Cosmology_Functions
    use Linear_Growth
    implicit none
    double precision                                      , intent(in   ) :: time
    class           (table1D                ), allocatable, intent(inout) :: deltaCritTable
    class           (cosmologyFunctionsClass)             , intent(inout) :: cosmologyFunctions_    
    class           (linearGrowthClass      )             , intent(inout) :: linearGrowth_    

    call Make_Table(time,deltaCritTable,calculationDeltaCrit,cosmologyFunctions_,linearGrowth_)
    return
  end subroutine Spherical_Collapse_Matter_Lambda_Critical_Overdensity_Tabulate

  subroutine Spherical_Collape_Matter_Lambda_Delta_Virial_Tabulate(time,deltaVirialTable,cosmologyFunctions_)
    !% Tabulate the virial density contrast for the spherical collapse model.
    use Tables
    use Cosmology_Functions
    implicit none
    double precision                                      , intent(in   ) :: time
    class           (table1D                ), allocatable, intent(inout) :: deltaVirialTable
    class           (cosmologyFunctionsClass)             , intent(inout) :: cosmologyFunctions_    

    call Make_Table(time,deltaVirialTable,calculationDeltaVirial,cosmologyFunctions_)
    return
  end subroutine Spherical_Collape_Matter_Lambda_Delta_Virial_Tabulate

  subroutine Make_Table(time,deltaTable,calculationType,cosmologyFunctions_,linearGrowth_)
    !% Tabulate $\delta_\mathrm{crit}$ or $\Delta_\mathrm{vir}$ vs. time.
    use Linear_Growth
    use Cosmology_Functions
    use Root_Finder
    use Tables
    use Kind_Numbers
    use Galacticus_Error
    implicit none
    double precision                                      , intent(in   ) :: time
    integer                                               , intent(in   ) :: calculationType
    class           (table1D                ), allocatable, intent(inout) :: deltaTable
    class           (cosmologyFunctionsClass)             , intent(inout) :: cosmologyFunctions_    
    class           (linearGrowthClass      ), optional   , intent(inout) :: linearGrowth_    
    double precision                         , parameter                  :: toleranceAbsolute         =0.0d+0, toleranceRelative         =1.0d-9
    type            (rootFinder             ), save                       :: finder
    !$omp threadprivate(finder)
    integer                                                               :: deltaTableNumberPoints           , iTime
    double precision                                                      :: aExpansionNow                    , epsilonPerturbation              , &
         &                                                                   epsilonPerturbationMaximum       , epsilonPerturbationMinimum       , &
         &                                                                   eta                              , normalization                    , &
         &                                                                   radiiRatio                       , radiusMaximum                    , &
         &                                                                   timeBigCrunch                    , deltaTableTimeMinimum            , &
         &                                                                   deltaTableTimeMaximum
    double complex                                                        :: a                                , b                                , &
         &                                                                   c                                , d                                , &
         &                                                                   Delta

    ! Validate input.
    if (calculationType == calculationDeltaCrit .and. .not.present(linearGrowth_)) call Galacticus_Error_Report('linearGrowth_ object must be provided for calcualtion of critical overdensity'//{introspection:location})
    ! Find minimum and maximum times to tabulate.
    if (allocated(deltaTable)) then
       ! Use currently tabulated range as the starting point.
       deltaTableTimeMinimum=deltaTable%x(+1)
       deltaTableTimeMaximum=deltaTable%x(-1)
    else
       ! Specify an initial default range.
       deltaTableTimeMinimum= 1.0d0
       deltaTableTimeMaximum=20.0d0
    end if
    ! Expand the range to ensure the requested time is included.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)
    timeBigCrunch=cosmologyFunctions_%timeBigCrunch()
    if (timeBigCrunch > 0.0d0) then
       ! A Big Crunch exists - avoid attempting to tabulate times beyond this epoch.
       if (deltaTableTimeMinimum > timeBigCrunch) deltaTableTimeMinimum= 0.5d0                                *timeBigCrunch
       if (deltaTableTimeMaximum > timeBigCrunch) deltaTableTimeMaximum=(1.0d0-timeToleranceRelativeBigCrunch)*timeBigCrunch
    end if
    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(log10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))
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
       do iTime=1,deltaTableNumberPoints
          tNow=deltaTable%x(iTime)

          ! Get the current expansion factor.
          aExpansionNow=cosmologyFunctions_%expansionFactor(tNow)
          ! Determine the largest (i.e. least negative) value of epsilonPerturbation for which a perturbation can collapse.
          if (cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=aExpansionNow)>0.0d0) then
             epsilonPerturbationMaximum=-(27.0d0*cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=aExpansionNow)*(cosmologyFunctions_%omegaMatterEpochal(expansionFactor=aExpansionNow)&
                  & **2)/4.0d0)**(1.0d0/3.0d0)
          else
             epsilonPerturbationMaximum=-1.0d-6
          end if

          ! Estimate a suitably negative minimum value for epsilon.
          epsilonPerturbationMinimum=-10.0d0

          OmegaM               =    cosmologyFunctions_%omegaMatterEpochal    (expansionFactor=aExpansionNow)
          OmegaDE              =    cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=aExpansionNow)
          hubbleParameterInvGyr=abs(cosmologyFunctions_%expansionRate         (                aExpansionNow))

          ! Find the value of epsilon for which the perturbation just collapses at this time.
          if (.not.finder%isInitialized()) then
             call finder%rootFunction(collapseRoot                       )
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
             call finder%rangeExpand (                                                             &
                  &                   rangeExpandUpward            =0.5d0                        , &
                  &                   rangeExpandDownward          =2.0d0                        , &
                  &                   rangeExpandType              =rangeExpandMultiplicative    , &
                  &                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
                  &                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
                  &                  )
          end if
          epsilonPerturbation=finder%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
          ! Compute the corresponding critical overdensity.
          select case (calculationType)
          case (calculationDeltaCrit)
             ! Critical linear overdensity.
             normalization=linearGrowth_%value(tNow,normalize=normalizeMatterDominated)/linearGrowth_%value(tNow)/aExpansionNow
             call deltaTable%populate(                                                                       &
                  &                   normalization*0.6d0*(1.0d0-OmegaM-OmegaDE-epsilonPerturbation)/OmegaM, &
                  &                   iTime                                                                  &
                  &                  )
             ! Check for non-monotonic decline.
             if (iTime > 1) then
                if (deltaTable%y(iTime) >= deltaTable%y(iTime-1)) call Galacticus_Error_Report('accuracy lost in tablulation of critical overdensity (usually results for computing critical overdensity for very large cosmic times)'//{introspection:location})
             end if
         case (calculationDeltaVirial)
             ! Compute the maximum radius of the perturbation.
             radiusMaximum=Perturbation_Maximum_Radius(epsilonPerturbation)
             ! Find the eta-factor (see Lahav et al. 1991) which measures the dark energy contribution to the energy of the
             ! perturbation.
             eta=2.0d0*(OmegaDE/OmegaM)*radiusMaximum**3
             ! Handle the open universe case separately.
             if (OmegaDE == 0.0d0) then
                radiiRatio=0.5d0
             else
                ! Coefficients of the cubic energy equation.
                a=cmplx(  2.0d0*eta ,0.0d0,kind=kind_dble)
                b=cmplx(  0.0d0     ,0.0d0,kind=kind_dble)
                c=cmplx(-(2.0d0+eta),0.0d0,kind=kind_dble)
                d=cmplx(  1.0d0     ,0.0d0,kind=kind_dble)
                ! Solve the cubic equation to get
                Delta=(sqrt(3.0d0)*sqrt(27.0d0*a**4*d**2+4.0d0*a**3*c**3)-9.0d0*a**2*d)**(1.0d0/3.0d0)
                radiiRatio=real(cmplx(1.0d0,-sqrt(3.0d0),kind=kind_dble)*c/2.0d0**(2.0d0/3.0d0)/3.0d0**(1.0d0/3.0d0)/Delta-cmplx(1.0d0,sqrt(3.0d0),kind=kind_dble)&
                     & *Delta /2.0d0/a /2.0d0**(1.0d0/3.0d0)/3.0d0**(2.0d0/3.0d0))
             end if
             call deltaTable%populate(                                                                        &
                  &                   1.0d0/(radiiRatio*Perturbation_Maximum_Radius(epsilonPerturbation))**3, &
                  &                   iTime                                                                   &
                  &                  )
          end select
       end do
    end select
    return
  end subroutine Make_Table

  double precision function collapseRoot(epsilonPerturbation)
    implicit none
    double precision, intent(in   ) :: epsilonPerturbation

    ! Evaluate the root function.
    collapseRoot=tCollapse(epsilonPerturbation)-tNow
    return
  end function collapseRoot

  double precision function Perturbation_Maximum_Radius(epsilonPerturbation)
    !% Find the maximum radius of a perturbation with initial curvature {\normalfont \ttfamily epsilonPerturbation}.
    use Root_Finder
    implicit none
    double precision            , intent(in   ) :: epsilonPerturbation
    double precision            , parameter     :: toleranceAbsolute  =0.0d0, toleranceRelative=1.0d-9
    type            (rootFinder), save          :: finder
    !$omp threadprivate(finder)
    double precision                            :: aMaximumHigh             , aMaximumLow


    if (OmegaDE == 0.0d0) then
       ! No cosmological constant - simple analytic solution.
       Perturbation_Maximum_Radius=-OmegaM/epsilonPerturbation
    else
       ! Cosmological constant - use root finder.
       epsilonPerturbationShared=epsilonPerturbation
       aMaximumLow=-OmegaM/epsilonPerturbationShared
       aMaximumHigh=(OmegaM/2.0d0/OmegaDE)**(1.0d0/3.0d0)
       if      (aMaximumRoot(aMaximumHigh) > 0.0d0) then
          ! If the root function is not negative at aMaximumHigh it is due to rounding errors in the calculation of aMaximumHigh
          ! which implies that aMaximumHigh is very close to the actual root.
          Perturbation_Maximum_Radius=aMaximumHigh
       else if (aMaximumRoot(aMaximumLow ) < 0.0d0)  then
          ! If the root function is not positive at aMaximumLow it is due to rounding errors in the calculation of aMaximumLow
          ! which implies that aMaximumLow is very close to the actual root.
          Perturbation_Maximum_Radius=aMaximumLow
       else

          if (.not.finder%isInitialized()) then
             call finder%rootFunction(aMaximumRoot                       )
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          end if
          Perturbation_Maximum_Radius=finder%find(rootRange=[aMaximumLow,aMaximumHigh])
       end if
    end if
    return
  end function Perturbation_Maximum_Radius

  double precision function aMaximumRoot(aMaximum)
    !% Root function for maximum expansion radius.
    double precision, intent(in   ) :: aMaximum

    ! Evaluate the root function.
    aMaximumRoot=OmegaM/aMaximum+epsilonPerturbationShared+OmegaDE*aMaximum**2
  end function aMaximumRoot

  double precision function tCollapse(epsilonPerturbation)
    use Numerical_Integration
    implicit none
    real   (kind=c_double             ), intent(in   )       :: epsilonPerturbation
    type   (fgsl_function             )               , save :: integrandFunction
    type   (fgsl_integration_workspace)               , save :: integrationWorkspace
    logical                                           , save :: integrationReset     =.true.
    !$omp threadprivate(integrandFunction,integrationWorkspace,integrationReset)
    real   (kind=c_double             ), parameter           :: aMinimum             =0.0d0
    real   (kind=c_double             ), parameter           :: numericalLimitEpsilon=1.0d-4
    real   (kind=c_double             )                      :: aMaximum                    , aUpperLimit, &
         &                                                      tMaximum

    ! Find the maximum radius of the perturbation. (This is the solution of a cubic equation.)
    aUpperLimit=Perturbation_Maximum_Radius(epsilonPerturbation)

    ! Compute maximum value of a for numerical integration.
    aMaximum=(1.0d0-numericalLimitEpsilon)*aUpperLimit

    ! Share the epsilon parameter.
    epsilonPerturbationShared=epsilonPerturbation

    ! Integrate the perturbation equation from size zero to maximum size to get the time to maximum expansion.
    tMaximum=Integrate(aMinimum,aMaximum,Perturbation_Integrand,integrandFunction,integrationWorkspace &
         &,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6,hasSingularities=.true.,reset=integrationReset)&
         &/hubbleParameterInvGyr
    ! Add on analytic correction for region close to aMaximum.
    tMaximum=tMaximum-2.0d0*sqrt(OmegaM/aMaximum+epsilonPerturbation+OmegaDE*aMaximum**2)/(2.0d0*OmegaDE*aMaximum-OmegaM&
         &/aMaximum**2)/hubbleParameterInvGyr
    ! Time to collapse is twice the time to maximum expansion.
    tCollapse=2.0d0*tMaximum
    return
  end function tCollapse

  double precision function Perturbation_Integrand(a)
    implicit none
    double precision, intent(in   ) :: a
    double precision                :: sqrtArgument

    ! Compute the integrand.
    sqrtArgument=OmegaM+epsilonPerturbationShared*a+OmegaDE*a**3
    if (sqrtArgument > 0.0d0) then
       Perturbation_Integrand=sqrt(a/sqrtArgument)
    else
       Perturbation_Integrand=0.0d0
    end if
    return
  end function Perturbation_Integrand

  subroutine Spherical_Collapse_Matter_Lambda_Nonlinear_Mapping(time,deltaTable,linearGrowth_,cosmologyFunctions_)
    !% Tabulate the critical overdensity for collapse for the spherical collapse model.
    use Tables
    use Cosmology_Functions
    use Linear_Growth
    use Table_Labels
    use Root_Finder
    use Numerical_Integration
    use Numerical_Ranges
    use Array_Utilities
    use Arrays_Search
    implicit none
    double precision                                         , intent(in   ) :: time
    class           (table2DLinLinLin          )             , intent(inout) :: deltaTable
    class           (cosmologyFunctionsClass   ), target     , intent(inout) :: cosmologyFunctions_    
    class           (linearGrowthClass         ), target     , intent(inout) :: linearGrowth_
    integer                                     , parameter                  :: tableIncrement          =100
    integer                                     , parameter                  :: timesPerDecade          = 10
    integer                                     , parameter                  :: overdensityLinearCount  =500
    double precision                            , parameter                  :: numericalLimitEpsilon   =  1.0d-4
    double precision                            , parameter                  :: toleranceAbsolute       =  0.0d+0, toleranceRelative         =1.0d-9
    double precision                            , parameter                  :: expansionFactorMinimum  =  1.0d-2, expansionFactorMaximum    =1.0d+0
    double precision                            , parameter                  :: overdensityLinearMinimum= -5.0d+0, overdensityLinearMaximum  =2.0d+0
    double precision                            , allocatable, dimension(:)  :: overdensityLinear                , overdensityNonlinear             , &
         &                                                                      overdensityLinearTmp             , overdensityNonLinearTmp          , &
         &                                                                      times                            , overdensitiesLinear
    double precision                                                         :: expansionFactor                  , epsilonPerturbationMaximum       , &
         &                                                                      epsilonPerturbationCollapsed     , radiusNow                        , &
         &                                                                      epsilonPerturbation              , epsilonPerturbationMinimum       , &
         &                                                                      timeMaximum                      , radiusUpperLimit                 , &
         &                                                                      normalization                    , overdensityNonlinear_            , &
         &                                                                      timesMinimum                     , timesMaximum
    type            (fgsl_function             )                             :: integrandFunction
    type            (fgsl_integration_workspace)                             :: integrationWorkspace
    type            (rootFinder                )                             :: finderPerturbation               , finderRadius
    logical                                                                  :: integrationReset
    integer                                                                  :: i                                , timeCount                        , &
         &                                                                      iOverdensityLinear               , iOverdensity                     , &
         &                                                                      iTime


    ! Find a suitable range of times to tabulate, and generate an array of times.
    timesMinimum      =min(0.5d0*time,cosmologyFunctions_%cosmicTime(expansionFactorMinimum))
    timesMaximum      =max(2.0d0*time,cosmologyFunctions_%cosmicTime(expansionFactorMaximum))
    timeCount         =int(log10(timesMaximum/timesMinimum)*dble(timesPerDecade))+1
    times             =Make_Range(timesMinimum,timesMaximum,timeCount,rangeTypeLogarithmic)
    ! Generate a range of linear overdensities.
    overdensitiesLinear=Make_Range(overdensityLinearMinimum,overdensityLinearMaximum,overdensityLinearCount,rangeTypeLinear)
    ! Create the table.
    call deltaTable%create(overdensitiesLinear,times)
    ! Iterate over times.
    do iTime=1,timeCount
       ! Get the current expansion factor.
       tNow           =times(iTime)
       expansionFactor=cosmologyFunctions_%expansionFactor(tNow)
       ! Determine the largest (i.e. least negative) value of epsilonPerturbation for which a perturbation can collapse.
       if (cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor) > 0.0d0) then
          epsilonPerturbationMaximum=-(                                                                                &
               &                       +27.0d0                                                                         &
               &                       / 4.0d0                                                                         &
               &                       *cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor)    &
               &                       *cosmologyFunctions_%omegaMatterEpochal    (expansionFactor=expansionFactor)**2 &
               &                      )**(1.0d0/3.0d0)
       else
          epsilonPerturbationMaximum=-1.0d-6
       end if
       ! Estimate a suitably negative minimum value for epsilon.
       epsilonPerturbationMinimum=-10.0d0
       ! Compute cosmological parametrers at this epoch.
       OmegaM               =cosmologyFunctions_%omegaMatterEpochal    (expansionFactor=expansionFactor)
       OmegaDE              =cosmologyFunctions_%omegaDarkEnergyEpochal(expansionFactor=expansionFactor)
       hubbleParameterInvGyr=cosmologyFunctions_%expansionRate         (                expansionFactor)
       ! Find the value of epsilon for which the perturbation just collapses at this time.
       if (.not.finderPerturbation%isInitialized()) then
          call finderPerturbation%rootFunction(collapseRoot                       )
          call finderPerturbation%tolerance   (toleranceAbsolute,toleranceRelative)
          call finderPerturbation%rangeExpand (                                                             &
               &                               rangeExpandUpward            =0.5d0                        , &
               &                               rangeExpandDownward          =2.0d0                        , &
               &                               rangeExpandType              =rangeExpandMultiplicative    , &
               &                               rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                               rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
               &                              )
       end if
       epsilonPerturbationCollapsed=finderPerturbation%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
       ! For non-collapsed regions, epsilon will be greater then that for a collapsed perturbation. Step through values until
       ! sufficiently low non-linear overdensity is reached.
       epsilonPerturbation=epsilonPerturbationCollapsed
       i=0
       do while (.true.)
          i                  =i                  +1
          epsilonPerturbation=epsilonPerturbation+1.0d-2*abs(epsilonPerturbationCollapsed)
          ! Share the epsilon parameter.
          epsilonPerturbationShared=epsilonPerturbation
          ! For collapsing perturbations, find the time of maximum radius.
          if (epsilonPerturbation > epsilonPerturbationMaximum) then
             ! This perturbation will not collapse. Maximum radius is reached at infinite time.
             radiusMaximum=huge(1.0d0)
             timeTarget   =     tNow
          else
             ! This perturbation will collapse. Find the maximum radius.
             radiusMaximum=Perturbation_Maximum_Radius(epsilonPerturbation)
             ! Compute maximum value of a for numerical integration.
             radiusUpperLimit=(1.0d0-numericalLimitEpsilon)*radiusMaximum
             ! Integrate the perturbation equation from size zero to maximum size to get the time to maximum expansion, adding on the
             ! analytic correction for the region close to maximum expansion.
             integrationReset=.true.
             timeMaximum     =+Integrate(                                          &
                  &                                        0.0d0                 , &
                  &                                        radiusUpperLimit      , &
                  &                                        Perturbation_Integrand, &
                  &                                        integrandFunction     , &
                  &                                        integrationWorkspace  , &
                  &                      toleranceAbsolute=0.0d0                 , &
                  &                      toleranceRelative=1.0d-6                , &
                  &                      hasSingularities =.true.                , &
                  &                      reset            =integrationReset        &
                  &                     )                                          &
                  &           /hubbleParameterInvGyr                               &
                  &           -2.0d0                                               &
                  &           *sqrt(                                               &
                  &                 +OmegaM                                        &
                  &                 /radiusUpperLimit                              &
                  &                 +epsilonPerturbation                           &
                  &                 +OmegaDE                                       &
                  &                 *radiusUpperLimit   **2                        &
                  &                )                                               &
                  &           /(                                                   &
                  &             +2.0d0                                             &
                  &             *OmegaDE                                           &
                  &             *radiusUpperLimit                                  &
                  &             -OmegaM                                            &
                  &             /radiusUpperLimit**2                               &
                  &            )&
                  &           /hubbleParameterInvGyr
             call Integrate_Done(integrandFunction,integrationWorkspace)
             ! Set the target time
             if (timeMaximum > tNow) then
                ! Expanding phase.
                timeTarget=+      tNow
             else
                ! Collapsing phase.
                timeTarget=+2.0d0*timeMaximum &
                     &     -      tNow
             end if
          end if
          ! Solve for the radius at the present time.
          if (.not.finderRadius%isInitialized()) then
             call finderRadius%rootFunction(Radius_Root                       )
             call finderRadius%tolerance   (toleranceAbsolute,toleranceRelative)
             call finderRadius%rangeExpand (                                                             &
                  &                         rangeExpandDownward          =0.5d0                        , &
                  &                         rangeExpandUpward            =2.0d0                        , &
                  &                         rangeExpandType              =rangeExpandMultiplicative    , &
                  &                         rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
                  &                         rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive  &
                  &                        )
          end if
          if (epsilonPerturbation <= epsilonPerturbationMaximum .and. Radius_Root(radiusMaximum) < 0.0d0) then
             ! Perturbation is close to maximum radius. Adopt this as the solution.
             radiusNow=radiusMaximum
          else
             ! Find the current radius.
             radiusNow=finderRadius%find(rootGuess=1.0d0)
          end if
          normalization=+linearGrowth_%value(tNow           ,normalize=normalizeMatterDominated) &
               &        /                    expansionFactor
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
          overdensityLinear   (i)=+normalization         &
               &                  *0.6d0                 &
               &                  *(                     &
               &                    +1.0d0               &
               &                    -OmegaM              &
               &                    -OmegaDE             &
               &                    -epsilonPerturbation &
               &                   )                     &
               &                  /OmegaM
          overdensityNonLinear(i)=+1.0d0                 &
               &                  /radiusNow**3          &
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
             ! Tabulated overdensity exceeds any we've computed, so this overdensity is already collapsed. Use highest nonlinear overdensity.
             overdensityNonLinear_=overdensityNonLinear(i)
          else
             ! Find the tabulated in those computed and interpolate.
             iOverdensityLinear=int(Search_Array(overdensityLinear,overdensitiesLinear(iOverdensity)))
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
          call deltaTable%populate(overdensityNonLinear_,iOverdensity,iTime)
       end do
    end do
    return
  end subroutine Spherical_Collapse_Matter_Lambda_Nonlinear_Mapping

  double precision function Radius_Root(radiusNow)
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: radiusNow
    double precision                            , parameter     :: numericalLimitEpsilon=1.0d-4
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    logical                                                     :: integrationReset 
    double precision                                            :: radiusUpperLimit

    radiusUpperLimit=min(                                              &
         &               +(1.0d0-numericalLimitEpsilon)*radiusMaximum, &
         &               +                              radiusNow      &
         &              )
    integrationReset=.true.
    Radius_Root     =+Integrate(                                          &
         &                                        0.0d0                 , &
         &                                        radiusUpperLimit      , &
         &                                        Perturbation_Integrand, &
         &                                        integrandFunction     , &
         &                                        integrationWorkspace  , &
         &                      toleranceAbsolute=0.0d+0                , &
         &                      toleranceRelative=1.0d-6                , &
         &                      hasSingularities =.true.                , &
         &                      reset            =integrationReset        &
         &                      )                                         &
         &           /hubbleParameterInvGyr                               &
         &           -timeTarget
    call Integrate_Done(integrandFunction,integrationWorkspace)
    if (radiusUpperLimit < radiusNow) then
       Radius_Root       =+Radius_Root                                                                                                                                                                 &
            &             -2.0d0*sqrt(OmegaM/radiusUpperLimit+epsilonPerturbationShared+OmegaDE*radiusUpperLimit**2)/(2.0d0*OmegaDE*radiusUpperLimit-OmegaM/radiusUpperLimit**2)/hubbleParameterInvGyr
       if (radiusNow < radiusMaximum)                                                                                                                                                                  &
            & Radius_Root=+Radius_Root                                                                                                                                                                 &
            &             +2.0d0*sqrt(OmegaM/radiusNow       +epsilonPerturbationShared+OmegaDE*radiusNow       **2)/(2.0d0*OmegaDE*radiusNow       -OmegaM/radiusNow       **2)/hubbleParameterInvGyr 
    end if
    return
  end function Radius_Root
  
end module Spherical_Collapse_Matter_Lambda
