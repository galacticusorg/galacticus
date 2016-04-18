!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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
  use FGSL
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Spherical_Collapse_Matter_Lambda_Critical_Overdensity_Tabulate,&
       & Spherical_Collape_Matter_Lambda_Delta_Virial_Tabulate,&
       & Spherical_Collapse_Matter_Lambda_State_Store, Spherical_Collapse_Matter_Lambda_State_Retrieve

  ! Variables to hold the tabulated critical overdensity data.
  double precision            :: deltaTableTimeMaximum     =20.0d0, deltaTableTimeMinimum =1.0d0
  integer         , parameter :: deltaTableNPointsPerDecade=1000
  !$omp threadprivate(deltaTableTimeMaximum,deltaTableTimeMinimum)
  ! Variables used in root finding.
  double precision            :: OmegaDE                          , OmegaM                      , &
       &                         epsilonPerturbationShared        , hubbleParameterInvGyr       , &
       &                         tNow
  !$omp threadprivate(OmegaDE,OmegaM,epsilonPerturbationShared,hubbleParameterInvGyr,tNow)
  
  ! Calculation types.
  integer         , parameter :: calculationDeltaCrit      =0     , calculationDeltaVirial=1

contains

  subroutine Spherical_Collapse_Matter_Lambda_Critical_Overdensity_Tabulate(time,deltaCritTable,linearGrowth_,cosmologyFunctions_)
    !% Tabulate the critical overdensity for collapse for the spherical collapse model.
    use Tables
    use Cosmology_Functions
    use Linear_Growth
    implicit none
    double precision                                      , intent(in   ) :: time
    class           (table1D                ), allocatable, intent(inout) :: deltaCritTable
    class           (cosmologyFunctionsClass), target     , intent(in   ) :: cosmologyFunctions_    
    class           (linearGrowthClass      ), target     , intent(in   ) :: linearGrowth_    

    call Make_Table(time,deltaCritTable,calculationDeltaCrit,linearGrowth_,cosmologyFunctions_)
    return
  end subroutine Spherical_Collapse_Matter_Lambda_Critical_Overdensity_Tabulate

  subroutine Spherical_Collape_Matter_Lambda_Delta_Virial_Tabulate(time,deltaVirialTable)
    !% Tabulate the virial density contrast for the spherical collapse model.
    use Tables
    implicit none
    double precision                      , intent(in   ) :: time
    class           (table1D), allocatable, intent(inout) :: deltaVirialTable

    call Make_Table(time,deltaVirialTable,calculationDeltaVirial)
    return
  end subroutine Spherical_Collape_Matter_Lambda_Delta_Virial_Tabulate

  subroutine Make_Table(time,deltaTable,calculationType,linearGrowth_,cosmologyFunctions_)
    !% Tabulate $\delta_{\mathrm crit}$ or $\Delta_{\mathrm vir}$ vs. time.
    use Linear_Growth
    use Cosmology_Functions
    use Root_Finder
    use Tables
    use Kind_Numbers
    use Galacticus_Error
    implicit none
    double precision                                      , intent(in   )           :: time
    integer                                               , intent(in   )           :: calculationType
    class           (table1D                ), allocatable, intent(inout)           :: deltaTable
    class           (cosmologyFunctionsClass), target     , intent(in   ), optional :: cosmologyFunctions_    
    class           (linearGrowthClass      ), target     , intent(in   ), optional :: linearGrowth_    
    double precision                         , parameter                            :: toleranceAbsolute         =0.0d0, toleranceRelative         =1.0d-9
    type            (rootFinder             ), save                                 :: finder
    !$omp threadprivate(finder)
    class           (cosmologyFunctionsClass), pointer                              :: cosmologyFunctions__
    class           (linearGrowthClass      ), pointer                              :: linearGrowth__
    integer                                                                         :: deltaTableNumberPoints          , iTime
    double precision                                                                :: aExpansionNow                   , epsilonPerturbation              , &
         &                                                                             epsilonPerturbationMaximum      , epsilonPerturbationMinimum       , &
         &                                                                             eta                             , normalization                    , &
         &                                                                             radiiRatio                      , radiusMaximum
    double complex                                                                  :: a,b,c,d,Delta

    ! Get required objects.
    if (present(cosmologyFunctions_)) then
       cosmologyFunctions__ => cosmologyFunctions_
    else
       cosmologyFunctions__ => cosmologyFunctions()
    end if
    if (present(linearGrowth_      )) then
       linearGrowth__       => linearGrowth_      
    else
       linearGrowth__       => linearGrowth      ()
    end if
    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)

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

          ! Get the current expansion factor.
          aExpansionNow=cosmologyFunctions__%expansionFactor(deltaTable%x(iTime))
          ! Determine the largest (i.e. least negative) value of epsilonPerturbation for which a perturbation can collapse.
          if (cosmologyFunctions__%omegaDarkEnergyEpochal(expansionFactor=aExpansionNow)>0.0d0) then
             epsilonPerturbationMaximum=-(27.0d0*cosmologyFunctions__%omegaDarkEnergyEpochal(expansionFactor=aExpansionNow)*(cosmologyFunctions__%omegaMatterEpochal(expansionFactor=aExpansionNow)&
                  & **2)/4.0d0)**(1.0d0/3.0d0)
          else
             epsilonPerturbationMaximum=-1.0d-6
          end if

          ! Estimate a suitably negative minimum value for epsilon.
          epsilonPerturbationMinimum=-10.0d0

          OmegaM               =cosmologyFunctions__%omegaMatterEpochal    (expansionFactor=aExpansionNow)
          OmegaDE              =cosmologyFunctions__%omegaDarkEnergyEpochal(expansionFactor=aExpansionNow)
          hubbleParameterInvGyr=cosmologyFunctions__%expansionRate         (                aExpansionNow)
          tNow                 =deltaTable%x(iTime)

          ! Find the value of epsilon for which the perturbation just collapses at this time.
          if (.not.finder%isInitialized()) then
             call finder%rootFunction(collapseRoot                       )
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
             call finder%rangeExpand (                                                           &
                  &                   rangeExpandUpward          =0.5d0                        , &
                  &                   rangeExpandType            =rangeExpandMultiplicative    , &
                  &                   rangeExpandUpwardSignExpect=rangeExpandSignExpectPositive  &
                  &                  )
          end if
          epsilonPerturbation=finder%find(rootRange=[epsilonPerturbationMinimum,epsilonPerturbationMaximum])
          ! Compute the corresponding critical overdensity.
          normalization=linearGrowth__%value(tNow,normalize=normalizeMatterDominated)/linearGrowth__%value(tNow)/aExpansionNow
          select case (calculationType)
          case (calculationDeltaCrit)
             ! Critical linear overdensity.
             call deltaTable%populate(                                                                       &
                  &                   normalization*0.6d0*(1.0d0-OmegaM-OmegaDE-epsilonPerturbation)/OmegaM, &
                  &                   iTime                                                                  &
                  &                  )
             ! Check for non-monotonic decline.
             if (iTime > 1) then
                if (deltaTable%y(iTime) >= deltaTable%y(iTime-1)) call Galacticus_Error_Report('Make_Table','accuracy lost in tablulation of critical overdensity (usually results for computing critical overdensity for very large cosmic times)')
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
    type   (c_ptr                     )                      :: parameterPointer
    real   (kind=c_double             )                      :: aMaximum                    , aUpperLimit, &
         &                                                      tMaximum

    ! Find the maximum radius of the perturbation. (This is the solution of a cubic equation.)
    aUpperLimit=Perturbation_Maximum_Radius(epsilonPerturbation)

    ! Compute maximum value of a for numerical integration.
    aMaximum=(1.0d0-numericalLimitEpsilon)*aUpperLimit

    ! Share the epsilon parameter.
    epsilonPerturbationShared=epsilonPerturbation

    ! Integrate the perturbation equation from size zero to maximum size to get the time to maximum expansion.
    tMaximum=Integrate(aMinimum,aMaximum,Perturbation_Integrand,parameterPointer,integrandFunction,integrationWorkspace &
         &,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-6,hasSingularities=.true.,reset=integrationReset)&
         &/hubbleParameterInvGyr
    ! Add on analytic correction for region close to aMaximum.
    tMaximum=tMaximum-2.0d0*sqrt(OmegaM/aMaximum+epsilonPerturbation+OmegaDE*aMaximum**2)/(2.0d0*OmegaDE*aMaximum-OmegaM&
         &/aMaximum**2)/hubbleParameterInvGyr
    ! Time to collapse is twice the time to maximum expansion.
    tCollapse=2.0d0*tMaximum
    return
  end function tCollapse

  function Perturbation_Integrand(a,parameterPointer) bind(c)
    implicit none
    real(kind=c_double)        :: Perturbation_Integrand
    real(kind=c_double), value :: a
    type(c_ptr        ), value :: parameterPointer
    real(kind=c_double)        :: sqrtArgument

    ! Compute the integrand.
    sqrtArgument=OmegaM+epsilonPerturbationShared*a+OmegaDE*a**3
    if (sqrtArgument>0.0d0) then
       Perturbation_Integrand=sqrt(a/sqrtArgument)
    else
       Perturbation_Integrand=0.0d0
    end if
    return
  end function Perturbation_Integrand

  !# <galacticusStateStoreTask>
  !#  <unitName>Spherical_Collapse_Matter_Lambda_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Spherical_Collapse_Matter_Lambda_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile
    !GCC$ attributes unused :: fgslStateFile
    
    write (stateFile) deltaTableTimeMinimum,deltaTableTimeMaximum
    return
  end subroutine Spherical_Collapse_Matter_Lambda_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Spherical_Collapse_Matter_Lambda_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Spherical_Collapse_Matter_Lambda_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile
    !GCC$ attributes unused :: fgslStateFile

    read (stateFile) deltaTableTimeMinimum,deltaTableTimeMaximum
    return
  end subroutine Spherical_Collapse_Matter_Lambda_State_Retrieve

end module Spherical_Collapse_Matter_Lambda
