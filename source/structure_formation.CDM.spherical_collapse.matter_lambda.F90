!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements calculations of spherical top hat collapse in cosmologies containing matter and a
!% cosmological constant.

module Spherical_Collapse_Matter_Lambda
  use FGSL
  use, intrinsic :: ISO_C_Binding
  private
  public :: Spherical_Collape_Delta_Critical_Initialize, Spherical_Collapse_Critical_Overdensity,&
       & Spherical_Collape_Delta_Virial_Initialize, Spherical_Collapse_Virial_Density_Contrast,&
       & Spherical_Collapse_Matter_Lambda_State_Store, Spherical_Collapse_Matter_Lambda_State_Retrieve

  ! Variables to hold the tabulated critical overdensity data.
  double precision               :: deltaTableTimeMinimum=1.0d0, deltaTableTimeMaximum=20.0d0
  integer,             parameter :: deltaTableNPointsPerDecade=300

  ! Variables used in root finding.
  double precision               :: OmegaDE,OmegaM,tNow,hubbleParameterInvGyr,epsilonPerturbationShared

  ! Calculation types.
  integer,             parameter :: calculationDeltaCrit=0, calculationDeltaVirial=1

contains

  !# <criticalOverdensityMethod>
  !#  <unitName>Spherical_Collape_Delta_Critical_Initialize</unitName>
  !# </criticalOverdensityMethod>
  subroutine Spherical_Collape_Delta_Critical_Initialize(criticalOverdensityMethod,Critical_Overdensity_Tabulate)
    !% Initializes the $\delta_{\rm crit}$ calculation for the spherical collapse module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: criticalOverdensityMethod
    procedure(),          pointer, intent(inout) :: Critical_Overdensity_Tabulate
    
    if (criticalOverdensityMethod.eq.'spherical top hat') Critical_Overdensity_Tabulate => Spherical_Collapse_Critical_Overdensity
    return
  end subroutine Spherical_Collape_Delta_Critical_Initialize

  !# <virialDensityContrastMethod>
  !#  <unitName>Spherical_Collape_Delta_Virial_Initialize</unitName>
  !# </virialDensityContrastMethod>
  subroutine Spherical_Collape_Delta_Virial_Initialize(virialDensityContrastMethod,Virial_Density_Contrast_Tabulate)
    !% Initializes the $\Delta_{\rm vir}$ calculation for the spherical collapse module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: virialDensityContrastMethod
    procedure(),          pointer, intent(inout) :: Virial_Density_Contrast_Tabulate
    
    if (virialDensityContrastMethod.eq.'spherical top hat') Virial_Density_Contrast_Tabulate =>&
         & Spherical_Collapse_Virial_Density_Contrast
    return
  end subroutine Spherical_Collape_Delta_Virial_Initialize

  subroutine Spherical_Collapse_Critical_Overdensity(time,deltaCritNumberPoints,deltaCritTime,deltaCritDeltaCrit)
    !% Tabulate the critical overdensity for collapse for the spherical collapse model.
    implicit none
    double precision, intent(in)                               :: time
    integer,          intent(out)                              :: deltaCritNumberPoints
    double precision, intent(inout), allocatable, dimension(:) :: deltaCritTime,deltaCritDeltaCrit

    !$omp critical(Spherical_Collapse_Make_Table)
    call Make_Table(time,deltaCritNumberPoints,deltaCritTime,deltaCritDeltaCrit,calculationDeltaCrit)
    !$omp end critical(Spherical_Collapse_Make_Table)

    return
  end subroutine Spherical_Collapse_Critical_Overdensity

  subroutine Spherical_Collapse_Virial_Density_Contrast(time,deltaVirialNumberPoints,deltaVirialTime,deltaVirialDeltaVirial)
    !% Tabulate the virial density contrast for the spherical collapse model.
    implicit none
    double precision, intent(in)                               :: time
    integer,          intent(out)                              :: deltaVirialNumberPoints
    double precision, intent(inout), allocatable, dimension(:) :: deltaVirialTime,deltaVirialDeltaVirial

    !$omp critical(Spherical_Collapse_Make_Table)
    call Make_Table(time,deltaVirialNumberPoints,deltaVirialTime,deltaVirialDeltaVirial,calculationDeltaVirial)
    !$omp end critical(Spherical_Collapse_Make_Table)

    return
  end subroutine Spherical_Collapse_Virial_Density_Contrast

  subroutine Make_Table(time,deltaTableNumberPoints,deltaTableTime,deltaTableDelta,calculationType)
    !% Tabulate $\delta_{\rm crit}$ or $\Delta_{\rm vir}$ vs. time.
    use Linear_Growth
    use Cosmology_Functions
    use Root_Finder
    use Numerical_Ranges
    use Numerical_Interpolation
    use Memory_Management
    implicit none
    double precision,        intent(in)                               :: time
    integer,                 intent(out)                              :: deltaTableNumberPoints
    double precision,        intent(inout), allocatable, dimension(:) :: deltaTableTime,deltaTableDelta
    integer,                 intent(in)                               :: calculationType
    type(fgsl_function),     save                                     :: rootFunction
    type(fgsl_root_fsolver), save                                     :: rootFunctionSolver
    double precision,        parameter                                :: toleranceAbsolute=1.0d-9,toleranceRelative=1.0d-9
    integer                                                           :: iTime
    type(c_ptr)                                                       :: parameterPointer
    double precision                                                  :: epsilonPerturbation,epsilonPerturbationMinimum &
         &,epsilonPerturbationMaximum,aExpansionNow,normalization

    ! Find minimum and maximum times to tabulate.
    deltaTableTimeMinimum=min(deltaTableTimeMinimum,time/2.0d0)
    deltaTableTimeMaximum=max(deltaTableTimeMaximum,time*2.0d0)
    
    ! Determine number of points to tabulate.
    deltaTableNumberPoints=int(dlog10(deltaTableTimeMaximum/deltaTableTimeMinimum)&
         &*dble(deltaTableNPointsPerDecade))
    
    ! Deallocate arrays if currently allocated.
    if (allocated(deltaTableTime )) call Dealloc_Array(deltaTableTime )
    if (allocated(deltaTableDelta)) call Dealloc_Array(deltaTableDelta)
    ! Allocate the arrays to current required size.
    call Alloc_Array(deltaTableTime ,[deltaTableNumberPoints])
    call Alloc_Array(deltaTableDelta,[deltaTableNumberPoints])
    
    ! Create set of grid points in time variable.
    deltaTableTime=Make_Range(deltaTableTimeMinimum,deltaTableTimeMaximum,deltaTableNumberPoints,rangeTypeLogarithmic)
    
    ! Solve ODE to get corresponding overdensities.
    do iTime=1,deltaTableNumberPoints
       ! Get the current expansion factor.
       aExpansionNow=Expansion_Factor(deltaTableTime(iTime))
       ! Determine the largest (i.e. least negative) value of epsilonPerturbation for which a perturbation can collapse.
       if (Omega_Dark_Energy(aExpansion=aExpansionNow)>0.0d0) then
          epsilonPerturbationMaximum=-(27.0d0*Omega_Dark_Energy(aExpansion=aExpansionNow)*(Omega_Matter(aExpansion=aExpansionNow)&
               & **2)/4.0d0)**(1.0d0/3.0d0)
       else
          epsilonPerturbationMaximum=-1.0d-6
       end if

       ! Estimate a suitably negative minimum value for epsilon.
       epsilonPerturbationMinimum=-10.0d0

       OmegaM               =Omega_Matter     (aExpansion=aExpansionNow)
       OmegaDE              =Omega_Dark_Energy(aExpansion=aExpansionNow)
       hubbleParameterInvGyr=Expansion_Rate   (           aExpansionNow)
       tNow                 =deltaTableTime(iTime)
       
       ! Find the value of epsilon for which the perturbation just collapses at this time.
       epsilonPerturbation=Root_Find(epsilonPerturbationMinimum,epsilonPerturbationMaximum,collapseRoot,parameterPointer &
            &,rootFunction,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
       
       ! Compute the corresponding critical overdensity.
       normalization=Linear_Growth_Factor(tNow,normalize=normalizeMatterDominated)/Linear_Growth_Factor(tNow)/aExpansionNow
       select case (calculationType)
       case (calculationDeltaCrit)
          ! Critical linear overdensity.
          deltaTableDelta(iTime)=normalization*0.6d0*(1.0d0-OmegaM-OmegaDE-epsilonPerturbation)/OmegaM
       case (calculationDeltaVirial)
          ! Virial density contrast at collapse.
          deltaTableDelta(iTime)=8.0d0/(Perturbation_Maximum_Radius(epsilonPerturbation)**3)
       end select
    end do
    return
  end subroutine Make_Table

  function collapseRoot(epsilonPerturbation,parameterPointer) bind(c)
    real(c_double), value   :: epsilonPerturbation
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: collapseRoot
   
    ! Evaluate the root function.
    collapseRoot=tCollapse(epsilonPerturbation)-tNow
  end function collapseRoot

  double precision function Perturbation_Maximum_Radius(epsilonPerturbation)
    !% Find the maximum radius of a perturbation with initial curvature {\tt epsilonPerturbation}.
    use Root_Finder
    implicit none
    double precision,        intent(in) :: epsilonPerturbation
    type(fgsl_function),     save       :: rootFunction
    type(fgsl_root_fsolver), save       :: rootFunctionSolver
    double precision,        parameter  :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-9
    type(c_ptr)                         :: parameterPointer
    double precision                    :: aMaximumLow,aMaximumHigh

    if (OmegaDE.eq.0.0d0) then
       ! No cosmological constant - simple analytic solution.
       Perturbation_Maximum_Radius=-OmegaM/epsilonPerturbation
    else
       ! Cosmological constant - use root finder.
       epsilonPerturbationShared=epsilonPerturbation
       aMaximumLow=-OmegaM/epsilonPerturbationShared
       aMaximumHigh=(OmegaM/2.0d0/OmegaDE)**(1.0d0/3.0d0)
       if (aMaximumRoot(aMaximumHigh,parameterPointer)>0.0d0) then
          ! If the root function is not negative at aMaximumHigh it is due to rounding errors in the calculation of aMaximumHigh
          ! which implies that aMaximumHigh is very close to the actual root.
          Perturbation_Maximum_Radius=aMaximumHigh          
       else
          Perturbation_Maximum_Radius=Root_Find(aMaximumLow,aMaximumHigh,aMaximumRoot,parameterPointer &
               &,rootFunction,rootFunctionSolver,toleranceAbsolute,toleranceRelative)
       end if
    end if
    return
  end function Perturbation_Maximum_Radius

  function aMaximumRoot(aMaximum,parameterPointer) bind(c)
    !% Root function for maximum expansion radius.
    real(c_double), value   :: aMaximum
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: aMaximumRoot
   
    ! Evaluate the root function.
    aMaximumRoot=OmegaM/aMaximum+epsilonPerturbationShared+OmegaDE*aMaximum**2
  end function aMaximumRoot

  double precision function tCollapse(epsilonPerturbation)
    use Numerical_Integration
    implicit none
    real(c_double),                   intent(in) :: epsilonPerturbation
    type(fgsl_function),              save       :: integrandFunction
    type(fgsl_integration_workspace), save       :: integrationWorkspace
    logical,                          save       :: integrationReset=.true.
    real(c_double),                   parameter  :: aMinimum=0.0d0
    real(c_double),                   parameter  :: numericalLimitEpsilon=1.0d-4
    type(c_ptr)                                  :: parameterPointer
    real(c_double)                               :: aMaximum,aUpperLimit,tMaximum

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
    tMaximum=tMaximum-2.0d0*dsqrt(OmegaM/aMaximum+epsilonPerturbation+OmegaDE*aMaximum**2)/(2.0d0*OmegaDE*aMaximum-OmegaM&
         &/aMaximum**2)/hubbleParameterInvGyr
    ! Time to collapse is twice the time to maximum expansion.
    tCollapse=2.0d0*tMaximum
    return
  end function tCollapse

  function Perturbation_Integrand(a,parameterPointer) bind(c)
    implicit none
    real(c_double)          :: Perturbation_Integrand
    real(c_double), value   :: a
    type(c_ptr),    value   :: parameterPointer
    real(c_double)          :: sqrtArgument

    ! Compute the integrand.
    sqrtArgument=OmegaM+epsilonPerturbationShared*a+OmegaDE*a**3
    if (sqrtArgument>0.0d0) then
       Perturbation_Integrand=dsqrt(a/sqrtArgument)
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
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) deltaTableTimeMinimum,deltaTableTimeMaximum
    return
  end subroutine Spherical_Collapse_Matter_Lambda_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Spherical_Collapse_Matter_Lambda_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Spherical_Collapse_Matter_Lambda_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) deltaTableTimeMinimum,deltaTableTimeMaximum
    return
  end subroutine Spherical_Collapse_Matter_Lambda_State_Retrieve
  
end module Spherical_Collapse_Matter_Lambda
