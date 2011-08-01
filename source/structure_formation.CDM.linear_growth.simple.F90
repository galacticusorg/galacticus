!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations of linear growth factor in simple cosomologies. Ignores pressure terms for the
!% growth of baryons and has no wavenumber dependence. Also assumes no growth of radiation perturbations.

module Linear_Growth_Simple
  !% Implements calculations of linear growth factor in simple cosmologies. Ignores pressure terms for the growth of baryons and
  !% has no wavenumber dependence. Also assumes no growth of radiation perturbations.
  use Cosmological_Parameters
  use Cosmology_Functions
  use, intrinsic :: ISO_C_Binding
  use FGSL
  implicit none
  private
  public :: Growth_Factor_Simple_Initialize, Linear_Growth_Simple_State_Store,&
       & Linear_Growth_Simple_State_Retrieve
  
  ! Variables to hold table of growth factor vs. cosmic time.
  double precision                    :: growthTableTimeMinimum=1.0d0, growthTableTimeMaximum=20.0d0
  double precision                    :: growthTableExpansionFactorMinimum
  integer,                  parameter :: growthTableNPointsPerDecade=1000

  ! Variables used in the ODE solver.
  type(fgsl_odeiv_step)               :: odeStepper
  type(fgsl_odeiv_control)            :: odeController
  type(fgsl_odeiv_evolve)             :: odeEvolver
  type(fgsl_odeiv_system)             :: odeSystem
  logical                             :: odeReset=.true. ! Ensure ODE variables will be reset on first call.

contains

  !# <linearGrowthMethod>
  !#  <unitName>Growth_Factor_Simple_Initialize</unitName>
  !# </linearGrowthMethod>
  subroutine Growth_Factor_Simple_Initialize(linearGrowthMethod,Linear_Growth_Tabulate)
    !% Initializes the simple growth factor module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: linearGrowthMethod
    procedure(),          pointer, intent(inout) :: Linear_Growth_Tabulate
    
    if (linearGrowthMethod == 'simple') Linear_Growth_Tabulate => Linear_Growth_Factor_Simple_Tabulate
    return
  end subroutine Growth_Factor_Simple_Initialize

  subroutine Linear_Growth_Factor_Simple_Tabulate(time,growthTableNumberPoints,growthTableTime,growthTableWavenumber&
       &,growthTableGrowthFactor ,normalizationMatterDominated)
    !% Returns the linear growth factor $D(a)$ for expansion factor {\tt aExpansion}, normalized such that $D(1)=1$ for a simple
    !% matter plus cosmological constant cosmology.
    use Numerical_Interpolation
    use Galacticus_Error
    use Memory_Management
    use Numerical_Ranges
    use ODE_Solver
    use Cosmology_Functions_Parameters
    implicit none
    double precision, intent(in)                                   :: time
    integer,          intent(out)                                  :: growthTableNumberPoints
    double precision, intent(inout), allocatable, dimension(:)     :: growthTableTime,growthTableWavenumber
    double precision, intent(inout), allocatable, dimension(:,:,:) :: growthTableGrowthFactor
    double precision, intent(out),                dimension(3)     :: normalizationMatterDominated
    double precision, parameter                                    :: dominateFactor=1.0d4
    double precision, parameter                                    :: odeToleranceAbsolute=1.0d-10, odeToleranceRelative=1.0d-10
    integer                                                        :: iTime,iComponent
    double precision                                               :: tMatterDominant,timeNow,tPresent,linearGrowthFactorPresent&
         &,growthFactorODEVariables(2),growthFactorDerivative,aMatterDominant
    type(c_ptr)                                                    :: parameterPointer
    type(fgsl_interp)                                              :: interpolationObject
    type(fgsl_interp_accel)                                        :: interpolationAccelerator
    logical                                                        :: resetInterpolation=.true.

    ! Find epoch of matter-dark energy equality.
    aMatterDominant=min(Expansion_Factor(growthTableTimeMinimum),Epoch_of_Matter_Domination(dominateFactor))
    tMatterDominant=Cosmology_Age(aMatterDominant)

    ! Find minimum and maximum times to tabulate.
    growthTableTimeMinimum=min(growthTableTimeMinimum,min(time/2.0,tMatterDominant))
    growthTableTimeMaximum=max(growthTableTimeMaximum,max(time,tMatterDominant)*2.0d0)
    
    ! Determine number of points to tabulate.
    growthTableNumberPoints=int(dlog10(growthTableTimeMaximum/growthTableTimeMinimum)*dble(growthTableNPointsPerDecade))

    ! Deallocate arrays if currently allocated.
    if (allocated(growthTableTime        )) call Dealloc_Array(growthTableTime        )
    if (allocated(growthTableWavenumber  )) call Dealloc_Array(growthTableWavenumber  )
    if (allocated(growthTableGrowthFactor)) call Dealloc_Array(growthTableGrowthFactor)
    ! Allocate the arrays to current required size.
    call Alloc_Array(growthTableTime        ,[    growthTableNumberPoints])
    call Alloc_Array(growthTableWavenumber  ,[  1                        ])
    call Alloc_Array(growthTableGrowthFactor,[3,1,growthTableNumberPoints])
    
    ! Set an arbitrary wavenumber (we do not tabulate as a function of wavenumber so the value does not matter).
    growthTableWavenumber=1.0d0

    ! Create set of grid points in time variable.
    growthTableTime=Make_Range(growthTableTimeMinimum,growthTableTimeMaximum,growthTableNumberPoints,rangeTypeLogarithmic)
    
    ! Solve ODE to get corresponding expansion factors. Initialize with solution for matter dominated phase.
    growthTableGrowthFactor(:,:,1)=1.0d0
    growthFactorDerivative=Expansion_Rate(Expansion_Factor(growthTableTime(1)))
    do iTime=2,growthTableNumberPoints
       timeNow=growthTableTime(iTime-1)
       growthFactorODEVariables(1)=growthTableGrowthFactor(1,1,iTime-1)
       growthFactorODEVariables(2)=growthFactorDerivative
       call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,timeNow,growthTableTime(iTime),2,growthFactorODEVariables &
            &,growthTableODEs,parameterPointer,odeToleranceAbsolute,odeToleranceRelative,reset=odeReset)
       growthTableGrowthFactor(1:2,1,iTime)=growthFactorODEVariables(1)
       growthFactorDerivative              =growthFactorODEVariables(2)
    end do
    growthTableGrowthFactor(3,1,:)=1.0d0 ! Assume no growth in radiation.
    
    ! Flag that interpolation must be reset.
    call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
    resetInterpolation=.true.
    
    ! Normalize to growth factor of unity at present day.
    tPresent=Cosmology_Age(1.0d0)
    do iComponent=1,3
       linearGrowthFactorPresent=Interpolate(growthTableNumberPoints,growthTableTime,growthTableGrowthFactor(iComponent,1,:)&
            &,interpolationObject,interpolationAccelerator,tPresent,reset=resetInterpolation)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
       growthTableGrowthFactor(iComponent,1,:)=growthTableGrowthFactor(iComponent,1,:)/linearGrowthFactorPresent
    end do

    ! Compute relative normalization factor such that growth factor behaves as expansion factor at early times.
    normalizationMatterDominated(:)=((9.0d0*Omega_0()/4.0d0)**(1.0d0/3.0d0))*((H_0_invGyr()*growthTableTime(1))**(2.0d0/3.0d0)) &
         &/growthTableGrowthFactor(:,1,1)
    return
  end subroutine Linear_Growth_Factor_Simple_Tabulate
  
  function growthTableODEs(t,D,dDdt,parameterPointer) bind(c)
    !% System of differential equations to solve for the growth factor.
    integer(c_int)                           :: growthTableODEs
    real(c_double), value                    :: t
    real(c_double), dimension(2), intent(in) :: D
    real(c_double), dimension(2)             :: dDdt
    type(c_ptr),    value                    :: parameterPointer
    real(c_double)                           :: aExpansion

    aExpansion=Expansion_Factor(t)
    dDdt(1)=D(2)
    dDdt(2)=1.5d0*(Expansion_Rate(aExpansion)**2)*Omega_Matter(aExpansion=aExpansion)*D(1)-2.0d0*Expansion_Rate(aExpansion)*D(2)
    growthTableODEs=FGSL_Success
  end function growthTableODEs
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Linear_Growth_Simple_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Linear_Growth_Simple_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) growthTableTimeMinimum,growthTableTimeMaximum
    return
  end subroutine Linear_Growth_Simple_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Linear_Growth_Simple_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Linear_Growth_Simple_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) growthTableTimeMinimum,growthTableTimeMaximum
    return
  end subroutine Linear_Growth_Simple_State_Retrieve
  
end module Linear_Growth_Simple
