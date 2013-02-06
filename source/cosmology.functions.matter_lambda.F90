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

!% Contains a module which implements useful cosmological functions. This implementation assumes a Universe filled with
!% collisionless matter and a cosmological constant.

module Cosmology_Functions_Matter_Lambda
  !% Implements useful cosmological functions. This implementation assumes a Universe filled with collisionless matter
  !% and a cosmological constant.
  use Cosmological_Parameters
  use FGSL
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Cosmology_Functions_Matter_Lambda_Initialize, Cosmology_Matter_Lambda_State_Store,&
       & Cosmology_Matter_Lambda_State_Retrieve

  ! Flag indicating if the module has been initialized.
  logical                                     :: moduleInitialized=.false.
  
  ! Variables used to track critical times in collapsing Universes.
  logical                                     :: collapsingUniverse=.false.
  integer                                     :: iTableTurnaround
  double precision                            :: tCosmologicalMax,tCosmologicalTurnaround,aExpansionMax
  !$omp threadprivate(collapsingUniverse,iTableTurnaround,aExpansionMax,tCosmologicalMax,tCosmologicalTurnaround)
  
  ! Factor by which one component of Universe must dominate others such that we can ignore the others.
  double precision, parameter                 :: dominateFactor=100.0d0

  ! Variables to hold table of expansion factor vs. cosmic time.
  logical                                     :: ageTableInitialized=.false.
  integer                                     :: ageTableNumberPoints
  double precision                            :: ageTableTimeMinimum=1.0d-4, ageTableTimeMaximum=20.0d0
  integer,          parameter                 :: ageTableNPointsPerDecade=300
  double precision, parameter                 :: ageTableNPointsPerOctave=dble(ageTableNPointsPerDecade)*dlog(2.0d0)/dlog(10.0d0)
  double precision, parameter                 :: ageTableIncrementFactor =dexp(int(ageTableNPointsPerOctave+1.0d0)*dlog(10.0d0)&
       &/dble(ageTableNPointsPerDecade))
  double precision, allocatable, dimension(:) :: ageTableTime, ageTableExpansionFactor
  type(fgsl_interp)                           :: interpolationObject     ,interpolationObjectInverse
  type(fgsl_interp_accel)                     :: interpolationAccelerator,interpolationAcceleratorInverse
  logical                                     :: resetInterpolation       =.true.
  logical                                     :: resetInterpolationInverse=.true.
  !$omp threadprivate(ageTableInitialized,ageTableNumberPoints,ageTableTimeMinimum,ageTableTimeMaximum,ageTableTime)
  !$omp threadprivate(ageTableExpansionFactor,interpolationObject,interpolationObjectInverse,interpolationAccelerator)
  !$omp threadprivate(interpolationAcceleratorInverse,resetInterpolation,resetInterpolationInverse)

  ! Variables to hold table of distance vs. cosmic time.
  logical                                     :: distanceTableInitialized=.false.
  integer                                     :: distanceTableNumberPoints
  double precision                            :: distanceTableTimeMinimum=1.0d-4, distanceTableTimeMaximum
  integer,          parameter                 :: distanceTableNPointsPerDecade=100
  double precision, allocatable, dimension(:) :: distanceTableTime,distanceTableComovingDistance&
       &,distanceTableLuminosityDistanceNegated ,distanceTableComovingDistanceNegated
  type(fgsl_interp)                           :: interpolationObjectDistance      ,interpolationObjectDistanceInverse      , &
       & interpolationObjectLuminosityDistance
  type(fgsl_interp_accel)                     :: interpolationAcceleratorDistance ,interpolationAcceleratorDistanceInverse , &
       & interpolationAcceleratorLuminosityDistance
  logical                                     :: resetInterpolationDistance=.true.,resetInterpolationDistanceInverse=.true., &
       & resetInterpolationLuminosityDistance=.true.

  ! Variables used in the ODE solver.
  type(fgsl_odeiv_step)                       :: odeStepper
  type(fgsl_odeiv_control)                    :: odeController
  type(fgsl_odeiv_evolve)                     :: odeEvolver
  type(fgsl_odeiv_system)                     :: odeSystem
  logical                                     :: odeReset=.true. ! Ensure ODE variables will be reset on first call.
  !$omp threadprivate (odeStepper,odeController,odeEvolver,odeSystem,odeReset)

contains

  !# <cosmologyMethod>
  !#  <unitName>Cosmology_Functions_Matter_Lambda_Initialize</unitName>
  !# </cosmologyMethod>
  subroutine Cosmology_Functions_Matter_Lambda_Initialize(cosmologyMethod,Expansion_Factor_Is_Valid_Get,Cosmic_Time_Is_Valid_Get &
       &,Cosmology_Age_Get,Expansion_Factor_Get,Hubble_Parameter_Get,Early_Time_Density_Scaling_Get,Omega_Matter_Total_Get &
       &,Omega_Dark_Energy_Get,Expansion_Rate_Get,Epoch_of_Matter_Dark_Energy_Equality_Get,Epoch_of_Matter_Domination_Get &
       &,Epoch_of_Matter_Curvature_Equality_Get,CMB_Temperature_Get,Comoving_Distance_Get,Time_From_Comoving_Distance_Get&
       &,Comoving_Distance_Conversion_Get)
    !% Initialize the module.
    use Numerical_Comparison
    use ISO_Varying_String
    use ODE_Solver
    implicit none
    type(varying_string),                 intent(in)    :: cosmologyMethod
    procedure(),                 pointer, intent(inout) :: Early_Time_Density_Scaling_Get
    procedure(logical),          pointer, intent(inout) :: Expansion_Factor_Is_Valid_Get,Cosmic_Time_Is_Valid_Get
    procedure(double precision), pointer, intent(inout) :: Cosmology_Age_Get ,Expansion_Factor_Get,Hubble_Parameter_Get &
         &,Omega_Matter_Total_Get,Omega_Dark_Energy_Get ,Expansion_Rate_Get,Epoch_of_Matter_Dark_Energy_Equality_Get &
         &,Epoch_of_Matter_Domination_Get ,Epoch_of_Matter_Curvature_Equality_Get,CMB_Temperature_Get,Comoving_Distance_Get&
         &,Time_From_Comoving_Distance_Get,Comoving_Distance_Conversion_Get
    double precision,            parameter              :: odeToleranceAbsolute=1.0d-9, odeToleranceRelative=1.0d-9
    double precision,            parameter              :: omegaTolerance=1.0d-9
    double precision                                    :: cubicTerm1,cubicTerm5,cubicTerm9,cubicTerm21Squared,cubicTerm21 &
         &,cubicTerm25Cubed,cubicTerm25,aMaximum,aDominant,timeMaximum(1),densityPower,Omega_Dominant
    type(c_ptr)                                         :: parameterPointer

    ! Check if our method is selected.
    if (cosmologyMethod == 'matter-lambda') then
       ! Set up procedure pointers.
       Expansion_Factor_Is_Valid_Get            => Expansion_Factor_Is_Valid_Matter_Lambda
       Cosmic_Time_Is_Valid_Get                 => Cosmic_Time_Is_Valid_Matter_Lambda
       Cosmology_Age_Get                        => Cosmology_Age_Matter_Lambda
       Expansion_Factor_Get                     => Expansion_Factor_Matter_Lambda
       Hubble_Parameter_Get                     => Hubble_Parameter_Matter_Lambda
       Early_Time_Density_Scaling_Get           => Early_Time_Density_Scaling_Matter_Lambda
       Omega_Matter_Total_Get                   => Omega_Matter_Total_Matter_Lambda
       Omega_Dark_Energy_Get                    => Omega_Dark_Energy_Matter_Lambda
       Expansion_Rate_Get                       => Expansion_Rate_Matter_Lambda
       Epoch_of_Matter_Dark_Energy_Equality_Get => Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda
       Epoch_of_Matter_Curvature_Equality_Get   => Epoch_of_Matter_Curvature_Equality_Matter_Lambda
       Epoch_of_Matter_Domination_Get           => Epoch_of_Matter_Domination_Matter_Lambda
       CMB_Temperature_Get                      => CMB_Temperature_Matter_Lambda
       Comoving_Distance_Get                    => Comoving_Distance_Matter_Lambda
       Time_From_Comoving_Distance_Get          => Time_From_Comoving_Distance_Matter_Lambda
       Comoving_Distance_Conversion_Get         => Comoving_Distance_Conversion_Matter_Lambda

       ! Determine if this universe will collapse. We take the Friedmann equation, which gives H^2 as a function of expansion factor,
       ! a, and solve for where H^2=0. If this has a real solution, then we have a collapsing universe.
       collapsingUniverse=.false.
       if (Values_Agree(Omega_K(),0.0d0,absTol=omegaTolerance)) then
          if (Values_Agree(Omega_DE(),0.0d0,absTol=omegaTolerance)) then
             ! Einstein-de Sitter case. Always expands to infinity.
             collapsingUniverse=.false.
          else
             ! Flat Universe with cosmological constant.
             if (Omega_DE() > 0.0d0) then
                ! Never collapses.
                collapsingUniverse=.false.
             else
                collapsingUniverse=.true.
                aExpansionMax=-(Omega_Matter()*Omega_DE()**2)**(1.0d0/3.0d0)/Omega_DE()
             end if
          end if
       else
          if (Values_Agree(Omega_DE(),0.0d0,absTol=omegaTolerance)) then
             ! Simple case for a matter-only universe.
             collapsingUniverse=Omega_Matter() > 1.0d0
             if (collapsingUniverse) aExpansionMax=Omega_Matter()/(Omega_Matter()-1.0d0)
          else
             ! Case of matter plus dark energy.
             cubicTerm1 =1.0d0/Omega_DE()
             cubicTerm5 =Omega_Matter()**2
             cubicTerm9 =Omega_DE()**2
             cubicTerm21Squared=-(-0.12d2+0.36d2*Omega_Matter()+36.0d0*Omega_DE()-0.36d2*cubicTerm5-0.72d2*Omega_Matter()*Omega_DE()-36.0d0&
                  &*cubicTerm9+0.12d2*cubicTerm5*Omega_Matter()-0.45d2*cubicTerm5*Omega_DE()+0.36d2*Omega_Matter()*cubicTerm9+12.0d0*cubicTerm9&
                  &*Omega_DE())*cubicTerm1
             if (cubicTerm21Squared > 0.0d0) then
                cubicTerm21=dsqrt(cubicTerm21Squared)
                cubicTerm25Cubed=(-0.108d3*Omega_Matter()+0.12d2*cubicTerm21)*cubicTerm9
                if (cubicTerm25Cubed >= 0.0d0) then
                   cubicTerm25=cubicTerm25Cubed**(1.0d0/3.0d0)
                else
                   cubicTerm25=-dabs(cubicTerm25Cubed)**(1.0d0/3.0d0)
                end if
                aMaximum=cubicTerm1*cubicTerm25/0.6d1+0.2d1*(-0.1d1+Omega_Matter()+Omega_DE())/cubicTerm25
                collapsingUniverse=aMaximum > 0.0d0
                if (collapsingUniverse) aExpansionMax=aMaximum
             end if
          end if
       end if
       
       ! If we have a collapsing Universe, find time of turnaround, and maximum time.
       if (collapsingUniverse) then
          ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
          call Early_Time_Density_Scaling_Matter_Lambda(dominateFactor,densityPower,aDominant,Omega_Dominant)       
          ! Find the corresponding time.
          timeMaximum(1)=1.0d0/H_0_invGyr()/dsqrt(Omega_Dominant)/aDominant**(0.5d0*densityPower)
          ! Solve Friedmann equation to get time at turnaround.
          odeReset=.true.
          call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,aDominant,aExpansionMax*(1.0d0-1.0d-4),1,timeMaximum&
               &,collapseODEs,parameterPointer,odeToleranceAbsolute,odeToleranceRelative,reset=odeReset)
          call ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
          odeReset=.true.
          ! Extract turnaround time from ODE variables and set maximum time to twice turnaround time.
          tCosmologicalTurnaround=timeMaximum(1)
          tCosmologicalMax       =2.0d0*tCosmologicalTurnaround
       end if       
    end if
    return
  end subroutine Cosmology_Functions_Matter_Lambda_Initialize
  
  function collapseODEs(a,t,dtda,parameterPointer) bind(c)
    !% System of differential equations to solve for age vs. expansion factor.
    integer(c_int)                           :: collapseODEs
    real(c_double), value                    :: a
    real(c_double), dimension(1), intent(in) :: t
    real(c_double), dimension(1)             :: dtda
    type(c_ptr),    value                    :: parameterPointer

    dtda(1)=1.0d0/a/Expansion_Rate_Matter_Lambda(a)
    collapseODEs=FGSL_Success
  end function collapseODEs
  
  logical function Expansion_Factor_Is_Valid_Matter_Lambda(aExpansion)
    !% Checks that the expansion factor falls within allowed ranges.
    implicit none
    double precision, intent(in) :: aExpansion

    Expansion_Factor_Is_Valid_Matter_Lambda=aExpansion>0.0d0 .and. (aExpansion<aExpansionMax .or. .not.collapsingUniverse)
    return
  end function Expansion_Factor_Is_Valid_Matter_Lambda
  
  logical function Cosmic_Time_Is_Valid_Matter_Lambda(time)
    !% Checks that the time falls within allowed ranges.
    implicit none
    double precision, intent(in) :: time
    
    Cosmic_Time_Is_Valid_Matter_Lambda=time>0.0d0 .and. (time<tCosmologicalMax .or. .not.collapsingUniverse)
    return
  end function Cosmic_Time_Is_Valid_Matter_Lambda
  
  double precision function Cosmology_Age_Matter_Lambda(aExpansion,collapsingPhase)
    use Galacticus_Error
    use Numerical_Interpolation
    implicit none
    double precision, intent(in)           :: aExpansion
    logical,          intent(in), optional :: collapsingPhase
    logical                                :: collapsingPhaseActual

    ! Initialize the module if necessary.
    ! Validate the input.
    if (.not.Expansion_Factor_Is_Valid_Matter_Lambda(aExpansion)) call Galacticus_Error_Report('Cosmology_Age_Matter_Lambda'&
         &,'expansion factor is invalid')

    ! Determine if we are in the expanding or collapsing phase for this universe.
    if (present(collapsingPhase)) then
       collapsingPhaseActual=collapsingPhase
    else
       collapsingPhaseActual=.false. ! Assume expanding phase by default.
    end if

    ! Ensure tabulation is initialized.
    if (.not.ageTableInitialized) call Make_Expansion_Factor_Table(ageTableTimeMinimum)

    ! Ensure that the tabulation spans a sufficient range of expansion factors.
    if (collapsingPhaseActual) then
       ! In collapsing phase just ensure that a sufficiently large expansion factor has been reached.
       do while (ageTableExpansionFactor(ageTableNumberPoints)<aExpansion)
          ageTableTimeMaximum=min(ageTableTimeMaximum*ageTableIncrementFactor,tCosmologicalTurnaround)
          call Make_Expansion_Factor_Table
       end do
    else
       ! In expanding phase ensure that sufficiently small and large expansion factors have been reached.
       do while (ageTableExpansionFactor(1)>aExpansion)
          ageTableTimeMinimum=ageTableTimeMinimum/ageTableIncrementFactor
          call Make_Expansion_Factor_Table
       end do
       do while (ageTableExpansionFactor(iTableTurnaround)<aExpansion)
          ageTableTimeMaximum=max(ageTableTimeMaximum*ageTableIncrementFactor,tCosmologicalTurnaround)
          call Make_Expansion_Factor_Table
       end do
    end if

    ! Interpolate to get cosmic time.
    Cosmology_Age_Matter_Lambda=Interpolate(ageTableNumberPoints,ageTableExpansionFactor&
         &,ageTableTime,interpolationObject,interpolationAccelerator&
         &,aExpansion,reset=resetInterpolation)
    if (collapsingPhaseActual) Cosmology_Age_Matter_Lambda=tCosmologicalMax-Cosmology_Age_Matter_Lambda

    return
  end function Cosmology_Age_Matter_Lambda

  double precision function Expansion_Factor_Matter_Lambda(tCosmological)
    !% Returns the expansion factor at cosmological time {\tt tCosmological}.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    double precision, intent(in) :: tCosmological
    double precision, save       :: tCosmologicalPrevious,expansionFactorPrevious
    !$omp threadprivate(tCosmologicalPrevious,expansionFactorPrevious)
    double precision             :: tEffective
    logical                      :: remakeTable

    ! Check if the time differs from the previous time.
    if (tCosmological /= tCosmologicalPrevious) then
       
       ! Quit on invalid input.
       if (tCosmological<0.0d0) call Galacticus_Error_Report('Expansion_Factor','cosmological time must be positive')
       
       ! Check if we need to recompute our table.
       if (ageTableInitialized) then
          remakeTable=(tCosmological<ageTableTime(1).or.tCosmological>ageTableTime(ageTableNumberPoints))
       else
          remakeTable=.true.
       end if
       if (remakeTable) call Make_Expansion_Factor_Table(tCosmological)
       
       ! Quit on invalid input.
       if (collapsingUniverse.and.tCosmological>tCosmologicalMax) call Galacticus_Error_Report('Expansion_Factor','cosmological time&
            & exceeds that at the Big Crunch')
       
       ! Interpolate to get the expansion factor.
       if (collapsingUniverse) then
          if (tCosmological <= tCosmologicalTurnaround) then
             tEffective=tCosmological
          else
             tEffective=tCosmologicalMax-tCosmological
          end if
       else
          tEffective=tCosmological
       end if
       expansionFactorPrevious=Interpolate(ageTableNumberPoints,ageTableTime,ageTableExpansionFactor,interpolationObjectInverse&
            &,interpolationAcceleratorInverse,tEffective,reset=resetInterpolationInverse)
       tCosmologicalPrevious=tCosmological
    end if
    
    ! Return the stored expansion factor.
    Expansion_Factor_Matter_Lambda=expansionFactorPrevious
    return
  end function Expansion_Factor_Matter_Lambda
  
  subroutine Make_Expansion_Factor_Table(tCosmological)
    !% Builds a table of expansion factor vs. time.
    use Numerical_Interpolation
    use Numerical_Ranges
    use ODE_Solver
    use Memory_Management
    use Array_Utilities
    implicit none    
    double precision, intent(in),  optional     :: tCosmological
    double precision, parameter                 :: odeToleranceAbsolute=1.0d-9, odeToleranceRelative=1.0d-9
    double precision, allocatable, dimension(:) :: ageTableTimeTemporary,ageTableExpansionFactorTemporary
    integer                                     :: iTime,prefixPointCount
    double precision                            :: densityPower,aDominant,Omega_Dominant,tDominant,time,aExpansion(1)
    type(c_ptr)                                 :: parameterPointer

    ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
    call Early_Time_Density_Scaling_Matter_Lambda(dominateFactor,densityPower,aDominant,Omega_Dominant)

    ! Find the corresponding time.
    tDominant=-2.0d0/densityPower/H_0_invGyr()/dsqrt(Omega_Dominant)/aDominant**(0.5d0*densityPower)

    ! Find minimum and maximum times to tabulate.
    if (present(tCosmological)) then
       time=tCosmological
       do while (ageTableTimeMinimum > min(time,tDominant)/2.0d0)
          ageTableTimeMinimum=ageTableTimeMinimum/ageTableIncrementFactor
       end do
       do while (ageTableTimeMaximum < max(time,tDominant)*2.0d0)
          ageTableTimeMaximum=ageTableTimeMaximum*ageTableIncrementFactor
       end do
    else
       do while (ageTableTimeMinimum > tDominant/2.0d0)
          ageTableTimeMinimum=ageTableTimeMinimum/ageTableIncrementFactor
       end do
       do while (ageTableTimeMaximum < tDominant*2.0d0)
          ageTableTimeMaximum=ageTableTimeMaximum*ageTableIncrementFactor
       end do
    end if
    if (collapsingUniverse) ageTableTimeMaximum=min(ageTableTimeMaximum,tCosmologicalTurnaround)

    ! Determine number of points to tabulate.
    ageTableNumberPoints=int(dlog10(ageTableTimeMaximum/ageTableTimeMinimum)*dble(ageTableNPointsPerDecade))+1
    ageTableTimeMaximum=ageTableTimeMinimum*10.0d0**(dble(ageTableNumberPoints)/dble(ageTableNPointsPerDecade))
    if (collapsingUniverse) ageTableTimeMaximum=min(ageTableTimeMaximum,tCosmologicalTurnaround)

    ! Deallocate arrays if currently allocated.
    if (allocated(ageTableTime)) then
       ! Determine number of points that are being added at the start of the array.
       prefixPointCount=int(dlog10(ageTableTime(1)/ageTableTimeMinimum)*dble(ageTableNPointsPerDecade)+0.5d0)
       call Move_Alloc(ageTableTime           ,ageTableTimeTemporary           )
       call Move_Alloc(ageTableExpansionFactor,ageTableExpansionFactorTemporary)
       ! Allocate the arrays to current required size.
       call Alloc_Array(ageTableTime,           [ageTableNumberPoints])
       call Alloc_Array(ageTableExpansionFactor,[ageTableNumberPoints])
       ! Create set of grid points in time variable.
       ageTableTime=Make_Range(ageTableTimeMinimum,ageTableTimeMaximum,ageTableNumberPoints,rangeTypeLogarithmic)
       ! Set the expansion factors to a negative value to indicate they are not yet computed.
       ageTableExpansionFactor=-1.0d0
       ! Paste in the previously computed regions.
       ageTableTime           (prefixPointCount+1:prefixPointCount+size(ageTableTimeTemporary))=ageTableTimeTemporary
       ageTableExpansionFactor(prefixPointCount+1:prefixPointCount+size(ageTableTimeTemporary))=ageTableExpansionFactorTemporary
       ! Deallocate the temporary arrays.
       call Dealloc_Array(ageTableTimeTemporary           )
       call Dealloc_Array(ageTableExpansionFactorTemporary)
    else    
       ! Allocate the arrays to current required size.
       call Alloc_Array(ageTableTime,           [ageTableNumberPoints])
       call Alloc_Array(ageTableExpansionFactor,[ageTableNumberPoints])
       ! Create set of grid points in time variable.
       ageTableTime=Make_Range(ageTableTimeMinimum,ageTableTimeMaximum,ageTableNumberPoints,rangeTypeLogarithmic)
       ! Set the expansion factors to a negative value to indicate they are not yet computed.
       ageTableExpansionFactor=-1.0d0
    end if

    ! For the initial time, we approximate that we are at sufficiently early times that a single component dominates the
    ! Universe and use the appropriate analytic solution.
    if (ageTableExpansionFactor(1) < 0.0d0) ageTableExpansionFactor(1)=(-0.5d0*densityPower*ageTableTime(1)*H_0_invGyr()&
         &*dsqrt(Omega_Dominant))**(-2.0d0/densityPower)

    ! Solve ODE to get corresponding expansion factors.
    iTableTurnaround=ageTableNumberPoints
    do iTime=2,ageTableNumberPoints
       ! Find the position in the table corresponding to turn around if we have a collapsing Universe.
       if (collapsingUniverse.and.ageTableTime(iTime-1)<tCosmologicalTurnaround.and.ageTableTime(iTime)&
            &>=tCosmologicalTurnaround) iTableTurnaround=iTime
       ! Compute the expansion factor if it is not already computed.
       if (ageTableExpansionFactor(iTime) < 0.0d0) then
          time=ageTableTime(iTime-1)
          aExpansion(1)=ageTableExpansionFactor(iTime-1)
          call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,time,ageTableTime(iTime),1,aExpansion&
               &,ageTableODEs,parameterPointer,odeToleranceAbsolute,odeToleranceRelative,reset=odeReset)
          ageTableExpansionFactor(iTime)=aExpansion(1)
       end if
    end do

    call Interpolate_Done(interpolationObject       ,interpolationAccelerator       ,resetInterpolation       )
    call Interpolate_Done(interpolationObjectInverse,interpolationAcceleratorInverse,resetInterpolationInverse)
    resetInterpolation       =.true.
    resetInterpolationInverse=.true.
    ! Flag that the table is now initialized.
    ageTableInitialized=.true.
    return
  end subroutine Make_Expansion_Factor_Table
  
  function ageTableODEs(t,a,dadt,parameterPointer) bind(c)
    !% System of differential equations to solve for expansion factor vs. age.
    integer(c_int)                           :: ageTableODEs
    real(c_double), value                    :: t
    real(c_double), dimension(1), intent(in) :: a
    real(c_double), dimension(1)             :: dadt
    type(c_ptr),    value                    :: parameterPointer

    dadt(1)=a(1)*Expansion_Rate_Matter_Lambda(a(1))
    ageTableODEs=FGSL_Success
  end function ageTableODEs
  
  double precision function Expansion_Rate_Matter_Lambda(aExpansion)
    !% Returns the cosmological expansion rate, $\dot{a}/a$ at expansion factor {\tt aExpansion}.
    implicit none
    double precision, intent(in) :: aExpansion
 
    ! Required value is simply the Hubble parameter but expressed in units of inverse Gyr.
    Expansion_Rate_Matter_Lambda=Hubble_Parameter_Matter_Lambda(aExpansion=aExpansion)*H_0_invGyr()/H_0()
    return
  end function Expansion_Rate_Matter_Lambda

  double precision function Hubble_Parameter_Matter_Lambda(tCosmological,aExpansion,collapsingPhase)
    !% Returns the Hubble parameter at the request cosmological time, {\tt tCosmological}, or expansion factor, {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: tCosmological,aExpansion
    logical,          intent(in), optional :: collapsingPhase
    double precision                       :: aExpansionActual,sqrtArgument
 
    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Hubble_Parameter_Matter_Lambda','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Lambda(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('Hubble_Parameter_Matter_Lambda','either a time or expansion factor must be specified')
       end if
    end if
    ! Compute the Hubble parameter at the specified expansion factor.
    sqrtArgument=max(Omega_Matter()/aExpansionActual**3+Omega_DE()+Omega_K()/aExpansionActual**2,0.0d0)
    Hubble_Parameter_Matter_Lambda=H_0()*dsqrt(sqrtArgument)
    ! Make the Hubble parameter negative if we are in the collapsing phase of the Universe.
    if (collapsingUniverse) then
       if (present(tCosmological)) then
          if (tCosmological>tCosmologicalTurnaround) Hubble_Parameter_Matter_Lambda=-Hubble_Parameter_Matter_Lambda
       else
          if (present(collapsingPhase)) then
             if (collapsingPhase) Hubble_Parameter_Matter_Lambda=-Hubble_Parameter_Matter_Lambda
          end if
       end if
    end if
    return
  end function Hubble_Parameter_Matter_Lambda
  
  subroutine Early_Time_Density_Scaling_Matter_Lambda(dominateFactor,densityPower,aDominant,Omega_Dominant)
    use Cosmological_Parameters
    implicit none
    double precision, intent(in)            :: dominateFactor
    double precision, intent(out)           :: densityPower,aDominant
    double precision, intent(out), optional :: Omega_Dominant

    ! For matter and cosmological constant, matter always dominates at early times.
    densityPower=-3.0d0 ! Power-law scaling of matter density with expansion factor.

    ! Choose present day as default - will be used if no other densities present (i.e. Einsetin-de Sitter).
    aDominant=Epoch_of_Matter_Domination_Matter_Lambda(dominateFactor)

    ! Return the density parameter in the dominant species if required.
    if (present(Omega_Dominant)) Omega_Dominant=Omega_Matter()
    return
  end subroutine Early_Time_Density_Scaling_Matter_Lambda

  double precision function Epoch_of_Matter_Domination_Matter_Lambda(dominateFactor)
    use Cosmological_Parameters
    use Cosmology_Functions_Parameters
    implicit none
    double precision, intent(in)            :: dominateFactor
    double precision                        :: aMatterEquality,aDominantDarkEnergy,aDominantCurvature

    ! Choose present day as default - will be used if no other densities present (i.e. Einsetin-de Sitter).
    Epoch_of_Matter_Domination_Matter_Lambda=1.0d0

    if (Omega_DE()/=0.0d0) then
       ! Find the expansion factor of matter-dark energy equality.
       aMatterEquality=Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda(requestTypeExpansionFactor)
       
       ! Find the earlier expansion factor at which matter dominates by the specified amount (ratio of matter
       ! to dark energy density scales as the cube of expansion factor).
       aDominantDarkEnergy=aMatterEquality/dominateFactor**(1.0d0/3.0d0)

       ! Choose earliest expansion factor.
       Epoch_of_Matter_Domination_Matter_Lambda=min(Epoch_of_Matter_Domination_Matter_Lambda,aDominantDarkEnergy)
    end if

    if (Omega_K()/=0.0d0) then
       ! Find the expansion factor of matter-curvature equality.
       aMatterEquality=Epoch_of_Matter_Curvature_Equality_Matter_Lambda(requestTypeExpansionFactor)
       
       ! Find the earlier expansion factor at which matter dominates by the specified amount (ratio of matter
       ! to curvature density scales as the expansion factor).
       aDominantCurvature=aMatterEquality/dominateFactor

       ! Choose earliest expansion factor.
       Epoch_of_Matter_Domination_Matter_Lambda=min(Epoch_of_Matter_Domination_Matter_Lambda,aDominantCurvature)
    end if
    return
  end function Epoch_of_Matter_Domination_Matter_Lambda

  double precision function Omega_Matter_Total_Matter_Lambda(tCosmological,aExpansion,collapsingPhase)
    !% Return the matter density parameter at expansion factor {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: tCosmological,aExpansion
    logical,          intent(in), optional :: collapsingPhase
    double precision                       :: aExpansionActual

    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Omega_Matter_Total_Matter_Lambda','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Lambda(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('Omega_Matter_Total_Matter_Lambda','either a time or expansion factor must be specified')
       end if
    end if
    Omega_Matter_Total_Matter_Lambda=Omega_Matter()*((H_0()/Hubble_Parameter_Matter_Lambda(aExpansion=aExpansionActual))**2)/(aExpansionActual**3)
    return
  end function Omega_Matter_Total_Matter_Lambda

  double precision function Omega_Dark_Energy_Matter_Lambda(tCosmological,aExpansion,collapsingPhase)
    !% Return the dark energy density parameter at expansion factor {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: tCosmological,aExpansion
    logical,          intent(in), optional :: collapsingPhase
    double precision                       :: aExpansionActual
    
    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Omega_Dark_Energy_Matter_Lambda','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Lambda(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('Omega_Dark_Energy_Matter_Lambda','either a time or expansion factor must be specified')
       end if
    end if 
    Omega_Dark_Energy_Matter_Lambda=Omega_DE()*((H_0()/Hubble_Parameter_Matter_Lambda(aExpansion=aExpansionActual))**2)
    return
  end function Omega_Dark_Energy_Matter_Lambda

  double precision function CMB_Temperature_Matter_Lambda(tCosmological,aExpansion,collapsingPhase)
    !% Return the temperature of the CMB at expansion factor {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: tCosmological,aExpansion
    logical,          intent(in), optional :: collapsingPhase
    double precision                       :: aExpansionActual

    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('CMB_Temperature_Matter_Lambda','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Lambda(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('CMB_Temperature_Matter_Lambda','either a time or expansion factor must be specified')
       end if
    end if
    CMB_Temperature_Matter_Lambda=T_CMB()/aExpansionActual
    return
  end function CMB_Temperature_Matter_Lambda

  double precision function Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda(requestType)
    !% Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    use Cosmology_Functions_Parameters
    implicit none
    integer, intent(in), optional :: requestType
    integer                       :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if

    Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda=(Omega_Matter()/dabs(Omega_DE()))**(1.0d0/3.0d0)
    if (requestType == requestTypeTime) Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda&
         &=Cosmology_Age_Matter_Lambda(Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda)
    return
  end function Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda

  double precision function Epoch_of_Matter_Curvature_Equality_Matter_Lambda(requestType)
    !% Return the epoch of matter-curvature magnitude equality (either expansion factor or cosmic time).
    use Cosmology_Functions_Parameters
    implicit none
    integer, intent(in), optional :: requestType
    integer                       :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if

    Epoch_of_Matter_Curvature_Equality_Matter_Lambda=Omega_Matter()/dabs(Omega_K())
    if (requestType == requestTypeTime) Epoch_of_Matter_Curvature_Equality_Matter_Lambda&
         &=Cosmology_Age_Matter_Lambda(Epoch_of_Matter_Curvature_Equality_Matter_Lambda)
    return
  end function Epoch_of_Matter_Curvature_Equality_Matter_Lambda

  double precision function Time_From_Comoving_Distance_Matter_Lambda(comovingDistance)
    !% Returns the cosmological time corresponding to given {\tt comovingDistance}.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    double precision, intent(in) :: comovingDistance
    double precision             :: tCosmological
    logical                      :: remakeTable

    ! Quit on invalid input.
    if (comovingDistance < 0.0d0) call Galacticus_Error_Report('Time_From_Comoving_Distance_Matter_Lambda','comoving distance must be positive')

    !$omp critical(Cosmology_Functions_Matter_Lambda_Distance_Initialize)
    ! Check if we need to recompute our table.
    remakeTable=.true.
    do while (remakeTable)
       if (distanceTableInitialized) then
          remakeTable=distanceTableComovingDistance(1) < comovingDistance
          tCosmological=0.5d0*distanceTableTime(1)
       else
          remakeTable=.true.
          tCosmological=distanceTableTimeMinimum
       end if
       ! Remake table if necessary.
       if (remakeTable) call Make_Distance_Table(tCosmological)
    end do
    !$omp end critical(Cosmology_Functions_Matter_Lambda_Distance_Initialize)

    ! Interpolate to get the comoving distance.
    Time_From_Comoving_Distance_Matter_Lambda=Interpolate(distanceTableNumberPoints,distanceTableComovingDistanceNegated&
         &,distanceTableTime ,interpolationObjectDistanceInverse ,interpolationAcceleratorDistanceInverse,-comovingDistance,reset &
         &=resetInterpolationDistanceInverse)
    return
  end function Time_From_Comoving_Distance_Matter_Lambda

  double precision function Comoving_Distance_Matter_Lambda(tCosmological)
    !% Returns the comoving distance to cosmological time {\tt tCosmological}.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    double precision, intent(in) :: tCosmological
    logical                      :: remakeTable

    ! Quit on invalid input.
    if (tCosmological < 0.0d0                             ) call Galacticus_Error_Report('Comoving_Distance_Matter_Lambda','cosmological time must be positive'   )
    if (tCosmological > Cosmology_Age_Matter_Lambda(1.0d0)) call Galacticus_Error_Report('Comoving_Distance_Matter_Lambda','cosmological time must be in the past')

    !$omp critical(Cosmology_Functions_Matter_Lambda_Distance_Initialize)
    ! Check if we need to recompute our table.
    if (distanceTableInitialized) then
       remakeTable=(tCosmological<distanceTableTime(1).or.tCosmological>distanceTableTime(distanceTableNumberPoints))
    else
       remakeTable=.true.
    end if
    if (remakeTable) call Make_Distance_Table(tCosmological)
    !$omp end critical(Cosmology_Functions_Matter_Lambda_Distance_Initialize)

    ! Quit on invalid input.
    if (collapsingUniverse.and.tCosmological>tCosmologicalMax) call Galacticus_Error_Report('Expansion_Factor','cosmological time&
         & exceeds that at the Big Crunch')

    ! Interpolate to get the comoving distance.
    Comoving_Distance_Matter_Lambda=Interpolate(distanceTableNumberPoints,distanceTableTime,distanceTableComovingDistance,interpolationObjectDistance&
         &,interpolationAcceleratorDistance,tCosmological,reset=resetInterpolationDistance)
    return
  end function Comoving_Distance_Matter_Lambda
  
  double precision function Comoving_Distance_Conversion_Matter_Lambda(output,distanceModulus,redshift)
    !% Convert bewteen different measures of distance.
    use Numerical_Interpolation
    use Galacticus_Error
    use Cosmology_Functions_Options
    implicit none
    integer         , intent(in)           :: output
    double precision, intent(in), optional :: distanceModulus,redshift
    logical                                :: remakeTable,gotComovingDistance
    double precision                       :: luminosityDistance,comovingDistance

    !$omp critical(Cosmology_Functions_Matter_Lambda_Distance_Initialize)
    ! Check if we need to recompute our table.
    if (.not.distanceTableInitialized) call Make_Distance_Table(Cosmology_Age_Matter_Lambda(1.0d0))
    !$omp end critical(Cosmology_Functions_Matter_Lambda_Distance_Initialize)

    ! Convert to comoving distance from whatever was supplied.
    gotComovingDistance=.false.
    if (present(distanceModulus)) then
       luminosityDistance =10.0d0**((distanceModulus-25.0d0)/5.0d0)
       !$omp critical(Cosmology_Functions_Matter_Lambda_Distance_Initialize)
       do while (luminosityDistance > -distanceTableLuminosityDistanceNegated(1))
          call Make_Distance_Table(0.5d0*distanceTableTimeMinimum) 
       end do
       !$omp end critical(Cosmology_Functions_Matter_Lambda_Distance_Initialize)
       comovingDistance   =-Interpolate(distanceTableNumberPoints,distanceTableLuminosityDistanceNegated &
            &,distanceTableComovingDistanceNegated,interpolationObjectLuminosityDistance&
            &,interpolationAcceleratorLuminosityDistance,-luminosityDistance,reset =resetInterpolationLuminosityDistance)
       gotComovingDistance=.true.
    end if
    if (present(redshift)) then
       comovingDistance   =Comoving_Distance_Matter_Lambda(Cosmology_Age_Matter_Lambda(1.0d0/(1.0d0+redshift)))
       gotComovingDistance=.true.
    end if
    if (.not.gotComovingDistance) call Galacticus_Error_Report('Comoving_Distance_Conversion_Matter_Lambda','no distance measure&
         & provided')
    
    ! Convert to required distance measure.
    select case (output)
    case (distanceTypeComoving)
       Comoving_Distance_Conversion_Matter_Lambda=comovingDistance
    case default
       call Galacticus_Error_Report('Comoving_Distance_Conversion_Matter_Lambda','unrecognized output option')
    end select

    return
  end function Comoving_Distance_Conversion_Matter_Lambda

  subroutine Make_Distance_Table(tCosmological)
    !% Builds a table of distance vs. time.
    use Numerical_Interpolation
    use Numerical_Ranges
    use Numerical_Integration
    use Memory_Management
    use Array_Utilities
    implicit none
    double precision,                intent(in) :: tCosmological
    double precision,                parameter  :: toleranceAbsolute=1.0d-5, toleranceRelative=1.0d-5
    integer                                     :: iTime
    logical                                     :: resetIntegration
    type(c_ptr)                                 :: parameterPointer
    type(fgsl_function)                         :: integrandFunction
    type(fgsl_integration_workspace)            :: integrationWorkspace

    ! Find minimum and maximum times to tabulate.
    distanceTableTimeMinimum=min(distanceTableTimeMinimum,0.5d0*tCosmological)
    distanceTableTimeMaximum=Cosmology_Age_Matter_Lambda(1.0d0)
 
    ! Determine number of points to tabulate.
    distanceTableNumberPoints=int(dlog10(distanceTableTimeMaximum/distanceTableTimeMinimum)*dble(distanceTableNPointsPerDecade))+1
 
    ! Deallocate arrays if currently allocated.
    if (allocated(distanceTableTime                     )) call Dealloc_Array(distanceTableTime                     )
    if (allocated(distanceTableComovingDistance         )) call Dealloc_Array(distanceTableComovingDistance         )
    if (allocated(distanceTableComovingDistanceNegated  )) call Dealloc_Array(distanceTableComovingDistanceNegated  )
    if (allocated(distanceTableLuminosityDistanceNegated)) call Dealloc_Array(distanceTableLuminosityDistanceNegated)
    ! Allocate the arrays to current required size.
    call Alloc_Array(distanceTableTime                     ,[distanceTableNumberPoints])
    call Alloc_Array(distanceTableComovingDistance         ,[distanceTableNumberPoints])
    call Alloc_Array(distanceTableComovingDistanceNegated  ,[distanceTableNumberPoints])
    call Alloc_Array(distanceTableLuminosityDistanceNegated,[distanceTableNumberPoints])
    
    ! Create the range of times.
    distanceTableTime=Make_Range(distanceTableTimeMinimum,distanceTableTimeMaximum,distanceTableNumberPoints,rangeTypeLogarithmic)

    ! Integrate to get the comoving distance.
    resetIntegration=.true.
    do iTime=1,distanceTableNumberPoints
       distanceTableComovingDistance(iTime)=Integrate(distanceTableTime(iTime),distanceTableTime(distanceTableNumberPoints) &
            &,Comoving_Distance_Integrand,parameterPointer,integrandFunction ,integrationWorkspace,toleranceAbsolute&
            &=toleranceAbsolute ,toleranceRelative=toleranceRelative,reset=resetIntegration)
       distanceTableLuminosityDistanceNegated(iTime)=distanceTableComovingDistance(iTime)&
            &/Expansion_Factor_Matter_Lambda(distanceTableTime(iTime))
    end do
    ! Make a negated copy of the distances so that we have an increasing array for use in interpolation routines.
    distanceTableComovingDistanceNegated  =-distanceTableComovingDistance
    distanceTableLuminosityDistanceNegated=-distanceTableLuminosityDistanceNegated
    ! Reset interpolators.
    call Interpolate_Done(interpolationObjectDistance          ,interpolationAcceleratorDistance          ,resetInterpolationDistance          )
    call Interpolate_Done(interpolationObjectDistanceInverse   ,interpolationAcceleratorDistanceInverse   ,resetInterpolationDistanceInverse   )
    call Interpolate_Done(interpolationObjectLuminosityDistance,interpolationAcceleratorLuminosityDistance,resetInterpolationLuminosityDistance)
    resetInterpolationDistance          =.true.
    resetInterpolationDistanceInverse   =.true.
    resetInterpolationLuminosityDistance=.true.

    ! Flag that the table is now initialized.
    distanceTableInitialized=.true.
    return
  end subroutine Make_Distance_Table
  
  function Comoving_Distance_Integrand(time,parameterPointer) bind(c)
    !% Integrand function used in computing the comoving distance.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    real(c_double)          :: Comoving_Distance_Integrand
    real(c_double), value   :: time
    type(c_ptr),    value   :: parameterPointer

    Comoving_Distance_Integrand=speedLight*gigaYear/megaParsec/Expansion_Factor_Matter_Lambda(time)
    return
  end function Comoving_Distance_Integrand

  !# <galacticusStateStoreTask>
  !#  <unitName>Cosmology_Matter_Lambda_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Cosmology_Matter_Lambda_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Store the full tables, as they are hysteretic and cannot be reconstructed precisely without knowing the path by which they
    ! were originally constructed.
    write (stateFile) ageTableNumberPoints,ageTableTimeMinimum,ageTableTimeMaximum
    write (stateFile) ageTableTime,ageTableExpansionFactor
    write (stateFile) distanceTableNumberPoints,distanceTableTimeMinimum,distanceTableTimeMaximum
    write (stateFile) distanceTableTime,distanceTableComovingDistance,distanceTableComovingDistanceNegated
    return
  end subroutine Cosmology_Matter_Lambda_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Cosmology_Matter_Lambda_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Cosmology_Matter_Lambda_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use Memory_Management
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the tabulations.
    read (stateFile) ageTableNumberPoints,ageTableTimeMinimum,ageTableTimeMaximum
    if (allocated(ageTableTime           )) call Dealloc_Array(ageTableTime           )
    if (allocated(ageTableExpansionFactor)) call Dealloc_Array(ageTableExpansionFactor)
    call Alloc_Array(ageTableTime           ,[ageTableNumberPoints])
    call Alloc_Array(ageTableExpansionFactor,[ageTableNumberPoints])
    read (stateFile) ageTableTime,ageTableExpansionFactor
    read (stateFile) distanceTableNumberPoints,distanceTableTimeMinimum,distanceTableTimeMaximum
    if (allocated(distanceTableTime                   )) call Dealloc_Array(distanceTableTime                   )
    if (allocated(distanceTableComovingDistance       )) call Dealloc_Array(distanceTableComovingDistance       )
    if (allocated(distanceTableComovingDistanceNegated)) call Dealloc_Array(distanceTableComovingDistanceNegated)
    call Alloc_Array(distanceTableTime                   ,[distanceTableNumberPoints])
    call Alloc_Array(distanceTableComovingDistance       ,[distanceTableNumberPoints])
    call Alloc_Array(distanceTableComovingDistanceNegated,[distanceTableNumberPoints])
    read (stateFile) distanceTableTime,distanceTableComovingDistance,distanceTableComovingDistanceNegated
    
    ! Ensure that interpolation objects will get reset.
    resetInterpolation               =.true.
    resetInterpolationInverse        =.true.
    resetInterpolationDistance       =.true.
    resetInterpolationDistanceInverse=.true.

    return
  end subroutine Cosmology_Matter_Lambda_State_Retrieve
  
end module Cosmology_Functions_Matter_Lambda
