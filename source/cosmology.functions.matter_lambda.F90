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


!% Contains a module which implements useful cosmological functions. This implementation assumes a Universe filled with
!% collisionless matter and a cosmological constant.

module Cosmology_Functions_Matter_Lambda
  !% Implements useful cosmological functions. This implementation assumes a Universe filled with collisionless matter
  !% and a cosmological constant.
  use Cosmological_Parameters
  use FGSL
  use, intrinsic :: ISO_C_Binding
  private
  public :: Cosmology_Functions_Matter_Lambda_Initialize, Cosmology_Matter_Lambda_State_Store,&
       & Cosmology_Matter_Lambda_State_Retrieve

  ! Flag indicating if the module has been initialized.
  logical                                     :: moduleInitialized=.false.
  
  ! Variables used to track critical times in collapsing Universes.
  logical                                     :: collapsingUniverse=.false.
  integer                                     :: iTableTurnaround
  double precision                            :: tCosmologicalMax,tCosmologicalTurnaround,aExpansionMax

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
  type(fgsl_interp)                           :: interpolationObject
  type(fgsl_interp_accel)                     :: interpolationAccelerator
  logical                                     :: resetInterpolation=.true.

  ! Variables to hold table of distance vs. cosmic time.
  logical                                     :: distanceTableInitialized=.false.
  integer                                     :: distanceTableNumberPoints
  double precision                            :: distanceTableTimeMinimum=1.0d-4, distanceTableTimeMaximum
  integer,          parameter                 :: distanceTableNPointsPerDecade=100
  double precision, allocatable, dimension(:) :: distanceTableTime,distanceTableComovingDistance,distanceTableComovingDistanceNegated
  type(fgsl_interp)                           :: interpolationObjectDistance      ,interpolationObjectDistanceInverse
  type(fgsl_interp_accel)                     :: interpolationAcceleratorDistance ,interpolationAcceleratorDistanceInverse
  logical                                     :: resetInterpolationDistance=.true.,resetInterpolationDistanceInverse=.true.

  ! Variables used in the ODE solver.
  type(fgsl_odeiv_step)                       :: odeStepper
  type(fgsl_odeiv_control)                    :: odeController
  type(fgsl_odeiv_evolve)                     :: odeEvolver
  type(fgsl_odeiv_system)                     :: odeSystem
  logical                                     :: odeReset=.true. ! Ensure ODE variables will be reset on first call.
  
contains

  !# <cosmologyMethod>
  !#  <unitName>Cosmology_Functions_Matter_Lambda_Initialize</unitName>
  !# </cosmologyMethod>
  subroutine Cosmology_Functions_Matter_Lambda_Initialize(cosmologyMethod,Expansion_Factor_Is_Valid_Get,Cosmic_Time_Is_Valid_Get &
       &,Cosmology_Age_Get,Expansion_Factor_Get,Hubble_Parameter_Get,Early_Time_Density_Scaling_Get,Omega_Matter_Get &
       &,Omega_Dark_Energy_Get,Expansion_Rate_Get,Epoch_of_Matter_Dark_Energy_Equality_Get,Epoch_of_Matter_Domination_Get&
       &,Epoch_of_Matter_Curvature_Equality_Get,CMB_Temperature_Get,Comoving_Distance_Get,Time_From_Comoving_Distance_Get)
    !% Initialize the module.
    use Numerical_Comparison
    use ISO_Varying_String
    use ODE_Solver
    implicit none
    type(varying_string),                 intent(in)    :: cosmologyMethod
    procedure(),                 pointer, intent(inout) :: Early_Time_Density_Scaling_Get
    procedure(logical),          pointer, intent(inout) :: Expansion_Factor_Is_Valid_Get,Cosmic_Time_Is_Valid_Get
    procedure(double precision), pointer, intent(inout) :: Cosmology_Age_Get ,Expansion_Factor_Get,Hubble_Parameter_Get &
         &,Omega_Matter_Get,Omega_Dark_Energy_Get ,Expansion_Rate_Get,Epoch_of_Matter_Dark_Energy_Equality_Get &
         &,Epoch_of_Matter_Domination_Get ,Epoch_of_Matter_Curvature_Equality_Get,CMB_Temperature_Get,Comoving_Distance_Get&
         &,Time_From_Comoving_Distance_Get
    double precision,            parameter              :: odeToleranceAbsolute=1.0d-9, odeToleranceRelative=1.0d-9
    double precision,            parameter              :: omegaTolerance=1.0d-9
    double precision                                    :: cubicTerm1,cubicTerm5,cubicTerm9,cubicTerm21Squared,cubicTerm21 &
         &,cubicTerm25Cubed,cubicTerm25,aMaximum,aDominant,timeMaximum(1),densityPower,Omega_Dominant
    type(c_ptr)                                         :: parameterPointer

    ! Check if our method is selected.
    if (cosmologyMethod == 'matter + lambda') then
       ! Set up procedure pointers.
       Expansion_Factor_Is_Valid_Get            => Expansion_Factor_Is_Valid_Matter_Lambda
       Cosmic_Time_Is_Valid_Get                 => Cosmic_Time_Is_Valid_Matter_Lambda
       Cosmology_Age_Get                        => Cosmology_Age_Matter_Lambda
       Expansion_Factor_Get                     => Expansion_Factor_Matter_Lambda
       Hubble_Parameter_Get                     => Hubble_Parameter_Matter_Lambda
       Early_Time_Density_Scaling_Get           => Early_Time_Density_Scaling_Matter_Lambda
       Omega_Matter_Get                         => Omega_Matter_Matter_Lambda
       Omega_Dark_Energy_Get                    => Omega_Dark_Energy_Matter_Lambda
       Expansion_Rate_Get                       => Expansion_Rate_Matter_Lambda
       Epoch_of_Matter_Dark_Energy_Equality_Get => Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda
       Epoch_of_Matter_Curvature_Equality_Get   => Epoch_of_Matter_Curvature_Equality_Matter_Lambda
       Epoch_of_Matter_Domination_Get           => Epoch_of_Matter_Domination_Matter_Lambda
       CMB_Temperature_Get                      => CMB_Temperature_Matter_Lambda
       Comoving_Distance_Get                    => Comoving_Distance_Matter_Lambda
       Time_From_Comoving_Distance_Get          => Time_From_Comoving_Distance_Matter_Lambda

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
                aExpansionMax=-(Omega_0()*Omega_DE()**2)**(1.0d0/3.0d0)/Omega_DE()
             end if
          end if
       else
          if (Values_Agree(Omega_DE(),0.0d0,absTol=omegaTolerance)) then
             ! Simple case for a matter-only universe.
             collapsingUniverse=Omega_0() > 1.0d0
             if (collapsingUniverse) aExpansionMax=Omega_0()/(Omega_0()-1.0d0)
          else
             ! Case of matter plus dark energy.
             cubicTerm1 =1.0d0/Omega_DE()
             cubicTerm5 =Omega_0()**2
             cubicTerm9 =Omega_DE()**2
             cubicTerm21Squared=-(-0.12d2+0.36d2*Omega_0()+36.0d0*Omega_DE()-0.36d2*cubicTerm5-0.72d2*Omega_0()*Omega_DE()-36.0d0&
                  &*cubicTerm9+0.12d2*cubicTerm5*Omega_0()-0.45d2*cubicTerm5*Omega_DE()+0.36d2*Omega_0()*cubicTerm9+12.0d0*cubicTerm9&
                  &*Omega_DE())*cubicTerm1
             if (cubicTerm21Squared > 0.0d0) then
                cubicTerm21=dsqrt(cubicTerm21Squared)
                cubicTerm25Cubed=(-0.108d3*Omega_0()+0.12d2*cubicTerm21)*cubicTerm9
                if (cubicTerm25Cubed >= 0.0d0) then
                   cubicTerm25=cubicTerm25Cubed**(1.0d0/3.0d0)
                else
                   cubicTerm25=-dabs(cubicTerm25Cubed)**(1.0d0/3.0d0)
                end if
                aMaximum=cubicTerm1*cubicTerm25/0.6d1+0.2d1*(-0.1d1+Omega_0()+Omega_DE())/cubicTerm25
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
          timeMaximum(1)=1.0/H_0_invGyr()/dsqrt(Omega_Dominant)/aDominant**(0.5d0*densityPower)
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
    !$omp critical(Cosmology_Functions_Initialize)
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
    !$omp end critical(Cosmology_Functions_Initialize)

    ! Interpolate to get cosmic time.
    !$omp critical(Cosmology_Functions_Interpolate)
    Cosmology_Age_Matter_Lambda=Interpolate(ageTableNumberPoints,ageTableExpansionFactor&
         &,ageTableTime,interpolationObject,interpolationAccelerator&
         &,aExpansion,reset=resetInterpolation)
    !$omp end critical(Cosmology_Functions_Interpolate)
    if (collapsingPhaseActual) Cosmology_Age_Matter_Lambda=tCosmologicalMax-Cosmology_Age_Matter_Lambda

    return
  end function Cosmology_Age_Matter_Lambda

  double precision function Expansion_Factor_Matter_Lambda(tCosmological)
    !% Returns the expansion factor at cosmological time {\tt tCosmological}.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    double precision, intent(in) :: tCosmological
    double precision             :: tEffective
    logical                      :: remakeTable

    ! Quit on invalid input.
    if (tCosmological<0.0d0) call Galacticus_Error_Report('Expansion_Factor','cosmological time must be positive')

    !$omp critical(Cosmology_Functions_Initialize)
    ! Check if we need to recompute our table.
    if (ageTableInitialized) then
       remakeTable=(tCosmological<ageTableTime(1).or.tCosmological>ageTableTime(ageTableNumberPoints))
    else
       remakeTable=.true.
    end if
    if (remakeTable) call Make_Expansion_Factor_Table(tCosmological)
    !$omp end critical(Cosmology_Functions_Initialize)

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
    !$omp critical(Cosmology_Functions_Interpolate)
    Expansion_Factor_Matter_Lambda=Interpolate(ageTableNumberPoints,ageTableTime,ageTableExpansionFactor,interpolationObject&
         &,interpolationAccelerator,tEffective,reset=resetInterpolation)
    !$omp end critical(Cosmology_Functions_Interpolate)
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

    !$omp critical(Cosmology_Functions_Interpolate)
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
    call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
    resetInterpolation=.true.
    
    ! Flag that the table is now initialized.
    ageTableInitialized=.true.
    !$omp end critical(Cosmology_Functions_Interpolate)
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
    sqrtArgument=max(Omega_0()/aExpansionActual**3+Omega_DE()+Omega_K()/aExpansionActual**2,0.0d0)
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
    if (present(Omega_Dominant)) Omega_Dominant=Omega_0()
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

  double precision function Omega_Matter_Matter_Lambda(tCosmological,aExpansion,collapsingPhase)
    !% Return the matter density parameter at expansion factor {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: tCosmological,aExpansion
    logical,          intent(in), optional :: collapsingPhase
    double precision                       :: aExpansionActual

    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Omega_Matter_Matter_Lambda','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Lambda(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('Omega_Matter_Matter_Lambda','either a time or expansion factor must be specified')
       end if
    end if
    Omega_Matter_Matter_Lambda=Omega_0()*((H_0()/Hubble_Parameter_Matter_Lambda(aExpansion=aExpansionActual))**2)/(aExpansionActual**3)
    return
  end function Omega_Matter_Matter_Lambda

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

    Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda=(Omega_0()/dabs(Omega_DE()))**(1.0d0/3.0d0)
    if (requestType.eq.requestTypeTime) Epoch_of_Matter_Dark_Energy_Equality_Matter_Lambda&
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

    Epoch_of_Matter_Curvature_Equality_Matter_Lambda=Omega_0()/dabs(Omega_K())
    if (requestType.eq.requestTypeTime) Epoch_of_Matter_Curvature_Equality_Matter_Lambda&
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
    if (allocated(distanceTableTime                   )) call Dealloc_Array(distanceTableTime                   )
    if (allocated(distanceTableComovingDistance       )) call Dealloc_Array(distanceTableComovingDistance       )
    if (allocated(distanceTableComovingDistanceNegated)) call Dealloc_Array(distanceTableComovingDistanceNegated)
    ! Allocate the arrays to current required size.
    call Alloc_Array(distanceTableTime                   ,[distanceTableNumberPoints])
    call Alloc_Array(distanceTableComovingDistance       ,[distanceTableNumberPoints])
    call Alloc_Array(distanceTableComovingDistanceNegated,[distanceTableNumberPoints])
    
    ! Create the range of times.
    distanceTableTime=Make_Range(distanceTableTimeMinimum,distanceTableTimeMaximum,distanceTableNumberPoints,rangeTypeLogarithmic)

    ! Integrate to get the comoving distance.
    resetIntegration=.true.
    do iTime=1,distanceTableNumberPoints
       distanceTableComovingDistance(iTime)=Integrate(distanceTableTime(iTime),distanceTableTime(distanceTableNumberPoints) &
            &,Comoving_Distance_Integrand,parameterPointer,integrandFunction ,integrationWorkspace,toleranceAbsolute&
            &=toleranceAbsolute ,toleranceRelative=toleranceRelative,reset=resetIntegration)
    end do
    ! Make a negated copy of the distances so that we have an increasing array for use in interpolation routines.
    distanceTableComovingDistanceNegated=-distanceTableComovingDistance
    ! Reset interpolators.
    call Interpolate_Done(interpolationObjectDistance       ,interpolationAcceleratorDistance       ,resetInterpolationDistance       )
    call Interpolate_Done(interpolationObjectDistanceInverse,interpolationAcceleratorDistanceInverse,resetInterpolationDistanceInverse)
    resetInterpolationDistance       =.true.
    resetInterpolationDistanceInverse=.true.
    
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

    write (stateFile) ageTableTimeMinimum,ageTableTimeMaximum
    return
  end subroutine Cosmology_Matter_Lambda_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Cosmology_Matter_Lambda_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Cosmology_Matter_Lambda_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile
    double precision            :: tCosmological

    ! Read the minimum and maximum tabulated times.
    read (stateFile) ageTableTimeMinimum,ageTableTimeMaximum
    ! Rebuild the table, choosing a time that won't affect the minimum/maximum tabulated.
    tCosmological=ageTableTimeMinimum*2.0d0
    call Make_Expansion_Factor_Table(tCosmological)
    return
  end subroutine Cosmology_Matter_Lambda_State_Retrieve
  
end module Cosmology_Functions_Matter_Lambda
