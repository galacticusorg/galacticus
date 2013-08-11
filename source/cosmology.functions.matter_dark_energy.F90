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
!% collisionless matter and dark energy with an equation of state of the form: $P=\rho^w$ with $w(a)=w_0+w_1 a
!% (1-a)$.

module Cosmology_Functions_Matter_Dark_Energy
  !% Implements useful cosmological functions. This implementation assumes a Universe filled with collisionless matter
  !% and dark energy with an equation of state of the form: $P=\rho^w$ with $w(a)=w_0+w_1 a (1-a)$.
  use Cosmological_Parameters
  use FGSL
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: Cosmology_Functions_Matter_Dark_Energy_Initialize, Cosmology_Matter_Dark_Energy_State_Store,&
       & Cosmology_Matter_Dark_Energy_State_Retrieve

  ! Variables used to track critical times in collapsing Universes.
  logical                                                         :: collapsingUniverse         =.false.
  integer                                                         :: iTableTurnaround
  double precision                                                :: aExpansionMax                                                                                                  , tCosmologicalMax                                                                            , &
       &                                                             tCosmologicalTurnaround
  !$omp threadprivate(collapsingUniverse,iTableTurnaround,aExpansionMax,tCosmologicalMax,tCosmologicalTurnaround)
  ! Dark energy equation of state.
  double precision                                                :: darkEnergyEquationOfStateW0                                                                                    , darkEnergyEquationOfStateW1

  ! Variables used in root finding.
  double precision                                                :: dominateFactorCurrent

  ! Variables to hold table of expansion factor vs. cosmic time.
  logical                                                         :: ageTableInitialized        =.false.
  integer                                                         :: ageTableNumberPoints
  double precision                                                :: ageTableTimeMaximum        =20.0d0                                                                             , ageTableTimeMinimum                                                                  =1.0d-4
  integer                             , parameter                 :: ageTableNPointsPerDecade   =300
  double precision                    , parameter                 :: ageTableNPointsPerOctave   =dble(ageTableNPointsPerDecade)*log(2.0d0)/log(10.0d0)
  double precision                    , parameter                 :: ageTableIncrementFactor    =exp(int(ageTableNPointsPerOctave+1.0d0)*log(10.0d0)/dble(ageTableNPointsPerDecade))
  double precision                    , allocatable, dimension(:) :: ageTableExpansionFactor                                                                                        , ageTableTime
  type            (fgsl_interp       )                            :: interpolationObject                                                                                            , interpolationObjectInverse
  type            (fgsl_interp_accel )                            :: interpolationAccelerator                                                                                       , interpolationAcceleratorInverse
  logical                                                         :: resetInterpolation         =.true.
  logical                                                         :: resetInterpolationInverse  =.true.
  !$omp threadprivate(ageTableInitialized,ageTableNumberPoints,ageTableTimeMinimum,ageTableTimeMaximum,ageTableTime)
  !$omp threadprivate(ageTableExpansionFactor,interpolationObject,interpolationObjectInverse,interpolationAccelerator)
  !$omp threadprivate(interpolationAcceleratorInverse,resetInterpolation,resetInterpolationInverse)
  ! Variables used in the ODE solver.
  type            (fgsl_odeiv_step   )                            :: odeStepper
  type            (fgsl_odeiv_control)                            :: odeController
  type            (fgsl_odeiv_evolve )                            :: odeEvolver
  type            (fgsl_odeiv_system )                            :: odeSystem
  logical                                                         :: odeReset                   =.true.                                                                                                               !   Ensure ODE variables will be reset on first call.
  !$omp threadprivate (odeStepper,odeController,odeEvolver,odeSystem,odeReset)
contains

  !# <cosmologyMethod>
  !#  <unitName>Cosmology_Functions_Matter_Dark_Energy_Initialize</unitName>
  !# </cosmologyMethod>
  subroutine Cosmology_Functions_Matter_Dark_Energy_Initialize(cosmologyMethod,Expansion_Factor_Is_Valid_Get,Cosmic_Time_Is_Valid_Get &
       &,Cosmology_Age_Get,Expansion_Factor_Get,Hubble_Parameter_Get,Early_Time_Density_Scaling_Get,Omega_Matter_Total_Get &
       &,Omega_Dark_Energy_Get,Expansion_Rate_Get,Epoch_of_Matter_Dark_Energy_Equality_Get,Epoch_of_Matter_Domination_Get &
       &,Epoch_of_Matter_Curvature_Equality_Get,CMB_Temperature_Get,Comoving_Distance_Get,Time_From_Comoving_Distance_Get&
       &,Comoving_Distance_Conversion_Get,Cosmology_Dark_Energy_Equation_Of_State_Get,Cosmology_Dark_Energy_Exponent_Get)
    !% Initialize the module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                                         ), intent(in   )          :: cosmologyMethod
    procedure(Early_Time_Density_Scaling_Matter_Dark_Energy          ), intent(inout), pointer :: Early_Time_Density_Scaling_Get
    procedure(Expansion_Factor_Is_Valid_Matter_Dark_Energy           ), intent(inout), pointer :: Expansion_Factor_Is_Valid_Get
    procedure(Cosmic_Time_Is_Valid_Matter_Dark_Energy                ), intent(inout), pointer :: Cosmic_Time_Is_Valid_Get
    procedure(Cosmology_Age_Matter_Dark_Energy                       ), intent(inout), pointer :: Cosmology_Age_Get
    procedure(Expansion_Factor_Matter_Dark_Energy                    ), intent(inout), pointer :: Expansion_Factor_Get
    procedure(Hubble_Parameter_Matter_Dark_Energy                    ), intent(inout), pointer :: Hubble_Parameter_Get
    procedure(Omega_Matter_Total_Matter_Dark_Energy                  ), intent(inout), pointer :: Omega_Matter_Total_Get
    procedure(Omega_Dark_Energy_Matter_Dark_Energy                   ), intent(inout), pointer :: Omega_Dark_Energy_Get
    procedure(Expansion_Rate_Matter_Dark_Energy                      ), intent(inout), pointer :: Expansion_Rate_Get
    procedure(Epoch_of_Matter_Dark_Energy_Equality_Matter_Dark_Energy), intent(inout), pointer :: Epoch_of_Matter_Dark_Energy_Equality_Get
    procedure(Epoch_of_Matter_Domination_Matter_Dark_Energy          ), intent(inout), pointer :: Epoch_of_Matter_Domination_Get
    procedure(Epoch_of_Matter_Curvature_Equality_Matter_Dark_Energy  ), intent(inout), pointer :: Epoch_of_Matter_Curvature_Equality_Get
    procedure(CMB_Temperature_Matter_Dark_Energy                     ), intent(inout), pointer :: CMB_Temperature_Get
    procedure(Comoving_Distance_Matter_Dark_Energy                   ), intent(inout), pointer :: Comoving_Distance_Get
    procedure(Time_From_Comoving_Distance_Matter_Dark_Energy         ), intent(inout), pointer :: Time_From_Comoving_Distance_Get
    procedure(Comoving_Distance_Conversion_Matter_Dark_Energy        ), intent(inout), pointer :: Comoving_Distance_Conversion_Get
    procedure(Cosmology_Dark_Energy_Equation_Of_State_Dark_Energy    ), intent(inout), pointer :: Cosmology_Dark_Energy_Equation_Of_State_Get
    procedure(Cosmology_Dark_Energy_Exponent_Dark_Energy             ), intent(inout), pointer :: Cosmology_Dark_Energy_Exponent_Get

    ! Check if our method is selected.
    if (cosmologyMethod == 'matter-darkEnergy') then
       ! Set up procedure pointers.
       Expansion_Factor_Is_Valid_Get               => Expansion_Factor_Is_Valid_Matter_Dark_Energy
       Cosmic_Time_Is_Valid_Get                    => Cosmic_Time_Is_Valid_Matter_Dark_Energy
       Cosmology_Age_Get                           => Cosmology_Age_Matter_Dark_Energy
       Expansion_Factor_Get                        => Expansion_Factor_Matter_Dark_Energy
       Hubble_Parameter_Get                        => Hubble_Parameter_Matter_Dark_Energy
       Early_Time_Density_Scaling_Get              => Early_Time_Density_Scaling_Matter_Dark_Energy
       Omega_Matter_Total_Get                      => Omega_Matter_Total_Matter_Dark_Energy
       Omega_Dark_Energy_Get                       => Omega_Dark_Energy_Matter_Dark_Energy
       Expansion_Rate_Get                          => Expansion_Rate_Matter_Dark_Energy
       Epoch_of_Matter_Dark_Energy_Equality_Get    => Epoch_of_Matter_Dark_Energy_Equality_Matter_Dark_Energy
       Epoch_of_Matter_Curvature_Equality_Get      => Epoch_of_Matter_Curvature_Equality_Matter_Dark_Energy
       Epoch_of_Matter_Domination_Get              => Epoch_of_Matter_Domination_Matter_Dark_Energy
       CMB_Temperature_Get                         => CMB_Temperature_Matter_Dark_Energy
       Comoving_Distance_Get                       => Comoving_Distance_Matter_Dark_Energy
       Time_From_Comoving_Distance_Get             => Time_From_Comoving_Distance_Matter_Dark_Energy
       Comoving_Distance_Conversion_Get            => Comoving_Distance_Conversion_Matter_Dark_Energy
       Cosmology_Dark_Energy_Equation_Of_State_Get => Cosmology_Dark_Energy_Equation_Of_State_Dark_Energy
       Cosmology_Dark_Energy_Exponent_Get          => Cosmology_Dark_Energy_Exponent_Dark_Energy
       ! Read the dark energy equation of state.
       !@ <inputParameter>
       !@   <name>darkEnergyEquationOfStateW0</name>
       !@   <defaultValue>-1 (cosmological constant)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The equation of state parameter for dark energy, $w_0$, defined such that $P=\rho^w$ with $w(a)=w_0+w_1 a (1-a)$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>cosmology</group>
       !@ </inputParameter>
       call Get_Input_Parameter('darkEnergyEquationOfStateW0',darkEnergyEquationOfStateW0,defaultValue=-1.0d0)
       !@ <inputParameter>
       !@   <name>darkEnergyEquationOfStateW1</name>
       !@   <defaultValue>0 (constant equation of state)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The equation of state parameter for dark energy, $w_1$, defined such that $P=\rho^w$ with $w(a)=w_0+w_1 a (1-a)$.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>cosmology</group>
       !@ </inputParameter>
       call Get_Input_Parameter('darkEnergyEquationOfStateW1',darkEnergyEquationOfStateW1,defaultValue=0.0d0)
    end if
    return
  end subroutine Cosmology_Functions_Matter_Dark_Energy_Initialize

  logical function Expansion_Factor_Is_Valid_Matter_Dark_Energy(aExpansion)
    !% Checks that the expansion factor falls within allowed ranges.
    implicit none
    double precision, intent(in   ) :: aExpansion

    Expansion_Factor_Is_Valid_Matter_Dark_Energy=(aExpansion > 0.0d0 .and. (aExpansion<aExpansionMax .or. .not.collapsingUniverse))
    return
  end function Expansion_Factor_Is_Valid_Matter_Dark_Energy

  logical function Cosmic_Time_Is_Valid_Matter_Dark_Energy(time)
    !% Checks that the time falls within allowed ranges.
    implicit none
    double precision, intent(in   ) :: time

    Cosmic_Time_Is_Valid_Matter_Dark_Energy=(time > 0.0d0 .and. (time<tCosmologicalMax .or. .not.collapsingUniverse))
    return
  end function Cosmic_Time_Is_Valid_Matter_Dark_Energy

  double precision function Cosmology_Age_Matter_Dark_Energy(aExpansion,collapsingPhase)
    use Galacticus_Error
    use Numerical_Interpolation
    implicit none
    double precision, intent(in   )           :: aExpansion
    logical         , intent(in   ), optional :: collapsingPhase
    logical                                   :: collapsingPhaseActual

    ! Initialize the module if necessary.
    ! Validate the input.
    if (.not.Expansion_Factor_Is_Valid_Matter_Dark_Energy(aExpansion)) call Galacticus_Error_Report('Cosmology_Age_Matter_Dark_Energy'&
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
       ! We assume that the universe does not collapse.
       call Galacticus_Error_Report('Cosmology_Age_Matter_Dark_Energy','non-collapsing universe assumed')
    else
       ! In expanding phase ensure that sufficiently small and large expansion factors have been reached.
       do while (ageTableExpansionFactor(1)>aExpansion)
          ageTableTimeMinimum=ageTableTimeMinimum/ageTableIncrementFactor
          call Make_Expansion_Factor_Table
       end do
       do while (ageTableExpansionFactor(ageTableNumberPoints)<aExpansion)
          ageTableTimeMaximum=ageTableTimeMaximum*ageTableIncrementFactor
          call Make_Expansion_Factor_Table
       end do
    end if

    ! Interpolate to get cosmic time.
    Cosmology_Age_Matter_Dark_Energy=Interpolate(ageTableNumberPoints,ageTableExpansionFactor&
         &,ageTableTime,interpolationObject,interpolationAccelerator&
         &,aExpansion,reset=resetInterpolation)
    return
  end function Cosmology_Age_Matter_Dark_Energy

  double precision function Expansion_Factor_Matter_Dark_Energy(tCosmological)
    !% Returns the expansion factor at cosmological time {\tt tCosmological}.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: tCosmological
    double precision, save          :: expansionFactorPrevious, tCosmologicalPrevious
    !$omp threadprivate(tCosmologicalPrevious,expansionFactorPrevious)
    logical                         :: remakeTable
    double precision                :: tEffective

    ! Check if the time differs from the previous time.
    if (tCosmological /= tCosmologicalPrevious) then

       ! Quit on invalid input.
       if (tCosmological<0.0d0) call Galacticus_Error_Report('Expansion_Factor_Matter_Dark_Energy','cosmological time must be positive')

       ! Check if we need to recompute our table.
       if (ageTableInitialized) then
          remakeTable=(tCosmological<ageTableTime(1).or.tCosmological>ageTableTime(ageTableNumberPoints))
       else
          remakeTable=.true.
       end if
       if (remakeTable) call Make_Expansion_Factor_Table(tCosmological)

       ! Quit on invalid input.
       if (collapsingUniverse.and.tCosmological>tCosmologicalMax) call Galacticus_Error_Report('Expansion_Factor_Matter_Dark_Energy','cosmological time&
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
    Expansion_Factor_Matter_Dark_Energy=expansionFactorPrevious
    return
  end function Expansion_Factor_Matter_Dark_Energy

  subroutine Make_Expansion_Factor_Table(tCosmological)
    !% Builds a table of expansion factor vs. time.
    use Numerical_Interpolation
    use Numerical_Ranges
    use Memory_Management
    implicit none
    double precision, intent(in   ), optional     :: tCosmological
    double precision, allocatable  , dimension(:) :: ageTableExpansionFactorTemporary           , ageTableTimeTemporary
    double precision, parameter                   :: turnaroundTimeTolerance            =1.0d-12
    double precision, parameter                   :: dominateFactor                     =100.0d0
    integer                                       :: iTime                                      , prefixPointCount
    double precision                              :: Omega_Dominant                             , aDominant            , &
         &                                           time                                       , deltaTime            , &
         &                                           densityPower                               , tDominant
    logical                                       :: solutionFound                              , timeExceeded

    ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
    call Early_Time_Density_Scaling_Matter_Dark_Energy(dominateFactor,densityPower,aDominant,Omega_Dominant)

    ! Find the corresponding time.
    tDominant=-2.0d0/densityPower/H_0_invGyr()/sqrt(Omega_Dominant)/aDominant**(0.5d0*densityPower)

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
    ageTableNumberPoints=int(log10(ageTableTimeMaximum/ageTableTimeMinimum)*dble(ageTableNPointsPerDecade))+1
    ageTableTimeMaximum=ageTableTimeMinimum*10.0d0**(dble(ageTableNumberPoints)/dble(ageTableNPointsPerDecade))

    ! Deallocate arrays if currently allocated.
    if (allocated(ageTableTime)) then
       ! Determine number of points that are being added at the start of the array.
       prefixPointCount=int(log10(ageTableTime(1)/ageTableTimeMinimum)*dble(ageTableNPointsPerDecade)+0.5d0)
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
         &*sqrt(Omega_Dominant))**(-2.0d0/densityPower)
    ! Solve ODE to get corresponding expansion factors.
    do iTime=2,ageTableNumberPoints
       ! Compute the expansion factor if it is not already computed.
       if (ageTableExpansionFactor(iTime) < 0.0d0) then
          ageTableExpansionFactor(iTime)=Expansion_Factor_Change(                                  &
               &                                                 ageTableTime           (iTime-1), &
               &                                                 ageTableTime           (iTime  ), &
               &                                                 ageTableExpansionFactor(iTime-1)  &
               &                                                )
          ! Check for a universe which is no longer expanding (i.e. has reached its maximum expansion).
          if (ageTableExpansionFactor(iTime) == ageTableExpansionFactor(iTime-1)) then
             ! Record that we have a collapsing Universe.
             collapsingUniverse=.true.
             ! Record the maximum expansion factor.
             aExpansionMax          =ageTableExpansionFactor(iTime-1)
             ! Find the time of maximum expansion by bisection.
             tCosmologicalTurnaround=(ageTableTime(iTime-1)+ageTableTime(iTime-2))/2.0d0
             deltaTime              =(ageTableTime(iTime-1)-ageTableTime(iTime-2))/2.0d0
             solutionFound=.false.
             do while (.not.solutionFound)
                timeExceeded=(                                                           &
                &              Expansion_Factor_Change(                                  &
                &                                      ageTableTime           (iTime-2), &
                &                                      tCosmologicalTurnaround         , &
                &                                      ageTableExpansionFactor(iTime-2)  &
                &                                     )                                  &
                &             >=                                                         &
                &              ageTableExpansionFactor(iTime-1)                          &
                &            )
                solutionFound=timeExceeded .and. deltaTime < turnaroundTimeTolerance*tCosmologicalTurnaround
                if (.not.solutionFound) then
                   deltaTime=0.5d0*deltaTime
                   if (timeExceeded) then
                      tCosmologicalTurnaround=tCosmologicalTurnaround-deltaTime
                   else
                      tCosmologicalTurnaround=tCosmologicalTurnaround+deltaTime
                   end if
                end if
             end do
             tCosmologicalMax       =2.0d0*tCosmologicalTurnaround
             ! Limit the tables to the expanding part of the evolution.
             iTableTurnaround       =iTime-2
             call Move_Alloc(ageTableTime           ,ageTableTimeTemporary           )
             call Move_Alloc(ageTableExpansionFactor,ageTableExpansionFactorTemporary)
             ageTableNumberPoints=iTableTurnaround
             call Alloc_Array(ageTableTime,           [ageTableNumberPoints])
             call Alloc_Array(ageTableExpansionFactor,[ageTableNumberPoints])
             ageTableTime           =ageTableTimeTemporary           (1:ageTableNumberPoints)
             ageTableExpansionFactor=ageTableExpansionFactorTemporary(1:ageTableNumberPoints)
             exit
          end if
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

  double precision function Expansion_Factor_Change(timeStart,timeEnd,expansionFactorStart)
    !% Compute the expansion factor at time {\tt timeEnd} given an initial value {\tt expansionFactorStart} at time {\tt
    !% timeStart}.
    use ODE_Solver
    implicit none
    double precision       , intent(in   ) :: expansionFactorStart       , timeEnd                     , &
         &                                    timeStart
    double precision       , dimension(1)  :: y
    double precision       , parameter     :: odeToleranceAbsolute=1.0d-9, odeToleranceRelative=1.0d-12
    double precision                       :: time
    type            (c_ptr)                :: parameterPointer

    time=timeStart
    y(1)=expansionFactorStart
    call ODE_Solve(odeStepper,odeController,odeEvolver,odeSystem,time,timeEnd,1,y&
         &,ageTableODEsDarkEnergy,parameterPointer,odeToleranceAbsolute,odeToleranceRelative,reset=odeReset)
    Expansion_Factor_Change=y(1)
    return
  end function Expansion_Factor_Change

  function ageTableODEsDarkEnergy(t,a,dadt,parameterPointer) bind(c)
    !% System of differential equations to solve for expansion factor vs. age.
    integer(kind=c_int   )                              :: ageTableODEsDarkEnergy
    real   (kind=c_double)              , value         :: t
    real   (kind=c_double), dimension(1), intent(in   ) :: a
    real   (kind=c_double), dimension(1)                :: dadt
    type   (c_ptr        )              , value         :: parameterPointer

    dadt(1)=a(1)*Expansion_Rate_Matter_Dark_Energy(a(1))
    ageTableODEsDarkEnergy=FGSL_Success
  end function ageTableODEsDarkEnergy

  double precision function Expansion_Rate_Matter_Dark_Energy(aExpansion)
    !% Returns the cosmological expansion rate, $\dot{a}/a$ at expansion factor {\tt aExpansion}.
    implicit none
    double precision, intent(in   ) :: aExpansion

    ! Required value is simply the Hubble parameter but expressed in units of inverse Gyr.
    Expansion_Rate_Matter_Dark_Energy=Hubble_Parameter_Matter_Dark_Energy(aExpansion=aExpansion)*H_0_invGyr()/H_0()
    return
  end function Expansion_Rate_Matter_Dark_Energy

  double precision function Hubble_Parameter_Matter_Dark_Energy(tCosmological,aExpansion,collapsingPhase)
    !% Returns the Hubble parameter at the request cosmological time, {\tt tCosmological}, or expansion factor, {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: aExpansion      , tCosmological
    logical         , intent(in   ), optional :: collapsingPhase
    double precision                          :: aExpansionActual, sqrtArgument

    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Hubble_Parameter_Matter_Dark_Energy','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Dark_Energy(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('Hubble_Parameter_Matter_Dark_Energy','either a time or expansion factor must be specified')
       end if
    end if
    ! Compute the Hubble parameter at the specified expansion factor.
    sqrtArgument=max(Omega_Matter()/aExpansionActual**3+Omega_DE()*aExpansionActual**Cosmology_Dark_Energy_Exponent_Dark_Energy(expansionFactor=aExpansionActual)+Omega_K()/aExpansionActual**2,0.0d0)
    Hubble_Parameter_Matter_Dark_Energy=H_0()*sqrt(sqrtArgument)
    ! Make the Hubble parameter negative if we are in the collapsing phase of the Universe.
    if (collapsingUniverse) then
       if (present(tCosmological)) then
          if (tCosmological>tCosmologicalTurnaround) Hubble_Parameter_Matter_Dark_Energy=-Hubble_Parameter_Matter_Dark_Energy
       else
          if (present(collapsingPhase)) then
             if (collapsingPhase) Hubble_Parameter_Matter_Dark_Energy=-Hubble_Parameter_Matter_Dark_Energy
          end if
       end if
    end if
    return
  end function Hubble_Parameter_Matter_Dark_Energy

  subroutine Early_Time_Density_Scaling_Matter_Dark_Energy(dominateFactor,densityPower,aDominant,Omega_Dominant)
    implicit none
    double precision, intent(in   )           :: dominateFactor
    double precision, intent(  out)           :: aDominant     , densityPower
    double precision, intent(  out), optional :: Omega_Dominant

    ! For matter and cosmological constant, matter always dominates at early times.
    densityPower=-3.0d0 ! Power-law scaling of matter density with expansion factor.

    ! Choose present day as default - will be used if no other densities present (i.e. Einsetin-de Sitter).
    aDominant=Epoch_of_Matter_Domination_Matter_Dark_Energy(dominateFactor)

    ! Return the density parameter in the dominant species if required.
    if (present(Omega_Dominant)) Omega_Dominant=Omega_Matter()
    return
  end subroutine Early_Time_Density_Scaling_Matter_Dark_Energy

  double precision function Epoch_of_Matter_Domination_Matter_Dark_Energy(dominateFactor)
    use Cosmology_Functions_Parameters
    use Root_Finder
    implicit none
    double precision                , intent(in   ) :: dominateFactor
    type            (rootFinder    ), save          :: finder
    !$omp threadprivate(finder)
    double precision                                :: aDominantCurvature , aDominantDarkEnergy      , &
         &                                             aMatterEquality    , darkEnergyExponentCurrent, &
         &                                             rangeExpandDownward, rangeExpandUpward

    ! Choose present day as default - will be used if no other densities present (i.e. Einsetin-de Sitter).
    Epoch_of_Matter_Domination_Matter_Dark_Energy=1.0d0

    ! Case where dark energy is present.
    if (Omega_DE() /= 0.0d0) then
       if (.not.finder%isInitialized()) then
          darkEnergyExponentCurrent=Cosmology_Dark_Energy_Exponent_Dark_Energy(expansionFactor=1.0d0)
          if (darkEnergyExponentCurrent > -3.0d0) then
             ! Dark energy density is decaying less rapidly than matter.
             if (Omega_Matter() < dominateFactor*Omega_DE()) then
                ! Matter density is less than dominated dark energy density. Search backward for epoch of domination.
                rangeExpandUpward  =1.0d0
                rangeExpandDownward=0.5d0
             else
                ! Matter density is greater than dominated dark energy density. Search forward for epoch of domination.
                rangeExpandUpward  =2.0d0
                rangeExpandDownward=1.0d0
             end if
          else
             ! Dark energy density is decaying more rapidly than matter.
             if (Omega_Matter() < dominateFactor*Omega_DE()) then
                ! Matter density is less than dominated dark energy density. Search forward for epoch of domination.
                rangeExpandUpward  =2.0d0
                rangeExpandDownward=1.0d0
             else
                ! Matter density is greater than dominated dark energy density. Search backward for epoch of domination.
                rangeExpandUpward  =1.0d0
                rangeExpandDownward=0.50d0
             end if
          end if
          call finder%rootFunction(                                               &
               &                   Matter_Dark_Energy_Domination                  &
               &                  )
          call finder%tolerance   (                                               &
               &                   toleranceAbsolute  =0.0d0                    , &
               &                   toleranceRelative  =1.0d-6                     &
               &                  )
          call finder%rangeExpand (                                               &
               &                   rangeExpandUpward  =rangeExpandUpward        , &
               &                   rangeExpandDownward=rangeExpandDownward      , &
               &                   rangeExpandType    =rangeExpandMultiplicative  &
               &                  )
       end if
       dominateFactorCurrent=dominateFactor
       aDominantDarkEnergy  =finder%find(rootGuess=1.0d0)

       ! Choose earliest expansion factor.
       Epoch_of_Matter_Domination_Matter_Dark_Energy=min(Epoch_of_Matter_Domination_Matter_Dark_Energy,aDominantDarkEnergy)
    end if

    if (Omega_K() /= 0.0d0) then
       ! Find the expansion factor of matter-curvature equality.
       aMatterEquality=Epoch_of_Matter_Curvature_Equality_Matter_Dark_Energy(requestTypeExpansionFactor)

       ! Find the earlier expansion factor at which matter dominates by the specified amount (ratio of matter
       ! to curvature density scales as the expansion factor).
       aDominantCurvature=aMatterEquality/dominateFactor

       ! Choose earliest expansion factor.
       Epoch_of_Matter_Domination_Matter_Dark_Energy=min(Epoch_of_Matter_Domination_Matter_Dark_Energy,aDominantCurvature)
    end if
    return
  end function Epoch_of_Matter_Domination_Matter_Dark_Energy

  double precision function Matter_Dark_Energy_Domination(expansionFactor)
    implicit none
    double precision, intent(in   ) :: expansionFactor

    Matter_Dark_Energy_Domination=Omega_Matter()/expansionFactor**3-dominateFactorCurrent*Omega_DE()*expansionFactor**Cosmology_Dark_Energy_Exponent_Dark_Energy(expansionFactor=expansionFactor)
    return
  end function Matter_Dark_Energy_Domination

  double precision function Omega_Matter_Total_Matter_Dark_Energy(tCosmological,aExpansion,collapsingPhase)
    !% Return the matter density parameter at expansion factor {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: aExpansion      , tCosmological
    logical         , intent(in   ), optional :: collapsingPhase
    double precision                          :: aExpansionActual

    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Omega_Matter_Total_Matter_Dark_Energy','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Dark_Energy(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('Omega_Matter_Total_Matter_Dark_Energy','either a time or expansion factor must be specified')
       end if
    end if
    Omega_Matter_Total_Matter_Dark_Energy=Omega_Matter()*((H_0()/Hubble_Parameter_Matter_Dark_Energy(aExpansion=aExpansionActual))**2)/(aExpansionActual**3)
    return
  end function Omega_Matter_Total_Matter_Dark_Energy

  double precision function Omega_Dark_Energy_Matter_Dark_Energy(tCosmological,aExpansion,collapsingPhase)
    !% Return the dark energy density parameter at expansion factor {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: aExpansion      , tCosmological
    logical         , intent(in   ), optional :: collapsingPhase
    double precision                          :: aExpansionActual

    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Omega_Dark_Energy_Matter_Dark_Energy','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Dark_Energy(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('Omega_Dark_Energy_Matter_Dark_Energy','either a time or expansion factor must be specified')
       end if
    end if
    Omega_Dark_Energy_Matter_Dark_Energy=Omega_DE()*aExpansionActual**Cosmology_Dark_Energy_Exponent_Dark_Energy(expansionFactor=aExpansionActual)*((H_0()/Hubble_Parameter_Matter_Dark_Energy(aExpansion=aExpansionActual))**2)
    return
  end function Omega_Dark_Energy_Matter_Dark_Energy

  double precision function CMB_Temperature_Matter_Dark_Energy(tCosmological,aExpansion,collapsingPhase)
    !% Return the temperature of the CMB at expansion factor {\tt aExpansion}.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: aExpansion      , tCosmological
    logical         , intent(in   ), optional :: collapsingPhase
    double precision                          :: aExpansionActual

    ! Determine the actual expansion factor to use.
    if (present(tCosmological)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('CMB_Temperature_Matter_Dark_Energy','only one of time or expansion factor can be specified')
       else
          aExpansionActual=Expansion_Factor_Matter_Dark_Energy(tCosmological)
       end if
    else
       if (present(aExpansion)) then
          aExpansionActual=aExpansion
       else
          call Galacticus_Error_Report('CMB_Temperature_Matter_Dark_Energy','either a time or expansion factor must be specified')
       end if
    end if
    CMB_Temperature_Matter_Dark_Energy=T_CMB()/aExpansionActual
    return
  end function CMB_Temperature_Matter_Dark_Energy

  double precision function Epoch_of_Matter_Dark_Energy_Equality_Matter_Dark_Energy(requestType)
    !% Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    use Cosmology_Functions_Parameters
    use Root_Finder
    implicit none
    integer                     , intent(in   ), optional :: requestType
    type            (rootFinder), save                    :: finder
    !$omp threadprivate(finder)
    integer                                               :: requestTypeActual
    double precision                                      :: darkEnergyExponentCurrent, rangeExpandDownward, &
         &                                                   rangeExpandUpward

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if
    if (.not.finder%isInitialized()) then
       darkEnergyExponentCurrent=Cosmology_Dark_Energy_Exponent_Dark_Energy(expansionFactor=1.0d0)
       if (darkEnergyExponentCurrent > -3.0d0) then
          ! Dark energy density is decaying less rapidly than matter.
          if (Omega_Matter() < Omega_DE()) then
             ! Matter density is less than dark energy density. Search backward for epoch of domination.
             rangeExpandUpward  =1.0d0
             rangeExpandDownward=0.5d0
          else
             ! Matter density is greater than dark energy density. Search forward for epoch of domination.
             rangeExpandUpward  =2.0d0
             rangeExpandDownward=1.0d0
          end if
       else
          ! Dark energy density is decaying more rapidly than matter.
          if (Omega_Matter() < Omega_DE()) then
             ! Matter density is less than dark energy density. Search forward for epoch of domination.
             rangeExpandUpward  =2.0d0
             rangeExpandDownward=1.0d0
          else
             ! Matter density is greater than dark energy density. Search backward for epoch of domination.
             rangeExpandUpward  =1.0d0
             rangeExpandDownward=0.50d0
          end if
       end if
       call finder%rootFunction(                                               &
            &                   Matter_Dark_Energy_Domination                  &
            &                  )
       call finder%tolerance   (                                               &
            &                   toleranceAbsolute  =0.0d0                    , &
            &                   toleranceRelative  =1.0d-6                     &
            &                  )
       call finder%rangeExpand (                                               &
            &                   rangeExpandUpward  =rangeExpandUpward        , &
            &                   rangeExpandDownward=rangeExpandDownward      , &
            &                   rangeExpandType    =rangeExpandMultiplicative  &
            &                  )
    end if
    dominateFactorCurrent=1.0d0
    Epoch_of_Matter_Dark_Energy_Equality_Matter_Dark_Energy=finder%find(rootGuess=1.0d0)
    if (present(requestType)) then
       if (requestType == requestTypeTime) Epoch_of_Matter_Dark_Energy_Equality_Matter_Dark_Energy &
            &=Cosmology_Age_Matter_Dark_Energy(Epoch_of_Matter_Dark_Energy_Equality_Matter_Dark_Energy)
    end if
    return
  end function Epoch_of_Matter_Dark_Energy_Equality_Matter_Dark_Energy

  double precision function Epoch_of_Matter_Curvature_Equality_Matter_Dark_Energy(requestType)
    !% Return the epoch of matter-curvature magnitude equality (either expansion factor or cosmic time).
    use Cosmology_Functions_Parameters
    implicit none
    integer, intent(in   ), optional :: requestType
    integer                          :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if

    Epoch_of_Matter_Curvature_Equality_Matter_Dark_Energy=Omega_Matter()/abs(Omega_K())
    if (requestType == requestTypeTime) Epoch_of_Matter_Curvature_Equality_Matter_Dark_Energy&
         &=Cosmology_Age_Matter_Dark_Energy(Epoch_of_Matter_Curvature_Equality_Matter_Dark_Energy)
    return
  end function Epoch_of_Matter_Curvature_Equality_Matter_Dark_Energy

  double precision function Time_From_Comoving_Distance_Matter_Dark_Energy(comovingDistance)
    !% Returns the cosmological time corresponding to given {\tt comovingDistance}.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: comovingDistance

    call Galacticus_Error_Report('Time_From_Comoving_Distance_Matter_Dark_Energy','functionality not implemented')
    return
  end function Time_From_Comoving_Distance_Matter_Dark_Energy

  double precision function Comoving_Distance_Matter_Dark_Energy(tCosmological)
    !% Returns the comoving distance to cosmological time {\tt tCosmological}.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ) :: tCosmological

    call Galacticus_Error_Report('Comoving_Distance_Matter_Dark_Energy','functionality not implemented')
   return
  end function Comoving_Distance_Matter_Dark_Energy

  double precision function Comoving_Distance_Conversion_Matter_Dark_Energy(output,distanceModulus,redshift)
    !% Convert bewteen different measures of distance.
    use Galacticus_Error
    implicit none
    integer         , intent(in   )           :: output
    double precision, intent(in   ), optional :: distanceModulus, redshift

    call Galacticus_Error_Report('Comoving_Distance_Conversion_Matter_Dark_Energy','functionality not implemented')
    return
  end function Comoving_Distance_Conversion_Matter_Dark_Energy

  double precision function Cosmology_Dark_Energy_Equation_Of_State_Dark_Energy(time,expansionFactor)
    !% Return the dark energy equation of state.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: expansionFactor      , time
    double precision                          :: expansionFactorActual

    if (present(expansionFactor)) then
       expansionFactorActual=expansionFactor
    else if (present(time)) then
       expansionFactorActual=Expansion_Factor_Matter_Dark_Energy(time)
    else
       if (darkEnergyEquationOfStateW1 /= 0.0d0) call Galacticus_Error_Report('Cosmology_Dark_Energy_Equation_Of_State_Dark_Energy','equation of state is time dependent, but no time given')
       expansionFactorActual=1.0d0
    end if
    Cosmology_Dark_Energy_Equation_Of_State_Dark_Energy=darkEnergyEquationOfStateW0+darkEnergyEquationOfStateW1&
         &*expansionFactorActual*(1.0d0-expansionFactorActual)
    return
  end function Cosmology_Dark_Energy_Equation_Of_State_Dark_Energy

  double precision function Cosmology_Dark_Energy_Exponent_Dark_Energy(time,expansionFactor)
    !% Return the dark energy equation of state.
    use Galacticus_Error
    implicit none
    double precision, intent(in   ), optional :: expansionFactor      , time
    double precision                          :: expansionFactorActual

    if      (present(expansionFactor)) then
       expansionFactorActual=expansionFactor
    else if (present(time           )) then
       expansionFactorActual=Expansion_Factor_Matter_Dark_Energy(time)
    else
       if (darkEnergyEquationOfStateW1 /= 0.0d0) call Galacticus_Error_Report('Cosmology_Dark_Energy_Exponent_Dark_Energy','equation of state is time dependent, but no time given')
       expansionFactorActual=1.0d0
    end if
    if (expansionFactorActual == 1.0d0) then
       Cosmology_Dark_Energy_Exponent_Dark_Energy=-3.0d0*(1.0d0+darkEnergyEquationOfStateW0)
    else
       Cosmology_Dark_Energy_Exponent_Dark_Energy=-3.0d0*(1.0d0+darkEnergyEquationOfStateW0)+3.0d0*darkEnergyEquationOfStateW1*(1.0d0-expansionFactorActual)**2/2.0d0/log(expansionFactorActual)
    end if
    return
  end function Cosmology_Dark_Energy_Exponent_Dark_Energy

  !# <galacticusStateStoreTask>
  !#  <unitName>Cosmology_Matter_Dark_Energy_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Cosmology_Matter_Dark_Energy_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    ! Store the full tables, as they are hysteretic and cannot be reconstructed precisely without knowing the path by which they
    ! were originally constructed.
    write (stateFile) ageTableNumberPoints,ageTableTimeMinimum,ageTableTimeMaximum
    write (stateFile) ageTableTime,ageTableExpansionFactor
    return
  end subroutine Cosmology_Matter_Dark_Energy_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Cosmology_Matter_Dark_Energy_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Cosmology_Matter_Dark_Energy_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use Memory_Management
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    ! Read the tabulations.
    read (stateFile) ageTableNumberPoints,ageTableTimeMinimum,ageTableTimeMaximum
    if (allocated(ageTableTime           )) call Dealloc_Array(ageTableTime           )
    if (allocated(ageTableExpansionFactor)) call Dealloc_Array(ageTableExpansionFactor)
    call Alloc_Array(ageTableTime           ,[ageTableNumberPoints])
    call Alloc_Array(ageTableExpansionFactor,[ageTableNumberPoints])
    read (stateFile) ageTableTime,ageTableExpansionFactor

    ! Ensure that interpolation objects will get reset.
    resetInterpolation       =.true.
    resetInterpolationInverse=.true.
    return
  end subroutine Cosmology_Matter_Dark_Energy_State_Retrieve

end module Cosmology_Functions_Matter_Dark_Energy
