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

!% Contains a module that implements calculations of the intergalactic medium thermal and ionization state.

module Intergalactic_Medium_State
  !% Implements calculations of the intergalactic medium thermal and ionization state.
  use Tables
  implicit none
  private
  public :: Intergalactic_Medium_Electron_Fraction, Intergalactic_Medium_Temperature,&
       & Intergalactic_Medium_Electron_Scattering_Optical_Depth, Intergalactic_Medium_Electron_Scattering_Time,&
       & Intergalactic_Medium_State_State_Retrieve, Intergalactic_Medium_State_State_Store

  ! Flag to indicate if this module has been initialized.
  logical                                                     :: igmStateInitialized                       =.false.

  ! Pointer to the function that actually does the calculation.
  procedure(Intergalactic_Medium_State_Get_Template), pointer :: Intergalactic_Medium_Electron_Fraction_Get=>null()
  procedure(Intergalactic_Medium_State_Get_Template), pointer :: Intergalactic_Medium_Temperature_Get      =>null()
  abstract interface
     double precision function Intergalactic_Medium_State_Get_Template(time)
       double precision, intent(in   ) :: time
     end function Intergalactic_Medium_State_Get_Template
  end interface

  ! Electron scattering optical depth tables.
  integer                                   , parameter   :: electronScatteringTablePointsPerDecade=100
  logical                                                 :: electronScatteringTableInitialized    =.false.
  integer                                                 :: electronScatteringTableNumberPoints
  double precision                                        :: electronScatteringTableTimeMaximum            , electronScatteringTableTimeMinimum
  type            (table1DLogarithmicLinear)              :: electronScattering                            , electronScatteringFullyIonized
  class           (table1D                 ), allocatable :: electronScatteringInverse                     , electronScatteringFullyIonizedInverse

  ! Option controlling whether electron scattering optical depth calculations should assume a fully ionized universe.
  logical                                                 :: fullyIonized

contains

  subroutine Intergalactic_Medium_State_Initialize
    !% Initialize the intergalactic medium state module.
    use ISO_Varying_String
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="intergalaticMediumStateMethod" type="moduleUse">
    include 'intergalactic_medium.state.modules.inc'
    !# </include>
    implicit none
    type(varying_string) :: intergalaticMediumStateMethod

    ! Initialize if necessary.
    if (.not.igmStateInitialized) then
       !$omp critical(Intergalactic_Medium_State_Initialization)
       if (.not.igmStateInitialized) then
          ! Get the cooling function method parameter.
          !@ <inputParameter>
          !@   <name>intergalaticMediumStateMethod</name>
          !@   <defaultValue>RecFast</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the state of the intergalactic medium.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('intergalaticMediumStateMethod',intergalaticMediumStateMethod,defaultValue='RecFast')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="intergalaticMediumStateMethod" type="functionCall" functionType="void">
          !#  <functionArgs>intergalaticMediumStateMethod,Intergalactic_Medium_Electron_Fraction_Get,Intergalactic_Medium_Temperature_Get</functionArgs>
          include 'intergalactic_medium.state.inc'
          !# </include>
          if (.not.(associated(Intergalactic_Medium_Electron_Fraction_Get).and.associated(Intergalactic_Medium_Temperature_Get))) call&
               & Galacticus_Error_Report('Intergalactic_Medium_State_Initialize','method ' //char(intergalaticMediumStateMethod)//' is unrecognized')

          igmStateInitialized=.true.
       end if
       !$omp end critical(Intergalactic_Medium_State_Initialization)
    end if
    return
  end subroutine Intergalactic_Medium_State_Initialize

  double precision function Intergalactic_Medium_Electron_Fraction(time)
    !% Return the electron fraction in the intergalactic medium at the specified {\tt time}.
    implicit none
    double precision, intent(in   ) :: time

    ! Initialize the module.
    call Intergalactic_Medium_State_Initialize

    Intergalactic_Medium_Electron_Fraction=Intergalactic_Medium_Electron_Fraction_Get(time)
    return
  end function Intergalactic_Medium_Electron_Fraction

  double precision function Intergalactic_Medium_Temperature(time)
    !% Return the temperature of the intergalactic medium at the specified {\tt time}.
    implicit none
    double precision, intent(in   ) :: time

    ! Initialize the module.
    call Intergalactic_Medium_State_Initialize

    Intergalactic_Medium_Temperature=Intergalactic_Medium_Temperature_Get(time)
    return
  end function Intergalactic_Medium_Temperature

  double precision function Intergalactic_Medium_Electron_Scattering_Optical_Depth(time,assumeFullyIonized)
    !% Return the electron scattering optical depth from the present day back to the given {\tt time} in the intergalactic medium.
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    double precision, intent(in   )           :: time
    logical         , intent(in   ), optional :: assumeFullyIonized
    logical                                   :: assumeFullyIonizedActual

    ! Ensure that the table is initialized.
    call IGM_State_Electron_Scattering_Tabulate(time)

    ! Check for invalid input.
    if (time > electronScatteringTableTimeMaximum) call Galacticus_Error_Report('Intergalactic_Medium_Electron_Scattering_Optical_Depth','time exceeds present age of the universe')

    !$omp critical     (IGM_State_Electron_Scattering_Interpolation)
    assumeFullyIonizedActual=.false.
    if (present(assumeFullyIonized)) assumeFullyIonizedActual=assumeFullyIonized
    if (assumeFullyIonizedActual) then
       Intergalactic_Medium_Electron_Scattering_Optical_Depth=-electronScatteringFullyIonized%interpolate(time)
    else
       Intergalactic_Medium_Electron_Scattering_Optical_Depth=-electronScattering            %interpolate(time)
    end if
    !$omp end critical (IGM_State_Electron_Scattering_Interpolation)
    return
  end function Intergalactic_Medium_Electron_Scattering_Optical_Depth

  double precision function Intergalactic_Medium_Electron_Scattering_Time(opticalDepth,assumeFullyIonized)
    !% Return the cosmological time at which the given electron scattering {\tt opticalDepth} is reached (integrating from the
    !% present day) in the intergalactic medium.
    use Numerical_Interpolation
    use Galacticus_Error
    use Cosmology_Functions
    implicit none
    double precision, intent(in   )           :: opticalDepth
    logical         , intent(in   ), optional :: assumeFullyIonized
    logical                                   :: assumeFullyIonizedActual
    double precision                          :: time

    ! Check for invalid input.
    if (opticalDepth < 0.0d0) call Galacticus_Error_Report('Intergalactic_Medium_Electron_Scattering_Time','optical depth must be non-negative')

    ! Determine which optical depth to use.
    assumeFullyIonizedActual=.true.!false.
    if (present(assumeFullyIonized)) assumeFullyIonizedActual=assumeFullyIonized

    ! Ensure that the table is initialized.
    time=Cosmology_Age(1.0d0)
    call IGM_State_Electron_Scattering_Tabulate(time)
    do while (                                                                                           &
         &     (.not.assumeFullyIonizedActual .and. electronScattering            %y(1) > -opticalDepth) &
         &    .or.                                                                                       &
         &     (     assumeFullyIonizedActual .and. electronScatteringFullyIonized%y(1) > -opticalDepth) &
         &   )
       time=time/2.0d0
       call IGM_State_Electron_Scattering_Tabulate(time)
    end do

    !$omp critical     (IGM_State_Electron_Scattering_Interpolation)
    if (assumeFullyIonizedActual) then
       Intergalactic_Medium_Electron_Scattering_Time=electronScatteringFullyIonizedInverse%interpolate(-opticalDepth)
    else
       Intergalactic_Medium_Electron_Scattering_Time=electronScatteringInverse            %interpolate(-opticalDepth)
    end if
    !$omp end critical (IGM_State_Electron_Scattering_Interpolation)
    return
  end function Intergalactic_Medium_Electron_Scattering_Time

  subroutine IGM_State_Electron_Scattering_Tabulate(time)
    !% Construct a table of electron scattering optical depth as a function of cosmological time.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Numerical_Interpolation
    use Cosmology_Functions
    use Memory_Management
    use Numerical_Ranges
    use FGSL
    implicit none
    double precision                            , intent(in   ) :: time
    type            (c_ptr                     )                :: parameterPointer
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    integer                                                     :: iTime

    !$omp critical (IGM_State_Electron_Scattering_Interpolation)
    if (.not.electronScatteringTableInitialized.or.time < electronScatteringTableTimeMinimum) then
       ! Find minimum and maximum times to tabulate.
       electronScatteringTableTimeMaximum=    Cosmology_Age(1.0d0)
       electronScatteringTableTimeMinimum=min(Cosmology_Age(1.0d0),time)/2.0d0
       ! Decide how many points to tabulate and allocate table arrays.
       electronScatteringTableNumberPoints=int(log10(electronScatteringTableTimeMaximum/electronScatteringTableTimeMinimum)&
            &*dble(electronScatteringTablePointsPerDecade))+1
       ! Create the tables.
       call electronScattering            %destroy()
       call electronScatteringFullyIonized%destroy()
       call electronScattering            %create (                                     &
            &                                      electronScatteringTableTimeMinimum , &
            &                                      electronScatteringTableTimeMaximum , &
            &                                      electronScatteringTableNumberPoints  &
            &                                     )
       call electronScatteringFullyIonized%create (                                     &
            &                                      electronScatteringTableTimeMinimum , &
            &                                      electronScatteringTableTimeMaximum , &
            &                                      electronScatteringTableNumberPoints  &
            &                                     )
       ! Loop over times and populate tables.
       do iTime=1,electronScatteringTableNumberPoints-1
          fullyIonized=.false.
          call electronScattering%populate(                                                 &
               &                 -Integrate(                                                &
               &                            electronScattering%x(iTime)                   , &
               &                            electronScatteringTableTimeMaximum            , &
               &                            IGM_State_Electron_Scattering_Integrand       , &
               &                            parameterPointer                              , &
               &                            integrandFunction                             , &
               &                            integrationWorkspace                          , &
               &                            toleranceAbsolute                      =0.0d+0, &
               &                            toleranceRelative                      =1.0d-3  &
               &                           )                                              , &
               &                  iTime                                                     &
               &                 )
          call Integrate_Done(integrandFunction,integrationWorkspace)
          fullyIonized=.true.
          call electronScatteringFullyIonized%populate(                                                           &
               &                                       -Integrate(                                                &
               &                                                  electronScatteringFullyIonized%x(iTime)       , &
               &                                                  electronScatteringTableTimeMaximum            , &
               &                                                  IGM_State_Electron_Scattering_Integrand       , &
               &                                                  parameterPointer                              , &
               &                                                  integrandFunction                             , &
               &                                                  integrationWorkspace                          , &
               &                                                  toleranceAbsolute                      =0.0d+0, &
               &                                                  toleranceRelative                      =1.0d-3  &
               &                                                 )                                              , &
               &                  iTime                                                                           &
               &                 )
          call Integrate_Done(integrandFunction,integrationWorkspace)
       end do
       call electronScattering            %populate(0.0d0,electronScatteringTableNumberPoints)
       call electronScatteringFullyIonized%populate(0.0d0,electronScatteringTableNumberPoints)
       call electronScattering            %reverse (electronScatteringInverse                )
       call electronScatteringFullyIonized%reverse (electronScatteringFullyIonizedInverse    )
       ! Specify that tabulation has been made.
       electronScatteringTableInitialized=.true.
    end if
    !$omp end critical (IGM_State_Electron_Scattering_Interpolation)
    return
  end subroutine IGM_State_Electron_Scattering_Tabulate

  function IGM_State_Electron_Scattering_Integrand(time,parameterPointer) bind(c)
    !% Integrand for electron scattering optical depth calculations.
    use, intrinsic :: ISO_C_Binding
    use Cosmology_Functions
    use Cosmological_Parameters
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    real            (kind=c_double)        :: IGM_State_Electron_Scattering_Integrand
    real            (kind=c_double), value :: time
    type            (c_ptr        ), value :: parameterPointer
    double precision                       :: electronFraction                       , expansionFactor

    expansionFactor=Expansion_Factor(time)
    if (fullyIonized) then
       electronFraction=hydrogenByMassPrimordial/atomicMassHydrogen+2.0d0*heliumByMassPrimordial/atomicMassHelium
    else
       electronFraction=hydrogenByMassPrimordial*Intergalactic_Medium_Electron_Fraction(time)/atomicMassHydrogen
    end if
    IGM_State_Electron_Scattering_Integrand=speedLight*gigaYear*thomsonCrossSection*Omega_b()*Critical_Density()*massSolar*electronFraction/atomicMassUnit/megaParsec**3/expansionFactor**3
    return
  end function IGM_State_Electron_Scattering_Integrand

  !# <galacticusStateStoreTask>
  !#  <unitName>Intergalactic_Medium_State_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Intergalactic_Medium_State_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    write (stateFile)
    return
  end subroutine Intergalactic_Medium_State_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Intergalactic_Medium_State_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Intergalactic_Medium_State_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile
    type   (fgsl_file), intent(in   ) :: fgslStateFile

    ! Read the table state.
    read (stateFile)
    ! Force retabulation.
    electronScatteringTableInitialized=.false.
    return
  end subroutine Intergalactic_Medium_State_State_Retrieve

end module Intergalactic_Medium_State
