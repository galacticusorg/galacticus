!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the virial overdensity for halos.

module Virial_Density_Contrast
  !% Implements the virial overdensity for halos.
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Halo_Virial_Density_Contrast, Halo_Virial_Density_Contrast_Rate_of_Change, Virial_Density_Contrast_State_Retrieve

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: deltaVirialInitialized=.false.

  ! Variables to hold the tabulated virial overdensity data.
  integer                                        :: deltaVirialTableNumberPoints
  double precision,    allocatable, dimension(:) :: deltaVirialTableTime,deltaVirialTableDeltaVirial
  type(fgsl_interp)                              :: interpolationObject
  type(fgsl_interp_accel)                        :: interpolationAccelerator
  logical                                        :: resetInterpolation, tablesInitialized=.false.
  !$omp threadprivate(deltaVirialTableNumberPoints,deltaVirialTableTime,deltaVirialTableDeltaVirial)
  !$omp threadprivate(interpolationObject,interpolationAccelerator,resetInterpolation,tablesInitialized)

  ! Name of virial overdensity method used.
  type(varying_string)                           :: virialDensityContrastMethod

  ! Pointer to the subroutine that tabulates the virial overdensity and template interface for that subroutine.
  procedure(Virial_Density_Contrast_Tabulate_Template), pointer :: Virial_Density_Contrast_Tabulate => null()
  abstract interface
     subroutine Virial_Density_Contrast_Tabulate_Template(time,deltaVirialTableNumberPoints,deltaVirialTime,deltaVirialDeltaVirial)
       double precision,                            intent(in)    :: time
       double precision, allocatable, dimension(:), intent(inout) :: deltaVirialTime,deltaVirialDeltaVirial
       integer,                                     intent(out)   :: deltaVirialTableNumberPoints
     end subroutine Virial_Density_Contrast_Tabulate_Template
  end interface
  
contains

  subroutine Virial_Density_Contrast_Initialize(time)
    !% Initializes the virial overdensity module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="virialDensityContrastMethod" type="moduleUse">
    include 'structure_formation.CDM.virial_overdensity.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: time

    ! Initialize the module if necessary.
    if (.not.deltaVirialInitialized) then
       !$omp critical (Virial_Density_Contrast_Initialize)
       if (.not.deltaVirialInitialized) then
          ! Get the virial overdensity method parameter.
          !@ <inputParameter>
          !@   <name>virialDensityContrastMethod</name>
          !@   <defaultValue>sphericalTopHat</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Selects the method to be used for computing halo virial density contrasts.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('virialDensityContrastMethod',virialDensityContrastMethod,defaultValue='sphericalTopHat')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="virialDensityContrastMethod" type="code" action="subroutine">
          !#  <subroutineArgs>virialDensityContrastMethod,Virial_Density_Contrast_Tabulate</subroutineArgs>
          include 'structure_formation.CDM.virial_overdensity.inc'
          !# </include>
          if (.not.associated(Virial_Density_Contrast_Tabulate)) call Galacticus_Error_Report('Virial_Density_Contrast_Initialize','method ' &
               &//char(virialDensityContrastMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          deltaVirialInitialized=.true.
       end if
       !$omp end critical (Virial_Density_Contrast_Initialize)
    end if
    return
  end subroutine Virial_Density_Contrast_Initialize
  
  subroutine Virial_Density_Contrast_Retabulate(time)
    !% Recompute the look-up tables for virial density contrast.
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: time
    logical                      :: remakeTable

    ! Ensure that the module is initialized.
    call Virial_Density_Contrast_Initialize(time)
    
    ! Check if we need to recompute our table.
    if (tablesInitialized) then
       remakeTable=(time<deltaVirialTableTime(1).or.time>deltaVirialTableTime(deltaVirialTableNumberPoints))
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Virial_Density_Contrast_Tabulate(time,deltaVirialTableNumberPoints,deltaVirialTableTime,deltaVirialTableDeltaVirial)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
       tablesInitialized =.true.
    end if
    return
  end subroutine Virial_Density_Contrast_Retabulate

  double precision function Halo_Virial_Density_Contrast(time,aExpansion,collapsing)
    !% Return the halo virial overdensity.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: aExpansion,time
    logical,          intent(in), optional :: collapsing
    logical                                :: collapsingActual
    double precision                       :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Halo_Virial_Density_Contrast','only one argument can be specified')
       else
          timeActual=time
       end if
    else
       if (present(aExpansion)) then
          if (present(collapsing)) then
             collapsingActual=collapsing
          else
             collapsingActual=.false.
          end if
          timeActual=Cosmology_Age(aExpansion,collapsingActual)
       else
          call Galacticus_Error_Report('Halo_Virial_Density_Contrast','at least one argument must be given')
       end if
    end if

    ! Remake the table if necessary.
    call Virial_Density_Contrast_Retabulate(timeActual)

    ! Interpolate to get the expansion factor.
    Halo_Virial_Density_Contrast=Interpolate(deltaVirialTableNumberPoints,deltaVirialTableTime,deltaVirialTableDeltaVirial &
         &,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    return
  end function Halo_Virial_Density_Contrast

  double precision function Halo_Virial_Density_Contrast_Rate_of_Change(time,aExpansion,collapsing)
    !% Return the halo virial overdensity.
    use Numerical_Interpolation
    use Cosmology_Functions
    use Galacticus_Error
    implicit none
    double precision, intent(in), optional :: aExpansion,time
    logical,          intent(in), optional :: collapsing
    logical                                :: collapsingActual
    double precision                       :: timeActual

    ! Determine which type of input we have.
    if (present(time)) then
       if (present(aExpansion)) then
          call Galacticus_Error_Report('Halo_Virial_Density_Contrast_Rate_of_Change','only one argument can be specified')
       else
          timeActual=time
       end if
    else
       if (present(aExpansion)) then
          if (present(collapsing)) then
             collapsingActual=collapsing
          else
             collapsingActual=.false.
          end if
          timeActual=Cosmology_Age(aExpansion,collapsingActual)
       else
          call Galacticus_Error_Report('Halo_Virial_Density_Contrast_Rate_of_Change','at least one argument must be given')
       end if
    end if

    ! Remake the table if necessary.
    call Virial_Density_Contrast_Retabulate(timeActual)

    ! Interpolate to get the expansion factor.
    Halo_Virial_Density_Contrast_Rate_of_Change=Interpolate_Derivative(deltaVirialTableNumberPoints,deltaVirialTableTime&
         &,deltaVirialTableDeltaVirial ,interpolationObject,interpolationAccelerator,timeActual,reset=resetInterpolation)
    return
  end function Halo_Virial_Density_Contrast_Rate_of_Change

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Virial_Density_Contrast_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Virial_Density_Contrast_State_Retrieve(stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use Memory_Management
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    deltaVirialTableNumberPoints=0
    if (allocated(deltaVirialTableTime       )) call Dealloc_Array(deltaVirialTableTime       )
    if (allocated(deltaVirialTableDeltaVirial)) call Dealloc_Array(deltaVirialTableDeltaVirial)
    tablesInitialized=.false.
    return
  end subroutine Virial_Density_Contrast_State_Retrieve
  
end module Virial_Density_Contrast
