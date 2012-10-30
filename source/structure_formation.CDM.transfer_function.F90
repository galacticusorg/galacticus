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

module CDM_Transfer_Function
  use ISO_Varying_String
  use FGSL
  implicit none
  private
  public :: Transfer_Function_CDM, CDM_Transfer_Function_State_Retrieve

  ! Flag to indicate if this module has been initialized.  
  logical                                        :: transferFunctionInitialized=.false., tablesInitialized=.false.

  ! Variables to hold the tabulated transfer function data.
  integer                                        :: transferFunctionNumberPoints=-1
  double precision,    allocatable, dimension(:) :: transferFunctionLogWavenumber,transferFunctionLogT
  type(fgsl_interp)                              :: interpolationObject
  type(fgsl_interp_accel)                        :: interpolationAccelerator
  logical                                        :: resetInterpolation=.true.

  ! Name of transfer function method used.
  type(varying_string)                           :: transferFunctionMethod

  ! Pointer to the subroutine that tabulates the transfer function and template interface for that subroutine.
  procedure(Transfer_Function_Tabulate_Template), pointer :: Transfer_Function_Tabulate => null()
  interface Transfer_Function_Tabulate_Template
     subroutine Transfer_Function_Tabulate_Template(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber &
          &,transferFunctionLogT)
    double precision,                            intent(in)    :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                     intent(out)   :: transferFunctionNumberPoints
  end subroutine Transfer_Function_Tabulate_Template
 end interface
  
contains

  double precision function Transfer_Function_CDM(wavenumber)
    !% Return the CDM transfer function for $k=${\tt wavenumber} [Mpc$^{-1}$].
    use Numerical_Interpolation
    implicit none
    double precision, intent(in) :: wavenumber
    double precision             :: logWavenumber

    ! Get logarithm of wavenumber.
    logWavenumber=dlog(wavenumber)

    !$omp critical(Transfer_Function_Initialization) 
    ! Initialize if necessary.
    if (.not.(transferFunctionInitialized.and.tablesInitialized)) then
       call Transfer_Function_Initialize(logWavenumber)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
    end if

    ! If wavenumber is out of range, attempt to remake the table.
    if (logWavenumber<transferFunctionLogWavenumber(1) .or. logWavenumber&
         &>transferFunctionLogWavenumber(transferFunctionNumberPoints) ) then
       call Transfer_Function_Tabulate(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
            &,transferFunctionLogT)
       call Interpolate_Done(interpolationObject,interpolationAccelerator,resetInterpolation)
       resetInterpolation=.true.
    end if
    !$omp end critical(Transfer_Function_Initialization)

    ! Interpolate in the tabulated function and return a value.
    Transfer_Function_CDM=dexp(Interpolate(transferFunctionNumberPoints,transferFunctionLogWavenumber,transferFunctionLogT &
         &,interpolationObject,interpolationAccelerator,logWavenumber,reset=resetInterpolation,interpolationType=fgsl_interp_cspline))

    return
  end function Transfer_Function_CDM

  subroutine Transfer_Function_Initialize(logWavenumber)
    !% Initializes the transfer function module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="transferFunctionMethod" type="moduleUse">
    include 'structure_formation.CDM.transfer_function.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: logWavenumber

    if (.not.transferFunctionInitialized) then
       ! Get the transfer function method parameter.
       !@ <inputParameter>
       !@   <name>transferFunctionMethod</name>
       !@   <defaultValue>Eisenstein-Hu1999</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing the transfer function.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('transferFunctionMethod',transferFunctionMethod,defaultValue='Eisenstein-Hu1999')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="transferFunctionMethod" type="functionCall" functionType="void">
       !#  <functionArgs>transferFunctionMethod,Transfer_Function_Tabulate</functionArgs>
       include 'structure_formation.CDM.transfer_function.inc'
       !# </include>
       if (.not.associated(Transfer_Function_Tabulate)) call Galacticus_Error_Report('Transfer_Function_Initialize','method '&
            &//char(transferFunctionMethod)//' is unrecognized')
    end if
    ! Call routine to initialize the transfer function.
    call Transfer_Function_Tabulate(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber,transferFunctionLogT)
    tablesInitialized=.true.
    ! Flag that the module is now initialized.
    transferFunctionInitialized=.true.
    return
  end subroutine Transfer_Function_Initialize

  !# <galacticusStateRetrieveTask>
  !#  <unitName>CDM_Transfer_Function_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine CDM_Transfer_Function_State_Retrieve(stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved. This will force tables to be rebuilt.
    use Memory_Management
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    transferFunctionNumberPoints=0
    if (allocated(transferFunctionLogWavenumber)) call Dealloc_Array(transferFunctionLogWavenumber)
    if (allocated(transferFunctionLogT         )) call Dealloc_Array(transferFunctionLogT         )
    tablesInitialized=.false.
    return
  end subroutine CDM_Transfer_Function_State_Retrieve
  
end module CDM_Transfer_Function
