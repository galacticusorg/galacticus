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








!% Contains a module which generates a tabulated transfer function using the BBKS fitting formula.

module Transfer_Function_BBKS
  !% Implements generation of a tabulated transfer function using the BBKS fitting formula.
  use ISO_Varying_String
  private
  public :: Transfer_Function_BBKS_Initialize, Transfer_Function_BBKS_State_Store, Transfer_Function_BBKS_State_Retrieve
  
  ! Flag to indicate if this module has been initialized.
  logical              :: transferFunctionInitialized=.false.

  ! avenumber range and fineness of gridding.
  double precision            :: logWavenumberMaximum=dlog(10.0d0)
  double precision            :: logWavenumberMinimum=dlog(1.0d-5)
  integer,          parameter :: numberPointsPerDecade=1000

contains
  
  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_BBKS_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_BBKS_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``transfer function from BBKS'' module.
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: transferFunctionMethod
    procedure(),          pointer, intent(inout) :: Transfer_Function_Tabulate
    
    if (transferFunctionMethod.eq.'BBKS') Transfer_Function_Tabulate => Transfer_Function_BBKS_Make
    return
  end subroutine Transfer_Function_BBKS_Initialize

  subroutine Transfer_Function_BBKS_Make(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Build a transfer function using the BBKS fitting formula.
    use Memory_Management
    use Cosmological_Parameters
    use Numerical_Ranges
    use Numerical_Constants_Math
    implicit none
    double precision,                            intent(in)    :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                     intent(out)   :: transferFunctionNumberPoints
    integer                                                    :: iWavenumber
    double precision                                           :: Gamma,q,wavenumberHUnits
 
    ! Set wavenumber range and number of points in table.
    logWavenumberMinimum=min(logWavenumberMinimum,logWavenumber-ln10)
    logWavenumberMaximum=max(logWavenumberMaximum,logWavenumber+ln10)
    transferFunctionNumberPoints=int((logWavenumberMaximum-logWavenumberMinimum)*dble(numberPointsPerDecade)/ln10)
    ! Deallocate arrays if currently allocated.
    if (allocated(transferFunctionLogWavenumber)) call Dealloc_Array(transferFunctionLogWavenumber)
    if (allocated(transferFunctionLogT))          call Dealloc_Array(transferFunctionLogT         )
    ! Allocate the arrays to current required size.
    call Alloc_Array(transferFunctionLogWavenumber,transferFunctionNumberPoints,'transferFunctionLogWavenumber')
    call Alloc_Array(transferFunctionLogT         ,transferFunctionNumberPoints,'transferFunctionLogT'         )
    ! Create range of wavenumbers.
    transferFunctionLogWavenumber=Make_Range(logWavenumberMinimum,logWavenumberMaximum,transferFunctionNumberPoints&
         &,rangeTypeLinear)
    ! Create transfer function.
    Gamma=Omega_0()*Little_H_0()*dexp(-Omega_b()*(1.0d0+dsqrt(2.0d0*Little_H_0())/Omega_0()))/((T_CMB()/2.7d0)**2)
    do iWavenumber=1,transferFunctionNumberPoints
       wavenumberHUnits=dexp(transferFunctionLogWavenumber(iWavenumber))/Little_H_0()
       q=wavenumberHUnits/Gamma
       transferFunctionLogT(iWavenumber)=dlog((dlog(1.0+2.34d0*q)/2.34d0/q)/(1.0d0+3.89d0*q+(16.1d0*q)**2+(5.46d0*q)**3+(6.71d0&
            &*q)**4)**0.25d0)
    end do
    return
  end subroutine Transfer_Function_BBKS_Make
  
  !# <galacticusStateStoreTask>
  !#  <unitName>Transfer_Function_BBKS_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Transfer_Function_BBKS_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) logWavenumberMinimum,logWavenumberMaximum
    return
  end subroutine Transfer_Function_BBKS_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Transfer_Function_BBKS_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Transfer_Function_BBKS_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) logWavenumberMinimum,logWavenumberMaximum
    return
  end subroutine Transfer_Function_BBKS_State_Retrieve
    
end module Transfer_Function_BBKS
