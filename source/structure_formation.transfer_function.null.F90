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

!% Contains a module which generates a null transfer function.

module Transfer_Function_Null
  !% Implements a null transfer function.
  implicit none
  private
  public :: Transfer_Function_Null_Initialize

  ! Wavenumber range and fineness of gridding.
  double precision            :: logWavenumberMaximum=log(10.0d0)
  double precision            :: logWavenumberMinimum=log(1.0d-5)
  integer         , parameter :: numberPointsPerDecade=1000

contains
  
  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_Null_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_Null_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``null transfer function'' module.
    use ISO_Varying_String
    implicit none
    type     (varying_string             ),          intent(in   ) :: transferFunctionMethod
    procedure(Transfer_Function_Null_Make), pointer, intent(inout) :: Transfer_Function_Tabulate
    
    if (transferFunctionMethod == 'null') Transfer_Function_Tabulate => Transfer_Function_Null_Make
    return
  end subroutine Transfer_Function_Null_Initialize

  subroutine Transfer_Function_Null_Make(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Build a null transfer function.
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Constants_Math
    implicit none
    double precision,                            intent(in   ) :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                     intent(  out) :: transferFunctionNumberPoints
 
    ! Set wavenumber range and number of points in table.
    logWavenumberMinimum=min(logWavenumberMinimum,logWavenumber-ln10)
    logWavenumberMaximum=max(logWavenumberMaximum,logWavenumber+ln10)
    transferFunctionNumberPoints=int((logWavenumberMaximum-logWavenumberMinimum)*dble(numberPointsPerDecade)/ln10)
    ! Deallocate arrays if currently allocated.
    if (allocated(transferFunctionLogWavenumber)) call Dealloc_Array(transferFunctionLogWavenumber)
    if (allocated(transferFunctionLogT))          call Dealloc_Array(transferFunctionLogT         )
    ! Allocate the arrays to current required size.
    call Alloc_Array(transferFunctionLogWavenumber,[transferFunctionNumberPoints])
    call Alloc_Array(transferFunctionLogT         ,[transferFunctionNumberPoints])
    ! Create range of wavenumbers.
    transferFunctionLogWavenumber=Make_Range(logWavenumberMinimum,logWavenumberMaximum,transferFunctionNumberPoints&
         &,rangeTypeLinear)
    ! Create transfer function.
    transferFunctionLogT=0.0d0
    return
  end subroutine Transfer_Function_Null_Make
  
end module Transfer_Function_Null
