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


!% Contains a module which generates a tabulated transfer function using the BBKS fitting formula.

module Transfer_Function_BBKS
  !% Implements generation of a tabulated transfer function using the BBKS fitting formula.
  use ISO_Varying_String
  private
  public :: Transfer_Function_BBKS_Initialize, Transfer_Function_BBKS_State_Store, Transfer_Function_BBKS_State_Retrieve
  
  ! Flag to indicate if this module has been initialized.
  logical                     :: transferFunctionInitialized=.false.

  ! Wavenumber range and fineness of gridding.
  double precision            :: logWavenumberMaximum=dlog(10.0d0)
  double precision            :: logWavenumberMinimum=dlog(1.0d-5)
  integer,          parameter :: numberPointsPerDecade=1000

  ! Warm dark matter free-streaming length.
  double precision            :: transferFunctionWDMFreeStreamingLength

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
    
    if (transferFunctionMethod == 'BBKS') then
       ! Return a pointer to our tabulation function.
       Transfer_Function_Tabulate => Transfer_Function_BBKS_Make

       ! Get input parameters.
       !@ <inputParameter>
       !@   <name>transferFunctionWDMFreeStreamingLength</name>
       !@   <defaultValue>0</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The warm dark matter free streaming length (in Mpc).
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('transferFunctionWDMFreeStreamingLength',transferFunctionWDMFreeStreamingLength,defaultValue=0.0d0)
    end if
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
    double precision                                           :: Gamma,q,wavenumberHUnits,wavenumber,wavenumberScaleFree
 
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
    Gamma=Omega_0()*Little_H_0()*dexp(-Omega_b()*(1.0d0+dsqrt(2.0d0*Little_H_0())/Omega_0()))/((T_CMB()/2.7d0)**2)
    do iWavenumber=1,transferFunctionNumberPoints
       wavenumber         =dexp(transferFunctionLogWavenumber(iWavenumber))
       wavenumberHUnits   =wavenumber/Little_H_0()
       wavenumberScaleFree=wavenumber*transferFunctionWDMFreeStreamingLength
       q                  =wavenumberHUnits/Gamma
       transferFunctionLogT(iWavenumber)=dlog((dlog(1.0+2.34d0*q)/2.34d0/q)/(1.0d0+3.89d0*q+(16.1d0*q)**2+(5.46d0*q)**3+(6.71d0&
            &*q)**4)**0.25d0)-0.5d0*wavenumberScaleFree*(1.0d0+wavenumberScaleFree)
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
