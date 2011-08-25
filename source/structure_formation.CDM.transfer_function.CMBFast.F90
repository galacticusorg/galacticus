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


!% Contains a module which generates a tabulated transfer function using {\sc CMBFast}.

module Transfer_Function_CMBFast
  !% Implements generation of a tabulated transfer function using {\sc CMBFast}.
  use ISO_Varying_String
  implicit none
  private
  public :: Transfer_Function_CMBFast_Initialize
  
  ! Flag to indicate if this module has been initialized.
  logical              :: transferFunctionInitialized=.false.

  ! File name for the transfer function data.
  type(varying_string) :: transferFunctionFile

  ! Smallest maximum wavenumber to tabulate.
  double precision, parameter :: logWavenumberMaximumDefault=dlog(10.0d0)

contains
  
  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_CMBFast_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_CMBFast_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``transfer function from CMBFast'' module.
    implicit none
    type(varying_string),          intent(in)    :: transferFunctionMethod
    procedure(),          pointer, intent(inout) :: Transfer_Function_Tabulate
    
    if (transferFunctionMethod == 'CMBFast') Transfer_Function_Tabulate => Transfer_Function_CMBFast_Make
    return
  end subroutine Transfer_Function_CMBFast_Initialize

  subroutine Transfer_Function_CMBFast_Make(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Build a transfer function using {\sc CMBFast}.
    use FoX_wxml
    use System_Command
    use Transfer_Function_File
    use Input_Parameters
    use Cosmological_Parameters
    use Numerical_Constants_Astronomical
    implicit none
    double precision,                            intent(in)    :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                     intent(out)   :: transferFunctionNumberPoints
    logical                                                    :: makeFile
    character(len=32)                                          :: parameterLabel,wavenumberLabel
    type(varying_string)                                       :: parameterFile,command
    type(xmlf_t)                                               :: parameterDoc

    ! Generate the name of the data file and an XML input parameter file.
    transferFunctionFile='data/transfer_function_CMBFast'
    parameterFile='data/transfer_function_parameters.xml'
    call xml_OpenFile(char(parameterFile),parameterDoc)
    call xml_NewElement(parameterDoc,"parameters")
    write (parameterLabel,'(f5.3)') Omega_Matter()
    transferFunctionFile=transferFunctionFile//'_OmegaMatter'//trim(parameterLabel)
    call Write_Parameter(parameterDoc,"Omega_Matter",parameterLabel)
    write (parameterLabel,'(f5.3)') Omega_DE()
    transferFunctionFile=transferFunctionFile//'_OmegaDE'//trim(parameterLabel)
    call Write_Parameter(parameterDoc,"Omega_DE",parameterLabel)
    write (parameterLabel,'(f6.4)') Omega_b()
    transferFunctionFile=transferFunctionFile//'_Omegab'//trim(parameterLabel)
    call Write_Parameter(parameterDoc,"Omega_b",parameterLabel)
    write (parameterLabel,'(f4.1)') H_0()
    transferFunctionFile=transferFunctionFile//'_H0'//trim(parameterLabel)
    call Write_Parameter(parameterDoc,"H_0",parameterLabel)
    write (parameterLabel,'(f5.3)') T_CMB()
    transferFunctionFile=transferFunctionFile//'_TCMB'//trim(parameterLabel)
    call Write_Parameter(parameterDoc,"T_CMB",parameterLabel)
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    transferFunctionFile=transferFunctionFile//'_YHe'//trim(parameterLabel)
    call Write_Parameter(parameterDoc,"Y_He",parameterLabel)
    transferFunctionFile=transferFunctionFile//'.xml'
    call xml_Close(parameterDoc)
    ! Determine if we need to reinitialize this module.
    if (.not.transferFunctionInitialized) then
       makeFile=.true.
    else
       makeFile=min(logWavenumber,logWavenumberMaximumDefault) > transferFunctionLogWavenumber(transferFunctionNumberPoints)
       if (makeFile) then
          ! Remove the transfer function file so that a new one will be created.
          command='rm -f '//transferFunctionFile
          call System_Command_Do(command)
       end if
    end if
    ! Read the file if this module has not been initialized or if the wavenumber is out of range.
    if (makeFile) then
       ! Run CMBFast wrapper script.
       write (wavenumberLabel,'(e12.6)') dexp(max(logWavenumber+1.0d0,logWavenumberMaximumDefault))
       command='./scripts/aux/CMBFast_Driver.pl '//parameterFile//' '//transferFunctionFile//' '//trim(wavenumberLabel)
       call System_Command_Do(command)

       ! Flag that transfer function is now initialized.
       transferFunctionInitialized=.true.
    end if

    ! Call routine to read in the tabulated data.
    call Transfer_Function_Named_File_Read(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber &
         &,transferFunctionLogT,transferFunctionFile)

    ! Remove the parameter file.
    command='rm -f '//parameterFile
    call System_Command_Do(command)
    return
  end subroutine Transfer_Function_CMBFast_Make
  
end module Transfer_Function_CMBFast
