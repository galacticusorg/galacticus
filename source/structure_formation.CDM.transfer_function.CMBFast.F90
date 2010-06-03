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








!% Contains a module which generates a tabulated transfer function using {\sc CMBFast}.

module Transfer_Function_CMBFast
  !% Implements generation of a tabulated transfer function using {\sc CMBFast}.
  use ISO_Varying_String
  private
  public :: Transfer_Function_CMBFast_Initialize
  
  ! Flag to indicate if this module has been initialized.
  logical              :: transferFunctionInitialized=.false.

  ! File name for the transfer function data.
  type(varying_string) :: transferFunctionFile

  ! Helium mass fraction.
  double precision     :: Y_He

  ! Smallest maximum wavenumber to tabulate.
  double precision, parameter :: logWavenumberMaximumDefault=dlog(10.0d0)

contains
  
  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_CMBFast_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_CMBFast_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``transfer function from CMBFast'' module.
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: transferFunctionMethod
    procedure(),          pointer, intent(inout) :: Transfer_Function_Tabulate
    
    if (transferFunctionMethod == 'CMBFast') then
       Transfer_Function_Tabulate => Transfer_Function_CMBFast_Make
       !@ <inputParameter>
       !@   <name>Y_He</name>
       !@   <attachedTo>module</attachedTo>
       !@   <defaultValue>0.2477 \citep{peimbert_primordial_2008}</defaultValue>
       !@   <description>
       !@     The mass fraction of helium in the primordial plasma for calculations of the transfer function using {\sc CMBFast}.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('Y_He',Y_He,defaultValue=0.2477d0)
    end if
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
    write (parameterLabel,'(f5.3)') Omega_0()
    transferFunctionFile=transferFunctionFile//'_Omega0'//trim(parameterLabel)
    call Write_Parameter(parameterDoc,"Omega_0",parameterLabel)
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
    write (parameterLabel,'(f4.2)') Y_He
    transferFunctionFile=transferFunctionFile//'_YHe'//trim(parameterLabel)
    call Write_Parameter(parameterDoc,"Y_He",parameterLabel)
    transferFunctionFile=transferFunctionFile//'.xml'
    call xml_Close(parameterDoc)
    ! Determine if we need to reinitialize this module.
    if (.not.transferFunctionInitialized) then
       makeFile=.true.
    else
       makeFile=(logWavenumber < transferFunctionLogWavenumber(1)) .or. (logWavenumber >&
            & transferFunctionLogWavenumber(transferFunctionNumberPoints))
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

       ! Call routine to read in the tabulated data.
       call Transfer_Function_Named_File_Read(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber &
            &,transferFunctionLogT,transferFunctionFile)

       ! Flag that transfer function is now initialized.
       transferFunctionInitialized=.true.
    end if

    ! Remove the parameter file.
    command='rm -f '//parameterFile
    call System_Command_Do(command)
    return
  end subroutine Transfer_Function_CMBFast_Make
  
end module Transfer_Function_CMBFast
