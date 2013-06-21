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

!% Contains a module which generates a tabulated transfer function using {\sc CMBFast}.

module Transfer_Function_CMBFast
  !% Implements generation of a tabulated transfer function using {\sc CMBFast}.
  use ISO_Varying_String
  implicit none
  private
  public :: Transfer_Function_CMBFast_Initialize
  
  ! Flag to indicate if this module has been initialized.
  logical                                     :: transferFunctionInitialized=.false.     
  
  ! File name for the transfer function data.
  type            (varying_string)            :: transferFunctionFile                    
  
  ! Smallest maximum wavenumber to tabulate.
  double precision                , parameter :: logWavenumberMaximumDefault=log(10.0d0) 
  
contains
  
  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_CMBFast_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_CMBFast_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``transfer function from CMBFast'' module.
    implicit none
    type     (varying_string                ), intent(in   )          :: transferFunctionMethod     
    procedure(Transfer_Function_CMBFast_Make), intent(inout), pointer :: Transfer_Function_Tabulate 
    
    if (transferFunctionMethod == 'CMBFast') Transfer_Function_Tabulate => Transfer_Function_CMBFast_Make
    return
  end subroutine Transfer_Function_CMBFast_Initialize

  subroutine Transfer_Function_CMBFast_Make(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Build a transfer function using {\sc CMBFast}.
    use FoX_wxml
    use System_Command
    use Transfer_Functions_File
    use Input_Parameters
    use Cosmological_Parameters
    use Numerical_Constants_Astronomical
    use Galacticus_Input_Paths
    use String_Handling
    implicit none
    double precision                                           , intent(in   ) :: logWavenumber                                               
    double precision                , allocatable, dimension(:), intent(inout) :: transferFunctionLogT        , transferFunctionLogWavenumber 
    integer                                                    , intent(  out) :: transferFunctionNumberPoints                                
    logical                                                                    :: makeFile                                                    
    character       (len=32        )                                           :: parameterLabel              , wavenumberLabel               
    type            (varying_string)                                           :: command                     , parameterFile                 
    type            (xmlf_t        )                                           :: parameterDoc                                                
    
    ! Generate the name of the data file and an XML input parameter file.    !# <uniqueLabel>    !#  <function>Transfer_Function_CMBFast_Label</function>    !#  <ignore>transferFunctionFile</ignore>    !# </uniqueLabel>
    transferFunctionFile=char(Galacticus_Input_Path())//'data/largeScaleStructure/transfer_function_CMBFast_'//Transfer_Function_CMBFast_Label(includeVersion=.true.,asHash=.true.)//".xml"
    parameterFile=char(Galacticus_Input_Path())//'data/transfer_function_parameters.xml'
    call xml_OpenFile(char(parameterFile),parameterDoc)
    call xml_NewElement(parameterDoc,"parameters")
    call xml_NewElement(parameterDoc,"uniqueLabel")
    call xml_AddCharacters(parameterDoc,char(Transfer_Function_CMBFast_Label(includeVersion=.true.)))
    call xml_EndElement(parameterDoc,"uniqueLabel")
    write (parameterLabel,'(f5.3)') Omega_Matter()
    call Write_Parameter(parameterDoc,"Omega_Matter",parameterLabel)
    write (parameterLabel,'(f5.3)') Omega_DE()
    call Write_Parameter(parameterDoc,"Omega_DE",parameterLabel)
    write (parameterLabel,'(f6.4)') Omega_b()
    call Write_Parameter(parameterDoc,"Omega_b",parameterLabel)
    write (parameterLabel,'(f4.1)') H_0()
    call Write_Parameter(parameterDoc,"H_0",parameterLabel)
    write (parameterLabel,'(f5.3)') T_CMB()
    call Write_Parameter(parameterDoc,"T_CMB",parameterLabel)
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    call Write_Parameter(parameterDoc,"Y_He",parameterLabel)
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
       write (wavenumberLabel,'(e12.6)') exp(max(logWavenumber+1.0d0,logWavenumberMaximumDefault))
       command=char(Galacticus_Input_Path())//'scripts/aux/CMBFast_Driver.pl '//parameterFile//' '//transferFunctionFile//' '//trim(wavenumberLabel)//' '//Transfer_Function_Named_File_Format_Version()
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
