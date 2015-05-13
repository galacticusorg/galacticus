!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which generates a tabulated transfer function using \href{http://camb.info}{\normalfont \scshape CAMB}.

module Transfer_Function_CAMB
  !% Implements generation of a tabulated transfer function using \href{http://camb.info}{\normalfont \scshape CAMB}.
  use ISO_Varying_String
  implicit none
  private
  public :: Transfer_Function_CAMB_Initialize

  ! Flag to indicate if this module has been initialized.
  logical                                     :: transferFunctionInitialized=.false.

  ! File name for the transfer function data.
  type            (varying_string)            :: transferFunctionFile

  ! Smallest maximum wavenumber to tabulate.
  double precision                , parameter :: logWavenumberMaximumDefault=log(10.0d0)

contains

  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_CAMB_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_CAMB_Initialize(transferFunctionMethod,Transfer_Function_Tabulate&
       &,Transfer_Function_Half_Mode_Mass)
    implicit none
    type     (varying_string                       ), intent(in   )          :: transferFunctionMethod
    procedure(Transfer_Function_CAMB_Make       ), intent(inout), pointer :: Transfer_Function_Tabulate
    procedure(Transfer_Function_Half_Mode_Mass_Null), intent(inout), pointer :: Transfer_Function_Half_Mode_Mass

    if (transferFunctionMethod == 'CAMB') then
       Transfer_Function_Tabulate       => Transfer_Function_CAMB_Make
       Transfer_Function_Half_Mode_Mass => Transfer_Function_Half_Mode_Mass_Null
    end if
    return
  end subroutine Transfer_Function_CAMB_Initialize

  subroutine Transfer_Function_CAMB_Make(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Build a transfer function using \href{http://camb.info}{\normalfont \scshape CAMB}.
    use FoX_wxml
    use System_Command
    use Transfer_Functions_File
    use Input_Parameters
    use Input_Parameters2
    use Cosmology_Parameters
    use Numerical_Constants_Astronomical
    use Galacticus_Input_Paths
    use String_Handling
    implicit none
    double precision                                               , intent(in   ) :: logWavenumber
    double precision                    , allocatable, dimension(:), intent(inout) :: transferFunctionLogT        , transferFunctionLogWavenumber
    integer                                                        , intent(  out) :: transferFunctionNumberPoints
    logical                                                                        :: makeFile
    character       (len=32            )                                           :: parameterLabel              , wavenumberLabel
    type            (varying_string    )                                           :: command                     , parameterFile
    type            (xmlf_t            )                                           :: parameterDoc
    type            (inputParameterList)                                           :: dependentParameters

    ! Generate the name of the data file and an XML input parameter file.
    !# <uniqueLabel>
    !#  <function>Transfer_Function_CAMB_Label</function>
    !#  <ignore>transferFunctionFile</ignore>
    !# </uniqueLabel>
    transferFunctionFile=char(Galacticus_Input_Path())//'data/largeScaleStructure/transfer_function_CAMB_'//Transfer_Function_CAMB_Label(includeSourceDigest=.true.,asHash=.true.,parameters=dependentParameters)//".xml"
    parameterFile=char(Galacticus_Input_Path())//'data/transfer_function_parameters.xml'
    call xml_OpenFile(char(parameterFile),parameterDoc)
    call xml_NewElement(parameterDoc,"parameters")
    call xml_NewElement(parameterDoc,"uniqueLabel")
    call xml_AddCharacters(parameterDoc,char(Transfer_Function_CAMB_Label(includeSourceDigest=.true.)))
    call xml_EndElement(parameterDoc,"uniqueLabel")
    write (parameterLabel,'(f4.2)') heliumByMassPrimordial
    call dependentParameters%add("Y_He",parameterLabel)
    call dependentParameters%serializeToXML(parameterDoc)
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
       ! Run CAMB wrapper script.
       write (wavenumberLabel,'(e12.6)') exp(max(logWavenumber+1.0d0,logWavenumberMaximumDefault))
       command=char(Galacticus_Input_Path())//'scripts/aux/CAMB_Driver.pl '//parameterFile//' '//transferFunctionFile//' '//trim(wavenumberLabel)//' '//Transfer_Function_Named_File_Format_Version()
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
  end subroutine Transfer_Function_CAMB_Make

  double precision function Transfer_Function_Half_Mode_Mass_Null()
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function. Not supported in this implementation.
    use Galacticus_Error
    implicit none

    call Galacticus_Error_Report('Transfer_Function_Half_Mode_Mass_Null','not supported by this implementation')
    return
  end function Transfer_Function_Half_Mode_Mass_Null

end module Transfer_Function_CAMB
