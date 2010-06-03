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








!% Contains a module which reads a tabulated transfer function from a file.

module Transfer_Function_File
  !% Implements reading of a tabulated transfer function from a file.
  use ISO_Varying_String
  private
  public :: Transfer_Function_File_Initialize, Transfer_Function_Named_File_Read
  
  ! Flag to indicate if this module has been initialized.
  logical              :: transferFunctionInitialized=.false.

  ! File name for the transfer function data.
  type(varying_string) :: transferFunctionFile

contains
  
  !# <transferFunctionMethod>
  !#  <unitName>Transfer_Function_File_Initialize</unitName>
  !# </transferFunctionMethod>
  subroutine Transfer_Function_File_Initialize(transferFunctionMethod,Transfer_Function_Tabulate)
    !% Initializes the ``transfer function from file'' module.
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: transferFunctionMethod
    procedure(),          pointer, intent(inout) :: Transfer_Function_Tabulate
    
    if (transferFunctionMethod.eq.'file') then
       Transfer_Function_Tabulate => Transfer_Function_File_Read
       !@ <inputParameter>
       !@   <name>transferFunctionFile</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of a file containing a tabulation of the transfer function for the ``file'' transfer function method.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('transferFunctionFile',transferFunctionFile)
    end if
    return
  end subroutine Transfer_Function_File_Initialize

  subroutine Transfer_Function_Named_File_Read(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT,fileName)
    !% Read the transfer function from a file specified by {\tt fileName}.
    implicit none
    double precision,                                intent(in)    :: logWavenumber
    double precision,     allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                         intent(out)   :: transferFunctionNumberPoints
    type(varying_string),                            intent(in)    :: fileName

    ! Set the filename to that specified.
    transferFunctionFile=fileName
    ! Flag that module is uninitialized again.
    transferFunctionInitialized=.false.
    ! Call routine to read in the data.
    call Transfer_Function_File_Read(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
         &,transferFunctionLogT)
   return
  end subroutine Transfer_Function_Named_File_Read

  subroutine Transfer_Function_File_Read(logWavenumber,transferFunctionNumberPoints,transferFunctionLogWavenumber&
       &,transferFunctionLogT)
    !% Reads a transfer function from an XML file.
    use FoX_dom
    use Numerical_Comparison
    use Memory_Management
    use Galacticus_Error
    use Galacticus_Display
    use Cosmological_Parameters
    implicit none
    double precision,                            intent(in)    :: logWavenumber
    double precision, allocatable, dimension(:), intent(inout) :: transferFunctionLogWavenumber,transferFunctionLogT
    integer,                                     intent(out)   :: transferFunctionNumberPoints
    type(Node),       pointer                                  :: doc,datum,thisParameter,nameElement,valueElement
    type(NodeList),   pointer                                  :: datumList,parameterList
    integer                                                    :: iDatum,ioErr,iParameter
    double precision                                           :: datumValues(2),parameterValue
    
    ! Read the file if this module has not been initialized.
    if (.not.transferFunctionInitialized) then
       ! Open and parse the data file.
       !$omp critical (FoX_DOM_Access)
       doc => parseFile(char(transferFunctionFile),iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Transfer_Function_File_Read','Unable to find transfer function file')
       ! Check that parameters match if any are present.
       parameterList => getElementsByTagname(doc,"parameter")
       do iParameter=0,getLength(parameterList)-1
          thisParameter => item(parameterList,iParameter)
          nameElement => item(getElementsByTagname(thisParameter,"name"),0)
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          call extractDataContent(valueElement,parameterValue)
          select case (getTextContent(nameElement))
          case ("Omega_b")
             if (Values_Differ(parameterValue,Omega_b(),absTol=1.0d-3)) call Galacticus_Display_Message('Omega_b from transfer &
                  & function file does not match internal value')
          case ("Omega_0")
             if (Values_Differ(parameterValue,Omega_0(),absTol=1.0d-3)) call Galacticus_Display_Message('Omega_0 from transfer &
                  & function file does not match internal value')
          case ("Omega_DE")
             if (Values_Differ(parameterValue,Omega_DE(),absTol=1.0d-3)) call Galacticus_Display_Message('Omega_DE from transfer &
                  & function file does not match internal value')
          case ("H_0")
             if (Values_Differ(parameterValue,H_0(),relTol=1.0d-3)) call Galacticus_Display_Message('H_0 from transfer &
                  & function file does not match internal value')
          case ("T_CMB")
             if (Values_Differ(parameterValue,T_CMB(),relTol=1.0d-3)) call Galacticus_Display_Message('T_CMB from transfer &
                  & function file does not match internal value')
          end select
       end do
       ! Get list of datum elements.
       datumList => getElementsByTagname(doc,"datum")
       ! Allocate transfer function arrays.
       transferFunctionNumberPoints=getLength(datumList)
       ! Deallocate arrays if currently allocated.
       if (allocated(transferFunctionLogWavenumber)) call Dealloc_Array(transferFunctionLogWavenumber)
       if (allocated(transferFunctionLogT))          call Dealloc_Array(transferFunctionLogT         )
       ! Allocate the arrays to current required size.
       call Alloc_Array(transferFunctionLogWavenumber,transferFunctionNumberPoints,'transferFunctionLogWavenumber')
       call Alloc_Array(transferFunctionLogT         ,transferFunctionNumberPoints,'transferFunctionLogT'         )
       ! Extract data from elements.
       do iDatum=0,getLength(datumList)-1
          datum => item(datumList,iDatum)
          call extractDataContent(datum,datumValues)
          transferFunctionLogWavenumber(iDatum+1)=dlog(datumValues(1))
          transferFunctionLogT         (iDatum+1)=dlog(datumValues(2))
       end do
       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)
       ! Flag that transfer function is now initialized.
       transferFunctionInitialized=.true.
    end if
    ! Check that the input wavenumber is within range.
    if (logWavenumber < transferFunctionLogWavenumber(1) .or. logWavenumber >&
         & transferFunctionLogWavenumber(transferFunctionNumberPoints)) call&
         & Galacticus_Error_Report('Transfer_Function_File_Read','wavenumber is not within tabulated range') 
    return
  end subroutine Transfer_Function_File_Read
  
end module Transfer_Function_File
