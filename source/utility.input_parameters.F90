!! Copyright 2009, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements reading of parameters from an XML data file.

module Input_Parameters
  !% Implements reading of parameters from an XML data file.
  use FoX_dom
  use HDF5
  use ISO_Varying_String
  use Galacticus_Error
  use Galacticus_HDF5_Groups
  private
  public :: Input_Parameters_File_Open, Input_Parameters_File_Close, Get_Input_Parameter, Get_Input_Parameter_Array_Size, Write_Parameter

  ! Node to hold the parameter document.
  type(Node),     pointer :: parameterDoc => null()
  type(NodeList), pointer :: parameterList
  integer                 :: parameterCount

  ! Generic interface to parameter value functions.
  interface Get_Input_Parameter
     module procedure Get_Input_Parameter_Char
     module procedure Get_Input_Parameter_Char_Array
     module procedure Get_Input_Parameter_VarString
     module procedure Get_Input_Parameter_VarString_Array
     module procedure Get_Input_Parameter_Double
     module procedure Get_Input_Parameter_Double_Array
     module procedure Get_Input_Parameter_Integer
     module procedure Get_Input_Parameter_Integer_Array
     module procedure Get_Input_Parameter_Logical
     module procedure Get_Input_Parameter_Logical_Array
  end interface

  ! Parameters group identifier in the output file.
  integer :: parametersGroupID=0
  
contains

  subroutine Input_Parameters_File_Open(parameterFile)
    !% Open an XML data file containing parameter values and parse it. The file should be structured as follows:
    !% \begin{verbatim}
    !% <parameters>
    !%   <parameter>
    !%     <name>parameter1Name</name>
    !%     <value>parameter1Value</value>
    !%   </parameter>
    !%   <parameter>
    !%     <name>parameter1Name</name>
    !%     <value>parameter1Value</value>
    !%   </parameter>
    !%   .
    !%   .
    !%   .
    !% </parameters>
    !% \end{verbatim}
    implicit none
    type(varying_string), intent(in) :: parameterFile
    integer                          :: ioErr

    ! Open and parse the data file.
    !$omp critical (FoX_DOM_Access)
    parameterDoc => parseFile(char(parameterFile),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Input_Parameters_File_Open','Unable to find parameter file')
    parameterList => getElementsByTagname(parameterDoc,"parameter")
    parameterCount=getLength(parameterList)
    !$omp end critical (FoX_DOM_Access)
    return
  end subroutine Input_Parameters_File_Open

  subroutine Input_Parameters_File_Close
    !% Close the parameter file (actually just destroy the internal record of it and clean up memory).
    implicit none

    call destroy(parameterDoc)
    return
  end subroutine Input_Parameters_File_Close

  integer function Get_Input_Parameter_Array_Size(parameterName)
    !% Get the number of elements in the parameter specified by parameter name is specified by {\tt parameterName}.
    use String_Handling
    implicit none
    character(len=*),     intent(in)           :: parameterName
    type(Node),           pointer              :: thisParameter,nameElement,valueElement
    integer                                    :: iParameter
    logical                                    :: foundMatch
    type(varying_string)                       :: parameterText

    ! If no parameter file has been read return zero.
    if (.not.associated(parameterDoc)) then
      Get_Input_Parameter_Array_Size=0
      return
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          parameterText=getTextContent(valueElement)
          Get_Input_Parameter_Array_Size=String_Count_Words(char(parameterText))
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       Get_Input_Parameter_Array_Size=0
       return
    end if
    return
  end function Get_Input_Parameter_Array_Size

  subroutine Get_Input_Parameter_Char(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt varying\_string} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*),     intent(out)          :: parameterValue
    character(len=*),     intent(in)           :: parameterName
    character(len=*),     intent(in), optional :: defaultValue
    logical,              intent(in), optional :: writeOutput
    type(Node),           pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                        :: datasetID
    integer                                    :: iParameter
    logical                                    :: foundMatch,writeOutputActual
    character(len=len(parameterValue))         :: datasetChar(1)

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Char','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          parameterValue=getTextContent(valueElement)
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Char','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       datasetChar=parameterValue
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',datasetChar)
    end if
    return
  end subroutine Get_Input_Parameter_Char

  subroutine Get_Input_Parameter_Char_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt varying\_string} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    use String_Handling
    implicit none
    character(len=*),     intent(out)          :: parameterValue(:)
    character(len=*),     intent(in)           :: parameterName
    character(len=*),     intent(in), optional :: defaultValue(:)
    logical,              intent(in), optional :: writeOutput
    type(Node),           pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                        :: datasetID
    integer                                    :: iParameter,nEntries
    logical                                    :: foundMatch,writeOutputActual
    type(varying_string)                       :: parameterText

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Char_Array','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          parameterText=getTextContent(valueElement)
          nEntries=String_Count_Words(char(parameterText))
          if (nEntries > size(parameterValue)) then
             call Galacticus_Error_Report('Get_Input_Parameter_Char_Array','array parameter has too many entries')
          else
             call String_Split_Words(parameterValue(1:nEntries),char(parameterText))
          end if
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Char_Array','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',parameterValue)
    end if
    return
  end subroutine Get_Input_Parameter_Char_Array

  subroutine Get_Input_Parameter_VarString(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt varying\_string} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    type(varying_string), intent(out)          :: parameterValue
    character(len=*),     intent(in)           :: parameterName
    character(len=*),     intent(in), optional :: defaultValue
    logical,              intent(in), optional :: writeOutput
    type(Node),           pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                        :: datasetID
    integer                                    :: iParameter
    logical                                    :: foundMatch,writeOutputActual
    type(varying_string)                       :: datasetVarString(1)

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_VarString','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          parameterValue=getTextContent(valueElement)
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_VarString','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       datasetVarString=parameterValue
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',datasetVarString)
    end if
    return
  end subroutine Get_Input_Parameter_VarString

  subroutine Get_Input_Parameter_VarString_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt varying\_string} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    use String_Handling
    implicit none
    type(varying_string), intent(out)          :: parameterValue(:)
    character(len=*),     intent(in)           :: parameterName
    character(len=*),     intent(in), optional :: defaultValue(:)
    logical,              intent(in), optional :: writeOutput
    type(Node),           pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                        :: datasetID
    integer                                    :: iParameter,nEntries
    logical                                    :: foundMatch,writeOutputActual
    type(varying_string)                       ::parameterText

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_VarString_Array','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          parameterText=getTextContent(valueElement)
          nEntries=String_Count_Words(char(parameterText))
          if (nEntries > size(parameterValue)) then
             call Galacticus_Error_Report('Get_Input_Parameter_VarString_Array','array parameter has too many entries')
          else
             call String_Split_Words(parameterValue(1:nEntries),char(parameterText))
          end if
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_VarString_Array','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',parameterValue)
    end if
    return
  end subroutine Get_Input_Parameter_VarString_Array

  subroutine Get_Input_Parameter_Double(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt double precision} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in)           :: parameterName
    double precision, intent(out)          :: parameterValue
    double precision, intent(in), optional :: defaultValue
    logical,          intent(in), optional :: writeOutput
    type(Node),       pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                    :: datasetID
    integer                                :: iParameter
    logical                                :: foundMatch,writeOutputActual
    double precision                       :: datasetValue(1)

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Double','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          call extractDataContent(valueElement,parameterValue)
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Double','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       datasetValue=[parameterValue]
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',datasetValue)
    end if
    return
  end subroutine Get_Input_Parameter_Double

  subroutine Get_Input_Parameter_Double_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt double precision} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in)           :: parameterName
    double precision, intent(out)          :: parameterValue(:)
    double precision, intent(in), optional :: defaultValue(:)
    logical,          intent(in), optional :: writeOutput
    type(Node),       pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                    :: datasetID
    integer                                :: iParameter
    logical                                :: foundMatch,writeOutputActual

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Double','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          call extractDataContent(valueElement,parameterValue)
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Double_Array','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',parameterValue)
    end if
    return
  end subroutine Get_Input_Parameter_Double_Array

  subroutine Get_Input_Parameter_Integer(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt integer} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in)           :: parameterName
    integer,          intent(out)          :: parameterValue
    integer,          intent(in), optional :: defaultValue
    logical,          intent(in), optional :: writeOutput
    type(Node),       pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                    :: datasetID
    integer                                :: iParameter
    logical                                :: foundMatch,writeOutputActual
    integer                                :: datasetValue(1)

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Integer','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          call extractDataContent(valueElement,parameterValue)
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Integer','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       datasetValue=[parameterValue]
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',datasetValue)
    end if
    return
  end subroutine Get_Input_Parameter_Integer

  subroutine Get_Input_Parameter_Integer_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read an {\tt integer} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in)           :: parameterName
    integer,          intent(out)          :: parameterValue(:)
    integer,          intent(in), optional :: defaultValue(:)
    logical,          intent(in), optional :: writeOutput
    type(Node),       pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                    :: datasetID
    integer                                :: iParameter
    logical                                :: foundMatch,writeOutputActual

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Integer','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          call extractDataContent(valueElement,parameterValue)
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Integer_Array','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',parameterValue)
    end if
    return
  end subroutine Get_Input_Parameter_Integer_Array

  subroutine Get_Input_Parameter_Logical(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt logical} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in)           :: parameterName
    logical,          intent(out)          :: parameterValue
    logical,          intent(in), optional :: defaultValue
    logical,          intent(in), optional :: writeOutput
    type(Node),       pointer              :: thisParameter,nameElement,valueElement
    integer(kind=HID_T)                    :: datasetID
    integer                                :: iParameter
    logical                                :: foundMatch,writeOutputActual
    character(len=5)                       :: datasetValue(1)

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Logical','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          call extractDataContent(valueElement,parameterValue)
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Logical','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       select case (parameterValue)
       case (.false.)
          datasetValue=['false']
       case (.true.)
          datasetValue=['true']
       end select
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',datasetValue)
    end if
    return
  end subroutine Get_Input_Parameter_Logical

  subroutine Get_Input_Parameter_Logical_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read an {\tt logical} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in)           :: parameterName
    logical,          intent(out)          :: parameterValue(:)
    logical,          intent(in), optional :: defaultValue(:)
    logical,          intent(in), optional :: writeOutput
    type(Node),       pointer              :: thisParameter,nameElement,valueElement
    character(len=5)                       :: datasetValue(size(parameterValue))
    integer(kind=HID_T)                    :: datasetID
    integer                                :: iParameter
    logical                                :: foundMatch,writeOutputActual

    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Logical','parameter file has not been parsed.')
       end if
    end if

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    foundMatch=.false.
    do while (.not.foundMatch.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) then
          valueElement => item(getElementsByTagname(thisParameter,"value"),0)
          call extractDataContent(valueElement,parameterValue)
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Logical_Array','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=.true.
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       datasetID=0
       where (parameterValue)
          datasetValue='true'
       elsewhere
          datasetValue='false'
       end where
       call Galacticus_Output_Dataset(parametersGroupID,datasetID,parameterName,'',datasetValue)
    end if
    return
  end subroutine Get_Input_Parameter_Logical_Array

  subroutine Write_Parameter(parameterDoc,parameterName,parameterValue)
    !% Add a parameter to the specified XML file.
    use FoX_wxml
    implicit none
    type(xmlf_t),     intent(inout) :: parameterDoc
    character(len=*), intent(in)    :: parameterName,parameterValue

    call xml_NewElement(parameterDoc,"parameter")
    call xml_NewElement(parameterDoc,"name")
    call xml_AddCharacters(parameterDoc,trim(parameterName))
    call xml_EndElement(parameterDoc,"name")
    call xml_NewElement(parameterDoc,"value")
    call xml_AddCharacters(parameterDoc,trim(parameterValue))
    call xml_EndElement(parameterDoc,"value")
    call xml_EndElement(parameterDoc,"parameter")
    return
  end subroutine Write_Parameter

  subroutine Make_Parameters_Group
    !% Create a group in the \glc\ output file in which to store parameters.
    use ISO_Varying_String
    implicit none
    type(varying_string) :: groupName,commentText

    if (parametersGroupID == 0) then
       groupName='Parameters'
       commentText='Contains values for Galacticus parameters'
       parametersGroupID=Galacticus_Output_Make_Group(groupName,commentText)
    end if
    return
  end subroutine Make_Parameters_Group

end module Input_Parameters
