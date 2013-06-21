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

!% Contains a module which implements reading of parameters from an XML data file.

module Input_Parameters
  !% Implements reading of parameters from an XML data file.
  use, intrinsic :: ISO_C_Binding
  use FoX_dom
  use HDF5
  use IO_HDF5
  use ISO_Varying_String
  use Galacticus_Error
  use Hashes_Cryptographic
  use Galacticus_Versioning
  implicit none
  private
  public :: Input_Parameters_File_Open, Close_Parameters_Group, Input_Parameters_File_Close, Get_Input_Parameter, Get_Input_Parameter_Array_Size,&
       & Write_Parameter, Input_Parameter_Is_Present

  ! Include public specifiers for functions that will generate unique labels for modules.
  include 'utility.input_parameters.unique_labels.visibilities.inc'

  ! Node to hold the parameter document.
  type   (Node    ), pointer :: parameterDoc  =>null() 
  type   (NodeList), pointer :: parameterList          
  integer                    :: parameterCount         
  
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
     module procedure Get_Input_Parameter_Integer_Long
     module procedure Get_Input_Parameter_Integer_Long_Array
     module procedure Get_Input_Parameter_Logical
     module procedure Get_Input_Parameter_Logical_Array
  end interface

  ! Maximum length for parameters that must be extracted as text.
  integer            , parameter :: parameterLengthMaximum=1024    
  
  ! Parameters group identifier in the output file.
  logical                        :: parametersGroupCreated=.false. 
  type   (hdf5Object)            :: parametersGroup                
  
  ! Local pointer to the main output file object.
  logical                        :: haveOutputFile                 
  type   (hdf5Object), pointer   :: outputFileObject               
  
contains

  subroutine Input_Parameters_File_Open(parameterFile,outputFileObjectTarget,allowedParametersFile)
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
    use Galacticus_Input_Paths
    use String_Handling
    !$ use OMP_Lib
    implicit none
    type     (varying_string)         , intent(in   )                   :: parameterFile                                       
    type     (hdf5Object    )         , intent(in   ), optional, target :: outputFileObjectTarget                              
    character(len=*         )         , intent(in   ), optional         :: allowedParametersFile                               
    type     (Node          ), pointer                                  :: allowedParameterDoc   , nameElement             , & 
         &                                                                 thisParameter                                       
    type     (NodeList      ), pointer                                  :: allowedParameterList                                
    logical                                                             :: parameterMatched      , unknownParametersPresent    
    integer                                                             :: allowedParameterCount , distance                , & 
         &                                                                 iParameter            , ioErr                   , & 
         &                                                                 jParameter            , minimumDistance             
    type     (varying_string)                                           :: possibleMatch         , unknownParameter            
    
    ! Open and parse the data file.
    !$omp critical (FoX_DOM_Access)
    parameterDoc => parseFile(char(parameterFile),iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Input_Parameters_File_Open','Unable to find or parse parameter file')
    parameterList => getElementsByTagname(parameterDoc,"parameter")
    parameterCount=getLength(parameterList)
    !$omp end critical (FoX_DOM_Access)

    ! Create a pointer to the output object if given.
    if (present(outputFileObjectTarget)) then
       outputFileObject => outputFileObjectTarget
       haveOutputFile=.true.
    else
       haveOutputFile=.false.
    end if

    ! Open and parse the allowed parameters file if present.
    if (present(allowedParametersFile)) then
       !$omp critical (FoX_DOM_Access)
       ! Parse the file.
       allowedParameterDoc => parseFile(char(Galacticus_Input_Path())//'work/build/'//allowedParametersFile,iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Input_Parameters_File_Open','Unable to find or parse allowed parameters file')
       ! Extract the list of parameters.
       allowedParameterList => getElementsByTagname(allowedParameterDoc,"parameter")
       allowedParameterCount=getLength(allowedParameterList)
       ! Loop over all parameters in the input file.
       unknownParametersPresent=.false.
       do iParameter=0,parameterCount-1
          thisParameter => item(parameterList,iParameter)
          nameElement   => item(getElementsByTagname(thisParameter,"name"),0)
          jParameter    =  0
          parameterMatched=.false.
          do while (.not.parameterMatched .and. jParameter < allowedParameterCount)
             thisParameter   => item(allowedParameterList,jParameter)
             parameterMatched=(getTextContent(nameElement) == getTextContent(thisParameter))
             jParameter=jParameter+1
          end do
          if (.not.parameterMatched) then
             if (.not.unknownParametersPresent) then
                !$ if (omp_in_parallel()) then
                !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
                !$ else
                !$    write (0,'(a2,a2,$)') "MM",": "
                !$ end if
                write (0,'(a)') '-> WARNING: unknown parameters present:'
             end if
             unknownParametersPresent=.true.
             !$ if (omp_in_parallel()) then
             !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
             !$ else
             !$    write (0,'(a2,a2,$)') "MM",": "
             !$ end if
            unknownParameter=getTextContent(nameElement)
            minimumDistance=10000
            do jParameter=0,allowedParameterCount-1
               thisParameter => item(allowedParameterList,jParameter)
               distance=String_Levenshtein_Distance(char(unknownParameter),getTextContent(thisParameter))
               if (distance < minimumDistance) then
                  minimumDistance=distance
                  possibleMatch=getTextContent(thisParameter)
               end if
            end do
            write (0,'(5a)') '    ',char(unknownParameter),' [did you mean "',char(possibleMatch),'"?]'
          end if
       end do
       ! Destroy the document.
       call destroy(allowedParameterDoc)
       !$omp end critical (FoX_DOM_Access)
       if (unknownParametersPresent) then
          !$ if (omp_in_parallel()) then
          !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
          !$ else
          !$    write (0,'(a2,a2,$)') "MM",": "
          !$ end if
          write (0,'(a)') '<-'
       end if
    end if

    return
  end subroutine Input_Parameters_File_Open

  subroutine Input_Parameters_File_Close
    !% Close the parameter file (actually just destroy the internal record of it and clean up memory).
    implicit none

    call Close_Parameters_Group()
    call destroy(parameterDoc)
    return
  end subroutine Input_Parameters_File_Close

  logical function Input_Parameter_Is_Present(parameterName)
    !% Return true if {\tt parameterName} is present in the input file.
    implicit none
    character(len=*), intent(in   ) :: parameterName                
    type     (Node ), pointer       :: nameElement  , thisParameter 
    integer                         :: iParameter                   
    
    ! If no parameter file has been read, stop with an error message.
    if (.not.associated(parameterDoc)) call Galacticus_Error_Report('Input_Parameter_Is_Present','parameter file has not been parsed.')

    !$omp critical (FoX_DOM_Access)
    iParameter=0
    Input_Parameter_Is_Present=.false.
    do while (.not.Input_Parameter_Is_Present.and.iParameter<parameterCount)
       thisParameter => item(parameterList, iParameter)
       nameElement => item(getElementsByTagname(thisParameter,"name"),0)
       if (parameterName == getTextContent(nameElement)) Input_Parameter_Is_Present=.true.
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    return
  end function Input_Parameter_Is_Present

  integer function Get_Input_Parameter_Array_Size(parameterName)
    !% Get the number of elements in the parameter specified by parameter name is specified by {\tt parameterName}.
    use String_Handling
    implicit none
    character(len=*         ), intent(in   ) :: parameterName                              
    type     (Node          ), pointer       :: nameElement  , thisParameter, valueElement 
    integer                                  :: iParameter                                 
    logical                                  :: foundMatch                                 
    type     (varying_string)                :: parameterText                              
    
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
    character(len=*), intent(  out)           :: parameterValue                                  
    character(len=*), intent(in   )           :: parameterName                                   
    character(len=*), intent(in   ), optional :: defaultValue                                    
    logical         , intent(in   ), optional :: writeOutput                                     
    type     (Node ), pointer                 :: nameElement   , thisParameter    , valueElement 
    integer                                   :: iParameter                                      
    logical                                   :: foundMatch    , writeOutputActual               
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(trim(parameterValue),trim(parameterName))
       !$omp end critical(HDF5_Access)
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
    character(len=*         ), intent(  out)           :: parameterValue(:)                                  
    character(len=*         ), intent(in   )           :: parameterName                                      
    character(len=*         ), intent(in   ), optional :: defaultValue  (:)                                  
    logical                  , intent(in   ), optional :: writeOutput                                        
    type     (Node          ), pointer                 :: nameElement      , thisParameter    , valueElement 
    integer                                            :: iParameter       , nEntries                        
    logical                                            :: foundMatch       , writeOutputActual               
    type     (varying_string)                          :: parameterText                                      
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Get_Input_Parameter_Char_Array

  subroutine Get_Input_Parameter_VarString(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt varying\_string} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    type     (varying_string), intent(  out)           :: parameterValue                                  
    character(len=*         ), intent(in   )           :: parameterName                                   
    character(len=*         ), intent(in   ), optional :: defaultValue                                    
    logical                  , intent(in   ), optional :: writeOutput                                     
    type     (Node          ), pointer                 :: nameElement   , thisParameter    , valueElement 
    integer                                            :: iParameter                                      
    logical                                            :: foundMatch    , writeOutputActual               
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))    
       !$omp end critical(HDF5_Access)
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
    type     (varying_string), intent(  out)           :: parameterValue(:)                                  
    character(len=*         ), intent(in   )           :: parameterName                                      
    character(len=*         ), intent(in   ), optional :: defaultValue  (:)                                  
    logical                  , intent(in   ), optional :: writeOutput                                        
    type     (Node          ), pointer                 :: nameElement      , thisParameter    , valueElement 
    integer                                            :: iParameter       , nEntries                        
    logical                                            :: foundMatch       , writeOutputActual               
    type     (varying_string)                          :: parameterText                                      
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Get_Input_Parameter_VarString_Array

  subroutine Get_Input_Parameter_Double(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt double precision} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character       (len=*), intent(in   )           :: parameterName                                   
    double precision       , intent(  out)           :: parameterValue                                  
    double precision       , intent(in   ), optional :: defaultValue                                    
    logical                , intent(in   ), optional :: writeOutput                                     
    type            (Node ), pointer                 :: nameElement   , thisParameter    , valueElement 
    integer                                          :: iParameter                                      
    logical                                          :: foundMatch    , writeOutputActual               
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Get_Input_Parameter_Double

  subroutine Get_Input_Parameter_Double_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt double precision} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character       (len=*), intent(in   )           :: parameterName                                      
    double precision       , intent(  out)           :: parameterValue(:)                                  
    double precision       , intent(in   ), optional :: defaultValue  (:)                                  
    logical                , intent(in   ), optional :: writeOutput                                        
    type            (Node ), pointer                 :: nameElement      , thisParameter    , valueElement 
    integer                                          :: iParameter                                         
    logical                                          :: foundMatch       , writeOutputActual               
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Get_Input_Parameter_Double_Array

  subroutine Get_Input_Parameter_Integer(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt integer} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in   )           :: parameterName                                   
    integer         , intent(  out)           :: parameterValue                                  
    integer         , intent(in   ), optional :: defaultValue                                    
    logical         , intent(in   ), optional :: writeOutput                                     
    type     (Node ), pointer                 :: nameElement   , thisParameter    , valueElement 
    integer                                   :: iParameter                                      
    logical                                   :: foundMatch    , writeOutputActual               
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Get_Input_Parameter_Integer

  subroutine Get_Input_Parameter_Integer_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read an {\tt integer} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in   )           :: parameterName                                      
    integer         , intent(  out)           :: parameterValue(:)                                  
    integer         , intent(in   ), optional :: defaultValue  (:)                                  
    logical         , intent(in   ), optional :: writeOutput                                        
    type     (Node ), pointer                 :: nameElement      , thisParameter    , valueElement 
    integer                                   :: iParameter                                         
    logical                                   :: foundMatch       , writeOutputActual               
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Get_Input_Parameter_Integer_Array

  subroutine Get_Input_Parameter_Logical(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\tt logical} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in   )           :: parameterName                                   
    logical         , intent(  out)           :: parameterValue                                  
    logical         , intent(in   ), optional :: defaultValue                                    
    logical         , intent(in   ), optional :: writeOutput                                     
    type     (Node ), pointer                 :: nameElement   , thisParameter    , valueElement 
    integer                                   :: iParameter                                      
    logical                                   :: foundMatch    , writeOutputActual               
    character(len=5)                          :: datasetValue                                    
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       select case (parameterValue)
       case (.false.)
          datasetValue='false'
       case (.true.)
          datasetValue='true'
       end select
       call parametersGroup%writeAttribute(datasetValue,trim(parameterName))
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Get_Input_Parameter_Logical

  subroutine Get_Input_Parameter_Logical_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read an {\tt logical} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    implicit none
    character(len=*), intent(in   )           :: parameterName                                              
    logical         , intent(  out)           :: parameterValue(:                   )                       
    logical         , intent(in   ), optional :: defaultValue  (:                   )                       
    logical         , intent(in   ), optional :: writeOutput                                                
    type     (Node ), pointer                 :: nameElement                         , thisParameter    , & 
         &                                       valueElement                                               
    character(len=5)                          :: datasetValue  (size(parameterValue))                       
    integer                                   :: iParameter                                                 
    logical                                   :: foundMatch                          , writeOutputActual    
    
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
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       !$omp critical(HDF5_Access)
       call Make_Parameters_Group
       where (parameterValue)
          datasetValue='true'
       elsewhere
          datasetValue='false'
       end where
       call parametersGroup%writeAttribute(datasetValue,trim(parameterName))
       !$omp end critical(HDF5_Access)
    end if
    return
  end subroutine Get_Input_Parameter_Logical_Array

  subroutine Get_Input_Parameter_Integer_Long(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a long {\tt integer} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    use Kind_Numbers
    implicit none
    character(len=*                     ), intent(in   )           :: parameterName                                   
    integer  (kind=kind_int8            ), intent(  out)           :: parameterValue                                  
    integer  (kind=kind_int8            ), intent(in   ), optional :: defaultValue                                    
    logical                              , intent(in   ), optional :: writeOutput                                     
    type     (Node                      ), pointer                 :: nameElement   , thisParameter    , valueElement 
    integer                                                        :: iParameter                                      
    logical                                                        :: foundMatch    , writeOutputActual               
    character(len=parameterLengthMaximum)                          :: parameterText                                   
    
    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Integer_Long','parameter file has not been parsed.')
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
          read (parameterText,*) parameterValue
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Integer_Long','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))
    end if
    return
  end subroutine Get_Input_Parameter_Integer_Long

  subroutine Get_Input_Parameter_Integer_Long_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a long {\tt integer} parameter from the parameter file. The parameter name is specified by {\tt parameterName} and
    !% its value is returned in {\tt parameterValue}. If no parameter file has been opened by
    !% \hyperlink{utility.input_parameters.F90:input_parameters:input_parameters_file_open}{{\tt Input\_Parameters\_File\_Open}} or no matching parameter is found, the
    !%  default value (if any) given by {\tt defaultValue} is returned. (If no default value is present an error occurs instead.)
    use Kind_Numbers
    implicit none
    character(len=*                     ), intent(in   )           :: parameterName                                      
    integer  (kind=kind_int8            ), intent(  out)           :: parameterValue(:)                                  
    integer  (kind=kind_int8            ), intent(in   ), optional :: defaultValue  (:)                                  
    logical                              , intent(in   ), optional :: writeOutput                                        
    type     (Node                      ), pointer                 :: nameElement      , thisParameter    , valueElement 
    integer                                                        :: iParameter                                         
    logical                                                        :: foundMatch       , writeOutputActual               
    character(len=parameterLengthMaximum)                          :: parameterText                                      
    
    ! If no parameter file has been read, either return the default or stop with an error message.
    if (.not.associated(parameterDoc)) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Integer_Long_Array','parameter file has not been parsed.')
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
          read (parameterText,*) parameterValue
          foundMatch=.true.
       end if
       iParameter=iParameter+1
    end do
    !$omp end critical (FoX_DOM_Access)
    if (.not.foundMatch) then
       if (present(defaultValue)) then
          parameterValue=defaultValue
       else
          call Galacticus_Error_Report('Get_Input_Parameter_Integer_Long_Array','parameter '//trim(parameterName)//' can not be found')
       end if
    end if

    ! Write the parameter to the output file.
    if (present(writeOutput)) then
       writeOutputActual=writeOutput
    else
       writeOutputActual=haveOutputFile .and. outputFileObject%isOpen()
    end if
    if (writeOutputActual) then
       call Make_Parameters_Group
       call parametersGroup%writeAttribute(parameterValue,trim(parameterName))
    end if
    return
  end subroutine Get_Input_Parameter_Integer_Long_Array

  subroutine Write_Parameter(parameterDoc,parameterName,parameterValue)
    !% Add a parameter to the specified XML file.
    use FoX_wxml
    implicit none
    type     (xmlf_t), intent(inout) :: parameterDoc                  
    character(len=* ), intent(in   ) :: parameterName, parameterValue 
    
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
    implicit none

    if (.not.parametersGroupCreated) then
       parametersGroup=outputFileObject%openGroup('Parameters','Contains values for Galacticus parameters')
       parametersGroupCreated=.true.
    end if
    return
  end subroutine Make_Parameters_Group

  !# <hdfPreCloseTask>
  !# <unitName>Close_Parameters_Group</unitName>
  !# </hdfPreCloseTask>
  subroutine Close_Parameters_Group() 
    implicit none
    
    !$omp critical (HDF5_Access)
    if (parametersGroup%isOpen()) call parametersGroup%close()
    !$omp end critical (HDF5_Access)
    return
  end subroutine Close_Parameters_Group
  
  subroutine Get_Input_Parameter_Double_C(parameterNameLength,parameterName,parameterValue,defaultValue) bind(c,name="Get_Input_Parameter_Double")
    !% C-bound wrapper function for getting {\tt double precision} parameter values.
    use ISO_Varying_String
    use String_Handling
    implicit none
    integer  (kind=c_int    ), value :: parameterNameLength                      
    character(kind=c_char   )        :: parameterName      (parameterNameLength) 
    real     (kind=c_double )        :: parameterValue                           
    real     (kind=c_double )        :: defaultValue                             
    type     (varying_string)        :: parameterNameF                           
    
    parameterNameF=String_C_to_Fortran(parameterName)
    call Get_Input_Parameter_Double(char(parameterNameF),parameterValue,defaultValue)
    return
  end subroutine Get_Input_Parameter_Double_C

  subroutine Get_Input_Parameter_Integer_C(parameterNameLength,parameterName,parameterValue,defaultValue) bind(c,name="Get_Input_Parameter_Integer")
    !% C-bound wrapper function for getting {\tt integer} parameter values.
    use ISO_Varying_String
    use String_Handling
    implicit none
    integer  (kind=c_int    ), value :: parameterNameLength                      
    character(kind=c_char   )        :: parameterName      (parameterNameLength) 
    integer  (kind=c_int    )        :: parameterValue                           
    integer  (kind=c_int    )        :: defaultValue                             
    type     (varying_string)        :: parameterNameF                           
    
    parameterNameF=String_C_to_Fortran(parameterName)
    call Get_Input_Parameter_Integer(char(parameterNameF),parameterValue,defaultValue)
    return
  end subroutine Get_Input_Parameter_Integer_C

  ! Include functions that generate unique labels for modules.
  include 'utility.input_parameters.unique_labels.inc'

end module Input_Parameters
