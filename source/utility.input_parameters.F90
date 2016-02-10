!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use IO_XML
  use IO_HDF5
  use ISO_Varying_String
  use Galacticus_Error
  use Hashes_Cryptographic
  use Galacticus_Versioning
  use Galacticus_Build
  use Input_Parameters2
  implicit none
  private
  public :: Input_Parameters_File_Open, Input_Parameters_File_Close, Get_Input_Parameter, Get_Input_Parameter_Array_Size,&
       & Input_Parameter_Is_Present

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
  end interface Get_Input_Parameter

  type(inputParameters), target :: globalParameters1
  
contains

  subroutine Input_Parameters_File_Open(parameterFile,outputFileObject,allowedParametersFile)
    !% Open an XML data file containing parameter values and parse it.
    implicit none
    type     (varying_string)         , intent(in   )           :: parameterFile
    type     (hdf5Object    )         , intent(in   ), optional :: outputFileObject
    character(len=*         )         , intent(in   ), optional :: allowedParametersFile

    ! Wrap the new input parameters code.
    if (present(outputFileObject)) then
       globalParameters1=inputParameters(parameterFile,allowedParametersFile=allowedParametersFile,outputParametersGroup=outputFileObject)
    else
       globalParameters1=inputParameters(parameterFile,allowedParametersFile=allowedParametersFile                                       )
    end if
    call globalParameters1%markGlobal()
    globalParameters => globalParameters1
    return
  end subroutine Input_Parameters_File_Open

  !# <hdfPreCloseTask>
  !# <unitName>Input_Parameters_File_Close</unitName>
  !# </hdfPreCloseTask>
  subroutine Input_Parameters_File_Close()
    !% Close the parameter file.
    implicit none

    call globalParameters1%destroy()
    return
  end subroutine Input_Parameters_File_Close

  logical function Input_Parameter_Is_Present(parameterName)
    !% Return true if {\normalfont \ttfamily parameterName} is present in the input file.
    implicit none
    character(len=*), intent(in   ) :: parameterName

    Input_Parameter_Is_Present=globalParameters1%isPresent(parameterName)
    return
  end function Input_Parameter_Is_Present

  integer function Get_Input_Parameter_Array_Size(parameterName)
    !% Get the number of elements in the parameter specified by parameter name is specified by {\normalfont \ttfamily parameterName}.
    use String_Handling
    implicit none
    character(len=*         ), intent(in   ) :: parameterName

    Get_Input_Parameter_Array_Size=globalParameters1%count(parameterName,zeroIfNotPresent=.true.)
    return
  end function Get_Input_Parameter_Array_Size

  subroutine Get_Input_Parameter_Char(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\normalfont \ttfamily varying\_string} parameter from the parameter file.
    implicit none
    character(len=*       ), intent(  out)           :: parameterValue
    character(len=*       ), intent(in   )           :: parameterName
    character(len=*       ), intent(in   ), optional :: defaultValue
    logical                , intent(in   ), optional :: writeOutput
 
    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Char

  subroutine Get_Input_Parameter_Char_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\normalfont \ttfamily varying\_string} parameter from the parameter file.
    use String_Handling
    implicit none
    character(len=*         ), intent(  out)           :: parameterValue(:)
    character(len=*         ), intent(in   )           :: parameterName
    character(len=*         ), intent(in   ), optional :: defaultValue  (:)
    logical                  , intent(in   ), optional :: writeOutput

    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Char_Array

  subroutine Get_Input_Parameter_VarString(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\normalfont \ttfamily varying\_string} parameter from the parameter file.
    implicit none
    type     (varying_string), intent(  out)           :: parameterValue
    character(len=*         ), intent(in   )           :: parameterName
    character(len=*         ), intent(in   ), optional :: defaultValue
    logical                  , intent(in   ), optional :: writeOutput

    if (present(defaultValue)) then
       call globalParameters1%value(parameterName,parameterValue,defaultValue=var_str(defaultValue),writeOutput=writeOutput)    
    else
       call globalParameters1%value(parameterName,parameterValue                                   ,writeOutput=writeOutput)    
    end if
    return
  end subroutine Get_Input_Parameter_VarString

  subroutine Get_Input_Parameter_VarString_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\normalfont \ttfamily varying\_string} parameter from the parameter file.
    implicit none
    type     (varying_string), intent(  out)           :: parameterValue(:)
    character(len=*         ), intent(in   )           :: parameterName
    character(len=*         ), intent(in   ), optional :: defaultValue  (:)
    logical                  , intent(in   ), optional :: writeOutput

    if (present(defaultValue)) then
       call globalParameters1%value(parameterName,parameterValue,defaultValue=var_str(defaultValue),writeOutput=writeOutput)
    else
       call globalParameters1%value(parameterName,parameterValue                                   ,writeOutput=writeOutput)
    end if
    return
  end subroutine Get_Input_Parameter_VarString_Array
  
  subroutine Get_Input_Parameter_Double(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\normalfont \ttfamily double precision} parameter from the parameter file.
    implicit none
    character       (len=*       ), intent(in   )           :: parameterName
    double precision              , intent(  out)           :: parameterValue
    double precision              , intent(in   ), optional :: defaultValue
    logical                       , intent(in   ), optional :: writeOutput

    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Double

  subroutine Get_Input_Parameter_Double_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\normalfont \ttfamily double precision} parameter from the parameter file.
    implicit none
    character       (len=*       ), intent(in   )           :: parameterName
    double precision              , intent(  out)           :: parameterValue(:)
    double precision              , intent(in   ), optional :: defaultValue  (:)
    logical                       , intent(in   ), optional :: writeOutput
    
    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Double_Array

  subroutine Get_Input_Parameter_Integer(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\normalfont \ttfamily integer} parameter from the parameter file.
    implicit none
    character(len=*       ), intent(in   )           :: parameterName
    integer                , intent(  out)           :: parameterValue
    integer                , intent(in   ), optional :: defaultValue
    logical                , intent(in   ), optional :: writeOutput

    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Integer

  subroutine Get_Input_Parameter_Integer_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read an {\normalfont \ttfamily integer} parameter from the parameter file.
    implicit none
    character(len=*       ), intent(in   )           :: parameterName
    integer                , intent(  out)           :: parameterValue(:)
    integer                , intent(in   ), optional :: defaultValue  (:)
    logical                , intent(in   ), optional :: writeOutput

    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Integer_Array

  subroutine Get_Input_Parameter_Logical(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a {\normalfont \ttfamily logical} parameter from the parameter file.
    implicit none
    character(len=*       ), intent(in   )           :: parameterName
    logical                , intent(  out)           :: parameterValue
    logical                , intent(in   ), optional :: defaultValue
    logical                , intent(in   ), optional :: writeOutput

    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Logical

  subroutine Get_Input_Parameter_Logical_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read an {\normalfont \ttfamily logical} parameter from the parameter file.
    implicit none
    character(len=*       ), intent(in   )           :: parameterName
    logical                , intent(  out)           :: parameterValue(:)
    logical                , intent(in   ), optional :: defaultValue  (:)
    logical                , intent(in   ), optional :: writeOutput

    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Logical_Array

  subroutine Get_Input_Parameter_Integer_Long(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a long {\normalfont \ttfamily integer} parameter from the parameter file.
    use Kind_Numbers
    implicit none
    character(len=*                     ), intent(in   )           :: parameterName
    integer  (kind=kind_int8            ), intent(  out)           :: parameterValue
    integer  (kind=kind_int8            ), intent(in   ), optional :: defaultValue
    logical                              , intent(in   ), optional :: writeOutput

    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Integer_Long

  subroutine Get_Input_Parameter_Integer_Long_Array(parameterName,parameterValue,defaultValue,writeOutput)
    !% Read a long {\normalfont \ttfamily integer} parameter from the parameter file.
    use Kind_Numbers
    implicit none
    character(len=*                     ), intent(in   )           :: parameterName
    integer  (kind=kind_int8            ), intent(  out)           :: parameterValue(:)
    integer  (kind=kind_int8            ), intent(in   ), optional :: defaultValue  (:)
    logical                              , intent(in   ), optional :: writeOutput

    call globalParameters1%value(parameterName,parameterValue,defaultValue=defaultValue,writeOutput=writeOutput)    
    return
  end subroutine Get_Input_Parameter_Integer_Long_Array

  subroutine Get_Input_Parameter_Double_C(parameterNameLength,parameterName,parameterValue,defaultValue) bind(c,name="Get_Input_Parameter_Double")
    !% C-bound wrapper function for getting {\normalfont \ttfamily double precision} parameter values.
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
    !% C-bound wrapper function for getting {\normalfont \ttfamily integer} parameter values.
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

end module Input_Parameters
