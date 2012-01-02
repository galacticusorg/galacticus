!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements various useful functionality for manipulating character strings.

module String_Handling
  !% Implements various useful functionality for manipulating character strings.
  use ISO_Varying_String
  implicit none
  private
  public :: operator(//), String_Split_Words, String_Count_Words, String_Upper_Case, String_Lower_Case, String_Upper_Case_First,&
       & Convert_VarString_To_Char, String_C_to_Fortran, String_Subscript, String_Superscript

  interface operator(//)
     module procedure Concatenate_VarStr_Integer
     module procedure Concatenate_VarStr_Integer8
  end interface

  interface String_Split_Words
     module procedure String_Split_Words_VarString
     module procedure String_Split_Words_Char
  end interface

  ! Maximum length of string needed to hold integer values.
  integer,          parameter :: maxIntegerSize  =   20
  character(len=5), parameter :: maxIntegerFormat='(i20)'

  ! Character strings used in converting upper to lower case and vice-versa.
  character(len=*), parameter :: charactersLowerCase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter :: charactersUpperCase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  ! Character strings used in whitespace detection.
  character(len=*), parameter :: charactersWhiteSpace=' '//char(0)//char(9)//char(10)//char(13)

  ! Character strings used in converting to subscripts and superscripts.
  character(len=*), parameter :: charactersScript     ='0123456789+-=()'
  character(len=*), parameter :: charactersSubscript  ='₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎'
  character(len=*), parameter :: charactersSuperscript='⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾'

contains

  integer function String_Count_Words(inputString,separator)
    !% Return a count of the number of space separated words in {\tt inputString}.
    implicit none
    character(len=*), intent(in)           :: inputString
    character(len=*), intent(in), optional :: separator
    logical                                :: inWord
    integer                                :: iCharacter
    type(varying_string)                   :: separatorActual

    ! Decide what separator to use.
    if (present(separator)) then
       separatorActual=separator
    else
       separatorActual=charactersWhiteSpace
    end if

    String_Count_Words=0
    inWord=.false.
    do iCharacter=1,len_trim(inputString)
       if (index(separatorActual,inputString(iCharacter:iCharacter)) /= 0) then
          if (inWord) String_Count_Words=String_Count_Words+1
          inWord=.false.
       else
          inWord=.true.
       end if
    end do
    if (inWord) String_Count_Words=String_Count_Words+1
    return
  end function String_Count_Words

  subroutine String_Split_Words_VarString(words,inputString,separator)
    !% Split {\tt inputString} into words and return as an array.
    implicit none
    type(varying_string), intent(out), dimension(:) :: words
    character(len=*),     intent(in)                :: inputString
    character(len=*),     intent(in),  optional     :: separator
    logical                                         :: inWord
    integer                                         :: iCharacterStart,iCharacter,iWord
    type(varying_string)                            :: separatorActual

    ! Decide what separator to use.
    if (present(separator)) then
       separatorActual=separator
    else
       separatorActual=charactersWhiteSpace
    end if

    inWord=.false.
    iWord=0
    do iCharacter=1,len_trim(inputString)
       if (index(separatorActual,inputString(iCharacter:iCharacter)) /= 0) then
          if (inWord) then
             iWord=iWord+1
             words(iWord)=inputString(iCharacterStart:iCharacter-1)
             if (iWord == size(words)) return
          end if
          inWord=.false.
       else
          if (.not.inWord) iCharacterStart=iCharacter
          inWord=.true.
       end if
    end do
    if (inWord) then
       iWord=iWord+1
       words(iWord)=inputString(iCharacterStart:len_trim(inputString))
    end if
    return
  end subroutine String_Split_Words_VarString

  subroutine String_Split_Words_Char(words,inputString,separator)
    !% Split {\tt inputString} into words and return as an array.
    implicit none
    character(len=*), intent(out), dimension(:) :: words
    character(len=*), intent(in)                :: inputString
    character(len=*), intent(in),  optional     :: separator
    logical                                     :: inWord
    integer                                     :: iCharacterStart,iCharacter,iWord
    type(varying_string)                        :: separatorActual

    ! Decide what separator to use.
    if (present(separator)) then
       separatorActual=separator
    else
       separatorActual=charactersWhiteSpace
    end if

    inWord=.false.
    iWord=0
    do iCharacter=1,len_trim(inputString)
       if (index(separatorActual,inputString(iCharacter:iCharacter)) /= 0) then
          if (inWord) then
             iWord=iWord+1
             words(iWord)=inputString(iCharacterStart:iCharacter-1)
             if (iWord == size(words)) return
          end if
          inWord=.false.
       else
          if (.not.inWord) iCharacterStart=iCharacter
          inWord=.true.
       end if
    end do
    if (inWord) then
       iWord=iWord+1
       words(iWord)=inputString(iCharacterStart:len_trim(inputString))
    end if
    return
  end subroutine String_Split_Words_Char

  function Concatenate_VarStr_Integer(varStrVariable,intVariable)
    !% Provides a concatenation operator to append an integer number to a {\tt varying\_string}.
    implicit none
    type(varying_string),         intent(in) :: varStrVariable
    integer,                      intent(in) :: intVariable
    type(varying_string)                     :: Concatenate_VarStr_Integer
    character(len=maxIntegerSize)            :: intString

    write (intString,maxIntegerFormat) intVariable
    Concatenate_VarStr_Integer=varStrVariable//trim(adjustl(intString))
    return
  end function Concatenate_VarStr_Integer

  function Concatenate_VarStr_Integer8(varStrVariable,intVariable)
    !% Provides a concatenation operator to append an integer number to a {\tt varying\_string}.
    use Kind_Numbers
    implicit none
    type(varying_string),         intent(in) :: varStrVariable
    integer(kind=kind_int8),      intent(in) :: intVariable
    type(varying_string)                     :: Concatenate_VarStr_Integer8
    character(len=maxIntegerSize)            :: intString

    write (intString,maxIntegerFormat) intVariable
    Concatenate_VarStr_Integer8=varStrVariable//trim(adjustl(intString))
    return
  end function Concatenate_VarStr_Integer8

  function String_Upper_Case(stringInput) result (stringOutput)
    !% Converts an input string to upper case.
    character(len=*),           intent(in) :: stringInput
    character(len(stringInput))            :: stringOutput
    integer                                :: iString,iCharacter

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Loop through each character in the string.
    do iString=1,len(stringOutput)
       ! Find position of current character in string in list of lower case characters.
       iCharacter=index(charactersLowerCase,stringOutput(iString:iString))
       ! If a match is found, repace with the upper case equivalent.
       if (iCharacter /= 0) stringOutput(iString:iString)=charactersUpperCase(iCharacter:iCharacter)
    end do
    return
  end function String_Upper_Case

  function String_Lower_Case(stringInput) result (stringOutput)
    !% Converts an input string to lower case.
    character(len=*),           intent(in) :: stringInput
    character(len(stringInput))            :: stringOutput
    integer                                :: iString,iCharacter

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Loop through each character in the string.
    do iString=1,len(stringOutput)
       ! Find position of current character in string in list of upper case characters.
       iCharacter=index(charactersUpperCase,stringOutput(iString:iString))
       ! If a match is found, repace with the lower case equivalent.
       if (iCharacter /= 0) stringOutput(iString:iString)=charactersLowerCase(iCharacter:iCharacter)
    end do
    return
  end function String_Lower_Case

  function String_Upper_Case_First(stringInput) result (stringOutput)
    !% Converts an input string to upper case.
    character(len=*),           intent(in) :: stringInput
    character(len(stringInput))            :: stringOutput
    integer                                :: iCharacter

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Find position of first character in string in list of lower case characters.
    iCharacter=index(charactersLowerCase,stringOutput(1:1))
    ! If a match is found, repace with the upper case equivalent.
    if (iCharacter /= 0) stringOutput(1:1)=charactersUpperCase(iCharacter:iCharacter)
    return
  end function String_Upper_Case_First

  function Convert_VarString_To_Char(varStrings)
    !% Convert an array of varying strings into an array of characters.
    implicit none
    type(varying_string),                   intent(in), dimension(:)                :: varStrings
    character(len=maxval(len(varStrings))),             dimension(size(varStrings)) :: Convert_VarString_To_Char
    integer                                                                         :: iString

    do iString=1,size(varStrings)
       Convert_VarString_To_Char(iString)=varStrings(iString)
    end do
    return
  end function Convert_VarString_To_Char

  function String_C_to_Fortran(charArray)
    !% Convert a C-style character array into a Fortran varying string variable.
    use, intrinsic :: ISO_C_Binding
    implicit none
    type(varying_string)          :: String_C_to_Fortran
    character(c_char), intent(in) :: charArray(:)
    integer                       :: i

    String_C_to_Fortran=""
    do i=1,size(charArray)
       String_C_to_Fortran=String_C_to_Fortran//charArray(i)
    end do
    return
  end function String_C_to_Fortran

  function String_Subscript(stringInput) result (stringOutput)
    !% Converts an input string to Unicode subscripts.
    character(len=*),           intent(in) :: stringInput
    character(len(stringInput))            :: stringOutput
    integer                                :: iString,iCharacter

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Loop through each character in the string.
    do iString=1,len(stringOutput)
       ! Find position of current character in string in list of upper case characters.
       iCharacter=index(charactersScript,stringOutput(iString:iString))
       ! If a match is found, repace with the lower case equivalent.
       if (iCharacter /= 0) stringOutput(iString:iString)=charactersSubscript(iCharacter:iCharacter)
    end do
    return
  end function String_Subscript

  function String_Superscript(stringInput) result (stringOutput)
    !% Converts an input string to Unicode superscripts.
    character(len=*),           intent(in) :: stringInput
    character(len(stringInput))            :: stringOutput
    integer                                :: iString,iCharacter

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Loop through each character in the string.
    do iString=1,len(stringOutput)
       ! Find position of current character in string in list of upper case characters.
       iCharacter=index(charactersScript,stringOutput(iString:iString))
       ! If a match is found, repace with the lower case equivalent.
       if (iCharacter /= 0) stringOutput(iString:iString)=charactersSuperscript(iCharacter:iCharacter)
    end do
    return
  end function String_Superscript

end module String_Handling
