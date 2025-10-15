!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!!{
Contains a module which implements various useful functionality for manipulating character strings.
!!}

module String_Handling
  !!{
  Implements various useful functionality for manipulating character strings.
  !!}
  implicit none
  private
  public :: operator(//)              , char                        , String_Split_Words                 , String_Count_Words         , &
       &    String_Upper_Case         , String_Lower_Case           , String_Upper_Case_First            , Convert_VarString_To_Char  , &
       &    String_C_to_Fortran       , String_Subscript            , String_Superscript                 , String_Levenshtein_Distance, &
       &    String_Join               , String_Strip                , String_Lower_Case_First            , String_Value_Type          , &
       &    String_Value_Extract_Float, String_Value_Extract_Integer, String_Value_Extract_Integer_Size_T, stringSubstitute           , &
       &    stringXMLFormat

  interface operator(//)
     module procedure Concatenate_VarStr_Integer
     module procedure Concatenate_VarStr_Integer8
  end interface operator(//)

  interface String_Split_Words
     module procedure String_Split_Words_VarString
     module procedure String_Split_Words_Char
  end interface String_Split_Words

  interface char
     module procedure Char_Logical
  end interface char

  interface stringXMLFormat
     module procedure stringXMLFormatVarStr
     module procedure stringXMLFormatChar
  end interface stringXMLFormat

  ! Maximum length of string needed to hold integer values.
  integer         , parameter :: maxIntegerSize=20
  character(len=5), parameter :: maxIntegerFormat='(i20)'

  ! Character strings used in converting upper to lower case and vice-versa.
  character(len=*), parameter :: charactersLowerCase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter :: charactersUpperCase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  ! Character strings used in whitespace detection.
  character(len=*), parameter :: charactersWhiteSpace=' '//char(0)//char(9)//char(10)//char(13)
  
  !![
  <enumeration>
   <name>valueType</name>
   <description>Enumerates possible value types for strings.</description>
   <visibility>public</visibility>
   <entry label="floating"/>
   <entry label="integer" />
   <entry label="other"   />
  </enumeration>
  !!]

contains

  integer function String_Count_Words(inputString,separator,bracketing)
    !!{
    Return a count of the number of space separated words in {\normalfont \ttfamily inputString}.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=), index
    implicit none
    character(len=*         ), intent(in   )           :: inputString
    character(len=*         ), intent(in   ), optional :: separator
    character(len=2         ), intent(in   ), optional :: bracketing
    logical                                            :: inWord
    integer                                            :: iCharacter     , inBracket
    type     (varying_string)                          :: separatorActual

    ! Decide what separator to use.
    if (present(separator)) then
       separatorActual=separator
    else
       separatorActual=charactersWhiteSpace
    end if

    String_Count_Words=0
    inWord            =.false.
    inBracket         = 0
    do iCharacter=1,len_trim(inputString)
       if (present(bracketing)) then
          if (inputString(iCharacter:iCharacter) == bracketing(1:1)) inBracket=inBracket+1
          if (inputString(iCharacter:iCharacter) == bracketing(2:2)) inBracket=inBracket-1
       end if
       if (index(separatorActual,inputString(iCharacter:iCharacter)) /= 0) then
          if (inBracket == 0) then
             if (inWord) String_Count_Words=String_Count_Words+1
             inWord=.false.
          end if
       else
          inWord=.true.
       end if
    end do
    if (inWord) String_Count_Words=String_Count_Words+1
    return
  end function String_Count_Words

  subroutine String_Split_Words_VarString(words,inputString,separator,bracketing)
    !!{
    Split {\normalfont \ttfamily inputString} into words and return as an array.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=), index
    implicit none
    type     (varying_string), dimension(:), intent(  out)           :: words
    character(len=*         )              , intent(in   )           :: inputString
    character(len=*         )              , intent(in   ), optional :: separator
    character(len=2         )              , intent(in   ), optional :: bracketing
    logical                                                          :: inWord
    integer                                                          :: iCharacter     , iCharacterStart, iWord, inBracket
    type     (varying_string)                                        :: separatorActual

    ! Decide what separator to use.
    if (present(separator)) then
       separatorActual=separator
    else
       separatorActual=charactersWhiteSpace
    end if

    words          =""
    inWord         =.false.
    iWord          = 0
    inBracket      = 0
    iCharacterStart=-1
    do iCharacter=1,len_trim(inputString)
       if (present(bracketing)) then
          if (inputString(iCharacter:iCharacter) == bracketing(1:1)) inBracket=inBracket+1
          if (inputString(iCharacter:iCharacter) == bracketing(2:2)) inBracket=inBracket-1
       end if
       if (index(separatorActual,inputString(iCharacter:iCharacter)) /= 0) then
          if (inBracket == 0) then
             if (inWord) then
                iWord=iWord+1
                words(iWord)=inputString(iCharacterStart:iCharacter-1)
                if (iWord == size(words)) return
             end if
             inWord=.false.
          end if
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

  subroutine String_Split_Words_Char(words,inputString,separator,bracketing)
    !!{
    Split {\normalfont \ttfamily inputString} into words and return as an array.
    !!}
    use :: ISO_Varying_String, only : varying_string, index, assignment(=)
    implicit none
    character(len=*         ), dimension(:), intent(  out)           :: words
    character(len=*         )              , intent(in   )           :: inputString
    character(len=*         )              , intent(in   ), optional :: separator
    character(len=2         )              , intent(in   ), optional :: bracketing
    logical                                                          :: inWord
    integer                                                          :: iCharacter     , iCharacterStart, iWord, inBracket
    type     (varying_string)                                        :: separatorActual

    ! Decide what separator to use.
    if (present(separator)) then
       separatorActual=separator
    else
       separatorActual=charactersWhiteSpace
    end if

    words          =""
    inWord         =.false.
    iWord          = 0
    inBracket      = 0
    iCharacterStart=-1
    do iCharacter=1,len_trim(inputString)
       if (present(bracketing)) then
          if (inputString(iCharacter:iCharacter) == bracketing(1:1)) inBracket=inBracket+1
          if (inputString(iCharacter:iCharacter) == bracketing(2:2)) inBracket=inBracket-1
       end if
       if (index(separatorActual,inputString(iCharacter:iCharacter)) /= 0) then
          if (inBracket == 0) then
             if (inWord) then
                iWord=iWord+1
                words(iWord)=inputString(iCharacterStart:iCharacter-1)
                if (iWord == size(words)) return
             end if
             inWord=.false.
          end if
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
    !!{
    Provides a concatenation operator to append an integer number to a {\normalfont \ttfamily varying\_string}.
    !!}
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    type     (varying_string    ), intent(in   ) :: varStrVariable
    integer                      , intent(in   ) :: intVariable
    type     (varying_string    )                :: Concatenate_VarStr_Integer
    character(len=maxIntegerSize)                :: intString

    write (intString,maxIntegerFormat) intVariable
    Concatenate_VarStr_Integer=varStrVariable//trim(adjustl(intString))
    return
  end function Concatenate_VarStr_Integer

  function Concatenate_VarStr_Integer8(varStrVariable,intVariable)
    !!{
    Provides a concatenation operator to append an integer number to a {\normalfont \ttfamily varying\_string}.
    !!}
    use :: Kind_Numbers      , only : kind_int8
    use :: ISO_Varying_String, only : varying_string, operator(//)
    implicit none
    type     (varying_string    ), intent(in   ) :: varStrVariable
    integer  (kind=kind_int8    ), intent(in   ) :: intVariable
    type     (varying_string    )                :: Concatenate_VarStr_Integer8
    character(len=maxIntegerSize)                :: intString

    write (intString,maxIntegerFormat) intVariable
    Concatenate_VarStr_Integer8=varStrVariable//trim(adjustl(intString))
    return
  end function Concatenate_VarStr_Integer8

  function String_Upper_Case(stringInput) result (stringOutput)
    !!{
    Converts an input string to upper case.
    !!}
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter  , iString

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Loop through each character in the string.
    do iString=1,len(stringOutput)
       ! Find position of current character in string in list of lower case characters.
       iCharacter=index(charactersLowerCase,stringOutput(iString:iString))
       ! If a match is found, replace with the upper case equivalent.
       if (iCharacter /= 0) stringOutput(iString:iString)=charactersUpperCase(iCharacter:iCharacter)
    end do
    return
  end function String_Upper_Case

  elemental function String_Lower_Case(stringInput) result (stringOutput)
    !!{
    Converts an input string to lower case.
    !!}
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter  , iString

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Loop through each character in the string.
    do iString=1,len(stringOutput)
       ! Find position of current character in string in list of upper case characters.
       iCharacter=index(charactersUpperCase,stringOutput(iString:iString))
       ! If a match is found, replace with the lower case equivalent.
       if (iCharacter /= 0) stringOutput(iString:iString)=charactersLowerCase(iCharacter:iCharacter)
    end do
    return
  end function String_Lower_Case

  function String_Upper_Case_First(stringInput) result (stringOutput)
    !!{
    Converts an input string to upper case.
    !!}
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Find position of first character in string in list of lower case characters.
    iCharacter=index(charactersLowerCase,stringOutput(1:1))
    ! If a match is found, replace with the upper case equivalent.
    if (iCharacter /= 0) stringOutput(1:1)=charactersUpperCase(iCharacter:iCharacter)
    return
  end function String_Upper_Case_First

  function String_Lower_Case_First(stringInput) result (stringOutput)
    !!{
    Converts the first character of an input string to lower case.
    !!}
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter

    ! Transfer input string to output string.
    stringOutput=stringInput
    ! Find position of first character in string in list of upper case characters.
    iCharacter=index(charactersUpperCase,stringOutput(1:1))
    ! If a match is found, replace with the lower case equivalent.
    if (iCharacter /= 0) stringOutput(1:1)=charactersLowerCase(iCharacter:iCharacter)
    return
  end function String_Lower_Case_First

  function Convert_VarString_To_Char(varStrings)
    !!{
    Convert an array of varying strings into an array of characters.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=), len
    implicit none
    type     (varying_string             ), dimension(:)               , intent(in   ) :: varStrings
    character(len=maxval(len(varStrings))), dimension(size(varStrings))                :: Convert_VarString_To_Char
    integer                                                                            :: iString

    do iString=1,size(varStrings)
       Convert_VarString_To_Char(iString)=varStrings(iString)
    end do
    return
  end function Convert_VarString_To_Char

  function String_C_to_Fortran(charArray)
    !!{
    Convert a C-style character array into a Fortran varying string variable.
    !!}
    use, intrinsic :: ISO_C_Binding     , only : c_char, c_null_char
    use            :: ISO_Varying_String, only : varying_string, assignment(=), operator(//)
    implicit none
    type     (varying_string)                              :: String_C_to_Fortran
    character(kind=c_char   ), intent(in   ), dimension(*) :: charArray
    integer                                                :: i

    String_C_to_Fortran=""
    i=1
    do while (charArray(i) /= c_null_char)
       String_C_to_Fortran=String_C_to_Fortran//charArray(i)
       i=i+1
    end do
    return
  end function String_C_to_Fortran

  function String_Subscript(stringInput) result (stringOutput)
    !!{
    Converts an input string to Unicode subscripts.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=), operator(//)
    character(len=*         ), intent(in   ) :: stringInput
    type     (varying_string)                :: stringOutput
    integer                                  :: iString
    
    stringOutput=""
    do iString=1,len(stringInput)      
       select case (stringInput(iString:iString))
       case ("0")
          stringOutput=stringOutput//"₀"
       case ("1")
          stringOutput=stringOutput//"₁"
       case ("2")
          stringOutput=stringOutput//"₂"
       case ("3")
          stringOutput=stringOutput//"₃"
       case ("4")
          stringOutput=stringOutput//"₄"
       case ("5")
          stringOutput=stringOutput//"₅"
       case ("6")
          stringOutput=stringOutput//"₆"
       case ("7")
          stringOutput=stringOutput//"₇"
       case ("8")
          stringOutput=stringOutput//"₈"
       case ("9")
          stringOutput=stringOutput//"₉"
       case ("+")
          stringOutput=stringOutput//"₊"
       case ("-")
          stringOutput=stringOutput//"₋"
       case ("=")
          stringOutput=stringOutput//"₌"
       case ("(")
          stringOutput=stringOutput//"₍"
       case (")")
          stringOutput=stringOutput//"₎"
       case default
          stringOutput=stringOutput//stringInput(iString:iString)
       end select
    end do
    return
  end function String_Subscript

  function String_Superscript(stringInput) result (stringOutput)
    !!{
    Converts an input string to Unicode superscripts.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=), operator(//)
    character(len=*         ), intent(in   ) :: stringInput
    type     (varying_string)                :: stringOutput
    integer                                  :: iString
    
    stringOutput=""
    do iString=1,len(stringInput)      
       select case (stringInput(iString:iString))
       case ("0")
          stringOutput=stringOutput//"⁰"
       case ("1")
          stringOutput=stringOutput//"¹"
       case ("2")
          stringOutput=stringOutput//"²"
       case ("3")
          stringOutput=stringOutput//"³"
       case ("4")
          stringOutput=stringOutput//"⁴"
       case ("5")
          stringOutput=stringOutput//"⁵"
       case ("6")
          stringOutput=stringOutput//"⁶"
       case ("7")
          stringOutput=stringOutput//"⁷"
       case ("8")
          stringOutput=stringOutput//"⁸"
       case ("9")
          stringOutput=stringOutput//"⁹"
       case ("+")
          stringOutput=stringOutput//"⁺"
       case ("-")
          stringOutput=stringOutput//"⁻"
       case ("=")
          stringOutput=stringOutput//"⁼"
       case ("(")
          stringOutput=stringOutput//"⁽"
       case (")")
          stringOutput=stringOutput//"⁾"
       case default
          stringOutput=stringOutput//stringInput(iString:iString)
       end select
    end do
    return
  end function String_Superscript

  integer function String_Levenshtein_Distance(s,t)
    !!{
    Compute the \href{http://en.wikipedia.org/wiki/Levenshtein_distance}{Levenshtein distance} between strings {\normalfont \ttfamily a} and {\normalfont \ttfamily
    b}.
    !!}
    implicit none
    character(len=*), intent(in   )                :: s, t
    integer         , dimension(0:len(s),0:len(t)) :: d
    integer                                        :: i, j, m, n

    m=len(s)
    n=len(t)
    do i=0,m
       d(i,0)=i ! The distance of any first string to an empty second string.
    end do
    do j=0,n
       d(0,j)=j ! The distance of any second string to an empty first string.
    end do
    do j=1,n
       do i=1,m
          if (s(i:i) == t(j:j)) then
             d(i,j)=d(i-1,j-1)       ! No operation required.
          else
             d(i,j)=minval(               &
                  &        [              &
                  &         d(i-1,j  )+1, & ! A deletion.
                  &         d(i  ,j-1)+1, & ! An insertion.
                  &         d(i-1,j-1)+1  & ! A substitution.
                  &        ]              &
                  &       )
                end if
             end do
          end do
          String_Levenshtein_Distance=d(m,n)
    return
  end function String_Levenshtein_Distance
  
  function String_Join(strings,separator)
    !!{
    Joins an array of strings into one long string with the given separator.
    !!}
    use :: ISO_Varying_String, only : varying_string, operator(//), assignment(=)
    implicit none
    type     (varying_string)                              :: String_Join
    type     (varying_string), dimension(:), intent(in   ) :: strings
    character(len=*         )              , intent(in   ) :: separator
    integer                                                :: i

    String_Join=""
    do i=1,size(strings)
       String_Join=String_Join//strings(i)
       if (i < size(strings)) String_Join=String_Join//separator
    end do
    return
  end function String_Join

  function String_Strip(string)
    !!{
    Strips a string of leading and trailing whitespace, including tabs.
    !!}
    use :: ISO_Varying_String, only : varying_string, assignment(=), index
    implicit none
    type     (varying_string)                :: String_Strip
    character(len=*         ), intent(in   ) :: string
    integer                                  :: iBegin      , iEnd

    ! Remove leading whitespace.
    String_Strip=""
    iBegin=1
    do while (index(charactersWhiteSpace,string(iBegin:iBegin)) /= 0)
       iBegin=iBegin+1
       if (iBegin > len(string)) return
    end do
    if (iBegin == len(string)) then
       String_Strip=string(iBegin:iBegin)
       return
    end if
    ! Remove trailing whitespace.
    iEnd=len(string)
    do while (index(charactersWhiteSpace,string(iEnd:iEnd)) /= 0)
       iEnd=iEnd-1
       if (iEnd == iBegin) exit
    end do
    String_Strip=string(iBegin:iEnd)
    return
  end function String_Strip

  elemental function Char_Logical(input)
    !!{
    Convert a logical to a string.
    !!}
    implicit none
    character(len=5)                :: Char_Logical
    logical         , intent(in   ) :: input

    if (input) then
       Char_Logical="true"
    else
       Char_Logical="false"
    end if
    return
  end function Char_Logical

 function String_Value_Type(input) result(valueType)
    !!{
    Attempt to detect whether a string corresponds to a floating point number, integer, or other.
    !!}
    implicit none
    type            (enumerationValueTypeType)                :: valueType
    character       (len=*                   ), intent(in   ) :: input
    double precision                                          :: valueFloating
    integer                                                   :: valueInteger , status

    read (input,*,ioStat=status) valueInteger
    if (status == 0) then
       valueType=valueTypeInteger
       return
    end if
    read (input,*,ioStat=status) valueFloating
    if (status == 0) then
       valueType=valueTypeFloating
       return
    end if
    valueType=valueTypeOther
    return
  end function String_Value_Type

  double precision function String_Value_Extract_Float(input,status) result(valueFloat)
    !!{
    Extract a floating-point value from a string.
    !!}
    implicit none
    character(len=*), intent(in   )           :: input
    integer         , intent(  out), optional :: status
    integer                                   :: status_

    read (input,*,ioStat=status_) valueFloat
    if (present(status)) status=status_
    return
  end function String_Value_Extract_Float
  
  integer function String_Value_Extract_Integer(input,status) result(valueInteger)
    !!{
    Extract an integer value from a string.
    !!}
    implicit none
    character(len=*), intent(in   )           :: input
    integer         , intent(  out), optional :: status
    integer                                   :: status_

    read (input,*,ioStat=status_) valueInteger
    if (present(status)) status=status_
    return
  end function String_Value_Extract_Integer

  function String_Value_Extract_Integer_Size_T(input,status) result(valueInteger)
    !!{
    Extract a {\normalfont \ttfamily size\_t} integer value from a string.
    !!}
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    implicit none
    integer  (c_size_t)                          :: valueInteger
    character(len=*   ), intent(in   )           :: input
    integer            , intent(  out), optional :: status
    integer                                      :: status_
    
    read (input,*,ioStat=status_) valueInteger
    if (present(status)) status=status_
    return
  end function String_Value_Extract_Integer_Size_T

  function stringSubstitute(string,find,replace)
    use :: ISO_Varying_String, only : varying_string, assignment(=), len, extract, &
         &                            operator(//)  , operator(==)
    implicit none
    type     (varying_string)                :: stringSubstitute
    type     (varying_string), intent(in   ) :: string
    character(len=1         )                :: find
    character(len=*         )                :: replace
    integer                                  :: i

    stringSubstitute=""
    do i=1,len(string)
       if (extract(string,i,i) == find) then
          stringSubstitute=stringSubstitute//replace
       else
          stringSubstitute=stringSubstitute//extract(string,i,i)
       end if
    end do
    return
  end function stringSubstitute

  function stringXMLFormatVarStr(stringIn,indentStep,indentInitial,forceColor) result (stringOut)
    !!{
    Format an XML string with pretty indentation and coloring. Valid XML strings will be automatically pretty-formatted with one
    element per line, automatic indenting (an initial indent, if required, can be specified via the optional {\normalfont
    \ttfamily indentInitial} argument). Some special formatting codes are supported:
    \begin{description}
      \item[{\normalfont \ttfamily **B}:] Highlight the remainder of the line using bold.
      \item[{\normalfont \ttfamily **C}:] Display a continuation line (to indicate arbitrary additional content), {\normalfont \ttfamily ......}.
    \end{description}
    !!}
    use :: ISO_Varying_String, only : varying_string, char
    implicit none
    type   (varying_string)                          :: stringOut
    type   (varying_string), intent(in   )           :: stringIn
    integer                , intent(in   ), optional :: indentStep, indentInitial
    logical                , intent(in   ), optional :: forceColor

    stringOut=stringXMLFormat(char(stringIn),indentStep,indentInitial,forceColor)
    return
  end function stringXMLFormatVarStr
  
  function stringXMLFormatChar(stringIn,indentStep,indentInitial,forceColor) result (stringOut)
    !!{
    Format an XML string with pretty indentation and coloring. Valid XML strings will be automatically pretty-formatted with one
    element per line, automatic indenting (an initial indent, if required, can be specified via the optional {\normalfont
    \ttfamily indentInitial} argument). Some special formatting codes are supported:
    \begin{description}
      \item[{\normalfont \ttfamily **B}:] Highlight the remainder of the line using bold.
      \item[{\normalfont \ttfamily **C}:] Display a continuation line (to indicate arbitrary additional content), {\normalfont \ttfamily ......}.
    \end{description}
    !!}
    use :: ISO_Varying_String, only : varying_string, len         , assignment(=), operator(//), &
         &                            extract       , operator(==)
    use :: System_Output     , only : stdOutIsATTY
    implicit none
    type     (varying_string)                          :: stringOut
    character(len=*         ), intent(in   )           :: stringIn
    integer                  , intent(in   ), optional :: indentStep            , indentInitial
    logical                  , intent(in   ), optional :: forceColor
    character(len=*         ), parameter               :: ESC         =achar(27)
    integer                                            :: lenStringIn           , i            , &
         &                                                indent
    logical                                            :: inElement             , inTagName    , &
         &                                                startTagName          , inAttribute  , &
         &                                                startValue            , inValue      , &
         &                                                useColor
    character(len=1         )                          :: c
    character(len=8         )                          :: reset
    !![
    <optionalArgument name="indentStep"    defaultsTo="2"      />
    <optionalArgument name="indentInitial" defaultsTo="0"      />
    <optionalArgument name="forceColor"    defaultsTo=".false."/>
    !!]

    ! Determine if color is to be added.
    if (forceColor_) then
       useColor=.true.
    else
#ifdef USEMPI
       useColor=.false.
#else
       useColor=stdOutIsATTY()
#endif
    end if
    ! Begin parsing.
    stringOut   =""
    lenStringIn =len(stringIn)
    inElement   =.false.
    inTagName   =.false.
    inAttribute =.false.
    inValue     =.false.
    startTagName=.false.
    startValue  =.false.
    indent      =indentInitial_
    reset       =ESC//"[0m"
    i           =0
    do while (i < lenStringIn)
       ! Move to the next character.
       i=i+1
       c=extract(stringIn,i,i)
       ! Skip whitespace if not in an element.
       ! Handle white space.
       if (c == " " .and. .not.inElement) cycle
       ! Reset at end of tag names.
       if (inTagName) then
          if (c == " " .or. c == "/" .or. c == ">") then
             if (useColor) stringOut=stringOut//trim(reset)
             inTagName=.false.
          end if
       end if
       ! Skip newlines - we do our own formatting.
       if (c == char(10)) cycle
       ! Detect formatting.
       if (c == "*") then
          if (extract(stringIn,i+1,i+2) == "*B") then
             ! Switch on bold formatting.
             if (useColor) stringOut=stringOut//ESC//"[1m"
             reset=ESC//"[0m"//ESC//"[1m"
             i=i+2
             cycle
          end if
          if (extract(stringIn,i+1,i+2) == "*C") then
             ! Add a continuation line.
             stringOut=stringOut//repeat(" ",indent)//"......"//char(10)
             i=i+2
             cycle
          end if
       end if
       ! Detect element opening.
       if (c == "<") then
          inElement   =.true.
          inTagName   =.true.
          startTagName=.true.
          if (extract(stringIn,i+1,i+1) == "/") then
             ! Closing element.
             indent   =indent-indentStep_
             stringOut=stringOut//repeat(" ",indent)
             stringOut=stringOut//"<"
             c        ="/"
             i        =i+1
          else
             ! Opening element.
             stringOut=stringOut//repeat(" ",indent)
             indent   =indent+indentStep_
          end if
       end if
       ! Detect element closing.
       if (c == "/") then
          if (extract(stringIn,i+1,i+1) == ">") indent=indent-indentStep_
       end if
       if (c == ">") inElement=.false.
       ! Detect attribute name.
       if (inElement .and. .not.inTagName .and. .not.inAttribute .and. c /= " " .and. c /= "/" .and. c /= ">") then
          if (useColor) stringOut=stringOut//ESC//"[33m"
          inAttribute=.true.
       end if
       ! Detect end of attribute name.
       if (inAttribute .and. c == "=") then
          if (useColor) stringOut=stringOut//trim(reset)
       end if
       ! Detect start/end of attribute value.
       if (inAttribute .and. c == '"') then
          if (inValue) then
             ! End of the value.
             if (useColor) stringOut=stringOut//trim(reset)
             inValue    =.false.
             inAttribute=.false.
          else
             ! Start of the value.
             inValue   =.true.
             startValue=.true.
          end if
       end if
       ! Append the current character.
       stringOut=stringOut//c
       ! Add new line.
       if (c == ">") then
          reset=ESC//"[0m"
          if (useColor) stringOut=stringOut//trim(reset)
          if (i < lenStringIn) stringOut=stringOut//char(10)
       end if
       ! Color tag names.
       if (startTagName) then
          if (useColor) stringOut=stringOut//ESC//"[34m"
          startTagName=.false.
       end if
       ! Color attribute values names.
       if (startValue) then
          if (useColor) stringOut=stringOut//ESC//"[32m"
          startValue=.false.
       end if
    end do
    return
  end function stringXMLFormatChar

end module String_Handling
