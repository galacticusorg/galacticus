!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements various useful functionality for manipulating character strings.

module String_Handling
  !% Implements various useful functionality for manipulating character strings.
  use ISO_Varying_String
  implicit none
  private
  public :: operator(//), String_Split_Words, String_Count_Words, String_Upper_Case, String_Lower_Case, String_Upper_Case_First,&
       & Convert_VarString_To_Char, String_C_to_Fortran, String_Subscript, String_Superscript, String_Levenshtein_Distance,&
       & String_Join

  interface operator(//)
     module procedure Concatenate_VarStr_Integer
     module procedure Concatenate_VarStr_Integer8
  end interface

  interface String_Split_Words
     module procedure String_Split_Words_VarString
     module procedure String_Split_Words_Char
  end interface

  ! Maximum length of string needed to hold integer values.
  integer, parameter :: maxIntegerSize=20
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
    character(len=*         ), intent(in   )           :: inputString
    character(len=*         ), intent(in   ), optional :: separator
    logical                                            :: inWord
    integer                                            :: iCharacter
    type     (varying_string)                          :: separatorActual

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

  subroutine String_Split_Words_VarString(words,inputString,separator,bracketing)
    !% Split {\tt inputString} into words and return as an array.
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

    words=""
    inWord=.false.
    iWord=0
    inBracket=0
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
    !% Split {\tt inputString} into words and return as an array.
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

    words=""
    inWord=.false.
    iWord=0
    inBracket=0
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
    !% Provides a concatenation operator to append an integer number to a {\tt varying\_string}.
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
    !% Provides a concatenation operator to append an integer number to a {\tt varying\_string}.
    use Kind_Numbers
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
    !% Converts an input string to upper case.
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter  , iString

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

  elemental function String_Lower_Case(stringInput) result (stringOutput)
    !% Converts an input string to lower case.
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter  , iString

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
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter

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
    type     (varying_string             ), dimension(:)               , intent(in   ) :: varStrings
    character(len=maxval(len(varStrings))), dimension(size(varStrings))                :: Convert_VarString_To_Char
    integer                                                                            :: iString

    do iString=1,size(varStrings)
       Convert_VarString_To_Char(iString)=varStrings(iString)
    end do
    return
  end function Convert_VarString_To_Char

  function String_C_to_Fortran(charArray)
    !% Convert a C-style character array into a Fortran varying string variable.
    use, intrinsic :: ISO_C_Binding
    implicit none
    type     (varying_string)                :: String_C_to_Fortran
    character(kind=c_char   ), intent(in   ) :: charArray          (:)
    integer                                  :: i

    String_C_to_Fortran=""
    do i=1,size(charArray)
       String_C_to_Fortran=String_C_to_Fortran//charArray(i)
    end do
    return
  end function String_C_to_Fortran

  function String_Subscript(stringInput) result (stringOutput)
    !% Converts an input string to Unicode subscripts.
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter  , iString

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
    character(len=*               ), intent(in   ) :: stringInput
    character(len=len(stringInput))                :: stringOutput
    integer                                        :: iCharacter  , iString

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

  integer function String_Levenshtein_Distance(s,t)
    !% Compute the \href{http://en.wikipedia.org/wiki/Levenshtein_distance}{Levenshtein distance} between strings {\tt a} and {\tt
    !% b}.
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
    !% Joins an array of strings into one long string with the given separator.
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

end module String_Handling
