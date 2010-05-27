!% Contains a module which implements various useful functionality for manipulating character strings.

module String_Handling
  !% Implements various useful functionality for manipulating character strings.
  use ISO_Varying_String
  private
  public :: operator(//), String_Split_Words, String_Count_Words, String_Upper_Case, String_Lower_Case

  interface operator(//)
     module procedure Concatenate_VarStr_Integer
  end interface

  interface String_Split_Words
     module procedure String_Split_Words_VarString
     module procedure String_Split_Words_Char
  end interface

  ! Maximum length of string needed to hold integer values.
  integer, parameter  :: maxIntegerSize=20

  ! Character strings used in converting upper to lower case and vice-versa.
  character(len=*), parameter :: charactersLowerCase='abcdefghijklmnopqrstuvwxyz'
  character(len=*), parameter :: charactersUpperCase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

contains

  integer function String_Count_Words(inputString)
    !% Return a count of the number of space separated words in {\tt inputString}.
    implicit none
    character(len=*), intent(in) :: inputString
    logical                      :: inWord
    integer                      :: iCharacter

    String_Count_Words=0
    inWord=.false.
    do iCharacter=1,len_trim(inputString)
       if (inputString(iCharacter:iCharacter) == " ") then
          if (inWord) String_Count_Words=String_Count_Words+1
          inWord=.false.
       else
          inWord=.true.
       end if
    end do
    if (inWord) String_Count_Words=String_Count_Words+1
    return
  end function String_Count_Words

  subroutine String_Split_Words_VarString(words,inputString)
    !% Split {\tt inputString} into words and return as an array.
    implicit none
    type(varying_string), intent(out), dimension(:) :: words
    character(len=*),     intent(in)                :: inputString
    logical                                         :: inWord
    integer                                         :: iCharacterStart,iCharacter,iWord

    inWord=.false.
    iWord=0
    do iCharacter=1,len_trim(inputString)
       if (inputString(iCharacter:iCharacter) == " ") then
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

  subroutine String_Split_Words_Char(words,inputString)
    !% Split {\tt inputString} into words and return as an array.
    implicit none
    character(len=*), intent(out), dimension(:) :: words
    character(len=*), intent(in)                :: inputString
    logical                                     :: inWord
    integer                                     :: iCharacterStart,iCharacter,iWord

    inWord=.false.
    iWord=0
    do iCharacter=1,len_trim(inputString)
       if (inputString(iCharacter:iCharacter) == " ") then
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

    write (intString,*) intVariable
    Concatenate_VarStr_Integer=varStrVariable//trim(adjustl(intString))
    return
  end function Concatenate_VarStr_Integer

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

end module String_Handling
