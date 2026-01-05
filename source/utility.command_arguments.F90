!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Contains a module which provides an interface to read command line arguments of arbitrary type.
!!}

module Command_Arguments
  !!{
  Provides an interface to read command line arguments of arbitrary type.
  !!}
  implicit none
  private
  public :: Get_Argument

  interface Get_Argument
     !!{
     Generic interface to routines that read command line arguments.
     !!}
     module procedure Get_Argument_Varying_String
     module procedure Get_Argument_Character
     module procedure Get_Argument_Real
     module procedure Get_Argument_Integer
     module procedure Get_Argument_Double
     module procedure Get_Argument_Logical
  end interface Get_Argument

contains

  subroutine Get_Argument_Varying_String(argumentNumber,varStrArgument)
    !!{
    Reads a varying string command line argument.
    !!}
    use :: ISO_Varying_String, only : varying_string
    implicit none
    integer                  , intent(in   ) :: argumentNumber
    type     (varying_string), intent(  out) :: varStrArgument
    integer                                  :: argumentLength
    character(len=10        )                :: shortCharacter

    call Get_Command_Argument(argumentNumber,shortCharacter,length=argumentLength)
    call Get_Temporary_String(argumentNumber,varStrArgument,argumentLength)
    return
  end subroutine Get_Argument_Varying_String

  subroutine Get_Temporary_String(argumentNumber,varStrArgument,argumentLength)
    !!{
    Reads a command line argument into a temporary string of the correct length, and returns it as a varying string.
    !!}
    use :: ISO_Varying_String, only : varying_string
    implicit none
    integer                      , intent(in   ) :: argumentLength    , argumentNumber
    type     (varying_string    ), intent(  out) :: varStrArgument
    character(len=argumentLength)                :: characterTemporary

    call Get_Command_Argument(argumentNumber,characterTemporary)
    varStrArgument=trim(characterTemporary)
    return
  end subroutine Get_Temporary_String

  subroutine Get_Argument_Character(argumentNumber,characterArgument)
    !!{
    Reads a character command line argument.
    !!}
    implicit none
    integer         , intent(in   ) :: argumentNumber
    character(len=*), intent(  out) :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    return
  end subroutine Get_Argument_Character

  subroutine Get_Argument_Integer(argumentNumber,integerArgument)
    !!{
    Reads a integer command line argument.
    !!}
    implicit none
    integer          , intent(in   ) :: argumentNumber
    integer          , intent(  out) :: integerArgument
    character(len=10)                :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    read (characterArgument,*) integerArgument
    return
  end subroutine Get_Argument_Integer

  subroutine Get_Argument_Real(argumentNumber,realArgument)
    !!{
    Reads a real command line argument.
    !!}
    implicit none
    integer          , intent(in   ) :: argumentNumber
    real             , intent(  out) :: realArgument
    character(len=20)                :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    read (characterArgument,*) realArgument
    return
  end subroutine Get_Argument_Real

  subroutine Get_Argument_Double(argumentNumber,doubleArgument)
    !!{
    Reads a double command line argument.
    !!}
    implicit none
    integer                 , intent(in   ) :: argumentNumber
    double precision        , intent(  out) :: doubleArgument
    character       (len=20)                :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    read (characterArgument,*) doubleArgument
    return
  end subroutine Get_Argument_Double

  subroutine Get_Argument_Logical(argumentNumber,logicalArgument)
    !!{
    Reads a logical command line argument.
    !!}
    implicit none
    integer          , intent(in   ) :: argumentNumber
    logical          , intent(  out) :: logicalArgument
    character(len=10)                :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    read (characterArgument,*) logicalArgument
    return
  end subroutine Get_Argument_Logical

end module Command_Arguments

