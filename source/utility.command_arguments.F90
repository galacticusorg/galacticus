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


!% Contains a module which provides an interface to read command line arguments of arbitrary type.

module Command_Arguments
  !% Provides an interface to read command line arguments of arbitrary type.
  implicit none
  private
  public :: Get_Argument

  interface Get_Argument
     !% Generic interface to routines that read command line arguments.
     module procedure Get_Argument_Varying_String
     module procedure Get_Argument_Character
     module procedure Get_Argument_Real
     module procedure Get_Argument_Integer
     module procedure Get_Argument_Double
     module procedure Get_Argument_Logical
  end interface

contains

  subroutine Get_Argument_Varying_String(argumentNumber,varStrArgument)
    !% Reads a varying string command line argument.
    use ISO_Varying_String
    implicit none
    integer,               intent(in)  :: argumentNumber    
    type (varying_string), intent(out) :: varStrArgument
    integer                            :: argumentLength
    character(len=10)                  :: shortCharacter

    call Get_Command_Argument(argumentNumber,shortCharacter,length=argumentLength)
    call Get_Temporary_String(argumentNumber,varStrArgument,argumentLength)
    return
  end subroutine Get_Argument_Varying_String

  subroutine Get_Temporary_String(argumentNumber,varStrArgument,argumentLength)
    !% Reads a command line argument into a temporary string of the correct length, and returns it as a varying string.
    use ISO_Varying_String
    implicit none
    integer,                       intent(in)  :: argumentNumber,argumentLength
    type (varying_string),         intent(out) :: varStrArgument
    character(len=argumentLength)              :: characterTemporary

    call Get_Command_Argument(argumentNumber,characterTemporary)
    varStrArgument=trim(characterTemporary)
    return
  end subroutine Get_Temporary_String

  subroutine Get_Argument_Character(argumentNumber,characterArgument)
    !% Reads a character command line argument.
    implicit none
    integer,          intent(in)  :: argumentNumber
    character(len=*), intent(out) :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    return
  end subroutine Get_Argument_Character

  subroutine Get_Argument_Integer(argumentNumber,integerArgument)
    !% Reads a integer command line argument.
    implicit none
    integer,           intent(in)  :: argumentNumber
    integer,           intent(out) :: integerArgument
    character(len=10)              :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    read (characterArgument,*) integerArgument
    return
  end subroutine Get_Argument_Integer

  subroutine Get_Argument_Real(argumentNumber,realArgument)
    !% Reads a real command line argument.
    implicit none
    integer,           intent(in)  :: argumentNumber
    real,              intent(out) :: realArgument
    character(len=20)              :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    read (characterArgument,*) realArgument
    return
  end subroutine Get_Argument_Real

  subroutine Get_Argument_Double(argumentNumber,doubleArgument)
    !% Reads a double command line argument.
    implicit none
    integer,           intent(in)  :: argumentNumber
    double precision,  intent(out) :: doubleArgument
    character(len=20)              :: characterArgument

    call Get_Command_Argument(argumentNumber,characterArgument)
    read (characterArgument,*) doubleArgument
    return
  end subroutine Get_Argument_Double

  subroutine Get_Argument_Logical(argumentNumber,logicalArgument)
    !% Reads a logical command line argument.
    implicit none
    integer,           intent(in)  :: argumentNumber
    logical,           intent(out) :: logicalArgument
    character(len=10)              :: characterArgument
     
    call Get_Command_Argument(argumentNumber,characterArgument)
    read (characterArgument,*) logicalArgument
    return
  end subroutine Get_Argument_Logical

end module Command_Arguments

