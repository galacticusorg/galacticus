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


!% Contains a module which stores file units and finds available file units.

module File_Utilities
  !% Contains a function which returns an available file unit. Also stores the name of the output directory and unit numbers for
  !% various files which remain open throughout.
  use iso_varying_string
  private
  public :: File_Units_Get,Count_Lines_in_File,File_Exists

  interface Count_Lines_in_File
     !% Generic interface for {\tt Count\_Lines\_in\_File} function.
     module procedure Count_Lines_in_File_Char
     module procedure Count_Lines_in_File_VarStr
  end interface

  interface File_Exists
     !% Generic interface for functions that check for a files existance.
     module procedure File_Exists_Char
     module procedure File_Exists_VarStr
  end interface

contains

  integer function File_Units_Get()
    !% Returns the number of an unused file unit. It tries unit numbers from 100 upwards and returns the first free unit.
    !% It aborts if it runs out of valid unit numbers.
    use Galacticus_Error
    implicit none
    integer :: iUnit
    logical :: open, exists

    iUnit=100
    open  =.true.
    exists=.true.
    do while (open.and.exists)
       iUnit=iUnit+1
       inquire(unit=iUnit,opened=open,exist=exists)
    end do

    if (.not.exists) call Galacticus_Error_Report('File_Units_Get','ran out of valid unit numbers')

    File_Units_Get=iUnit

    return
  end function File_Units_Get

  logical function File_Exists_VarStr(FileName)
    !% Checks for existance of file {\tt FileName} (version for varying string argument).
    implicit none
    type(varying_string), intent(in) :: FileName

    File_Exists_VarStr=File_Exists_Char(char(FileName))
    return
  end function File_Exists_VarStr

  logical function File_Exists_Char(FileName)
    !% Checks for existance of file {\tt FileName} (version for character argument).
    implicit none
    character(len=*), intent(in) :: FileName

    inquire(file=FileName,exist=File_Exists_Char)
    return
  end function File_Exists_Char

  integer function Count_Lines_in_File_VarStr(in_file,comment_char)
    !% Returns the number of lines in the file {\tt in\_file} (version for varying string argument).
    use Galacticus_Error
    implicit none
    type (varying_string), intent(in)           :: in_file
    character,             intent(in), optional :: comment_char*(1)

    if (present(comment_char)) then
       Count_Lines_in_File_VarStr=Count_Lines_in_File_Char(char(in_file),comment_char)
    else
       Count_Lines_in_File_VarStr=Count_Lines_in_File_Char(char(in_file))
    end if
    return
  end function Count_Lines_in_File_VarStr

  integer function Count_Lines_in_File_Char(in_file,comment_char)
    !% Returns the number of lines in the file {\tt in\_file} (version for character argument).
    use Galacticus_Error
    implicit none
    character, intent(in)           :: in_file*(*)
    character, intent(in), optional :: comment_char*(1)
    character                       :: first_char*(1)
    integer,   save                 :: io_status,i_unit
    !$omp threadprivate(io_status,i_unit)

    i_unit=File_Units_Get()
    open(i_unit,file=in_file,status='old',form='formatted',iostat=io_status)
    if (io_status.ne.0) then
       write (0,*) 'Count_Lines_in_File(): FATAL - cannot open file ',trim(in_file)
       call Galacticus_Error_Report
    end if
    Count_Lines_in_File_Char=0
    do while (io_status.eq.0)
       read (i_unit,*,iostat=io_status) first_char
       if (io_status.eq.0) then
          if (present(comment_char)) then
             if (first_char.ne.comment_char) Count_Lines_in_File_Char=Count_Lines_in_File_Char+1
          else
             Count_Lines_in_File_Char=Count_Lines_in_File_Char+1
          end if
       end if
    end do
    close(i_unit)
    return
  end function Count_Lines_in_File_Char

end module File_Utilities

