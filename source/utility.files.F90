!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which stores file units and finds available file units.

! Specify an explicit dependence on the flock.o object file.
!: $(BUILDPATH)/flock.o

module File_Utilities
  !% Contains a function which returns an available file unit. Also stores the name of the output directory and unit numbers for
  !% various files which remain open throughout.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  !$ use OMP_Lib
  implicit none
  private
  public :: Count_Lines_in_File, File_Exists, File_Lock_Initialize, File_Lock, File_Unlock, Executable_Find

  interface Count_Lines_in_File
     !% Generic interface for {\normalfont \ttfamily Count\_Lines\_in\_File} function.
     module procedure Count_Lines_in_File_Char
     module procedure Count_Lines_in_File_VarStr
  end interface Count_Lines_in_File

  interface File_Exists
     !% Generic interface for functions that check for a files existance.
     module procedure File_Exists_Char
     module procedure File_Exists_VarStr
  end interface File_Exists

  interface
     subroutine flock_C(name,ld,lockIsShared) bind(c,name='flock_C')
       !% Template for a C function that calls {\normalfont \ttfamily flock()} to lock a file.
       import
       character(c_char)        :: name
       type     (c_ptr )        :: ld
       integer  (c_int ), value :: lockIsShared
     end subroutine flock_C
  end interface

  interface
     subroutine funlock_C(ld) bind(c,name='funlock_C')
       !% Template for a C function that calls {\normalfont \ttfamily flock()} to unlock a file.
       import
       type(c_ptr) :: ld
     end subroutine funlock_C
  end interface

  ! Declare the interface for POSIX fsync function.
  interface
     function fsync (fd) bind(c,name="fsync")
       import
       integer(c_int), value :: fd
       integer(c_int) :: fsync
     end function fsync
  end interface

  type, public :: lockDescriptor
     !% Type used to store file lock descriptors.
     private
     type      (c_ptr         ) :: lockDescriptorC
     !$ integer(omp_lock_kind ) :: threadLock
     type      (varying_string) :: fileName
  end type lockDescriptor

contains

  logical function File_Exists_VarStr(FileName)
    !% Checks for existance of file {\normalfont \ttfamily FileName} (version for varying string argument).
    implicit none
    type(varying_string), intent(in   ) :: FileName

    File_Exists_VarStr=File_Exists_Char(char(FileName))
    return
  end function File_Exists_VarStr

  logical function File_Exists_Char(FileName)
    !% Checks for existance of file {\normalfont \ttfamily FileName} (version for character argument).
    implicit none
    character(len=*), intent(in   ) :: FileName

    inquire(file=FileName,exist=File_Exists_Char)
    return
  end function File_Exists_Char

  integer function Count_Lines_in_File_VarStr(in_file,comment_char)
    !% Returns the number of lines in the file {\normalfont \ttfamily in\_file} (version for varying string argument).
    implicit none
    type     (varying_string), intent(in   )           :: in_file
    character(len=1         ), intent(in   ), optional :: comment_char

    if (present(comment_char)) then
       Count_Lines_in_File_VarStr=Count_Lines_in_File_Char(char(in_file),comment_char)
    else
       Count_Lines_in_File_VarStr=Count_Lines_in_File_Char(char(in_file))
    end if
    return
  end function Count_Lines_in_File_VarStr

  integer function Count_Lines_in_File_Char(in_file,comment_char)
    !% Returns the number of lines in the file {\normalfont \ttfamily in\_file} (version for character argument).
    use Galacticus_Error
    implicit none
    character(len=*), intent(in   )           :: in_file
    character(len=1), intent(in   ), optional :: comment_char
    character(len=1)                          :: first_char
    integer         , save                    :: i_unit      , io_status
    !$omp threadprivate(io_status,i_unit)
    open(newunit=i_unit,file=in_file,status='old',form='formatted',iostat=io_status)
    if (io_status /= 0) then
       call Galacticus_Error_Report('Count_Lines_in_File(): FATAL - cannot open file ',trim(in_file))
    end if
    Count_Lines_in_File_Char=0
    do while (io_status == 0)
       read (i_unit,*,iostat=io_status) first_char
       if (io_status == 0) then
          if (present(comment_char)) then
             if (first_char /= comment_char) Count_Lines_in_File_Char=Count_Lines_in_File_Char+1
          else
             Count_Lines_in_File_Char=Count_Lines_in_File_Char+1
          end if
       end if
    end do
    close(i_unit)
    return
  end function Count_Lines_in_File_Char

  subroutine File_Lock_Initialize(lock)
    !% Initilialize the per thread lock in a file lock object.
    implicit none
    type(lockDescriptor), intent(inout) :: lock

    !$ call OMP_Init_Lock(lock%threadLock)
    return
  end subroutine File_Lock_Initialize

  subroutine File_Lock(fileName,lock,lockIsShared)
    !% Place a lock on a file.
    implicit none
    character(len=*         ), intent(in   )           :: fileName
    type     (lockDescriptor), intent(inout)           :: lock
    logical                  , intent(in   ), optional :: lockIsShared
    integer  (c_int)                                   :: lockIsShared_

    lockIsShared_=0
    if (present(lockIsShared).and.lockIsShared) lockIsShared_=1
    ! First obtain a per-thread lock, since POSIX SETLKW locks only per process.
    !$ call OMP_Set_Lock(lock%threadLock)
    ! Now obtain the lock on the file.
    call flock_C(trim(fileName)//".lock"//char(0),lock%lockDescriptorC,lockIsShared_)
    lock%fileName=trim(fileName)
    return
  end subroutine File_Lock

  subroutine File_Unlock(lock)
    !% Remove a lock from a file.
    use Galacticus_Error
    implicit none
    type   (lockDescriptor), intent(inout) :: lock
    integer                                :: fileUnit

    ! First unlock the file.
    call funlock_C(lock%lockDescriptorC)
    open(newUnit=fileUnit,file=char(lock%fileName),status='unknown')
    if (fsync(fnum(fileUnit)) /= 0) call Galacticus_Error_Report('File_Unlock','error syncing file at unlock')
    close(fileUnit)
    ! Then release the per-thread lock.
    !$ call OMP_Unset_Lock(lock%threadLock)
    return
  end subroutine File_Unlock

  function Executable_Find(executableName)
    !% Return the full path to the executable of the given name.
    use ISO_Varying_String
    use String_Handling
    implicit none
    type     (varying_string)                            :: Executable_Find
    character(len=*         ), intent(in   )             :: executableName
    type     (varying_string), allocatable, dimension(:) :: path
    integer                                              :: pathsLength    , pathsStatus, i
    type     (varying_string)                            :: paths

    call Get_Environment_Variable('PATH',length=pathsLength,status=pathsStatus)
    call Get_Paths               (              pathsLength                   )
    allocate(path(String_Count_Words(char(paths),":")))
    call String_Split_Words(path,char(paths),":")
    do i=1,size(path)
       Executable_Find=path(i)//"/"//executableName
       if (File_Exists(Executable_Find)) return
    end do
    Executable_Find=""
    deallocate(path)
    return

  contains
    
    subroutine Get_Paths(pathsLength)
      !% Retrieve the {\normalfont \ttfamily PATH} environment variable.
      implicit none
      integer                     , intent(in   ) :: pathsLength
      character(len=pathsLength+1)                :: pathsName
      
      ! Get the paths.
      call Get_Environment_Variable("PATH",value=pathsName)
      ! Store the paths.
      paths=pathsName
      return
    end subroutine Get_Paths

  end function Executable_Find
  
end module File_Utilities
