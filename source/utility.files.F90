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
Contains a module which implements various file-related utilities.
!!}

! Specify an explicit dependence on the C interface files.
!: $(BUILDPATH)/flock.o $(BUILDPATH)/mkdir.o $(BUILDPATH)/unlink.o $(BUILDPATH)/rmdir.o $(BUILDPATH)/rename.o $(BUILDPATH)/access.o

module File_Utilities
  !!{
  Implements various file-related utilities.
  !!}
  use, intrinsic :: ISO_C_Binding     , only : c_char        , c_int, c_ptr
  use            :: ISO_Varying_String, only : varying_string
  use            :: Locks             , only : ompLock
  implicit none
  private
  public :: Count_Lines_in_File, File_Exists           , File_Rename   , File_Lock       , &
       &    File_Unlock        , Executable_Find       , File_Path     , File_Name       , &
       &    File_Name_Temporary, File_Remove           , Directory_Make, File_Name_Expand, &
       &    Directory_Remove   , File_Modification_Time

  interface Count_Lines_in_File
     !!{
     Generic interface for {\normalfont \ttfamily Count\_Lines\_in\_File} function.
     !!}
     module procedure Count_Lines_in_File_Char
     module procedure Count_Lines_in_File_VarStr
  end interface Count_Lines_in_File

  interface File_Exists
     !!{
     Generic interface for functions that check for a files existence.
     !!}
     module procedure File_Exists_Char
     module procedure File_Exists_VarStr
  end interface File_Exists

  interface File_Lock
     !!{
     Generic interface for functions that lock a file
     !!}
     module procedure File_Lock_Char
     module procedure File_Lock_VarStr
  end interface File_Lock

  interface File_Modification_Time
     !!{
     Generic interface for file modification functions.
     !!}
     module procedure File_Modification_Time_Char
     module procedure File_Modification_Time_VarStr
  end interface File_Modification_Time

  interface File_Path
     !!{
     Generic interface for functions that return the path to a file.
     !!}
     module procedure File_Path_Char
     module procedure File_Path_VarStr
  end interface File_Path

  interface File_Name
     !!{
     Generic interface for functions that return the name of a file.
     !!}
     module procedure File_Name_Char
     module procedure File_Name_VarStr
  end interface File_Name

  interface File_Remove
     !!{
     Generic interface for functions that remove a file.
     !!}
     module procedure File_Remove_Char
     module procedure File_Remove_VarStr
  end interface File_Remove

  interface Directory_Remove
     !!{
     Generic interface for functions that remove a directory.
     !!}
     module procedure Directory_Remove_Char
     module procedure Directory_Remove_VarStr
  end interface Directory_Remove

  interface Directory_Make
     !!{
     Generic interface for functions that create a directory.
     !!}
     module procedure Directory_Make_Char
     module procedure Directory_Make_VarStr
  end interface Directory_Make

  interface
     function mkdir_C(name) bind(c,name='mkdir_C')
       !!{
       Template for a C function that calls {\normalfont \ttfamily mkdir()} to make a directory.
       !!}
       import
       integer  (c_int ) :: mkdir_C
       character(c_char) :: name
     end function mkdir_C
  end interface

  interface
     function rmdir_C(name) bind(c,name='rmdir_C')
       !!{
       Template for a C function that calls {\normalfont \ttfamily rmdir()} to remove a directory.
       !!}
       import
       integer  (c_int ) :: rmdir_C
       character(c_char) :: name
     end function rmdir_C
  end interface

  interface
     function unlink_C(name) bind(c,name='unlink_C')
       !!{
       Template for a C function that calls {\normalfont \ttfamily unlink()} to remove a file.
       !!}
       import
       integer  (c_int ) :: unlink_C
       character(c_char) :: name
     end function unlink_C
  end interface

  interface
     function rename_C(nameOld,nameNew) bind(c,name='rename_C')
       !!{
       Template for a C function that calls {\normalfont \ttfamily rename()} to rename a file.
       !!}
       import
       integer  (c_int ) :: rename_C
       character(c_char) :: nameOld , nameNew
     end function rename_C
  end interface

  interface
     function flock_C(name,ld,lockIsShared,timeSleep,countAttempts) bind(c,name='flock_C')
       !!{
       Template for a C function that calls {\normalfont \ttfamily flock()} to lock a file.
       !!}
       import
       integer  (c_int )        :: flock_C
       character(c_char)        :: name
       type     (c_ptr )        :: ld
       integer  (c_int ), value :: lockIsShared , timeSleep, &
            &                      countAttempts
     end function flock_C
  end interface

  interface
     subroutine funlock_C(ld) bind(c,name='funlock_C')
       !!{
       Template for a C function that calls {\normalfont \ttfamily flock()} to unlock a file.
       !!}
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

  ! Declare interface for the POSIX access() function.
  interface
     function access_C(name) bind(c,name='access_C')
       !!{
       Template for a C function that calls {\normalfont \ttfamily access()} to check for file existence.
       !!}
       import
       integer  (c_int ) :: access_C
       character(c_char) :: name
     end function access_C
  end interface

  ! Declare interface for a directory sync function.
  interface
     subroutine syncdir_C(name) bind(c,name='syncdir_C')
       !!{
       Template for a C function that syncs a directory.
       !!}
       import
       character(c_char) :: name
     end subroutine syncdir_C
  end interface

  type, public :: lockDescriptor
     !!{
     Type used to store file lock descriptors.
     !!}
     private
     type(c_ptr         ) :: lockDescriptorC
     type(varying_string) :: fileName
  end type lockDescriptor

  type   (ompLock) :: posixOpenMPFileLock
  logical          :: posixOpenMPFileLockInitialized=.false.
  integer          :: posixOpenMPFileLockCount      =0
  !$omp threadprivate(posixOpenMPFileLockCount)
  
contains

  logical function File_Exists_VarStr(fileName)
    !!{
    Checks for existence of file {\normalfont \ttfamily fileName} (version for varying string argument).
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    type(varying_string), intent(in   ) :: fileName

    File_Exists_VarStr=File_Exists_Char(char(fileName))
    return
  end function File_Exists_VarStr

  logical function File_Exists_Char(fileName)
    !!{
    Checks for existence of file {\normalfont \ttfamily fileName} (version for character argument).
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    character(len=*            ), intent(in   ) :: fileName
    character(len=len(fileName))                :: parentName
    integer                                     :: lengthParent
    !$GLC attributes initialized ::  parentName

    ! Handle empty file names.
    if (fileName == "") then
       File_Exists_Char=.false.
       return
    end if
    ! Find the name of the parent directory.
    parentName  =      fileName
    lengthParent=len(parentName)
    if (parentName(lengthParent:lengthParent) == "/") &
         & parentName=parentName(1:lengthParent-1)
    parentName       =char(File_Path(parentName))
    ! Sync the parent directory to ensure that any file system caching is updated.
    call syncdir_C(trim(parentName)//char(0))
    ! Test for file existence.
    !![
    <workaround type="gfortran" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;ml&#x2F;fortran&#x2F;2019-12&#x2F;msg00012.html">
     <description>Segfault triggered by inquire when running multiple OpenMP threads and a large number of MPI processes. Cause unknown. To workaround this we use the POSIX access() function to test for file existence.</description>
    </workaround>
    !!]
    !! inquire(file=fileName,exist=File_Exists_Char)
    File_Exists_Char=access_C(trim(fileName)//char(0)) == 0
    return
  end function File_Exists_Char

  integer function Count_Lines_in_File_VarStr(in_file,comment_char)
    !!{
    Returns the number of lines in the file {\normalfont \ttfamily in\_file} (version for varying string argument).
    !!}
    use :: ISO_Varying_String, only : char
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
    !!{
    Returns the number of lines in the file {\normalfont \ttfamily in\_file} (version for character argument).
    !!}
    use :: Error, only : Error_Report
    implicit none
    character(len=*), intent(in   )           :: in_file
    character(len=1), intent(in   ), optional :: comment_char
    character(len=1)                          :: first_char
    integer         , save                    :: i_unit      , io_status
    !$omp threadprivate(io_status,i_unit)
    open(newunit=i_unit,file=in_file,status='old',form='formatted',iostat=io_status)
    if (io_status /= 0) then
       call Error_Report("cannot open file '"//trim(in_file)//{introspection:location})
    end if
    Count_Lines_in_File_Char=0
    do while (io_status == 0)
       read (i_unit,'(a1)',iostat=io_status) first_char
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

  subroutine File_Lock_VarStr(fileName,lock,lockIsShared,timeSleep,countAttempts)
    !!{
    Place a lock on a file.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    type   (varying_string), intent(in   )           :: fileName
    type   (lockDescriptor), intent(inout)           :: lock
    logical                , intent(in   ), optional :: lockIsShared
    integer(c_int         ), intent(in   ), optional :: timeSleep   , countAttempts

    call File_Lock(char(fileName),lock,lockIsShared,timeSleep,countAttempts)
    return
  end subroutine File_Lock_VarStr
  
  subroutine File_Lock_Char(fileName,lock,lockIsShared,timeSleep,countAttempts)
    !!{
    Place a lock on a file.
    !!}
    use :: Display                     , only : displayMessage, displayReset, displayMagenta 
    use :: Error                       , only : Error_Report  , Warn
    use :: ISO_Varying_String          , only : assignment(=) , trim
    use :: Numerical_Constants_Prefixes, only : siFormat
    implicit none
    character(len=*         ), intent(in   )           :: fileName
    type     (lockDescriptor), intent(inout)           :: lock
    logical                  , intent(in   ), optional :: lockIsShared
    integer  (c_int         ), intent(in   ), optional :: timeSleep                         , countAttempts
    logical                  , save                    :: fileLockNotAvailableWarned=.false.
    integer  (c_int         )                          :: lockIsShared_                     , status       , &
         &                                                timeWait
    !![
    <optionalArgument name="timeSleep"     defaultsTo="1_c_int" />
    <optionalArgument name="countAttempts" defaultsTo="60_c_int"/>
    !!]
    
    lockIsShared_=0
    if (present(lockIsShared).and.lockIsShared) lockIsShared_=1
    ! Even if using OFD locks we must lock per OpenMP thread since POSIX file locks, since we may have the file open multiple
    ! times with different "descriptions" (see https://www.gnu.org/software/libc/manual/html_node/Open-File-Description-Locks.html).
    if (.not.posixOpenMPFileLockInitialized) then
       !$omp critical(ofdOpenMPInitialize)
       if (.not.posixOpenMPFileLockInitialized) then
          posixOpenMPFileLock           =ompLock()
          posixOpenMPFileLockInitialized=.true.
       end if
       !$omp end critical(ofdOpenMPInitialize)
    end if
    ! Obtain an OpenMP lock, if this thread has no other file locks active.
    if (posixOpenMPFileLockCount == 0) call posixOpenMPFileLock%set()
    posixOpenMPFileLockCount=posixOpenMPFileLockCount+1
    ! Now obtain the lock on the file.
    timeWait     =      0_c_int
    status       =-huge(0_c_int)
    do while (status /= 0)
       status=flock_C(trim(fileName)//".lock"//char(0),lock%lockDescriptorC,lockIsShared_,timeSleep_,countAttempts_)
       select case (status)
       case ( 0_c_int)
          ! Success.
       case (-1_c_int)
          ! File locking is not supported.
          !$omp critical(fileLockNotAvailableWarned)
          if (.not.fileLockNotAvailableWarned) then
             call Warn('file locking is not available - proceeding without locking files')
             fileLockNotAvailableWarned=.true.
          end if
          !$omp end critical(fileLockNotAvailableWarned)
       case (-2_c_int)
          ! File lock was not obtained.
          timeWait=+timeWait       &
               &   +timeSleep_     &
               &   *countAttempts_
          call displayMessage(displayMagenta()//"WARNING:"//displayReset()//" failed to obtain lock on file '"//fileName//"' after "//trim(adjustl(siFormat(dble(timeWait),'f5.1,1x')))//"s - trying again")
       case default
          ! Unknown return code.
          call Error_Report('unknown error code from flock()'//{introspection:location})
       end select
    end do
    lock%fileName=trim(fileName)
    return
  end subroutine File_Lock_Char

  subroutine File_Unlock(lock,sync)
    !!{
    Remove a lock from a file.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    implicit none
    type   (lockDescriptor), intent(inout)           :: lock
    logical                , intent(in   ), optional :: sync
    integer                                          :: fileUnit      , errorStatus, &
         &                                              fileDescriptor
    !![
    <optionalArgument name="sync" defaultsTo=".true." />
    !!]

    if (sync_) then
       open(newUnit=fileUnit,file=char(lock%fileName),status='unknown',iostat=errorStatus)
       if (errorStatus == 0) then
          fileDescriptor=fnum(fileUnit)
          if (fsync(fileDescriptor) /= 0) call Error_Report('error syncing file at unlock'//{introspection:location})
          close(fileUnit)
      end if
    end if
    ! First unlock the file.
    call funlock_C(lock%lockDescriptorC)
    ! Even if using OFD locks we must lock per OpenMP thread since POSIX file locks, since we may have the file open multiple
    ! times with different "descriptions" (see https://www.gnu.org/software/libc/manual/html_node/Open-File-Description-Locks.html).
    ! Release that lock now, if this thread has no more file locks.
    posixOpenMPFileLockCount=posixOpenMPFileLockCount-1
    if (posixOpenMPFileLockCount == 0) call posixOpenMPFileLock%unset()
    return
  end subroutine File_Unlock

  function Executable_Find(executableName)
    !!{
    Return the full path to the executable of the given name.
    !!}
    use :: ISO_Varying_String, only : assignment(=)     , char              , operator(//)
    use :: String_Handling   , only : String_Count_Words, String_Split_Words
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
      !!{
      Retrieve the {\normalfont \ttfamily PATH} environment variable.
      !!}
      use :: ISO_Varying_String, only : assignment(=)
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

  subroutine Directory_Make_VarStr(pathName)
    !!{
    Make the given directory path. Will create intermediate directories in the path if necessary.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    type(varying_string), intent(in   ) :: pathName

    call Directory_Make(char(pathName))
    return
  end subroutine Directory_Make_VarStr

  subroutine Directory_Make_Char(pathName)
    !!{
    Make the given directory path. Will create intermediate directories in the path if necessary.
    !!}
    use :: Error             , only : Error_Report       , Kernel_EACCES, Kernel_ELOOP , Kernel_EMLINK , &
         &                            Kernel_ENAMETOOLONG, Kernel_ENOENT, Kernel_ENOSPC, Kernel_ENOTDIR, &
         &                            Kernel_EROFS
    use :: ISO_Varying_String, only : trim               , var_str
    implicit none
    character(len=* ), intent(in   ) :: pathName
    integer  (c_int )                :: status
    integer                          :: i
    character(len=24)                :: reason

    ! Quick return if the path exists.
    if (File_Exists(var_str(pathName))) return
    do i=2,len(trim(pathName))
       if (pathName(i:i) == "/") then
          if (File_Exists(var_str(pathName(1:i-1)))) cycle
          !$omp critical(mkdir)
          status=mkdir_C(pathName(1:i-1)//char(0))
          if (status /= 0) then
             select case (status)
             case (Kernel_EACCES)
                reason="permission denied"
             case (Kernel_ELOOP       )
                reason="loop in symbolic links"
             case (Kernel_EMLINK      )
                reason="too many links"
             case (Kernel_ENAMETOOLONG)
                reason="name too long"
             case (Kernel_ENOENT      )
                reason="non-existant directory"
             case (Kernel_ENOSPC      )
                reason="no space on file system"
             case (Kernel_ENOTDIR     )
                reason="not a directory"
             case (Kernel_EROFS       )
                reason="read-only file system"
             case default
                reason="unknown"
             end select
             call Error_Report('failed to make intermediate directory (reason: '//trim(reason)//') "'//pathName(1:i-1)//'"'//{introspection:location})
          end if
          !$omp end critical(mkdir)
       end if
    end do
    !$omp critical(mkdir)
    status=mkdir_C(trim(pathName)//char(0))
    if (status /= 0) call Error_Report('failed to make directory "'//trim(pathName)//'"'//{introspection:location})
    !$omp end critical(mkdir)
    return
  end subroutine Directory_Make_Char

  function File_Path_VarStr(fileName)
    !!{
    Returns the path to the file.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    type(varying_string)                :: File_Path_VarStr
    type(varying_string), intent(in   ) :: fileName

    File_Path_VarStr=File_Path(char(fileName))
    return
  end function File_Path_VarStr

  function File_Path_Char(fileName)
    !!{
    Returns the path to the file.
    !!}
    use :: ISO_Varying_String, only : assignment(=), char, extract, index
    implicit none
    type     (varying_string)                :: File_Path_Char
    character(len=*         ), intent(in   ) :: fileName

    if (index(fileName,"/",back=.true.) > 0) then
       File_Path_Char=extract(fileName,1,index(fileName,"/",back=.true.))
    else
       File_Path_Char="./"
    end if
    return
  end function File_Path_Char

  function File_Name_VarStr(fileName)
    !!{
    Returns the path to the file.
    !!}
    use :: ISO_Varying_String, only : char, varying_string
    implicit none
    type(varying_string)                :: File_Name_VarStr
    type(varying_string), intent(in   ) :: fileName

    File_Name_VarStr=File_Name(char(fileName))
    return
  end function File_Name_VarStr

  function File_Name_Char(fileName)
    !!{
    Returns the path to the file.
    !!}
    use :: ISO_Varying_String, only : assignment(=), extract, index, varying_string
    implicit none
    type     (varying_string)                :: File_Name_Char
    character(len=*         ), intent(in   ) :: fileName

    if (index(fileName,"/",back=.true.) > 0) then
       File_Name_Char=extract(fileName,index(fileName,"/",back=.true.)+1,len(fileName))
    else
       File_Name_Char=fileName
    end if
    return
  end function File_Name_Char

  function File_Name_Temporary(fileRootName,path) result(fileName)
    !!{
    Returns the path to the file.
    !!}
#ifdef USEMPI
    use    :: MPI_Utilities     , only : mpiSelf
#endif
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num, OMP_In_Parallel
    use    :: ISO_Varying_String, only : assignment(=)     , operator(//)   , var_str, varying_string
    use    :: String_Handling   , only : operator(//)
    implicit none
    type     (varying_string)                          :: fileName
    character(len=*         ), intent(in   )           :: fileRootName
    character(len=*         ), intent(in   ), optional :: path
    integer                                            :: i

    i=0
    do while (i == 0 .or. File_Exists(fileName))
       i=i+1
       if (present(path)) then
          fileName=path//"/"
       else
          fileName=var_str("/tmp/")
       end if
       fileName=fileName//trim(fileRootName)//"."//GetPID()
       !$ if (OMP_In_Parallel()) fileName=fileName//"."//OMP_Get_Thread_Num()
#ifdef USEMPI
       fileName=fileName//"."//mpiSelf%rankLabel()
#endif
       fileName=fileName//"."//i
    end do
    return
  end function File_Name_Temporary

  subroutine File_Remove_VarStr(fileName)
    !!{
    Remove a file.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    type(varying_string), intent(in   ) :: fileName

    call File_Remove_Char(char(fileName))
    return
  end subroutine File_Remove_VarStr

  subroutine File_Remove_Char(fileName)
    !!{
    Remove a file.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    implicit none
    character(len=*), intent(in   ) :: fileName
    integer  (c_int)                :: status

    if (File_Exists(fileName)) then
       status=unlink_C(trim(fileName)//char(0))
       if (status /= 0) call Error_Report('failed to remove file "'//trim(fileName)//'"'//{introspection:location})
    end if
    return
  end subroutine File_Remove_Char

  subroutine Directory_Remove_VarStr(directoryName)
    !!{
    Remove a directory.
    !!}
    use :: ISO_Varying_String, only : char
    implicit none
    type(varying_string), intent(in   ) :: directoryName

    call Directory_Remove_Char(char(directoryName))
    return
  end subroutine Directory_Remove_VarStr

  subroutine Directory_Remove_Char(directoryName)
    !!{
    Remove a file.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    implicit none
    character(len=*), intent(in   ) :: directoryName
    integer  (c_int)                :: status

    if (File_Exists(directoryName)) then
       status=rmdir_C(trim(directoryName)//char(0))
       if (status /= 0) call Error_Report('failed to remove directory "'//trim(directoryName)//'"'//{introspection:location})
    end if
    return
  end subroutine Directory_Remove_Char
  
  subroutine File_Rename(nameOld,nameNew,overwrite)
    !!{
    Remove a file.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char        , operator(//), trim
    implicit none
    type   (varying_string), intent(in   )           :: nameOld  , nameNew
    logical                , intent(in   ), optional :: overwrite
    integer(c_int         )                          :: status
    !![
    <optionalArgument name="overwrite" defaultsTo=".false."/>
    !!]

    if (overwrite_ .and. File_Exists(nameNew)) call File_Remove(nameNew)
    status=rename_C(char(nameOld)//char(0),char(nameNew)//char(0))
    if (status /= 0) call Error_Report('failed to rename file "'//trim(nameOld)//'" to "'//trim(nameNew)//'"'//{introspection:location})
    return
  end subroutine File_Rename

  function File_Name_Expand(fileNameIn) result(fileNameOut)
    !!{
    Expands placeholders for Galacticus paths in file names.
    !!}
    use :: Input_Paths       , only : inputPath    , pathTypeDataDynamic, pathTypeDataStatic, pathTypeExec
    use :: ISO_Varying_String, only : assignment(=), replace            , index             , extract     , &
         &                            char         , var_str            , operator(//)
    use :: Error             , only : Error_Report
    implicit none
    type     (varying_string)                :: fileNameOut
    character(len=*         ), intent(in   ) :: fileNameIn
    integer                  , parameter     :: iterationMaximum=10000
    type     (varying_string)                :: variableName          , variableContent
    integer                                  :: indexStart            , indexEnd       , &
         &                                      variableLength        , status         , &
         &                                      iteration

    fileNameOut=fileNameIn
    ! Handle custom paths.
    fileNameOut=replace(fileNameOut,"%EXECPATH%"       ,inputPath(pathTypeExec       ),every=.true.)
    fileNameOut=replace(fileNameOut,"%DATASTATICPATH%" ,inputPath(pathTypeDataStatic ),every=.true.)
    fileNameOut=replace(fileNameOut,"%DATADYNAMICPATH%",inputPath(pathTypeDataDynamic),every=.true.)
    ! Handle generic environment variables.
    iteration=0
    do while (index(fileNameOut,"%") /= 0)
       iteration =iteration+1
       if (iteration > iterationMaximum) call Error_Report("infinite loop expanding file name '"//fileNameIn//"':"//char(10)//"   currently: '"//char(fileNameOut)//"'"//{introspection:location})
       indexStart=index(        fileNameOut              ,"%")
       indexEnd  =index(extract(fileNameOut,indexStart+1),"%")
       if (indexEnd == 0) exit
       variableName=extract(fileNameOut,indexStart+1,indexStart+indexEnd-1)
       call Get_Environment_Variable(char(variableName),length=variableLength,status=status)
       if (status == 0) then
          call variableRetrieve(variableContent,char(variableName),variableLength)
       else
          call Error_Report(var_str('environment variable `')//variableName//'` is not defined'//{introspection:location})
       end if
       fileNameOut=extract(fileNameOut,1,indexStart-1)//variableContent//extract(fileNameOut,indexStart+indexEnd+1)
    end do
    return

  contains

    subroutine variableRetrieve(variableContent,variableName,variableLength)
      !!{
      Retrieve an environment variable
      !!}
      use :: ISO_Varying_String, only : assignment(=)
      implicit none
      type     (varying_string    ), intent(  out) :: variableContent
      integer                      , intent(in   ) :: variableLength
      character(len=*             ), intent(in   ) :: variableName
      character(len=variableLength)                :: variableContent_
      
      call Get_Environment_Variable(variableName,value=variableContent_)
      variableContent=variableContent_
      return
    end subroutine variableRetrieve
    
  end function File_Name_Expand

  function File_Modification_Time_VarStr(fileName,status) result(timeModification)
    !!{
    Return the modification time of the named file.
    !!}
    use :: ISO_Varying_String, only : varying_string, char
    implicit none
    character(len=30        )                          :: timeModification
    type     (varying_string), intent(in   )           :: fileName
    integer                  , intent(  out), optional :: status
    
    timeModification=File_Modification_Time(char(fileName),status)
    return
  end function File_Modification_Time_VarStr
  
  function File_Modification_Time_Char(fileName,status) result(timeModification)
    !!{
    Return the modification time of the named file.
    !!}
    use :: Error, only : Error_Report, errorStatusSuccess, errorStatusNotExist
    implicit none
    character(len=30)                          :: timeModification
    character(len=* ), intent(in   )           :: fileName
    integer          , intent(  out), optional :: status
    integer          , dimension(13)           :: statValues

    if (present(status)) status=errorStatusSuccess
    if (File_Exists(fileName)) then       
       call stat(fileName,statValues)
       timeModification=ctime(statValues(10))
    else
       timeModification=""
       if (present(status)) then
          status=errorStatusNotExist
       else
          call Error_Report('file "'//fileName//'" does not exist'//{introspection:location})
       end if
    end if
    return
  end function File_Modification_Time_Char
  
end module File_Utilities
