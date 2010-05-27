!! Copyright 2009, Andrew Benson <abenson@caltech.edu>
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




!% Contains a module for storing and reporting memory usage by the code.

module Memory_Management
  !% Routines and data type for storing and reporting on memory usage. Also contains routines for allocating and deallocating
  !% arrays with automatic error checking and deallocation at program termination and memory usage reporting.
  !$ use omp_lib
  include 'utility.memory_management.use.inc'
  private
  public :: Memory_Management_Type,Report_Memory_Usage,Memory_Type,Free_Memory,Code_Memory_Usage,Alloc_Array,Dealloc_Array&
       &,Memory_Usage_Record

  integer          :: successive_decreases=0
  double precision :: last_report_usage=0.0d0

  type Memory_Type
     !% Dervied type variable for storing the properties of a single class of memory storage (memory usage, divider for outputting
     !% and suffix for outputting)
     double precision   :: usage,divider
     character (len=20) :: name
     character (len=2)  :: suffix
  end type Memory_Type

  ! Identifiers for the various memory types used
  integer, public, parameter :: MEM_TYPE_CODE     =1
  integer, public, parameter :: MEM_TYPE_TREE     =2
  integer, public, parameter :: MEM_TYPE_GALAXY   =3
  integer, public, parameter :: MEM_TYPE_SATELLITE=4
  integer, public, parameter :: MEM_TYPE_SSP      =5
  integer, public, parameter :: MEM_TYPE_LEVELS   =6
  integer, public, parameter :: MEM_TYPE_SAMPLE   =7
  integer, public, parameter :: MEM_TYPE_MISC     =8
  integer, public, parameter :: MEM_TYPE_TOTAL    =9

  type Memory_Management_Type
     !% Dervied type variable for storing all memory usage in the code.
     type (Memory_Type) :: mem_type(9)
  end type Memory_Management_Type

  type (Memory_Management_Type) :: Memory_Used=Memory_management_Type( (/Memory_Type(-1.0,1.0,'code'     ,'??'),&
       & Memory_Type( 0.0 ,1.0,'tree'     ,'??'), Memory_Type( 0.0,1.0,'galaxy'   ,'??'), Memory_Type( 0.0,1.0,'satellite','??'),&
       & Memory_Type( 0.0 ,1.0,'ssp'      ,'??'), Memory_Type( 0.0,1.0,'levels'   ,'??'), Memory_Type( 0.0,1.0,'sample'   ,'??'),&
       & Memory_Type( 0.0 ,1.0,'misc'     ,'??'), Memory_Type( 0.0,1.0,'total'    ,'??')/))

  include 'utility.memory_management.precontain.inc'

contains
  subroutine Report_Memory_Usage
    !% Writes a report on the current memory usage. The total memory use is evaluated and all usages are scaled into convenient
    !% units prior to output.
    use ISO_Varying_String
    use Galacticus_Display
    implicit none
    double precision, parameter :: total_change_factor=1.2d0
    logical                     :: report_now
    type (Memory_Type)          :: galaxy_total
    type (varying_string)       :: str1,str2
    character(len=1)            :: join

    !$omp critical (MemAdd)
    report_now=.false.
    galaxy_total%usage=Memory_Used%mem_type(MEM_TYPE_GALAXY)%usage+Memory_Used%mem_type(MEM_TYPE_SATELLITE)%usage
    galaxy_total%name='galaxy'
    Memory_Used%mem_type(MEM_TYPE_TOTAL)%usage=0.0
    Memory_Used%mem_type(MEM_TYPE_TOTAL)%usage=sum(Memory_Used%mem_type(:)%usage)
    if (last_report_usage.le.0.0) then ! First call, so always report memory usage.
       report_now=.true.
    else if (Memory_Used%mem_type(MEM_TYPE_TOTAL)%usage.gt.0.0) then
       ! Only report memory usage if the total usage has changed by more than total_change_factor.
       if (abs(log(Memory_Used%mem_type(MEM_TYPE_TOTAL)%usage/last_report_usage)).gt.log(total_change_factor)) then
          if (Memory_Used%mem_type(MEM_TYPE_TOTAL)%usage.lt.last_report_usage) then
             successive_decreases=successive_decreases+1
             if (successive_decreases.ge.2) then
                report_now=.true.
                successive_decreases=0
             end if
          else
             successive_decreases=0
             report_now=.true.
          end if
       end if
    endif
    if (report_now) then
       last_report_usage=Memory_Used%mem_type(MEM_TYPE_TOTAL)%usage
       call Set_Memory_Prefix(galaxy_total)
       call Set_Memory_Prefix(Memory_Used%mem_type(MEM_TYPE_CODE  )) ! Set dividers and suffixes for output.
       call Set_Memory_Prefix(Memory_Used%mem_type(MEM_TYPE_TREE  ))
       call Set_Memory_Prefix(Memory_Used%mem_type(MEM_TYPE_SSP   ))
       call Set_Memory_Prefix(Memory_Used%mem_type(MEM_TYPE_LEVELS))
       call Set_Memory_Prefix(Memory_Used%mem_type(MEM_TYPE_SAMPLE))
       call Set_Memory_Prefix(Memory_Used%mem_type(MEM_TYPE_MISC  ))
       call Set_Memory_Prefix(Memory_Used%mem_type(MEM_TYPE_TOTAL ))
       ! Output the memory usage.
       join=' '
       str1='Memory: '
       str2='       ' ! Note that these are non-breaking spaces.
       call Add_Memory_Component(Memory_Used%mem_type(MEM_TYPE_CODE)  ,str1,str2,join)
       call Add_Memory_Component(Memory_Used%mem_type(MEM_TYPE_TREE)  ,str1,str2,join)
       call Add_Memory_Component(galaxy_total                         ,str1,str2,join)
       call Add_Memory_Component(Memory_Used%mem_type(MEM_TYPE_SSP)   ,str1,str2,join)
       call Add_Memory_Component(Memory_Used%mem_type(MEM_TYPE_LEVELS),str1,str2,join)
       call Add_Memory_Component(Memory_Used%mem_type(MEM_TYPE_SAMPLE),str1,str2,join)
       call Add_Memory_Component(Memory_Used%mem_type(MEM_TYPE_MISC)  ,str1,str2,join)
       join='='
       call Add_Memory_Component(Memory_Used%mem_type(MEM_TYPE_TOTAL) ,str1,str2,join)
       call Galacticus_Display_Message(str1)
       call Galacticus_Display_Message(str2)
    end if
    !$omp end critical (MemAdd)
  end subroutine Report_Memory_Usage

  subroutine Add_Memory_Component(Memory,str1,str2,join)
    !% Add a memory type to the memory reporting strings.
    use iso_varying_string
    implicit none
    type (Memory_Type),    intent(in)    :: Memory
    type (varying_string), intent(inout) :: str1,str2
    character(len=1),      intent(inout) :: join
    integer                              :: ispace
    character(len=20)                    :: format_string
    character(len=len(str1)+40)          :: tmp_string

    if (Memory%usage.gt.0) then
       ispace=max(0,10-len_trim(Memory%name))
       write (format_string,'(a,i1,a)') '(a,1x,a1,',ispace,'x,a)'
       write (tmp_string,format_string) trim(char(str1)),join,trim(Memory%name)
       str1=trim(tmp_string)
       write (tmp_string,'(1x,a1,1x,f7.3,a2)') join,Memory%usage/Memory%divider,Memory%suffix
       str2=trim(str2)//trim(tmp_string)
       join='+'
    end if
    return
  end subroutine Add_Memory_Component

  subroutine Set_Memory_Prefix(This_Memory)
    !% Given a memory variable, sets the divider and suffix required to put the memory usage into convenient units for output.
    implicit none
    type (Memory_Type), intent(inout)   :: This_Memory
    double precision,   parameter       :: log10_Kilo=3.0103d0 ! log_10 of 1024
    double precision,   parameter       :: kilo=1024.0d0
    integer                             :: usage_base

    if (This_Memory%usage.gt.0.0) then
       usage_base=int(log10(This_Memory%usage)/log10_Kilo+0.01)
       select case (usage_base)
       case (:0)
          This_Memory%divider=1.0
          This_Memory%suffix='b '
       case (1)
          This_Memory%divider=kilo
          This_Memory%suffix='kb'
       case (2)
          This_Memory%divider=kilo**2
          This_Memory%suffix='Mb'
       case (3:)
          This_Memory%divider=kilo**3
          This_Memory%suffix='Gb'
       end select
    else
       This_Memory%divider=1.0
       This_Memory%suffix='b '
    end if
    return
  end subroutine Set_Memory_Prefix

  subroutine Code_Memory_Usage(sizefile)
    !% If present reads the file {\tt $\langle$executable$\rangle$.size} to determine the amount
    !% of memory the {\tt $\langle$executable$\rangle$.exe} code needs before other memory is
    !% allocated. This is stored to allow an accurate calculation of the
    !% memory used by the code.
    !%
    !% The {\tt $\langle$executable$\rangle$.size} file is made by running the Perl script {\tt Find\_Executable\_Size.pl} (which is done
    !% automatically when the executable is built by {\tt make}).
    use File_Utilities
    use ISO_Varying_String
    implicit none
    character(len=*), intent(in) :: sizefile
    integer                      :: file_err,iunit
    double precision             :: dummy
    character                    :: line*80
    type (varying_string)        :: sizefile_ext

    Memory_Used%mem_type(MEM_TYPE_CODE)%usage=0.0  ! Default value in case size file is unreadable.
    iunit=File_Units_Get()
    sizefile_ext='./work/build/'//trim(sizefile)
    open (iunit,file=char(sizefile_ext),iostat=file_err,status='old',form='formatted')
    read (iunit,'(a80)',iostat=file_err) line ! Read header line.
    if (file_err.eq.0) then
       read (iunit,*) dummy,dummy,dummy,Memory_Used%mem_type(MEM_TYPE_CODE)%usage
       close (iunit)
       return
    else
#ifdef WARN
       write (0,*) 'Code_Memory_Usage(): WARNING - ',trim(char(sizefile_ext)),' not present or unreadable'
       write (0,*) '                     this file can be made using Find_Executable_Size.pl'
#endif
       close (iunit)
       return
    endif
  end subroutine Code_Memory_Usage

  subroutine Memory_Usage_Record(elementsUsed,memoryType,addRemove)
    !% Record a change in memory usage.
    implicit none
    integer, intent(in)           :: elementsUsed
    integer, intent(in), optional :: memoryType,addRemove
    integer                       :: memoryTypeActual,addRemoveActual

    if (present(memoryType)) then
       memoryTypeActual=memoryType
    else
       memoryTypeActual=MEM_TYPE_MISC
    end if
    if (present(addRemove)) then
       addRemoveActual=addRemove
    else
       addRemoveActual=1
    end if
    !$omp atomic
    Memory_Used%mem_type(memoryTypeActual)%usage=Memory_Used%mem_type(memoryTypeActual)%usage+dble(elementsUsed*addRemoveActual)
    return
  end subroutine Memory_Usage_Record

  include 'utility.memory_management.postcontain.inc'

end module Memory_Management
