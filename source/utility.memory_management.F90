!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module for storing and reporting memory usage by the code.

module Memory_Management
  !% Routines and data type for storing and reporting on memory usage. Also contains routines for allocating and deallocating
  !% arrays with automatic error checking and deallocation at program termination and memory usage reporting.
  use, intrinsic :: ISO_C_Binding
  use Kind_Numbers
  implicit none
  private
  public :: Memory_Usage_Report,Code_Memory_Usage,Alloc_Array,Dealloc_Array,Memory_Usage_Record
#ifdef PROCPS
  public :: Memory_Usage_Get
#endif

  ! Count of number of successive decreases in memory usage.
  integer                 :: successiveDecreaseCount=0

  ! Record of memory usage at the last time it was reported.
  integer(kind=kind_int8) :: usageAtPreviousReport=0

  type memoryUsage
     !% Dervied type variable for storing the properties of a single class of memory storage (memory usage, divisor for outputting
     !% and suffix for outputting)
     integer(kind=kind_int8) :: usage,divisor
     character(len=20)       :: name
     character(len=3)        :: suffix
  end type memoryUsage

  ! Identifiers for the various memory types used.
  integer, public, parameter :: memoryTypeCode =1
  integer, public, parameter :: memoryTypeNodes=2
  integer, public, parameter :: memoryTypeMisc =3
  integer, public, parameter :: memoryTypeTotal=4

  type memoryUsageList
     !% Dervied type variable for storing all memory usage in the code.
     type (memoryUsage) :: memoryType(4)
  end type memoryUsageList

  ! List of all types of memory usage.
  type (memoryUsageList) :: usedMemory=memoryUsageList([memoryUsage(-1,1,'code' ,'??'), &
       &                                                memoryUsage( 0,1,'nodes','??'), &
       &                                                memoryUsage( 0,1,'misc' ,'??'), &
       &                                                memoryUsage( 0,1,'total','??')  &
       &                                               ]                                &
       &                                              )

  ! Overhead memory (in bytes) per allocation.
  integer(kind=kind_int8) :: allocationOverhead=8

#ifdef PROCPS
  interface
     !: ./work/build/utility.memory_usage.o
     function Memory_Usage_Get_C() bind(c,name='Memory_Usage_Get_C')
       !% Template for a C function that returns the current memory usage.
       import
       integer(c_long) :: Memory_Usage_Get_C
     end function Memory_Usage_Get_C
  end interface
#endif

  ! Include automatically generated inferfaces.
  include 'utility.memory_management.precontain.inc'

contains

  subroutine Memory_Usage_Report
    !% Writes a report on the current memory usage. The total memory use is evaluated and all usages are scaled into convenient
    !% units prior to output.
    use ISO_Varying_String
    use Galacticus_Display
    implicit none
    double precision, parameter :: newReportChangeFactor=1.2d0
    logical                     :: issueNewReport
    type(varying_string)        :: headerText,usageText
    character(len=1)            :: join

    !$omp critical (MemAdd)
    issueNewReport=.false.
    usedMemory%memoryType(memoryTypeTotal)%usage=0
    usedMemory%memoryType(memoryTypeTotal)%usage=sum(usedMemory%memoryType(:)%usage)
    if (usageAtPreviousReport <= 0) then ! First call, so always report memory usage.
       issueNewReport=.true.
    else if (usedMemory%memoryType(memoryTypeTotal)%usage > 0) then
       ! Only report memory usage if the total usage has changed by more than newReportChangeFactor.
       if (dabs(dlog(dble(usedMemory%memoryType(memoryTypeTotal)%usage)/dble(usageAtPreviousReport))) > dlog(newReportChangeFactor)) then
          if (usedMemory%memoryType(memoryTypeTotal)%usage < usageAtPreviousReport) then
             successiveDecreaseCount=successiveDecreaseCount+1
             if (successiveDecreaseCount >= 2) then
                issueNewReport         =.true.
                successiveDecreaseCount=0
             end if
          else
             successiveDecreaseCount=0
             issueNewReport         =.true.
          end if
       end if
    end if
    if (issueNewReport) then
       ! Record new memory usage reported.
       usageAtPreviousReport=usedMemory%memoryType(memoryTypeTotal)%usage
       ! Set divisors and suffixes for output.
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeCode ))
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeNodes))
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeMisc ))
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeTotal))
       ! Output the memory usage.
       join=' '
       headerText='Memory: '
       usageText='       ' ! Note that these are non-breaking spaces.
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeCode) ,headerText,usageText,join)
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeNodes),headerText,usageText,join)
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeMisc) ,headerText,usageText,join)
       join='='
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeTotal),headerText,usageText,join)
       call Galacticus_Display_Message(headerText)
       call Galacticus_Display_Message(usageText )
    end if
    !$omp end critical (MemAdd)
    return
  end subroutine Memory_Usage_Report

  subroutine Add_Memory_Component(thisMemoryUsage,headerText,usageText,join)
    !% Add a memory type to the memory reporting strings.
    use ISO_Varying_String
    implicit none
    type(memoryUsage),    intent(in)    :: thisMemoryUsage
    type(varying_string), intent(inout) :: headerText,usageText
    character(len=1),     intent(inout) :: join
    integer                             :: spaceCount
    character(len=20)                   :: formatString
    character(len=len(headerText)+40)   :: temporaryString
    character(len=13)                   :: usageString

    if (thisMemoryUsage%usage > 0) then
       spaceCount=max(0,11-len_trim(thisMemoryUsage%name))
       write (formatString,'(a,i1,a)') '(a,1x,a1,',spaceCount,'x,a)'
       write (temporaryString,formatString) trim(char(headerText)),join,trim(thisMemoryUsage%name)
       headerText=trim(temporaryString)
       write (usageString,'(1x,a1,1x,f7.3,a3)') join,dble(thisMemoryUsage%usage)/dble(thisMemoryUsage%divisor),thisMemoryUsage%suffix
       usageText=trim(usageText)//usageString
       join='+'
    end if
    return
  end subroutine Add_Memory_Component

  subroutine Set_Memory_Prefix(thisMemoryUsage)
    !% Given a memory variable, sets the divisor and suffix required to put the memory usage into convenient units for output.
    use ISO_Varying_String
    implicit none
    type(memoryUsage),       intent(inout) :: thisMemoryUsage
    integer(kind=kind_int8), parameter     :: kilo     =1024
    double precision,        parameter     :: log10kilo=dlog10(dble(kilo))
    integer                                :: usageDecade

    if (thisMemoryUsage%usage > 0) then
       usageDecade=int(dlog10(dble(thisMemoryUsage%usage))/log10kilo+0.01d0)
       select case (usageDecade)
       case (:0)
          thisMemoryUsage%divisor=1
          thisMemoryUsage%suffix='  b'
       case (1)
          thisMemoryUsage%divisor=kilo
          thisMemoryUsage%suffix='kib'
       case (2)
          thisMemoryUsage%divisor=kilo**2
          thisMemoryUsage%suffix='Mib'
       case (3:)
          thisMemoryUsage%divisor=kilo**3
          thisMemoryUsage%suffix='Gib'
       end select
    else
       thisMemoryUsage%divisor=1
       thisMemoryUsage%suffix='  b'
    end if
    return
  end subroutine Set_Memory_Prefix

  subroutine Code_Memory_Usage(codeSizeFile)
    !% If present reads the file {\tt $\langle$executable$\rangle$.size} to determine the amount
    !% of memory the {\tt $\langle$executable$\rangle$.exe} code needs before other memory is
    !% allocated. This is stored to allow an accurate calculation of the
    !% memory used by the code.
    !%
    !% The {\tt $\langle$executable$\rangle$.size} file is made by running the Perl script {\tt Find\_Executable\_Size.pl} (which is done
    !% automatically when the executable is built by {\tt make}).
    use ISO_Varying_String
    use Galacticus_Display
    use Galacticus_Input_Paths
    implicit none
    character(len=*),    intent(in) :: codeSizeFile
    integer                         :: ioError,unitNumber
    double precision                :: dummy
    character(len=80)               :: line
    type(varying_string)            :: codeSizeFileExtension

    usedMemory%memoryType(memoryTypeCode)%usage=0  ! Default value in case size file is unreadable.
    codeSizeFileExtension=char(Galacticus_Input_Path())//'work/build/'//trim(codeSizeFile)
    open (newunit=unitNumber,file=char(codeSizeFileExtension),iostat=ioError,status='old',form='formatted')
    read (unitNumber,'(a80)',iostat=ioError) line ! Read header line.
    if (ioError == 0) then
       read (unitNumber,*) dummy,dummy,dummy,usedMemory%memoryType(memoryTypeCode)%usage
       close (unitNumber)
       return
    else
       call Galacticus_Display_Message('Code size file ['//codeSizeFileExtension//'] not present or unreadable: can be made using Find_Executable_Size.pl')
       close (unitNumber)
       return
    endif
  end subroutine Code_Memory_Usage

  subroutine Memory_Usage_Record(elementsUsed,memoryType,addRemove,blockCount)
    !% Record a change in memory usage.
    use, intrinsic :: ISO_C_Binding
    implicit none
    integer(kind=C_SIZE_T), intent(in)           :: elementsUsed
    integer,                intent(in), optional :: addRemove,memoryType,blockCount
    integer                                      :: memoryTypeActual,blockCountActual,addRemoveActual

    if (present(memoryType)) then
       memoryTypeActual=memoryType
    else
       memoryTypeActual=memoryTypeMisc
    end if
    if (present(addRemove)) then
       addRemoveActual=addRemove
    else
       addRemoveActual=1
    end if
    if (present(blockCount)) then
       blockCountActual=blockCount
    else
       blockCountActual=1
    end if
    !$omp critical(Memory_Management_Usage)
    usedMemory%memoryType(memoryTypeActual)%usage=usedMemory%memoryType(memoryTypeActual)%usage+elementsUsed*addRemoveActual
    usedMemory%memoryType(memoryTypeActual)%usage=usedMemory%memoryType(memoryTypeActual)%usage+sign(blockCountActual,addRemoveActual)*allocationOverhead
    !$omp end critical(Memory_Management_Usage)
    return
  end subroutine Memory_Usage_Record

  include 'utility.memory_management.postcontain.inc'

#ifdef PROCPS
  function Memory_Usage_Get()
    implicit none
    integer(kind_int8) :: Memory_Usage_Get(2)

    usedMemory%memoryType(memoryTypeTotal)%usage=0
    usedMemory%memoryType(memoryTypeTotal)%usage=sum(usedMemory%memoryType(:)%usage)
    Memory_Usage_Get(1)=Memory_Usage_Get_C()*4096_kind_int8
    Memory_Usage_Get(2)=usedMemory%memoryType(memoryTypeTotal)%usage
    return
  end function Memory_Usage_Get
#endif

end module Memory_Management
