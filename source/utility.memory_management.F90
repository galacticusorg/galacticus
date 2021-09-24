!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module for storing and reporting memory usage by the code.
!!}

!: $(BUILDPATH)/utility.memory_usage.o

module Memory_Management
  !!{
  Routines and data type for storing and reporting on memory usage. Also contains routines for allocating and deallocating
  arrays with automatic error checking and deallocation at program termination and memory usage reporting.
  !!}
  use            :: Galacticus_Error, only : Galacticus_Error_Report
  use, intrinsic :: ISO_C_Binding   , only : c_double_complex       , c_int    , c_long
  use            :: Kind_Numbers    , only : kind_int8              , kind_quad
  implicit none
  private
  public :: Memory_Usage_Report,Code_Memory_Usage,allocateArray,deallocateArray,Memory_Usage_Record
#ifdef PROCPS
  public :: Memory_Usage_Get
#endif

  ! Code memory size initialization status.
  logical                 :: codeMemoryUsageInitialized=.false.

  ! Count of number of successive decreases in memory usage.
  integer                 :: successiveDecreaseCount=0

  ! Record of memory usage at the last time it was reported.
  integer(kind=kind_int8) :: usageAtPreviousReport  =0

  type memoryUsage
     !!{
     Dervied type variable for storing the properties of a single class of memory storage (memory usage, divisor for outputting
     and suffix for outputting)
     !!}
     integer  (kind=kind_int8) :: divisor, usage
     character(len=20        ) :: name
     character(len=3         ) :: suffix
  end type memoryUsage

  ! Identifiers for the various memory types used.
  integer, parameter, public :: memoryTypeCode   =1
  integer, parameter, public :: memoryTypeNodes  =2
  integer, parameter, public :: memoryTypeMisc   =3
  integer, parameter, public :: memoryTypeUnknown=4
  integer, parameter, public :: memoryTypeTotal  =5

  type memoryUsageList
     !!{
     Dervied type variable for storing all memory usage in the code.
     !!}
     type(memoryUsage) :: memoryType(5)
  end type memoryUsageList

  ! List of all types of memory usage.
  type (memoryUsageList) :: usedMemory=memoryUsageList([memoryUsage(-1,1,'code'   ,'??'), &
       &                                                memoryUsage( 0,1,'nodes'  ,'??'), &
       &                                                memoryUsage( 0,1,'misc'   ,'??'), &
       &                                                memoryUsage( 0,1,'unknown','??'), &
       &                                                memoryUsage( 0,1,'total'  ,'??')  &
       &                                               ]                                  &
       &                                              )

  ! Overhead memory (in bytes) per allocation.
  integer(kind=kind_int8) :: allocationOverhead=8

#ifdef PROCPS
  interface
     function Memory_Usage_Get_C() bind(c,name='Memory_Usage_Get_C')
       !!{
       Template for a C function that returns the current memory usage.
       !!}
       import
       integer(kind=c_long) :: Memory_Usage_Get_C
     end function Memory_Usage_Get_C
  end interface
#endif

  ! Interface to getpagesize() function.
  interface
     function getPageSize() bind(c)
       import c_int
       integer(c_int) :: getPageSize
     end function getPageSize
  end interface

  ! Include automatically generated inferfaces.
  include 'utility.memory_management.preContain.inc'

contains

  subroutine Memory_Usage_Report()
    !!{
    Writes a report on the current memory usage. The total memory use is evaluated and all usages are scaled into convenient
    units prior to output.
    !!}
    use :: Display           , only : displayMessage
    use :: ISO_Varying_String, only : assignment(=) , varying_string
    implicit none
    double precision                , parameter    :: newReportChangeFactor=1.2d0
    logical                                        :: issueNewReport
    type            (varying_string)               :: headerText                 , usageText
    character       (len =1        )               :: join
#ifdef PROCPS
    integer         (kind=kind_int8), dimension(2) :: memoryUsage
#endif
    
    !$omp critical(memoryUsageReport)
    call Code_Memory_Usage()
    issueNewReport=.false.
    usedMemory%memoryType(memoryTypeUnknown)%usage=0_kind_int8
#ifdef PROCPS
    memoryUsage=Memory_Usage_Get()
    usedMemory%memoryType(memoryTypeUnknown)%usage=+memoryUsage(1) &
         &                                         -memoryUsage(2)
#endif
    usedMemory%memoryType(memoryTypeTotal)%usage=sum(usedMemory%memoryType(1:memoryTypeTotal-1)%usage)
    if (usageAtPreviousReport <= 0) then ! First call, so always report memory usage.
       issueNewReport=.true.
    else if (usedMemory%memoryType(memoryTypeTotal)%usage > 0) then
       ! Only report memory usage if the total usage has changed by more than newReportChangeFactor.
       if (abs(log(dble(usedMemory%memoryType(memoryTypeTotal)%usage)/dble(usageAtPreviousReport))) > log(newReportChangeFactor)) then
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
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeCode   ))
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeNodes  ))
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeMisc   ))
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeUnknown))
       call Set_Memory_Prefix(usedMemory%memoryType(memoryTypeTotal  ))
       ! Output the memory usage.
       join=' '
       headerText='Memory: '
       usageText='       ' ! Note that these are non-breaking spaces.
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeCode   ),headerText,usageText,join)
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeNodes  ),headerText,usageText,join)
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeMisc   ),headerText,usageText,join)
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeUnknown),headerText,usageText,join)
       join='='
       call Add_Memory_Component(usedMemory%memoryType(memoryTypeTotal  ),headerText,usageText,join)
       call displayMessage(headerText)
       call displayMessage(usageText )
    end if
    !$omp end critical(memoryUsageReport)
    return
  end subroutine Memory_Usage_Report

  subroutine Add_Memory_Component(memoryUsage_,headerText,usageText,join)
    !!{
    Add a memory type to the memory reporting strings.
    !!}
    use :: ISO_Varying_String, only : assignment(=), char          , len, operator(//), &
          &                           trim         , varying_string
    implicit none
    type     (memoryUsage           ), intent(in   ) :: memoryUsage_
    type     (varying_string        ), intent(inout) :: headerText     , usageText
    character(len=1                 ), intent(inout) :: join
    integer                                          :: spaceCount
    character(len=20                )                :: formatString
    character(len=len(headerText)+40)                :: temporaryString
    character(len=13                )                :: usageString
    
    if (memoryUsage_%usage > 0) then
       !![
       <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
        <description>Internal file I/O in gfortran can be non-thread safe.</description>
       </workaround>
       !!]
#ifdef THREADSAFEIO
       !$omp critical(gfortranInternalIO)
#endif
       spaceCount=max(0,11-len_trim(memoryUsage_%name))
       write (formatString,'(a,i1,a)') '(a,1x,a1,',spaceCount,'x,a)'
       write (temporaryString,formatString) trim(char(headerText)),join,trim(memoryUsage_%name)
       headerText=trim(temporaryString)
       write (usageString,'(1x,a1,1x,f7.3,a3)') join,dble(memoryUsage_%usage)/dble(memoryUsage_%divisor),memoryUsage_%suffix
       usageText=trim(usageText)//usageString
       join='+'
#ifdef THREADSAFEIO
       !$omp end critical(gfortranInternalIO)
#endif
    end if
    return
  end subroutine Add_Memory_Component

  subroutine Set_Memory_Prefix(memoryUsage_)
    !!{
    Given a memory variable, sets the divisor and suffix required to put the memory usage into convenient units for output.
    !!}
    implicit none
    type            (memoryUsage   ), intent(inout) :: memoryUsage_
    integer         (kind=kind_int8), parameter     :: kilo        =1024
    double precision                , parameter     :: log10kilo   =log10(dble(kilo))
    integer                                         :: usageDecade

    if (memoryUsage_%usage > 0) then
       usageDecade=int(log10(dble(memoryUsage_%usage))/log10kilo+0.01d0)
       select case (usageDecade)
       case (:0)
          memoryUsage_%divisor=1
          memoryUsage_%suffix='  b'
       case (1)
          memoryUsage_%divisor=kilo
          memoryUsage_%suffix='kib'
       case (2)
          memoryUsage_%divisor=kilo**2
          memoryUsage_%suffix='Mib'
       case (3:)
          memoryUsage_%divisor=kilo**3
          memoryUsage_%suffix='Gib'
       end select
    else
       memoryUsage_%divisor=1
       memoryUsage_%suffix='  b'
    end if
    return
  end subroutine Set_Memory_Prefix

  subroutine Code_Memory_Usage()
    !!{
    Determines the size of the ``text'' (i.e. code) size.
    !!}
    implicit none
    character(len=80) :: procFileName, pid
    integer           :: ioStatus    , proc         , &
         &               pagesTotal  , pagesResident, &
         &               pagesShared , pagesText    , &
         &               pagesData   , pagesLibrary , &
         &               pagesDirty

    if (.not.codeMemoryUsageInitialized) then
       usedMemory%memoryType(memoryTypeCode)%usage=0  ! Default value in case size file is unreadable.
       write (pid         ,'(i12)'  ) getPID()
       write (procFileName,'(a,a,a)') '/proc/',trim(adjustl(pid)),'/statm'
       open (newUnit=proc,file=procFileName,status='old',form='formatted',ioStat=ioStatus)
       if (ioStatus == 0) then
          read (proc,*) pagesTotal,pagesResident,pagesShared,pagesText,pagesData,pagesLibrary,pagesDirty
          if (ioStatus == 0) usedMemory%memoryType(memoryTypeCode)%usage=pagesText*getPageSize()
       end if
       close(proc)
       codeMemoryUsageInitialized=.true.
    end if
    return
  end subroutine Code_Memory_Usage

  subroutine Memory_Usage_Record(elementsUsed,memoryType,addRemove,blockCount,file,line)
    !!{
    Record a change in memory usage.
    !!}
    use            :: Display           , only : displayMessage, displayVerbosity, verbosityLevelDebug
    use, intrinsic :: ISO_C_Binding     , only : c_size_t
    use            :: ISO_Varying_String, only : assignment(=) , operator(//)    , varying_string
    use            :: String_Handling   , only : operator(//)
    implicit none
    integer  (kind=c_size_t ), intent(in   )           :: elementsUsed
    integer                  , intent(in   ), optional :: addRemove      , blockCount      , line            , memoryType
    character(len=*         ), intent(in   ), optional :: file
    integer                                            :: addRemoveActual, blockCountActual, memoryTypeActual
    integer  (kind_int8     )                          :: accumulation
    type     (varying_string)                          :: message

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
    accumulation=elementsUsed*addRemoveActual+sign(blockCountActual,addRemoveActual)*allocationOverhead
    !$omp atomic
    usedMemory%memoryType(memoryTypeActual)%usage=usedMemory%memoryType(memoryTypeActual)%usage+accumulation
    if (displayVerbosity() >= verbosityLevelDebug) then
       if (present(file).and.present(line)) then
          message='memory record: '
          message=message//elementsUsed*addRemoveActual+sign(blockCountActual,addRemoveActual)*allocationOverhead
          message=message//' ['//file//':'//line//']'
          call displayMessage(message)
       end if
    end if
    return
  end subroutine Memory_Usage_Record

  include 'utility.memory_management.postContain.inc'

#ifdef PROCPS
  function Memory_Usage_Get()
    implicit none
    integer(kind=kind_int8) :: Memory_Usage_Get(2)

    usedMemory%memoryType(memoryTypeTotal)%usage=0
    usedMemory%memoryType(memoryTypeTotal)%usage=sum(usedMemory%memoryType(:)%usage)
    Memory_Usage_Get(1)=Memory_Usage_Get_C()
    Memory_Usage_Get(2)=usedMemory%memoryType(memoryTypeTotal)%usage
    return
  end function Memory_Usage_Get
#endif

end module Memory_Management
