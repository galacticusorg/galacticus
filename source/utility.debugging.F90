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
Contains a module that implements useful debugging utilities.
!!}

! Specify an explicit dependence on the backtrace.o object file.
!: $(BUILDPATH)/backtrace.o
  
module Debugging
  !!{
  Implements useful debugging utilities.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_int, c_char
  private
  public :: debugLog , debugOn, debugOff, isDebugging, &
       &    getCaller

  ! Debugging status.
  logical :: debugging_          =.false.
  !$omp threadprivate(debugging_)
  
  ! File unit for output of debugging information.
  logical :: debugFileInitialized=.false.
  integer :: debugFile
  !$omp threadprivate(debugFile,debugFileInitialized)

  interface
     subroutine getCallerC(callerSize,caller) bind(c,name='getCallerC')
       !!{
       Template for a C function that returns the name of the caller function.
       !!}
       import
       integer  (kind=c_int ), value :: callerSize
       character(kind=c_char)        :: caller    (callerSize)
     end subroutine getCallerC
  end interface

contains

  subroutine initialize()
    !!{
    Initialize debug output.
    !!}
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num
    implicit none
    character(len=32) :: fileName
    
    if (.not.debugFileInitialized) then
       write    (fileName,'(a)'       ) 'debug.log'
       !$ write (fileName,'(a,i4.4,a)') 'debug_thread',OMP_Get_Thread_Num(),'.log'
       open(newUnit=debugFile,file='debug_thread.log',status='unknown',form='formatted')
       debugFileInitialized=.true.
    end if
    return
  end subroutine initialize

  subroutine debugOn()
    !!{
    Switch debugging on.
    !!}

    debugging_=.true.
    return
  end subroutine debugOn

  subroutine debugOff()
    !!{
    Switch debugging off.
    !!}

    debugging_=.false.
    return
  end subroutine debugOff

  logical function isDebugging()
    !!{
    Return true if debugging is on.
    !!}

    isDebugging=debugging_
    return
  end function isDebugging
  
  subroutine debugLog(message)
    !!{
    Write a message to the debug log.
    !!}
    use :: ISO_Varying_String, only : varying_string, char
    implicit none
    type(varying_string), intent(in   ) :: message

    call initialize()
    write (debugFile,'(a)') char(message)
    return
  end subroutine debugLog
  
  function getCaller()
    !!{
    Return the name of the calling function.
    !!}
    use :: ISO_Varying_String, only : varying_string     , index, extract
    use :: String_Handling   , only : String_C_to_Fortran
    implicit none
    type     (varying_string)                 :: getCaller , caller_
    character(kind=c_char   ), dimension(256) :: caller
    integer                                   :: indexStart, indexEnd

    call getCallerC(256,caller)
    caller_   =String_C_to_Fortran(caller)
    indexStart=index(caller_,"(")
    indexEnd  =index(caller_,"+")
    getCaller=extract(caller_,indexStart+1,indexEnd-1)
    return
  end function getCaller
  
end module Debugging
