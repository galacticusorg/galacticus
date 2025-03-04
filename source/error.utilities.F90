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
Contains a module which provides error handling utilities.
!!}

module Error_Utilities
  !!{
  Provides error handling utilities.
  !!}
  implicit none
  private
  public :: Error_Wait_Set_From_Parameters

  ! Enumeration of signal numbers.
  !![
  <enumeration>
    <name>signal</name>
    <description>Enumeration of UNIX signal numbers.</description>
    <decodeFunction>yes</decodeFunction>
    <visibility>public</visibility>
    <indexing>1</indexing>
    <entry label="SIGHUP"   /> <!--  1 -->
    <entry label="SIGINT"   /> <!--  2 -->
    <entry label="SIGQUIT"  /> <!--  3 -->
    <entry label="SIGILL"   /> <!--  4 -->
    <entry label="SIGTRAP"  /> <!--  5 -->
    <entry label="SIGABRT"  /> <!--  6 -->
    <entry label="SIGBUS"   /> <!--  7 -->
    <entry label="SIGFPE"   /> <!--  8 -->
    <entry label="SIGKILL"  /> <!--  9 -->
    <entry label="SIGUSR1"  /> <!-- 10 -->
    <entry label="SIGSEGV"  /> <!-- 11 -->
    <entry label="SIGUSR2"  /> <!-- 12 -->
    <entry label="SIGPIPE"  /> <!-- 13 -->
    <entry label="SIGALRM"  /> <!-- 14 -->
    <entry label="SIGTERM"  /> <!-- 15 -->
    <entry label="SIGSTKFLT"/> <!-- 16 -->
    <entry label="SIGCHLD"  /> <!-- 17 -->
    <entry label="SIGCONT"  /> <!-- 18 -->
    <entry label="SIGSTOP"  /> <!-- 19 -->
    <entry label="SIGTSTP"  /> <!-- 20 -->
    <entry label="SIGTTIN"  /> <!-- 21 -->
    <entry label="SIGTTOU"  /> <!-- 22 -->
    <entry label="SIGURG"   /> <!-- 23 -->
    <entry label="SIGXCPU"  /> <!-- 24 -->
  </enumeration>
  !!]
  
contains

  subroutine Error_Wait_Set_From_Parameters(parameters)
    !!{
    Read the parameter that controls the verbosity level, and set that level.
    !!}
    use :: Error           , only : Error_Wait_Set
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (inputParameters), intent(inout) :: parameters
    integer                                 :: errorWaitTime

    !![
    <inputParameter>
      <name>errorWaitTime</name>
      <defaultValue>0</defaultValue>
      <description>The time, in seconds, for which \glc\ should sleep after a fatal error when running under MPI.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    call Error_Wait_Set(errorWaitTime)
    return
  end subroutine Error_Wait_Set_From_Parameters

end module Error_Utilities
