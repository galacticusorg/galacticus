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
Contains a module which sets system resource limits.
!!}

! Specify an explicit dependence on the rlimit.o object file.
!: $(BUILDPATH)/rlimit.o

module System_Limits
  !!{
  Set resource limits.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_int, c_long
  implicit none
  private
  public :: System_Limits_Set

  interface
     function setResourceLimit(resource,limit) bind(c,name='setResourceLimit_C')
       !!{
       Template for a C function that calls {\normalfont \ttfamily setrlimit()} to set a resource limit.
       !!}
       import
       integer(c_int )        :: setResourceLimit
       integer(c_int ), value :: resource
       integer(c_long), value :: limit
     end function setResourceLimit
  end interface

  integer(c_int), bind(C) :: rLimitCPU

contains

  subroutine System_Limits_Set(parameters)
    !!{
    Set system resource limits.
    !!}
    use    :: Error           , only : Error_Report
    use    :: Input_Parameters, only : inputParameter     , inputParameters
    !$ use :: OMP_Lib         , only : OMP_Get_Max_Threads
    implicit none
    type   (inputParameters), intent(inout) :: parameters
    integer(c_int          )                :: status
    integer(c_long         )                :: cpuLimit

    !![
    <inputParameter>
      <name>cpuLimit</name>
      <defaultValue>0_c_long</defaultValue>
      <description>The CPU time limit for the run, in seconds.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (cpuLimit > 0_c_long) then
       !$ cpuLimit=cpuLimit*OMP_Get_Max_Threads()
       status=setResourceLimit(rLimitCPU,cpuLimit)
       if (status /= 0) call Error_Report('failed to set resource limit: CPU'//{introspection:location})
    end if
    return
  end subroutine System_Limits_Set

end module System_Limits
