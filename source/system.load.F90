!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which reports system load averages.

module System_Load
  !% Reports system load averages.
  implicit none
  private
  public :: System_Load_Get, System_Processor_Count

contains

  subroutine System_Load_Get(loadAverage,activeTasks,totalTasks)
    !% Reports on system load via {\tt /proc/loadavg}.
    implicit none
    double precision         , intent(  out) :: loadAverage       (3)                 
    integer                  , intent(  out) :: activeTasks          , totalTasks     
    integer                                  :: lUnit                , slashPosition  
    character       (len=256)                :: loadAverageEncoded                    
    
    !$omp critical (System_Load_Read)                                                                               
    open(newUnit=lUnit,file="/proc/loadavg",status='old',form='formatted')
    read (lUnit,'(a)') loadAverageEncoded
    close(lUnit)
    !$omp end critical (System_Load_Read)
    slashPosition=index(loadAverageEncoded,"/")
    read (loadAverageEncoded(:slashPosition-1 ),*) loadAverage,activeTasks
    read (loadAverageEncoded( slashPosition+1:),*) totalTasks
    return
  end subroutine System_Load_Get

  integer function System_Processor_Count()
    !% Return a count of the number of available processors.
    implicit none
    integer           :: ioError , lUnit  
    character(len=32) :: infoLine         
    
    !$omp critical (System_CPU_Read)                                   
    System_Processor_Count=0
    open(newUnit=lUnit,file="/proc/cpuinfo",status='old',form='formatted',iostat=ioError)
    if (ioError == 0) then
       do while (ioError == 0)
          read (lUnit,'(a)',iostat=ioError) infoLine
          if (ioError == 0 .and. infoLine(1:9) == "processor") System_Processor_Count=System_Processor_Count+1
       end do
       close(lUnit)
    end if
    !$omp end critical (System_CPU_Read)
    return
  end function System_Processor_Count

end module System_Load
