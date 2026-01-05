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
Contains a module which provides a class that implements general tasks to be performed by \glc.
!!}

module Task_Evolve_Forests_Work_Shares
  !!{
  Provides a class that implements general tasks to be performed by \glc.
  !!}
  use   , intrinsic :: ISO_C_Binding, only : c_size_t
  !$ use            :: OMP_Lib      , only : OMP_Get_Thread_Num
  private

  !![
  <functionClass>
   <name>evolveForestsWorkShare</name>
   <descriptiveName>Evolve Forests Work Share</descriptiveName>
   <description>Class providing work sharing for the evolve forests task.</description>
   <default>FCFS</default>
   <data>integer :: workerIDOffset_=-1, workerCount_=-1</data>
   <method name="forestNumber" >
    <description>Return the number of the forest to process.</description>
    <type>integer(c_size_t)</type>
    <pass>yes</pass>
    <argument>logical, intent(in   ) :: utilizeOpenMPThreads</argument>
   </method>
   <method name="workerID" >
    <description>Return a unique worker ID.</description>
    <type>integer</type>
    <pass>yes</pass>
    <argument>logical, intent(in   ) :: utilizeOpenMPThreads</argument>
    <code>
     if (self%workerIDOffset_ &lt; 0) call evolveForestsWorkerIDs(self,utilizeOpenMPThreads)
     evolveForestsWorkShareWorkerID=self%workerIDOffset_
     !$ if (utilizeOpenMPThreads) evolveForestsWorkShareWorkerID=evolveForestsWorkShareWorkerID+OMP_Get_Thread_Num()
    </code>
   </method>
   <method name="workerCount" >
    <description>Return the count of workers.</description>
    <type>integer</type>
    <pass>yes</pass>
    <argument>logical, intent(in   ) :: utilizeOpenMPThreads</argument>
    <code>
     if (self%workerIDOffset_ &lt; 0) call evolveForestsWorkerIDs(self,utilizeOpenMPThreads)
     evolveForestsWorkShareWorkerCount=self%workerCount_
    </code>
   </method>
  </functionClass>
  !!]

contains

  subroutine evolveForestsWorkerIDs(self,utilizeOpenMPThreads)
    !!{
    Constructs an ID for this worker which is unique across all MPI and OpenMP threads.
    !!}
#ifdef USEMPI
    use    :: MPI_Utilities, only : mpiSelf
#endif
    !$ use :: OMP_Lib      , only : OMP_Get_Max_Threads
    implicit none
    class  (evolveForestsWorkShareClass), intent(inout)               :: self
    logical                             , intent(in   )               :: utilizeOpenMPThreads
#ifdef USEMPI
    integer                             , allocatable  , dimension(:) :: mpiProcessOpenMPThreadCount, mpiProcessIDOffset
    integer                                                           :: i                          , openMPThreadCount

    ! Determine the number of OpenMP threads available to each MPI process.
    if (utilizeOpenMPThreads) then
       ! Both MPI and OpenMP threads are being used - must account for both types.
       allocate(mpiProcessOpenMPThreadCount(0:mpiSelf%count()-1))
       allocate(mpiProcessIDOffset         (0:mpiSelf%count()-1))
       openMPThreadCount   =1
       !$ openMPThreadCount=OMP_Get_Max_Threads()
       mpiProcessOpenMPThreadCount=mpiSelf%gather(openMPThreadCount)
       ! Compute the MPI process offset.
       mpiProcessIDOffset(0)=0
       do i=1,mpiSelf%count()-1
          mpiProcessIDOffset(i)=+mpiProcessIDOffset         (i-1) &
               &                +mpiProcessOpenMPThreadCount(i-1)
       end do
       ! Compute the ID for this thread.
       self%workerIDOffset_=+mpiProcessIDOffset(mpiSelf%rank())
       ! Compute total number of workers.
       self%workerCount_=sum(mpiProcessOpenMPThreadCount)
    else
       ! OpenMP threads are not being utilized, so ID and count depend only on MPI process.
       ! Compute the ID for this thread.
       self%workerIDOffset_=mpiSelf%rank ()
       ! Compute total number of workers.
       self%workerCount_   =mpiSelf%count()
    end if
#else
    ! No MPI, so ID and count are based just on OpenMP threads (if those are being utilized).
    self   %workerIDOffset_=0
    self   %workerCount_   =1
    if (utilizeOpenMPThreads) then
       !$ self%workerIDOffset_=0
       !$ self%workerCount_   =OMP_Get_Max_Threads()
    end if
#endif
    return
  end subroutine evolveForestsWorkerIDs

end module Task_Evolve_Forests_Work_Shares
