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
Contains a module for storing and reporting memory usage by the code.
!!}

  ! Add a dependence on the file providing the interface to mallinfo2().
  !: $(BUILDPATH)/utility.memory_reporting.mallinfo2.o

module Memory_Reporting
  !!{
  Provide reporting functions for memory usage.
  !!}
  use            :: Error        , only : Error_Report
  use, intrinsic :: ISO_C_Binding, only : c_size_t    , c_int
  use            :: Locks        , only : ompLock
  implicit none
  private
  public :: reportMemoryUsage

  ! Code memory size initialization status.
  logical           :: codeMemoryUsageInitialized =.false.

  ! Lock used to coordinate memory reporting.
  type   (ompLock ) :: memoryUsageLock
  logical           :: memoryUsageLockInitialized =.false.

  ! Count of number of successive decreases in memory usage.
  integer           :: successiveDecreaseCount    =0

  ! Record of memory usage at the last time it was reported.
  integer(c_size_t) :: memoryUsageAtPreviousReport=0

  ! Record of maximum memory usage.
  integer(c_size_t) :: memoryUsageMaximum         =0

  ! Record of code size and available memory.
  integer(c_size_t) :: memoryUsageCode                    , memoryAvailable

  ! Interface to getpagesize() function.
  interface
     function getPageSize() bind(c)
       import c_int
       integer(c_int) :: getPageSize
     end function getPageSize
  end interface

  ! Interface to mallinfo2_c() function.
  interface
     function mallinfo2_c() bind(c)
       import c_size_t
       integer(c_size_t) :: mallinfo2_c
     end function mallinfo2_c
  end interface

  ! Interface to getTotalSystemMemory() function.
  interface
     function getTotalSystemMemory() bind(c)
       import c_size_t
       integer(c_size_t) :: getTotalSystemMemory
     end function getTotalSystemMemory
  end interface

contains

  subroutine reportMemoryUsage()
    !!{
    Writes a report on the current memory usage.
    !!}
    use :: Display           , only : displayMessage, displayRed  , displayYellow, displayGreen, &
         &                            displayBold   , displayReset
    use :: ISO_Varying_String, only : varying_string, operator(//), assignment(=)
    implicit none
    double precision                , parameter :: newReportChangeFactor=1.2d0
    logical                                     :: issueNewReport
    type            (varying_string)            :: usageText
    integer         (c_size_t      )            :: memoryUsage                , divisor
    character       (len =2        )            :: suffix
    character       (len =7        )            :: label
    double precision                            :: memoryFraction

    if (.not.memoryUsageLockInitialized) then
       !$omp critical(memoryUsageReport)
       if (.not.memoryUsageLockInitialized) then
          memoryUsageLock           =ompLock()
          memoryUsageLockInitialized=.true.
       end if
       !$omp end critical(memoryUsageReport)
    end if
    ! Attempt to get a lock to coordinate memory usage reporting - if the lock is held by another thread, just return (no need for
    ! us to report memory also).
    if (.not.memoryUsageLock%setNonBlocking()) return
    ! Ensure that we have the code memory usage.
    call codeUsageGet()
    ! Get the current memory usage.
    memoryUsage=+memoryUsageCode   &
         &      +mallinfo2_C    ()
    ! Record the maximum memory usage.
    memoryUsageMaximum=max(memoryUsageMaximum,memoryUsage)
    ! Decide whether to report.
    issueNewReport=.false.
    if (memoryUsageAtPreviousReport <= 0) then ! First call, so always report memory usage.
       issueNewReport=.true.
    else if (memoryUsage > 0) then
       ! Only report memory usage if the total usage has changed by more than newReportChangeFactor.
       if (abs(log(dble(memoryUsage)/dble(memoryUsageAtPreviousReport))) > log(newReportChangeFactor)) then
          if (memoryUsage < memoryUsageAtPreviousReport) then
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
       memoryUsageAtPreviousReport=memoryUsage
       ! Generate the report.
       usageText     ='Memory usage: '
       memoryFraction=dble(memoryUsage)/dble(memoryAvailable)
       if      (memoryFraction < 0.25d0) then
          usageText=usageText//displayGreen ()
       else if (memoryFraction < 0.50d0) then
          usageText=usageText//displayYellow()
       else if (memoryFraction < 0.75d0) then
          usageText=usageText//displayRed   ()
       else
          usageText=usageText//displayRed   ()//displayBold()
       end if
       call getSuffix(memoryUsage    ,divisor,suffix)
       write (label,'(f7.3)') dble(memoryUsage    )/dble(divisor)
       usageText=usageText//trim(adjustl(label))//' '//trim(adjustl(suffix))//' / '
       call getSuffix(memoryAvailable,divisor,suffix)
       write (label,'(f7.3)') dble(memoryAvailable)/dble(divisor)
       usageText=usageText//trim(adjustl(label))//' '//trim(adjustl(suffix))//displayReset()
       ! Display the report.
       call displayMessage(usageText)
    end if
    call memoryUsageLock%unset()
    return
  end subroutine reportMemoryUsage

  subroutine getSuffix(memoryUsage,divisor,suffix)
    !!{
    Compute an appropriate suffix (and divisor) for a given memory usage.
    !!}
    implicit none
    integer         (c_size_t), intent(in   ) :: memoryUsage
    integer         (c_size_t), intent(  out) :: divisor
    character       (len =2  ), intent(  out) :: suffix
    integer         (c_size_t), parameter     :: kilo       =1024
    double precision          , parameter     :: log10kilo  =log10(dble(kilo))
    integer                                   :: usageDecade

    if (memoryUsage > 0) then
       usageDecade=int(log10(dble(memoryUsage))/log10kilo+0.01d0)
       select case (usageDecade)
       case (:0)
          divisor=1
          suffix=' B'
       case (1)
          divisor=kilo
          suffix='KB'
       case (2)
          divisor=kilo**2
          suffix='MB'
       case (3)
          divisor=kilo**3
          suffix='GB'
       case (4)
          divisor=kilo**4
          suffix='TB'
       case (5)
          divisor=kilo**4
          suffix='PB'
       case (6)
          divisor=kilo**4
          suffix='EB'
       case (7)
          divisor=kilo**4
          suffix='ZB'
       case (8:)
          divisor=kilo**4
          suffix='YB'
       end select
    else
       divisor=1
       suffix=' B'
    end if
    return
  end subroutine getSuffix

  subroutine codeUsageGet()
    !!{
    Determines the size of the ``text'' (i.e. code) size.
    !!}
    use :: Display      , only : displayMessage, displayGreen , displayReset
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    implicit none
    integer  (c_size_t), parameter :: kilo             =1024
    integer  (c_size_t)            :: divisor
    character(len= 2  )            :: suffix
    character(len=80  )            :: procFileName          , pid                 , &
         &                            slurmMemoryString     , sourceMemory
    integer                        :: ioStatus              , proc                , &
         &                            pagesTotal            , pagesResident       , &
         &                            pagesShared           , pagesText           , &
         &                            pagesData             , pagesLibrary        , &
         &                            pagesDirty            , statusSlurmPerNode  , &
         &                            statusSlurmPerCPU     , statusSlurmCountCPUs, &
         &                            countCPUs

    if (.not.codeMemoryUsageInitialized) then
       ! Find memory used by code itself.
       memoryUsageCode=0_c_size_t  ! Default value in case size file is unreadable.
       write (pid         ,'(i12)'  ) getPID()
       write (procFileName,'(a,a,a)') '/proc/',trim(adjustl(pid)),'/statm'
       open (newUnit=proc,file=procFileName,status='old',form='formatted',ioStat=ioStatus)
       if (ioStatus == 0) then
          read (proc,*) pagesTotal,pagesResident,pagesShared,pagesText,pagesData,pagesLibrary,pagesDirty
          if (ioStatus == 0) memoryUsageCode=pagesText*getPageSize()
       end if
       close(proc)
       ! Find available system memory. We first check if limits have been imposed by SLURM - if so, we use them, otherwise we use
       ! the system total memory.
       call Get_Environment_Variable('SLURM_MEM_PER_NODE',status=statusSlurmPerNode)
       call Get_Environment_Variable('SLURM_MEM_PER_CPU' ,status=statusSlurmPerCPU )
       if      (statusSlurmPerNode == 0) then
          call Get_Environment_Variable('SLURM_MEM_PER_NODE',value=slurmMemoryString,status=statusSlurmPerNode)
          if (statusSlurmPerNode /= 0) call Error_Report("failed to read environment variable 'SLURM_MEM_PER_NODE'"//{introspection:location})
          read (slurmMemoryString,*) memoryAvailable
          memoryAvailable=memoryAvailable*kilo**2
          sourceMemory="SLURM"
       else if (statusSlurmPerCPU  == 0) then
          call Get_Environment_Variable('SLURM_MEM_PER_CPU',value=slurmMemoryString,status=statusSlurmPerCPU)
          if (statusSlurmPerCPU    /= 0) call Error_Report("failed to read environment variable 'SLURM_MEM_PER_CPU'"//{introspection:location})
          read (slurmMemoryString,*) memoryAvailable
          memoryAvailable=memoryAvailable*kilo**2
          call Get_Environment_Variable('SLURM_CPUS_ON_NODE',value=slurmMemoryString,status=statusSlurmCountCPUs)
          if (statusSlurmCountCPUs /= 0) call Error_Report("failed to read environment variable 'SLURM_CPUS_ON_NODE'"//{introspection:location})
          read (slurmMemoryString,*) countCPUs
          memoryAvailable=memoryAvailable*countCPUs
          sourceMemory   ="SLURM"
       else
          memoryAvailable=getTotalSystemMemory()
          sourceMemory   ="system"
       end if
#ifdef USEMPI
       memoryAvailable=+        memoryAvailable   &
            &          /mpiSelf%countOnNode    ()
#endif
       ! Report on source of memory limits.
       call getSuffix(memoryAvailable,divisor,suffix)
       write (slurmMemoryString,'(f7.3,a1,a2)')  dble(memoryAvailable)/dble(divisor)," ",trim(adjustl(suffix))
       call displayMessage(displayGreen()//"NOTE: "//displayReset()//"memory available ("//trim(adjustl(slurmMemoryString))//") determined from "//trim(sourceMemory))
       ! Mark that these results are now initialized.
       codeMemoryUsageInitialized=.true.
    end if
    return
  end subroutine codeUsageGet
  
  !![
  <outputFileClose function="outputMemoryUsageMaximum"/>
  !!]
  subroutine outputMemoryUsageMaximum()
    !!{
    Output maximum memory usage information to the main output file.
    !!}
    use :: IO_HDF5    , only : hdf5Object
    use :: HDF5_Access, only : hdf5Access
    use :: Output_HDF5, only : outputFile
    implicit none
    type(hdf5Object) :: versionGroup

    !$ call hdf5Access%set()
    versionGroup=outputFile%openGroup('Version')
    call versionGroup%writeAttribute(memoryUsageMaximum,'memoryUsageMaximum')
    !$ call hdf5Access%unset()
    return
  end subroutine outputMemoryUsageMaximum

end module Memory_Reporting
