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

!!{RST
Contains a program to test MPI functions.
!!}

program Test_MPI
  !!{RST
  Tests of MPI functions.
  !!}
  use   , intrinsic :: ISO_C_Binding     , only : c_size_t
  use               :: Unit_Tests        , only : Assert               , Unit_Tests_Begin_Group, Unit_Tests_End_Group  , Unit_Tests_Finish
  use               :: Display           , only : displayMessage       , displayVerbositySet   , verbosityLevelStandard
  !$ use            :: OMP_Lib           , only : OMP_Get_Max_Threads
#ifdef USEMPI
  use               :: Events_Hooks      , only : eventsHooksInitialize
  use               :: MPI_F08           , only : MPI_Thread_Multiple
  use               :: MPI_Utilities     , only : mpiBarrier           , mpiCounter            , mpiFinalize           , mpiInitialize    , &
       &                                          mpiSelf
  use               :: ISO_Varying_String, only : char
#endif
  implicit none
#ifdef USEMPI
  type            (mpiCounter)                                :: counter
  integer         (c_size_t  )                                :: i               , counterSum     , &
       &                                                         counterSumOver  , counters       , &
       &                                                         counterMaxOver  , j
  integer                                                     :: rankOriginal    , countOriginal
  logical                     , allocatable, dimension(:    ) :: maskExcludeFirst, maskExcludeLast
  double precision            , allocatable, dimension(:,:  ) :: dataDouble1D
  double precision            , allocatable, dimension(:,:,:) :: dataDouble2D
  integer                     , allocatable, dimension(:,:  ) :: dataInteger1D
  logical                     , allocatable, dimension(:,:  ) :: dataLogical1D
  integer                                                     :: rankNext
#endif
  
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

#ifdef USEMPI
  call mpiInitialize        (MPI_Thread_Multiple)
  call eventsHooksInitialize(                   )
  call Unit_Tests_Begin_Group("MPI: rank "//char(mpiSelf%rankLabel()))
  ! Test MPI/OpenMP counters.
  counters  =int(OMP_Get_Max_Threads()*mpiSelf%count(),kind=c_size_t)
  counter   =mpiCounter()
  counterSum=0_c_size_t
  j         =0_c_size_t
  call counter%reset()
  !$omp parallel private(i) reduction(max : j)
  do while (.true.)
     i=counter%increment()
     if (i >= counters) exit
     j=i
     !$omp atomic
     counterSum=counterSum+i
  end do
  !$omp end parallel
  call mpiBarrier()
  counterMaxOver=mpiSelf%maxval(         j)
  counterSumOver=mpiSelf%sum   (counterSum)
  if (mpiSelf%rank() == 0) then
     call Unit_Tests_Begin_Group("Counters")
     call Assert(                            &
          &      "MPI/OpenMP counter total", &
          &       counterMaxOver           , &
          &       counters-1                 &
          &     )
     call Assert(                            &
          &      "MPI/OpenMP counter sum"  , &
          &       counterSumOver           , &
          &       counters*(counters-1)/2    &
          &     )
     call Unit_Tests_End_Group()
  end if
  ! Test masked reductions. Each process contributes (rank+1), and the mask excludes either the first or the last process, so
  ! the minimum over the included processes is 2, and the maximum is the process count less one. For the logical reductions,
  ! only the excluded first process differs from the others, so it must not affect the result.
  if (mpiSelf%count() > 1) then
     allocate(maskExcludeFirst(0:mpiSelf%count()-1))
     allocate(maskExcludeLast (0:mpiSelf%count()-1))
     maskExcludeFirst                   =.true.
     maskExcludeFirst(0                )=.false.
     maskExcludeLast                    =.true.
     maskExcludeLast (mpiSelf%count()-1)=.false.
     call Unit_Tests_Begin_Group("Masked reductions")
     call Assert("minval (double array)"  ,mpiSelf%minval([dble(mpiSelf%rank()+1)              ],maskExcludeFirst),[2.0d0                           ])
     call Assert("minval (integer array)" ,mpiSelf%minval([     mpiSelf%rank()+1               ],maskExcludeFirst),[2                               ])
     call Assert("minval (double scalar)" ,mpiSelf%minval( dble(mpiSelf%rank()+1)               ,maskExcludeFirst), 2.0d0                            )
     call Assert("minval (integer scalar)",mpiSelf%minval(      mpiSelf%rank()+1                ,maskExcludeFirst), 2                                )
     call Assert("maxval (double array)"  ,mpiSelf%maxval([dble(mpiSelf%rank()+1)              ],maskExcludeLast ),[dble(mpiSelf%count()-1         )])
     call Assert("maxval (size_t array)"  ,mpiSelf%maxval([int (mpiSelf%rank()+1,kind=c_size_t)],maskExcludeLast ),[int (mpiSelf%count()-1,c_size_t)])
     call Assert("minloc (double array)"  ,mpiSelf%minloc([dble(mpiSelf%rank()+1)              ],maskExcludeFirst),[1                               ])
     call Assert("maxloc (double array)"  ,mpiSelf%maxloc([dble(mpiSelf%rank()+1)              ],maskExcludeLast ),[mpiSelf%count()-2               ])
     call Assert("all (logical scalar)"   ,mpiSelf%all   (      mpiSelf%rank() /= 0             ,maskExcludeFirst),.true.                            )
     call Assert("any (logical scalar)"   ,mpiSelf%any   (      mpiSelf%rank() == 0             ,maskExcludeFirst),.false.                           )
     call Unit_Tests_End_Group()
     deallocate(maskExcludeFirst,maskExcludeLast)
     ! Test data requests. Each process requests data from the next process (cyclically), which returns data identifying itself.
     rankNext     =modulo(mpiSelf%rank()+1,mpiSelf%count())
     dataDouble1D =mpiSelf%requestData([rankNext],[dble(mpiSelf%rank()),-dble(mpiSelf%rank())])
     dataDouble2D =mpiSelf%requestData([rankNext],reshape([dble(mpiSelf%rank()),1.0d0,2.0d0,-dble(mpiSelf%rank())],[2,2]))
     dataInteger1D=mpiSelf%requestData([rankNext],[mpiSelf%rank(),2*mpiSelf%rank()])
     dataLogical1D=mpiSelf%requestData([rankNext],[modulo(mpiSelf%rank(),2) == 0,.true.])
     call Unit_Tests_Begin_Group("Data requests")
     call Assert("requestData (double 1-D)" ,dataDouble1D (:,  1),[dble(rankNext),-dble(rankNext)]                           )
     call Assert("requestData (double 2-D)" ,dataDouble2D (:,:,1),reshape([dble(rankNext),1.0d0,2.0d0,-dble(rankNext)],[2,2]))
     call Assert("requestData (integer 1-D)",dataInteger1D(:,  1),[rankNext,2*rankNext]                                      )
     call Assert("requestData (logical 1-D)",dataLogical1D(:,  1),[modulo(rankNext,2) == 0,.true.]                           )
     call Unit_Tests_End_Group()
  end if
  ! Test communicators.
  call Unit_Tests_Begin_Group("Communicators")
  rankOriginal =mpiSelf%rank ()
  countOriginal=mpiSelf%count()
  call mpiSelf%communicatorPush(color=mpiSelf%rank())
  call Assert("Rank updated"  ,mpiSelf%rank (),            0)
  call Assert("Count updated" ,mpiSelf%count(),            1)
  call mpiSelf%communicatorPop()
  call Assert("Rank restored" ,mpiSelf%rank (), rankOriginal)
  call Assert("Count restored",mpiSelf%count(),countOriginal)
  call Unit_Tests_End_Group()
  ! Finished tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
  call mpiFinalize()
#else
  call displayMessage('SKIPPED: code was not compiled for MPI')
#endif
end program Test_MPI
