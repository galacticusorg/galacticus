!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a program to test MPI functions.
!!}

program Test_MPI
  !!{
  Tests of MPI functions.
  !!}
  use   , intrinsic :: ISO_C_Binding, only : c_size_t
  use               :: Unit_Tests   , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group  , Unit_Tests_Finish
  use               :: Display      , only : displayMessage     , displayVerbositySet   , verbosityLevelStandard
  !$ use            :: OMP_Lib      , only : OMP_Get_Max_Threads
#ifdef USEMPI
  use               :: MPI          , only : MPI_Thread_Multiple
  use               :: MPI_Utilities, only : mpiBarrier         , mpiCounter            , mpiFinalize           , mpiInitialize    , &
          &                                  mpiSelf
  implicit none
  type   (mpiCounter) :: counter
  integer(c_size_t  ) :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  call mpiInitialize(MPI_Thread_Multiple)
  if (mpiSelf%rank() == 0) call Unit_Tests_Begin_Group("MPI")
  ! Test MPI/OpenMP counters.
  counter=mpiCounter()
  !$omp parallel
  i=counter%increment()
  !$omp end parallel
  call mpiBarrier()
  if (mpiSelf%rank() == 0)                                                        &
       & call Assert(                                                             &
       &             "MPI/OpenMP counter"                                       , &
       &                                       counter%get   ()                 , &
       &             int(OMP_Get_Max_Threads()*mpiSelf %count()-1,kind=c_size_t)  &
       &            )
  ! Finished tests.
  if (mpiSelf%rank() == 0) then
     call Unit_Tests_End_Group()
     call Unit_Tests_Finish()
  end if
  call mpiFinalize()
#else
  call displayMessage('SKIPPED: code was not compiled for MPI')
#endif
end program Test_MPI
