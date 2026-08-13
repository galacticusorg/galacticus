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

!+    Contributions to this file made by: Andrew Benson. This test was drafted with assistance from Claude, and reviewed and
!+    verified by Andrew Benson.

!!{RST
Contains a program to test the accumulation of warning messages.
!!}

program Test_Warnings
  !!{RST
  Tests that the cost of issuing a warning does not grow with the number of warnings already
  issued. Warnings are recorded (so that they can be reviewed if the run subsequently fails) in a
  list which is maintained inside a global critical section. If the cost of adding to that list
  scales with the length of the list then a code path which warns frequently (for example, once per
  node) becomes quadratically expensive and, since the work is serialized, reduces a multi-threaded
  run to a single-threaded one.
  !!}
  use :: Display   , only : displayVerbositySet, verbosityLevelStandard
  use :: Error     , only : Warn
  use :: Unit_Tests, only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish, &
       &                    compareLessThan
  implicit none
  ! Number of warnings issued in each of two successive batches, and the maximum by which the second
  ! batch is permitted to be slower than the first. If the cost of recording a warning were to scale
  ! with the number already recorded the second batch would be approximately three times slower than
  ! the first.
  integer                     , parameter :: countBatch     =50000
  double precision            , parameter :: slowdownMaximum=2.0d0
  integer                                 :: i
  double precision                        :: timeBegin            , timeMiddle, &
       &                                     timeEnd
  character       (len=32    )            :: message

  ! Set the verbosity level below `warn`, so that warnings are recorded (for possible review on an
  ! error condition) rather than simply displayed.
  call displayVerbositySet    (verbosityLevelStandard)
  call Unit_Tests_Begin_Group ("Warnings"            )
  ! Issue two successive batches of warnings, all of which are unique - the worst case for the
  ! accumulated list of warnings - and time each batch.
  call CPU_Time(timeBegin )
  do i=1,countBatch
     write (message,'(a,i8)') "test warning ",i
     call Warn(message)
  end do
  call CPU_Time(timeMiddle)
  do i=countBatch+1,2*countBatch
     write (message,'(a,i8)') "test warning ",i
     call Warn(message)
  end do
  call CPU_Time(timeEnd   )
  call Assert              ("cost of issuing a warning is independent of the number already issued",(timeEnd-timeMiddle)/(timeMiddle-timeBegin),slowdownMaximum,compareLessThan)
  call Unit_Tests_End_Group(                                                                                                                                                   )
  call Unit_Tests_Finish   (                                                                                                                                                   )
end program Test_Warnings
