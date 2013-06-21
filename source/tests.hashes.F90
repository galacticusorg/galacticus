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

!% Contains a program to test features of the hashes (i.e. associative arrays) module.

program Test_Hashes
  !% Tests features of the hashes (i.e. associative arrays) module.
  use Unit_Tests
  use ISO_Varying_String
  use Memory_Management
  use Hashes
  use ISO_Varying_String
  implicit none
  type(integerScalarHash) :: myHash

  ! Read in basic code memory usage.
  call Code_Memory_Usage('tests.hashes.size')

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Hashes")

  ! Initialize the hash.
  call myHash%initialize()

  ! Create some entries in the hash.
  call myHash%set("dude"    , 34)
  call myHash%set("zuncular", 12)
  call myHash%set("munky"   ,-86)

  ! Assert the size of the hash.
  call Assert("hash size is correct",myHash%size(),3)

  ! Assert that the set entries exist.
  call Assert("hash entries exist",                                                      &
       &      [myHash%exists("dude") ,myHash%exists("munky"),myHash%exists("zuncular")], &
       &      [.true.                ,.true.                ,.true.                   ]  &
       &     )

  ! Assert that unset entries to not exist.
  call Assert("hash non-entries do not exist",                                           &
       &      [myHash%exists("buffy"),myHash%exists("giles"),myHash%exists("willow")  ], &
       &      [.false.               ,.false.               ,.false.                  ]  &
       &     )

  ! Assert values of set entries.
  call Assert("hash entries have correct values",                                        &
       &      [myHash%value ("dude") ,myHash%value ("munky"),myHash%value ("zuncular")], &
       &      [34                    ,-86                   ,12                       ]  &
       &     )

  ! Delete a hash entry.
  call myHash%delete("munky")

  ! Assert that deleted entry no longer exists.
  call Assert("deleted hash entries no longer exists",myHash%exists("munky"),.false.)

  ! Change the value of a hash entry.
  call myHash%set("dude",876)

  ! Assert that new value is correct.
  call Assert("changed hash entry has correct value",myHash%value("dude"),876)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Hashes
