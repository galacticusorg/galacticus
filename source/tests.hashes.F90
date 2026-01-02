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
Contains a program to test features of the hashes (i.e. associative arrays) module.
!!}

program Test_Hashes
  !!{
  Tests features of the hashes (i.e. associative arrays) module.
  !!}
  use :: Display         , only : displayVerbositySet, verbosityLevelStandard
  use :: Hashes          , only : genericHash        , integerHash
  use :: Input_Parameters, only : inputParameter     , inputParameters
  use :: Unit_Tests      , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type (integerHash    )              :: myHash
  type (genericHash    )              :: genericHash_
  class(inputParameters), allocatable :: inputParameters_
  type (inputParameter )              :: inputParameter_
  class(*              ), pointer     :: generic_

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Hashes")
  ! Tests of integer scalar hashes.
  !! Initialize the hash.
  call myHash%initialize()
  !! Create some entries in the hash. These are chosen to ensure that the hash must be re-ordered as they entries are added - this
  !! therefore tests that the re-ordering works correctly.
  call myHash%set("node092", 34)
  call myHash%set("node078", 12)
  call myHash%set("node077",-86)
  call Assert("hash entries exist",                                                         &
       &      [myHash%exists("node092"),myHash%exists("node078"),myHash%exists("node077")], &
       &      [.true.                  ,.true.                  ,.true.                  ]  &
       &     )
  !! Create some more entries in the hash. 
  call myHash%set("dude"    , 34)
  call myHash%set("zuncular", 12)
  call myHash%set("munky"   ,-86)
  !! Assert the size of the hash.
  call Assert("hash size is correct",myHash%size(),6)
  !! Assert that the set entries exist.
  call Assert("hash entries exist",                                                      &
       &      [myHash%exists("dude") ,myHash%exists("munky"),myHash%exists("zuncular")], &
       &      [.true.                ,.true.                ,.true.                   ]  &
       &     )
  !! Assert that unset entries to not exist.
  call Assert("hash non-entries do not exist",                                           &
       &      [myHash%exists("buffy"),myHash%exists("giles"),myHash%exists("willow")  ], &
       &      [.false.               ,.false.               ,.false.                  ]  &
       &     )
  !! Assert values of set entries.
  call Assert("hash entries have correct values",                                        &
       &      [myHash%value ("dude") ,myHash%value ("munky"),myHash%value ("zuncular")], &
       &      [34                    ,-86                   ,12                       ]  &
       &     )
  !! Delete a hash entry.
  call myHash%delete("munky")
  !! Assert that deleted entry no longer exists.
  call Assert("deleted hash entries no longer exists",myHash%exists("munky"),.false.)
  !! Change the value of a hash entry.
  call myHash%set("dude",876)
  !! Assert that new value is correct.
  call Assert("changed hash entry has correct value",myHash%value("dude"),876)
  !! Tests of generic scalar hashes.
  !! Initialize the hash.
  call genericHash_%initialize()
  !! Create some entries in the hash.
  allocate(inputParameters :: inputParameters_)
  call genericHash_%set("class",inputParameters_)
  call genericHash_%set("type" ,inputParameter_ )
  !! Assert types of entries.
  generic_ => genericHash_%value("class")
  call Assert("correct class for polymoprhic object"    ,same_type_as(generic_,inputParameters_),.true.)
  generic_ => genericHash_%value("type" )
  call Assert("correct class for non-polymoprhic object",same_type_as(generic_,inputParameter_ ),.true.)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Hashes
