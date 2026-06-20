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
Contains a program to test features of the dictionaries (i.e. associative arrays) module.
!!}

program Test_Dictionaries
  !!{RST
  Tests features of the dictionaries (i.e. associative arrays) module.
  !!}
  use :: Display         , only : displayVerbositySet, verbosityLevelStandard
  use :: Dictionaries    , only : genericDictionary  , integerDictionary
  use :: Input_Parameters, only : inputParameter     , inputParameters
  use :: Unit_Tests      , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type (integerDictionary)              :: myDictionary
  type (genericDictionary)              :: genericDictionary_
  class(inputParameters  ), allocatable :: inputParameters_
  type (inputParameter   )              :: inputParameter_
  class(*                ), pointer     :: generic_

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Dictionaries")
  ! Tests of integer scalar dictionaries.
  !! Initialize the dictionary.
  call myDictionary%initialize()
  !! Create some entries in the dictionary. These are chosen to ensure that the dictionary must be re-ordered as they entries are added - this
  !! therefore tests that the re-ordering works correctly.
  call myDictionary%set("node092", 34)
  call myDictionary%set("node078", 12)
  call myDictionary%set("node077",-86)
  call Assert("dictionary entries exist",                                                         &
       &      [myDictionary%exists("node092"),myDictionary%exists("node078"),myDictionary%exists("node077")], &
       &      [.true.                  ,.true.                  ,.true.                  ]  &
       &     )
  !! Create some more entries in the dictionary. 
  call myDictionary%set("dude"    , 34)
  call myDictionary%set("zuncular", 12)
  call myDictionary%set("munky"   ,-86)
  !! Assert the size of the dictionary.
  call Assert("dictionary size is correct",myDictionary%size(),6)
  !! Assert that the set entries exist.
  call Assert("dictionary entries exist",                                                      &
       &      [myDictionary%exists("dude") ,myDictionary%exists("munky"),myDictionary%exists("zuncular")], &
       &      [.true.                ,.true.                ,.true.                   ]  &
       &     )
  !! Assert that unset entries to not exist.
  call Assert("dictionary non-entries do not exist",                                           &
       &      [myDictionary%exists("buffy"),myDictionary%exists("giles"),myDictionary%exists("willow")  ], &
       &      [.false.               ,.false.               ,.false.                  ]  &
       &     )
  !! Assert values of set entries.
  call Assert("dictionary entries have correct values",                                        &
       &      [myDictionary%value ("dude") ,myDictionary%value ("munky"),myDictionary%value ("zuncular")], &
       &      [34                    ,-86                   ,12                       ]  &
       &     )
  !! Delete a dictionary entry.
  call myDictionary%delete("munky")
  !! Assert that deleted entry no longer exists.
  call Assert("deleted dictionary entries no longer exists",myDictionary%exists("munky"),.false.)
  !! Change the value of a dictionary entry.
  call myDictionary%set("dude",876)
  !! Assert that new value is correct.
  call Assert("changed dictionary entry has correct value",myDictionary%value("dude"),876)
  !! Tests of generic scalar dictionaries.
  !! Initialize the dictionary.
  call genericDictionary_%initialize()
  !! Create some entries in the dictionary.
  allocate(inputParameters :: inputParameters_)
  call genericDictionary_%set("class",inputParameters_)
  call genericDictionary_%set("type" ,inputParameter_ )
  !! Assert types of entries.
  generic_ => genericDictionary_%value("class")
  call Assert("correct class for polymoprhic object"    ,same_type_as(generic_,inputParameters_),.true.)
  generic_ => genericDictionary_%value("type" )
  call Assert("correct class for non-polymoprhic object",same_type_as(generic_,inputParameter_ ),.true.)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Dictionaries
