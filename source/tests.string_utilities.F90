!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program to test string handling utilities

program Test_String_Utilities
  !% Tests that numerical range making code works correctly.
  use Unit_Tests
  use ISO_Varying_String
  use String_Handling
  use Kind_Numbers
  implicit none
  character(len=10),    dimension(3) :: words
  type(varying_string), dimension(3) :: myStrings
  type(varying_string)               :: myString1,myString2

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("String handling utilities")
  
  ! Test word counting.
  call Assert("count words: empty"                            ,String_Count_Words("  "                                                   ),0)
  call Assert("count words: 'one'"                            ,String_Count_Words("one"                                                  ),1)
  call Assert("count words: 'one two three'"                  ,String_Count_Words("one two three"                                        ),3)
  call Assert("count words: 'one two three' [null separated]" ,String_Count_Words("one"//char( 0)//"two"//char( 0)//"three"              ),3)
  call Assert("count words: 'one two three' [LF separated]"   ,String_Count_Words("one"//char(10)//"two"//char(10)//"three"              ),3)
  call Assert("count words: 'one two three' [CR separated]"   ,String_Count_Words("one"//char(13)//"two"//char(13)//"three"              ),3)
  call Assert("count words: 'one two three' [tab separated]"  ,String_Count_Words("one"//char( 9)//"two"//char( 9)//"three"              ),3)
  call Assert("count words: 'one two three' [comma separated]",String_Count_Words("one,two,three"                          ,separator=","),3)

  ! Test word splitting.
  call String_Split_Words(words,"one"                                                  )
  call Assert("split words: 'one'"                            ,words(1:1),["one"                  ])
  call String_Split_Words(words,"one two three"                                        )
  call Assert("split words: 'one two three'"                  ,words(1:3),["one  ","two  ","three"])
  call String_Split_Words(words,"one"//char( 0)//"two"//char( 0)//"three"              )
  call Assert("split words: 'one two three' [null separated]" ,words(1:3),["one  ","two  ","three"])
  call String_Split_Words(words,"one"//char(10)//"two"//char(10)//"three"              )
  call Assert("split words: 'one two three' [null separated]" ,words(1:3),["one  ","two  ","three"])
  call String_Split_Words(words,"one"//char(13)//"two"//char(13)//"three"              )
  call Assert("split words: 'one two three' [null separated]" ,words(1:3),["one  ","two  ","three"])
  call String_Split_Words(words,"one"//char( 9)//"two"//char( 9)//"three"              )
  call Assert("split words: 'one two three' [null separated]" ,words(1:3),["one  ","two  ","three"])
  call String_Split_Words(words,"one,two,three"                          ,separator=",")
  call Assert("split words: 'one two three' [comma separated]",words(1:3),["one  ","two  ","three"])

  ! Test concatenation of varying strings.
  myString1='test'
  myString2='test1234'
  call Assert('concatenate varying string and integer',myString1//1234          ,myString2)
  call Assert('concatenate varying string and integer',myString1//1234_kind_int8,myString2)

  ! Test case changing.
  call Assert('string to upper case',String_Upper_Case("ConVerT mE to ALL uPpeR Case 8204!&#$"),"CONVERT ME TO ALL UPPER CASE 8204!&#$")
  call Assert('string to lower case',String_Lower_Case("ConVerT mE to ALL LoWeR Case 8204!&#$"),"convert me to all lower case 8204!&#$")
  call Assert('first lower case to upper case',String_Upper_Case_First("convert first character to upper case"),"Convert first character to upper case")

  ! Conversion of varying string to character array.
  myStrings=["one  ","two  ","three"]
  call Assert('convert varying string to character array',Convert_VarString_To_Char(myStrings),["one  ","two  ","three"])

  ! Test Levenshtein distance.
  call Assert(                                                  &
       &      'measure Levenshtein distance'                  , &
       &      [                                                 &
       &       String_Levenshtein_Distance('kitten','sitting'), &
       &       String_Levenshtein_Distance('grapes','wine'   ), &
       &       String_Levenshtein_Distance('monkey','human'  ), &
       &       String_Levenshtein_Distance('lead'  ,'gold'   )  &
       &      ]                                               , &
       &      [                                                 &
       &       3                                              , &
       &       5                                              , &
       &       6                                              , &
       &       3                                                &
       &      ]                                                 &
       &     )
  
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_String_Utilities
