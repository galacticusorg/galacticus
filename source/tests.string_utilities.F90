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
Contains a program to test string handling utilities
!!}

program Test_String_Utilities
  !!{
  Tests that numerical range making code works correctly.
  !!}
  use :: Display           , only : displayVerbositySet      , verbosityLevelStandard, displayMessage
  use :: ISO_Varying_String, only : assignment(=)            , varying_string        , char                       , var_str           , &
       &                            operator(//)             , operator(==)
  use :: Kind_Numbers      , only : kind_int8
  use :: String_Handling   , only : Convert_VarString_To_Char, String_Count_Words    , String_Levenshtein_Distance, String_Lower_Case , &
       &                            String_Split_Words       , String_Upper_Case     , String_Upper_Case_First    , String_Superscript, &
       &                            String_Subscript         , operator(//)          , stringXMLFormat
  use :: System_Command    , only : System_Command_Do
  use :: Unit_Tests        , only : Assert                   , Unit_Tests_Begin_Group, Unit_Tests_End_Group       , Unit_Tests_Finish
  implicit none
  character(len=20        ), dimension(3) :: words
  type     (varying_string), dimension(3) :: myStrings
  type     (varying_string)               :: myString1, myString2      , &
       &                                     xmlFormat, xmlFormatTarget

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("String handling utilities")

  ! Test word counting.
  call Assert("count words: empty"                            ,String_Count_Words("  "                                                     ),0)
  call Assert("count words: 'one'"                            ,String_Count_Words("one"                                                    ),1)
  call Assert("count words: 'one two three'"                  ,String_Count_Words("one two three"                                          ),3)
  call Assert("count words: 'one two three' [null separated]" ,String_Count_Words("one"//char( 0)//"two"//char( 0)//"three"                ),3)
  call Assert("count words: 'one two three' [LF separated]"   ,String_Count_Words("one"//char(10)//"two"//char(10)//"three"                ),3)
  call Assert("count words: 'one two three' [CR separated]"   ,String_Count_Words("one"//char(13)//"two"//char(13)//"three"                ),3)
  call Assert("count words: 'one two three' [tab separated]"  ,String_Count_Words("one"//char( 9)//"two"//char( 9)//"three"                ),3)
  call Assert("count words: 'one two three' [comma separated]",String_Count_Words("one,two,three"                          ,separator ="," ),3)
  call Assert("count words: 'one two three' [bracketed]"      ,String_Count_Words("one{1 2 3} two{5 6 7} three{9 10 11}"   ,bracketing="{}"),3)

  ! Test word splitting.
  call String_Split_Words(words,"one"                                                    )
  call Assert("split words: 'one'"                                         ,words(1:1),["one"                                                   ])
  call String_Split_Words(words,"one two three"                                          )
  call Assert("split words: 'one two three'"                               ,words(1:3),["one  "           ,"two  "           ,"three"           ])
  call String_Split_Words(words,"one"//char( 0)//"two"//char( 0)//"three"                )
  call Assert("split words: 'one two three' [null separated]"              ,words(1:3),["one  "           ,"two  "           ,"three"           ])
  call String_Split_Words(words,"one"//char(10)//"two"//char(10)//"three"                )
  call Assert("split words: 'one two three' [LF separated]"                ,words(1:3),["one  "           ,"two  "           ,"three"           ])
  call String_Split_Words(words,"one"//char(13)//"two"//char(13)//"three"                )
  call Assert("split words: 'one two three' [CR separated]"                ,words(1:3),["one  "           ,"two  "           ,"three"           ])
  call String_Split_Words(words,"one"//char( 9)//"two"//char( 9)//"three"                )
  call Assert("split words: 'one two three' [tab separated]"               ,words(1:3),["one  "           ,"two  "           ,"three"           ])
  call String_Split_Words(words,"one,two,three"                          ,separator=","  )
  call Assert("split words: 'one two three' [comma separated]"             ,words(1:3),["one  "           ,"two  "           ,"three"           ])
  call String_Split_Words(words,"one (two and a half) three"             ,bracketing="()")
  call Assert("split words: 'one (two and a half) three' [with bracketing]",words(1:3),["one             ","(two and a half)","three           "])

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

  ! Conversion to superscript.
  call Assert('convert to superscript',char(String_Superscript("Result=(+10293-84756)?")),"Result⁼⁽⁺¹⁰²⁹³⁻⁸⁴⁷⁵⁶⁾?")

  ! Conversion to superscript.
  call Assert('convert to subscript'  ,char(String_Subscript  ("Result=(+10293-84756)?")),"Result₌₍₊₁₀₂₉₃₋₈₄₇₅₆₎?")
  
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

  ! Pretty-print an XML string.
  xmlFormat      =stringXMLFormat('<make><this xml="look"/>**B<pretty value="please" and="thank you"/></make>',forceColor=.true.)
  xmlFormatTarget='<[34mmake[0m>[0m'//char(10)//'  <[34mthis[0m [33mxml[0m="[32mlook[0m"/>[0m'//char(10)//'[1m  <[34mpretty[0m[1m [33mvalue[0m[1m="[32mplease[0m[1m" [33mand[0m[1m="[32mthank you[0m[1m"/>[0m'//char(10)//'</[34mmake[0m>[0m'
  call Assert('XML formatter string #1',xmlFormat,xmlFormatTarget)
  call displayMessage(var_str('XML formatted string #1:')//char(10)//xmlFormat)

  xmlFormat      =stringXMLFormat('<mergerTreeEvolver value="doop">**B<allTreesExistAtFinalTime value="false" />**C</mergerTreeEvolver>',indentInitial=6,forceColor=.true.)
  xmlFormatTarget='      <[34mmergerTreeEvolver[0m [33mvalue[0m="[32mdoop[0m">[0m'//char(10)//'[1m        <[34mallTreesExistAtFinalTime[0m[1m [33mvalue[0m[1m="[32mfalse[0m[1m" />[0m'//char(10)//'        ......'//char(10)//'      </[34mmergerTreeEvolver[0m>[0m'
  call Assert('XML formatter string #2',xmlFormat,xmlFormatTarget)
  call displayMessage(var_str('XML formatted string #2:')//char(10)//xmlFormat)
  
  ! block
  !   ! This is an example of how to get the target string for these assertions. Once the formatter is producing the correct output,
  !   ! it can be written to file, and then included into this source code for comparison.
  !   integer :: i
  !   open(newUnit=i,file='xmlFormatted.txt',form='formatted',status='unknown')
  !   write (i,'(a)') char(xmlFormat)
  !   close(i)
  !   call System_Command_Do("sed -i~ -r ':a;N;$!ba;s/\n/'\''\/\/char\(10\)\/\/'\''/g' xmlFormatted.txt")
  ! end block

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_String_Utilities
