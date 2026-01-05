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
Contains a module which implements unit testing.
!!}

module Unit_Tests
  !!{
  Implements unit testing.
  !!}
  use :: ISO_Varying_String, only : varying_string
  implicit none

  !![
  <generic identifier="Type">
   <instance label="integerScalar"       intrinsic="integer"              rank=""                       format="i8"    reduction="regEx¦(.*)¦$1¦"     />
   <instance label="integerRank1"        intrinsic="integer"              rank=", dimension(:)"         format="i8"    reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="integerLongScalar"   intrinsic="integer(kind_int8)"   rank=""                       format="i16"   reduction="regEx¦(.*)¦$1¦"     />
   <instance label="integerLongRank1"    intrinsic="integer(kind_int8)"   rank=", dimension(:)"         format="i16"   reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="logicalScalar"       intrinsic="logical"              rank=""                       format="l1"    reduction="regEx¦(.*)¦$1¦"     />
   <instance label="logicalRank1"        intrinsic="logical"              rank=", dimension(:)"         format="l1"    reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="characterScalar"     intrinsic="character(len=*)"     rank=""                       format=""      reduction="regEx¦(.*)¦$1¦"     />
   <instance label="characterRank1"      intrinsic="character(len=*)"     rank=", dimension(:)"         format=""      reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="varyingStringScalar" intrinsic="type(varying_string)" rank=""                       format=""      reduction="regEx¦(.*)¦$1¦"     />
   <instance label="varyingStringRank1"  intrinsic="type(varying_string)" rank=", dimension(:)"         format=""      reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="doubleScalar"        intrinsic="double precision"     rank=""                       format="e16.8" reduction="regEx¦(.*)¦$1¦"     />
   <instance label="doubleRank1"         intrinsic="double precision"     rank=", dimension(:)"         format="e16.8" reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="doubleRank2"         intrinsic="double precision"     rank=", dimension(:,:)"       format="e16.8" reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="doubleRank3"         intrinsic="double precision"     rank=", dimension(:,:,:)"     format="e16.8" reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="doubleRank4"         intrinsic="double precision"     rank=", dimension(:,:,:,:)"   format="e16.8" reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="doubleRank5"         intrinsic="double precision"     rank=", dimension(:,:,:,:,:)" format="e16.8" reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="realRank1"           intrinsic="real"                 rank=", dimension(:)"         format="e16.8" reduction="regEx¦(.*)¦all($1)¦"/>
   <instance label="doubleComplexRank1"  intrinsic="double complex"       rank=", dimension(:)"         format="e16.8" reduction="regEx¦(.*)¦all($1)¦"/>
  </generic>
  !!]

  private
  public :: Assert, Skip, Unit_Tests_Finish, Unit_Tests_Begin_Group, Unit_Tests_End_Group

  ! Types of comparison.
  integer, parameter, public :: compareEquals            =0
  integer, parameter, public :: compareNotEqual          =1
  integer, parameter, public :: compareLessThan          =2
  integer, parameter, public :: compareGreaterThan       =3
  integer, parameter, public :: compareLessThanOrEqual   =4
  integer, parameter, public :: compareGreaterThanOrEqual=5

  !![
  <enumeration>
   <name>test</name>
   <description>Statuses for unit tests.</description>
   <entry label="passed"/>
   <entry label="failed" />
   <entry label="skipped"/>
  </enumeration>
  !!]

  ! Type for assert results.
  type assertResult
     !!{
     A derived type for storing results of asserts.
     !!}
     type    (enumerationTestType)          :: result
     logical                                :: beginGroup          , endGroup
     type   (varying_string      )          :: label               , note
     type   (assertResult        ), pointer :: nextResult => null()
  end type assertResult

  ! Results list.
  type   (assertResult), target  :: firstResult
  type   (assertResult), pointer :: currentResult
  logical                        :: firstAssert  =.true.

  ! Interface for assert routines.
  interface Assert
     !!{
     Generic interface for assert routines.
     !!}
     module procedure Assert{Type¦label}
  end interface Assert

contains

  function getStatus(passed)
    !!{
    Return the status code for a test on the basis of a boolean pass/fail.
    !!}
    implicit none
    type   (enumerationTestType)                :: getStatus
    logical                     , intent(in   ) :: passed

    if (passed) then
       getStatus=testPassed
    else
       getStatus=testFailed
    end if
    return
  end function getStatus

  subroutine Skip(testName,reason)
    !!{
    Record that a test was skipped.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*       ), intent(in   ) :: testName  , reason
    type     (assertResult), pointer       :: result

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()
    result%result=testSkipped
    result%label =trim(testName)
    result%note  =trim(reason  )
    return
  end subroutine Skip

  subroutine Assert{Type¦label}(testName,value1,value2,compare{Type¦match¦^(double|real)¦,absTol,relTol¦})
    !!{
    Assess and record an assertion.
    !!}
    use                                  :: Error               , only : Error_Report
    use                                  :: Numerical_Comparison, only : Values_Agree
    {Type¦match¦^integerLong¦        use :: Kind_Numbers        , only : kind_int8¦}
    {Type¦match¦^varyingString¦      use :: ISO_Varying_String  , only : operator(==) , operator(/=), operator(<)  , operator(>), operator(>=), operator(<=), assignment(=), operator(//), var_str¦}
    {Type¦match¦^(?!(varyingString))¦use :: ISO_Varying_String  , only : assignment(=), operator(//), var_str¦}
    character                                  (len=*         ), intent(in   )              :: testName
    {Type¦intrinsic}                                           , intent(in   )  {Type¦rank} :: value1    , value2
    integer                                                    , intent(in   ), optional    :: compare
    {Type¦match¦^(double|real)¦{Type¦intrinsic}                , intent(in   ), optional    :: absTol    , relTol¦}
    type                                       (assertResult  ), pointer                    :: result
    logical                                                    , allocatable    {Type¦rank} :: passed
    type                                       (varying_string)                             :: comparison
    character                                  (len=256       )                             :: label
    character                                  (len= 2        )                             :: separator    
    !![
    <optionalArgument name="compare" defaultsTo="compareEquals"/>
    !!]

    ! Allocate array for results.
    !![
    <allocate variable="passed" size="value1"/>
    !!]
    ! Perform the comparison.
    select case (compare_)
    case (compareEquals)
       {Type¦match¦^(double|real)¦passed=Values_Agree(value1,value2,absTol,relTol)¦}
       {Type¦match¦^logical¦passed=value1 .eqv. value2¦}
       {Type¦match¦^(?!(double|real|logical))¦passed=value1 == value2¦}
       comparison="≉"
    case (compareNotEqual          )
       passed={Type¦match¦^logical¦value1 .neqv. value2¦value1 /= value2}
       comparison="="
   case (compareLessThan          )
       passed={Type¦match¦^(logical|doubleComplex)¦.false.¦value1  < value2}
       comparison="≮"
       {Type¦match¦^(logical|doubleComplex)¦call Error_Report('unsupported comparison'//{introspection:location})¦}
    case (compareGreaterThan       )
       passed={Type¦match¦^(logical|doubleComplex)¦.false.¦value1  > value2}
       comparison="≯"
       {Type¦match¦^(logical|doubleComplex)¦call Error_Report('unsupported comparison'//{introspection:location})¦}
    case (compareLessThanOrEqual   )
       passed={Type¦match¦^(logical|doubleComplex)¦.false.¦value1 <= value2}
       comparison="≰"
       {Type¦match¦^(logical|doubleComplex)¦call Error_Report('unsupported comparison'//{introspection:location})¦}
    case (compareGreaterThanOrEqual)
       passed={Type¦match¦^(logical|doubleComplex)¦.false.¦value1 >= value2}
       comparison="≱"
       {Type¦match¦^(logical|doubleComplex)¦call Error_Report('unsupported comparison'//{introspection:location})¦}
    case default
       passed=.false.
       comparison=""
       call Error_Report('unknown comparison'//{introspection:location})
    end select
    ! Get an object to store the results in.
    result => Get_New_Assert_Result()
    ! Store the result.
    result%result=getStatus({Type¦reduction¦passed})
    result%label =trim(testName)
    if (.not.{Type¦reduction¦passed}) then
       result%note=var_str('')
       separator  =''
       !![
       <forEach variable="passed">
	 if (.not.passed{index}) then
	 if (separator /= "") result%note=result%note//separator
         separator='; '
         write (label,%index%) {{index}}; result%note=result%note//trim(adjustl(label))//': '
         {Type¦match¦^varyingString¦result%note=result%note//value1{index}¦}
         {Type¦match¦^character¦result%note=result%note//trim(value1{index})¦}
         {Type¦match¦^(?!(character|varyingString))¦write (label,'({Type¦format})') value1{index}; result%note=result%note//trim(adjustl(label))¦}
         result%note=result%note//" "//comparison//" "
         {Type¦match¦^varyingString¦result%note=result%note//value2{index}¦}
         {Type¦match¦^character¦result%note=result%note//trim(value2{index})¦}
         {Type¦match¦^(?!(character|varyingString))¦write (label,'({Type¦format})') value2{index}; result%note=result%note//trim(adjustl(label))¦}
        end if
       </forEach>
       !!]
     end if
    return
  end subroutine Assert{Type¦label}

  subroutine Unit_Tests_Begin_Group(groupName)
    !!{
    Marks that a unit test group has begun.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*       ), intent(in   ) :: groupName
    type     (assertResult), pointer       :: result

    result => Get_New_Assert_Result()
    result%beginGroup=.true.
    result%label     =trim(groupName)
    return
  end subroutine Unit_Tests_Begin_Group

  subroutine Unit_Tests_End_Group
    !!{
    Marks that a unit test group has ended.
    !!}
    implicit none
    type(assertResult), pointer :: result

    result => Get_New_Assert_Result()
    result%endGroup=.true.
    return
  end subroutine Unit_Tests_End_Group

  subroutine Unit_Tests_Finish
    !!{
    Write out the results of unit testing.
    !!}
    use :: Display           , only : displayIndent, displayMessage, displayUnindent
    use :: ISO_Varying_String, only : assignment(=), operator(//)  , operator(/=)
    use :: String_Handling   , only : operator(//)
    implicit none
    type   (assertResult  ), pointer :: nextResult, result
    integer                          :: failCount , passCount , percentage
    type   (varying_string)          :: message

    passCount=0
    failCount=0
    if (.not.firstAssert) then
       result => firstResult
       do while (associated(result))
          if (.not.(result%beginGroup.or.result%endGroup)) then
             select case (result%result%ID)
             case (testFailed%ID)
                failCount=failCount+1
                message=" FAILED: "//result%label
                if (result%note /= "") message=message//" {"//result%note//"}"
             case (testPassed%ID)
                passCount=passCount+1
                message=" passed: "//result%label
             case (testSkipped%ID)
                message="skipped: "//result%label//" {"//result%note//"}"
             end select
             call displayMessage(message)
          else
             if (result%beginGroup) call displayIndent  (result%label)
             if (result%endGroup  ) call displayUnindent("")
          end if
          result => result%nextResult
       end do
    end if
    if (passCount+failCount > 0) then
       percentage=int(100.0d0*dble(passCount)/dble(passCount+failCount)+0.5d0)
    else
       percentage=100
    end if
    message="Tests passed: "
    message=message//passCount//" ("//percentage//"%)"
    call displayMessage(message)
    if (passCount+failCount > 0) then
       percentage=int(100.0d0*dble(failCount)/dble(passCount+failCount)+0.5d0)
    else
       percentage=  0
    end if
    message="Tests failed: "
    message=message//failCount//" ("//percentage//"%)"
    call displayMessage(message)

    ! Destroy the list of results.
    result => firstResult%nextResult
    do while (associated(result))
       nextResult => result%nextResult
       call result%label%destroy()
       deallocate(result)
       result => nextResult
    end do
    return
  end subroutine Unit_Tests_Finish

  function Get_New_Assert_Result() result(newResult)
    !!{
    Get a new assert result object.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type(assertResult), pointer :: newResult

    ! Return the first result if this is the first assert, otherwise, allocate a new one.
    select case (firstAssert)
    case (.false.)
       allocate(currentResult%nextResult)
       newResult => currentResult%nextResult
    case (.true. )
       newResult => firstResult
       firstAssert=.false.
    end select
    currentResult => newResult
    newResult%beginGroup=.false.
    newResult%endGroup  =.false.
    newResult%result    =testFailed
    newResult%note      =""
    newResult%nextResult => null()
    return
  end function Get_New_Assert_Result

end module Unit_Tests
