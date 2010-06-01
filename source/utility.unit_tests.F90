!% Contains a module which implements unit testing.

module Unit_Tests
  !% Implements unit testing.
  use ISO_Varying_String
  private
  public :: Assert, Unit_Tests_Finish, Unit_Tests_Begin_Group, Unit_Tests_End_Group

  ! Types of comparison.
  integer, parameter, public :: compareEquals            =0
  integer, parameter, public :: compareNotEqual          =1
  integer, parameter, public :: compareLessThan          =2
  integer, parameter, public :: compareGreaterThan       =3
  integer, parameter, public :: compareLessThanOrEqual   =4
  integer, parameter, public :: compareGreaterThanOrEqual=5

  ! Type for assert results.
  type assertResult
     !% A derived type for storing results of asserts.
     logical                      :: beginGroup,endGroup,result
     type(varying_string)         :: label
     type(assertResult),  pointer :: nextResult
  end type assertResult

  ! Results list.
  type(assertResult), target  :: firstResult
  type(assertResult), pointer :: currentResult
  logical                     :: firstAssert=.true.

contains

  subroutine Assert(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion.
    use Numerical_Comparison
    implicit none
    character(len=*),   intent(in)           :: testName
    double precision,   intent(in)           :: value1,value2
    integer,            intent(in), optional :: compare
    double precision,   intent(in), optional :: absTol,relTol
    type(assertResult), pointer              :: thisResult
    integer                                  :: compareActual
    logical                                  :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if
    
    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals)
       passed=.not.Values_Differ(value1,value2,absTol,relTol)
    case (compareNotEqual          )
       passed=(value1 /= value2)
    case (compareLessThan          )
       passed=(value1  < value2)
    case (compareGreaterThan       )
       passed=(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=(value1 >= value2)
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert

  subroutine Unit_Tests_Begin_Group(groupName)
    !% Marks that a unit test group has begun.
    implicit none
    character(len=*),   intent(in) :: groupName
    type(assertResult), pointer    :: thisResult

    thisResult => Get_New_Assert_Result()
    thisResult%beginGroup=.true.
    thisResult%label     =trim(groupName)
    return
  end subroutine Unit_Tests_Begin_Group

  subroutine Unit_Tests_End_Group
    !% Marks that a unit test group has ended.
    implicit none
    type(assertResult), pointer :: thisResult

    thisResult => Get_New_Assert_Result()
    thisResult%endGroup=.true.
    return
  end subroutine Unit_Tests_End_Group

  subroutine Unit_Tests_Finish
    !% Write out the results of unit testing.
    use Galacticus_Display
    use String_Handling
    implicit none
    type(assertResult),   pointer :: thisResult,nextResult
    integer                       :: passCount,failCount,percentage
    type(varying_string)          :: message 

    passCount=0
    failCount=0
    thisResult => firstResult
    do while (associated(thisResult))
       if (.not.(thisResult%beginGroup.or.thisResult%endGroup)) then
          select case (thisResult%result)
          case (.false.)
             failCount=failCount+1
             message="FAILED: "//thisResult%label  
          case (.true. )
             passCount=passCount+1             
             message="passed: "//thisResult%label  
          end select
          call Galacticus_Display_Message(message)
       else
          if (thisResult%beginGroup) call Galacticus_Display_Indent  (thisResult%label)
          if (thisResult%endGroup  ) call Galacticus_Display_Unindent("")
       end if
       thisResult => thisResult%nextResult
    end do
    percentage=int(100.0d0*dble(passCount)/dble(passCount+failCount))
    message="Tests passed: "
    message=message//passCount//" ("//percentage//"%)"
    call Galacticus_Display_Message(message)
    percentage=int(100.0d0*dble(failCount)/dble(passCount+failCount))
    message="Tests failed: "
    message=message//failCount//" ("//percentage//"%)"
    call Galacticus_Display_Message(message)

    ! Destroy the list of results.
    thisResult => firstResult%nextResult
    do while (associated(thisResult))
       nextResult => thisResult%nextResult
       deallocate(thisResult)
       thisResult => nextResult
    end do
    return
  end subroutine Unit_Tests_Finish

  function Get_New_Assert_Result() result(newResult)
    !% Get a new assert result object.
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
    newResult%result    =.false.
    newResult%nextResult => null()
    return
  end function Get_New_Assert_Result

end module Unit_Tests
