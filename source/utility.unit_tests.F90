!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements unit testing.

module Unit_Tests
  !% Implements unit testing.
  use :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: Assert, Skip, Unit_Tests_Finish, Unit_Tests_Begin_Group, Unit_Tests_End_Group

  ! Types of comparison.
  integer, parameter, public :: compareEquals            =0
  integer, parameter, public :: compareNotEqual          =1
  integer, parameter, public :: compareLessThan          =2
  integer, parameter, public :: compareGreaterThan       =3
  integer, parameter, public :: compareLessThanOrEqual   =4
  integer, parameter, public :: compareGreaterThanOrEqual=5

  !# <enumeration>
  !#  <name>test</name>
  !#  <description>Statuses for unit tests.</description>
  !#  <entry label="passed"/>
  !#  <entry label="failed" />
  !#  <entry label="skipped"/>
  !# </enumeration>

  ! Type for assert results.
  type assertResult
     !% A derived type for storing results of asserts.
     integer                          :: result
     logical                          :: beginGroup, endGroup
     type   (varying_string)          :: label     , note
     type   (assertResult  ), pointer :: nextResult
  end type assertResult

  ! Results list.
  type   (assertResult), target  :: firstResult
  type   (assertResult), pointer :: currentResult
  logical                        :: firstAssert  =.true.

  ! Interface for assert routines.
  interface Assert
     !% Generic interface for assert routines.
     module procedure Assert_Real_Scalar
     module procedure Assert_Double_Scalar
     module procedure Assert_Integer_Scalar
     module procedure Assert_Integer8_Scalar
     module procedure Assert_Character_Scalar
     module procedure Assert_VarString_Scalar
     module procedure Assert_Logical_Scalar
     module procedure Assert_Real_1D_Array
     module procedure Assert_Double_1D_Array
     module procedure Assert_Double_Complex_1D_Array
     module procedure Assert_Integer_1D_Array
     module procedure Assert_Integer8_1D_Array
     module procedure Assert_Character_1D_Array
     module procedure Assert_VarString_1D_Array
     module procedure Assert_Logical_1D_Array
     module procedure Assert_Double_2D_Array
     module procedure Assert_Double_3D_Array
     module procedure Assert_Double_4D_Array
     module procedure Assert_Double_5D_Array
  end interface Assert

contains

  integer function getStatus(passed)
    !% Return the status code for a test on the basis of a boolean pass/fail.
    implicit none
    logical, intent(in   ) :: passed

    if (passed) then
       getStatus=testPassed
    else
       getStatus=testFailed
    end if
    return
  end function getStatus

  subroutine Skip(testName,reason)
    !% Record that a test was skipped.
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

  subroutine Assert_Real_Scalar(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about real arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : assignment(=)
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    character(len=*       ), intent(in   )           :: testName
    real                   , intent(in   )           :: value1       , value2
    integer                , intent(in   ), optional :: compare
    real                   , intent(in   ), optional :: absTol       , relTol
    type     (assertResult), pointer                 :: result
    integer                                          :: compareActual
    logical                                          :: passed

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
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Real_Scalar

  subroutine Assert_Double_Scalar(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : operator(//)           , assignment(=)
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    character       (len=*       ), intent(in   )           :: testName
    double precision              , intent(in   )           :: value1       , value2
    integer                       , intent(in   ), optional :: compare
    double precision              , intent(in   ), optional :: absTol       , relTol
    type            (assertResult), pointer                 :: result
    integer                                                 :: compareActual
    logical                                                 :: passed
    character       (len=16      )                          :: label        , comparison

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
       comparison="≉"
    case (compareNotEqual          )
       passed=(value1 /= value2)
       comparison="="
    case (compareLessThan          )
       passed=(value1  < value2)
       comparison="≮"
    case (compareGreaterThan       )
       passed=(value1  > value2)
       comparison="≯"
    case (compareLessThanOrEqual   )
       passed=(value1 <= value2)
       comparison="≰"
    case (compareGreaterThanOrEqual)
       passed=(value1 >= value2)
       comparison="≱"
    case default
       passed=.false.
       comparison=""
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)
    if (.not.passed) then
       write (label,'(e16.8)') value1
       result%note=                 trim(adjustl(label))
       result%note=result%note//" "//trim(comparison)//" "
       write (label,'(e16.8)') value2
       result%note=result%note//trim(adjustl(label))
    end if
    return
  end subroutine Assert_Double_Scalar

  subroutine Assert_Integer_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*       ), intent(in   )           :: testName
    integer                , intent(in   )           :: value1       , value2
    integer                , intent(in   ), optional :: compare
    type     (assertResult), pointer                 :: result
    integer                                          :: compareActual
    logical                                          :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=(value1 == value2)
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
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Integer_Scalar

  subroutine Assert_Integer8_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    use :: Kind_Numbers      , only : kind_int8
    implicit none
    character(len=*         ), intent(in   )           :: testName
    integer  (kind=kind_int8), intent(in   )           :: value1       , value2
    integer                  , intent(in   ), optional :: compare
    type     (assertResult  ), pointer                 :: result
    integer                                            :: compareActual
    logical                                            :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=(value1 == value2)
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
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Integer8_Scalar

  subroutine Assert_Character_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about character arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*       ), intent(in   )           :: testName
    character(len=*       ), intent(in   )           :: value1       , value2
    integer                , intent(in   ), optional :: compare
    type     (assertResult), pointer                 :: result
    integer                                          :: compareActual
    logical                                          :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=(value1 == value2)
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
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Character_Scalar

  subroutine Assert_VarString_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about character arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : operator(==)           , operator(/=), operator(<)  , operator(>), &
         &                            operator(>=)           , operator(<=), assignment(=)
    implicit none
    character(len=*         ), intent(in   )           :: testName
    type     (varying_string), intent(in   )           :: value1       , value2
    integer                  , intent(in   ), optional :: compare
    type     (assertResult  ), pointer                 :: result
    integer                                            :: compareActual
    logical                                            :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=(value1 == value2)
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
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_VarString_Scalar

  subroutine Assert_Real_1D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about real arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : assignment(=)
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    character(len=*       )              , intent(in   )           :: testName
    real                   , dimension(:), intent(in   )           :: value1       , value2
    integer                              , intent(in   ), optional :: compare
    real                                 , intent(in   ), optional :: absTol       , relTol
    type     (assertResult), pointer                               :: result
    integer                                                        :: compareActual, iTest
    logical                                                        :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals)
       passed=.true.
       do iTest=1,min(size(value1),size(value2))
          if (Values_Differ(value1(iTest),value2(iTest),absTol,relTol)) then
             passed=.false.
             exit
          end if
       end do
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Real_1D_Array

  subroutine Assert_Double_1D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : assignment(=)          , operator(//), var_str
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    character       (len=*       )                         , intent(in   )           :: testName
    double precision              , dimension(          : ), intent(in   )           :: value1       , value2
    integer                                                , intent(in   ), optional :: compare
    double precision                                       , intent(in   ), optional :: absTol       , relTol
    type            (assertResult), pointer                                          :: result
    integer                                                                          :: compareActual, iTest
    logical                       , dimension(size(value1))                          :: passed
    character       (len=22      )                                                   :: label        , comparison

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals)
       passed=.true.
       do iTest=1,min(size(value1),size(value2))
          if (.not.Values_Agree(value1(iTest),value2(iTest),absTol,relTol)) passed(iTest)=.false.
       end do
       comparison="≉"
   case (compareNotEqual          )
       passed=value1 /= value2
       comparison="="
    case (compareLessThan          )
       passed=value1  < value2
       comparison="≮"
    case (compareGreaterThan       )
       passed=value1  > value2
       comparison="≯"
    case (compareLessThanOrEqual   )
       passed=value1 <= value2
       comparison="≰"
    case (compareGreaterThanOrEqual)
       passed=value1 >= value2
       comparison="≱"
    case default
       passed=.false.
       comparison=""
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(all(passed))
    result%label =trim(testName)
    if (.not.all(passed)) then
       result%note=var_str('')
       do iTest=1,min(size(value1),size(value2))
          if (passed(iTest)) cycle
          write (label,'(i2,1x,a1,1x,e16.8)') iTest,':',value1(iTest)
          result%note=result%note//trim(adjustl(label))
          result%note=result%note//" "//trim(comparison)//" "
          write (label,'(e16.8)') value2(iTest)
          result%note=result%note//trim(adjustl(label))//"; "
       end do
    end if
    return
  end subroutine Assert_Double_1D_Array

  subroutine Assert_Double_Complex_1D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double complex arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : assignment(=)
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    character       (len=*       )              , intent(in   )           :: testName
    double complex                , dimension(:), intent(in   )           :: value1       , value2
    integer                                     , intent(in   ), optional :: compare
    double complex                              , intent(in   ), optional :: absTol       , relTol
    double complex                                                        :: absTolActual , relTolActual
    type            (assertResult), pointer                               :: result
    integer                                                               :: compareActual, iTest
    logical                                                               :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if
    ! Determine tolerances.
    if (present(absTol)) then
       absTolActual=absTol
    else
       absTolActual=dcmplx(huge(1.0d0),huge(1.0d0))
    end if
    if (present(relTol)) then
       relTolActual=relTol
    else
       relTolActual=dcmplx(huge(1.0d0),huge(1.0d0))
    end if
    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals)
       passed=.true.
       do iTest=1,min(size(value1),size(value2))
          if     (                                                                                              &
               &   Values_Differ(real(value1(iTest)),real(value2(iTest)),real(absTolActual),real(relTolActual)) &
               &  .or.                                                                                          &
               &   Values_Differ(imag(value1(iTest)),imag(value2(iTest)),imag(absTolActual),imag(relTolActual)) &
               & ) then
             passed=.false.
             exit
          end if
       end do
    case default
       passed=.false.
       call Galacticus_Error_Report('comparison not supported for complex values'//{introspection:location})
    end select
    ! Get an object to store the results in.
    result => Get_New_Assert_Result()
    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)
    return
  end subroutine Assert_Double_Complex_1D_Array

  subroutine Assert_Double_2D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : assignment(=)
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    character       (len=*       )                , intent(in   )           :: testName
    double precision              , dimension(:,:), intent(in   )           :: value1       , value2
    integer                                       , intent(in   ), optional :: compare
    double precision                              , intent(in   ), optional :: absTol       , relTol
    type            (assertResult), pointer                                 :: result
    integer                                                                 :: compareActual, iTest , jTest
    logical                                                                 :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals)
       passed=.true.
       do iTest=1,min(size(value1,dim=1),size(value2,dim=1))
          do jTest=1,min(size(value1,dim=2),size(value2,dim=2))
             if (.not.Values_Agree(value1(iTest,jTest),value2(iTest,jTest),absTol,relTol)) then
                passed=.false.
                exit
             end if
          end do
       end do
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Double_2D_Array

  subroutine Assert_Double_3D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : assignment(=)
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    character       (len=*       )                  , intent(in   )           :: testName
    double precision              , dimension(:,:,:), intent(in   )           :: value1       , value2
    integer                                         , intent(in   ), optional :: compare
    double precision                                , intent(in   ), optional :: absTol       , relTol
    type            (assertResult), pointer                                   :: result
    integer                                                                   :: compareActual, iTest , jTest, kTest
    logical                                                                   :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals)
       passed=.true.
       do iTest=1,min(size(value1,dim=1),size(value2,dim=1))
          do jTest=1,min(size(value1,dim=2),size(value2,dim=2))
             do kTest=1,min(size(value1,dim=3),size(value2,dim=3))
                if (.not.Values_Agree(value1(iTest,jTest,kTest),value2(iTest,jTest,kTest),absTol,relTol)) then
                   passed=.false.
                   exit
                end if
             end do
          end do
       end do
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Double_3D_Array

  subroutine Assert_Double_4D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : assignment(=)
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    character       (len=*       )                    , intent(in   )           :: testName
    double precision              , dimension(:,:,:,:), intent(in   )           :: value1       , value2
    integer                                           , intent(in   ), optional :: compare
    double precision                                  , intent(in   ), optional :: absTol       , relTol
    type            (assertResult), pointer                                     :: result
    integer                                                                     :: compareActual, iTest , jTest, kTest, &
         &                                                                         lTest
    logical                                                                     :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals)
       passed=.true.
       do iTest=1,min(size(value1,dim=1),size(value2,dim=1))
          do jTest=1,min(size(value1,dim=2),size(value2,dim=2))
             do kTest=1,min(size(value1,dim=3),size(value2,dim=3))
                do lTest=1,min(size(value1,dim=4),size(value2,dim=4))
                   if (.not.Values_Agree(value1(iTest,jTest,kTest,lTest),value2(iTest,jTest,kTest,lTest),absTol,relTol)) then
                      passed=.false.
                      exit
                   end if
                end do
             end do
          end do
       end do
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Double_4D_Array

  subroutine Assert_Double_5D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use :: Galacticus_Error    , only : Galacticus_Error_Report
    use :: ISO_Varying_String  , only : assignment(=)
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    character       (len=*       )                      , intent(in   )           :: testName
    double precision              , dimension(:,:,:,:,:), intent(in   )           :: value1       , value2
    integer                                             , intent(in   ), optional :: compare
    double precision                                    , intent(in   ), optional :: absTol       , relTol
    type            (assertResult), pointer                                       :: result
    integer                                                                       :: compareActual, iTest , jTest, kTest, &
         &                                                                           lTest        , mTest
    logical                                                                       :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals)
       passed=.true.
       do iTest=1,min(size(value1,dim=1),size(value2,dim=1))
          do jTest=1,min(size(value1,dim=2),size(value2,dim=2))
             do kTest=1,min(size(value1,dim=3),size(value2,dim=3))
                do lTest=1,min(size(value1,dim=4),size(value2,dim=4))
                   do mTest=1,min(size(value1,dim=5),size(value2,dim=5))
                      if (.not.Values_Agree(value1(iTest,jTest,kTest,lTest,mTest),value2(iTest,jTest,kTest,lTest,mTest),absTol,relTol)) then
                         passed=.false.
                         exit
                      end if
                   end do
                end do
             end do
          end do
       end do
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Double_5D_Array

  subroutine Assert_Integer_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*       )              , intent(in   )           :: testName
    integer                , dimension(:), intent(in   )           :: value1       , value2
    integer                              , intent(in   ), optional :: compare
    type     (assertResult), pointer                               :: result
    integer                                                        :: compareActual
    logical                                                        :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=all(value1 == value2)
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Integer_1D_Array

  subroutine Assert_Logical_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*       )              , intent(in   )           :: testName
    logical                , dimension(:), intent(in   )           :: value1       , value2
    integer                              , intent(in   ), optional :: compare
    type     (assertResult), pointer                               :: result
    integer                                                        :: compareActual
    logical                                                        :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals  )
       passed=all(value1 .eqv.  value2)
    case (compareNotEqual)
       passed=all(value1 .neqv. value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('assertions about logical variables can be equality or non-equality only'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Logical_1D_Array

  subroutine Assert_Integer8_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    use :: Kind_Numbers      , only : kind_int8
    implicit none
    character(len=*         )              , intent(in   )           :: testName
    integer  (kind=kind_int8), dimension(:), intent(in   )           :: value1       , value2
    integer                                , intent(in   ), optional :: compare
    type     (assertResult  ), pointer                               :: result
    integer                                                          :: compareActual
    logical                                                          :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=all(value1 == value2)
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Integer8_1D_Array

  subroutine Assert_Logical_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about logical arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*       ), intent(in   )           :: testName
    logical                , intent(in   )           :: value1       , value2
    integer                , intent(in   ), optional :: compare
    type     (assertResult), pointer                 :: result
    integer                                          :: compareActual
    logical                                          :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=value1 .eqv.  value2
    case (compareNotEqual          )
       passed=value1 .neqv. value2
    case default
       passed=.false.
       call Galacticus_Error_Report('assertions about logical variables can be equality or non-equality only'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Logical_Scalar

  subroutine Assert_Character_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about character arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*       )              , intent(in   )           :: testName
    character(len=*       ), dimension(:), intent(in   )           :: value1       , value2
    integer                              , intent(in   ), optional :: compare
    type     (assertResult), pointer                               :: result
    integer                                                        :: compareActual
    logical                                                        :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=all(value1 == value2)
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_Character_1D_Array

  subroutine Assert_VarString_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about character arguments.
    use :: Galacticus_Error  , only : Galacticus_Error_Report
    use :: ISO_Varying_String, only : operator(==)           , operator(/=), operator(<)  , operator(>), &
         &                            operator(>=)           , operator(<=), assignment(=)
    implicit none
    character(len=*         )              , intent(in   )           :: testName
    type     (varying_string), dimension(:), intent(in   )           :: value1       , value2
    integer                                , intent(in   ), optional :: compare
    type     (assertResult  ), pointer                               :: result
    integer                                                          :: compareActual
    logical                                                          :: passed

    ! Determine what type of comparison to perform.
    if (present(compare)) then
       compareActual=compare
    else
       compareActual=compareEquals
    end if

    ! Perform the comparison.
    select case (compareActual)
    case (compareEquals            )
       passed=all(value1 == value2)
    case (compareNotEqual          )
       passed=all(value1 /= value2)
    case (compareLessThan          )
       passed=all(value1  < value2)
    case (compareGreaterThan       )
       passed=all(value1  > value2)
    case (compareLessThanOrEqual   )
       passed=all(value1 <= value2)
    case (compareGreaterThanOrEqual)
       passed=all(value1 >= value2)
    case default
       passed=.false.
       call Galacticus_Error_Report('unknown comparison'//{introspection:location})
    end select

    ! Get an object to store the results in.
    result => Get_New_Assert_Result()

    ! Store the result.
    result%result=getStatus(passed)
    result%label =trim(testName)

    return
  end subroutine Assert_VarString_1D_Array

  subroutine Unit_Tests_Begin_Group(groupName)
    !% Marks that a unit test group has begun.
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
    !% Marks that a unit test group has ended.
    implicit none
    type(assertResult), pointer :: result

    result => Get_New_Assert_Result()
    result%endGroup=.true.
    return
  end subroutine Unit_Tests_End_Group

  subroutine Unit_Tests_Finish
    !% Write out the results of unit testing.
    use :: Galacticus_Display, only : Galacticus_Display_Indent, Galacticus_Display_Message, Galacticus_Display_Unindent
    use :: ISO_Varying_String, only : assignment(=)            , operator(//)              , operator(/=)
    use :: Memory_Management , only : Memory_Usage_Record
    use :: String_Handling   , only : operator(//)
    implicit none
    type   (assertResult  ), pointer :: nextResult, result
    integer                          :: failCount , passCount , percentage
    type   (varying_string)          :: message

    passCount=0
    failCount=0
    result => firstResult
    do while (associated(result))
       if (.not.(result%beginGroup.or.result%endGroup)) then
          select case (result%result)
          case (testFailed)
             failCount=failCount+1
             message=" FAILED: "//result%label
             if (result%note /= "") message=message//" ["//result%note//"]"
          case (testPassed)
             passCount=passCount+1
             message=" passed: "//result%label
          case (testSkipped)
             message="skipped: "//result%label//" ["//result%note//"]"
          end select
          call Galacticus_Display_Message(message)
       else
          if (result%beginGroup) call Galacticus_Display_Indent  (result%label)
          if (result%endGroup  ) call Galacticus_Display_Unindent("")
       end if
       result => result%nextResult
    end do
    if (passCount+failCount > 0) then
       percentage=int(100.0d0*dble(passCount)/dble(passCount+failCount))
    else
       percentage=100
    end if
    message="Tests passed: "
    message=message//passCount//" ("//percentage//"%)"
    call Galacticus_Display_Message(message)
    if (passCount+failCount > 0) then
       percentage=int(100.0d0*dble(failCount)/dble(passCount+failCount))
    else
       percentage=100
    end if
    message="Tests failed: "
    message=message//failCount//" ("//percentage//"%)"
    call Galacticus_Display_Message(message)

    ! Destroy the list of results.
    result => firstResult%nextResult
    do while (associated(result))
       nextResult => result%nextResult
       call result%label%destroy()
       deallocate(result)
       call Memory_Usage_Record(sizeof(result),addRemove=-1)
       result => nextResult
    end do
    return
  end subroutine Unit_Tests_Finish

  function Get_New_Assert_Result() result(newResult)
    !% Get a new assert result object.
    use :: ISO_Varying_String, only : assignment(=)
    use :: Memory_Management , only : Memory_Usage_Record
    implicit none
    type(assertResult), pointer :: newResult

    ! Return the first result if this is the first assert, otherwise, allocate a new one.
    select case (firstAssert)
    case (.false.)
       allocate(currentResult%nextResult)
       call Memory_Usage_Record(sizeof(currentResult%nextResult))
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
