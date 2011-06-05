!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


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

  ! Interface for assert routines.
  interface Assert
     !% Generic interface for assert routines.
     module procedure Assert_Double_Scalar
     module procedure Assert_Integer_Scalar
     module procedure Assert_Integer8_Scalar
     module procedure Assert_Character_Scalar
     module procedure Assert_VarString_Scalar
     module procedure Assert_Double_1D_Array
     module procedure Assert_Integer_1D_Array
     module procedure Assert_Integer8_1D_Array
     module procedure Assert_Character_1D_Array
     module procedure Assert_VarString_1D_Array
     module procedure Assert_Double_2D_Array
     module procedure Assert_Double_3D_Array
     module procedure Assert_Double_4D_Array
     module procedure Assert_Double_5D_Array
  end interface Assert

contains

  subroutine Assert_Double_Scalar(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
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
  end subroutine Assert_Double_Scalar

  subroutine Assert_Integer_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    implicit none
    character(len=*),   intent(in)           :: testName
    integer,            intent(in)           :: value1,value2
    integer,            intent(in), optional :: compare
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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Integer_Scalar

  subroutine Assert_Integer8_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    use Kind_Numbers
    implicit none
    character(len=*),        intent(in)           :: testName
    integer(kind=kind_int8), intent(in)           :: value1,value2
    integer,                 intent(in), optional :: compare
    type(assertResult),      pointer              :: thisResult
    integer                                       :: compareActual
    logical                                       :: passed

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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Integer8_Scalar

  subroutine Assert_Character_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about character arguments.
    implicit none
    character(len=*),   intent(in)           :: testName
    character(len=*),   intent(in)           :: value1,value2
    integer,            intent(in), optional :: compare
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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Character_Scalar

  subroutine Assert_VarString_Scalar(testName,value1,value2,compare)
    !% Assess and record an assertion about character arguments.
    implicit none
    character(len=*),     intent(in)           :: testName
    type(varying_string), intent(in)           :: value1,value2
    integer,              intent(in), optional :: compare
    type(assertResult),   pointer              :: thisResult
    integer                                    :: compareActual
    logical                                    :: passed

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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_VarString_Scalar

  subroutine Assert_Double_1D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use Numerical_Comparison
    implicit none
    character(len=*),   intent(in)               :: testName
    double precision,   intent(in), dimension(:) :: value1,value2
    integer,            intent(in), optional     :: compare
    double precision,   intent(in), optional     :: absTol,relTol
    type(assertResult), pointer                  :: thisResult
    integer                                      :: compareActual,iTest
    logical                                      :: passed

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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Double_1D_Array

  subroutine Assert_Double_2D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use Numerical_Comparison
    implicit none
    character(len=*),   intent(in)                 :: testName
    double precision,   intent(in), dimension(:,:) :: value1,value2
    integer,            intent(in), optional       :: compare
    double precision,   intent(in), optional       :: absTol,relTol
    type(assertResult), pointer                    :: thisResult
    integer                                        :: compareActual,iTest,jTest
    logical                                        :: passed

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
             if (Values_Differ(value1(iTest,jTest),value2(iTest,jTest),absTol,relTol)) then
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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Double_2D_Array

  subroutine Assert_Double_3D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use Numerical_Comparison
    implicit none
    character(len=*),   intent(in)                   :: testName
    double precision,   intent(in), dimension(:,:,:) :: value1,value2
    integer,            intent(in), optional         :: compare
    double precision,   intent(in), optional         :: absTol,relTol
    type(assertResult), pointer                      :: thisResult
    integer                                          :: compareActual,iTest,jTest,kTest
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
       passed=.true.
       do iTest=1,min(size(value1,dim=1),size(value2,dim=1))
          do jTest=1,min(size(value1,dim=2),size(value2,dim=2))
             do kTest=1,min(size(value1,dim=3),size(value2,dim=3))
                if (Values_Differ(value1(iTest,jTest,kTest),value2(iTest,jTest,kTest),absTol,relTol)) then
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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Double_3D_Array

  subroutine Assert_Double_4D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use Numerical_Comparison
    implicit none
    character(len=*),   intent(in)                     :: testName
    double precision,   intent(in), dimension(:,:,:,:) :: value1,value2
    integer,            intent(in), optional           :: compare
    double precision,   intent(in), optional           :: absTol,relTol
    type(assertResult), pointer                        :: thisResult
    integer                                            :: compareActual,iTest,jTest,kTest,lTest
    logical                                            :: passed

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
                   if (Values_Differ(value1(iTest,jTest,kTest,lTest),value2(iTest,jTest,kTest,lTest),absTol,relTol)) then
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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Double_4D_Array
  
  subroutine Assert_Double_5D_Array(testName,value1,value2,compare,absTol,relTol)
    !% Assess and record an assertion about double precision arguments.
    use Numerical_Comparison
    implicit none
    character(len=*),   intent(in)                       :: testName
    double precision,   intent(in), dimension(:,:,:,:,:) :: value1,value2
    integer,            intent(in), optional             :: compare
    double precision,   intent(in), optional             :: absTol,relTol
    type(assertResult), pointer                          :: thisResult
    integer                                              :: compareActual,iTest,jTest,kTest,lTest,mTest
    logical                                              :: passed

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
                      if (Values_Differ(value1(iTest,jTest,kTest,lTest,mTest),value2(iTest,jTest,kTest,lTest,mTest),absTol,relTol)) then
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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Double_5D_Array
  
  subroutine Assert_Integer_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    implicit none
    character(len=*),   intent(in)               :: testName
    integer,            intent(in), dimension(:) :: value1,value2
    integer,            intent(in), optional     :: compare
    type(assertResult), pointer                  :: thisResult
    integer                                      :: compareActual
    logical                                      :: passed

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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Integer_1D_Array

  subroutine Assert_Integer8_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about integer arguments.
    use Kind_Numbers
    implicit none
    character(len=*),        intent(in)               :: testName
    integer(kind=kind_int8), intent(in), dimension(:) :: value1,value2
    integer,                 intent(in), optional     :: compare
    type(assertResult),      pointer                  :: thisResult
    integer                                           :: compareActual
    logical                                           :: passed

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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Integer8_1D_Array

  subroutine Assert_Character_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about character arguments.
    implicit none
    character(len=*),   intent(in)               :: testName
    character(len=*),   intent(in), dimension(:) :: value1,value2
    integer,            intent(in), optional     :: compare
    type(assertResult), pointer                  :: thisResult
    integer                                      :: compareActual
    logical                                      :: passed

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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_Character_1D_Array

  subroutine Assert_VarString_1D_Array(testName,value1,value2,compare)
    !% Assess and record an assertion about character arguments.
    implicit none
    character(len=*),     intent(in)               :: testName
    type(varying_string), intent(in), dimension(:) :: value1,value2
    integer,              intent(in), optional     :: compare
    type(assertResult),   pointer                  :: thisResult
    integer                                        :: compareActual
    logical                                        :: passed

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
    end select

    ! Get an object to store the results in.
    thisResult => Get_New_Assert_Result()

    ! Store the result.
    thisResult%result=passed
    thisResult%label =trim(testName)

    return
  end subroutine Assert_VarString_1D_Array

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
    use Memory_Management
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
    if (passCount+FailCount > 0) then
       percentage=int(100.0d0*dble(passCount)/dble(passCount+failCount))
    else
       percentage=100.0d0
    end if
    message="Tests passed: "
    message=message//passCount//" ("//percentage//"%)"
    call Galacticus_Display_Message(message)
    if (passCount+FailCount > 0) then
       percentage=int(100.0d0*dble(failCount)/dble(passCount+failCount))
    else
       percentage=100.0d0
    end if
    message="Tests failed: "
    message=message//failCount//" ("//percentage//"%)"
    call Galacticus_Display_Message(message)

    ! Destroy the list of results.
    thisResult => firstResult%nextResult
    do while (associated(thisResult))
       nextResult => thisResult%nextResult
       call thisResult%label%destroy()
       deallocate(thisResult)
       call Memory_Usage_Record(sizeof(thisResult),addRemove=-1)
       thisResult => nextResult
    end do
    return
  end subroutine Unit_Tests_Finish

  function Get_New_Assert_Result() result(newResult)
    !% Get a new assert result object.
    use Memory_Management
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
    newResult%result    =.false.
    newResult%nextResult => null()
    return
  end function Get_New_Assert_Result

end module Unit_Tests
