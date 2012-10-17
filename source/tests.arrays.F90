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

!% Contains a program to test the array functions.

program Test_Array_Monotonicity
  !% Tests that array functions.
  use Unit_Tests
  use ISO_Varying_String
  use Array_Utilities
  use Kind_Numbers
  implicit none
  double precision,        target , dimension( 1,2) :: singleElementArrays       =reshape([                                                                         &
       &                                                                                     1.23d0                                                                 &
       &                                                                                   ,-2.31d0                                                                 &
       &                                                                                  ]                                                                         &
       &                                                                                  ,shape(singleElementArrays)                                               &
       &                                                                                 )
  integer(kind=kind_int8), target , dimension( 1,2) :: singleElementArraysInteger=reshape([                                                                         &
       &                                                                                     123                                                                    &
       &                                                                                   ,-231                                                                    &
       &                                                                                  ]                                                                         &
       &                                                                                  ,shape(singleElementArraysInteger)                                        &
       &                                                                                 )
  logical,                 target , dimension( 9,2) :: singleElementExpectations =reshape([                                                                         &
       &                                                                                    .true. ,.true. ,.true. ,.true. ,.true. ,.true. ,.true. ,.true. ,.true.  &
       &                                                                                   ,.true. ,.true. ,.true. ,.true. ,.true. ,.true. ,.true. ,.true. ,.true.  &
       &                                                                                  ]                                                                         &
       &                                                                                  ,shape(singleElementExpectations)                                         &
       &                                                                                 )
  character(len=128),      target , dimension(   2) :: singleElementNames        =        [                                                                         &
       &                                                                                    'Single element array (positive value)'                                 &
       &                                                                                   ,'Single element array (negative value)'                                 &
       &                                                                                  ]
  double precision,        target , dimension(10,6) :: tenElementArrays          =reshape([                                                                         &
       &                                                                                    1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0            &
       &                                                                                   ,3.0d0,2.0d0,1.0d0,0.0d0,-1.0d0,-2.0d0,-3.0d0,-4.0d0,-5.0d0,-6.0d0       &
       &                                                                                   ,1.0d0,0.0d0,3.0d0,4.0d0,-3.0d0,5.0d0,1.0d0,8.0d0,2.0d0,3.0d0            &
       &                                                                                   ,1.0d0,1.0d0,1.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0            &
       &                                                                                   ,3.0d0,2.0d0,1.0d0,1.0d0,-1.0d0,-2.0d0,-3.0d0,-4.0d0,-5.0d0,-6.0d0       &
       &                                                                                   ,1.0d0,0.0d0,3.0d0,3.0d0,-3.0d0,5.0d0,1.0d0,8.0d0,2.0d0,3.0d0            &
       &                                                                                  ]                                                                         &
       &                                                                                  ,shape(tenElementArrays)                                                  &
       &                                                                                 )
  integer(kind=kind_int8),            target , dimension(10,6) :: tenElementArraysInteger   =reshape([                                                                         &
       &                                                                                    10,20,30,40, 50, 60, 70, 80, 90, 100                                    &
       &                                                                                   ,30,20,10,00,-10,-20,-30,-40,-50,-60                                     &
       &                                                                                   ,10,00,30,40,-30, 50, 10, 80, 20, 30                                     &
       &                                                                                   ,10,10,10,40, 50, 60, 70, 80, 90,100                                     &
       &                                                                                   ,30,20,10,10,-10,-20,-30,-40,-50,-60                                     &
       &                                                                                   ,10,00,30,30,-30, 50, 10, 80, 20, 30                                     &
       &                                                                                  ]                                                                         &
       &                                                                                  ,shape(tenElementArrays)                                                  &
       &                                                                                 )
  logical,                 target , dimension( 9,6) :: tenElementExpectations    =reshape([                                                                         &
       &                                                                                    .true. ,.true. ,.false.,.true. ,.true. ,.false.,.true. ,.true. ,.false. &
       &                                                                                   ,.true. ,.false.,.true. ,.true. ,.false.,.true. ,.true. ,.false.,.true.  &
       &                                                                                   ,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false. &
       &                                                                                   ,.false.,.false.,.false.,.false.,.false.,.false.,.true. ,.true. ,.false. &
       &                                                                                   ,.false.,.false.,.false.,.false.,.false.,.false.,.true. ,.false.,.true.  &
       &                                                                                   ,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false. &
       &                                                                                  ]                                                                         &
       &                                                                                  ,shape(tenElementExpectations)                                            &
       &                                                                                 )
  character(len=128),      target , dimension(   6) :: tenElementNames           =        [                                                                         &
       &                                                                                    'Increasing array (no equalities)     '                                 &
       &                                                                                   ,'Decreasing array (no equalities)     '                                 &
       &                                                                                   ,'Non-monotinic array (no equalities)  '                                 &
       &                                                                                   ,'Increasing array (with equalities)   '                                 &
       &                                                                                   ,'Decreasing array (with equalities)   '                                 &
       &                                                                                   ,'Non-monotinic array (with equalities)'                                 &
       &                                                                                  ]
  double precision,                 dimension(  10) :: doubleArray         =[0.0d0,1.0d0,2.0d0,3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0]
  double precision,                 dimension(  10) :: doubleArrayReversed =[9.0d0,8.0d0,7.0d0,6.0d0, 5.0d0, 4.0d0, 3.0d0, 2.0d0, 1.0d0, 0.0d0]
  double precision,                 dimension(  10) :: doubleArrayCumulated=[0.0d0,1.0d0,3.0d0,6.0d0,10.0d0,15.0d0,21.0d0,28.0d0,36.0d0,45.0d0]
  real,                             dimension(  10) :: realArray           =[0.0e0,1.0e0,2.0e0,3.0e0, 4.0e0, 5.0e0, 6.0e0, 7.0e0, 8.0e0, 9.0e0]
  real,                             dimension(  10) :: realArrayReversed   =[9.0e0,8.0e0,7.0e0,6.0e0, 5.0e0, 4.0e0, 3.0e0, 2.0e0, 1.0e0, 0.0e0]
  double precision,        pointer, dimension( :,:) :: thisArraySet
  integer(kind=kind_int8), pointer, dimension( :,:) :: thisArraySetInteger
  logical,                 pointer, dimension( :,:) :: thisExpectations
  character(len=128),      pointer, dimension(   :) :: thisNames
  logical                                           :: isMonotonic
  integer                                           :: iArraySet,iArray
  type(varying_string)                              :: test

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Array functions")

  ! Loop over all sets of arrays.
  do iArraySet=1,2

     ! Create pointers to the specific array set to process.
     select case (iArraySet)
     case (1)
        thisArraySet     => singleElementArrays
        thisExpectations => singleElementExpectations
        thisNames        => singleElementNames
     case (2)
        thisArraySet     => tenElementArrays
        thisExpectations => tenElementExpectations
        thisNames        => tenElementNames
     end select

     ! Loop over all arrays in this array set.
     do iArray=1,size(thisArraySet,dim=2)

        ! Perform tests on this array for each of the nine permutations of options:
        !  direction=   increasing|decreasing|undefined
        !  allowEquals= true      |false     |undefined

        test=trim(thisNames(iArray))//' is monotonic [no direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray)                                                 )
        call Assert(char(test),isMonotonic,thisExpectations(1,iArray))

        test=trim(thisNames(iArray))//' is monotonic [increasing direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray),direction=directionIncreasing                   )
        call Assert(char(test),isMonotonic,thisExpectations(2,iArray))

        test=trim(thisNames(iArray))//' is monotonic [decreasing direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray),direction=directionDecreasing                   )
        call Assert(char(test),isMonotonic,thisExpectations(3,iArray))

        test=trim(thisNames(iArray))//' is monotonic [no direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray)                              ,allowEqual=.false.)
        call Assert(char(test),isMonotonic,thisExpectations(4,iArray))

        test=trim(thisNames(iArray))//' is monotonic [increasing direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray),direction=directionIncreasing,allowEqual=.false.)
        call Assert(char(test),isMonotonic,thisExpectations(5,iArray))

        test=trim(thisNames(iArray))//' is monotonic [decreasing direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray),direction=directionDecreasing,allowEqual=.false.)
        call Assert(char(test),isMonotonic,thisExpectations(6,iArray))

        test=trim(thisNames(iArray))//' is monotonic [no direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray)                              ,allowEqual=.true. )
        call Assert(char(test),isMonotonic,thisExpectations(7,iArray))

        test=trim(thisNames(iArray))//' is monotonic [increasing direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray),direction=directionIncreasing,allowEqual=.true. )
        call Assert(char(test),isMonotonic,thisExpectations(8,iArray))

        test=trim(thisNames(iArray))//' is monotonic [decreasing direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySet(:,iArray),direction=directionDecreasing,allowEqual=.true. )
        call Assert(char(test),isMonotonic,thisExpectations(9,iArray))

     end do
  end do
  
  ! Loop over all sets of integer arrays.
  do iArraySet=1,2
     
     ! Create pointers to the specific array set to process.
     select case (iArraySet)
     case (1)
        thisArraySetInteger => singleElementArraysInteger
        thisExpectations    => singleElementExpectations
        thisNames           => singleElementNames
     case (2)
        thisArraySetInteger => tenElementArraysInteger
        thisExpectations    => tenElementExpectations
        thisNames           => tenElementNames
     end select

     ! Loop over all arrays in this array set.
     do iArray=1,size(thisArraySetInteger,dim=2)

        test=trim(thisNames(iArray))//' [integer] is monotonic [no direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray)                                                 )
        call Assert(char(test),isMonotonic,thisExpectations(1,iArray))

        test=trim(thisNames(iArray))//' [integer] is monotonic [increasing direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray),direction=directionIncreasing                   )
        call Assert(char(test),isMonotonic,thisExpectations(2,iArray))

        test=trim(thisNames(iArray))//' [integer] is monotonic [decreasing direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray),direction=directionDecreasing                   )
        call Assert(char(test),isMonotonic,thisExpectations(3,iArray))

        test=trim(thisNames(iArray))//' [integer] is monotonic [no direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray)                              ,allowEqual=.false.)
        call Assert(char(test),isMonotonic,thisExpectations(4,iArray))

        test=trim(thisNames(iArray))//' [integer] is monotonic [increasing direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray),direction=directionIncreasing,allowEqual=.false.)
        call Assert(char(test),isMonotonic,thisExpectations(5,iArray))

        test=trim(thisNames(iArray))//' [integer] is monotonic [decreasing direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray),direction=directionDecreasing,allowEqual=.false.)
        call Assert(char(test),isMonotonic,thisExpectations(6,iArray))

        test=trim(thisNames(iArray))//' [integer] is monotonic [no direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray)                              ,allowEqual=.true. )
        call Assert(char(test),isMonotonic,thisExpectations(7,iArray))

        test=trim(thisNames(iArray))//' [integer] is monotonic [increasing direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray),direction=directionIncreasing,allowEqual=.true. )
        call Assert(char(test),isMonotonic,thisExpectations(8,iArray))

        test=trim(thisNames(iArray))//' [integer] is monotonic [decreasing direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(thisArraySetInteger(:,iArray),direction=directionDecreasing,allowEqual=.true. )
        call Assert(char(test),isMonotonic,thisExpectations(9,iArray))

     end do
  end do
  
  ! Test array reversal.
  call Assert('Reverse double precision array',Array_Reverse(doubleArray),doubleArrayReversed)
  call Assert('Reverse single precision array',Array_Reverse(realArray  ),realArrayReversed  )

  ! Test array cumulation.
  call Assert('Cumulate double precision array',Array_Cumulate(doubleArray),doubleArrayCumulated)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Array_Monotonicity
