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
Contains a program to test the array functions.
!!}

program Test_Array_Monotonicity
  !!{
  Tests that array functions.
  !!}
  use, intrinsic :: ISO_C_Binding     , only : c_size_t
  use            :: Array_Utilities   , only : Array_Cumulate     , Array_Is_Monotonic      , Array_Reverse       , directionDecreasing, &
          &                                    directionIncreasing, operator(.intersection.), slice5Dto3D
  use            :: Display           , only : displayVerbositySet, verbosityLevelStandard
  use            :: ISO_Varying_String, only : assignment(=)      , char                    , varying_string
  use            :: Kind_Numbers      , only : kind_int8
  use            :: Unit_Tests        , only : Assert             , Unit_Tests_Begin_Group  , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                , dimension( 1,2      ), target     :: singleElementArrays       =reshape([1.23d0,-2.31d0],shape(singleElementArrays))
  integer         (kind=kind_int8), dimension( 1,2      ), target     :: singleElementArraysInteger=reshape([123,-231],shape(singleElementArraysInteger))
  logical                         , dimension( 9,2      ), target     :: singleElementExpectations =reshape([.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.],shape(singleElementExpectations))
  character(len=128),      target , dimension(   2      )             :: singleElementNames=    [                                                                              &
       &                                                                                          'Single element array (positive value)'                                      &
       &                                                                                         ,'Single element array (negative value)'                                      &
       &                                                                                        ]
  double precision                , dimension(10,6      ), target     :: tenElementArrays       =reshape([1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0,3.0d0,2.0d0,1.0d0,0.0d0,-1.0d0,-2.0d0,-3.0d0,-4.0d0,-5.0d0,-6.0d0,1.0d0,0.0d0,3.0d0,4.0d0,-3.0d0,5.0d0,1.0d0,8.0d0,2.0d0,3.0d0,1.0d0,1.0d0,1.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0,10.0d0,3.0d0,2.0d0,1.0d0,1.0d0,-1.0d0,-2.0d0,-3.0d0,-4.0d0,-5.0d0,-6.0d0,1.0d0,0.0d0,3.0d0,3.0d0,-3.0d0,5.0d0,1.0d0,8.0d0,2.0d0,3.0d0],shape(tenElementArrays))
  integer         (kind=kind_int8), dimension(10,6      ), target     :: tenElementArraysInteger=reshape([10,20,30,40,50,60,70,80,90,100,30,20,10,00,-10,-20,-30,-40,-50,-60,10,00,30,40,-30,50,10,80,20,30,10,10,10,40,50,60,70,80,90,100,30,20,10,10,-10,-20,-30,-40,-50,-60,10,00,30,30,-30,50,10,80,20,30],shape(tenElementArrays))
  logical                         , dimension( 9,6      ), target     :: tenElementExpectations =reshape([.true.,.true.,.false.,.true.,.true.,.false.,.true.,.true.,.false.,.true.,.false.,.true.,.true.,.false.,.true.,.true.,.false.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.true.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.true.,.false.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.],shape(tenElementExpectations))
  character(len=128),      target , dimension(   6      )             :: tenElementNames=       [                                                                              &
       &                                                                                          'Increasing array (no equalities)     '                                      &
       &                                                                                         ,'Decreasing array (no equalities)     '                                      &
       &                                                                                         ,'Non-monotinic array (no equalities)  '                                      &
       &                                                                                         ,'Increasing array (with equalities)   '                                      &
       &                                                                                         ,'Decreasing array (with equalities)   '                                      &
       &                                                                                         ,'Non-monotinic array (with equalities)'                                      &
       &                                                                                        ]
  double precision                , dimension(10        )              :: doubleArray         =[0.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,7.0d0,8.0d0,9.0d0]
  double precision                , dimension(10        )              :: doubleArrayReversed =[9.0d0,8.0d0,7.0d0,6.0d0,5.0d0,4.0d0,3.0d0,2.0d0,1.0d0,0.0d0]
  double precision                , dimension(10        )              :: doubleArrayCumulated=[0.0d0,1.0d0,3.0d0,6.0d0,10.0d0,15.0d0,21.0d0,28.0d0,36.0d0,45.0d0]
  real                            , dimension(10        )              :: realArray           =[0.0e0,1.0e0,2.0e0,3.0e0,4.0e0,5.0e0,6.0e0,7.0e0,8.0e0,9.0e0]
  real                            , dimension(10        )              :: realArrayReversed   =[9.0e0,8.0e0,7.0e0,6.0e0,5.0e0,4.0e0,3.0e0,2.0e0,1.0e0,0.0e0]
  double precision                , dimension( :,:,:    ), allocatable :: array3D                                                                                 , array3DTarget
  double precision                , dimension( :,:,:,:,:), allocatable :: array5D
  integer         (c_size_t      ), dimension( 2        )              :: sliceDimension                                                                          , sliceIndex
  double precision                , dimension( :,:      ), pointer     :: arraySet
  integer         (kind=kind_int8), dimension( :,:      ), pointer     :: arraySetInteger
  logical                         , dimension( :,:      ), pointer     :: expectations
  character       (len=128       ), dimension(   :      ), pointer     :: names
  type            (varying_string), dimension(   :      ), allocatable :: set1                                                                                    , set2     , &
       &                                                                  set3
  logical                                                              :: isMonotonic
  integer                                                              :: iArray                                                                                  , iArraySet, &
       &                                                                  index1                                                                                  , index2   , &
       &                                                                  index3                                                                                  , index4   , &
       &                                                                  index5
  type            (varying_string)                                     :: test

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Array functions")

  ! Loop over all sets of arrays.
  do iArraySet=1,2

     ! Create pointers to the specific array set to process.
     select case (iArraySet)
     case (1)
        arraySet     => singleElementArrays
        expectations => singleElementExpectations
        names        => singleElementNames
     case (2)
        arraySet     => tenElementArrays
        expectations => tenElementExpectations
        names        => tenElementNames
     end select

     ! Loop over all arrays in this array set.
     do iArray=1,size(arraySet,dim=2)

        ! Perform tests on this array for each of the nine permutations of options:
        !  direction=   increasing|decreasing|undefined
        !  allowEquals= true      |false     |undefined

        test=trim(names(iArray))//' is monotonic [no direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray)                                                 )
        call Assert(char(test),isMonotonic,expectations(1,iArray))

        test=trim(names(iArray))//' is monotonic [increasing direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray),direction=directionIncreasing                   )
        call Assert(char(test),isMonotonic,expectations(2,iArray))

        test=trim(names(iArray))//' is monotonic [decreasing direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray),direction=directionDecreasing                   )
        call Assert(char(test),isMonotonic,expectations(3,iArray))

        test=trim(names(iArray))//' is monotonic [no direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray)                              ,allowEqual=.false.)
        call Assert(char(test),isMonotonic,expectations(4,iArray))

        test=trim(names(iArray))//' is monotonic [increasing direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray),direction=directionIncreasing,allowEqual=.false.)
        call Assert(char(test),isMonotonic,expectations(5,iArray))

        test=trim(names(iArray))//' is monotonic [decreasing direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray),direction=directionDecreasing,allowEqual=.false.)
        call Assert(char(test),isMonotonic,expectations(6,iArray))

        test=trim(names(iArray))//' is monotonic [no direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray)                              ,allowEqual=.true. )
        call Assert(char(test),isMonotonic,expectations(7,iArray))

        test=trim(names(iArray))//' is monotonic [increasing direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray),direction=directionIncreasing,allowEqual=.true. )
        call Assert(char(test),isMonotonic,expectations(8,iArray))

        test=trim(names(iArray))//' is monotonic [decreasing direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(arraySet(:,iArray),direction=directionDecreasing,allowEqual=.true. )
        call Assert(char(test),isMonotonic,expectations(9,iArray))

     end do
  end do

  ! Loop over all sets of integer arrays.
  do iArraySet=1,2

     ! Create pointers to the specific array set to process.
     select case (iArraySet)
     case (1)
        arraySetInteger => singleElementArraysInteger
        expectations    => singleElementExpectations
        names           => singleElementNames
     case (2)
        arraySetInteger => tenElementArraysInteger
        expectations    => tenElementExpectations
        names           => tenElementNames
     end select

     ! Loop over all arrays in this array set.
     do iArray=1,size(arraySetInteger,dim=2)

        test=trim(names(iArray))//' [integer] is monotonic [no direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray)                                                 )
        call Assert(char(test),isMonotonic,expectations(1,iArray))

        test=trim(names(iArray))//' [integer] is monotonic [increasing direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray),direction=directionIncreasing                   )
        call Assert(char(test),isMonotonic,expectations(2,iArray))

        test=trim(names(iArray))//' [integer] is monotonic [decreasing direction asserted, no equality condition]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray),direction=directionDecreasing                   )
        call Assert(char(test),isMonotonic,expectations(3,iArray))

        test=trim(names(iArray))//' [integer] is monotonic [no direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray)                              ,allowEqual=.false.)
        call Assert(char(test),isMonotonic,expectations(4,iArray))

        test=trim(names(iArray))//' [integer] is monotonic [increasing direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray),direction=directionIncreasing,allowEqual=.false.)
        call Assert(char(test),isMonotonic,expectations(5,iArray))

        test=trim(names(iArray))//' [integer] is monotonic [decreasing direction asserted, equality disallowed]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray),direction=directionDecreasing,allowEqual=.false.)
        call Assert(char(test),isMonotonic,expectations(6,iArray))

        test=trim(names(iArray))//' [integer] is monotonic [no direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray)                              ,allowEqual=.true. )
        call Assert(char(test),isMonotonic,expectations(7,iArray))

        test=trim(names(iArray))//' [integer] is monotonic [increasing direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray),direction=directionIncreasing,allowEqual=.true. )
        call Assert(char(test),isMonotonic,expectations(8,iArray))

        test=trim(names(iArray))//' [integer] is monotonic [decreasing direction asserted, equality allowed]'
        isMonotonic=Array_Is_Monotonic(arraySetInteger(:,iArray),direction=directionDecreasing,allowEqual=.true. )
        call Assert(char(test),isMonotonic,expectations(9,iArray))

     end do
  end do

  ! Test array reversal.
  call Assert('Reverse double precision array',Array_Reverse(doubleArray),doubleArrayReversed)
  call Assert('Reverse single precision array',Array_Reverse(realArray  ),realArrayReversed  )

  ! Test array cumulation.
  call Assert('Cumulate double precision array',Array_Cumulate(doubleArray),doubleArrayCumulated)

  ! Test array intersection.
  allocate(set1(3))
  allocate(set2(3))
  allocate(set3(1))
  set1=['abc','def','ghi']
  set2=['123','abc','789']
  set3=['abc'            ]
  call Assert('Array intersection',set1.intersection.set2,set3)
  deallocate(set1,set2,set3)

  ! Test array slicing.
  allocate(array5D      (3,3,3,3,3))
  allocate(array3DTarget(3,  3,  3))
  do index1=1,3
     do index2=1,3
        do index3=1,3
           do index4=1,3
              do index5=1,3
                 array5D(index1,index2,index3,index4,index5)=dble(index1)+dble(index2)+dble(index3)+dble(index4)+dble(index5)
                 if (index4 == 2 .and. index2 == 1) &
                      & array3DTarget(index1,index3,index5)=array5D(index1,index2,index3,index4,index5)
              end do
           end do
        end do
     end do
  end do
  sliceDimension=[4,2]
  sliceIndex    =[2,1]
  array3D       =slice5Dto3D(array5D,sliceDimension,sliceIndex)
  call Assert("Array slicing",array3D,array3DTarget,absTol=1.0d-30)
  
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Array_Monotonicity
