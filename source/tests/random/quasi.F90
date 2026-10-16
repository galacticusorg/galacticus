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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program to test quasi-random number functions.
!!}

program Test_Quasi_Random
  !!{RST
  Tests that quasi-random number functions work.
  !!}
  use :: Display                         , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Quasi_Random_Sequences, only : gsl_qrng_sobol     , quasiRandomNumberGenerator
  use :: Unit_Tests                      , only : Assert             , Unit_Tests_Begin_Group    , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                                     , parameter                               :: countDimensions       =40, countPoints          =1024
  type            (quasiRandomNumberGenerator)                                          :: quasiRandomSequence      , quasiRandomSequence3D     , &
       &                                                                                   quasiRandomSequence40D
  double precision                            , dimension(                7          ) :: r                         , rSobol
  double precision                            , dimension(3              ,7          ) :: r3                        , rSobol3
  double precision                            , dimension(countDimensions,countPoints) :: design
  integer                                     , dimension(                countPoints) :: countPerBin
  logical                                     , dimension(countDimensions            ) :: isBalanced
  integer                                                                              :: i                         , j                         , &
       &                                                                                  k

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  quasiRandomSequence=quasiRandomNumberGenerator(gsl_qrng_sobol)
  call Unit_Tests_Begin_Group("quasi-random number sequences")
  do i=1,7
     r(i)=quasiRandomSequence%get()
  end do
  rSobol=[1.0d0/2.0d0,3.0d0/4.0d0,1.0d0/4.0d0,3.0d0/8.0d0,7.0d0/8.0d0,5.0d0/8.0d0,1.0d0/8.0d0]
  call Assert('Sobol sequence (1D)',r,rSobol,absTol=1.0d-6)
  ! Three-dimensional sequence. The first two dimensions are common to all Sobol constructions (and agree with, for example,
  ! SciPy's implementation, after its leading origin). The third depends on the choice of direction numbers, and is GSL's.
  ! Each sequence is held in its own variable: assigning a new sequence onto a used one double-finalizes its resource manager
  ! (GCC PR 110626).
  quasiRandomSequence3D=quasiRandomNumberGenerator(gsl_qrng_sobol,countDimensions=3)
  call Assert('Sobol sequence (3D): dimension',quasiRandomSequence3D%dimensions(),3)
  do i=1,7
     call quasiRandomSequence3D%getVector(r3(:,i))
  end do
  rSobol3=reshape(                                       &
       &          [                                      &
       &           1.0d0/2.0d0,1.0d0/2.0d0,1.0d0/2.0d0,  &
       &           3.0d0/4.0d0,1.0d0/4.0d0,3.0d0/4.0d0,  &
       &           1.0d0/4.0d0,3.0d0/4.0d0,1.0d0/4.0d0,  &
       &           3.0d0/8.0d0,3.0d0/8.0d0,5.0d0/8.0d0,  &
       &           7.0d0/8.0d0,7.0d0/8.0d0,1.0d0/8.0d0,  &
       &           5.0d0/8.0d0,1.0d0/8.0d0,3.0d0/8.0d0,  &
       &           1.0d0/8.0d0,5.0d0/8.0d0,7.0d0/8.0d0   &
       &          ]                                    , &
       &          [3,7]                                  &
       &         )
  call Assert('Sobol sequence (3D): first points',r3,rSobol3,absTol=1.0d-6)
  ! Balance of a 2^m-point design in the maximum supported dimension. GSL omits the origin, so the origin plus the first
  ! 2^m-1 points form a (t,m,s)-net: in every one-dimensional projection, each interval [k/2^m,(k+1)/2^m) holds exactly one
  ! point. This is the property that makes 2^m-point Sobol designs balanced.
  quasiRandomSequence40D=quasiRandomNumberGenerator(gsl_qrng_sobol,countDimensions=countDimensions)
  design(:,1)=0.0d0
  do i=2,countPoints
     call quasiRandomSequence40D%getVector(design(:,i))
  end do
  do j=1,countDimensions
     countPerBin=0
     do i=1,countPoints
        k             =min(int(design(j,i)*dble(countPoints))+1,countPoints)
        countPerBin(k)=countPerBin(k)+1
     end do
     isBalanced(j)=all(countPerBin == 1)
  end do
  call Assert('Sobol sequence (40D): 1024-point design is balanced in every dimension',all(isBalanced),.true.)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Quasi_Random
