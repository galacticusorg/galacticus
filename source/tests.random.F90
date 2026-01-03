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
Contains a program to test random number functions.
!!}

program Test_Random
  !!{
  Tests that random number functions work.
  !!}
  use            :: Display                 , only : displayVerbositySet     , verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding           , only : c_long
  use            :: Input_Parameters        , only : inputParameters
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorGSL
  use            :: Unit_Tests              , only : Assert                  , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                                   , parameter             :: sampleCount    =10000000                                  , limitCount    =6
  double precision                          , dimension(limitCount) :: upperLimit     =[1.0d-5,1.0d-4,1.0d-3,1.0d-2,1.0d-1,0.5d0]
  double precision                          , dimension(limitCount) :: frequency                                                 , frequencyError  , &
       &                                                               deviation
  integer                                   , dimension(limitCount) :: upperLimitCount
  type            (randomNumberGeneratorGSL)                        :: randomSequence
  double precision                                                  :: x
  integer                                                           :: i                                                         , j
  type            (inputParameters         )                        :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize parameters.
  parameters=inputParameters()
  ! Begin unit tests.
  randomSequence=randomNumberGeneratorGSL(295_c_long)
  call Unit_Tests_Begin_Group("random numbers")
  upperLimitCount=0
  do i=1,sampleCount
     x=randomSequence%uniformSample()
     do j=1,size(upperLimit)
        if (x <= upperLimit(j)) upperLimitCount(j)=upperLimitCount(j)+1
     end do
  end do
  frequency     =dble(upperLimitCount)/dble(sampleCount)
  frequencyError=sqrt(dble(sampleCount)*upperLimit*(1.0d0-upperLimit))/dble(sampleCount)
  deviation     =abs(frequency-upperLimit)/frequencyError
  call Assert('uniformity of random numbers',all(deviation < 3.0d0),.true.)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Random
