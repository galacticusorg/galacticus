!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a program to test mathematical distributions.

program Test_Math_Distributions
  !% Tests of mathematical distributions.
  use, intrinsic :: ISO_C_Binding                      , only : c_long
  use            :: Galacticus_Display                 , only : Galacticus_Verbosity_Level_Set, verbosityStandard
  use            :: Input_Parameters                   , only : inputParameters
  use            :: Math_Distributions_Poisson_Binomial, only : Poisson_Binomial_Distribution , Poisson_Binomial_Distribution_Mean_Pairs
  use            :: Numerical_Random_Numbers           , only : randomNumberGeneratorGSL
  use            :: Statistics_Distributions           , only : distributionFunction1DGamma
  use            :: Unit_Tests                         , only : Assert                        , Unit_Tests_Begin_Group                  , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                             , dimension(  10) :: p                , x           , y
  integer                                      , dimension(0:10) :: trials
  double precision                             , dimension(0:10) :: Pk               , PkMonteCarlo, errorMonteCarlo
  integer                                      , parameter       :: trialCount=100000
  type            (distributionFunction1DGamma)                  :: distributionGamma_
  integer                                                        :: i                , j           , k
  type            (randomNumberGeneratorGSL   )                  :: prng
  type            (inputParameters            )                  :: parameters

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Math: distributions")
  ! Initialize parameters.
  parameters=inputParameters()
  call parameters%markGlobal()
  ! Test Poisson binomial distribution.
  p=[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0]
  do k=0,10
     Pk(k)=Poisson_Binomial_Distribution(k,p)
  end do
  call Assert("Poisson binomial: normalization",sum(Pk),1.0d0,absTol=1.0d-6)
  ! Generate a Monte Carlo realization of the Poisson binomial distribution for comparison.
  prng  =randomNumberGeneratorGSL(923_c_long)
  trials=0
  do i=1,trialCount
     k=0
     do j=1,size(p)
        if (prng%uniformSample() <= p(j)) k=k+1
     end do
     trials(k)=trials(k)+1
  end do
  PkMonteCarlo   =     dble(trials) /dble(trialCount)
  errorMonteCarlo=sqrt(dble(trials))/dble(trialCount)
  where(trials == 0)
     errorMonteCarlo=100.0d0
  end where
  call Assert("Poisson binomial: distribution",all(abs(Pk-PkMonteCarlo)/errorMonteCarlo < 3.0d0),.true.)
  ! Test Poisson binomial mean pairs in the trivial case.
  p=1.0d0
  call Assert("Poisson binomial: mean pairs",Poisson_Binomial_Distribution_Mean_Pairs(p),90.0d0,absTol=1.0d-4)

  ! Gamma distribution.
  distributionGamma_=distributionFunction1DGamma(2.0d0,1.2d0,prng,limitLower=0.3d0,limitUpper=6.0d0)
  x=[0.3d0,0.8d0,1.3d0,1.8d0,2.3d0,2.8d0,3.3d0,3.8d0,5.0d0,5.9d0]
  do i=1,10
     p(i)=distributionGamma_%cumulative(x(i))
     y(i)=distributionGamma_%inverse(p(i))
  end do
  call Assert("Gamma: inversion",x,y,relTol=1.0d-4)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Math_Distributions
