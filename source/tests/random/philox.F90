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
Contains a program to test the Philox counter-based random number generator.
!!}

program Test_Random_Philox
  !!{RST
  Tests the Philox counter-based random number generator. The raw sequence is compared with that of ``numpy.random.Philox``, and
  the distributions of uniform, standard normal, Poisson, and bounded integer deviates are checked.
  !!}
  use            :: Display                 , only : displayVerbositySet        , verbosityLevelStandard
  use, intrinsic :: ISO_C_Binding           , only : c_int64_t                  , c_long
  use            :: Numerical_Random_Numbers, only : randomNumberGeneratorPhilox
  use            :: Unit_Tests              , only : Assert                     , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                                      , parameter     :: countSamples         =4000000
  type            (randomNumberGeneratorPhilox)                :: generatorSequence                   , generatorCounter, &
       &                                                          generatorDistribution
  double precision                             , dimension(4)  :: uniform
  double precision                             , dimension(16) :: normalFirst                         , normalSecond
  double precision                                             :: x                                   , sum1            , &
       &                                                          sum2                                , sum4            , &
       &                                                          mean                                , variance
  integer                                                      :: i                                   , k
  integer         (c_long                     ), dimension(7)  :: countBins
  double precision                             , dimension(2)  :: meansPoisson         =[3.0d0,50.0d0]

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Philox random numbers")
  ! Compare the sequence with that of numpy.random.Philox. The expected uniform deviates are the upper 53 bits of the raw 64-bit
  ! output of numpy.random.Philox(key=[219,0],counter=[0,0,0,0]), multiplied by 2⁻⁵³, so are exactly representable.
  generatorSequence=randomNumberGeneratorPhilox(219_c_long)
  do i=1,4
     uniform(i)=generatorSequence%uniformSample()
  end do
  call Assert('sequence matches numpy.random.Philox'                          ,uniform,[0.99805879664803070d0,0.11367866618742983d0,0.70629886148656410d0,0.94630274906420550d0],relTol=1.0d-15)
  ! Positioning by counter, compared with numpy.random.Philox(key=[12345,0],counter=[0,5,2,0]), and with a negative label,
  ! numpy.random.Philox(key=[12345,0],counter=[0,2⁶⁴-3,6,0]), whose bit pattern is that of the signed integer -3.
  generatorCounter=randomNumberGeneratorPhilox(12345_c_long)
  call generatorCounter%counterSet([5_c_int64_t,2_c_int64_t,0_c_int64_t])
  do i=1,4
     uniform(i)=generatorCounter%uniformSample()
  end do
  call Assert('positioned stream matches numpy.random.Philox'                 ,uniform,[0.66020962438660540d0,0.71925321135806510d0,0.30926866951845644d0,0.51535218389321100d0],relTol=1.0d-15)
  call generatorCounter%counterSet([-3_c_int64_t,6_c_int64_t,0_c_int64_t])
  do i=1,4
     uniform(i)=generatorCounter%uniformSample()
  end do
  call Assert('positioned stream (negative label) matches numpy.random.Philox',uniform,[0.48509410970785083d0,0.15738237599677818d0,0.23576092109238445d0,0.27509801803624334d0],relTol=1.0d-15)
  ! Repositioning to the same counter regenerates the same deviates, regardless of what was drawn in between (including a
  ! partially used block of random bits and a stored Gaussian deviate).
  call generatorCounter%counterSet([7_c_int64_t,1_c_int64_t,0_c_int64_t])
  do i=1,16
     normalFirst(i)=generatorCounter%standardNormalSample()
  end do
  x=generatorCounter%standardNormalSample()
  x=generatorCounter%uniformSample       ()
  call generatorCounter%counterSet([7_c_int64_t,1_c_int64_t,0_c_int64_t])
  do i=1,16
     normalSecond(i)=generatorCounter%standardNormalSample()
  end do
  call Assert('repositioned stream is reproducible',normalSecond,normalFirst,relTol=0.0d0)
  call generatorCounter%counterSet([7_c_int64_t,2_c_int64_t,0_c_int64_t])
  do i=1,16
     normalSecond(i)=generatorCounter%standardNormalSample()
  end do
  call Assert('streams with different labels differ',any(normalSecond /= normalFirst),.true.)
  ! Uniform deviates: mean 1/2 and variance 1/12. The tolerances are about 5σ for this number of samples.
  generatorDistribution=randomNumberGeneratorPhilox(8231_c_long)
  sum1=0.0d0
  sum2=0.0d0
  do i=1,countSamples
     x   =generatorDistribution%uniformSample()
     sum1=sum1+x
     sum2=sum2+x**2
  end do
  mean    =sum1/dble(countSamples)
  variance=sum2/dble(countSamples)-mean**2
  call Assert('uniform deviates: mean'    ,mean    ,1.0d0/ 2.0d0,absTol=5.0d0*sqrt(1.0d0/ 12.0d0/dble(countSamples)))
  call Assert('uniform deviates: variance',variance,1.0d0/12.0d0,absTol=5.0d0*sqrt(1.0d0/180.0d0/dble(countSamples)))
  ! Standard normal deviates: mean 0, variance 1, and kurtosis 3.
  sum1=0.0d0
  sum2=0.0d0
  sum4=0.0d0
  do i=1,countSamples
     x   =generatorDistribution%standardNormalSample()
     sum1=sum1+x
     sum2=sum2+x**2
     sum4=sum4+x**4
  end do
  call Assert('normal deviates: mean'    ,sum1/dble(countSamples),0.0d0,absTol=5.0d0*sqrt( 1.0d0/dble(countSamples)))
  call Assert('normal deviates: variance',sum2/dble(countSamples),1.0d0,absTol=5.0d0*sqrt( 2.0d0/dble(countSamples)))
  call Assert('normal deviates: kurtosis',sum4/dble(countSamples),3.0d0,absTol=5.0d0*sqrt(96.0d0/dble(countSamples)))
  ! Poisson deviates, by inversion (mean 3) and by transformed rejection (mean 50): the mean and variance both equal the mean.
  do k=1,size(meansPoisson)
     sum1=0.0d0
     sum2=0.0d0
     do i=1,countSamples
        x   =dble(generatorDistribution%poissonSample(meansPoisson(k)))
        sum1=sum1+x
        sum2=sum2+x**2
     end do
     mean    =sum1/dble(countSamples)
     variance=sum2/dble(countSamples)-mean**2
     call Assert('Poisson deviates: mean'    ,mean    ,meansPoisson(k),absTol=5.0d0*sqrt(meansPoisson(k)                              /dble(countSamples)))
     call Assert('Poisson deviates: variance',variance,meansPoisson(k),absTol=5.0d0*sqrt(meansPoisson(k)*(1.0d0+2.0d0*meansPoisson(k))/dble(countSamples)))
  end do
  ! Bounded integers: each value in [0,6] is equally likely.
  countBins=0_c_long
  do i=1,countSamples
     k           =int(generatorDistribution%sample(7_c_long))+1
     countBins(k)=countBins(k)+1_c_long
  end do
  call Assert('bounded integers: uniform over [0,n-1]',dble(countBins)/dble(countSamples),spread(1.0d0/7.0d0,1,7),absTol=5.0d0*sqrt(6.0d0/49.0d0/dble(countSamples)))
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Test_Random_Philox
