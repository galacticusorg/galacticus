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

!!{RST
Contains a program to test the accuracy of the series-approximation branch of the NFW kinematics distribution against the exact (dilogarithm-based) form. Both branches evaluate the 1D velocity dispersion of a self-gravitating NFW profile via the Jeans equation, but the ``useSeriesApproximation=.true.`` path avoids the dilogarithm by expanding around five anchor points in scale-free radius. This test quantifies the relative error between the two paths across the radius range encountered in subhalo evolution.
!!}

program Test_Kinematic_Distributions_NFW
  !!{RST
  Tests the series-approximation branch of the NFW kinematics distribution against the exact form.
  !!}
  use :: Coordinates       , only : assignment(=)         , coordinateSpherical
  use :: Display           , only : displayMessage        , displayVerbositySet   , verbosityLevelStandard
  use :: Events_Hooks      , only : eventsHooksInitialize
  use :: IO_HDF5           , only : ioHDF5AccessInitialize
  use :: Mass_Distributions, only : massDistributionClass , massDistributionNFW   , kinematicsDistributionClass, kinematicsDistributionNFW
  use :: Numerical_Ranges  , only : Make_Range            , rangeTypeLogarithmic
  use :: Unit_Tests        , only : Assert                , Unit_Tests_Begin_Group, Unit_Tests_End_Group       , Unit_Tests_Finish
  implicit none
  integer                                                  , parameter              :: radiusCount             =301
  double precision                                         , parameter              :: radiusScaleFreeMinimum  =1.0d-4, radiusScaleFreeMaximum=1.0d+4
  class           (massDistributionClass      )                        , allocatable :: massDistribution_
  class           (kinematicsDistributionClass)                        , pointer     :: kinematicsExact_       , kinematicsSeries_
  type            (coordinateSpherical        )                                      :: coordinates
  double precision                                         , dimension(radiusCount)  :: radiusScaleFree        , velocityDispersionExact, &
       &                                                                                velocityDispersionSeries, errorRelative
  double precision                                                                   :: errorRelativeMaximum   , errorRelativeRMS       , &
       &                                                                                radiusErrorMaximum
  integer                                                                            :: i                      , iErrorMaximum
  character       (len=128                    )                                      :: message

  ! Set verbosity level.
  call displayVerbositySet  (verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize(                      )
  ! Initialize HDF5 lock.
  call ioHDF5AccessInitialize(                     )
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("NFW kinematics distribution: series vs exact")

  ! Construct a dimensionless self-gravitating NFW mass distribution. Both branches of velocityDispersion1D reduce to
  ! scale-free expressions in r/r_s, so this single distribution suffices.
  allocate(massDistributionNFW :: massDistribution_)
  select type (massDistribution_)
  type is (massDistributionNFW)
     massDistribution_=massDistributionNFW(scaleLength=1.0d0,densityNormalization=1.0d0,dimensionless=.true.)
  end select

  ! Construct the two kinematics distributions: one using the exact (dilogarithm-based) form, one using the series
  ! approximation.
  allocate(kinematicsDistributionNFW :: kinematicsExact_ )
  select type (kinematicsExact_ )
  type is (kinematicsDistributionNFW)
     kinematicsExact_ =kinematicsDistributionNFW(useSeriesApproximation=.false.)
  end select
  allocate(kinematicsDistributionNFW :: kinematicsSeries_)
  select type (kinematicsSeries_)
  type is (kinematicsDistributionNFW)
     kinematicsSeries_=kinematicsDistributionNFW(useSeriesApproximation=.true. )
  end select

  ! Build a logarithmic grid in scale-free radius covering the range likely to be encountered in subhalo evolution.
  radiusScaleFree=Make_Range(radiusScaleFreeMinimum,radiusScaleFreeMaximum,radiusCount,rangeTypeLogarithmic)

  ! Evaluate the 1D velocity dispersion from both forms. The third argument to velocityDispersion1D must be the same target
  ! as the second so that the self-gravitating NFW analytic branch is taken inside the kinematicsDistributionNFW.
  do i=1,radiusCount
     coordinates                =[radiusScaleFree(i),0.0d0,0.0d0]
     velocityDispersionExact (i)=kinematicsExact_ %velocityDispersion1D(coordinates,massDistribution_,massDistribution_)
     velocityDispersionSeries(i)=kinematicsSeries_%velocityDispersion1D(coordinates,massDistribution_,massDistribution_)
  end do

  ! Compute relative error per point and aggregate statistics.
  errorRelative       =abs(velocityDispersionSeries-velocityDispersionExact)/velocityDispersionExact
  errorRelativeMaximum=maxval(errorRelative       )
  iErrorMaximum       =maxloc(errorRelative,dim =1)
  radiusErrorMaximum  =       radiusScaleFree(iErrorMaximum)
  errorRelativeRMS    =sqrt(sum(errorRelative**2)/dble(radiusCount))

  ! Report worst-case statistics so the developer can judge whether the series approximation is acceptable as the production
  ! default. These messages are emitted regardless of the assertion outcome.
  write (message,'(a,e13.6,a,e13.6)') '  Maximum relative error: ',errorRelativeMaximum,' at r/r_s = ',radiusErrorMaximum
  call displayMessage(trim(message))
  write (message,'(a,e13.6           )') '  RMS     relative error: ',errorRelativeRMS
  call displayMessage(trim(message))

  ! Assert that the series approximation matches the exact form to better than 10⁻⁵ everywhere on the test grid. The observed
  ! maximum is ~5×10⁻⁶ near r/r_s≈0.3 (the join between the expansions around 0 and around 1/2), so this threshold gives a
  ! small headroom for compiler/optimization variation while still catching genuine accuracy regressions.
  call Assert("Series approximation matches exact (relative)",errorRelative,spread(0.0d0,1,radiusCount),absTol=1.0d-5)

  ! Clean up.
  deallocate(massDistribution_)
  deallocate(kinematicsExact_ )
  deallocate(kinematicsSeries_)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Test_Kinematic_Distributions_NFW
