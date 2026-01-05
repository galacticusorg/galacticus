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
Contains a program which tests Gaunt factor functions.
!!}

program Test_Gaunt_Factors
  !!{
  Tests Gaunt factor functions.
  !!}
  use :: Atomic_Ionization_Potentials  , only : atomicIonizationPotentialVerner
  use :: Atomic_Radiation_Gaunt_Factors, only : gauntFactorVanHoof2014
  use :: Display                       , only : displayVerbositySet            , verbosityLevelStandard
  use :: Input_Parameters              , only : inputParameters
  use :: Unit_Tests                    , only : Assert                         , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (atomicIonizationPotentialVerner)                          :: atomicIonizationPotential_
  type            (gauntFactorVanHoof2014         )                          :: gauntFactorVanHoof2014_
  double precision                                 , dimension(5)            :: gauntFactors
  ! Target Gaunt factors computed using http://data.nublado.org/gauntff/interpolate2.f
  double precision                                 , dimension(5), parameter :: temperatures              =[1.0000000000000000d4,1.0000000000000000d5,1.0000000000000000d6,1.0000000000000000d7,1.0000000000000000d8]
  double precision                                 , dimension(5), parameter :: gauntFactorsTarget        =[1.2651040532869235d0,1.4144240304642641d0,1.4005277750508995d0,1.2421867138942060d0,1.1500622299841357d0]
  integer                                                                    :: i
  type            (inputParameters                )                          :: parameters

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize parameters.
  parameters=inputParameters()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("van Hoof et al. (2014) fitting function:")
  atomicIonizationPotential_=atomicIonizationPotentialVerner(                          )
  gauntFactorVanHoof2014_   =gauntFactorVanHoof2014         (atomicIonizationPotential_)
  do i=1,size(temperatures)
     gauntFactors(i)=gauntFactorVanHoof2014_%total(1,1,temperatures(i))
  end do
  call Assert('total',gauntFactors,gauntFactorsTarget,relTol=1.0d-6)
  ! End unit tests.
  call Unit_Tests_End_Group       ()
  call Unit_Tests_Finish          ()

end program Test_Gaunt_Factors
