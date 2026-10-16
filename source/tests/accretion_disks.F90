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
Contains a program to test accretion disk functions.
!!}

!+ Contributions to this file made by: Andrew Benson, Claude.

program Test_Accretion_Disks
  !!{RST
  Tests of accretion disk functions. The expected :term:`ADAF` jet power efficiencies were computed using an independent
  implementation of the :cite:t:`benson_maximum_2009` model (``bensonBabul2009_galacticusAlpha.py`` in the `galacticusDevTools
  <https://github.com/galacticusorg/galacticusDevTools>`_ repository) with the :math:`\alpha(j)` fitting functions used by
  Galacticus and the Kerr metric factor :math:`\mathcal{A}=1+j^2/r^2+2j^2/r^3`.
  !!}
  use :: Accretion_Disks             , only : accretionDisksADAF            , adafEnergyPureADAF    , adafFieldEnhancementExponential, adafRadiativeEfficiencyTypeFixed, &
          &                                   adafViscosityFit
  use :: Display                     , only : displayVerbositySet           , verbosityLevelStandard
  use :: Galacticus_Nodes            , only : nodeComponentBlackHoleStandard
  use :: Numerical_Constants_Physical, only : speedLight
  use :: Numerical_Constants_Prefixes, only : kilo
  use :: Unit_Tests                  , only : Assert                        , Unit_Tests_Begin_Group, Unit_Tests_End_Group           , Unit_Tests_Finish               , &
          &                                   compareEquals
  implicit none
  double precision                                , dimension(6) :: spin            =[0.0000d+0,0.2000d+0,0.4000d+0,0.6000d+0,0.8000d+0,0.9500d+0]
  double precision                                , dimension(6) :: jetPowerExpected=[2.9928d-3,4.0715d-3,7.5336d-3,1.9602d-2,5.6359d-2,4.2037d-1]
  double precision                                , dimension(6) :: jetPower
  type            (accretionDisksADAF            )               :: accretionDisk
  type            (nodeComponentBlackHoleStandard)               :: blackHole
  integer                                                        :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Accretion disks")
  ! Build an accretion disk.
  accretionDisk=accretionDisksADAF(                                                       &
       &                                                adafEnergyPureADAF              , &
       &                                                adafFieldEnhancementExponential , &
       &                                                adafRadiativeEfficiencyTypeFixed, &
       &                                                adafViscosityFit                , &
       &                           efficiencyJetMaximum=2.000d0                         , &
       &                           efficiencyRadiation =0.010d0                         , &
       &                           adiabaticIndex      =1.444d0                           &
       &                          )
  ! Compute jet power efficiency.
  do i=1,size(spin)
     call blackHole%massSet(  1.0d0)
     call blackHole%spinSet(spin(i))
     jetPower(i)=accretionDisk%powerJet(blackHole,accretionRateMass=1.0d0)/(speedLight/kilo)**2
  end do
  ! Test jet power efficiency.
  call Assert("Jet power efficiency vs. spin [ADAF]" ,jetPower,jetPowerExpected,compareEquals,relTol=1.0d-3)
  ! End the unit testing.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
end program Test_Accretion_Disks
