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
Contains a program that tests primordial power spectrum calculations.
!!}

program Tests_Power_Spectrum_Primordial
  !!{
  Tests primordial power spectrum calculations.
  !!}
  use :: Display                 , only : displayVerbositySet                     , verbosityLevelStandard
  use :: Power_Spectra_Primordial, only : powerSpectrumPrimordialPiecewisePowerLaw
  use :: Unit_Tests              , only : Assert                                  , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (powerSpectrumPrimordialPiecewisePowerLaw), pointer   :: powerSpectrumPrimordial_
  double precision                                          , parameter :: tolerance               =1.0d-6

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Primordial power spectra")
  ! Test that index is correct.
  allocate(powerSpectrumPrimordial_)
  !![
  <referenceConstruct object="powerSpectrumPrimordial_">
   <constructor>
     powerSpectrumPrimordialPiecewisePowerLaw(                                                &amp;
      &amp;                                   index_             =[0.900d0,0.950d0, 1.000d0], &amp;
      &amp;                                   running            =[0.000d0,0.000d0, 0.000d0], &amp;
      &amp;                                   runningRunning     =[0.000d0,0.000d0, 0.000d0], &amp;
      &amp;                                   wavenumberReference=[        1.000d0,10.000d0]  &amp;
      &amp;                                  )
   </constructor>
  </referenceConstruct>
  !!]
  call Assert(                                                                                                                                                                         &
       &      'dlnPₚ/dlnk(k) = nₛ(k)'                                                                                                                                                , &
       &      [powerSpectrumPrimordial_%logarithmicDerivative(0.50d0),powerSpectrumPrimordial_%logarithmicDerivative(3.00d0),powerSpectrumPrimordial_%logarithmicDerivative(20.00d0)], &
       &      [                                               0.90d0 ,                                               0.95d0 ,                                                1.00d0 ], &
       &      relTol=1.0d-6                                                                                                                                                            &
       &     )
  !![
  <objectDestructor name="powerSpectrumPrimordial_"/>
  !!]
  ! Test that power spectrum is continuous.
  allocate(powerSpectrumPrimordial_)
  !![
  <referenceConstruct object="powerSpectrumPrimordial_">
   <constructor>
     powerSpectrumPrimordialPiecewisePowerLaw(                                                &amp;
      &amp;                                   index_             =[0.900d0,0.950d0, 1.000d0], &amp;
      &amp;                                   running            =[0.020d0,0.010d0, 0.000d0], &amp;
      &amp;                                   runningRunning     =[0.005d0,0.003d0, 0.000d0], &amp;
      &amp;                                   wavenumberReference=[        1.000d0,10.000d0]  &amp;
      &amp;                                  )
   </constructor>
  </referenceConstruct>
  !!]
  call Assert(                                                                                                                    &
       &      'Pₚ(k) continuous'                                                                                                , &
       &      [powerSpectrumPrimordial_%power(1.0d0*(1.0d0-tolerance)),powerSpectrumPrimordial_%power(10.0d0*(1.0d0-tolerance))], &
       &      [powerSpectrumPrimordial_%power(1.0d0*(1.0d0+tolerance)),powerSpectrumPrimordial_%power(10.0d0*(1.0d0+tolerance))], &
       &      relTol=2.0d0*tolerance                                                                                              &
       &     )
  !![
  <objectDestructor name="powerSpectrumPrimordial_"/>
  !!]
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Power_Spectrum_Primordial
