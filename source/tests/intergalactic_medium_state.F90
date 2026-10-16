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

program Tests_Intergalactic_Medium_State
  !!{RST
  Tests calculations of the state of the intergalactic medium---specifically, of the optical depth to electron scattering, which
  is found by tabulating an integral over cosmic time.

  The universe used is Einstein-de Sitter and is taken to have reionized at an epoch earlier than any tabulated here, so that the
  electron fraction is constant over the whole of the tabulation and the optical depth has a closed form: with
  :math:`a=(t/t_0)^{2/3}` the integral of :math:`a^{-3}` from :math:`t` to :math:`t_0` is
  :math:`t_0^2(1/t-1/t_0)`. Both the optical depth and its inverse are checked against that form.

  The remaining tests are of the *determinism* of the tabulation: the optical depths returned must not depend on the order in
  which they were asked for. The tabulation is extended as earlier times are requested, so an object asked for its optical depths
  from the present day backwards builds its table in several steps, where one asked for them in the opposite order builds it in a
  single step. Both must return identical values, which requires that the times tabulated lie on an absolute lattice and that the
  values already computed be carried over exactly when the tabulation is extended.
  !!}
  use :: Cosmology_Functions             , only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters            , only : cosmologyParametersSimple
  use :: Display                         , only : displayVerbositySet           , verbosityLevelStandard
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Intergalactic_Medium_State      , only : intergalacticMediumStateSimple
  use :: Numerical_Constants_Astronomical, only : gigaYear                      , heliumByMassPrimordial, hydrogenByMassPrimordial, massSolar, &
       &                                          megaParsec
  use :: Numerical_Constants_Atomic      , only : atomicMassHelium              , atomicMassHydrogen    , atomicMassUnit
  use :: Numerical_Constants_Physical    , only : speedLight                    , thomsonCrossSection
  use :: Unit_Tests                      , only : Assert                        , Unit_Tests_Begin_Group, Unit_Tests_End_Group    , Unit_Tests_Finish
  implicit none
  double precision                                , dimension(6), parameter :: redshift=[0.1d0,0.5d0,1.0d0,3.0d0,9.0d0,29.0d0]
  double precision                                , dimension(6)            :: time                               , opticalDepthAscending              , &
       &                                                                       opticalDepthDescending             , opticalDepthAnalytic               , &
       &                                                                       timeInverse
  type            (cosmologyParametersSimple     )                          :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda)                          :: cosmologyFunctions_
  ! Three identically-constructed objects: one asked for its optical depths from the earliest time to the latest, one from the
  ! latest to the earliest, and one asked to invert them.
  type            (intergalacticMediumStateSimple)                          :: intergalacticMediumStateAscending_, intergalacticMediumStateDescending_, &
       &                                                                       intergalacticMediumStateInverse_
  character       (len=1024                      )                          :: message
  integer                                                                   :: i
  double precision                                                          :: timePresent                       , opticalDepthNormalization

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Initialize event hooks.
  call eventsHooksInitialize()
  ! Construct an Einstein-de Sitter model.
  cosmologyParameters_               =cosmologyParametersSimple     (                                                                  &
       &                                                             OmegaMatter               = 1.00d0                              , &
       &                                                             OmegaBaryon               = 0.05d0                              , &
       &                                                             OmegaDarkEnergy           = 0.00d0                              , &
       &                                                             temperatureCMB            = 2.78d0                              , &
       &                                                             HubbleConstant            =70.00d0                                &
       &                                                            )
  cosmologyFunctions_                =cosmologyFunctionsMatterLambda(                                                                  &
       &                                                             cosmologyParameters_      =cosmologyParameters_                   &
       &                                                            )
  ! Reionization is placed at an epoch far earlier than any time tabulated below, so that the intergalactic medium is fully
  ! ionized throughout the range tabulated and the optical depth therefore has the closed form used here.
  intergalacticMediumStateAscending_ =intergalacticMediumStateSimple(                                                                  &
       &                                                             reionizationRedshift      = 1.00d3                              , &
       &                                                             reionizationTemperature   = 1.00d4                              , &
       &                                                             preReionizationTemperature= 1.00d4                              , &
       &                                                             cosmologyFunctions_       =cosmologyFunctions_                  , &
       &                                                             cosmologyParameters_      =cosmologyParameters_                   &
       &                                                            )
  intergalacticMediumStateDescending_=intergalacticMediumStateSimple(                                                                  &
       &                                                             reionizationRedshift      = 1.00d3                              , &
       &                                                             reionizationTemperature   = 1.00d4                              , &
       &                                                             preReionizationTemperature= 1.00d4                              , &
       &                                                             cosmologyFunctions_       =cosmologyFunctions_                  , &
       &                                                             cosmologyParameters_      =cosmologyParameters_                   &
       &                                                            )
  intergalacticMediumStateInverse_   =intergalacticMediumStateSimple(                                                                  &
       &                                                             reionizationRedshift      = 1.00d3                              , &
       &                                                             reionizationTemperature   = 1.00d4                              , &
       &                                                             preReionizationTemperature= 1.00d4                              , &
       &                                                             cosmologyFunctions_       =cosmologyFunctions_                  , &
       &                                                             cosmologyParameters_      =cosmologyParameters_                   &
       &                                                            )
  ! Find the times corresponding to our redshifts, and the analytic optical depths at those times.
  timePresent              =cosmologyFunctions_%cosmicTime(1.0d0)
  opticalDepthNormalization=+speedLight                             &
       &                    *gigaYear                               &
       &                    *thomsonCrossSection                    &
       &                    *massSolar                              &
       &                    /atomicMassUnit                         &
       &                    /megaParsec         **3                 &
       &                    *cosmologyParameters_%OmegaBaryon    () &
       &                    *cosmologyParameters_%densityCritical() &
       &                    *(                                      &
       &                      +      hydrogenByMassPrimordial       &
       &                      /atomicMassHydrogen                   &
       &                      +2.0d0*heliumByMassPrimordial         &
       &                      /atomicMassHelium                     &
       &                     )
  do i=1,size(redshift)
     time                (i)=+cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift(i)))
     opticalDepthAnalytic(i)=+opticalDepthNormalization &
          &                  *timePresent          **2  &
          &                  *(                         &
          &                    +1.0d0/time       (i)    &
          &                    -1.0d0/timePresent       &
          &                   )
  end do
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Intergalactic medium state")
  ! Evaluate the optical depth to electron scattering, in each direction of time. The redshifts are in increasing order, so the
  ! first object is asked for the earliest time first, and builds its tabulation in a single step, while the second is asked for
  ! the latest time first, and so extends its tabulation repeatedly.
  do i=size(redshift),1,-1
     opticalDepthAscending (i)=intergalacticMediumStateAscending_ %electronScatteringOpticalDepth(time(i),assumeFullyIonized=.true.)
  end do
  do i=1,size(redshift)
     opticalDepthDescending(i)=intergalacticMediumStateDescending_%electronScatteringOpticalDepth(time(i),assumeFullyIonized=.true.)
  end do
  do i=1,size(redshift)
     write (message,'(a,f6.1,a)') "optical depth to electron scattering [z = ",redshift(i),"]"
     call Assert(trim(message),opticalDepthAscending(i),opticalDepthAnalytic(i),relTol=1.0d-3)
  end do
  call Assert("optical depth is independent of the order in which times are requested",opticalDepthAscending,opticalDepthDescending,absTol=0.0d0)
  ! Invert the optical depths, which uses a second pair of tabulations built from the first. The times are requested in the
  ! opposite order to the optical depths themselves, so this object reaches its final range differently again.
  do i=1,size(redshift)
     timeInverse(i)=intergalacticMediumStateInverse_%electronScatteringTime(opticalDepthAnalytic(i),assumeFullyIonized=.true.)
  end do
  do i=1,size(redshift)
     write (message,'(a,f6.1,a)') "time at which a given optical depth is reached [z = ",redshift(i),"]"
     call Assert(trim(message),timeInverse(i),time(i),relTol=1.0d-3)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Intergalactic_Medium_State
