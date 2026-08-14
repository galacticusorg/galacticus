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

program Tests_Conditional_Mass_Functions
  !!{RST
  Tests the conditional mass function of :cite:t:`behroozi_comprehensive_2010`.

  That class tabulates the halo mass corresponding to each mass, and reverses the tabulation to obtain the mass corresponding to
  each halo mass. The range tabulated is extended whenever a halo mass outside it is requested, so the first test is of the
  *determinism* of that tabulation: the conditional mass function returned for a given halo mass must not depend on which other
  halo masses have been asked for. Three objects are used, differing only in what they are asked and in what order: one is asked
  for the smallest halo mass first, so that its tabulation is extended to its full extent immediately, one is asked in the
  opposite order, so that it is extended repeatedly, and a third is built afresh for each halo mass, so that its tabulation is
  extended only as far as that one halo mass requires. All three must agree exactly.

  The tabulation is then checked against the fitting function which it tabulates, which is reimplemented here: at the mass for
  which the median relation gives a particular halo mass, precisely half of halos of that mass should host a central galaxy more
  massive than it.
  !!}
  use :: Conditional_Mass_Functions, only : conditionalMassFunctionBehroozi2010, haloModelGalaxyTypeCentral
  use :: Display                   , only : displayVerbositySet                , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Unit_Tests                , only : Assert                             , Unit_Tests_Begin_Group    , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  ! Halo masses at which the conditional mass function is evaluated. The seeded tabulation spans masses of 10⁸ to 10¹³ M☉, which
  ! corresponds to halo masses above 4.9×10¹⁰ M☉, so the smaller of these require the tabulation to be extended - by five decades
  ! in the case of the smallest.
  double precision                                     , dimension(7), parameter :: massHalo                        =[1.0d8,1.0d9,1.0d10,1.0d11,1.0d12,1.0d13,1.0d14]
  ! Parameters of the fitting function, from the fit of Leauthaud et al. (2011).
  double precision                                                   , parameter :: alphaSatellite                  = 1.0000d0, log10M1                          =12.5200d0, &
       &                                                                            log10Mstar0                     =10.9160d0, beta                             = 0.4570d0, &
       &                                                                            delta                           = 0.5666d0, gamma                            = 1.5300d0, &
       &                                                                            sigmaLogMstar                   = 0.2060d0, BCut                             = 1.4700d0, &
       &                                                                            BSatellite                      =10.6200d0, betaCut                          =-0.1300d0, &
       &                                                                            betaSatellite                   = 0.8590d0
  double precision                                     , dimension(7)            :: massFunctionAscending                     , massFunctionDescending                     , &
       &                                                                            massFunctionSingle
  type            (conditionalMassFunctionBehroozi2010)                          :: conditionalMassFunctionAscending          , conditionalMassFunctionDescending
  ! One object per halo mass, each constructed once and asked a single question - these must not be re-assigned in a loop, as the
  ! assignment would finalize the previous object.
  type            (conditionalMassFunctionBehroozi2010), dimension(7)            :: conditionalMassFunctionSingle
  character       (len=1024                           )                          :: message
  integer                                                                        :: i
  double precision                                                               :: mass                                      , numberCentrals

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  call eventsHooksInitialize()
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Conditional mass functions")
  ! Evaluate the conditional mass function in each direction of halo mass, and for an object built afresh for each halo mass. A
  ! fixed mass is used throughout, so that only the halo mass - and so only the extent of the tabulation - differs between them.
  conditionalMassFunctionAscending =conditionalMassFunctionBehroozi2010(alphaSatellite,log10M1,log10Mstar0,beta,delta,gamma,sigmaLogMstar,BCut,BSatellite,betaCut,betaSatellite)
  conditionalMassFunctionDescending=conditionalMassFunctionBehroozi2010(alphaSatellite,log10M1,log10Mstar0,beta,delta,gamma,sigmaLogMstar,BCut,BSatellite,betaCut,betaSatellite)
  do i=1,size(massHalo)
     massFunctionAscending (i)=conditionalMassFunctionAscending %massFunction(massHalo(i),1.0d10)
  end do
  do i=size(massHalo),1,-1
     massFunctionDescending(i)=conditionalMassFunctionDescending%massFunction(massHalo(i),1.0d10)
  end do
  do i=1,size(massHalo)
     conditionalMassFunctionSingle(i)=conditionalMassFunctionBehroozi2010(alphaSatellite,log10M1,log10Mstar0,beta,delta,gamma,sigmaLogMstar,BCut,BSatellite,betaCut,betaSatellite)
     massFunctionSingle           (i)=conditionalMassFunctionSingle(i)%massFunction(massHalo(i),1.0d10)
  end do
  call Assert("conditional mass function is independent of the order in which halo masses are requested",massFunctionAscending,massFunctionDescending,absTol=0.0d0)
  call Assert("conditional mass function is independent of the range of masses tabulated"               ,massFunctionAscending,massFunctionSingle    ,absTol=0.0d0)
  ! Check the tabulated median relation against the fitting function which it tabulates. At the mass for which the median
  ! relation gives a halo mass, half of the halos of that mass should host a more massive central galaxy. The tolerance is set by
  ! the interpolation error of the tabulation, which is made at ten points per decade of mass: toward high halo masses the median
  ! relation steepens sharply, so that a decade of mass spans an ever larger range of halo mass and the reversed tabulation is
  ! correspondingly more coarsely sampled. That is a property of the tabulation's resolution, not of the range over which it is
  ! made, and is unchanged by the pinning of that range - the tests above are what constrain the latter.
  do i=3,size(massHalo)
     mass          =massStellar(massHalo(i))
     numberCentrals=conditionalMassFunctionAscending%massFunction(massHalo(i),mass,haloModelGalaxyTypeCentral)
     write (message,'(a,e9.3,a)') "half of halos host a central galaxy above the median mass [M = ",massHalo(i)," M☉]"
     call Assert(trim(message),numberCentrals,0.5d0,relTol=5.0d-2)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  double precision function massHaloMedian(mass)
    !!{RST
    The median halo mass vs. mass relation of :cite:t:`behroozi_comprehensive_2010`, reimplemented here to test the tabulation of
    it. The tests are limited to halo masses for which the logarithmic halo mass remains below the transition value above which
    the class being tested allows it to grow only logarithmically, so no such transition is needed here.
    !!}
    implicit none
    double precision, intent(in   ) :: mass
    double precision                :: massRelative

    massRelative  =+mass                            &
         &         /10.0d0**log10Mstar0
    massHaloMedian=+10.0d0**(                       &
         &                   +log10M1               &
         &                   +beta                  &
         &                   *log10(massRelative)   &
         &                   +massRelative**delta   &
         &                   /(                     &
         &                     +1.0d0               &
         &                     +1.0d0               &
         &                     /massRelative**gamma &
         &                    )                     &
         &                   -0.5d0                 &
         &                  )
    return
  end function massHaloMedian

  double precision function massStellar(massHaloTarget)
    !!{RST
    Return the mass for which the median relation gives the halo mass ``massHaloTarget``, found by bisection.
    !!}
    implicit none
    double precision, intent(in   ) :: massHaloTarget
    double precision                :: massLogarithmicLow , massLogarithmicHigh, &
         &                             massLogarithmicMean
    integer                         :: i

    massLogarithmicLow =-6.0d0
    massLogarithmicHigh=+14.0d0
    do i=1,200
       massLogarithmicMean=0.5d0*(massLogarithmicLow+massLogarithmicHigh)
       if (massHaloMedian(10.0d0**massLogarithmicMean) < massHaloTarget) then
          massLogarithmicLow =massLogarithmicMean
       else
          massLogarithmicHigh=massLogarithmicMean
       end if
    end do
    massStellar=10.0d0**(0.5d0*(massLogarithmicLow+massLogarithmicHigh))
    return
  end function massStellar

end program Tests_Conditional_Mass_Functions
