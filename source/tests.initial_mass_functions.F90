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
Contains a program to test stellar initial mass functions.
!!}

program Test_Initial_Mass_Functions
  !!{
  Tests of stellar initial mass functions.
  !!}
  use :: Display                                   , only : displayVerbositySet             , verbosityLevelStandard
  use :: Numerical_Integration2                    , only : integratorCompositeTrapezoidal1D
  use :: NUmerical_Constants_Astronomical          , only : metallicitySolar
  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionBPASS        , initialMassFunctionBaugh2005TopHeavy, initialMassFunctionChabrier2001        , initialMassFunctionClass            , &
          &                                                 initialMassFunctionKennicutt1983, initialMassFunctionKroupa2001       , initialMassFunctionMillerScalo1979     , initialMassFunctionPiecewisePowerLaw, &
          &                                                 initialMassFunctionSalpeter1955 , initialMassFunctionScalo1986
  use :: Supernovae_Type_Ia                        , only : supernovaeTypeIaNagashima2005   , supernovaeTypeIaPowerLawDTD         , supernovaeTypeIaPowerLawDTDDifferential, supernovaeTypeIaClass
  use :: Stellar_Astrophysics                      , only : stellarAstrophysicsFile
  use :: Unit_Tests                                , only : Assert                          , Unit_Tests_Begin_Group              , Unit_Tests_End_Group                   , Unit_Tests_Finish
  implicit none
  class           (initialMassFunctionClass               ), pointer      :: imf
  type            (initialMassFunctionChabrier2001        ), target       :: imfChabrier2001
  type            (initialMassFunctionPiecewisePowerLaw   ), target       :: imfPiecewisePowerLaw
  type            (initialMassFunctionSalpeter1955        ), target       :: imfSalpeter1955
  type            (initialMassFunctionBPASS               ), target       :: imfBPASS
  type            (initialMassFunctionBaugh2005TopHeavy   ), target       :: imfBaugh2005TopHeavy
  type            (initialMassFunctionKennicutt1983       ), target       :: imfKennicutt1983
  type            (initialMassFunctionKroupa2001          ), target       :: imfKroupa2001
  type            (initialMassFunctionMillerScalo1979     ), target       :: imfMillerScalo1979
  type            (initialMassFunctionScalo1986           ), target       :: imfScalo1986
  type            (stellarAstrophysicsFile                )               :: stellarAstrophysics_
  type            (supernovaeTypeIaNagashima2005          ), target       :: supernovaeTypeIaNagashima2005_
  type            (supernovaeTypeIaPowerLawDTD            ), target       :: supernovaeTypeIaPowerLawDTD_
  type            (supernovaeTypeIaPowerLawDTDDifferential), target       :: supernovaeTypeIaPowerLawDTDDifferential_
  class           (supernovaeTypeIaClass                  ), pointer      :: supernovaeTypeIa_
  type            (integratorCompositeTrapezoidal1D       )               :: integrator_
  ! Values of the cumulative number of type Ia SNe read from Figure 6 of Nagashima et al (2005; MNRAS; 363; 31;
  ! https://ui.adsabs.harvard.edu/abs/2005MNRAS.363L..31N) at 0.1, 1, and 10Gyr, and for a power-law delay time distribution.
  double precision                                      , dimension(3) :: ageTypeIa                     =[                                          &
       &                                                                                                  1.000000000000000d-1,                     &
       &                                                                                                  1.000000000000000d+0,                     &
       &                                                                                                  1.000000000000000d+1                      &
       &                                                                                                 ]                    ,                     &
       &                                                                  numberTypeIaSNeNagashima2005  =[                                          &
       &                                                                                                  6.000000000000000d-6,                     &
       &                                                                                                  6.600000000000000d-4,                     &
       &                                                                                                  2.200000000000000d-3                      &
       &                                                                                                 ]                    ,                     &
       &                                                                  numberTypeIaSNePowerLawDTD    =[                                          &
       &                                                                                                  2.334828201481421d-4,                     &
       &                                                                                                  7.581754850034585d-4,                     &
       &                                                                                                  1.204761370427590d-3                      &
       &                                                                                                 ]
  double precision                                                     :: massInInitialMassFunction                           , numberTypeIaSNe   , &
       &                                                                  massInitialMinimum                                  , massInitialMaximum
  integer                                                              :: i
  character       (len=12                              )               :: label
  
  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Stellar initial mass functions")
  call integrator_%initialize  (24           )
  call integrator_%toleranceSet(1.0d-7,1.0d-7)
  call integrator_%integrandSet(initialMassFunctionIntegrand)
  call Unit_Tests_Begin_Group("Normalization")
  imfChabrier2001     =initialMassFunctionChabrier2001     (                                              &
       &                                                    massLower         =+  0.10d0                , &
       &                                                    massUpper         =+125.00d0                , &
       &                                                    massTransition    =+  1.00d0                , &
       &                                                    massCharacteristic=+  0.08d0                , &
       &                                                    exponent          =-  2.30d0                , &
       &                                                    sigma             =+  0.69d0                  &
       &                                                   )
  ! This piecewise IMF is matched to the "Kennicutt IMF" used by Nagashima et al (2005; MNRAS; 363; 31;
  ! https://ui.adsabs.harvard.edu/abs/2005MNRAS.363L..31N) in their Type Ia SNe calculations.
  imfPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                              &
       &                                                    mass              =[+0.15d0,+1.0d0,+120.0d0], &
       &                                                    exponent          =[-1.40d0,-2.5d0         ]  &
       &                                                   )
  imfSalpeter1955     =initialMassFunctionSalpeter1955     (                                              &
       &                                                   )
  imfBPASS            =initialMassFunctionBPASS            (                                              &
       &                                                   )
  imfBaugh2005TopHeavy=initialMassFunctionBaugh2005TopHeavy(                                              &
       &                                                   )
  imfKennicutt1983    =initialMassFunctionKennicutt1983    (                                              &
       &                                                   )
  imfKroupa2001       =initialMassFunctionKroupa2001       (                                              &
       &                                                   )
  imfMillerScalo1979  =initialMassFunctionMillerScalo1979  (                                              &
       &                                                   )
  imfScalo1986        =initialMassFunctionScalo1986        (                                              &
       &                                                   )
  imf                       => imfChabrier2001
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('Chabrier (2001)'              ,massInInitialMassFunction,1.0d0,relTol=1.0d-6)
  imf                       => imfPiecewisePowerLaw
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('piecewise power-law'          ,massInInitialMassFunction,1.0d0,relTol=1.0d-6)
  imf                       => imfSalpeter1955
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('Salpeter (1955)'              ,massInInitialMassFunction,1.0d0,relTol=1.0d-6)
  imf                       => imfBPASS
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('BPASS'                        ,massInInitialMassFunction,1.0d0,relTol=1.0d-6)

  imf                       => imfBaugh2005TopHeavy
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('Baugh et al. (2005) top heavy',massInInitialMassFunction,1.0d0,relTol=1.0d-6)
  imf                       => imfKennicutt1983
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('Kennicutt (1983)'             ,massInInitialMassFunction,1.0d0,relTol=1.0d-6)
  imf                       => imfKroupa2001
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('Kroupa (2001)'                ,massInInitialMassFunction,1.0d0,relTol=1.0d-6)
  imf                       => imfMillerScalo1979
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('Miller & Scalo (1979)'        ,massInInitialMassFunction,1.0d0,relTol=1.0d-6)
  imf                       => imfScalo1986
  massInInitialMassFunction =  integrator_%evaluate(                   &
       &                                            imf%massMinimum(), &
       &                                            imf%massMaximum()  &
       &                                           )
  call Assert('Scalo (1986)'                 ,massInInitialMassFunction,1.0d0,relTol=1.0d-6)
  call Unit_Tests_End_Group()
  call Unit_Tests_Begin_Group('Type Ia SNe')
  call Unit_Tests_Begin_Group('Nagashima et al. (2005)')
  stellarAstrophysics_           =  stellarAstrophysicsFile      ('%DATASTATICPATH%/stellarAstrophysics/stellarPropertiesPortinariChiosiBressan1998.xml')
  supernovaeTypeIaNagashima2005_ =  supernovaeTypeIaNagashima2005(stellarAstrophysics_)
  supernovaeTypeIa_              => supernovaeTypeIaNagashima2005_
  call integrator_%toleranceSet(1.0d-6,1.0d-6)
  call integrator_%integrandSet(numberTypeIaSNeIntegrand)
  do i=1,size(ageTypeIa)
     call supernovaeTypeIaNagashima2005_%massInitialRange(imfPiecewisePowerLaw,ageTypeIa(i),metallicitySolar,massInitialMinimum,massInitialMaximum)
     massInitialMinimum=max(massInitialMinimum,imfPiecewisePowerLaw%massMinimum())
     massInitialMaximum=min(massInitialMaximum,imfPiecewisePowerLaw%massMaximum())
     numberTypeIaSNe=integrator_%evaluate(                    &
          &                               massInitialMinimum, &
          &                               massInitialMaximum  &
          &                              )  
     write (label,'(f4.1)') ageTypeIa(i)
     ! The tolerance for this test is (very) low. Nagashima et al. (2005) do not specify what they use for M(t) - the mass of a
     ! star leaving the main sequence at time t). Therefore, the best we can do is an approximate comparison.
     call Assert('Cumulative number of Type Ia SNe at age '//trim(label)//'Gyr',numberTypeIaSNe,numberTypeIaSNeNagashima2005(i),relTol=0.5d0)
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Begin_Group('Power law delay time distribution')
  supernovaeTypeIaPowerLawDTD_ =  supernovaeTypeIaPowerLawDTD(40.0d-3,-1.07d0,0.21d-3)
  supernovaeTypeIa_            => supernovaeTypeIaPowerLawDTD_
  call integrator_%toleranceSet(1.0d-6,1.0d-6)
  call integrator_%integrandSet(numberTypeIaSNeIntegrand)
  do i=1,size(ageTypeIa)
     call supernovaeTypeIaPowerLawDTD_%massInitialRange(imfPiecewisePowerLaw,ageTypeIa(i),metallicitySolar,massInitialMinimum,massInitialMaximum)
     massInitialMinimum=max(massInitialMinimum,imfPiecewisePowerLaw%massMinimum())
     massInitialMaximum=min(massInitialMaximum,imfPiecewisePowerLaw%massMaximum())
     numberTypeIaSNe=integrator_%evaluate(                    &
          &                               massInitialMinimum, &
          &                               massInitialMaximum  &
          &                              )  
     write (label,'(f4.1)') ageTypeIa(i)
     call Assert('Cumulative number of Type Ia SNe at age '//trim(label)//'Gyr',numberTypeIaSNe,numberTypeIaSNePowerLawDTD(i),relTol=1.0d-3)
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Begin_Group('Power law delay time distribution (differential)')
  supernovaeTypeIaPowerLawDTDDifferential_ =  supernovaeTypeIaPowerLawDTDDifferential(40.0d-3,-1.07d0,0.21d-3)
  supernovaeTypeIa_                        => supernovaeTypeIaPowerLawDTDDifferential_
  call integrator_%toleranceSet(1.0d-6,1.0d-6)
  call integrator_%integrandSet(numberTypeIaSNeIntegrand)
  do i=1,size(ageTypeIa)
     call supernovaeTypeIaPowerLawDTDDifferential_%massInitialRange(imfPiecewisePowerLaw,ageTypeIa(i),metallicitySolar,massInitialMinimum,massInitialMaximum)
     massInitialMinimum=max(massInitialMinimum,imfPiecewisePowerLaw%massMinimum())
     massInitialMaximum=min(massInitialMaximum,imfPiecewisePowerLaw%massMaximum())
     numberTypeIaSNe=integrator_%evaluate(                    &
          &                               massInitialMinimum, &
          &                               massInitialMaximum  &
          &                              )  
     write (label,'(f4.1)') ageTypeIa(i)
     call Assert('Cumulative number of Type Ia SNe at age '//trim(label)//'Gyr',numberTypeIaSNe,numberTypeIaSNePowerLawDTD(i),relTol=1.0d-3)
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_End_Group()
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

contains

  double precision function initialMassFunctionIntegrand(mass)
    !!{
    Integrand used to find the total mass in the initial mass function.
    !!}
    implicit none
    double precision, intent(in   ) :: mass

    initialMassFunctionIntegrand=+        mass  &
         &                       *imf%phi(mass)
    return
  end function initialMassFunctionIntegrand

  double precision function numberTypeIaSNeIntegrand(massSecondary)
    !!{
    Integrand used to find the cumulative number of Type Ia SNe.
    !!}
    implicit none
    double precision, intent(in   ) :: massSecondary

    numberTypeIaSNeIntegrand=+supernovaeTypeIa_   %number(imfPiecewisePowerLaw,massSecondary,ageTypeIa(i),metallicitySolar)  &
         &                   *imfPiecewisePowerLaw%phi   (                     massSecondary                              )
    return
  end function numberTypeIaSNeIntegrand

end program Test_Initial_Mass_Functions


