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
Contains a program to test stellar initial mass functions.
!!}

program Test_Initial_Mass_Functions
  !!{RST
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
  ! A second Chabrier (2001) initial mass function, identical to the first but with the transition between the log-normal and
  ! power-law branches moved away from its default of 1 M☉, used to measure the continuity of the branches at that transition.
  type            (initialMassFunctionChabrier2001        ), target       :: imfChabrier2001Shifted
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
  ! Masses at which the Chabrier (2001) initial mass function is evaluated, spanning both branches, and the values of the
  ! initial mass function there, computed independently by `chabrier2001IMF.py`.
  double precision                                      , dimension(8) :: massChabrier                  =[1.5000000000d-01,3.0000000000d-01,5.0000000000d-01,8.0000000000d-01,1.5000000000d+00,5.0000000000d+00,2.0000000000d+01,1.0000000000d+02]
  double precision                                      , dimension(8) :: phiChabrierReference          =[5.1331586852d+00,1.9636241594d+00,8.5624117148d-01,3.6415207035d-01,9.2614484754d-02,5.8084171733d-03,2.3950788779d-04,5.9113790861d-06]
  double precision                                      , dimension(8) :: phiChabrier
  ! Mass ranges over which the cumulative number of stars formed per unit mass is evaluated, and the reference values.
  double precision                                      , dimension(4) :: massChabrierLower             =[1.0000000000d-01,1.0000000000d+00,8.0000000000d+00,1.0000000000d-01]
  double precision                                      , dimension(4) :: massChabrierUpper             =[1.0000000000d+00,8.0000000000d+00,1.2500000000d+02,1.2500000000d+02]
  double precision                                      , dimension(4) :: numberChabrierReference       =[1.2873108407d+00,1.6890157242d-01,1.1786082627d-02,1.4679984958d+00]
  double precision                                      , dimension(4) :: numberChabrier
  ! A transition mass away from the default, used to check that the two branches join continuously and that the initial mass
  ! function remains normalized to unit mass there.
  double precision                                      , parameter    :: massTransitionShifted         =2.0000000000d+00
  double precision                                                     :: continuityRatio                                     , massInInitialMassFunctionShifted
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
  ! Shape of the Chabrier (2001) initial mass function. The normalization group above integrates each initial mass function
  ! between its own mass limits and requires unit total mass. That is a self-consistency check: it is satisfied by any shape
  ! which the constructor happens to normalize, and so cannot detect an incorrect functional form. The assertions below instead
  ! compare the initial mass function itself, and its cumulative number, against values computed independently by the
  ! `chabrier2001IMF.py` script in the galacticusDevTools repository, which derives the normalization of each branch in closed
  ! form (the log-normal branch via error functions, the power-law branch analytically) rather than transcribing the
  ! expressions used here.
  call Unit_Tests_Begin_Group('Chabrier (2001) shape')
  do i=1,size(massChabrier)
     phiChabrier   (i)=imfChabrier2001%phi             (massChabrier     (i)                     )
  end do
  do i=1,size(massChabrierLower)
     numberChabrier(i)=imfChabrier2001%numberCumulative(massChabrierLower(i),massChabrierUpper(i))
  end do
  call Assert('initial mass function'        ,phiChabrier              ,phiChabrierReference   ,relTol=1.0d-6)
  call Assert('cumulative number'            ,numberChabrier           ,numberChabrierReference,relTol=1.0d-6)
  ! Continuity and normalization away from the default transition mass. Requiring the log-normal and power-law branches to join
  ! continuously at the transition mass M_t fixes the coefficient of the power-law branch to exp(...)/M_t^(1+α), and the
  ! mass integral used to normalize the initial mass function is computed with exactly that coefficient. Both assertions below
  ! therefore hold for any M_t.
  !
  ! These were added because they did not: the coefficient used when *evaluating* the power-law branch omitted the exponent,
  ! reading exp(...)/M_t. The two agree only at M_t = 1 M☉, which is the default, so nothing shipped was affected - but for any
  ! other transition mass the branches were discontinuous by a factor M_t^α, and, because the normalization had been
  ! computed with the other coefficient, the initial mass function was not normalized to unit mass at all: the total came to
  ! 3.77 at M_t = 0.5 M☉ and 0.62 at M_t = 2 M☉. The normalization group above cannot detect this, as it exercises only the
  ! default transition mass.
  imfChabrier2001Shifted=initialMassFunctionChabrier2001(                                                  &
       &                                                 massLower         =imfChabrier2001%massMinimum(), &
       &                                                 massTransition    =massTransitionShifted        , &
       &                                                 massUpper         =imfChabrier2001%massMaximum(), &
       &                                                 exponent          =-2.3d0                       , &
       &                                                 massCharacteristic= 0.08d0                      , &
       &                                                 sigma             = 0.69d0                        &
       &                                                )
  continuityRatio=+imfChabrier2001Shifted%phi(massTransitionShifted*(1.0d0+1.0d-12)) &
       &          /imfChabrier2001Shifted%phi(massTransitionShifted*(1.0d0-1.0d-12))
  call Assert('branches join continuously at a shifted transition mass',continuityRatio,1.0d0,relTol=1.0d-6)
  imf                             => imfChabrier2001Shifted
  call integrator_%toleranceSet(1.0d-7,1.0d-7)
  call integrator_%integrandSet(initialMassFunctionIntegrand)
  massInInitialMassFunctionShifted=  integrator_%evaluate(                   &
       &                                                  imf%massMinimum(), &
       &                                                  imf%massMaximum()  &
       &                                                 )
  call Assert('unit mass at a shifted transition mass'               ,massInInitialMassFunctionShifted,1.0d0,relTol=1.0d-6)
  call Unit_Tests_End_Group()
  call Unit_Tests_Begin_Group('Type Ia SNe')
  call Unit_Tests_Begin_Group('Nagashima et al. (2005)')
  stellarAstrophysics_           =  stellarAstrophysicsFile      (                     '%DATASTATICPATH%/stellarAstrophysics/stellarPropertiesPortinariChiosiBressan1998.xml')
  supernovaeTypeIaNagashima2005_ =  supernovaeTypeIaNagashima2005(stellarAstrophysics_,'%DATASTATICPATH%/stellarAstrophysics/Supernovae_Type_Ia_Yields.xml'                  )
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
  supernovaeTypeIaPowerLawDTD_ =  supernovaeTypeIaPowerLawDTD(40.0d-3,-1.07d0,0.21d-3,'%DATASTATICPATH%/stellarAstrophysics/Supernovae_Type_Ia_Yields.xml')
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
  supernovaeTypeIaPowerLawDTDDifferential_ =  supernovaeTypeIaPowerLawDTDDifferential(40.0d-3,-1.07d0,0.21d-3,'%DATASTATICPATH%/stellarAstrophysics/Supernovae_Type_Ia_Yields.xml')
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
    !!{RST
    Integrand used to find the total mass in the initial mass function.
    !!}
    implicit none
    double precision, intent(in   ) :: mass

    initialMassFunctionIntegrand=+        mass  &
         &                       *imf%phi(mass)
    return
  end function initialMassFunctionIntegrand

  double precision function numberTypeIaSNeIntegrand(massSecondary)
    !!{RST
    Integrand used to find the cumulative number of Type Ia SNe.
    !!}
    implicit none
    double precision, intent(in   ) :: massSecondary

    numberTypeIaSNeIntegrand=+supernovaeTypeIa_   %number(imfPiecewisePowerLaw,massSecondary,ageTypeIa(i),metallicitySolar)  &
         &                   *imfPiecewisePowerLaw%phi   (                     massSecondary                              )
    return
  end function numberTypeIaSNeIntegrand

end program Test_Initial_Mass_Functions


