!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionBPASS        , initialMassFunctionBaugh2005TopHeavy, initialMassFunctionChabrier2001   , initialMassFunctionClass            , &
          &                                                 initialMassFunctionKennicutt1983, initialMassFunctionKroupa2001       , initialMassFunctionMillerScalo1979, initialMassFunctionPiecewisePowerLaw, &
          &                                                 initialMassFunctionSalpeter1955 , initialMassFunctionScalo1986
  use :: Unit_Tests                                , only : Assert                          , Unit_Tests_Begin_Group              , Unit_Tests_End_Group              , Unit_Tests_Finish
  implicit none
  class           (initialMassFunctionClass            ), pointer :: imf
  type            (initialMassFunctionChabrier2001     ), target  :: imfChabrier2001
  type            (initialMassFunctionPiecewisePowerLaw), target  :: imfPiecewisePowerLaw
  type            (initialMassFunctionSalpeter1955     ), target  :: imfSalpeter1955
  type            (initialMassFunctionBPASS            ), target  :: imfBPASS
  type            (initialMassFunctionBaugh2005TopHeavy), target  :: imfBaugh2005TopHeavy
  type            (initialMassFunctionKennicutt1983    ), target  :: imfKennicutt1983
  type            (initialMassFunctionKroupa2001       ), target  :: imfKroupa2001
  type            (initialMassFunctionMillerScalo1979  ), target  :: imfMillerScalo1979
  type            (initialMassFunctionScalo1986        ), target  :: imfScalo1986
  type            (integratorCompositeTrapezoidal1D    )          :: integrator_
  double precision                                                :: massInInitialMassFunction

  call displayVerbositySet(verbosityLevelStandard)
  call Unit_Tests_Begin_Group("Stellar initial mass functions")
  call integrator_%initialize  (24           )
  call integrator_%toleranceSet(1.0d-7,1.0d-7)
  call integrator_%integrandSet(initialMassFunctionIntegrand)
  call Unit_Tests_Begin_Group("Normalization")
  imfChabrier2001     =initialMassFunctionChabrier2001     (                                       &
       &                                                    massLower         =+  0.10d0         , &
       &                                                    massUpper         =+125.00d0         , &
       &                                                    massTransition    =+  1.00d0         , &
       &                                                    massCharacteristic=+  0.08d0         , &
       &                                                    exponent          =-  2.30d0         , &
       &                                                    sigma             =+  0.69d0           &
       &                                                   )
  imfPiecewisePowerLaw=initialMassFunctionPiecewisePowerLaw(                                       &
       &                                                    mass              =[+0.10d0,+125.0d0], &
       &                                                    exponent          =[-2.30d0         ]  &
       &                                                   )
  imfSalpeter1955     =initialMassFunctionSalpeter1955     (                                       &
       &                                                   )
  imfBPASS            =initialMassFunctionBPASS            (                                       &
       &                                                   )
  imfBaugh2005TopHeavy=initialMassFunctionBaugh2005TopHeavy(                                       &
       &                                                   )
  imfKennicutt1983    =initialMassFunctionKennicutt1983    (                                       &
       &                                                   )
  imfKroupa2001       =initialMassFunctionKroupa2001       (                                       &
       &                                                   )
  imfMillerScalo1979  =initialMassFunctionMillerScalo1979  (                                       &
       &                                                   )
  imfScalo1986        =initialMassFunctionScalo1986        (                                       &
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

end program Test_Initial_Mass_Functions


