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

program Tests_Cosmic_Age
  !!{
  Tests cosmic age calculations for various Universes. Ages calculated using Python
  \href{http://www.astro.ucla.edu/~wright/CC.python}{implementation} of Ned Wright's cosmology calculator.
  !!}
  use :: Cosmology_Functions     , only : cosmologyFunctions       , cosmologyFunctionsClass, cosmologyFunctionsMatterDarkEnergy, cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters    , only : cosmologyParametersSimple, hubbleUnitsTime
  use :: Display                 , only : displayVerbositySet      , verbosityLevelStandard
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                   , Unit_Tests_Begin_Group , Unit_Tests_End_Group              , Unit_Tests_Finish
  implicit none
  double precision                                    , dimension(8), parameter :: redshift                                      =[0.0000000000d+0,1.0000000000d+0,3.0000000000d+0,9.0000000000d+0,3.0000000000d+1,1.0000000000d+2,3.0000000000d+2,1.0000000000d+3]
  double precision                                    , dimension(8), parameter :: ageEdS                                        =[0.0651868207d+0,0.0230470058d+0,0.0081483367d+0,0.0020613766d+0,0.0003776671d+0,0.0000642171d+0,0.0000124804d+0,0.0000020571d+0]
  double precision                                    , dimension(8), parameter :: ageOpen                                       =[0.0790841462d+0,0.0327062977d+0,0.0128686687d+0,0.0035287732d+0,0.0006745641d+0,0.0001164482d+0,0.0000227374d+0,0.0000037552d+0]
  double precision                                    , dimension(8), parameter :: ageCosmologicalConstant                       =[0.0942699818d+0,0.0402619685d+0,0.0147878493d+0,0.0037621019d+0,0.0006895264d+0,0.0001172509d+0,0.0000227901d+0,0.0000037578d+0]
  double precision                                    , dimension(8), parameter :: ageClosed                                     =[3.4369560000d-2,8.6126540000d-3,2.7752830000d-3,6.7037040000d-4,1.2048780000d-4,2.0363050000d-5,3.9509380000d-6,6.5106720000d-7]
  type            (cosmologyParametersSimple         )                          :: cosmologyParametersClosed                                                                                                                                                       , cosmologyParametersCosmologicalConstant         , &
       &                                                                           cosmologyParametersOpen                                                                                                                                                         , cosmologyParametersEdS
  type            (cosmologyFunctionsMatterLambda    )                          :: cosmologyFunctionsCosmologicalConstant                                                                                                                                          , cosmologyFunctionsOpen                          , &
       &                                                                           cosmologyFunctionsEdS
  type            (cosmologyFunctionsMatterDarkEnergy)                          :: cosmologyFunctionsDarkEnergyClosed                                                                                                                                              , cosmologyFunctionsDarkEnergyCosmologicalConstant, &
       &                                                                           cosmologyFunctionsDarkEnergyOmegaMinusOneThird
  character       (len=1024                          )                          :: message
  integer                                                                       :: iExpansion
  double precision                                                              :: age                                                                                                                                                                             , expansionFactor                                 , &
       &                                                                           expansionFactorSymmetric                                                                                                                                                        , timeTurnaround

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cosmic age")
  ! Define cosmological parameters.
  !![
  <referenceConstruct object="cosmologyParametersEdS"                 >
   <constructor>
    cosmologyParametersSimple(                             &amp;
     &amp;                    OmegaMatter         =1.00d0, &amp;
     &amp;                    OmegaBaryon         =0.00d0, &amp;
     &amp;                    OmegaDarkEnergy     =0.00d0, &amp;
     &amp;                    temperatureCMB      =2.78d0, &amp;
     &amp;                    HubbleConstant      =1.00d4  &amp;
     &amp;                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParametersOpen"                >
   <constructor>
    cosmologyParametersSimple(                             &amp;
     &amp;                    OmegaMatter         =0.30d0, &amp;
     &amp;                    OmegaBaryon         =0.00d0, &amp;
     &amp;                    OmegaDarkEnergy     =0.00d0, &amp;
     &amp;                    temperatureCMB      =2.78d0, &amp;
     &amp;                    HubbleConstant      =1.00d4  &amp;
     &amp;                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParametersClosed">
   <constructor>
    cosmologyParametersSimple(                             &amp;
     &amp;                    OmegaMatter         =1.00d1, &amp;
     &amp;                    OmegaBaryon         =0.00d0, &amp;
     &amp;                    OmegaDarkEnergy     =0.00d0, &amp;
     &amp;                    temperatureCMB      =2.78d0, &amp;
     &amp;                    HubbleConstant      =1.00d4  &amp;
     &amp;                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParametersCosmologicalConstant">
   <constructor>
    cosmologyParametersSimple(                             &amp;
     &amp;                    OmegaMatter         =0.30d0, &amp;
     &amp;                    OmegaBaryon         =0.00d0, &amp;
     &amp;                    OmegaDarkEnergy     =0.70d0, &amp;
     &amp;                    temperatureCMB      =2.78d0, &amp;
     &amp;                    HubbleConstant      =1.00d4  &amp;
     &amp;                   )
   </constructor>
  </referenceConstruct>
  !!]
  ! Build cosmological function objects.
  !![
  <referenceConstruct object="cosmologyFunctionsEdS"                           >
   <constructor>
    cosmologyFunctionsMatterLambda    (                                                                     &amp;
     &amp;                             cosmologyParameters_       =cosmologyParametersEdS                   &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsOpen"                          >
   <constructor>
    cosmologyFunctionsMatterLambda    (                                                                     &amp;
     &amp;                             cosmologyParameters_       =cosmologyParametersOpen                  &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsCosmologicalConstant"          >
   <constructor>
    cosmologyFunctionsMatterLambda    (                                                                     &amp;
     &amp;                             cosmologyParameters_       =cosmologyParametersCosmologicalConstant  &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsDarkEnergyClosed"              >
   <constructor>
    cosmologyFunctionsMatterDarkEnergy(                                                                     &amp;
     &amp;                             cosmologyParameters_       =cosmologyParametersClosed              , &amp;
     &amp;                             darkEnergyEquationOfStateW0=+0.0d0                                 , &amp;
     &amp;                             darkEnergyEquationOfStateW1=+0.0d0                                   &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsDarkEnergyCosmologicalConstant">
   <constructor>
    cosmologyFunctionsMatterDarkEnergy(                                                                     &amp;
     &amp;                             cosmologyParameters_       =cosmologyParametersCosmologicalConstant, &amp;
     &amp;                             darkEnergyEquationOfStateW0=-1.0d0                                 , &amp;
     &amp;                             darkEnergyEquationOfStateW1=+0.0d0                                   &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctionsDarkEnergyOmegaMinusOneThird"  >
   <constructor>
    cosmologyFunctionsMatterDarkEnergy(                                                                     &amp;
     &amp;                             cosmologyParameters_       =cosmologyParametersCosmologicalConstant, &amp;
     &amp;                             darkEnergyEquationOfStateW0=-1.0d0/3.0d0                           , &amp;
     &amp;                             darkEnergyEquationOfStateW1=+0.0d0                                   &amp;
     &amp;                            )
   </constructor>
  </referenceConstruct>
  !!]
  ! Compute the time of maximum expansion for the Universe. In this simple, Ωₘ=10, Ω_DE=0 Universe this is analytically
  ! calculable and equals:
  timeTurnaround=(5.0d0/27.0d0)*Pi/cosmologyParametersClosed%HubbleConstant(hubbleUnitsTime)
  ! Evaluate ages for matter + cosmological constant universes.
  call Unit_Tests_Begin_Group("Matter + Cosmological Constant")
  do iExpansion=1,size(redshift)
     expansionFactor=cosmologyFunctionsEdS%expansionFactorFromRedshift(redshift       (iExpansion))
     age            =cosmologyFunctionsEdS%cosmicTime                 (expansionFactor            )
     write (message,'(a,f6.1,a)') "cosmic age: EdS                       [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageEdS                 (iExpansion),relTol=1.0d-3)
     age=cosmologyFunctionsOpen                %cosmicTime(expansionFactor)
     write (message,'(a,f6.1,a)') "cosmic age: open                      [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageOpen                (iExpansion),relTol=1.0d-3)
     age=cosmologyFunctionsCosmologicalConstant%cosmicTime(expansionFactor)
     write (message,'(a,f6.1,a)') "cosmic age: cosmological constant     [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageCosmologicalConstant(iExpansion),relTol=1.0d-3)
  end do
  call Unit_Tests_End_Group()
  ! Evaluate ages for a dark energy Universe. We consider a few simple cases for which external solutions are available:
  !  * A cosmological constant;
  !  * A closed universe with no dark energy (to test the solving of re-collapsing universes);
  !  * A flat universe with equation of state -1/3 which behaves like the equivalent open universe.
  call Unit_Tests_Begin_Group("Matter + Dark Energy")
  do iExpansion=1,size(redshift)
     expansionFactor         =cosmologyFunctionsDarkEnergyCosmologicalConstant%expansionFactorFromRedshift(redshift(iExpansion))
     age                     =cosmologyFunctionsDarkEnergyCosmologicalConstant%cosmicTime     (expansionFactor         )
     write (message,'(a,f6.1,a)') "cosmic age: cosmological constant     [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageCosmologicalConstant(iExpansion),relTol=1.0d-3)
     age                     =cosmologyFunctionsDarkEnergyClosed              %cosmicTime     (expansionFactor         )
     write (message,'(a,f6.1,a)') "cosmic age: closed                    [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageClosed              (iExpansion),relTol=1.0d-3)
     expansionFactorSymmetric=cosmologyFunctionsDarkEnergyClosed              %expansionFactor(2.0d0*timeTurnaround-age)
     write (message,'(a,f6.1,a)') "cosmic age: closed: collapse symmetry [z=",redshift(iExpansion),"]"
     call Assert(trim(message),expansionFactor,expansionFactorSymmetric,absTol=0.03d0)
     age                     =cosmologyFunctionsDarkEnergyOmegaMinusOneThird  %cosmicTime     (expansionFactor         )
     write (message,'(a,f6.1,a)') "cosmic age: w^{-1/3}                  [z=",redshift(iExpansion),"]"
     call Assert(trim(message),age,ageOpen                (iExpansion),relTol=1.0d-3)
  end do
  call Unit_Tests_End_Group()
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Cosmic_Age
