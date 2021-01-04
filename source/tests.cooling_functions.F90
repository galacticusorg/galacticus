!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a program which tests cooling function functionality.

program Test_Cooling_Functions
  !% Tests cooling function functionality.
  use :: Abundances_Structure            , only : abundances                             , metallicityTypeLinearByMassSolar
  use :: Chemical_Abundances_Structure   , only : chemicalAbundances                     , zeroChemicalAbundances
  use :: Chemical_States                 , only : chemicalState                          , chemicalStateClass
  use :: Cooling_Functions               , only : coolingFunction                        , coolingFunctionAtomicCIECloudy  , coolingFunctionCMBCompton, coolingFunctionClass
  use :: Cosmology_Functions             , only : cosmologyFunctions                     , cosmologyFunctionsClass
  use :: Galacticus_Display              , only : Galacticus_Verbosity_Level_Set         , verbosityStandard
  use :: ISO_Varying_String              , only : varying_string                         , assignment(=)
  use :: Input_Parameters                , only : inputParameters
  use :: Numerical_Constants_Astronomical, only : gigaYear
  use :: Numerical_Constants_Physical    , only : boltzmannsConstant
  use :: Numerical_Constants_Units       , only : ergs
  use :: Radiation_Fields                , only : radiationFieldCosmicMicrowaveBackground
  use :: Unit_Tests                      , only : Assert                                 , Unit_Tests_Begin_Group          , Unit_Tests_End_Group     , Unit_Tests_Finish
  implicit none
  type            (inputParameters                        ), target  :: testParameters
  class           (coolingFunctionClass                   ), pointer :: coolingFunction_
  class           (chemicalStateClass                     ), pointer :: chemicalState_
  class           (cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_
  type            (varying_string                         )          :: parameterFile
  type            (abundances                             )          :: gasAbundances
  type            (chemicalAbundances                     )          :: chemicalDensities
  type            (radiationFieldCosmicMicrowaveBackground)          :: radiation
  type            (coolingFunctionCMBCompton              )          :: coolingFunctionCMBCompton_
  type            (coolingFunctionAtomicCIECloudy)                   :: coolingFunctionAtomicCIECloudy_
  double precision                                                   :: numberDensityHydrogen          , temperature           , &
       &                                                                coolantSummation               , coolantCMBCompton     , &
       &                                                                timescaleCooling               , coolantAtomicCIECloudy

  ! Set verbosity level.
  call Galacticus_Verbosity_Level_Set(verbosityStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Cooling functions")
  ! Read in controlling parameters.
  parameterFile='testSuite/parameters/coolingFunctions.xml'
  testParameters=inputParameters(parameterFile)
  call testParameters%markGlobal()
  ! Get required objects.
  cosmologyFunctions_ => cosmologyFunctions()
  chemicalState_      => chemicalState     ()
  coolingFunction_    => coolingFunction   ()
  ! Define plasma conditions.
  numberDensityHydrogen=1.0d-4 ! cm^-3
  temperature          =1.0d+6 ! K
  chemicalDensities    =zeroChemicalAbundances
  radiation            =radiationFieldCosmicMicrowaveBackground(cosmologyFunctions_)
  call radiation    %       timeSet(1.0d0                                                  )
  call gasAbundances%metallicitySet(1.0d0,metallicityType= metallicityTypeLinearByMassSolar)
  ! Construct cooling functions.
  coolingFunctionCMBCompton_      =  coolingFunctionCMBCompton     (chemicalState_)
  coolingFunctionAtomicCIECloudy_ =  coolingFunctionAtomicCIECloudy(              )
  ! Summed cooling function should be twice the CMB Compton cooling function.
  coolantSummation       =  coolingFunction_          %coolingFunction                   (numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  coolantCMBCompton      =  coolingFunctionCMBCompton_%coolingFunction                   (numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  call Assert('summation'       ,coolantSummation,2.0d0*coolantCMBCompton,relTol=1.0d-6)
  ! Logarithmic derivatives of summed cooling function should equal those of the CMB Compton cooling function.
  coolantSummation       =  coolingFunction_          %coolingFunctionDensityLogSlope    (numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  coolantCMBCompton      =  coolingFunctionCMBCompton_%coolingFunctionDensityLogSlope    (numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  call Assert('density slope'    ,coolantSummation,      coolantCMBCompton,relTol=1.0d-6)
  ! Logarithmic derivatives of summed cooling function should equal those of the CMB Compton cooling function.
  coolantSummation       =  coolingFunction_          %coolingFunctionTemperatureLogSlope(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  coolantCMBCompton      =  coolingFunctionCMBCompton_%coolingFunctionTemperatureLogSlope(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  call Assert('temperature slope',coolantSummation,      coolantCMBCompton,relTol=1.0d-6)
  ! Begin repeatibility tests.
  call Unit_Tests_Begin_Group("Repeatibility")
  ! Compute timescale for Compton cooling off of CMB at the present day.
  call radiation%timeSet(13.8d0)
  coolantCMBCompton=+coolingFunctionCMBCompton_%coolingFunction(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  timescaleCooling =+1.5d0                 &
       &            *numberDensityHydrogen &
       &            *boltzmannsConstant    &
       &            *temperature           &
       &            /ergs                  &
       &            /coolantCMBCompton     &
       &            /gigaYear
  call Assert('CMB Compton cooling timescale at z=0',timescaleCooling,976.65342064763729d0,relTol=1.0d-6)
  ! Compute Cloudy cooling time.
  coolantAtomicCIECloudy=+coolingFunctionAtomicCIECloudy_%coolingFunction(numberDensityHydrogen,temperature,gasAbundances,chemicalDensities,radiation)
  call Assert('Cloudy CIE cooling function',coolantAtomicCIECloudy,1.4177100259738058d-30,relTol=1.0d-6)
  call Unit_Tests_End_Group       ()
  ! End unit tests.
  call Unit_Tests_End_Group       ()
  call Unit_Tests_Finish          ()
end program Test_Cooling_Functions
